#!/usr/bin/env python
"""Generate IF97 conformance and timing tables for the docs.

Emits Web/fluid_properties/IF97_conformance.rst.in, which the IF97 docs
page ``include``s.

Conformance protocol mirrors IAPWS G13-15 (Kunick et al., 2015), Tables
8-13: for each IF97 region we draw random :math:`(p, T)` samples
log-uniform in pressure and uniform in temperature, retain the ones the
IF97 region classifier assigns to that region, look up the enthalpy
:math:`h = h_\\mathrm{IF97}(p, T)` so the comparison can be expressed
with the published-style backward arguments :math:`(p, h)`, evaluate
each tested backend at :math:`(h, p)`, and report maximum and
root-mean-square deviation of the tested backend from IAPWS-IF97 for
:math:`T(p,h)`, :math:`v(p,h)`, :math:`s(p,h)`, :math:`w(p,h)`,
:math:`\\eta(p,h)`, and :math:`\\lambda(p,h)`.

Timing: per backend per property, mean wall time for ``PropsSI`` over
the same sample population, reported as ratio relative to IF97.

The script gracefully skips any backend whose factory string isn't
available (e.g. on a CoolProp build without the SVDSBTL backend
plumbed) — that way the docs build never breaks just because a
backend is missing.
"""
from __future__ import print_function

import math
import os
import random
import time

import CoolProp.CoolProp as CP

WEB_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_PATH = os.path.join(WEB_DIR, 'fluid_properties', 'IF97_conformance.rst.in')

# Truth source for the conformance comparison.
REFERENCE = 'IF97::Water'

# Backends to test against the reference.  Each entry is rendered as a
# block of tables; if the backend can't be instantiated the block is
# skipped (with a placeholder note in the RST).
TESTED = [
    ('SVDSBTL&IF97::Water', 'SVDSBTL with IF97 source'),
]

# Properties to compare.  Each entry: (CoolProp PropsSI key, math symbol,
# unit, kind, IAPWS perm-value, IAPWS perm-unit, G13-15 table id).
#   kind = 'rel'  -> relative deviation as %, reported in %
#   kind = 'abs_T'-> absolute deviation, reported in mK (with x1000)
#   kind = 'abs_s'-> absolute deviation, reported in 10^-3 J/(kg K)
PROPERTIES = [
    ('T', 'T',         'K',        'abs_T', 25.0,    'mK',                       8),
    ('D', r'\rho',     'kg/m^3',   'rel',   0.001,   r'\%',                      9),
    ('S', 's',         'J/(kg K)', 'abs_s', 1.0,     r'10^{-3}\ \mathrm{J/(kg\,K)}', 10),
    ('A', 'w',         'm/s',      'rel',   0.001,   r'\%',                     11),
    ('V', r'\eta',     'Pa s',     'rel',   0.001,   r'\%',                     12),
    ('L', r'\lambda',  'W/(m K)',  'rel',   0.001,   r'\%',                     13),
]

# Property keys we time (must be a subset of the comparison properties).
TIMING_PROPS = ['T', 'D', 'A', 'V']

# Per-region (T, p) sampling box [K, MPa]. Bounds come from IAPWS-IF97.
REGION_BOXES = {
    1: dict(T=(273.15, 623.15), p=(1.0e-3, 100.0)),
    2: dict(T=(273.15, 1073.15), p=(1.0e-3, 100.0)),
    3: dict(T=(623.15, 863.15), p=(16.5292, 100.0)),
    5: dict(T=(1073.15, 2273.15), p=(1.0e-3, 50.0)),
}
SAMPLES_PER_REGION = 2000
SEED = 0xC001DAD


# ---------------------------------------------------------------------------
# IF97 region classifier in (p, T). Bounds + B23 boundary from IF97 Rev.
# ---------------------------------------------------------------------------

def _p_B23(T_K):
    """B23 boundary curve p(T) in MPa."""
    n3 = 0.10192970039326e-2
    n4 = 0.57254459862746e3
    n5 = 0.1391883776670e2
    return n3 * (T_K - n4) ** 2 + n5


def _p_sat_MPa(T_K):
    try:
        return CP.PropsSI('P', 'T', T_K, 'Q', 0, REFERENCE) / 1e6
    except Exception:
        return None


def classify_region(T_K, p_MPa):
    """Return IF97 region index in {1, 2, 3, 5} or None if outside."""
    if not (273.15 <= T_K <= 2273.15):
        return None
    if p_MPa <= 0.0 or p_MPa > 100.0:
        return None
    if T_K > 1073.15:
        return 5 if p_MPa <= 50.0 else None
    if T_K <= 623.15:
        psat = _p_sat_MPa(T_K)
        if psat is None:
            return None
        return 1 if p_MPa > psat else 2
    if T_K <= 863.15 and p_MPa > _p_B23(T_K):
        return 3
    return 2


# ---------------------------------------------------------------------------
# Deviation accumulator + property comparison
# ---------------------------------------------------------------------------

class DevAcc(object):
    __slots__ = ('n', 'max_abs', 'sum_sq')

    def __init__(self):
        self.n = 0
        self.max_abs = 0.0
        self.sum_sq = 0.0

    def add(self, delta):
        if delta is None or not math.isfinite(delta):
            return
        a = abs(delta)
        if a > self.max_abs:
            self.max_abs = a
        self.sum_sq += a * a
        self.n += 1

    @property
    def rms(self):
        if self.n == 0:
            return float('nan')
        return math.sqrt(self.sum_sq / self.n)


def _deviation(ref, test, kind):
    if ref is None or test is None:
        return None
    if not math.isfinite(ref) or not math.isfinite(test):
        return None
    if kind == 'rel':
        if abs(ref) < 1e-30:
            return None
        # For density we report deviation of v = 1/rho; the table is
        # built on v so flip the sign convention to match G13-15.
        return 100.0 * (test - ref) / ref
    if kind == 'abs_T':
        return 1000.0 * (test - ref)  # K -> mK
    if kind == 'abs_s':
        return 1000.0 * (test - ref)  # J/(kg K) -> 10^-3 J/(kg K)
    return test - ref


def _safe_propssi(key, h, p, backend):
    try:
        return CP.PropsSI(key, 'H', h, 'P', p, backend)
    except Exception:
        return None


def _safe_truth_h(T, p):
    try:
        return CP.PropsSI('H', 'T', T, 'P', p, REFERENCE)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Conformance sweep
# ---------------------------------------------------------------------------

def backend_available(factory_str):
    try:
        CP.PropsSI('T', 'P', 1e5, 'Q', 0, factory_str)
        return True
    except Exception:
        try:
            CP.PropsSI('D', 'T', 500, 'P', 1e6, factory_str)
            return True
        except Exception:
            return False


def run_conformance(backend):
    """Return (results, counts, samples).

    results[region][prop_key] = DevAcc of (backend - reference) in the
    table's reported unit. counts[region] = in-region accepted samples.
    samples[region] = list of dicts {'h', 'p', 'devs': {prop: |dev|}} —
    raw per-sample record used by the failure-point figures.
    """
    rng = random.Random(SEED)
    results = {r: {p[0]: DevAcc() for p in PROPERTIES} for r in REGION_BOXES}
    counts = {r: 0 for r in REGION_BOXES}
    samples = {r: [] for r in REGION_BOXES}

    prop_keys = [p[0] for p in PROPERTIES]
    kinds = {p[0]: p[3] for p in PROPERTIES}

    for region, box in REGION_BOXES.items():
        T_lo, T_hi = box['T']
        p_lo, p_hi = box['p']
        for _ in range(SAMPLES_PER_REGION):
            T = rng.uniform(T_lo, T_hi)
            p_MPa = math.exp(rng.uniform(math.log(p_lo), math.log(p_hi)))
            if classify_region(T, p_MPa) != region:
                continue
            p_Pa = p_MPa * 1e6
            h = _safe_truth_h(T, p_Pa)
            if h is None or not math.isfinite(h):
                continue
            # Truth values for all comparison properties.
            ref_vals = {k: _safe_propssi(k, h, p_Pa, REFERENCE) for k in prop_keys}
            if any(v is None for v in ref_vals.values()):
                continue
            test_vals = {k: _safe_propssi(k, h, p_Pa, backend) for k in prop_keys}
            if any(v is None for v in test_vals.values()):
                continue
            counts[region] += 1
            per_sample_devs = {}
            for k in prop_keys:
                # v = 1/rho — flip so the table reports v deviation.
                if k == 'D':
                    ref_v = 1.0 / ref_vals[k]
                    test_v = 1.0 / test_vals[k]
                    dev = _deviation(ref_v, test_v, kinds[k])
                else:
                    dev = _deviation(ref_vals[k], test_vals[k], kinds[k])
                results[region][k].add(dev)
                per_sample_devs[k] = abs(dev) if (dev is not None and math.isfinite(dev)) else None
            samples[region].append({'h': h, 'p': p_Pa, 'devs': per_sample_devs})
    return results, counts, samples


def write_failure_figures(backend, samples_per_region, out_dir):
    """Emit one PNG per region showing sample population in PH coords
    with budget-violating points overlaid darker.

    Returns dict region → relative-path of the emitted PNG (relative to
    Web/fluid_properties/ for inclusion in the rst.in).  When the
    matplotlib import fails (headless CI without matplotlib), returns
    an empty dict and the rst renderer omits the figures section.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return {}

    # Slug for the backend name so the filename is filesystem-safe.
    slug = ''.join(c if c.isalnum() else '_' for c in backend)
    perms = {p[0]: p[4] for p in PROPERTIES}
    symbols = {p[0]: p[1] for p in PROPERTIES}
    paths = {}

    for region in sorted(samples_per_region):
        samples = samples_per_region[region]
        if not samples:
            continue
        # 2x3 grid: T, v, s, w, η, λ
        fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), sharex=True, sharey=True)
        h_all = [s['h'] / 1e3 for s in samples]   # h in kJ/kg
        p_all = [s['p'] / 1e6 for s in samples]   # p in MPa
        prop_order = [p[0] for p in PROPERTIES]
        for ax, prop_key in zip(axes.flat, prop_order):
            ax.scatter(h_all, p_all, c='lightgray', s=4, edgecolors='none', label='in-region samples')
            failing_h, failing_p = [], []
            for s in samples:
                d = s['devs'].get(prop_key)
                if d is not None and d > perms[prop_key]:
                    failing_h.append(s['h'] / 1e3)
                    failing_p.append(s['p'] / 1e6)
            if failing_h:
                ax.scatter(failing_h, failing_p, c='crimson', s=8, edgecolors='none', label='over budget')
            ax.set_yscale('log')
            ax.set_title(r'$|\Delta {0}|>\,{1:g}$ ({2}/{3})'.format(symbols[prop_key], perms[prop_key], len(failing_h), len(samples)))
            ax.grid(True, which='both', alpha=0.3)
        for ax in axes[-1]:
            ax.set_xlabel(r'$h$ [kJ/kg]')
        for ax in axes[:, 0]:
            ax.set_ylabel(r'$p$ [MPa]')
        fig.suptitle(r'IF97 Region {0} — `{1}` budget-violating samples ({2} drawn)'.format(region, backend, len(samples)))
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        out_path = os.path.join(out_dir, 'IF97_conformance_fails_{0}_R{1}.png'.format(slug, region))
        fig.savefig(out_path, dpi=110)
        plt.close(fig)
        paths[region] = os.path.basename(out_path)
    return paths


# ---------------------------------------------------------------------------
# Timing sweep
# ---------------------------------------------------------------------------

def run_timing(backends, rng_seed=0xCAFE):
    """Mean PropsSI time per (backend, property) over N timing calls.

    Returns timings[backend][prop_key] = ns/call.
    """
    rng = random.Random(rng_seed)
    # Build a single (T, p) sample population using the same per-region
    # sampling + classifier as run_conformance(): draw from each region
    # box, accept only points the IF97 region atlas assigns to that
    # region.  Then probe every backend in the timing set at each
    # candidate (with every property we'll time) and KEEP only samples
    # that succeed everywhere.  This is the right place to filter out
    # SVDSBTL's known rank-truncation gap near the critical singularity
    # (R3, !mayfail in the conformance tests): including those points
    # in the timed loop would either crash or require try/except inside
    # the inner per-call measurement.  Pre-filtering keeps the timed
    # loop a plain for-loop and timing data is the wall-time mean over
    # the SHARED-evaluable population for fair backend-to-backend
    # comparison.
    T_arr = []
    p_arr = []
    region_of = []
    for region, box in REGION_BOXES.items():
        T_lo, T_hi = box['T']
        p_lo, p_hi = box['p']
        attempts = 0
        accepted = 0
        target = SAMPLES_PER_REGION
        while accepted < target and attempts < 10 * target:
            T = rng.uniform(T_lo, T_hi)
            p_MPa = math.exp(rng.uniform(math.log(p_lo), math.log(p_hi)))
            attempts += 1
            if classify_region(T, p_MPa) != region:
                continue
            p_Pa = p_MPa * 1e6
            ok = True
            for backend in backends:
                for prop in TIMING_PROPS:
                    try:
                        v = CP.PropsSI(prop, 'T', T, 'P', p_Pa, backend)
                    except Exception:
                        ok = False
                        break
                    if not math.isfinite(v):
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue
            T_arr.append(T)
            p_arr.append(p_Pa)
            region_of.append(region)
            accepted += 1
    N = len(T_arr)
    # Pre-filter may reject everything if the backend set has zero
    # joint-domain overlap on the IF97 region boxes — leave timings as
    # NaN rather than divide by zero and abort the docs build.
    if N == 0:
        return {backend: {prop: float('nan') for prop in TIMING_PROPS}
                for backend in backends}
    timings = {}
    for backend in backends:
        timings[backend] = {}
        for prop in TIMING_PROPS:
            for k in range(min(64, N)):  # warmup
                CP.PropsSI(prop, 'T', T_arr[k], 'P', p_arr[k], backend)
            t0 = time.perf_counter()
            for k in range(N):
                CP.PropsSI(prop, 'T', T_arr[k], 'P', p_arr[k], backend)
            t1 = time.perf_counter()
            timings[backend][prop] = (t1 - t0) / N * 1e9
    return timings


# ---------------------------------------------------------------------------
# RST rendering
# ---------------------------------------------------------------------------

def _fmt_dev(value):
    if value is None or not math.isfinite(value):
        return '--'
    return '{:.3g}'.format(value)


def _fmt_perm(p):
    if p >= 1.0:
        return '{:g}'.format(p)
    return '{:.3g}'.format(p)


def _fmt_ns(ns):
    if not math.isfinite(ns):
        return '--'
    if ns >= 1e6:
        return '{:.2f} ms'.format(ns / 1e6)
    if ns >= 1e3:
        return '{:.2f} us'.format(ns / 1e3)
    return '{:.0f} ns'.format(ns)


def render_rst(results_per_backend, counts_per_backend, timings, figure_paths_per_backend=None):
    if figure_paths_per_backend is None:
        figure_paths_per_backend = {}
    lines = []
    lines.append('')
    lines.append('.. _IF97-Conformance:')
    lines.append('')
    lines.append('IF97 Conformance and Timing')
    lines.append('---------------------------')
    lines.append('')
    lines.append(
        'The tables in this section are regenerated by '
        '``Web/scripts/fluid_properties.IF97Conformance.py`` on every '
        'docs build. The protocol mirrors '
        ':ref:`IAPWS G13-15 Tables 8-13 <IAPWS-IF97>` '
        '(Kunick et al., 2015): for each IF97 region we draw '
        '{n} random :math:`(p, T)` samples log-uniform in pressure and '
        'uniform in temperature, classify them into the IF97 region atlas, '
        'compute :math:`h = h_\\mathrm{{IF97}}(p, T)`, evaluate each '
        'tested backend at :math:`(h, p)`, and tabulate the maximum and '
        'root-mean-square deviation from IAPWS-IF97 of the six properties '
        'in G13-15: :math:`T(p,h)`, :math:`v(p,h)`, :math:`s(p,h)`, '
        ':math:`w(p,h)`, :math:`\\eta(p,h)`, :math:`\\lambda(p,h)`.'
        .format(n=SAMPLES_PER_REGION))
    lines.append('')

    for backend, _label in TESTED:
        if backend not in results_per_backend:
            lines.append('Backend ``{}``: not available on this build.'.format(backend))
            lines.append('')
            continue
        results = results_per_backend[backend]
        fig_paths = figure_paths_per_backend.get(backend, {})

        section_title = 'Deviation of ``{}`` from IAPWS-IF97'.format(backend)
        lines.append(section_title)
        lines.append('~' * max(60, len(section_title)))
        lines.append('')

        # Failure-point figures (one per region) showing the
        # population in PH coords with budget-violating samples
        # overlaid in red.  Lets the reader see WHERE the SVDSBTL
        # truncation residual exceeds the IAPWS perm budget, not
        # just the bulk statistics.
        if fig_paths:
            lines.append('*Where the failures cluster (PH coords; light grey = all in-region samples, '
                         'red = budget-violating samples for that property; one panel per property; '
                         'pressure axis is log-scaled):*')
            lines.append('')
            for region in sorted(fig_paths):
                lines.append('.. figure:: {0}'.format(fig_paths[region]))
                lines.append('   :align: center')
                lines.append('   :width: 90%')
                lines.append('')
                lines.append('   IF97 Region {0} — {1}'.format(region, backend))
                lines.append('')

        for pk, sym, _unit, kind, perm, perm_unit, table_id in PROPERTIES:
            if kind == 'rel':
                unit_str = r'\%'
            elif kind == 'abs_T':
                unit_str = 'mK'
            elif kind == 'abs_s':
                unit_str = r'10^{-3}\ \mathrm{J/(kg\,K)}'
            else:
                unit_str = perm_unit

            prop_sym = sym if pk != 'D' else 'v'  # density → specific volume
            lines.append('')
            lines.append(
                '*G13-15 Table {tid}* — :math:`{psym}(p,h)`, '
                'permissible :math:`|\\Delta {psym}|_\\mathrm{{perm}} = {perm}\\ {pu}`.'
                .format(tid=table_id, psym=prop_sym, perm=_fmt_perm(perm), pu=unit_str))
            lines.append('')
            lines.append('.. list-table::')
            lines.append('   :header-rows: 1')
            lines.append('   :widths: 12 18 22 22')
            lines.append('')
            lines.append('   * - IF97 Region')
            lines.append('     - in-region samples')
            lines.append('     - :math:`|\\Delta {0}|_{{\\max}}` [{1}]'.format(prop_sym, unit_str))
            lines.append('     - :math:`(\\Delta {0})_{{\\mathrm{{RMS}}}}` [{1}]'.format(prop_sym, unit_str))
            for region in sorted(REGION_BOXES):
                acc = results[region][pk]
                lines.append('   * - {}'.format(region))
                lines.append('     - {}'.format(acc.n))
                lines.append('     - {}'.format(_fmt_dev(acc.max_abs)))
                lines.append('     - {}'.format(_fmt_dev(acc.rms)))
        lines.append('')

    # ------------------------------------------------------------------
    # Timing table
    # ------------------------------------------------------------------
    if timings:
        lines.append('Backend timing (PT inputs)')
        lines.append('~~~~~~~~~~~~~~~~~~~~~~~~~~')
        lines.append('')
        lines.append(
            'Mean wall time per ``PropsSI`` call for property look-ups at '
            'random :math:`(T, p)` samples drawn from the same per-IF97-region '
            'population the conformance sweep uses above (stratified across '
            'R1, R2, R3, R5; saturation cells filtered out by the IF97 '
            'region classifier so every timed call is single-phase, and '
            'samples that any tested backend cannot evaluate — e.g. cells '
            'inside SVDSBTL\'s known rank-truncation gap near the critical '
            'singularity — are pre-rejected during population so the timed '
            'loop is a plain for-loop). '
            'The ratio column is :math:`t_{\\mathrm{backend}} / '
            't_{\\mathrm{IF97}}`. Lower is faster than IF97. '
            'Absolute timings depend on hardware; the **ratios** are the '
            'load-bearing quantity. CI rebuilds these on a GitHub-hosted '
            'runner.')
        lines.append('')
        lines.append('.. list-table::')
        lines.append('   :header-rows: 1')
        lines.append('   :widths: 22 12 12 12 12 14')
        lines.append('')
        header = ['   * - Backend']
        for prop in TIMING_PROPS:
            header.append('     - {}'.format(prop))
        header.append('     - mean ratio vs IF97')
        lines.extend(header)

        baseline = timings.get(REFERENCE, {})
        ordered = [REFERENCE] + [b for b, _ in TESTED if b in timings]
        for backend in ordered:
            bt = timings.get(backend, {})
            if not bt:
                continue
            row = ['   * - ``{}``'.format(backend)]
            ratios = []
            for prop in TIMING_PROPS:
                row.append('     - {}'.format(_fmt_ns(bt.get(prop, float('nan')))))
                base = baseline.get(prop, float('nan'))
                this = bt.get(prop, float('nan'))
                if math.isfinite(base) and base > 0 and math.isfinite(this):
                    ratios.append(this / base)
            row.append('     - {:.2f}'.format(sum(ratios) / len(ratios)) if ratios else '     - --')
            lines.extend(row)
        lines.append('')
    return '\n'.join(lines) + '\n'


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    print('IF97 conformance script: starting')
    results_per_backend = {}
    counts_per_backend = {}
    figure_paths_per_backend = {}
    available = [REFERENCE]
    for backend, label in TESTED:
        if not backend_available(backend):
            print('  {} unavailable — skipping'.format(backend))
            continue
        print('  sweeping {} ({})'.format(backend, label))
        results, counts, samples = run_conformance(backend)
        results_per_backend[backend] = results
        counts_per_backend[backend] = counts
        available.append(backend)
        print('    in-region counts: {}'.format(counts))
        # Emit the per-region failure-point PH figures alongside the
        # .rst.in include.  Failures take their basename so the
        # rst-relative path is just the file name.
        fig_dir = os.path.dirname(OUT_PATH)
        paths = write_failure_figures(backend, samples, fig_dir)
        if paths:
            figure_paths_per_backend[backend] = paths
            print('    wrote {} failure figure(s)'.format(len(paths)))

    print('  timing pass over {}'.format(available))
    timings = run_timing(available)

    rst = render_rst(results_per_backend, counts_per_backend, timings, figure_paths_per_backend)
    with open(OUT_PATH, 'w') as fp:
        fp.write(rst)
    print('Wrote {}'.format(OUT_PATH))


if __name__ == '__main__':
    main()
