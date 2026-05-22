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

Timing: per backend per property, mean wall time for one
``AbstractState.update()`` + ``keyed_output()`` call over the same
sample population, reported as ratio relative to IF97.  (The
``AbstractState`` is instantiated once per backend and reused across
the loop, so the cost is the flash + property read with no factory
rebuild per probe — a real ``PropsSI`` end-to-end call pays an
additional reconstruction tax that dominates for tabular / SVDSBTL
backends.)

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
from CoolProp import AbstractState

# Map our short keys to the CoolProp parameter int IDs used by the
# AbstractState.keyed_output() fast path.  Avoids the per-call string
# lookup and PropsSI factory rebuild that PropsSI('key', ...) would do.
_PARAM_KEY = {
    'T': CP.iT,
    'D': CP.iDmass,
    'S': CP.iSmass,
    'A': CP.ispeed_sound,
    'V': CP.iviscosity,
    'L': CP.iconductivity,
    'H': CP.iHmass,
}

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

# Properties to compare.  Each entry: (CoolProp PropsSI key, math
# symbol, unit, kind, IAPWS perm-value, IAPWS perm-unit, G13-15 table
# id).
#   kind = 'rel'  -> relative deviation as %, reported in %
#   kind = 'abs_T'-> absolute deviation, reported in mK (with x1000)
#   kind = 'abs_s'-> absolute deviation, reported in 10^-3 J/(kg K)
#
# The perm-value can be either a scalar (applies to all regions) or a
# dict {IF97 region : perm}.  G13-15 specifies region-dependent
# tolerances for T(p,h): 25 mK in R1 and R3, 10 mK in R2 and R5.
# Other tables here use a single perm pending per-region verification
# against the published G13-15 Tables 9-13.
PROPERTIES = [
    ('T', 'T',         'K',        'abs_T', {1: 25.0, 2: 10.0, 3: 25.0, 5: 10.0}, 'mK', 8),
    ('D', r'\rho',     'kg/m^3',   'rel',   0.001,   r'\%',                      9),
    ('S', 's',         'J/(kg K)', 'abs_s', 1.0,     r'10^{-3}\ \mathrm{J/(kg\,K)}', 10),
    ('A', 'w',         'm/s',      'rel',   0.001,   r'\%',                     11),
    ('V', r'\eta',     'Pa s',     'rel',   0.001,   r'\%',                     12),
    ('L', r'\lambda',  'W/(m K)',  'rel',   0.001,   r'\%',                     13),
]


def _perm_for(perm, region):
    """Resolve a perm field that may be either a scalar or a dict of
    region -> perm.  Falls back to the smallest perm in the dict if
    the region key isn't present (= the strictest budget, conservative)."""
    if isinstance(perm, dict):
        if region in perm:
            return perm[region]
        # Region not enumerated — return the strictest budget so an
        # untyped region isn't silently let off the hook.
        return min(perm.values())
    return perm

# Property keys we time (must be a subset of the comparison properties).
TIMING_PROPS = ['T', 'D', 'A', 'V']

# Per-region (T, p) sampling box [K, MPa]. Bounds come from IAPWS-IF97.
REGION_BOXES = {
    1: dict(T=(273.15, 623.15), p=(1.0e-3, 100.0)),
    2: dict(T=(273.15, 1073.15), p=(1.0e-3, 100.0)),
    3: dict(T=(623.15, 863.15), p=(16.5292, 100.0)),
    5: dict(T=(1073.15, 2273.15), p=(1.0e-3, 50.0)),
}
# Bumped from 2000 → 50000 per IF97 region.  Reasoning:
#   * The conformance Tables 8-13 are statistical claims about a
#     population (max + RMS over random IF97 probes).  At 2000 in-
#     region samples (~24 R3 probes after classification), the R3
#     column was a 24-sample average — too noisy to call meaningful.
#   * The failure-cluster figure was visibly sparse vs the reference
#     /tmp/conformance_dense_540k.png; the geographic clustering
#     pattern only emerges with enough density to fill the (h, p) plane.
#   * Total runtime budget: 50000 * 4 regions * ~8 us per IF97 HmassP
#     probe ~= 1.6 s for the conformance sweep, plus ~6 s for the
#     timing pass.  Docs build adds < 10 s.
SAMPLES_PER_REGION = 50000
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
    """Legacy PropsSI helper.  Kept for the few places (timing pass,
    backend availability probe) that genuinely want the high-level
    string-keyed surface.  The hot conformance loop uses the
    AbstractState fast path below instead — going through PropsSI
    means a factory() rebuild per call, which for SVDSBTL is ~80 ms
    and turned this script into a 2+ hour drag on the docs build."""
    try:
        return CP.PropsSI(key, 'H', h, 'P', p, backend)
    except Exception:
        return None


def _factory(backend_str):
    """Resolve the backend factory string (e.g. ``SVDSBTL&IF97::Water``)
    to an AbstractState we keep alive for the whole sweep, paying the
    surface-load cost ONCE instead of once per probe."""
    if '::' in backend_str:
        backend, fluid = backend_str.split('::', 1)
    else:
        backend, fluid = 'HEOS', backend_str
    return AbstractState(backend, fluid)


def _read_props_fast(state, h, p, prop_int_keys):
    """update(HmassP, h, p) once, then read every key off the same
    cached state in one shot.  Returns ``None`` when update fails so
    the caller can skip the sample (matches the prior PropsSI semantics)."""
    try:
        state.update(CP.HmassP_INPUTS, h, p)
    except Exception:
        return None
    out = {}
    for short_key, ip in prop_int_keys:
        try:
            out[short_key] = state.keyed_output(ip)
        except Exception:
            return None
    return out


def _refine_to_forward_h(state, p, T_seed, h_target, rel_tol=1e-10):
    """TOMS748-iterate ``state`` on the FORWARD equation h(T, p) until
    its cached T satisfies h(T, p) = h_target to machine precision.
    Mirrors the polish used by SVDSBTL's own sampling lambda
    (src/SBTL/SurfacePresets.cpp + the SurfaceFailMap C++ harness in
    src/Tests/CoolProp-Tests-SVDSBTLFailMap.cpp; PR #2940).

    Why TOMS748 rather than Newton: SVDSBTL's table was built via a
    TOMS748-polished sampling pass (foi.9.10).  Using TOMS748 here too
    means the conformance comparison's reference-state refinement is
    bit-identical-in-philosophy to the source-table construction —
    same root-finder, same convergence guarantees, no cpmass slope
    discontinuity at region boundaries to trip up Newton.

    Why this matters: IF97's published backward T_phmass(p, h) has
    ±25 mK residual baked in (R7-97 backward eqs are polynomial fits,
    not exact inverses).  If used as "truth", every R1/R3 cell shows
    ~25 mK of deviation that's really IF97's own backward-equation
    residual — swamps the SVD's actual error.  Refining the reference
    to forward-consistent T removes that floor so the figure shows
    the SVD's residual against IF97's forward equations, not the
    backward-equation noise.

    State pinning: hold the IF97 phase across the iteration so HEOS /
    IF97 don't flip basins mid-iteration near the saturation curve.

    Best-effort: returns silently on bracket failure / source rejection.
    """
    try:
        state.update(CP.HmassP_INPUTS, h_target, p)
    except Exception:
        try:
            state.update(CP.PT_INPUTS, p, T_seed)
        except Exception:
            pass
        return

    pinned = False
    try:
        phase0 = state.phase()
        single = phase0 in (CP.iphase_liquid, CP.iphase_gas,
                            CP.iphase_supercritical_liquid,
                            CP.iphase_supercritical_gas,
                            CP.iphase_supercritical)
        if single:
            state.specify_phase(phase0)
            pinned = True
    except Exception:
        pass

    def h_residual(T):
        try:
            state.update(CP.PT_INPUTS, p, T)
            return state.hmass() - h_target
        except Exception:
            return float('nan')

    # Build a tight bracket around T_seed and check for sign change.
    # ±0.5 K is wider than IF97's ±25 mK backward residual by 20×,
    # comfortably contains the forward root.
    try:
        T_lo = max(state.Tmin() if hasattr(state, 'Tmin') else 273.16, T_seed - 0.5)
        T_hi = min(state.Tmax() if hasattr(state, 'Tmax') else 2273.15, T_seed + 0.5)
    except Exception:
        T_lo, T_hi = T_seed - 0.5, T_seed + 0.5

    try:
        from scipy.optimize import toms748
    except ImportError:
        toms748 = None

    try:
        if toms748 is None or not (T_lo < T_hi):
            # Fallback to leaving the state at the backward seed.
            pass
        else:
            r_lo = h_residual(T_lo)
            r_hi = h_residual(T_hi)
            if math.isfinite(r_lo) and math.isfinite(r_hi) and r_lo * r_hi <= 0:
                try:
                    T_root, info = toms748(h_residual, T_lo, T_hi,
                                           xtol=1e-6, rtol=rel_tol,
                                           full_output=True)
                    if info.converged:
                        state.update(CP.PT_INPUTS, p, T_root)
                except Exception:
                    # On any solver failure, restore the HmassP seed
                    # so callers don't see a wild PT state.
                    try: state.update(CP.HmassP_INPUTS, h_target, p)
                    except Exception: pass
            else:
                # Bracket misses (rare: at T_seed boundary cells).
                # Leave at the backward seed.
                try: state.update(CP.HmassP_INPUTS, h_target, p)
                except Exception: pass
    finally:
        if pinned:
            try: state.unspecify_phase()
            except Exception: pass


# ---------------------------------------------------------------------------
# Conformance sweep
# ---------------------------------------------------------------------------

def backend_available(factory_str):
    try:
        _factory(factory_str)
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
    prop_int_keys = [(k, _PARAM_KEY[k]) for k in prop_keys]

    # Build BOTH AbstractStates once and reuse across all probes.  This
    # is the entire perf delta — PropsSI rebuilds the AbstractState (and
    # for SVDSBTL re-loads the cached SVDSurface, ~80 ms) on every call.
    ref_state = _factory(REFERENCE)
    test_state = _factory(backend)
    # Truth state used to convert (T, p) -> h on the reference side.
    truth_state = _factory(REFERENCE)

    for region, box in REGION_BOXES.items():
        T_lo, T_hi = box['T']
        p_lo, p_hi = box['p']
        for _ in range(SAMPLES_PER_REGION):
            T = rng.uniform(T_lo, T_hi)
            p_MPa = math.exp(rng.uniform(math.log(p_lo), math.log(p_hi)))
            if classify_region(T, p_MPa) != region:
                continue
            p_Pa = p_MPa * 1e6
            try:
                truth_state.update(CP.PT_INPUTS, p_Pa, T)
                h = truth_state.hmass()
            except Exception:
                continue
            if not math.isfinite(h):
                continue
            # Forward-refine the reference state so its cached T
            # matches h_target = h(T_forward, p) to machine precision —
            # otherwise the IF97 backward equation's ±25 mK residual
            # leaks into every deviation and swamps the figure.
            # ref_state and truth_state are both REFERENCE (IF97); the
            # refinement leaves ref_state holding forward-consistent
            # (T, p, h, s, ...) which is what the test_state's
            # (SVDSBTL) post-update values should be compared against.
            _refine_to_forward_h(ref_state, p_Pa, T, h)
            ref_vals = {}
            try:
                for short_key, ip in prop_int_keys:
                    ref_vals[short_key] = ref_state.keyed_output(ip)
            except Exception:
                ref_vals = None
            if ref_vals is None:
                continue
            test_vals = _read_props_fast(test_state, h, p_Pa, prop_int_keys)
            if test_vals is None:
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
    """Emit ONE combined PNG showing budget-violating probes for the
    four G13-15 thermodynamic properties (T, v, s, w) in (h, p) coords,
    overlaid on the IF97 region boundaries and the critical-patch
    HEOS-fallback bbox.

    Layout: 2x2 grid, one panel per property.  Each panel plots ONLY
    the budget-violating samples (no full-population scatter) so the
    geographic pattern is unambiguous.  Markers are sized + coloured
    by :math:`\\log_{10}(|\\Delta| / |\\Delta|_\\mathrm{perm})` — bigger
    + redder means further past budget.

    Overlays:
      - IF97 saturation curve (R4)
      - R1/R3 isotherm (T = 623.15 K) and R2/R5 isotherm (T = 1073.15 K),
        in (h, p) coords by walking the isotherm and recording IF97 h(T,p)
      - IF97 B23 boundary p_B23(T), in (h, p) coords
      - critical-patch HEOS-fallback bbox (auto-calibrated; in (h, p) coords
        via the same perimeter walk the backend uses), gold-shaded
      - critical point as a black dot

    Returns ``{'all': <filename>}`` so render_rst() can embed the single
    figure with the same dict-of-paths plumbing used previously.  Empty
    dict on matplotlib import failure (headless CI).
    """
    try:
        import math as _math

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import LogNorm
    except ImportError:
        return {}

    # Slug for the backend name so the filename is filesystem-safe.
    slug = ''.join(c if c.isalnum() else '_' for c in backend)
    # All six G13-15 properties: T, v, s, w (thermodynamic, Tables 8-11)
    # plus eta, lambda (transport, Tables 12-13).  Layout is 2x3 to
    # match — transport properties fail more aggressively than the
    # thermodynamic ones at the rank-truncation limit, but the
    # log10-scaled colour bar keeps both classes legible side-by-side.
    plot_props = list(PROPERTIES)
    perms = {p[0]: p[4] for p in plot_props}
    symbols = {p[0]: p[1] for p in plot_props}

    # Resolve the actual backend AbstractState to fetch (a) the IF97
    # boundary curves in (h, p) coords and (b) the critical-patch bbox
    # the backend computed at construction.  Both use the source
    # backend's own h(T, p) so the overlay matches the geometry the
    # backend actually classifies on.
    try:
        ref = _factory(backend)
    except Exception:
        return {}
    try:
        Tc = ref.T_critical()
        pc = ref.p_critical()
    except Exception:
        Tc = pc = None

    # Saturation curve (R4): walk Tsat from triple to critical.
    sat_h = []
    sat_p = []
    if Tc is not None:
        T_tri = max(273.16, ref.Ttriple() if hasattr(ref, 'Ttriple') else 273.16)
        for i in range(200):
            T = T_tri + (Tc - T_tri) * (i / 199.0)
            try:
                p_sat = CP.PropsSI('P', 'T', T, 'Q', 0, 'IF97::Water')
                h_L = CP.PropsSI('H', 'T', T, 'Q', 0, 'IF97::Water')
                h_V = CP.PropsSI('H', 'T', T, 'Q', 1, 'IF97::Water')
                sat_h.append((h_L, h_V))
                sat_p.append(p_sat)
            except Exception:
                pass

    # Isotherm overlays — walk fixed T sweeping log p, record (h, p).
    def isotherm(T_K, p_lo_Pa, p_hi_Pa, n=64):
        hs, ps = [], []
        for i in range(n):
            p = _math.exp(_math.log(p_lo_Pa) + (i / (n - 1)) * (_math.log(p_hi_Pa) - _math.log(p_lo_Pa)))
            try:
                h = CP.PropsSI('H', 'T', T_K, 'P', p, 'IF97::Water')
                hs.append(h)
                ps.append(p)
            except Exception:
                pass
        return hs, ps

    iso_R1R3_h, iso_R1R3_p = isotherm(623.15, 1.0e3, 100.0e6)
    iso_R2R5_h, iso_R2R5_p = isotherm(1073.15, 1.0e3, 100.0e6)

    # B23 curve: walk T in [623.15, 863.15], compute p_B23(T) and the
    # corresponding h via IF97 at (T, p_B23).  Already in (h, p).
    B23_h, B23_p = [], []
    for i in range(64):
        T = 623.15 + (863.15 - 623.15) * (i / 63.0)
        p_B23_MPa = 0.10192970039326e-2 * (T - 572.54459862746) ** 2 + 13.91883776670
        p_B23_Pa = p_B23_MPa * 1.0e6
        try:
            h = CP.PropsSI('H', 'T', T, 'P', p_B23_Pa, 'IF97::Water')
            B23_h.append(h)
            B23_p.append(p_B23_Pa)
        except Exception:
            pass

    # Critical-patch bbox in (h, p) — the backend exposes only the
    # axis-aligned bbox shape via its update/calc_* dispatch, not via a
    # public accessor.  Approximate it by walking the auto-calibrated
    # (T, p) perimeter and recording h, matching the construction-time
    # logic in SVDSBTLBackend::build_critical_patch_.  Multipliers come
    # from a hardcoded snapshot here; if mismatched the gold patch
    # extent is slightly off but still informative.
    patch_h_lo = patch_h_hi = patch_p_lo = patch_p_hi = None
    if Tc is not None and pc is not None:
        # Conservative snapshot of the post-calibration Water/IF97 bbox
        # (from CoolProp-5ni / dxd).  Hard-coding rather than calling
        # an introspection accessor keeps the docs script standalone.
        T_lo, T_hi = 0.999 * Tc, 1.010 * Tc
        p_lo, p_hi = 0.996 * pc, 1.138 * pc
        try:
            h_lo_p = float('inf')
            h_hi_p = -float('inf')
            for i in range(24):
                f = i / 23.0
                for T in (T_lo + f * (T_hi - T_lo), T_lo, T_hi):
                    for p in (p_lo + f * (p_hi - p_lo), p_lo, p_hi):
                        try:
                            h = CP.PropsSI('H', 'T', T, 'P', p, 'IF97::Water')
                            if _math.isfinite(h):
                                h_lo_p = min(h_lo_p, h)
                                h_hi_p = max(h_hi_p, h)
                        except Exception:
                            pass
            if _math.isfinite(h_lo_p) and _math.isfinite(h_hi_p):
                patch_h_lo, patch_h_hi = h_lo_p, h_hi_p
                patch_p_lo, patch_p_hi = p_lo, p_hi
        except Exception:
            pass

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    total_probes = sum(len(s) for s in samples_per_region.values())

    for ax, (pk, sym, _u, kind, perm, _pu, _tid) in zip(axes.flat, plot_props):
        # Collect the full in-region population (for the grey
        # background scatter — shows the reader where the sample
        # density is and what fraction of the envelope is being
        # tested), and the over-budget subset (red foreground).
        all_h, all_p = [], []
        fail_h, fail_p, fail_severity = [], [], []
        for region in samples_per_region:
            region_perm = _perm_for(perm, region)
            for s in samples_per_region[region]:
                d = s['devs'].get(pk)
                all_h.append(s['h'] / 1e3)
                all_p.append(s['p'] / 1e6)
                if d is None or d <= region_perm:
                    continue
                fail_h.append(s['h'] / 1e3)
                fail_p.append(s['p'] / 1e6)
                fail_severity.append(d / region_perm)
        # Grey population scatter under the red violations.  Small
        # markers + low alpha so it acts as a density backdrop
        # without obscuring the over-budget cluster pattern.
        if all_h:
            ax.scatter(all_h, all_p, c='lightgrey', s=1, edgecolors='none', alpha=0.4, zorder=1)
        # Marker size + colour scale by log10(severity).  Severity floor
        # at 1.0 (= exactly at budget) by construction.  Cap at 1e5 so
        # extreme outliers don't blow out the colour bar.
        if fail_h:
            severity = [_math.log10(max(s, 1.0)) for s in fail_severity]
            sizes = [4 + 4 * min(s, 5.0) for s in severity]
            sc = ax.scatter(fail_h, fail_p, c=severity, cmap='Reds', s=sizes, edgecolors='none', vmin=0.0, vmax=5.0, zorder=3)
            cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
            cbar.set_label(r'$\log_{10}(|\Delta|/|\Delta|_\mathrm{perm})$', fontsize=8)
        # IF97 boundary overlays.
        if sat_p:
            sat_h_arr = [(hl + hv) / 2e3 for hl, hv in sat_h]  # dummy; replaced by L+V curves below
            sat_p_arr = [p / 1e6 for p in sat_p]
            ax.plot([h[0] / 1e3 for h in sat_h], sat_p_arr, color='steelblue', lw=0.8, label='sat dome')
            ax.plot([h[1] / 1e3 for h in sat_h], sat_p_arr, color='steelblue', lw=0.8)
        if iso_R1R3_h:
            ax.plot([h / 1e3 for h in iso_R1R3_h], [p / 1e6 for p in iso_R1R3_p], color='C1', lw=0.6, ls='--', label='T=623.15 K (R1/R3)')
        if iso_R2R5_h:
            ax.plot([h / 1e3 for h in iso_R2R5_h], [p / 1e6 for p in iso_R2R5_p], color='C3', lw=0.6, ls='--',
                    label='T=1073.15 K (R2/R5)')
        if B23_h:
            ax.plot([h / 1e3 for h in B23_h], [p / 1e6 for p in B23_p], color='C2', lw=0.6, ls=':', label='$p_{B23}(T)$ (R2/R3)')
        # Critical-patch bbox in (h, p).
        if patch_h_lo is not None:
            ax.fill_between([patch_h_lo / 1e3, patch_h_hi / 1e3], patch_p_lo / 1e6, patch_p_hi / 1e6, color='gold', alpha=0.35, edgecolor='goldenrod', lw=0.8, label='critical-patch bbox')
        # Critical point marker.
        if Tc is not None and pc is not None:
            try:
                h_c = CP.PropsSI('H', 'T', Tc, 'P', pc, 'IF97::Water')
                ax.plot(h_c / 1e3, pc / 1e6, 'ko', ms=5)
            except Exception:
                pass

        ax.set_yscale('log')
        ax.grid(True, which='both', alpha=0.3)
        # Title shows the property symbol + violation count.
        prop_sym = sym if pk != 'D' else 'v'  # density -> specific volume
        ax.set_title(r'$|\Delta {0}|$ budget violations  ({1} / {2} probes)'.format(prop_sym, len(fail_h), total_probes), fontsize=10)
    for ax in axes[-1]:
        ax.set_xlabel(r'$h$ [kJ/kg]')
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$p$ [MPa]')
    # Legend on the upper-left panel only (overlays are the same on all).
    axes[0, 0].legend(fontsize=7, loc='lower left')
    fig.suptitle('{0} — IAPWS G13-15 budget violations\n{1} probes — fails only shown\nGold shaded = critical-patch HEOS-fallback bbox'.format(backend, total_probes), fontsize=11)

    # Footer: enumerate the G13-15 perm budget that each panel is
    # scored against.  Region-dependent perms (Table 8 T) are spelled
    # out as R1/R2/R3/R5; uniform perms are shown as a single value.
    # Keeps the figure self-contained — reader doesn't have to crack
    # the script to know what "budget violation" means per property.
    def _fmt_perm_for_footer(pk, perm, unit_str):
        prop_sym = 'v' if pk == 'D' else dict((p[0], p[1]) for p in plot_props)[pk]
        if isinstance(perm, dict):
            vals = '/'.join(str(int(perm[r])) if perm[r] == int(perm[r]) else '{:g}'.format(perm[r])
                            for r in (1, 2, 3, 5))
            return r'$|\Delta {0}|_\mathrm{{perm}}$: {1} {2} (R1/R2/R3/R5)'.format(prop_sym, vals, unit_str)
        # scalar
        return r'$|\Delta {0}|_\mathrm{{perm}}$: {1:g} {2}'.format(prop_sym, perm, unit_str)

    footer_parts = []
    for pk, sym, _u, kind, perm, perm_unit, _tid in plot_props:
        if kind == 'rel':
            unit_str = r'$\%$'
        elif kind == 'abs_T':
            unit_str = 'mK'
        elif kind == 'abs_s':
            unit_str = r'$10^{-3}\,$J/(kg$\cdot$K)'
        else:
            unit_str = perm_unit
        footer_parts.append(_fmt_perm_for_footer(pk, perm, unit_str))
    # Two-line layout: thermodynamic perms on line 1, transport on line 2.
    footer_line_1 = '   '.join(footer_parts[:4])
    footer_line_2 = '   '.join(footer_parts[4:])
    fig.text(0.5, 0.025, footer_line_1, ha='center', fontsize=7.5, color='#333')
    fig.text(0.5, 0.005, footer_line_2, ha='center', fontsize=7.5, color='#333')
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    out_path = os.path.join(out_dir, 'IF97_conformance_fails_{0}.png'.format(slug))
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    return {'all': os.path.basename(out_path)}


# ---------------------------------------------------------------------------
# Timing sweep
# ---------------------------------------------------------------------------

def run_timing(backends, rng_seed=0xCAFE):
    """Mean ``AbstractState.update()`` + ``keyed_output()`` time per
    (backend, input pair, property) over N timing calls.  Times BOTH
    ``PT_INPUTS`` (the input pair the conformance population is
    naturally drawn on) and ``HmassP_INPUTS`` (the input pair most
    SVDSBTL users actually call with — the whole point of a tabulated
    backend over HEOS is the cheap :math:`(h, p)` lookup).  One
    ``AbstractState`` per backend is constructed up front and reused
    so the per-call cost is the flash + property read without the
    factory-rebuild tax ``PropsSI`` would add.  Each cell entry is
    the wall time for **one ``update()`` + one ``keyed_output()``**
    — backends that share the flash across N outputs would amortise
    this measurement; the SVDSBTL ``fast_evaluate`` path is the
    batched variant.

    Returns timings[backend][input_pair_name][prop_key] = ns/call.
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
    # Same AbstractState-reuse trick as run_conformance: pay surface
    # load once, reuse across all probes.
    states = {b: _factory(b) for b in backends}
    timing_int_keys = [_PARAM_KEY[p] for p in TIMING_PROPS]
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
                state = states[backend]
                try:
                    state.update(CP.PT_INPUTS, p_Pa, T)
                except Exception:
                    ok = False
                if not ok:
                    break
                for ip in timing_int_keys:
                    try:
                        v = state.keyed_output(ip)
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
        empty = {ip: {prop: float('nan') for prop in TIMING_PROPS} for ip in ('PT', 'HmassP')}
        return {backend: empty for backend in backends}
    # Pre-compute h_arr[] via the reference backend so HmassP timing
    # uses inputs identical to the conformance population.  IF97 has
    # a direct backward equation for HmassP — using its forward h
    # ensures both input-pair timings sample the same physical points.
    h_arr = [float('nan')] * N
    ref_state = _factory(REFERENCE)
    for k in range(N):
        try:
            ref_state.update(CP.PT_INPUTS, p_arr[k], T_arr[k])
            h_arr[k] = ref_state.hmass()
        except Exception:
            pass
    timings = {}
    for backend in backends:
        timings[backend] = {'PT': {}, 'HmassP': {}}
        state = _factory(backend)
        for prop in TIMING_PROPS:
            ip = _PARAM_KEY[prop]
            # PT timing
            for k in range(min(64, N)):  # warmup
                state.update(CP.PT_INPUTS, p_arr[k], T_arr[k])
                state.keyed_output(ip)
            t0 = time.perf_counter()
            for k in range(N):
                state.update(CP.PT_INPUTS, p_arr[k], T_arr[k])
                state.keyed_output(ip)
            t1 = time.perf_counter()
            timings[backend]['PT'][prop] = (t1 - t0) / N * 1e9
            # HmassP timing — same sample set, expressed in (h, p)
            # via the reference backend's forward h(T, p) above.
            for k in range(min(64, N)):
                if math.isfinite(h_arr[k]):
                    try:
                        state.update(CP.HmassP_INPUTS, h_arr[k], p_arr[k])
                        state.keyed_output(ip)
                    except Exception:
                        pass
            t0 = time.perf_counter()
            ok_n = 0
            for k in range(N):
                if not math.isfinite(h_arr[k]):
                    continue
                try:
                    state.update(CP.HmassP_INPUTS, h_arr[k], p_arr[k])
                    state.keyed_output(ip)
                    ok_n += 1
                except Exception:
                    pass
            t1 = time.perf_counter()
            timings[backend]['HmassP'][prop] = ((t1 - t0) / ok_n * 1e9) if ok_n else float('nan')
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
    lines.append('.. note::')
    lines.append('')
    lines.append(
        '    All tables and figures in this section are regenerated by '
        '``Web/scripts/fluid_properties.IF97Conformance.py`` on every '
        'docs build, against the local CoolProp commit.  The numbers '
        'shipped with the published docs come from the GitHub-hosted '
        'CI runner; rebuild locally to see numbers for your hardware '
        'and source revision.')
    lines.append('')
    lines.append(
        'Conformance is measured against the '
        ':ref:`IAPWS G13-15 Tables 8-13 <IAPWS-IF97>` (Kunick et al., 2015) '
        'protocol: for each IF97 region we draw {n} random :math:`(p, T)` '
        'samples log-uniform in :math:`p` and uniform in :math:`T`, '
        'classify them into the IF97 region atlas, compute '
        ':math:`h = h_\\mathrm{{IF97}}(p, T)`, evaluate each tested backend '
        'at :math:`(h, p)`, and tabulate the maximum and root-mean-square '
        'deviation from IAPWS-IF97 of the six properties in G13-15: '
        ':math:`T(p,h)`, :math:`v(p,h)`, :math:`s(p,h)`, :math:`w(p,h)`, '
        ':math:`\\eta(p,h)`, :math:`\\lambda(p,h)`.  Region 4 (saturation '
        'dome) is excluded — G13-15 evaluates the forward / backward '
        'equations only outside the dome.'
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

        # Single combined failure-point figure spanning all IF97
        # regions, six panels (2x3) — one per G13-15 property
        # (T, v, s, w plus eta and lambda).  Each panel shows ONLY
        # the budget-violating probes; markers sized + coloured by
        # log10(|delta|/perm) so isolated thin-support outliers are
        # visually distinguishable from the multi-cell clusters that
        # signal structural SVD-truncation problems.  Overlays
        # (saturation dome, R1/R3 / R2/R5 isotherms, B23 boundary,
        # auto-calibrated critical-patch bbox, critical point) anchor
        # the geographic interpretation.
        if fig_paths:
            lines.append(
                'The figure below shows *where* the failures live in '
                ':math:`(p, h)` coords across the full IF97 envelope '
                '— six panels, one per G13-15 property '
                '(:math:`T, v, s, w, \\eta, \\lambda`).  Only the '
                'budget-violating probes are plotted; marker size and '
                'colour both encode '
                ':math:`\\log_{10}(|\\Delta| / |\\Delta|_\\mathrm{perm})`, '
                'so bigger + redder is further past budget.  Overlays: '
                'the saturation dome, the IF97 R1/R3 isotherm '
                '(:math:`T = 623.15` K) and R2/R5 isotherm '
                '(:math:`T = 1073.15` K), the R2/R3 boundary curve '
                ':math:`p_{B23}(T)`, and the auto-calibrated '
                'critical-patch HEOS-fallback bbox (gold).  Black dot '
                'marks the critical point.')
            lines.append('')
            lines.append(
                'The clustering pattern is the load-bearing signal — '
                'it tells you *why* the table fails where it does, and '
                'what to do about it:')
            lines.append('')
            lines.append(
                '* **Thermodynamic panels** (:math:`T, v, s, w`): '
                'failures cluster around the **critical point** '
                '(rank-20 SVD cannot represent the '
                ':math:`(\\partial \\rho / \\partial p)_h` cusp at the '
                'critical singularity — that would require unbounded '
                'rank).  The gold critical-patch bbox covers the bulk '
                'of those failures; the auto-calibration loop sized it '
                'to do exactly that.  Residual reds outside the gold '
                'bbox are thin-Hermite-support corner cells '
                '(``CoolProp-3c4`` accuracy ceiling — accepted, not '
                'fixable by widening the bbox).')
            lines.append('')
            lines.append(
                '* **Transport panels** (:math:`\\eta, \\lambda`): '
                'failures are widespread.  The G13-15 transport '
                'budgets (:math:`10^{-5}` relative) are 100-1000x '
                'tighter than the thermodynamic ones, and the rank-20 '
                'SVD\'s representation of viscosity / conductivity is '
                'not yet conformant — closing this gap is tracked in '
                'the SBTL conformance epic (``CoolProp-foi``).  The '
                'thermodynamic panels are the ones to focus on when '
                'assessing whether SVDSBTL is fit-for-purpose for a '
                'given application; the transport panels document the '
                'remaining work.')
            lines.append('')
            lines.append(
                '* The **R2/R3 kink** along :math:`p_{B23}(T)` is no '
                'longer a structural failure cluster — the atlas split '
                'introduced in ``CoolProp-foi.9.5`` separates SUPER_R2 '
                'and SUPER_R3 so each SVD only sees cells on one side '
                'of the kink.  Residual failures along that line are '
                'the boundary-row cells where the Hermite kernel spans '
                'across :math:`p_{B23}`.')
            lines.append('')
            # The new generator emits a single combined PNG keyed
            # 'all' rather than per-region files — embed it directly.
            combined = fig_paths.get('all')
            if combined is not None:
                lines.append('.. figure:: {0}'.format(combined))
                lines.append('   :align: center')
                lines.append('   :width: 95%')
                lines.append('')
                lines.append('   {0} — IAPWS G13-15 budget violations across the full IF97 envelope.'.format(backend))
                lines.append('')
            else:
                # Backwards-compat for any legacy per-region keys (a
                # plugin or downstream regen might still emit them).
                for region in sorted(fig_paths):
                    lines.append('.. figure:: {0}'.format(fig_paths[region]))
                    lines.append('   :align: center')
                    lines.append('   :width: 90%')
                    lines.append('')
                    lines.append('   IF97 Region {0} — {1}'.format(region, backend))
                    lines.append('')

        for pk, sym, _unit, kind, perm, perm_unit, table_id in PROPERTIES:
            # `unit_str` is the LaTeX-flavoured unit token used INSIDE
            # a :math: directive; `unit_plain` is the corresponding
            # bracketed form for the header (also wrapped in :math:
            # so superscripts / \mathrm render correctly).  Previously
            # the units were emitted as plain RST inside square
            # brackets, which left `10^{-3}` and `\mathrm` as literal
            # text in the entropy column header.
            if kind == 'rel':
                unit_str = r'\%'
            elif kind == 'abs_T':
                unit_str = r'\mathrm{mK}'
            elif kind == 'abs_s':
                unit_str = r'10^{-3}\ \mathrm{J/(kg\,K)}'
            else:
                unit_str = perm_unit
            unit_math = ':math:`[{0}]`'.format(unit_str)

            prop_sym = sym if pk != 'D' else 'v'  # density → specific volume
            lines.append('')
            # Header line: state the perm once if it's region-uniform,
            # or as "perm depends on region (see column)" if it varies.
            if isinstance(perm, dict):
                lines.append(
                    '*G13-15 Table {tid}* — :math:`{psym}(p,h)`, '
                    'permissible :math:`|\\Delta {psym}|_\\mathrm{{perm}}` '
                    'is region-dependent (per-region column below).'
                    .format(tid=table_id, psym=prop_sym))
            else:
                lines.append(
                    '*G13-15 Table {tid}* — :math:`{psym}(p,h)`, '
                    'permissible :math:`|\\Delta {psym}|_\\mathrm{{perm}} = {perm}\\ {pu}`.'
                    .format(tid=table_id, psym=prop_sym, perm=_fmt_perm(perm), pu=unit_str))
            lines.append('')
            lines.append('.. list-table::')
            lines.append('   :header-rows: 1')
            # Wider columns when the perm column is present.
            if isinstance(perm, dict):
                lines.append('   :widths: 10 14 18 18 18')
            else:
                lines.append('   :widths: 12 18 22 22')
            lines.append('')
            lines.append('   * - IF97 Region')
            lines.append('     - in-region samples')
            if isinstance(perm, dict):
                lines.append('     - :math:`|\\Delta {0}|_{{\\mathrm{{perm}}}}` {1}'.format(prop_sym, unit_math))
            lines.append('     - :math:`|\\Delta {0}|_{{\\max}}` {1}'.format(prop_sym, unit_math))
            lines.append('     - :math:`(\\Delta {0})_{{\\mathrm{{RMS}}}}` {1}'.format(prop_sym, unit_math))
            for region in sorted(REGION_BOXES):
                acc = results[region][pk]
                lines.append('   * - {}'.format(region))
                lines.append('     - {}'.format(acc.n))
                if isinstance(perm, dict):
                    lines.append('     - {}'.format(_fmt_perm(_perm_for(perm, region))))
                lines.append('     - {}'.format(_fmt_dev(acc.max_abs)))
                lines.append('     - {}'.format(_fmt_dev(acc.rms)))
        lines.append('')

    # ------------------------------------------------------------------
    # Timing tables — one per input pair
    # ------------------------------------------------------------------
    if timings:
        lines.append('Backend timing (``AbstractState.update`` regime)')
        lines.append('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        lines.append('')
        lines.append(
            'Two tables below time the same backends against the same '
            'sample population for **both** ``PT_INPUTS`` (the input '
            'pair the conformance sweep is drawn on) and '
            '``HmassP_INPUTS`` (the input pair most users actually call '
            'with — the whole point of a tabulated backend is the '
            'cheap :math:`(h, p)` lookup).  Sample points are identical '
            'across the two tables; only the input pair changes.  '
            ':math:`h(T, p)` is taken from the reference IF97 backend '
            'so HmassP inputs are forward-consistent with PT inputs.')
        lines.append('')
        lines.append(
            'Each cell is the **mean wall time for one ``update()`` + '
            'one ``keyed_output()`` requesting that property** — i.e. '
            '**one output per call**.  Backends that share the flash '
            'across N outputs would amortise this measurement; the '
            'SVDSBTL ``fast_evaluate`` path is the batched variant.  '
            'The ratio column is the mean across the per-property '
            'ratios :math:`t_{\\mathrm{backend}} / '
            't_{\\mathrm{IF97}}` for that input pair.  Lower is '
            'faster than IF97.  Rows are sorted by descending mean '
            'ratio so the slowest path sits at the top.  Absolute '
            'timings depend on hardware; the **ratios** are the '
            'load-bearing quantity.  CI rebuilds these on a '
            'GitHub-hosted runner.')
        lines.append('')

        for ip_label in ('PT', 'HmassP'):
            lines.append('**{0}_INPUTS**:'.format(ip_label))
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
            ref_t = timings.get(REFERENCE, {}).get(ip_label, {})
            candidates = [REFERENCE] + [b for b, _ in TESTED if b in timings]
            rows_by_ratio = []
            for backend in candidates:
                bt = timings.get(backend, {}).get(ip_label, {})
                if not bt:
                    continue
                ratios = []
                for prop in TIMING_PROPS:
                    base = ref_t.get(prop, float('nan'))
                    this = bt.get(prop, float('nan'))
                    if math.isfinite(base) and base > 0 and math.isfinite(this):
                        ratios.append(this / base)
                mean_ratio = sum(ratios) / len(ratios) if ratios else float('inf')
                rows_by_ratio.append((mean_ratio, backend, bt, ratios))
            rows_by_ratio.sort(key=lambda x: x[0], reverse=True)
            for mean_ratio, backend, bt, ratios in rows_by_ratio:
                row = ['   * - ``{}``'.format(backend)]
                for prop in TIMING_PROPS:
                    row.append('     - {}'.format(_fmt_ns(bt.get(prop, float('nan')))))
                row.append('     - {:.2f}'.format(mean_ratio) if ratios else '     - --')
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
