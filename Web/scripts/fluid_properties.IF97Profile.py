#!/usr/bin/env python
"""Profile HmassP_INPUTS backends via batched timing.

Builds a sparse (T, p) probe grid via IF97 to get h_truth(T, p), then
times each available backend by *batches* — either via
:meth:`AbstractState.fast_evaluate` (the cache-bypassing zero-allocation
batch API) when the backend implements it, or via a manual same-probe
loop inside a single :func:`time.perf_counter_ns` block as the fallback.

Why batch
---------

Per-call timing in a Python loop bottoms out at the
``perf_counter_ns`` resolution and is dominated by the wrapper
marshalling cost (~1 us / call), so the fastest backends all look the
same. Timing whole batches at a time amortizes both effects: the
wrapper overhead becomes a one-time tax over ``N`` points instead of
``N`` taxes of one, and a many-tick batch is statistically robust.

``fast_evaluate`` exposes the zero-allocation C++ batch path directly;
for backends that don't implement it yet we time
``[update(...) for each point]`` inside one tick block as the closest
equivalent. Both numbers are reported as ``ns/point``.

Run it
------

::

    python Web/scripts/fluid_properties.IF97Profile.py

Output: ``Web/fluid_properties/IF97_profile.png`` + a stdout summary
listing each backend's per-point ns from N_RUNS independent batches.
"""
from __future__ import annotations

import os
import statistics
import time

import numpy as np
import CoolProp.CoolProp as CP
from CoolProp import AbstractState

WEB_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_PNG = os.path.join(WEB_DIR, 'fluid_properties', 'IF97_profile.png')

# Sparse-by-design probe grid; we want enough states to amortize the
# fast_evaluate batch and give the manual loop a few hundred us to
# stabilize, not statistical coverage.
NT = 30
NP = 30
T_RANGE = (290.0, 1000.0)   # K — covers R1/R2/R5 for Water under IF97
P_RANGE = (1.0e4, 5.0e7)    # Pa — 0.01 to 50 MPa

# Number of independent batch runs per backend.  Each run times the full
# valid probe set; we report the distribution (median, P10, P90) across
# runs.  Cheap to bump if the noise looks high.
N_RUNS = 200

# Order left-to-right: HEOS (the slow iterative baseline) first, then
# IF97 (closed-form backward equations), then the tabular backends — so
# the figure reads as a descending cost ladder.
BACKENDS = [
    ('HEOS',         'Water', 'HEOS'),
    ('IF97',         'Water', 'IF97'),
    ('SVDSBTL&HEOS', 'Water', 'SVDSBTL+HEOS'),
    ('SVDSBTL&IF97', 'Water', 'SVDSBTL+IF97'),
]


def make_probes() -> tuple[np.ndarray, np.ndarray]:
    """Return contiguous (h, p) arrays for every IF97-valid single-phase point.

    Probes that IF97 can't evaluate (above the temperature limit, etc.) are
    dropped silently — the count is reported by the caller via len() before
    timing.  Raises RuntimeError if the grid yielded zero valid points so
    the timing loop never divides by zero.
    """
    Ts = np.linspace(*T_RANGE, NT)
    ps = np.geomspace(*P_RANGE, NP)
    if97 = AbstractState('IF97', 'Water')
    hs, prs = [], []
    skipped = 0
    for T in Ts:
        for p in ps:
            try:
                if97.update(CP.PT_INPUTS, p, T)
            except Exception:
                # IF97 rejected this (T, p) — outside the steam-table envelope,
                # below T_triple, etc.  Drop silently; the grid spec is broad
                # by design.
                skipped += 1
                continue
            h = if97.hmass()
            if np.isfinite(h):
                hs.append(h)
                prs.append(p)
            else:
                skipped += 1
    if not hs:
        raise RuntimeError(
            f'make_probes(): grid {NT}x{NP} over T in {T_RANGE} K, p in {P_RANGE} Pa '
            f'yielded zero IF97-valid points (skipped {skipped}). Adjust ranges.')
    if skipped:
        print(f'  ({skipped} of {NT * NP} grid points dropped as IF97-invalid)')
    return np.ascontiguousarray(hs), np.ascontiguousarray(prs)


def try_factory(backend: str, fluid: str):
    try:
        return AbstractState(backend, fluid)
    except Exception as exc:
        print(f'  {backend}: skipped ({exc})')
        return None


def _backend_valid_subset(s, h: np.ndarray, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Filter (h, p) down to the points the backend can actually evaluate.

    Backends differ in coverage — IF97 takes everything in the IF97 envelope,
    HEOS takes a superset, SVDSBTL might reject points outside its surface
    bbox.  Both timing paths must measure on the same subset or the
    ns/point comparison is biased by varying NaN-fill counts.
    """
    keep = np.zeros(h.size, dtype=bool)
    for k in range(h.size):
        try:
            s.update(CP.HmassP_INPUTS, float(h[k]), float(p[k]))
            keep[k] = True
        except Exception:
            # Backend rejected this probe — leave it out of both fast_evaluate
            # and manual-loop timings so the per-point cost is comparable.
            pass
    return np.ascontiguousarray(h[keep]), np.ascontiguousarray(p[keep])


def time_fast_evaluate(s, h: np.ndarray, p: np.ndarray, n_runs: int) -> tuple[list[float], bool]:
    """Time ``fast_evaluate`` over the supplied probe set; returns ns/point per run.

    Caller is expected to have pre-filtered (h, p) via _backend_valid_subset
    so every probe contributes a real evaluation, not a NaN-fill.
    """
    if h.size == 0:
        return [], False
    outputs = np.array([CP.iT, CP.iDmass], dtype=np.int32)
    out = np.zeros((h.size, 2), dtype=np.float64)
    status = np.zeros(h.size, dtype=np.int32)
    try:
        s.fast_evaluate(CP.HmassP_INPUTS, h, p, outputs, out, status)  # warmup
    except AttributeError:
        # Older Python wrapper without the binding.
        return [], False
    except ValueError as exc:
        # Backend's default-impl NotImplementedError surfaces as ValueError
        # through the wrapper.  Caller falls back to the manual-loop path.
        if 'not implemented' in str(exc).lower():
            return [], False
        raise
    samples = []
    for _ in range(n_runs):
        t0 = time.perf_counter_ns()
        s.fast_evaluate(CP.HmassP_INPUTS, h, p, outputs, out, status)
        samples.append((time.perf_counter_ns() - t0) / h.size)
    return samples, True


def time_manual_batch(s, h: np.ndarray, p: np.ndarray, n_runs: int) -> list[float]:
    """Time ``[update(...) for each point]`` inside one tick block.

    Caller has already filtered (h, p) via _backend_valid_subset so the
    inner update() can run without try/except — one invalid probe in the
    timed loop would otherwise abort the whole batch.
    """
    if h.size == 0:
        return []
    # Warmup pass primes any backend-side caches.
    for k in range(h.size):
        s.update(CP.HmassP_INPUTS, float(h[k]), float(p[k]))
    samples = []
    for _ in range(n_runs):
        t0 = time.perf_counter_ns()
        for k in range(h.size):
            s.update(CP.HmassP_INPUTS, float(h[k]), float(p[k]))
        samples.append((time.perf_counter_ns() - t0) / h.size)
    return samples


def main() -> None:
    import matplotlib.pyplot as plt

    print(f'Building {NT}x{NP} probe grid via IF97...')
    h, p = make_probes()
    print(f'  {h.size} valid single-phase probes\n')

    results: list[tuple[str, list[float], str]] = []
    print(f'Timing {N_RUNS} batch runs per backend (ns/point):')
    print(f'  {"backend":18s} {"mechanism":18s} {"N":>5s} {"median":>8s} {"P10":>8s} {"P90":>8s}')
    for backend, fluid, label in BACKENDS:
        s = try_factory(backend, fluid)
        if s is None:
            continue
        # Filter the probe set down to states this backend can evaluate.
        # Both timing paths run on the SAME subset so the per-point cost
        # is comparable — otherwise fast_evaluate would silently absorb
        # NaN-fill points into its average while manual would either
        # crash or skip them.
        h_sub, p_sub = _backend_valid_subset(s, h, p)
        if h_sub.size == 0:
            print(f'  {label:18s} (no backend-valid probes; skipped)')
            continue
        fe_samples, fe_ok = time_fast_evaluate(s, h_sub, p_sub, N_RUNS)
        if fe_ok:
            results.append((label, fe_samples, 'fast_evaluate'))
            med = statistics.median(fe_samples)
            p10 = float(np.percentile(fe_samples, 10))
            p90 = float(np.percentile(fe_samples, 90))
            print(f'  {label:18s} {"fast_evaluate":18s} {h_sub.size:5d} {med:8.0f} {p10:8.0f} {p90:8.0f}')
        else:
            mb_samples = time_manual_batch(s, h_sub, p_sub, N_RUNS)
            results.append((label, mb_samples, 'manual update loop'))
            med = statistics.median(mb_samples)
            p10 = float(np.percentile(mb_samples, 10))
            p90 = float(np.percentile(mb_samples, 90))
            print(f'  {label:18s} {"manual update loop":18s} {h_sub.size:5d} {med:8.0f} {p10:8.0f} {p90:8.0f}')

    if not results:
        print('No backends available; nothing to plot.')
        return

    # Box-and-violin: each backend gets one column.  Log y because
    # SVDSBTL and HEOS differ by ~3 orders of magnitude.
    labels = [f'{lbl}\n({mech})' for lbl, _, mech in results]
    medians = [statistics.median(s) for _, s, _ in results]
    samples = [s for _, s, _ in results]

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    parts = ax.violinplot(samples, showmedians=True, widths=0.7)
    for body in parts['bodies']:
        body.set_alpha(0.55)
    ax.set_xticks(range(1, len(results) + 1))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_yscale('log')
    ax.set_ylabel('per-point batch wall time [ns]')
    ax.grid(True, axis='y', which='both', alpha=0.25)

    for i, med in enumerate(medians, start=1):
        ax.annotate(f'{med:.0f}', xy=(i, med), xytext=(8, 0),
                    textcoords='offset points', fontsize=9, va='center',
                    color='black', fontweight='bold')

    title = (f'HmassP_INPUTS batched timing — Water  ({h.size} probes, '
             f'{N_RUNS} batches per backend)\n'
             f'fast_evaluate where available; manual-update-loop fallback '
             f'(same probe set, single perf_counter block) otherwise')
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130, bbox_inches='tight')
    print(f'\nWrote {OUT_PNG}')


if __name__ == '__main__':
    main()
