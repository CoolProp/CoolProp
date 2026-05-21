#!/usr/bin/env python
"""Profile IF97 (and SVDSBTL, when available) per-call timing.

Walks a sparse (T, p) probe grid via IF97 to get h(T, p), then queries
each backend at ``HmassP_INPUTS(h, p)`` and times the per-call wall cost
via :func:`time.perf_counter_ns`.  Backends not in the build are
skipped without breaking the rest.

What this measures
------------------

* **per-call** — one ``AbstractState.update()`` call per probe, repeated
  ``REPEATS`` times; median ns/call reported.  Includes the Python-side
  marshalling cost (~1 us/call floor).  This is what your code sees if
  you ``update()`` once per state-point in a Python loop.
* **batch** — :meth:`AbstractState.fast_evaluate`, the cache-bypassing
  zero-allocation batch API.  Times the whole probe set as one C++ call,
  divides by point count.  This is the C++-native per-point cost — what
  you'd see in tight inner loops once Python overhead is amortized.

Run it
------

::

    python Web/scripts/fluid_properties.IF97Profile.py

Output: ``Web/fluid_properties/IF97_profile.png`` + stdout summary.
The figure pages itself to however many backends were actually
available; missing ones drop out cleanly.

Caveats
-------

* The first call against each backend warms the in-memory surface.
  Reported timings are the *median* of warm-cache calls.
* Per-call ``update()`` times are ~1 us higher than C++-native because
  of the wrapper marshalling.  ``fast_evaluate`` does not pay this.
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

# Sparse-by-design.  G13-15 conformance probes 8000+ points per region
# but profiling is about per-call cost, not statistical coverage; 30x30
# is enough to estimate medians and show spatial trends without making
# the docs build slow.
NT = 30
NP = 30
T_RANGE = (290.0, 1000.0)   # K — covers R1/R2/R5 for Water under IF97
P_RANGE = (1.0e4, 5.0e7)    # Pa — 0.01 to 50 MPa
REPEATS = 100               # per-point repeats; median absorbs jitter

# Backends to try.  Missing ones drop cleanly; the figure pages to the
# number that worked.
BACKENDS = [
    ('IF97',         'Water', 'IF97'),
    ('HEOS',         'Water', 'HEOS'),
    ('SVDSBTL&HEOS', 'Water', 'SVDSBTL+HEOS'),
    ('SVDSBTL&IF97', 'Water', 'SVDSBTL+IF97'),
]


def make_probes() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    Ts = np.linspace(*T_RANGE, NT)
    ps = np.geomspace(*P_RANGE, NP)
    TT, PP = np.meshgrid(Ts, ps, indexing='ij')
    HH = np.full_like(TT, np.nan)
    if97 = AbstractState('IF97', 'Water')
    for i in range(NT):
        for j in range(NP):
            try:
                if97.update(CP.PT_INPUTS, PP[i, j], TT[i, j])
                if np.isfinite(if97.hmass()):
                    HH[i, j] = if97.hmass()
            except Exception:
                pass
    return TT, PP, HH


def try_factory(backend: str, fluid: str):
    try:
        return AbstractState(backend, fluid)
    except Exception as exc:
        print(f'  {backend}: skipped ({exc})')
        return None


def time_per_call(s, h: np.ndarray, p: np.ndarray) -> np.ndarray:
    ns = np.full_like(h, np.nan)
    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            if not np.isfinite(h[i, j]):
                continue
            try:
                s.update(CP.HmassP_INPUTS, h[i, j], p[i, j])  # warmup
            except Exception:
                continue
            samples = []
            for _ in range(REPEATS):
                t0 = time.perf_counter_ns()
                s.update(CP.HmassP_INPUTS, h[i, j], p[i, j])
                samples.append(time.perf_counter_ns() - t0)
            ns[i, j] = statistics.median(samples)
    return ns


def time_fast_evaluate(s, h: np.ndarray, p: np.ndarray) -> float:
    h_flat = h.ravel()
    p_flat = p.ravel()
    valid = np.isfinite(h_flat)
    h_flat = np.ascontiguousarray(h_flat[valid])
    p_flat = np.ascontiguousarray(p_flat[valid])
    if h_flat.size == 0:
        return float('nan')
    outputs = np.array([CP.iT, CP.iDmass], dtype=np.int32)
    out = np.zeros((h_flat.size, 2), dtype=np.float64)
    status = np.zeros(h_flat.size, dtype=np.int32)
    try:
        s.fast_evaluate(CP.HmassP_INPUTS, h_flat, p_flat, outputs, out, status)
    except Exception:
        return float('nan')
    samples = []
    for _ in range(20):
        t0 = time.perf_counter_ns()
        s.fast_evaluate(CP.HmassP_INPUTS, h_flat, p_flat, outputs, out, status)
        samples.append((time.perf_counter_ns() - t0) / h_flat.size)
    return statistics.median(samples)


def main() -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    print(f'Building {NT}x{NP} probe grid via IF97...')
    TT, PP, HH = make_probes()
    n_valid = int(np.isfinite(HH).sum())
    print(f'  {n_valid} valid single-phase probes\n')

    states = []
    print('per-call update(HmassP_INPUTS, h, p):')
    for backend, fluid, label in BACKENDS:
        s = try_factory(backend, fluid)
        if s is None:
            continue
        ns = time_per_call(s, HH, PP)
        finite = ns[np.isfinite(ns)]
        if finite.size == 0:
            print(f'  {label}: no valid points')
            continue
        med = float(np.median(finite))
        print(f'  {label:14s} median {med:7.0f} ns/call  ({finite.size} pts)')
        states.append((label, ns, s, med))

    print('\nfast_evaluate (batch, high-level API):')
    fast_medians = {}
    for label, _, s, _ in states:
        per_pt = time_fast_evaluate(s, HH, PP)
        if np.isfinite(per_pt):
            print(f'  {label:14s} {per_pt:7.0f} ns/point (batch)')
            fast_medians[label] = per_pt
        else:
            print(f'  {label:14s} fast_evaluate not supported')

    if not states:
        print('\nNo backends available; nothing to plot.')
        return

    # Lay out N panels in a 2-up grid; trailing panels stay blank if odd.
    n = len(states)
    cols = 2 if n > 1 else 1
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.4 * cols, 4.0 * rows),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.flatten()

    pos = np.concatenate([ns[np.isfinite(ns)] for _, ns, _, _ in states if np.isfinite(ns).any()])
    vmin = float(np.percentile(pos, 5))
    vmax = float(np.percentile(pos, 95))
    norm = LogNorm(vmin=max(vmin, 1.0), vmax=max(vmax, vmin * 10))

    sc = None
    for ax, (label, ns, _, med) in zip(axes, states):
        h_kJ = HH / 1e3
        sc = ax.scatter(h_kJ.ravel(), PP.ravel() / 1e6,
                        c=ns.ravel(), s=22, cmap='viridis', norm=norm,
                        edgecolor='none')
        ax.set_yscale('log')
        title = f'{label}\nmedian {med:.0f} ns/call'
        if label in fast_medians:
            title += f'  |  fast_evaluate {fast_medians[label]:.0f} ns/pt'
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('h [kJ/kg]')
        ax.set_ylabel('p [MPa]')
        ax.grid(True, alpha=0.2)
    for ax in axes[n:]:
        ax.set_visible(False)
    if sc is not None:
        cbar = fig.colorbar(sc, ax=axes[:n].tolist(), fraction=0.025, pad=0.04)
        cbar.set_label('per-call wall time [ns]')

    fig.suptitle(f'CoolProp HmassP_INPUTS timing — Water  ({n_valid} probes, '
                 f'{REPEATS} repeats / point)', fontsize=12)
    fig.savefig(OUT_PNG, dpi=120, bbox_inches='tight')
    print(f'\nWrote {OUT_PNG}')


if __name__ == '__main__':
    main()
