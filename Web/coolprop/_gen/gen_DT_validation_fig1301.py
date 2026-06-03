"""Reproduce yufang67's (T, ρ) error scatter from issue #1301 for the
DT-indexed SVDSBTL preset (CoolProp-i7j), fully in Python — does its
own (T, ρ) sweep via CoolProp (no C++/CSV dependency) and renders the
two-panel figure in the ORIGINAL report's format.

Requires a CoolProp Python wrapper built from THIS branch (so the
SVDSBTL backend understands DmassT_INPUTS).  Build it with, e.g.:

    cmake --build build_catch --target CoolProp   # or the Python target
    pip install -e wrappers/Python                 # editable wrapper

Figure format (matches #1301):
  * x = T [K] over [220, 500], y = ρ [kg/m³] over [0, 1200] (linear)
  * color = |P_backend − P_EOS| / P_EOS · 100  [%], jet, log-normed
  * left  = BICUBIC&HEOS  (reproduces the near-saturation error band
            + negative-pressure cells)
  * right = SVDSBTL&HEOS DT-indexed (the fix)
  * P < 0 cells overplotted as black ×

Usage:
    python3 Web/coolprop/_gen/gen_DT_validation_fig1301.py
    python3 ... --n 200000        # denser scatter (default 20000)
    python3 ... --grid 400,1600,30  # bump SVDSBTL NT,NR,rank via options
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# yufang67's bounds from #1301.
T_LO, T_HI = 220.0, 500.0
RHO_LO, RHO_HI = 1.0e-2, 1200.0
OUT_DEFAULT = Path(__file__).resolve().parents[1] / "fig1301_dt_validation.png"


def _sweep(n_samples: int, grid: str | None):
    import CoolProp.CoolProp as CP

    svd_fluid = "CarbonDioxide"
    if grid:
        nt, nr, rank = (int(x) for x in grid.split(","))
        svd_fluid = f'CarbonDioxide?{{"grid":{{"NT":{nt},"NR":{nr},"rank":{rank}}}}}'

    heos = CP.AbstractState("HEOS", "CarbonDioxide")
    svd = CP.AbstractState("SVDSBTL&HEOS", svd_fluid)
    bic = CP.AbstractState("BICUBIC&HEOS", "CarbonDioxide")

    rng = np.random.default_rng(1301)
    T = rng.uniform(T_LO, T_HI, n_samples)
    rho = rng.uniform(RHO_LO, RHO_HI, n_samples)

    # SVDSBTL: vectorized via fast_evaluate (DmassT) — bypasses the
    # per-call cache overhead that dominates a 200k point-by-point
    # sweep.  val1 = D, val2 = T for DmassT_INPUTS.  out[:,0] = p.
    out = np.full((n_samples, 1), np.nan)
    status = np.zeros(n_samples, dtype=np.int32)
    svd.fast_evaluate(
        CP.DmassT_INPUTS,
        np.ascontiguousarray(rho), np.ascontiguousarray(T),
        np.array([CP.iP], dtype=np.int32),
        out, status,
    )
    p_svd = np.where(status == 0, out[:, 0], np.nan)

    # HEOS (truth) and BICUBIC have no DmassT fast path (HEOS has none;
    # tabular fast_evaluate only does PT / HmolarP), so they go
    # point-by-point.  Both are cheap per call relative to SVDSBTL's
    # cache check, so this is no longer the bottleneck.
    p_heos = np.full(n_samples, np.nan)
    p_bic = np.full(n_samples, np.nan)
    for i in range(n_samples):
        for state, arr in ((heos, p_heos), (bic, p_bic)):
            try:
                state.update(CP.DmassT_INPUTS, rho[i], T[i])
                arr[i] = state.p()
            except Exception:
                pass  # leave NaN — outside this backend's validity envelope
    return T, rho, p_heos, p_svd, p_bic


def _panel(ax, T, rho, p_heos, p_bk, title):
    valid = np.isfinite(p_heos) & np.isfinite(p_bk) & (p_heos != 0.0)
    pct = np.full_like(p_heos, np.nan)
    pct[valid] = np.abs(p_bk[valid] - p_heos[valid]) / np.abs(p_heos[valid]) * 100.0
    vals = np.maximum(pct, 1e-10)
    finite = np.isfinite(vals)
    sc = ax.scatter(
        T[finite], rho[finite], c=vals[finite], s=6, edgecolor="none",
        cmap="jet", norm=LogNorm(vmin=1e-8, vmax=1e2),
    )
    neg = np.isfinite(p_heos) & (p_heos > 0) & np.isfinite(p_bk) & (p_bk <= 0)
    if neg.any():
        ax.scatter(T[neg], rho[neg], s=22, marker="x", color="black", linewidths=1.1, label="P < 0 (sign flip)")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.set_xlim(T_LO, T_HI)
    ax.set_ylim(0.0, RHO_HI)
    ax.set_xlabel("Temperature  [K]")
    ax.set_ylabel("density  [kg/m³]")
    ax.set_title(f"{title}  ({int(neg.sum())} cells P<0)")
    plt.colorbar(sc, ax=ax, label="Error in pressure  [%]")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20000, help="number of random (T, ρ) draws")
    ap.add_argument("--grid", type=str, default=None, help="SVDSBTL grid as NT,NR,rank (e.g. 400,1600,30)")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args(argv)

    try:
        T, rho, p_heos, p_svd, p_bic = _sweep(args.n, args.grid)
    except ImportError:
        sys.stderr.write(
            "CoolProp Python wrapper not importable — build it from this branch first "
            "(SVDSBTL DmassT support is required).\n"
        )
        return 1

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
    fig.suptitle(
        f"CO₂ P(ρ, T) error vs HEOS — issue #1301 reproduction (CoolProp-i7j, N={args.n} random draws)",
        fontsize=13,
    )
    _panel(axes[0], T, rho, p_heos, p_bic, "BICUBIC&HEOS")
    _panel(axes[1], T, rho, p_heos, p_svd, "SVDSBTL&HEOS (DT-indexed)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140, bbox_inches="tight")
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
