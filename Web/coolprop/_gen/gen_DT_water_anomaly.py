"""Water density-anomaly validation for the DT-indexed SVDSBTL surface
(CoolProp-i7j).  Companion to gen_DT_validation_fig1301.py, focused on
the compressed-liquid band just above the saturation dome, bracketing
the 4 degC density maximum (T_anom ~277.13 K) and running up to the
melting / HEOS-validity boundary.

This is the region the LIQUID_LO / LIQUID_HI anomaly split exists for:
rho_sat,L(T) is non-monotone here (peaks at T_anom), so a single LIQUID
region would straddle the density-maximum hairpin.  The figure shows
the split keeps P(rho,T) clean across it.

Requires a CoolProp Python wrapper built from this branch (DmassT
SVDSBTL support).

Format mirrors the #1301 figure:
  * x = T [K], y = rho [kg/m3] (linear)
  * color = |P_backend - P_EOS| / P_EOS * 100  [%], jet, log-normed
  * left = BICUBIC&HEOS, right = SVDSBTL&HEOS (DT-indexed)
  * rho_sat,L(T) overlaid (white); T_anom marked (dashed)

Usage:
    python3 Web/coolprop/_gen/gen_DT_water_anomaly.py
    python3 ... --n 100000
    python3 ... --interactive      # live zoom/pan/inspect window
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

# Interactive mode must pick a GUI backend BEFORE importing pyplot.
# Detect the flag from argv directly (argparse runs later).  Only force
# "MacOSX" on darwin — that backend's Cocoa support isn't available on
# Linux/Windows; elsewhere let matplotlib auto-select a GUI backend.
_INTERACTIVE = "--interactive" in sys.argv
if _INTERACTIVE and sys.platform == "darwin":
    matplotlib.use("MacOSX")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

# Band brackets the anomaly (277.13 K) with compressed-liquid densities
# from just above rho_sat,L up to the melting/validity ceiling.
T_LO, T_HI = 273.2, 300.0
RHO_LO, RHO_HI = 996.0, 1090.0
OUT_DEFAULT = Path(__file__).resolve().parents[1] / "fig_dt_water_anomaly.png"


def _sweep(n_samples: int):
    import CoolProp.CoolProp as CP

    heos = CP.AbstractState("HEOS", "Water")
    svd = CP.AbstractState("SVDSBTL&HEOS", "Water")
    bic = CP.AbstractState("BICUBIC&HEOS", "Water")

    rng = np.random.default_rng(277)
    T = rng.uniform(T_LO, T_HI, n_samples)
    rho = rng.uniform(RHO_LO, RHO_HI, n_samples)

    # SVDSBTL via fast_evaluate (DmassT): val1=rho, val2=T, out=[iP].
    out = np.full((n_samples, 1), np.nan)
    status = np.zeros(n_samples, dtype=np.int32)
    svd.fast_evaluate(CP.DmassT_INPUTS, np.ascontiguousarray(rho), np.ascontiguousarray(T),
                      np.array([CP.iP], dtype=np.int32), out, status)
    p_svd = np.where(status == 0, out[:, 0], np.nan)

    p_heos = np.full(n_samples, np.nan)
    p_bic = np.full(n_samples, np.nan)
    for i in range(n_samples):
        for state, arr in ((heos, p_heos), (bic, p_bic)):
            try:
                state.update(CP.DmassT_INPUTS, rho[i], T[i])
                arr[i] = state.p()
            except Exception:  # noqa: BLE001
                # (rho, T) outside this backend's validity (below the
                # melting line, or off BICUBIC's rectangular grid) —
                # keep the pre-filled NaN; _panel masks non-finite cells.
                pass

    # rho_sat,L(T) overlay — shows the anomaly hairpin (peak at ~277 K).
    Tcurve = np.linspace(T_LO, T_HI, 200)
    rho_satL = np.full_like(Tcurve, np.nan)
    for i, Tc in enumerate(Tcurve):
        try:
            heos.update(CP.QT_INPUTS, 0.0, Tc)
            rho_satL[i] = heos.rhomass()
        except Exception:  # noqa: BLE001
            # No saturation point at this T (e.g. above Tc) — leave the
            # overlay NaN so the curve simply doesn't plot there.
            pass
    return T, rho, p_heos, p_svd, p_bic, Tcurve, rho_satL


def _panel(ax, T, rho, p_heos, p_bk, Tcurve, rho_satL, title):
    valid = np.isfinite(p_heos) & np.isfinite(p_bk) & (p_heos != 0.0)
    pct = np.full_like(p_heos, np.nan)
    pct[valid] = np.abs(p_bk[valid] - p_heos[valid]) / np.abs(p_heos[valid]) * 100.0
    vals = np.maximum(pct, 1e-10)
    finite = np.isfinite(vals)
    sc = ax.scatter(T[finite], rho[finite], c=vals[finite], s=6, edgecolor="none",
                    cmap="jet", norm=LogNorm(vmin=1e-8, vmax=1e2))
    neg = np.isfinite(p_heos) & (p_heos > 0) & np.isfinite(p_bk) & (p_bk <= 0)
    if neg.any():
        ax.scatter(T[neg], rho[neg], s=22, marker="x", color="black", linewidths=1.1, label="P < 0")
    ax.plot(Tcurve, rho_satL, color="white", lw=1.4, label=r"$\rho_{sat,L}(T)$")
    ax.axvline(277.13, color="white", ls="--", lw=0.8, alpha=0.7, label=r"$T_{anom}\approx277$ K")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.85)
    ax.set_xlim(T_LO, T_HI)
    ax.set_ylim(RHO_LO, RHO_HI)
    ax.set_xlabel("Temperature  [K]")
    ax.set_ylabel("density  [kg/m³]")
    ax.set_title(f"{title}  ({int(neg.sum())} cells P<0)")
    plt.colorbar(sc, ax=ax, label="Error in pressure  [%]")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=50000)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--interactive", action="store_true",
                    help="open a live GUI window (zoom/pan/inspect) instead of only writing the PNG")
    args = ap.parse_args(argv)
    try:
        T, rho, p_heos, p_svd, p_bic, Tcurve, rho_satL = _sweep(args.n)
    except ImportError:
        sys.stderr.write("CoolProp wrapper not importable — build from this branch first.\n")
        return 1
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
    fig.suptitle(
        f"Water P(ρ,T) error vs HEOS near the 4 °C density anomaly, sat → melting "
        f"(CoolProp-i7j, N={args.n})", fontsize=13)
    _panel(axes[0], T, rho, p_heos, p_bic, Tcurve, rho_satL, "BICUBIC&HEOS")
    _panel(axes[1], T, rho, p_heos, p_svd, Tcurve, rho_satL, "SVDSBTL&HEOS (DT-indexed)")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140, bbox_inches="tight")
    print(f"Wrote {args.out}")
    if args.interactive:
        # Live window: matplotlib's nav toolbar gives zoom/pan/cursor
        # readout.  Blocks until the window is closed.
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
