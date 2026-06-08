#!/usr/bin/env python3
"""Per-fluid error-scatter plots for Phase 2a validation CSVs.

Reads /tmp/svd_sbtl_e2e_*.csv produced by dev/svd_sbtl_e2e and writes
one (h, p) scatter PNG per fluid in the same jet-LogNorm style as the
Python PoC at dev/svd_sbtl_error_plot.py.  Also prints a per-fluid
summary and the multi-fluid aggregate to stdout.

Usage:
    python3 dev/svd_sbtl_e2e_plot.py [CSV_DIR]

CSV_DIR defaults to /tmp.  PNGs land alongside the CSVs.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
from typing import Optional

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

# CoolProp is optional at plot time -- needed only for the saturation-dome
# overlay.  If it's unavailable the script still produces the per-fluid
# error-scatter PNGs (just without the red dome curve).
try:
    import CoolProp.CoolProp as CP
except ImportError:
    CP = None


def fluid_from_path(path: str) -> Optional[str]:
    m = re.match(r"svd_sbtl_e2e_(.+)\.csv$", os.path.basename(path))
    return m.group(1) if m else None


def load_csv(path: str) -> dict:
    """Load a Phase 2a results CSV.  Raises ValueError on malformed input."""
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    if arr.size == 0:
        raise ValueError(f"CSV has no data rows: {path}")
    if arr.ndim == 1:  # single row
        arr = arr[np.newaxis, :]
    if arr.shape[1] < 6:
        raise ValueError(f"CSV must have at least 6 columns: {path}")
    return {
        "T": arr[:, 0],
        "p": arr[:, 1],
        "h": arr[:, 2],
        "rho_truth": arr[:, 3],
        "rho_pred": arr[:, 4],
        "rel_err": arr[:, 5],
    }


def plot_fluid(fluid: str, csv_path: str) -> dict:
    data = load_csv(csv_path)
    p = data["p"]
    h = data["h"]
    rel = np.maximum(data["rel_err"], 1e-12) * 100.0  # percent, with floor

    max_pct = float(rel.max())
    p99 = float(np.percentile(rel, 99))
    p50 = float(np.percentile(rel, 50))

    # Saturation curve for context.  Skip cleanly if CoolProp isn't
    # importable or if the fluid is unknown to PropsSI -- the plot is
    # still useful without the dome overlay.
    psat = None
    hLs = None
    hVs = None
    if CP is not None:
        try:
            Tt = CP.PropsSI("Ttriple", fluid)
            Tc = CP.PropsSI("Tcrit", fluid)
            Tsat = np.linspace(Tt + 0.1, Tc - 0.01, 400)
            psat = CP.PropsSI("P", "T", Tsat, "Q", 0, fluid)
            hLs = CP.PropsSI("Hmass", "T", Tsat, "Q", 0, fluid)
            hVs = CP.PropsSI("Hmass", "T", Tsat, "Q", 1, fluid)
        except Exception:
            psat = None
            hLs = None
            hVs = None

    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    cnorm = colors.LogNorm(vmin=1e-10, vmax=10)
    ax.scatter(h / 1e3, p, c=rel, s=4, cmap="jet", norm=cnorm, edgecolors="none")
    if psat is not None:
        ax.plot(hLs / 1e3, psat, "k-", lw=1.0, alpha=0.8)
        ax.plot(hVs / 1e3, psat, "k-", lw=1.0, alpha=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("h [kJ/kg]")
    ax.set_ylabel("p [Pa]")
    ax.set_title(
        f"{fluid} HmassP ρ rel error: SVDSBTL e2e (C++ port)\n"
        f"N={len(rel)}  med={p50:.1e}  p99={p99:.1e}  max={max_pct:.1e} %"
    )
    plt.colorbar(
        plt.cm.ScalarMappable(norm=cnorm, cmap="jet"),
        ax=ax,
        label=r"$|\rho_{\rm pred}/\rho_{\rm EOS} - 1| \times 100\,[\%]$",
        shrink=0.85,
    )
    plt.tight_layout()
    png_path = csv_path.replace(".csv", ".png")
    fig.savefig(png_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return {
        "fluid": fluid,
        "n": int(len(rel)),
        "max_pct": max_pct,
        "p99_pct": p99,
        "p50_pct": p50,
        "png_path": png_path,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_dir", nargs="?", default="/tmp")
    args = ap.parse_args()

    csv_glob = os.path.join(args.csv_dir, "svd_sbtl_e2e_*.csv")
    csvs = sorted(glob.glob(csv_glob))
    if not csvs:
        raise SystemExit(f"no CSVs match {csv_glob}")

    rows = []
    for csv in csvs:
        fluid = fluid_from_path(csv)
        if fluid is None or os.path.getsize(csv) == 0:
            continue
        try:
            rows.append(plot_fluid(fluid, csv))
        except Exception as exc:
            # One malformed CSV must not abort the whole batch.
            print(f"skip {csv}: {exc}")

    print(f"\n{'fluid':<12s} {'N':>6s} {'max %':>12s} {'p99 %':>12s} {'p50 %':>12s}  png")
    print("-" * 90)
    for r in rows:
        print(
            f"{r['fluid']:<12s} {r['n']:>6d} {r['max_pct']:>12.3e} "
            f"{r['p99_pct']:>12.3e} {r['p50_pct']:>12.3e}  {r['png_path']}"
        )


if __name__ == "__main__":
    main()
