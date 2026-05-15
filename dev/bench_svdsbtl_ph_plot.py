"""Plot SVDSBTL PH-benchmark CSVs as state-point heatmaps.

Reads /tmp/bench_svdsbtl_ph_<fluid>.csv produced by
dev/bench_svdsbtl_ph.cpp and emits a single PNG per fluid laid out as:

  Row 0 (accuracy):  rho rel-err     T rel-err     s rel-err
  Row 1 (timing):    SVDSBTL ns      HEOS ns       REFPROP ns

The three timing panels share a single log color scale so the eye
can directly compare backends.  REFPROP columns are blank if the
benchmark was run without BENCH_REFPROP_PATH.

Usage:
    python3 dev/bench_svdsbtl_ph_plot.py
    python3 dev/bench_svdsbtl_ph_plot.py Water Propane
    python3 dev/bench_svdsbtl_ph_plot.py /tmp/bench_svdsbtl_ph_Water.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm


def _err_panel(ax, h_kJ, p_MPa, vals, title, vmin=1e-7, vmax=1e-1):
    vals = np.maximum(vals, 1e-16)
    sc = ax.scatter(h_kJ, p_MPa, c=vals, s=6, edgecolor="none", cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax))
    ax.set_yscale("log")
    ax.set_xlabel("h  [kJ/kg]")
    ax.set_ylabel("p  [MPa]")
    ax.set_title(title)
    plt.colorbar(sc, ax=ax, label="relative error")


def _time_panel(ax, h_kJ, p_MPa, ns, title, norm):
    if np.all(np.isnan(ns)):
        ax.set_axis_off()
        ax.set_title(title + "  (n/a)")
        return
    valid = ~np.isnan(ns)
    sc = ax.scatter(h_kJ[valid], p_MPa[valid], c=ns[valid], s=6, edgecolor="none", cmap="magma", norm=norm)
    ax.set_yscale("log")
    ax.set_xlabel("h  [kJ/kg]")
    ax.set_ylabel("p  [MPa]")
    ax.set_title(title)
    plt.colorbar(sc, ax=ax, label="ns / call")


def plot_one(csv_path: Path, out_path: Path) -> None:
    df = pd.read_csv(csv_path)
    if df.empty:
        print(f"  {csv_path.name}: no rows, skipping")
        return

    fluid = csv_path.stem.replace("bench_svdsbtl_ph_", "")
    has_refprop = df["ns_per_call_refprop"].notna().any()
    has_refprop_direct = "ns_per_call_refprop_direct" in df.columns and df["ns_per_call_refprop_direct"].notna().any()

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    fig.suptitle(
        f"SVDSBTL benchmark in PH coordinates -- {fluid}   N={len(df)}",
        fontsize=13,
        y=1.02,
    )

    h_kJ = df["h"].to_numpy() / 1e3
    p_MPa = df["p"].to_numpy() / 1e6

    # Row 0 -- accuracy panels (SVDSBTL vs HEOS).
    _err_panel(axes[0, 0], h_kJ, p_MPa, df["rel_err_rho"].to_numpy(), r"$\rho$  rel. error  vs HEOS")
    _err_panel(axes[0, 1], h_kJ, p_MPa, df["rel_err_T"].to_numpy(), "T  rel. error  vs HEOS")
    _err_panel(axes[0, 2], h_kJ, p_MPa, df["rel_err_s"].to_numpy(), "s  rel. error  vs HEOS")

    # Row 1 -- timing panels on a SHARED log color scale so the
    # cross-backend speed delta is visually unambiguous.  Prefer the
    # direct PHFLSHdll measurement for the REFPROP panel when we have
    # it; the via-AbstractState run for REFPROP only adds a couple
    # hundred ns of wrapper overhead on top, which is invisible at
    # this scale.
    refprop_ns = df["ns_per_call_refprop_direct"].to_numpy() if has_refprop_direct else (df["ns_per_call_refprop"].to_numpy() if has_refprop else np.full(len(df), np.nan))
    timings = [df["ns_per_call_svd"].to_numpy(), df["ns_per_call_heos"].to_numpy(), refprop_ns]
    all_ns = np.concatenate([t[~np.isnan(t)] for t in timings])
    vmin = max(50.0, float(np.percentile(all_ns, 1)))
    vmax = float(np.percentile(all_ns, 99))
    shared = LogNorm(vmin=vmin, vmax=vmax)

    _time_panel(axes[1, 0], h_kJ, p_MPa, df["ns_per_call_svd"].to_numpy(), "SVDSBTL  ns / (update + rhomass)", shared)
    _time_panel(axes[1, 1], h_kJ, p_MPa, df["ns_per_call_heos"].to_numpy(), "HEOS  ns / (update + rhomass)", shared)
    refprop_title = "REFPROP direct PHFLSHdll  ns / call" if has_refprop_direct else "REFPROP  ns / (update + rhomass via AbstractState)"
    _time_panel(axes[1, 2], h_kJ, p_MPa, refprop_ns, refprop_title, shared)

    # Summary footer.
    max_rho = df["rel_err_rho"].max()
    p99_rho = df["rel_err_rho"].quantile(0.99)
    med_rho = df["rel_err_rho"].median()
    mean_ns_svd = df["ns_per_call_svd"].mean()
    mean_ns_heos = df["ns_per_call_heos"].mean()
    parts = [
        f"rho rel-err: max={max_rho:.2e}  p99={p99_rho:.2e}  median={med_rho:.2e}",
        f"mean ns/call: SVDSBTL={mean_ns_svd:.0f}  HEOS={mean_ns_heos:.0f}",
        f"speedup vs HEOS={mean_ns_heos / mean_ns_svd:.1f}x",
    ]
    if has_refprop_direct:
        mean_ns_refprop_direct = df["ns_per_call_refprop_direct"].mean()
        mean_ns_refprop_as = df["ns_per_call_refprop"].mean() if has_refprop else float("nan")
        parts[1] += f"  REFPROP(direct PHFLSHdll)={mean_ns_refprop_direct:.0f}"
        if not np.isnan(mean_ns_refprop_as):
            parts[1] += f"  REFPROP(via AS)={mean_ns_refprop_as:.0f}"
        parts.append(f"speedup vs REFPROP(direct)={mean_ns_refprop_direct / mean_ns_svd:.1f}x")
    elif has_refprop:
        mean_ns_refprop = df["ns_per_call_refprop"].mean()
        parts[1] += f"  REFPROP={mean_ns_refprop:.0f}"
        parts.append(f"speedup vs REFPROP={mean_ns_refprop / mean_ns_svd:.1f}x")
    fig.text(0.5, -0.03, "    ".join(parts), ha="center", fontsize=10, family="monospace")

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  {csv_path.name}  ->  {out_path}")


def main() -> int:
    args = sys.argv[1:]
    if not args:
        args = ["Water"]

    paths: list[Path] = []
    for a in args:
        p = Path(a)
        if p.suffix == ".csv" and p.exists():
            paths.append(p)
        else:
            paths.append(Path(f"/tmp/bench_svdsbtl_ph_{a}.csv"))

    out_dir = Path("dev")
    for p in paths:
        if not p.exists():
            print(f"  MISSING: {p}")
            continue
        out_path = out_dir / f"bench_svdsbtl_ph_{p.stem.replace('bench_svdsbtl_ph_', '')}.png"
        plot_one(p, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
