"""Plot SVDSBTL PH-benchmark CSVs as state-point heatmaps.

Reads /tmp/bench_svdsbtl_ph_<fluid>.csv produced by
dev/bench_svdsbtl_ph.cpp and emits a single PNG per fluid laid out as:

  Row 0 (accuracy):  rho rel-err  T rel-err  s rel-err  u rel-err
  Row 1 (timing):    SVDSBTL ns   IF97 ns    HEOS ns    REFPROP ns

The four timing panels share a single log color scale so the eye
can directly compare backends.  IF97 / REFPROP panels are blank if
the run didn't have them available.

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


def _err_panel(ax, h_kJ, p_MPa, vals, title, vmin=1e-7, vmax=1e-1, dome=None, p_crit_MPa=None):
    vals = np.maximum(np.asarray(vals, dtype=float), 1e-16)
    valid = ~np.isnan(vals)
    sc = ax.scatter(h_kJ[valid], p_MPa[valid], c=vals[valid], s=4, edgecolor="none", cmap="viridis", norm=LogNorm(vmin=vmin, vmax=vmax))
    ax.set_yscale("log")
    ax.set_xlabel("h  [kJ/kg]")
    ax.set_ylabel("p  [MPa]")
    ax.set_title(title)
    if dome is not None:
        ax.plot(dome["hL"] / 1e3, dome["p"] / 1e6, color="red", lw=0.8, alpha=0.8)
        ax.plot(dome["hV"] / 1e3, dome["p"] / 1e6, color="red", lw=0.8, alpha=0.8)
    if p_crit_MPa is not None:
        ax.axhline(p_crit_MPa, color="gray", ls="--", lw=0.6, alpha=0.6)
    plt.colorbar(sc, ax=ax, label="relative error")


def _time_panel(ax, h_kJ, p_MPa, ns, title, norm, dome=None, p_crit_MPa=None):
    if np.all(np.isnan(ns)):
        ax.set_axis_off()
        ax.set_title(title + "  (n/a)")
        return
    valid = ~np.isnan(ns)
    sc = ax.scatter(h_kJ[valid], p_MPa[valid], c=ns[valid], s=4, edgecolor="none", cmap="magma", norm=norm)
    ax.set_yscale("log")
    ax.set_xlabel("h  [kJ/kg]")
    ax.set_ylabel("p  [MPa]")
    ax.set_title(title)
    if dome is not None:
        ax.plot(dome["hL"] / 1e3, dome["p"] / 1e6, color="cyan", lw=0.8, alpha=0.8)
        ax.plot(dome["hV"] / 1e3, dome["p"] / 1e6, color="cyan", lw=0.8, alpha=0.8)
    if p_crit_MPa is not None:
        ax.axhline(p_crit_MPa, color="gray", ls="--", lw=0.6, alpha=0.6)
    plt.colorbar(sc, ax=ax, label="ns / call")


def _saturation_dome(fluid: str) -> tuple[dict | None, float | None]:
    """Walk the saturation curve via CoolProp.HEOS for the dome overlay.

    Returns (dict with arrays p, hL, hV, or None if HEOS doesn't have
    this fluid; p_crit in Pa or None).  Catches all exceptions so the
    plot still produces something useful if CoolProp isn't importable.
    """
    try:
        import CoolProp.CoolProp as cp
    except Exception:
        return None, None

    try:
        p_crit = cp.PropsSI("pcrit", fluid)
        p_trip = cp.PropsSI("ptriple", fluid)
    except Exception:
        return None, None

    p_arr = np.geomspace(p_trip * 1.01, p_crit * 0.999, 200)
    hL = np.full(p_arr.size, np.nan)
    hV = np.full(p_arr.size, np.nan)
    for i, p in enumerate(p_arr):
        try:
            hL[i] = cp.PropsSI("H", "P", p, "Q", 0.0, fluid)
            hV[i] = cp.PropsSI("H", "P", p, "Q", 1.0, fluid)
        except Exception:
            pass
    return {"p": p_arr, "hL": hL, "hV": hV}, p_crit


def plot_one(csv_path: Path, out_path: Path) -> None:
    df = pd.read_csv(csv_path)
    if df.empty:
        print(f"  {csv_path.name}: no rows, skipping")
        return

    fluid = csv_path.stem.replace("bench_svdsbtl_ph_", "")
    dome, p_crit_Pa = _saturation_dome(fluid)
    p_crit_MPa = p_crit_Pa / 1e6 if p_crit_Pa is not None else None
    has_refprop = df["ns_per_call_refprop"].notna().any()
    has_refprop_direct = "ns_per_call_refprop_direct" in df.columns and df["ns_per_call_refprop_direct"].notna().any()
    has_refprop_phfl1 = "ns_per_call_refprop_phfl1" in df.columns and df["ns_per_call_refprop_phfl1"].notna().any()
    has_if97 = "ns_per_call_if97" in df.columns and df["ns_per_call_if97"].notna().any()
    has_if97_direct = "ns_per_call_if97_direct" in df.columns and df["ns_per_call_if97_direct"].notna().any()

    fig, axes = plt.subplots(2, 4, figsize=(19, 9), constrained_layout=True)
    fig.suptitle(
        f"SVDSBTL benchmark in PH coordinates -- {fluid}   N={len(df)}",
        fontsize=13,
        y=1.02,
    )

    h_kJ = df["h"].to_numpy() / 1e3
    p_MPa = df["p"].to_numpy() / 1e6

    # Row 0 -- accuracy panels (SVDSBTL vs HEOS).  Sat dome (red) +
    # critical-pressure line (gray dashed) overlaid so the eye can
    # see where coverage stops and where the surface is fighting
    # boundary geometry.
    common = dict(dome=dome, p_crit_MPa=p_crit_MPa)
    _err_panel(axes[0, 0], h_kJ, p_MPa, df["rel_err_rho"].to_numpy(), r"$\rho$  rel. error  vs HEOS", **common)
    _err_panel(axes[0, 1], h_kJ, p_MPa, df["rel_err_T"].to_numpy(), "T  rel. error  vs HEOS", **common)
    _err_panel(axes[0, 2], h_kJ, p_MPa, df["rel_err_s"].to_numpy(), "s  rel. error  vs HEOS", **common)
    _err_panel(axes[0, 3], h_kJ, p_MPa, df["rel_err_u"].to_numpy(), "u  rel. error  vs HEOS", **common)

    # Row 1 -- timing panels on a SHARED log color scale so the
    # cross-backend speed delta is visually unambiguous.  For IF97
    # and REFPROP we prefer the direct-call measurement (header-only
    # IF97::rhomass_phmass / dlsym'd PHFLSHdll) when available;
    # via-AbstractState measurements get pulled in as fallback.
    if97_ns = (
        df["ns_per_call_if97_direct"].to_numpy() if has_if97_direct
        else (df["ns_per_call_if97"].to_numpy() if has_if97 else np.full(len(df), np.nan))
    )
    # Prefer PHFL1 (single-phase fast path, apples-to-apples with
    # SVDSBTL which is also single-phase-only) over PHFLSH for the
    # REFPROP timing panel.  Fall back to PHFLSH then via-AS.
    if has_refprop_phfl1:
        refprop_ns = df["ns_per_call_refprop_phfl1"].to_numpy()
        refprop_panel_title = "REFPROP direct PHFL1dll  ns / call"
    elif has_refprop_direct:
        refprop_ns = df["ns_per_call_refprop_direct"].to_numpy()
        refprop_panel_title = "REFPROP direct PHFLSHdll  ns / call"
    elif has_refprop:
        refprop_ns = df["ns_per_call_refprop"].to_numpy()
        refprop_panel_title = "REFPROP  ns / (update + rhomass via AS)"
    else:
        refprop_ns = np.full(len(df), np.nan)
        refprop_panel_title = "REFPROP"
    timings = [df["ns_per_call_svd"].to_numpy(), df["ns_per_call_heos"].to_numpy(), if97_ns, refprop_ns]
    valid_arrays = [t[~np.isnan(t)] for t in timings if not np.all(np.isnan(t))]
    all_ns = np.concatenate(valid_arrays) if valid_arrays else np.array([100.0, 10000.0])
    vmin = max(50.0, float(np.percentile(all_ns, 1)))
    vmax = float(np.percentile(all_ns, 99))
    shared = LogNorm(vmin=vmin, vmax=vmax)

    tcommon = dict(dome=dome, p_crit_MPa=p_crit_MPa)
    _time_panel(axes[1, 0], h_kJ, p_MPa, df["ns_per_call_svd"].to_numpy(), "SVDSBTL  ns / (update + rhomass)", shared, **tcommon)
    if97_title = "IF97 direct rhomass_phmass  ns / call" if has_if97_direct else "IF97  ns / (update + rhomass via AS)"
    _time_panel(axes[1, 1], h_kJ, p_MPa, if97_ns, if97_title, shared, **tcommon)
    _time_panel(axes[1, 2], h_kJ, p_MPa, df["ns_per_call_heos"].to_numpy(), "HEOS  ns / (update + rhomass)", shared, **tcommon)
    _time_panel(axes[1, 3], h_kJ, p_MPa, refprop_ns, refprop_panel_title, shared, **tcommon)

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
    if has_if97_direct:
        mean_ns_if97_direct = df["ns_per_call_if97_direct"].mean()
        mean_ns_if97_as = df["ns_per_call_if97"].mean() if has_if97 else float("nan")
        parts[1] += f"  IF97(direct)={mean_ns_if97_direct:.0f}"
        if not np.isnan(mean_ns_if97_as):
            parts[1] += f"  IF97(via AS)={mean_ns_if97_as:.0f}"
        parts.append(f"speedup vs IF97(direct)={mean_ns_if97_direct / mean_ns_svd:.2f}x")
    elif has_if97:
        mean_ns_if97 = df["ns_per_call_if97"].mean()
        parts[1] += f"  IF97={mean_ns_if97:.0f}"
        parts.append(f"speedup vs IF97={mean_ns_if97 / mean_ns_svd:.2f}x")
    if has_refprop_direct or has_refprop_phfl1:
        mean_ns_refprop_direct = df["ns_per_call_refprop_direct"].mean() if has_refprop_direct else float("nan")
        mean_ns_refprop_phfl1 = df["ns_per_call_refprop_phfl1"].mean() if has_refprop_phfl1 else float("nan")
        bits = []
        if not np.isnan(mean_ns_refprop_phfl1):
            bits.append(f"PHFL1={mean_ns_refprop_phfl1:.0f}")
        if not np.isnan(mean_ns_refprop_direct):
            bits.append(f"PHFLSH={mean_ns_refprop_direct:.0f}")
        parts[1] += "  REFPROP(" + ", ".join(bits) + ")"
        ref_mean = mean_ns_refprop_phfl1 if not np.isnan(mean_ns_refprop_phfl1) else mean_ns_refprop_direct
        ref_label = "PHFL1" if not np.isnan(mean_ns_refprop_phfl1) else "PHFLSH"
        parts.append(f"speedup vs REFPROP({ref_label})={ref_mean / mean_ns_svd:.1f}x")
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
