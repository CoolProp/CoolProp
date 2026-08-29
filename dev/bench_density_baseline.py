"""
Plot and summarize bench_density_baseline.cpp output across systems.

Reads:
  dev/bench_density_<system>.csv

Writes:
  dev/bench_density_<system>.png   — 2x2 panel: phase envelope, flash time,
                                     solver_rho_Tp time, tier markers
  stdout                           — markdown failure-mode summary

Run:
    python dev/bench_density_baseline.py [system1 system2 ...]
"""

import os
import sys
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import CoolProp

DEV = os.path.dirname(os.path.abspath(__file__))

# Must match SYSTEMS in bench_density_baseline.cpp
SYSTEMS = {
    "amarillo": dict(
        fluids=("Methane&Nitrogen&CarbonDioxide&Ethane&Propane"
                "&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane"),
        z=[0.906724, 0.031284, 0.004676, 0.045279, 0.008280,
           0.001037, 0.001563, 0.000321, 0.000443, 0.000393],
        NT=60, NP=60,
        T_min=140.0, T_max=260.0,
        p_min=0.5e5, p_max=80e5,
        title="Amarillo AGA-8 natural gas (10-component)",
    ),
    "ch4h2s": dict(
        fluids="Methane&HydrogenSulfide",
        z=[0.70, 0.30],
        NT=60, NP=60,
        T_min=200.0, T_max=420.0,
        p_min=1e5, p_max=150e5,
        title="CH4/H2S 70/30 (sour gas, H2S near-critical at 373 K)",
    ),
    "co2n2": dict(
        fluids="CarbonDioxide&Nitrogen",
        z=[0.80, 0.20],
        NT=60, NP=60,
        T_min=250.0, T_max=400.0,
        p_min=10e5, p_max=250e5,
        title="CO2/N2 80/20 (cuts CO2 critical 304 K, 73.8 bar)",
    ),
    "c10c1": dict(
        fluids="n-Decane&Methane",
        z=[0.30, 0.70],
        NT=60, NP=60,
        T_min=250.0, T_max=500.0,
        p_min=1e5, p_max=500e5,
        title="n-Decane/Methane 30/70 (low-T heavy + light, retrograde)",
    ),
}


def trace_envelope(fluids, z, T_min, T_max, n=120, p_clip=(1e3, 1e9)):
    Ts = np.linspace(T_min, T_max, n)
    out = {0: ([], []), 1: ([], [])}
    for Q in (0, 1):
        for T in Ts:
            try:
                AS = CoolProp.AbstractState("HEOS", fluids)
                AS.set_mole_fractions(z)
                AS.update(CoolProp.QT_INPUTS, float(Q), float(T))
                p = AS.p()
                if p_clip[0] < p < p_clip[1]:
                    out[Q][0].append(T)
                    out[Q][1].append(p)
            except Exception:
                pass
    return ((np.asarray(out[0][0]), np.asarray(out[0][1])),
            (np.asarray(out[1][0]), np.asarray(out[1][1])))


def load_grid(path, NT, NP):
    raw = np.genfromtxt(path, delimiter=",", names=True)
    return {name: raw[name].reshape(NP, NT) for name in raw.dtype.names}


def plot_system(key, spec, csv_path, png_path):
    g = load_grid(csv_path, spec["NT"], spec["NP"])
    TT, PP = g["T_K"], g["p_Pa"]
    P_bar = PP / 1e5
    flash_ms = g["flash_ms"]
    solver_us = g["solver_us"]
    phase = g["phase"].astype(int)
    solver_ok = g["solver_ok"].astype(int)

    err_flash = phase == -1

    valid_flash = flash_ms[~err_flash]
    valid_solver = solver_us[(solver_ok == 1) & np.isfinite(solver_us)]

    (T_bub, p_bub), (T_dew, p_dew) = trace_envelope(
        spec["fluids"], spec["z"], spec["T_min"], spec["T_max"]
    )

    # Stats
    flash_fail = int(err_flash.sum())
    solver_fail = int(((solver_ok == 0) & (~err_flash)).sum())
    n_total = phase.size
    summary = {
        "system": key,
        "title": spec["title"],
        "n_total": n_total,
        "flash_fail": flash_fail,
        "solver_fail": solver_fail,
        "flash_ms_p50": float(np.percentile(valid_flash, 50)) if valid_flash.size else float("nan"),
        "flash_ms_p95": float(np.percentile(valid_flash, 95)) if valid_flash.size else float("nan"),
        "flash_ms_max": float(valid_flash.max()) if valid_flash.size else float("nan"),
        "solver_us_p50": float(np.percentile(valid_solver, 50)) if valid_solver.size else float("nan"),
        "solver_us_p95": float(np.percentile(valid_solver, 95)) if valid_solver.size else float("nan"),
        "solver_us_max": float(valid_solver.max()) if valid_solver.size else float("nan"),
    }

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(spec["title"], fontsize=12)

    def add_envelope(ax):
        if T_bub.size: ax.plot(T_bub, p_bub / 1e5, "w--", lw=1.2, alpha=0.85)
        if T_dew.size: ax.plot(T_dew, p_dew / 1e5, "w-", lw=1.2, alpha=0.85)

    def style(ax, title, txt):
        ax.set_yscale("log")
        ax.set_ylim(spec["p_min"] / 1e5, spec["p_max"] / 1e5)
        ax.set_xlim(spec["T_min"], spec["T_max"])
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Pressure (bar)")
        ax.set_title(title)
        ax.grid(True, alpha=0.2, which="both")
        ax.text(0.02, 0.97, txt, transform=ax.transAxes, fontsize=8,
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.75))

    # [0,0] Phase envelope
    ax = axes[0, 0]
    if T_bub.size: ax.plot(T_bub, p_bub / 1e5, "r--", lw=2, label="Bubble")
    if T_dew.size: ax.plot(T_dew, p_dew / 1e5, "b-", lw=2, label="Dew")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_xlim(spec["T_min"], spec["T_max"])
    if T_bub.size or T_dew.size:
        ax.set_yscale("log")
        ax.set_ylim(spec["p_min"] / 1e5, spec["p_max"] / 1e5)
    ax.set_title("Phase envelope (QT trace)")
    ax.grid(True, alpha=0.3, which="both")
    if T_bub.size or T_dew.size: ax.legend(fontsize=9)

    # [0,1] Full flash time
    ax = axes[0, 1]
    if valid_flash.size:
        norm = mcolors.LogNorm(vmin=max(1e-3, np.percentile(valid_flash, 2)),
                                vmax=np.percentile(valid_flash, 98))
        im = ax.pcolormesh(TT, P_bar, flash_ms, norm=norm, cmap="plasma", shading="auto")
        fig.colorbar(im, ax=ax, label="ms (log)")
    if err_flash.any():
        ax.pcolormesh(TT, P_bar, np.where(err_flash, 1.0, np.nan),
                      cmap=mcolors.ListedColormap(["#ff000088"]),
                      shading="auto", vmin=0, vmax=1)
    add_envelope(ax)
    style(ax, "Full PT flash time",
          f"min  {valid_flash.min():.2f} ms\n"
          f"med  {summary['flash_ms_p50']:.2f} ms\n"
          f"p95  {summary['flash_ms_p95']:.2f} ms\n"
          f"max  {summary['flash_ms_max']:.1f} ms\n"
          f"errors  {flash_fail}/{n_total}")

    # [1,0] Solver_rho_Tp time
    ax = axes[1, 0]
    if valid_solver.size:
        norm = mcolors.LogNorm(vmin=max(1e-2, np.percentile(valid_solver, 2)),
                                vmax=np.percentile(valid_solver, 98))
        # Mask out cells where solver did not run / failed
        plot_us = np.where((solver_ok == 1) & np.isfinite(solver_us), solver_us, np.nan)
        im = ax.pcolormesh(TT, P_bar, plot_us, norm=norm, cmap="viridis", shading="auto")
        fig.colorbar(im, ax=ax, label="μs (log)")
    fail_mask = (solver_ok == 0) & ~err_flash
    if fail_mask.any():
        ax.pcolormesh(TT, P_bar, np.where(fail_mask, 1.0, np.nan),
                      cmap=mcolors.ListedColormap(["#ff000088"]),
                      shading="auto", vmin=0, vmax=1)
    add_envelope(ax)
    style(ax, "Isolated solver_rho_Tp time (per call)",
          f"min  {valid_solver.min():.2f} μs\n"
          f"med  {summary['solver_us_p50']:.2f} μs\n"
          f"p95  {summary['solver_us_p95']:.2f} μs\n"
          f"max  {summary['solver_us_max']:.1f} μs\n"
          f"failures  {solver_fail}/{n_total - flash_fail}")

    # [1,1] Phase classification
    ax = axes[1, 1]
    phase_for_plot = np.where(err_flash, np.nan, phase.astype(float))
    im = ax.pcolormesh(TT, P_bar, phase_for_plot, cmap="tab20",
                       shading="auto", vmin=0, vmax=10)
    fig.colorbar(im, ax=ax, label="CoolProp phase index")
    add_envelope(ax)
    phase_counts = Counter(phase[~err_flash].tolist())
    phase_txt = "  ".join(f"{p}:{c}" for p, c in sorted(phase_counts.items()))
    style(ax, "Phase classification", phase_txt or "(no data)")

    plt.tight_layout()
    plt.savefig(png_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return summary


def emit_markdown(summaries):
    print("# Baseline characterization — solver_rho_Tp\n")
    print("Issue: CoolProp-aor\n")
    print("| System | Cells | Flash err | Solver fail | Flash p50/p95/max (ms) | Solver p50/p95/max (μs) |")
    print("|---|---:|---:|---:|---|---|")
    for s in summaries:
        print(f"| {s['system']} | {s['n_total']} | {s['flash_fail']} | {s['solver_fail']} "
              f"| {s['flash_ms_p50']:.2f} / {s['flash_ms_p95']:.2f} / {s['flash_ms_max']:.1f} "
              f"| {s['solver_us_p50']:.2f} / {s['solver_us_p95']:.2f} / {s['solver_us_max']:.1f} |")
    print()
    for s in summaries:
        print(f"## {s['system']} — {s['title']}")
        print(f"- PNG: `dev/bench_density_{s['system']}.png`")
        print(f"- Cells: {s['n_total']}, flash errors: {s['flash_fail']}, "
              f"solver failures (within successful flashes): {s['solver_fail']}")
        print()


def main():
    keys = sys.argv[1:] if len(sys.argv) > 1 else list(SYSTEMS.keys())
    summaries = []
    for key in keys:
        if key not in SYSTEMS:
            print(f"# Unknown system: {key}", file=sys.stderr)
            continue
        csv_path = os.path.join(DEV, f"bench_density_{key}.csv")
        if not os.path.exists(csv_path):
            print(f"# Missing CSV: {csv_path} — skip", file=sys.stderr)
            continue
        png_path = os.path.join(DEV, f"bench_density_{key}.png")
        summaries.append(plot_system(key, SYSTEMS[key], csv_path, png_path))
    if summaries:
        emit_markdown(summaries)


if __name__ == "__main__":
    main()
