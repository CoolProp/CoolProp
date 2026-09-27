#!/usr/bin/env python3
"""Where in (T, p) the updated CoolProp PT flash is faster or slower than REFPROP 10 TPFLSH.

Reads the per-state CSV(s) from bench_ratio.cpp; one panel per mixture, each state colored by
t_CoolProp / t_TPFLSH on a log diverging scale (blue: CoolProp faster, red: slower).  States CoolProp
publishes as two-phase are outlined dark, so the phase envelope shows.

    python3 plot_ratio_pT.py out.png ratio_on.csv [ratio_more.csv ...]
"""
import csv
import math
import sys
from collections import OrderedDict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

# dataviz reference palette: diverging blue <-> red, neutral gray midpoint, light chart surface
BLUE, RED, MID, SURFACE, INK, MUTED, GRID = "#2a78d6", "#e34948", "#f0efec", "#fcfcfb", "#1f1f1e", "#6b6a64", "#e6e5df"
CMAP = LinearSegmentedColormap.from_list("div", [BLUE, MID, RED])
LIM = math.log10(30.0)  # clip the color scale at 30x either way


def load(files):
    mixes = OrderedDict()
    for fn in files:
        for r in csv.DictReader(open(fn)):
            mixes.setdefault(r["mixture"], []).append(r)
    return mixes


def main(out, files):
    mixes = load(files)
    n = len(mixes)
    nstates = max(len(v) for v in mixes.values())
    ncol = 3
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.1 * nrow + 0.9), squeeze=False)
    fig.subplots_adjust(hspace=0.5, wspace=0.18, bottom=0.17, top=0.86)
    fig.patch.set_facecolor(SURFACE)
    norm = Normalize(-LIM, LIM)
    for k, (name, rows) in enumerate(mixes.items()):
        ax = axes[k // ncol][k % ncol]
        ax.set_facecolor(SURFACE)
        pts, bad_cp, bad_rp = [], 0, 0
        for r in rows:
            tcp, trp = float(r["t_cp_us"]), float(r["t_rp_us"])
            fcp, frp = int(r["cp_fail"]) != 0, int(r["ierr_rp"]) > 0
            bad_cp += fcp
            bad_rp += frp
            if fcp or frp or not (tcp > 0 and trp > 0):
                continue
            Q = float(r["Q_cp"])
            pts.append((float(r["T"]), float(r["p"]) / 1e6, math.log10(tcp / trp), 0 < Q < 1))
        # marker size shrinks with density so a 10k-state map still reads as points
        sc = math.sqrt(2000.0 / max(len(rows), 1))
        # single phase first, two-phase on top
        for two in (False, True):
            sel = [q for q in pts if q[3] == two]
            if not sel:
                continue
            ax.scatter([q[0] for q in sel], [q[1] for q in sel], c=[max(-LIM, min(LIM, q[2])) for q in sel], cmap=CMAP, norm=norm,
                       s=(11 if two else 9) * sc, linewidths=(0.6 if two else 0.25) * math.sqrt(sc), edgecolors=INK if two else "#b8b7b0", zorder=3 if two else 2)
        lr = sorted(q[2] for q in pts)
        med = 10 ** lr[len(lr) // 2]
        faster = sum(1 for v in lr if v < 0) / len(lr)
        medtxt = f"{1 / med:.1f}× faster" if med < 1 else f"{med:.1f}× slower"
        ax.set_title(f"{name}\nmedian {medtxt} · faster at {100 * faster:.0f}% of states", fontsize=9.5, color=INK, loc="left")
        ax.set_yscale("log")
        ax.set_ylim(0.01, 30)
        ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
        for s in ax.spines.values():
            s.set_color("#d4d3cc")
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.set_xlabel("T / K", fontsize=9, color=INK)
        if k % ncol == 0:
            ax.set_ylabel("p / MPa", fontsize=9, color=INK)
        if bad_cp or bad_rp:
            ax.text(0.98, 0.03, f"not plotted: REFPROP failed {bad_rp}, CoolProp failed {bad_cp}", transform=ax.transAxes, ha="right",
                    fontsize=7, color=MUTED, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5))
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    cax = fig.add_axes([0.2, 0.06, 0.6, 0.022])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    ticks = [-math.log10(30), -1, -math.log10(3), 0, math.log10(3), 1, math.log10(30)]
    cb.set_ticks(ticks)
    cb.set_ticklabels(["30× faster", "10× faster", "3× faster", "equal", "3× slower", "10× slower", "30× slower"])
    cb.ax.tick_params(labelsize=8, colors=INK)
    cb.outline.set_edgecolor("#d4d3cc")
    cb.set_label("CoolProp PT flash time / REFPROP 10 TPFLSH time", fontsize=9, color=INK)
    fig.suptitle(f"Updated CoolProp PT flash vs REFPROP 10 TPFLSH, per state ({nstates} states per mixture, min of 3 timings each{', mixtures run as parallel processes' if nstates > 2000 else ''})\n"
                 "CoolProp: master + #3427 SS-skip + Chebyshev density kernel.  Dark outline: CoolProp publishes a two-phase state.\n"
                 "GERG-2008 on both sides, except R454B (CoolProp HEOS vs REFPROP default mixture model).",
                 fontsize=10, color=INK, x=0.06, ha="left", y=0.975)
    fig.savefig(out, dpi=140, bbox_inches="tight", facecolor=SURFACE)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
