#!/usr/bin/env python3
"""Where in (T, p) the updated CoolProp PT flash is faster or slower than REFPROP 10 TPFLSH.

Reads the per-state CSV(s) from bench_ratio.cpp; one panel per mixture, each state colored by
t_CoolProp / t_TPFLSH on a log diverging scale (blue: CoolProp faster, red: slower).  States CoolProp
publishes as two-phase are outlined dark, so the phase envelope shows.

    python3 plot_ratio_pT.py out.png ratio_on.csv [ratio_more.csv ...]
    python3 plot_ratio_pT.py out.png a.csv b.csv -- c.csv d.csv      # '--' starts a new row of panels
    ... v:verdicts_a.csv v:verdicts_b.csv                             # overlay per-state verdicts (tools/verdict.cpp)
    ... build:"<one-line description of the CoolProp build>"          # names the build in the title (required)
"""
import csv
import math
from math import lcm
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


# Verdict markers: CoolProp wrong / knife-edge in ink; TPFLSH wrong / minor / out of scope smaller and lighter.
# Every verdict is checked with REFPROP's own GERG routines (no CoolProp code), see tools/README.md.
VMARK = {
    "cp_wrong": dict(marker="X", s=70, c=INK, edgecolors=SURFACE, linewidths=0.8, zorder=6, label="CoolProp wrong: misses or mis-splits a verified split"),
    "cp_history": dict(marker="D", s=34, facecolors="none", edgecolors=INK, linewidths=1.1, zorder=6, label="CoolProp verdict flips under a ~1e-10 input change"),
    "rp_wrong": dict(marker="^", s=9, facecolors="none", edgecolors="#5c5b55", linewidths=0.45, zorder=5, label="TPFLSH wrong: split off equilibrium by >=1e-3 (its own ln f, p)"),
    "rp_minor": dict(marker="o", s=6, facecolors="none", edgecolors="#a3a29b", linewidths=0.35, zorder=4, label="TPFLSH minor: phase label only, loose (1e-5..1e-3), or knife-edge"),
    "rp_scope": dict(marker="v", s=7, facecolors="none", edgecolors="#8a4fb3", linewidths=0.35, alpha=0.75, zorder=5, label="TPFLSH misses a verified split (LLE, water condensation)"),
}


def load_verdicts(files):
    v = {}
    for fn in files:
        for r in csv.DictReader(open(fn)):
            v[(r["mixture"], r["i"])] = r["verdict"]
    return v


def main(out, args):
    verdicts = load_verdicts([a[2:] for a in args if a.startswith("v:")])
    build = next((a[6:] for a in args if a.startswith("build:")), None)
    if not build:
        sys.exit("build:<description> is required - the figure must say which CoolProp build it shows")
    args = [a for a in args if not a.startswith("v:") and not a.startswith("build:")]
    # Rows: groups of files separated by '--'; without '--', 3 panels per row.
    groups, cur = [], []
    for a in args:
        if a == "--":
            groups.append(cur)
            cur = []
        else:
            cur.append(a)
    groups.append(cur)
    files = [f for g in groups for f in g]
    mixes = load(files)
    n = len(mixes)
    nstates = max(len(v) for v in mixes.values())
    if len(groups) == 1:  # default: 3 per row
        names = list(mixes)
        rows = [names[i:i + 3] for i in range(0, n, 3)]
    else:
        rows = [list(load(g)) for g in groups]
    ncol = max(len(r) for r in rows)
    grid = lcm(*[len(r) for r in rows]) if all(rows) else ncol
    nrow = len(rows)
    fig = plt.figure(figsize=(4.6 * ncol, 4.1 * nrow + 0.9))
    extra = 0.055 if verdicts else 0.0  # room for the marker legend under the colorbar
    gs = fig.add_gridspec(nrow, grid, hspace=0.62, wspace=0.9 if grid > ncol else 0.18, bottom=(0.17 if nrow <= 2 else 0.1) + extra, top=0.86 if nrow <= 2 else 0.91)
    placement = {}
    for r, names in enumerate(rows):
        w = grid // len(names)
        for c, name in enumerate(names):
            placement[name] = (r, c * w, (c + 1) * w, c == 0)
    fig.patch.set_facecolor(SURFACE)
    norm = Normalize(-LIM, LIM)
    for k, (name, rows) in enumerate(mixes.items()):
        r0, c0, c1, first_in_row = placement[name]
        ax = fig.add_subplot(gs[r0, c0:c1])
        ax.set_facecolor(SURFACE)
        pts, bad_cp, bad_rp, sum_cp, sum_rp = [], 0, 0, 0.0, 0.0
        for r in rows:
            tcp, trp = float(r["t_cp_us"]), float(r["t_rp_us"])
            fcp, frp = int(r["cp_fail"]) != 0, int(r["ierr_rp"]) > 0
            bad_cp += fcp
            bad_rp += frp
            if fcp or frp or not (tcp > 0 and trp > 0):
                continue
            sum_cp += tcp
            sum_rp += trp
            Q = float(r["Q_cp"])
            pts.append((float(r["T"]), float(r["p"]) / 1e6, math.log10(tcp / trp), 0 < Q < 1, verdicts.get((name, r["i"]))))
        # marker size shrinks with density so a 10k-state map still reads as points
        sc = math.sqrt(2000.0 / max(len(rows), 1))
        # single phase first, two-phase on top
        for two in (False, True):
            sel = [q for q in pts if q[3] == two]
            if not sel:
                continue
            ax.scatter([q[0] for q in sel], [q[1] for q in sel], c=[max(-LIM, min(LIM, q[2])) for q in sel], cmap=CMAP, norm=norm,
                       s=(11 if two else 9) * sc, linewidths=(0.6 if two else 0.25) * math.sqrt(sc), edgecolors=INK if two else "#b8b7b0", zorder=3 if two else 2)
        for key, style in VMARK.items():
            sel = [q for q in pts if q[4] == key]
            if sel:
                st = {k: v for k, v in style.items() if k != "label"}
                ax.scatter([q[0] for q in sel], [q[1] for q in sel], **st)
        lr = sorted(q[2] for q in pts)
        med = 10 ** lr[len(lr) // 2]
        faster = sum(1 for v in lr if v < 0) / len(lr)
        def fs(x):
            return f"{1 / x:.1f}× faster" if x < 1 else f"{x:.1f}× slower"
        tot = sum_cp / sum_rp
        ax.set_title(f"{name}\nmedian state {fs(med)} · total time {fs(tot)}\nfaster at {100 * faster:.0f}% of states", fontsize=9, color=INK, loc="left")
        ax.set_yscale("log")
        ax.set_ylim(0.01, 30)
        ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
        for s in ax.spines.values():
            s.set_color("#d4d3cc")
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.set_xlabel("T / K", fontsize=9, color=INK)
        if first_in_row:
            ax.set_ylabel("p / MPa", fontsize=9, color=INK)
        if bad_cp or bad_rp:
            ax.text(0.98, 0.03, f"not plotted (error raised): TPFLSH ierr>0 {bad_rp}, CoolProp threw {bad_cp}", transform=ax.transAxes, ha="right",
                    fontsize=7, color=MUTED, zorder=1, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5, alpha=0.8))  # under the verdict markers
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    cax = fig.add_axes([0.2, (0.06 if nrow <= 2 else 0.035) + extra, 0.6, 0.022 if nrow <= 2 else 0.012])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    ticks = [-math.log10(30), -1, -math.log10(3), 0, math.log10(3), 1, math.log10(30)]
    cb.set_ticks(ticks)
    cb.set_ticklabels(["30× faster", "10× faster", "3× faster", "equal", "3× slower", "10× slower", "30× slower"])
    cb.ax.tick_params(labelsize=8, colors=INK)
    cb.outline.set_edgecolor("#d4d3cc")
    cb.set_label("CoolProp PT flash time / REFPROP 10 TPFLSH time", fontsize=9, color=INK)
    if verdicts:
        from matplotlib.lines import Line2D
        handles = []
        vcount = {}
        for vv in verdicts.values():
            vcount[vv] = vcount.get(vv, 0) + 1
        for key, st in VMARK.items():
            if not vcount.get(key):
                continue
            face = st.get("c", st.get("facecolors"))
            handles.append(Line2D([], [], linestyle="", marker=st["marker"], markersize=7 if key.startswith("cp") else 5.5,
                                  markerfacecolor=face if face != "none" else "none", markeredgecolor=st["edgecolors"] if key != "cp_wrong" else INK,
                                  markeredgewidth=1.0, label=f"{st['label']} ({vcount[key]})"))
        fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.0),
                   labelcolor=INK)
    fig.suptitle(f"CoolProp PT flash vs REFPROP 10 TPFLSH, per state ({nstates} states per mixture, min of 3 timings each{', mixtures run as parallel processes' if nstates > 2000 else ''})\n"
                 f"CoolProp build: {build}.  Dark outline: CoolProp publishes a two-phase state.\n"
                 "GERG-2008 on both sides, except R454B (CoolProp HEOS vs REFPROP default mixture model; its disagreements are not judged).",
                 fontsize=10, color=INK, x=0.06, ha="left", y=0.975 if nrow <= 2 else 0.985)
    fig.savefig(out, dpi=140, bbox_inches="tight", facecolor=SURFACE)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
