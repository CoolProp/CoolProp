#!/usr/bin/env python3
"""Failure maps from bench_reliab.cpp's CSV: where each solver fails, or returns a root that is not
the stable one, in (T, p).  One figure per model; rows = mixtures, columns = solvers.
Our solver has no failures, so it has no column.

    python3 plot_reliab.py reliab.csv figs/
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# validated categorical pair (dataviz reference palette, slots 1-2; light surface #fcfcfb)
C_FAIL, C_NONSTABLE, SURFACE, INK, MUTED = "#eb6834", "#2a78d6", "#fcfcfb", "#1f1f1e", "#6b6a64"
SOLVERS = [("CoolProp", "CoolProp solver_rho_Tp"), ("REFPROP_kph2", "REFPROP TPRHO kph=2"), ("REFPROP_kph1", "REFPROP TPRHO kph=1")]


def main(csvname, outdir):
    rows = defaultdict(list)
    for r in csv.DictReader(open(csvname)):
        rows[(r["model"], r["mixture"])].append(r)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    for model in ("HEOS", "GERG2008"):
        mixes = [m for (mo, m) in rows if mo == model and not m.startswith("humid")]
        mixes = sorted(set(mixes), key=lambda s: [k for (mo, k) in rows].index(s))
        fig, axes = plt.subplots(len(mixes), 3, figsize=(11, 2.3 * len(mixes)), sharex=False, sharey=True, squeeze=False)
        fig.patch.set_facecolor(SURFACE)
        for i, mix in enumerate(mixes):
            data = rows[(model, mix)]
            for j, (key, title) in enumerate(SOLVERS):
                ax = axes[i][j]
                ax.set_facecolor(SURFACE)
                d = [r for r in data if r["solver"] == key]
                for outcome, color, marker, face, label in [("nonstable", C_NONSTABLE, "o", "none", "returned a non-stable root"),
                                                            ("fail", C_FAIL, "x", C_FAIL, "failed (no result)"),
                                                            ("notroot", C_FAIL, "^", "none", "not a root of this model")]:
                    pts = [(float(r["T"]), float(r["p"]) / 1e6) for r in d if r["outcome"] == outcome]
                    if pts:
                        ax.scatter([a for a, _ in pts], [b for _, b in pts], s=9, marker=marker, facecolors=face if marker != "x" else None,
                                   edgecolors=color if marker != "x" else None, c=color if marker == "x" else None, linewidths=0.8, label=label)
                ax.set_yscale("log")
                ax.set_ylim(0.01, 30)
                ax.grid(True, color="#e6e5df", linewidth=0.6)
                for s in ax.spines.values():
                    s.set_color("#d4d3cc")
                ax.tick_params(colors=MUTED, labelsize=8)
                if i == 0:
                    ax.set_title(title, fontsize=10, color=INK)
                if j == 0:
                    ax.set_ylabel(f"{mix}\np / MPa", fontsize=9, color=INK)
                if i == len(mixes) - 1:
                    ax.set_xlabel("T / K", fontsize=9, color=INK)
                n_f = sum(r["outcome"] == "fail" for r in d)
                n_n = sum(r["outcome"] == "nonstable" for r in d)
                ax.text(0.02, 0.04, f"{n_f} failed · {n_n} non-stable", transform=ax.transAxes, fontsize=7.5, color=MUTED)
        handles = [plt.Line2D([], [], marker="x", color=C_FAIL, linestyle="", label="failed (no result)"),
                   plt.Line2D([], [], marker="o", markerfacecolor="none", markeredgecolor=C_NONSTABLE, linestyle="", label="returned a non-stable root"),
                   plt.Line2D([], [], marker="^", markerfacecolor="none", markeredgecolor=C_FAIL, linestyle="", label="not a root of this model")]
        fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 1.0))
        mname = "reference EOS (CoolProp HEOS; REFPROP default)" if model == "HEOS" else "GERG-2008 (CoolProp GERG2008; REFPROP GERG mode)"
        fig.suptitle(f"Where density solvers fail or pick a non-stable root — {mname}\n5000 states per mixture; the Chebyshev solver had no failures",
                     fontsize=11, color=INK, y=1.04)
        fig.tight_layout()
        f = out / f"reliab_map_{model}.png"
        fig.savefig(f, dpi=130, bbox_inches="tight", facecolor=SURFACE)
        print("wrote", f)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
