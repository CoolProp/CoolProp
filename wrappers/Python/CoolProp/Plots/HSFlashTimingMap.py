# -*- coding: utf-8 -*-
"""
Interactive map of the single-phase H,S-flash cascade cost over the full fluid
domain (water by default), in p-T and/or h-s coordinates.

The data is produced by the C++ Catch2 test ``[HS_watermap]`` (see
``src/Tests/CoolProp-Tests-HS.cpp``), which sweeps a dense (log p, linear T) grid
-- single-phase points only, exactly as the ConsistencyPlot suite samples -- and
writes one CSV row per point with the EOS-evaluation count, wall time, and the
cascade leg that solved it:

    ./build_rel/CatchTestRunner "[HS_watermap]"        # writes water_hs_map.csv

Then render the interactive HTML (hover shows every field; pan/zoom/box-select):

    python -m CoolProp.Plots.HSFlashTimingMap water_hs_map.csv             # both coords
    python -m CoolProp.Plots.HSFlashTimingMap water_hs_map.csv --coords hs

Each coordinate system gets a 3-panel figure (EOS evals, solve time, winning
cascade leg) with the phase boundary overlaid -- the p-T saturation curve, or the
h-s two-phase dome.  Requires plotly for the interactive HTML; falls back to a
static matplotlib figure if plotly is unavailable.
"""
from __future__ import print_function, division, absolute_import
import os
import argparse
import warnings

import numpy as np
import pandas


# Cascade leg labels (must match solve_HS_cascade in the C++ source).
METHOD_LABEL = {0: "unsolved", 1: "saturation-anchor", 2: "supercritical-isentrope", 3: "ideal-gas-departure"}
METHOD_COLOR = {0: "#000000", 1: "#1f77b4", 2: "#d62728", 3: "#2ca02c"}

# Coordinate systems: (df x-column, df y-column, x-label, y-label).
COORDS = {
    "pt": ("T_K", "log10p", "T [K]", "log10(p [Pa])"),
    "hs": ("s_JmolK", "h_Jmol", "s [J/mol/K]", "h [J/mol]"),
}


def load(csv_path):
    df = pandas.read_csv(csv_path)
    df["log10p"] = np.log10(df["p_Pa"])
    df["method_label"] = df["method"].map(METHOD_LABEL)
    return df


def phase_boundary(fluid, n=240):
    """Phase boundary in both coordinate systems, computed from CoolProp.

    Returns a dict with the p-T saturation curve and the h-s two-phase dome
    (liquid branch up + vapor branch back down, a closed outline) plus the
    critical point in each system.  Returns {} if CoolProp is unavailable."""
    try:
        import CoolProp.CoolProp as CP
    except ImportError:
        return {}
    try:
        Tt = CP.PropsSI("Ttriple", fluid)
        Tc = CP.PropsSI("Tcrit", fluid)
    except Exception as exc:
        # A bad --fluid (or any CoolProp lookup failure) should not be swallowed:
        # warn with the fluid name and underlying error so the omitted overlay is
        # explained rather than silently dropped, then skip the overlay.
        warnings.warn("phase_boundary: cannot resolve triple/critical point for fluid {0!r}: {1}".format(fluid, exc))
        return {}
    T, p, hL, sL, hV, sV = [], [], [], [], [], []
    for Ti in np.linspace(Tt, Tc, n):
        Ti = float(Ti)
        try:
            row = (CP.PropsSI("P", "T", Ti, "Q", 0, fluid),
                   CP.PropsSI("Hmolar", "T", Ti, "Q", 0, fluid), CP.PropsSI("Smolar", "T", Ti, "Q", 0, fluid),
                   CP.PropsSI("Hmolar", "T", Ti, "Q", 1, fluid), CP.PropsSI("Smolar", "T", Ti, "Q", 1, fluid))
        except Exception:
            continue
        if not all(np.isfinite(v) for v in row):
            continue
        T.append(Ti)
        p.append(row[0]); hL.append(row[1]); sL.append(row[2]); hV.append(row[3]); sV.append(row[4])
    if not T:
        return {}
    T, p = np.array(T), np.array(p)
    hL, sL, hV, sV = map(np.array, (hL, sL, hV, sV))
    # Closed h-s dome: up the liquid branch, back down the vapor branch.
    dome_s = np.concatenate([sL, sV[::-1]])
    dome_h = np.concatenate([hL, hV[::-1]])
    return {
        "pt": dict(lines=[(T, np.log10(p), "saturation")], crit=(T[-1], np.log10(p[-1]))),
        "hs": dict(lines=[(dome_s, dome_h, "two-phase dome")], crit=(0.5 * (sL[-1] + sV[-1]), 0.5 * (hL[-1] + hV[-1]))),
    }


def log_color(values):
    """Map positive values to log10 for color, plus 1-2-5 decade colorbar ticks
    labelled with the ORIGINAL values (so the bulk variation is visible instead of
    being washed out by a few high-cost outliers)."""
    v = np.asarray(values, dtype=float)
    vp = v[v > 0]
    vmin = vp.min() if len(vp) else 1.0
    c = np.log10(np.clip(v, vmin, None))
    ticks = [m * 10.0 ** e for e in range(int(np.floor(np.log10(vmin))), int(np.ceil(np.log10(v.max()))) + 1)
             for m in (1, 2, 5) if vmin <= m * 10.0 ** e <= v.max()]
    return c, np.log10(ticks), ["{0:g}".format(t) for t in ticks]


def make_plotly(df, title, html_path, xcol, ycol, xlabel, ylabel, overlay):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    hover = (
        "T = %{customdata[0]:.3f} K<br>"
        "p = %{customdata[1]:.4g} Pa<br>"
        "rho = %{customdata[2]:.4g} mol/dm^3<br>"
        "h = %{customdata[6]:.4g} J/mol<br>"
        "s = %{customdata[7]:.4g} J/mol/K<br>"
        "evals = %{customdata[3]:d}<br>"
        "time = %{customdata[4]:.2f} us<br>"
        "leg = %{customdata[5]}<extra></extra>"
    )

    def cd(d):
        return np.stack([d["T_K"], d["p_Pa"], d["rho_moldm3"], d["evals"], d["microseconds"], d["method_label"],
                         d["h_Jmol"], d["s_JmolK"]], axis=-1)

    fig = make_subplots(rows=1, cols=3, shared_yaxes=True,
                        subplot_titles=("EOS evaluations (log)", "Solve time [us] (log)", "Cascade leg"), horizontal_spacing=0.06)
    ev_c, ev_tv, ev_tt = log_color(df["evals"])
    us_c, us_tv, us_tt = log_color(df["microseconds"])
    fig.add_trace(go.Scattergl(x=df[xcol], y=df[ycol], mode="markers",
                  marker=dict(size=4, color=ev_c, colorscale="Viridis",
                              colorbar=dict(title="evals", x=0.28, len=0.9, tickvals=ev_tv, ticktext=ev_tt)),
                  customdata=cd(df), hovertemplate=hover, name="evals"), row=1, col=1)
    fig.add_trace(go.Scattergl(x=df[xcol], y=df[ycol], mode="markers",
                  marker=dict(size=4, color=us_c, colorscale="Inferno",
                              colorbar=dict(title="us", x=0.63, len=0.9, tickvals=us_tv, ticktext=us_tt)),
                  customdata=cd(df), hovertemplate=hover, name="time"), row=1, col=2)
    for m, label in METHOD_LABEL.items():
        sub = df[df["method"] == m]
        if len(sub):
            fig.add_trace(go.Scattergl(x=sub[xcol], y=sub[ycol], mode="markers", marker=dict(size=4, color=METHOD_COLOR[m]),
                          customdata=cd(sub), hovertemplate=hover, name=label), row=1, col=3)

    # Phase boundary overlaid on every panel.  Scattergl (same WebGL layer as the
    # markers, added last -> on top) drawn as a black halo under a yellow line so
    # it stays visible over Viridis, Inferno, AND the white categorical panel.
    if overlay:
        for c in (1, 2, 3):
            for (lx, ly, lname) in overlay["lines"]:
                fig.add_trace(go.Scattergl(x=lx, y=ly, mode="lines", line=dict(color="black", width=6),
                              legendgroup="bound", showlegend=False, hoverinfo="skip"), row=1, col=c)
                fig.add_trace(go.Scattergl(x=lx, y=ly, mode="lines", line=dict(color="yellow", width=2.5),
                              name=lname, legendgroup="bound", showlegend=(c == 1),
                              hovertemplate=lname + ": x=%{x:.4g}, y=%{y:.4g}<extra></extra>"), row=1, col=c)
            cx, cy = overlay["crit"]
            fig.add_trace(go.Scattergl(x=[cx], y=[cy], mode="markers",
                          marker=dict(color="yellow", size=11, symbol="x", line=dict(color="black", width=1.5)),
                          name="critical point", legendgroup="cp", showlegend=(c == 1),
                          hovertemplate="critical point<extra></extra>"), row=1, col=c)

    for c in (1, 2, 3):
        fig.update_xaxes(title_text=xlabel, row=1, col=c)
    fig.update_yaxes(title_text=ylabel, row=1, col=1)
    fig.update_layout(title=title, width=1500, height=620, legend=dict(x=1.0, y=1.0))
    fig.write_html(html_path)
    print("Wrote interactive plot:", html_path)


def make_matplotlib(df, title, png_path, xcol, ycol, xlabel, ylabel, overlay):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    sc0 = axes[0].scatter(df[xcol], df[ycol], c=df["evals"].clip(lower=1), s=6, cmap="viridis", norm=LogNorm())
    fig.colorbar(sc0, ax=axes[0], label="EOS evals (log)")
    axes[0].set_title("EOS evaluations")
    sc1 = axes[1].scatter(df[xcol], df[ycol], c=df["microseconds"].clip(lower=1e-3), s=6, cmap="inferno", norm=LogNorm())
    fig.colorbar(sc1, ax=axes[1], label="solve time [us] (log)")
    axes[1].set_title("Solve time")
    for m, label in METHOD_LABEL.items():
        sub = df[df["method"] == m]
        if len(sub):
            axes[2].scatter(sub[xcol], sub[ycol], s=6, c=METHOD_COLOR[m], label=label)
    axes[2].legend(loc="best", fontsize=8)
    axes[2].set_title("Cascade leg")
    if overlay:
        for ax in axes:
            for (lx, ly, _) in overlay["lines"]:
                ax.plot(lx, ly, "k-", lw=3)
                ax.plot(lx, ly, "y-", lw=1.4)
            ax.plot(*overlay["crit"], "x", color="yellow", mec="black", ms=9)
    for ax in axes:
        ax.set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(png_path, dpi=130)
    print("Wrote static plot:", png_path)


def render(df, title, base, coord, overlay):
    xcol, ycol, xlabel, ylabel = COORDS[coord]
    try:
        out = base + "_" + coord + ".html"
        make_plotly(df, title + "  [" + coord + "]", out, xcol, ycol, xlabel, ylabel, overlay)
    except ImportError:
        out = base + "_" + coord + ".png"
        print("plotly not available; falling back to matplotlib ->", out)
        make_matplotlib(df, title + "  [" + coord + "]", out, xcol, ycol, xlabel, ylabel, overlay)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="?", default="water_hs_map.csv", help="CSV produced by the [HS_watermap] test")
    parser.add_argument("-o", "--output", default=None, help="output basename (coord + extension appended)")
    parser.add_argument("--coords", choices=["pt", "hs", "both"], default="both", help="coordinate system(s) to render")
    parser.add_argument("--fluid", default="Water", help="fluid name (for the title and phase boundary)")
    args = parser.parse_args(argv)

    if not os.path.exists(args.csv):
        parser.error("CSV not found: {0} (run the [HS_watermap] Catch2 test first)".format(args.csv))
    df = load(args.csv)
    n_unsolved = int((df["ok"] == 0).sum())
    title = "{0}: single-phase H,S cascade -- {1} pts, {2} unsolved, mean {3:.1f} evals / {4:.2f} us".format(
        args.fluid, len(df), n_unsolved, df["evals"].mean(), df["microseconds"].mean())
    boundary = phase_boundary(args.fluid)

    base = os.path.splitext(args.output or args.csv)[0]
    for coord in (["pt", "hs"] if args.coords == "both" else [args.coords]):
        render(df, title, base, coord, boundary.get(coord))


if __name__ == "__main__":
    main()
