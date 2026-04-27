#!/usr/bin/env python3
"""Regenerate the Tauri icon set from the canonical CoolProp logo.

Renders the water phase-diagram logo (orange melt + saturation curves over a
blue density colormap) at 1024x1024, then runs `tauri icon` to produce the
full per-platform icon set under src-tauri/icons/.

Run from wrappers/GUI/. Requires CoolProp Python bindings, matplotlib, scipy,
Pillow, and a working `npx tauri` (i.e. `npm install` already done).
"""
import os
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

import CoolProp as CP


def render_logo(out_path: str, size_px: int = 1024) -> None:
    Water = CP.AbstractState("HEOS", "Water")
    Tc = Water.keyed_output(CP.iT_critical)
    T_min, T_max = 200, 1000
    p_max = Water.keyed_output(CP.iP_max)
    p_triple = 611.657
    T_triple = 273.16

    steps = 2000
    PP = list(np.logspace(np.log10(p_triple), np.log10(p_max), steps))
    TT = [Water.melting_line(CP.iT, CP.iP, p) for p in PP]

    for T in np.linspace(max(TT), 355, int(steps / 10)):
        TT.append(T)
        theta = T / 273.31
        PP.append((1 - 1.07476 * (1 - theta**4.6)) * 632.4e6)

    for T in np.linspace(355, 715, int(steps / 10)):
        TT.append(T)
        theta = T / 355
        lnpi = (
            0.173683e1 * (1 - 1 / theta)
            - 0.544606e-1 * (1 - theta**5)
            + 0.806106e-7 * (1 - theta**22)
        )
        PP.append(np.exp(lnpi) * 2216e6)

    steps_small = int(steps / 10)
    p_melt = np.logspace(np.log10(np.min(PP)), np.log10(np.max(PP)), steps_small)
    T_melt = scipy.interpolate.interp1d(np.log10(PP), TT)(np.log10(p_melt))

    T_sat = np.linspace(T_triple, Tc, len(T_melt))
    p_sat = CP.CoolProp.PropsSI("P", "T", T_sat, "Q", [0] * len(T_sat), "Water")

    Tg, Dg, Pg = [], [], []
    for T in np.linspace(T_min, T_max, steps_small):
        for p in np.logspace(np.log10(np.min(p_melt)), np.log10(np.max(p_melt)), steps_small):
            Tm = scipy.interpolate.interp1d(p_melt, T_melt)(p)
            if T < Tm:
                continue
            D = CP.CoolProp.PropsSI("D", "T", T, "P", min(p, p_max), "Water")
            Tg.append(T)
            Dg.append(np.log10(D))
            Pg.append(p)

    melt_args = dict(color="orange", lw=6, solid_capstyle="round")
    nm = matplotlib.colors.Normalize(min(Dg), max(Dg))
    rho_args = dict(cmap=plt.get_cmap("Blues"), norm=nm)

    inches = size_px / 100.0
    fig = plt.figure(figsize=(inches, inches), dpi=100)
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.plot(T_melt, p_melt, **melt_args)
    ax.plot(T_sat, p_sat, **melt_args)
    ax.scatter(Tg, Pg, c=Dg, edgecolor="none", s=60, **rho_args)
    dx = np.min(T_melt) * 0.01
    ax.set_xlim([np.min(T_melt) - dx, np.max(T_melt) + dx])
    ax.set_ylim([np.min(p_melt) * 0.875, np.max(p_melt) * 1.125])
    ax.set_yscale("log")
    ax.axis("off")
    fig.savefig(out_path, transparent=True, dpi=100)
    plt.close(fig)


def main() -> int:
    here = os.path.abspath(os.path.dirname(__file__))
    gui_dir = os.path.dirname(here)
    if not os.path.isdir(os.path.join(gui_dir, "src-tauri")):
        print(f"error: expected src-tauri/ next to scripts/, got {gui_dir}", file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, "CoolPropLogo_1024.png")
        render_logo(src, size_px=1024)
        subprocess.run(["npx", "tauri", "icon", src], cwd=gui_dir, check=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
