#!/usr/bin/env python3
"""SBTL/BICUBIC accuracy maps in (T, P) coordinates for any pure fluid.

Probes density relative error vs HEOS over a (T, P) grid using
PT_INPUTS for both BICUBIC and SBTL.  SBTL routes pure-fluid PT through
the coordinate-aligned NormalizedPTTable (LIQUID / VAPOR / SUPER) plus
a small HEOS-fallback box around (T_c, p_c); BICUBIC uses Hermite-bicubic
on the legacy logpT table.

Three windows: A) around critical, B) subcritical, C) supercritical.
Solid-region points (T < T_melt) are filtered.

Usage:  python3 sbtl_acc_fig_pt.py [FluidName] [/path/output.png]

Note: 200×200 default × 3 backends × ~10 ms/HEOS ≈ 25 min/fluid.  Lower
NT/NP in the probe() call for faster turnaround.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import CoolProp.CoolProp as CP

ref = sys.argv[1] if len(sys.argv) > 1 else 'CO2'
output_path = sys.argv[2] if len(sys.argv) > 2 else f'/tmp/sbtl_acc_fig_pt_{ref.lower()}.png'

heos = CP.AbstractState('HEOS', ref)
T_triple = heos.Ttriple()
T_crit = heos.T_critical()
p_crit = heos.p_critical()
T_min_eos = max(T_triple, heos.Tmin())
T_max_eos = heos.Tmax()
p_max_eos = heos.pmax()
heos.update(CP.QT_INPUTS, 0.0, T_triple * 1.0001)
p_triple = heos.p()

print(f'Fluid: {ref}', flush=True)
print(f'  Tc = {T_crit:.3f} K   pc = {p_crit / 1e5:.4f} bar', flush=True)

# Auto-windows
def T_range(lo_frac, hi_frac):
    return (max(T_triple * 1.001, lo_frac * T_crit), hi_frac * T_crit)


# Window A: tight box around the critical point
A_T = (T_crit * 0.95, T_crit * 1.05)
A_p = (0.95 * p_crit, 1.36 * p_crit)

# Window B: subcritical.  Stay clear of the triple point — SBTL's
# logpT cell selection can hang for queries right on the triple-line
# corner, where every neighboring cell straddles either the sat dome
# or the EOS validity boundary.
B_T = (max(T_triple * 1.05, 0.5 * T_crit), 0.99 * T_crit)
B_p = (max(p_triple * 1.5, p_crit * 1e-3), 0.99 * p_crit)

# Window C: supercritical (and supercritical-T at high P)
C_T = (1.05 * T_crit, min(T_max_eos * 0.7, 3 * T_crit))
C_p = (1.02 * p_crit, min(p_max_eos, 4.0 * p_crit))

# Backends
prime_T = max(T_triple * 1.5, 1.5 * T_min_eos)
prime_p = max(p_triple * 2, 1e5)
print('Building BICUBIC...', flush=True)
bicu = CP.AbstractState('BICUBIC&HEOS', ref); bicu.update(CP.PT_INPUTS, prime_p, prime_T)
print('Building SBTL...', flush=True)
sbtl = CP.AbstractState('SBTL&HEOS', ref); sbtl.update(CP.PT_INPUTS, prime_p, prime_T)


def probe(BCK, T_lim, p_lim, NT=200, NP=200):
    """ρ relative error vs HEOS on a (T, p) grid via PT_INPUTS."""
    has_melt = heos.has_melting_line()
    Tg = np.linspace(*T_lim, NT)
    Pg = np.geomspace(*p_lim, NP)
    TT, PP = np.meshgrid(Tg, Pg)
    err = np.full_like(TT, np.nan)
    for i in range(NP):
        for j in range(NT):
            T = TT[i, j]; P = PP[i, j]
            try:
                heos.update(CP.PT_INPUTS, P, T)
                if 0.0 < heos.Q() < 1.0: continue
                rho = heos.rhomass()
            except Exception:
                continue
            if T < T_triple:
                continue
            if has_melt:
                try:
                    if T < heos.melting_line(CP.iT, CP.iP, P): continue
                except Exception:
                    pass
            try:
                BCK.update(CP.PT_INPUTS, P, T)
                err[i, j] = abs(BCK.rhomass() - rho) / rho
            except Exception:
                pass
    return TT, PP, err


def saturation_curve_pt():
    Ts = np.linspace(T_triple + 0.01, T_crit - 0.001, 400)
    p_arr = []
    for t in Ts:
        try:
            heos.update(CP.QT_INPUTS, 0.0, t); p_arr.append(heos.p())
        except Exception:
            p_arr.append(np.nan)
    return Ts, np.array(p_arr)


windows = [
    ('A) Critical region', A_T, A_p),
    ('B) Subcritical',     B_T, B_p),
    ('C) Supercritical',   C_T, C_p),
]

cNorm = colors.LogNorm(vmin=1e-7, vmax=1e0)
fig, axes = plt.subplots(3, 2, figsize=(12, 14))
Ts_sat, p_sat_arr = saturation_curve_pt()

for row_i, (w_label, T_lim, p_lim) in enumerate(windows):
    label_suffix = (
        f'(T={T_lim[0]:.0f}-{T_lim[1]:.0f} K, '
        f'P={p_lim[0] / 1e5:.2f}-{p_lim[1] / 1e5:.0f} bar)'
    )
    for col_i, (BCK, label) in enumerate([(bicu, 'BICUBIC PT'), (sbtl, 'SBTL PT')]):
        ax = axes[row_i, col_i]
        TT, PP, err = probe(BCK, T_lim, p_lim)
        im = ax.pcolormesh(TT, PP / 1e5, np.maximum(err, 1e-12),
                           cmap='plasma', shading='auto', norm=cNorm)
        plt.colorbar(im, ax=ax, label=r'$|\rho_{\mathrm{tab}} / \rho_{\mathrm{HEOS}} - 1|$',
                     shrink=0.85)
        if Ts_sat.size:
            mask = (Ts_sat >= T_lim[0]) & (Ts_sat <= T_lim[1])
            ax.plot(Ts_sat[mask], p_sat_arr[mask] / 1e5, 'w-', lw=1.5, alpha=0.85)
        ax.axhline(p_crit / 1e5, color='w', ls=':', lw=1, alpha=0.6)
        ax.axvline(T_crit, color='w', ls=':', lw=1, alpha=0.6)
        ax.set_xlim(*T_lim)
        ax.set_ylim(p_lim[0] / 1e5, p_lim[1] / 1e5)
        ax.set_yscale('log')
        med = np.nanmedian(err); p99 = np.nanpercentile(err, 99); mx = np.nanmax(err)
        ax.set_title(
            f'{w_label} {label_suffix}: {label}\n'
            f'med={med:.1e}  99th={p99:.1e}  max={mx:.1e}',
            fontsize=10)
        if col_i == 0:
            ax.set_ylabel('P / bar')
        if row_i == 2:
            ax.set_xlabel('T / K')

plt.suptitle(
    f'{ref}  ρ relative error vs HEOS via PT_INPUTS  (default 200×200 tables)\n'
    f'Saturation curve in white; T_crit and p_crit dotted',
    fontsize=12, y=1.0)
plt.tight_layout()
plt.savefig(output_path, dpi=110, bbox_inches='tight')
print(f'Saved {output_path}')
