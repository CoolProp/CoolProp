#!/usr/bin/env python3
"""Diagram of the Kunick SBTL (v,u) region scheme for any pure fluid.

Mirrors IAPWS G13-15 Fig. 9 for water and generalizes to other fluids:
  - Region L  (liquid)               : sat-L  →  p_max isobar, bounded above by u_c
  - Region G  (gas + supercritical-low-D) : sat-V  →  p_triple, bounded below by u_c, above by T_HT
  - Region HT (high-temperature)     : T > T_HT  (1073.15 K for water)
  - Region TP (two-phase dome)       : x ∈ (0, 1)

The u_c = u(D_crit, T_crit) iso-energy line splits the supercritical
region into L (high D) and G (low D) — Kunick's choice instead of a
constant-D split.

Usage:  python3 sbtl_du_diagram.py [FluidName] [/path/output.png]
"""
import sys

import numpy as np
import matplotlib.pyplot as plt

import CoolProp.CoolProp as CP

ref = sys.argv[1] if len(sys.argv) > 1 else 'Water'
output = sys.argv[2] if len(sys.argv) > 2 else f'/tmp/sbtl_du_diagram_{ref.lower()}.png'

heos = CP.AbstractState('HEOS', ref)
T_triple = max(heos.Ttriple(), heos.Tmin())
T_crit = heos.T_critical()
p_crit = heos.p_critical()
D_crit = heos.rhomolar_critical()
T_max_eos = heos.Tmax()
p_max_eos = heos.pmax()

# u_c = u at the critical point (Kunick's L/G boundary for supercritical)
heos.update(CP.DmolarT_INPUTS, D_crit, T_crit)
u_c_molar = heos.umolar()
u_c_mass = heos.umass()

# Triple-point sat properties (Kunick's lower-T boundary of L and G)
heos.update(CP.QT_INPUTS, 0, T_triple)
v_sat_L_triple = 1.0 / heos.rhomass()
u_sat_L_triple = heos.umass()
heos.update(CP.QT_INPUTS, 1, T_triple)
v_sat_V_triple = 1.0 / heos.rhomass()
u_sat_V_triple = heos.umass()
p_sat_triple = heos.p()

# Kunick uses T_HT = 1073.15 K for water (= T_crit * 1.658).  For other
# fluids, scale proportionally — a fixed isotherm well above T_crit.
T_HT = T_crit * 1073.15 / 647.096

# p_max for the L region — Kunick uses 100 MPa for water.  Generalize as
# min(p_max_eos, ~5 * p_crit), capped at the EOS limit.
p_max_L = min(p_max_eos, 100e6 if ref == 'Water' else 5 * p_crit)

# p_min — Kunick uses 611.212 Pa for water (the triple-point pressure).
p_min = p_sat_triple

print(f'Fluid: {ref}', flush=True)
print(f'  Tc={T_crit:.2f}K  pc={p_crit / 1e6:.3f}MPa  Dc={D_crit:.2f}mol/m³  '
      f'u_c={u_c_mass / 1e3:.2f} kJ/kg', flush=True)
print(f'  Kunick boundaries:  T_HT={T_HT:.1f} K  p_max_L={p_max_L / 1e6:.1f} MPa  '
      f'p_min={p_min:.3f} Pa', flush=True)


# Sat dome
def sat_curve():
    Ts = np.linspace(T_triple + 1e-3, T_crit - 1e-4, 400)
    vL, uL, vV, uV = [], [], [], []
    for T in Ts:
        try:
            heos.update(CP.QT_INPUTS, 0, T); vL.append(1.0 / heos.rhomass()); uL.append(heos.umass())
            heos.update(CP.QT_INPUTS, 1, T); vV.append(1.0 / heos.rhomass()); uV.append(heos.umass())
        except Exception:
            # query landed outside table envelope; skip this point
            pass
    return np.array(vL), np.array(uL), np.array(vV), np.array(uV)


def isobar(P, T_lo, T_hi, N=200):
    Ts = np.linspace(T_lo, T_hi, N)
    vs, us = [], []
    for T in Ts:
        try:
            heos.update(CP.PT_INPUTS, P, T)
            if 0 < heos.Q() < 1: continue
            vs.append(1.0 / heos.rhomass()); us.append(heos.umass())
        except Exception:
            # query landed outside table envelope; skip this point
            pass
    return np.array(vs), np.array(us)


def isotherm(T, p_lo, p_hi, N=200):
    Ps = np.geomspace(max(p_lo, 1e-4), p_hi, N)
    vs, us = [], []
    for P in Ps:
        try:
            heos.update(CP.PT_INPUTS, P, T)
            if 0 < heos.Q() < 1: continue
            vs.append(1.0 / heos.rhomass()); us.append(heos.umass())
        except Exception:
            # query landed outside table envelope; skip this point
            pass
    return np.array(vs), np.array(us)


vL, uL, vV, uV = sat_curve()

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# ---------------------------------------------------------------------------
# Panel 1: (v, u) — Kunick's native coords, log v
# ---------------------------------------------------------------------------
ax = axes[0]
ax.set_xscale('log')

# Region boundaries
v_p_max, u_p_max = isobar(p_max_L, T_triple, T_HT)              # p=p_max_L isobar
v_p_min, u_p_min = isobar(p_min, T_triple, T_HT)                # p=p_min isobar
v_T_HT, u_T_HT = isotherm(T_HT, p_min, p_max_L)                 # T=T_HT isotherm
v_T_min, u_T_min = isotherm(T_triple, p_min, p_max_L)           # T=T_triple isotherm

# u_c iso-energy line (L/G supercritical boundary)
v_uc = np.geomspace(min(vL.min(), v_p_max.min()) * 0.5, max(v_p_min.max(), 1e2), 200)
u_uc_line = np.full_like(v_uc, u_c_mass)

# Region fills (approximate, for visual cue)
# L region: between sat-L (right) and p_max_L isobar (left), bounded above by u_c (supercritical part)
# G region: between sat-V (left) and p_min isobar (right), bounded below by u_c (supercritical), above by T_HT
# HT region: T > T_HT
# TP region: inside the dome

# Sat dome
ax.plot(vL, uL / 1e3, 'k-', lw=2, label='sat-L (x=0)')
ax.plot(vV, uV / 1e3, 'k-', lw=2)
ax.plot([1.0 / heos.rhomass_critical()] if hasattr(heos, 'rhomass_critical') else [1.0 / (D_crit * heos.molar_mass())],
        [u_c_mass / 1e3], 'r*', ms=20, label=f'critical pt (u_c={u_c_mass / 1e3:.0f} kJ/kg)')

# Region boundaries
ax.plot(v_p_max, u_p_max / 1e3, 'b-', lw=1.5, label=f'p = p_max ({p_max_L / 1e6:.0f} MPa)')
ax.plot(v_p_min, u_p_min / 1e3, 'g-', lw=1.5, label=f'p = p_min ({p_min:.0f} Pa)')
ax.plot(v_T_HT, u_T_HT / 1e3, 'm-', lw=1.5, label=f'T = T_HT ({T_HT:.0f} K)')
ax.plot(v_T_min, u_T_min / 1e3, 'c--', lw=1, alpha=0.7, label=f'T = T_min ({T_triple:.0f} K)')
ax.plot(v_uc, u_uc_line / 1e3, 'r--', lw=1.5, alpha=0.8,
        label=f'u = u_c (L/G supercritical split)')

# Annotations
v_min_global = min(vL.min(), v_p_max.min()) * 0.5
v_max_global = max(v_p_min.max(), 1e3)
ax.text(2e-3, 1500, 'L', fontsize=22, fontweight='bold', ha='center', va='center', color='C0')
ax.text(1e1, 2900, 'G', fontsize=22, fontweight='bold', ha='center', va='center', color='C2')
ax.text(1e1, 3700, 'HT', fontsize=22, fontweight='bold', ha='center', va='center', color='m')
ax.text(0.7, 1200, 'TP', fontsize=22, fontweight='bold', ha='center', va='center', color='gray')

ax.set_xlim(v_min_global, v_max_global)
ax.set_ylim(0, max(uL.max(), uV.max(), u_p_max.max() if u_p_max.size else 0,
                   u_T_HT.max() if u_T_HT.size else 0) / 1e3 * 1.1)
ax.set_xlabel('v / (m³/kg)  — log scale')
ax.set_ylabel('u / (kJ/kg)')
ax.set_title(f'{ref}  —  Kunick SBTL (v, u) region scheme   [mirror of IAPWS G13-15 Fig. 9]')
ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
ax.grid(True, which='both', alpha=0.3)

# ---------------------------------------------------------------------------
# Panel 2: (D, U) — CoolProp's native coords, log D
# ---------------------------------------------------------------------------
ax = axes[1]
ax.set_xscale('log')

DL = 1.0 / vL; DV = 1.0 / vV
ax.plot(DL, uL / 1e3, 'k-', lw=2, label='sat-L')
ax.plot(DV, uV / 1e3, 'k-', lw=2, label='sat-V')

D_p_max = 1.0 / v_p_max; D_p_min = 1.0 / v_p_min
D_T_HT = 1.0 / v_T_HT; D_T_min = 1.0 / v_T_min
ax.plot(D_p_max, u_p_max / 1e3, 'b-', lw=1.5, label=f'p_max = {p_max_L / 1e6:.0f} MPa')
ax.plot(D_p_min, u_p_min / 1e3, 'g-', lw=1.5, label=f'p_min = {p_min:.0f} Pa')
ax.plot(D_T_HT, u_T_HT / 1e3, 'm-', lw=1.5, label=f'T_HT = {T_HT:.0f} K')
ax.plot(D_T_min, u_T_min / 1e3, 'c--', lw=1, alpha=0.7, label=f'T_min = {T_triple:.0f} K')

# u_c boundary (horizontal in (D, u))
ax.axhline(u_c_mass / 1e3, color='r', ls='--', lw=1.5, alpha=0.8,
           label=f'u = u_c (L/G supercritical split)')
ax.axvline(D_crit * heos.molar_mass(), color='r', ls=':', lw=1, alpha=0.5,
           label=f'D = D_crit')

D_crit_mass = D_crit * heos.molar_mass()
ax.plot([D_crit_mass], [u_c_mass / 1e3], 'r*', ms=20, label='critical point')

# Region annotations
ax.text(8e2, 1500, 'L', fontsize=22, fontweight='bold', ha='center', va='center', color='C0')
ax.text(1e0, 2900, 'G', fontsize=22, fontweight='bold', ha='center', va='center', color='C2')
ax.text(1e0, 3700, 'HT', fontsize=22, fontweight='bold', ha='center', va='center', color='m')
ax.text(20, 1200, 'TP', fontsize=22, fontweight='bold', ha='center', va='center', color='gray')

ax.set_xlim(min(DL.min(), DV.min()) * 0.5, max(DL.max(), DV.max()) * 1.5)
ax.set_ylim(0, max(uL.max(), uV.max(), u_p_max.max() if u_p_max.size else 0,
                   u_T_HT.max() if u_T_HT.size else 0) / 1e3 * 1.1)
ax.set_xlabel(r'$\rho$ / (kg/m³)  — log scale')
ax.set_ylabel('u / (kJ/kg)')
ax.set_title(f'{ref}  —  same scheme in (ρ, u) coords used by CoolProp')
ax.legend(loc='lower left', fontsize=8, framealpha=0.9)
ax.grid(True, which='both', alpha=0.3)

fig.suptitle(
    f'Kunick SBTL region subdivision in (v, u) for {ref}\n'
    f'L · G · HT · TP — split by sat dome (x=0/1), p_max isobar, p_min isobar, T_HT isotherm, '
    f'and u_c (= {u_c_mass / 1e3:.1f} kJ/kg)',
    fontsize=11, y=1.00)
plt.tight_layout()
plt.savefig(output, dpi=110, bbox_inches='tight')
print(f'Saved {output}')
