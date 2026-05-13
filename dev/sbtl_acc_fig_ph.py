#!/usr/bin/env python3
"""SBTL/BICUBIC accuracy maps in (h, P) coordinates for any pure fluid.

Usage:
  python3 sbtl_acc_fig_ph.py [FluidName] [/path/to/output.png]

Three rows: A) around critical, B) subcritical, C) supercritical.
Three columns: BICUBIC HmassP, BICUBIC PT, SBTL HmassP.  Window (h, P)
ranges auto-scale from the fluid's critical properties.  Solid-region
points (T < T_melt) are filtered so the accuracy maps reflect the
actual fluid region, not extrapolation noise.

For pure fluids, SBTL HmassP routes through the coordinate-aligned
NormalizedPHTable (LIQUID / VAPOR / SUPER) plus a small HEOS-fallback
box around (h_crit, p_crit).  See dev/sbtl_normalized_ph_design.md for
the design.  Sibling PT figure: dev/sbtl_acc_fig_pt.py.
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import CoolProp.CoolProp as CP

ref = sys.argv[1] if len(sys.argv) > 1 else 'CO2'
output_path = sys.argv[2] if len(sys.argv) > 2 else f'/tmp/sbtl_acc_fig_ph_{ref.lower()}.png'

# ---------------------------------------------------------------------------
# Fluid metadata
# ---------------------------------------------------------------------------
heos = CP.AbstractState('HEOS', ref)
T_triple = heos.Ttriple()
T_crit = heos.T_critical()
p_crit = heos.p_critical()
T_min = max(T_triple, heos.Tmin())
heos.update(CP.DmolarT_INPUTS, heos.rhomolar_critical(), T_crit)
h_crit = heos.hmass()
heos.update(CP.QT_INPUTS, 0.0, T_triple * 1.0001)
p_triple = heos.p()
p_max_eos = heos.pmax()

print(f'Fluid: {ref}', flush=True)
print(f'  Tc = {T_crit:.3f} K   pc = {p_crit / 1e5:.4f} bar   h_crit = {h_crit / 1e3:.3f} kJ/kg', flush=True)
print(f'  T_triple = {T_triple:.3f} K   p_triple ≈ {p_triple / 1e5:.6f} bar   pmax = {p_max_eos / 1e6:.0f} MPa', flush=True)

# ---------------------------------------------------------------------------
# Auto-windows scaled from critical properties.  For fluids with h_crit close
# to zero (rare reference choices), fall back to absolute scaling on |h_crit|
# plus an offset so the windows have non-degenerate width.
# ---------------------------------------------------------------------------
def h_range(lo_frac, hi_frac):
    """Return (h_lo, h_hi) scaled around h_crit."""
    if h_crit > 0:
        return (lo_frac * h_crit, hi_frac * h_crit)
    span = max(abs(h_crit), 100e3)
    return (h_crit - lo_frac * span, h_crit + hi_frac * span)


# Window A: tight box around the critical point.  ~±40% in h, ±35% in P.
A_h = h_range(0.6, 1.4)
A_p = (0.95 * p_crit, 1.36 * p_crit)

# Window B: full subcritical PH plane up to (just below) p_crit.
B_h = h_range(0.0, 1.5)
B_p = (max(p_triple * 1.001, p_crit * 1e-4), 0.99 * p_crit)

# Window C: supercritical, span ~ pcrit to a few pcrit (clamped to pmax).
C_h = h_range(1.0, 3.5)
C_p = (1.02 * p_crit, min(p_max_eos, 4.0 * p_crit))

# ---------------------------------------------------------------------------
# Backends — primed with a known-good state.
# ---------------------------------------------------------------------------
prime_T = max(T_triple * 1.5, 1.5 * T_min)
prime_p = max(p_triple * 2, 1e5)
print('Building BICUBIC...', flush=True)
bicu = CP.AbstractState('BICUBIC&HEOS', ref); bicu.update(CP.PT_INPUTS, prime_p, prime_T)
print('Building SBTL...', flush=True)
sbtl = CP.AbstractState('SBTL&HEOS', ref); sbtl.update(CP.PT_INPUTS, prime_p, prime_T)


def probe(BCK, h_lim, p_lim, mode, NH=200, NP=200):
    """mode = 'HP' uses HmassP; 'PT' uses HEOS to get T then PT_INPUTS.

    Skips solid-region points (T < T_triple or T < T_melt(p) where the
    melting ancillary is defined) so the maps don't show extrapolation
    artifacts in the solid phase as 'errors'.
    """
    has_melt = heos.has_melting_line()
    hg = np.linspace(*h_lim, NH)
    pg = np.geomspace(*p_lim, NP)
    HH, PP = np.meshgrid(hg, pg)
    err = np.full_like(HH, np.nan)
    for i in range(NP):
        for j in range(NH):
            h = HH[i, j]; P = PP[i, j]
            try:
                heos.update(CP.HmassP_INPUTS, h, P)
                if 0.0 < heos.Q() < 1.0:  # skip 2-phase
                    continue
                rho = heos.rhomass()
                T_eos = heos.T()
            except Exception:
                continue
            if T_eos < T_triple:
                continue
            if has_melt:
                try:
                    if T_eos < heos.melting_line(CP.iT, CP.iP, P):
                        continue
                except Exception:
                    pass
            try:
                if mode == 'HP':
                    BCK.update(CP.HmassP_INPUTS, h, P)
                else:
                    BCK.update(CP.PT_INPUTS, P, T_eos)
                err[i, j] = abs(BCK.rhomass() - rho) / rho
            except Exception:
                # query landed outside table envelope; skip this point
                pass
    return HH, PP, err


def saturation_curve_ph():
    """Return (h_L, h_V, P) along the sat curve in PH coords."""
    Ts = np.linspace(T_triple + 0.01, T_crit - 0.001, 400)
    h_L, h_V, p_arr = [], [], []
    for t in Ts:
        try:
            heos.update(CP.QT_INPUTS, 0.0, t); h_L.append(heos.hmass()); p_arr.append(heos.p())
            heos.update(CP.QT_INPUTS, 1.0, t); h_V.append(heos.hmass())
        except Exception:
            # query landed outside table envelope; skip this point
            pass
    return np.array(h_L), np.array(h_V), np.array(p_arr)


# ---------------------------------------------------------------------------
# Draw the figure
# ---------------------------------------------------------------------------
windows = [
    ('A) Critical region', A_h, A_p, True),
    ('B) Subcritical',     B_h, B_p, False),
    ('C) Supercritical',   C_h, C_p, True),
]

cNorm = colors.LogNorm(vmin=1e-7, vmax=1e0)
fig, axes = plt.subplots(3, 3, figsize=(16, 14))
h_L_full, h_V_full, p_full = saturation_curve_ph()

for row_i, (w_label, h_lim, p_lim, allow_PT) in enumerate(windows):
    h_lo, h_hi = h_lim
    p_lo, p_hi = p_lim
    label_suffix = (
        f'(h~{h_lo / 1e3:.0f}-{h_hi / 1e3:.0f} kJ/kg, '
        f'P={p_lo / 1e5:.2f}-{p_hi / 1e5:.0f} bar)'
    )
    for col_i, (BCK, label, mode) in enumerate([
        (bicu, 'BICUBIC HmassP', 'HP'),
        (bicu, 'BICUBIC PT',     'PT'),
        (sbtl, 'SBTL HmassP',    'HP'),
    ]):
        ax = axes[row_i, col_i]
        skip_PT = (BCK is sbtl and mode == 'PT') or (mode == 'PT' and not allow_PT)
        if skip_PT:
            ax.set_axis_off()
            ax.set_title(f'{w_label} {label_suffix}: {label}\n(skipped)', fontsize=10)
            continue
        HH, PP, err = probe(BCK, h_lim, p_lim, mode)
        im = ax.pcolormesh(HH / 1e3, PP / 1e5, np.maximum(err, 1e-12),
                           cmap='plasma', shading='auto', norm=cNorm)
        plt.colorbar(im, ax=ax, label=r'$|\rho_{\mathrm{tab}}/\rho_{\mathrm{HEOS}} - 1|$',
                     shrink=0.85)
        if h_L_full.size:
            mask = (p_full >= p_lim[0]) & (p_full <= p_lim[1])
            ax.plot(h_L_full[mask] / 1e3, p_full[mask] / 1e5, 'w-', lw=1.5, alpha=0.85)
            ax.plot(h_V_full[mask] / 1e3, p_full[mask] / 1e5, 'w-', lw=1.5, alpha=0.85)
        ax.axhline(p_crit / 1e5, color='w', ls=':', lw=1, alpha=0.6)
        # Critical-region HEOS-fallback box (Box 1: pc ± 0.1 MPa, h_crit ± 30 %)
        rect_p_lo = (p_crit - 0.1e6) / 1e5
        rect_p_hi = (p_crit + 0.1e6) / 1e5
        rect_h_lo = h_crit * 0.70 / 1e3
        rect_h_hi = h_crit * 1.30 / 1e3
        ax.plot([rect_h_lo, rect_h_hi, rect_h_hi, rect_h_lo, rect_h_lo],
                [rect_p_lo, rect_p_lo, rect_p_hi, rect_p_hi, rect_p_lo],
                'c-', lw=1.5, alpha=0.9)
        ax.set_xlim(h_lo / 1e3, h_hi / 1e3)
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
            ax.set_xlabel('h / (kJ/kg)')

plt.suptitle(
    f'{ref}  ρ relative error  vs HEOS  (default 200×200 tables)\n'
    f'Probed in (h, P) coordinates; saturation dome white; cyan rectangle = SBTL HEOS-fallback box',
    fontsize=12, y=1.0)
plt.tight_layout()
plt.savefig(output_path, dpi=110, bbox_inches='tight')
print(f'Saved {output_path}')
