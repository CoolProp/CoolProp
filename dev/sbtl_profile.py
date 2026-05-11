#!/usr/bin/env python3
"""SBTL / BICUBIC / HEOS timing + accuracy profile in three boxes.

  Box 1: queries inside the SBTL HEOS-fallback box around the critical point
         (T_crit ± 5 % T_c, p_crit ± 0.4 MPa).
  Box 2: strict supercritical (T well above T_crit, p > p_crit + a margin).
  Box 3: strict subcritical (T below T_crit, p < p_crit − a margin).

Reports for each box:
  * accuracy vs HEOS for HmassP_INPUTS lookups (median, 99th, max)
  * per-call timing for HEOS, BICUBIC, SBTL on HmassP_INPUTS

Usage:  python3 sbtl_profile.py [FluidName]
"""
import sys
import time

import numpy as np

import CoolProp.CoolProp as CP

ref = sys.argv[1] if len(sys.argv) > 1 else 'CO2'

heos = CP.AbstractState('HEOS', ref)
T_crit = heos.T_critical()
p_crit = heos.p_critical()

# Prime backends with a known-good single-phase state.
prime_T = max(0.5 * T_crit, heos.Ttriple() + 5)
prime_p = max(p_crit * 0.1, 1e5)
heos.update(CP.PT_INPUTS, prime_p, prime_T)
print(f'Building BICUBIC backend...', flush=True)
bicu = CP.AbstractState('BICUBIC&HEOS', ref); bicu.update(CP.PT_INPUTS, prime_p, prime_T)
print(f'Building SBTL backend (~minutes for first build)...', flush=True)
sbtl = CP.AbstractState('SBTL&HEOS', ref); sbtl.update(CP.PT_INPUTS, prime_p, prime_T)

print(f'\n{ref}:  Tc={T_crit:.2f} K   pc={p_crit / 1e5:.3f} bar', flush=True)


def accuracy_box(T_lim, p_lim, label, N=80, verbose_worst=False):
    """ρ relative error vs HEOS on a (T, p) grid, single-phase only."""
    Tg = np.linspace(*T_lim, N)
    Pg = np.geomspace(*p_lim, N)
    err_b, err_s = [], []
    worst = []
    for T in Tg:
        for P in Pg:
            try:
                heos.update(CP.PT_INPUTS, P, T)
                if 0 < heos.Q() < 1: continue
                rho = heos.rhomass(); h = heos.hmass()
            except Exception:
                continue
            try:
                bicu.update(CP.HmassP_INPUTS, h, P)
                err_b.append(abs(bicu.rhomass() - rho) / rho)
            except Exception:
                pass
            try:
                sbtl.update(CP.HmassP_INPUTS, h, P)
                e = abs(sbtl.rhomass() - rho) / rho
                err_s.append(e)
                if verbose_worst and e > 1e-3:
                    worst.append((e, T, P, h, rho, sbtl.rhomass()))
            except Exception:
                pass
    eb = np.array(err_b); es = np.array(err_s)
    if eb.size:
        print(f'  BICUBIC HP: med={np.median(eb):.2e}  99th={np.percentile(eb, 99):.2e}  max={eb.max():.2e}  N={len(eb)}',
              flush=True)
    if es.size:
        print(f'  SBTL    HP: med={np.median(es):.2e}  99th={np.percentile(es, 99):.2e}  max={es.max():.2e}  N={len(es)}',
              flush=True)
    if verbose_worst and worst:
        worst.sort(key=lambda x: -x[0])
        print('  Top-10 worst SBTL points:', flush=True)
        for w in worst[:10]:
            print(f'    err={w[0]:.3e} T={w[1]:.1f} P={w[2] / 1e5:.2f}bar h={w[3] / 1e3:.1f} eos={w[4]:.2f} sbtl={w[5]:.2f}',
                  flush=True)


def timing_box(T_lim, p_lim, label, N=2000):
    np.random.seed(0)
    Tarr = np.random.uniform(*T_lim, N)
    Parr = np.exp(np.random.uniform(np.log(p_lim[0]), np.log(p_lim[1]), N))
    h_arr = np.zeros(N); valid = np.zeros(N, dtype=bool)
    for k in range(N):
        try:
            heos.update(CP.PT_INPUTS, Parr[k], Tarr[k])
            if 0 < heos.Q() < 1: continue
            h_arr[k] = heos.hmass(); valid[k] = True
        except Exception:
            pass
    Tarr = Tarr[valid]; Parr = Parr[valid]; h_arr = h_arr[valid]
    M = len(Tarr)
    print(f'\n  {label}  (M={M} valid points)', flush=True)
    for name, BCK in [('BICUBIC', bicu), ('SBTL', sbtl), ('HEOS', heos)]:
        t0 = time.perf_counter()
        for k in range(M):
            try: BCK.update(CP.HmassP_INPUTS, h_arr[k], Parr[k])
            except Exception: pass
        dt = (time.perf_counter() - t0) * 1e6 / M
        print(f'    {name:<8} HmassP {dt:7.3f} µs/call ({1e6 / dt:>9,.0f} calls/s)', flush=True)


# Three boxes scaled from the fluid's critical properties.
boxes = [
    ('Box 1: INSIDE HEOS-fallback (around critical)',
     (T_crit * 0.95, T_crit * 1.05),
     (p_crit - 0.4e6, p_crit + 0.4e6)),
    ('Box 2: STRICT supercritical',
     (T_crit * 1.15, T_crit * 1.65),
     (p_crit + 2e6, min(p_crit * 9, heos.pmax() * 0.5))),
    ('Box 3: STRICT subcritical',
     (max(heos.Ttriple() + 1, T_crit * 0.4), T_crit * 0.97),
     (max(p_crit * 0.005, 1e4), p_crit - 2e6)),
]

for box_label, T_lim, p_lim in boxes:
    print('\n' + '=' * 72)
    print(box_label, flush=True)
    print(f'  T ∈ [{T_lim[0]:.1f}, {T_lim[1]:.1f}] K   '
          f'p ∈ [{p_lim[0] / 1e5:.3f}, {p_lim[1] / 1e5:.1f}] bar', flush=True)
    print('=' * 72, flush=True)
    accuracy_box(T_lim, p_lim, box_label, verbose_worst=box_label.startswith('Box 3'))
    timing_box(T_lim, p_lim, box_label)
