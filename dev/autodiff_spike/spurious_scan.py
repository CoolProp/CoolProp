#!/usr/bin/env python3
"""
Where do spurious (>3) density roots of PC-SAFT appear, relative to the triple point, and do
the Liang et al. universal constants remove them?  Brute-force root count of
eta Z(eta) = p q/(RT) on eta in (0, 0.74) over a (T, p) grid.

Caveat: the Liang constant sets were published with refitted pure-component parameters; here
they are paired with the Gross & Sadowski 2001 parameters, so this is only a qualitative check
of the dispersion polynomials' high-eta behaviour.
"""
import numpy as np
import pcsaft_cheb as M

GS = (M.A.copy(), M.Bc.copy())
LIANG2012 = (np.array([[0.836215101666, 2.201683842453, -11.25210310939, 37.841836899902, -68.035304263396, 69.952369867326, -42.828905226651],
                       [-0.411727190935913, 2.37426400571265, -21.0620603419144, 105.65718855671, -298.665894225644, 468.695983731173, -316.673589664169],
                       [0.0319867672916212, -2.75137756155503, 25.7581175334397, -103.082737163044, 239.569856365622, -320.622085430506, 186.50494276364]]),
             np.array([[0.627209841336118, 4.02517816132384, -22.6660554011051, 70.0153445172765, -129.548046066679, 150.197401680241, -86.9664928046989],
                       [-0.622507280536237, 1.86478114256654, -22.5660800172653, 153.069895818654, -321.710866534244, 453.083735030445, -233.232187026907],
                       [-0.0303061275320169, 2.59554209415371, 11.5803588817289, -145.915305288352, 354.901071174401, -151.675970200232, -282.837925415568]]))
LIANG2014 = (np.array([[0.961597, 0.414449, 0.689253, -7.43899, 31.8755, -54.8833, 27.3613],
                       [-0.333416, 0.358440, -0.219088, 0.285586, 5.88256, -22.4931, 21.2109],
                       [-0.0632483, 0.287134, -0.105309, 7.16607, -37.5008, 68.8002, -40.3177]]),
             np.array([[0.548398, 2.07176, -1.84013, -29.9683, 160.445, -206.106, 51.6201],
                       [-0.277226, 0.0621105, 4.29640, -39.9453, 214.911, -232.563, 60.9734],
                       [-0.150684, -0.507877, 1.24227, -48.7052, 63.1155, 205.380, -262.437]]))

FLUIDS = {  # (params, x, triple point of the pure fluid / heaviest component [K])
    "methane": (dict(m=[1.0], sigma=[3.7039], eps=[150.03]), [1.0], 90.7),
    "propane": (dict(m=[2.0020], sigma=[3.6184], eps=[208.11]), [1.0], 85.5),
    "n-decane": (dict(m=[4.6627], sigma=[3.8384], eps=[243.87]), [1.0], 243.5),
    "C1C2C3": (dict(m=[1.0, 1.6069, 2.0020], sigma=[3.7039, 3.5206, 3.6184], eps=[150.03, 191.42, 208.11]), [0.5, 0.3, 0.2], None),
    "C1-nC10 70/30": (dict(m=[1.0, 4.6627], sigma=[3.7039, 3.8384], eps=[150.03, 243.87]), [0.7, 0.3], None),
}


def nroots(f, target, e):
    g = f(e) - target
    return int(np.sum(np.sign(g[:-1]) != np.sign(g[1:])))


def scan(consts):
    M.A[:], M.Bc[:] = consts
    e = np.linspace(1e-9, 0.74, 40001)
    ps = np.logspace(-4, 3, 22)  # MPa
    res = {}
    for name, (p, x, Ttp) in FLUIDS.items():
        Tmax4, etas = None, []
        for T in np.arange(20.0, 400.0, 2.0):
            S = M.scalars(p, T, x)
            Fe = M.F_decomposed(S, e)
            f = lambda ee, Fe=Fe: Fe  # precomputed on the grid
            for pM in ps:
                t = pM * 1e6 * S["q"] / (M.R * T)
                g = Fe - t
                k = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
                if len(k) > 3:
                    Tmax4 = T
                    etas.append(e[k[-1]])
        res[name] = (Tmax4, (min(etas), max(etas)) if etas else None, Ttp)
    M.A[:], M.Bc[:] = GS
    return res


for label, consts in (("Gross & Sadowski 2001", GS), ("Liang 2012", LIANG2012), ("Liang 2014", LIANG2014)):
    print(f"\n{label}: highest T [K] with >3 roots (p = 1e-4..1e3 MPa), eta range of the extra root, triple point")
    for name, (Tmax, er, Ttp) in scan(consts).items():
        tp = f"{Ttp:.1f}" if Ttp else "  -  "
        if Tmax is None:
            print(f"   {name:14s}  none (20-400 K)                         T_tp {tp}")
        else:
            print(f"   {name:14s}  T <= {Tmax:5.0f}   extra root eta {er[0]:.3f}-{er[1]:.3f}   T_tp {tp}   T/T_tp = {Tmax / Ttp if Ttp else float('nan'):.2f}")
