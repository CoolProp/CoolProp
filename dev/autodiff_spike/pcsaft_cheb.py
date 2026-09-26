#!/usr/bin/env python3
"""
Experiment 4: can PC-SAFT use Chebyshev density rootfinding (Bell & Alpert 2018 style)?

Idea.  With eta = zeta_3 = rho q(T) and zeta_n = eta r_n(T), at fixed (T, x):
    F(eta) := eta Z = eta + eta^2 d(alphar)/d(eta)
and the pressure equation p = rho R T Z becomes
    F(eta) = p q(T) / (R T)                                     (target only shifts c_0)
F is a weighted sum of a few *universal* functions of eta (no fluid parameters), plus
two small parameter families:
    alphar = mbar [A1 phi1 + A2 phi2 + A3 phi3] + B1 eta I1(eta; mbar) + B2 eta C1 I2(eta; mbar)
             - sum_i x_i (m_i - 1) ln g(eta; a_i),      g = 1/(1-e) + a e/(1-e)^2 + (2a^2/9) e^2/(1-e)^3
I1 is linear in (1, c1m, c2m) -> 3 universal polynomials; C1 I2 needs a 1-D family in mbar;
ln g needs a 1-D family in a_i = 3 D_i r2 (b_i = 2 a_i^2 / 9 identically).
So the tables live on bounded rectangles: eta in [0, eta_max], mbar in [1, m_max], a in [a_lo, a_hi].

This script checks, in order:
  1. the decomposition reproduces Z from the direct model
  2. where the complex singularities sit -> Chebyshev degree per piece
  3. piecewise Chebyshev rootfinding finds every root a brute-force scan finds, incl. the
     spurious high-density roots (Privat et al. 2010)
"""
import numpy as np
from numpy.polynomial import chebyshev as C

NA, kB = 6.02214076e23, 1.380649e-23
R = NA * kB
A = np.array([[0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084],
              [-0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930],
              [-0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368]])
Bc = np.array([[0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612],
               [-0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346],
               [0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585]])
ETA_CP = np.pi / (3 * np.sqrt(2))  # close packing, 0.7405


# ------------------------------------------------------------------ direct model (complex-step Z)
def alphar_direct(p, T, rho, x):
    m, sig, eps = (np.asarray(p[k], dtype=float) for k in ("m", "sigma", "eps"))
    x = np.asarray(x, dtype=float)
    d = sig * (1 - 0.12 * np.exp(-3 * eps / T))
    rhoN = rho * NA * 1e-30
    z = [np.pi / 6 * rhoN * np.sum(x * m * d ** n) for n in range(4)]
    mbar = np.sum(x * m)
    om = 1 - z[3]
    ahs = (3 * z[1] * z[2] / om + z[2] ** 3 / (z[3] * om ** 2) + (z[2] ** 3 / z[3] ** 2 - z[0]) * np.log(om)) / z[0]
    ahc = mbar * ahs
    for i in range(len(m)):
        D = d[i] / 2
        g = 1 / om + D * 3 * z[2] / om ** 2 + D ** 2 * 2 * z[2] ** 2 / om ** 3
        ahc = ahc - x[i] * (m[i] - 1) * np.log(g)
    eta = z[3]
    c1, c2 = (mbar - 1) / mbar, (mbar - 1) / mbar * (mbar - 2) / mbar
    I1 = sum((A[0, k] + c1 * A[1, k] + c2 * A[2, k]) * eta ** k for k in range(7))
    I2 = sum((Bc[0, k] + c1 * Bc[1, k] + c2 * Bc[2, k]) * eta ** k for k in range(7))
    sij = 0.5 * (sig[:, None] + sig[None, :])
    eij = np.sqrt(eps[:, None] * eps[None, :])
    pre = np.outer(x * m, x * m) * sij ** 3
    s1, s2 = np.sum(pre * eij) / T, np.sum(pre * eij ** 2) / T ** 2
    C1 = 1 / (1 + mbar * (8 * eta - 2 * eta ** 2) / om ** 4
              + (1 - mbar) * (20 * eta - 27 * eta ** 2 + 12 * eta ** 3 - 2 * eta ** 4) / (om * (2 - eta)) ** 2)
    return ahc - 2 * np.pi * rhoN * I1 * s1 - np.pi * rhoN * mbar * C1 * I2 * s2


def Z_direct(p, T, rho, x):
    h = 1e-30 * rho
    return 1 + rho * alphar_direct(p, T, rho + 1j * h, x).imag / h


# ------------------------------------------------------------------ eta-decomposition
def scalars(p, T, x):
    """Everything that depends on (T, x) only."""
    m, sig, eps = (np.asarray(p[k], dtype=float) for k in ("m", "sigma", "eps"))
    x = np.asarray(x, dtype=float)
    d = sig * (1 - 0.12 * np.exp(-3 * eps / T))
    s = [np.sum(x * m * d ** n) for n in range(4)]
    r0, r1, r2 = s[0] / s[3], s[1] / s[3], s[2] / s[3]
    mbar = np.sum(x * m)
    sij = 0.5 * (sig[:, None] + sig[None, :])
    eij = np.sqrt(eps[:, None] * eps[None, :])
    pre = np.outer(x * m, x * m) * sij ** 3
    E1, E2 = np.sum(pre * eij), np.sum(pre * eij ** 2)
    return dict(q=np.pi / 6 * NA * 1e-30 * s[3], mbar=mbar, A1=3 * r1 * r2 / r0, A2=r2 ** 3 / r0, A3=r2 ** 3 / r0 - 1,
                B1=-12 * E1 / (T * s[3]), B2=-6 * mbar * E2 / (T ** 2 * s[3]),
                a=3 * (d / 2) * r2, w=x * (m - 1), c1m=(mbar - 1) / mbar, c2m=(mbar - 1) / mbar * (mbar - 2) / mbar)


# universal eta-functions, as eta^2 * d/deta of the alphar pieces (complex step in eta)
def _cs(f, e):
    h = 1e-30
    return (e ** 2) * f(e + 1j * h).imag / h


phi1 = lambda e: e / (1 - e)
phi2 = lambda e: e / (1 - e) ** 2
phi3 = lambda e: np.log(1 - e)
I1r = [lambda e, r=r: e * sum(A[r, k] * e ** k for k in range(7)) for r in range(3)]  # eta * I1 basis


def C1I2(e, mbar):
    c1, c2 = (mbar - 1) / mbar, (mbar - 1) / mbar * (mbar - 2) / mbar
    I2 = sum((Bc[0, k] + c1 * Bc[1, k] + c2 * Bc[2, k]) * e ** k for k in range(7))
    om = 1 - e
    C1 = 1 / (1 + mbar * (8 * e - 2 * e ** 2) / om ** 4 + (1 - mbar) * (20 * e - 27 * e ** 2 + 12 * e ** 3 - 2 * e ** 4) / (om * (2 - e)) ** 2)
    return e * C1 * I2


def lng(e, a):
    om = 1 - e
    return np.log(1 / om + a * e / om ** 2 + 2 * a * a / 9 * e ** 2 / om ** 3)


def F_decomposed(S, e):
    """F(eta) = eta Z, assembled from the universal pieces with (T, x) weights."""
    F = e.astype(float) + S["mbar"] * (S["A1"] * _cs(phi1, e) + S["A2"] * _cs(phi2, e) + S["A3"] * _cs(phi3, e))
    F += S["B1"] * (_cs(I1r[0], e) + S["c1m"] * _cs(I1r[1], e) + S["c2m"] * _cs(I1r[2], e))
    F += S["B2"] * _cs(lambda z: C1I2(z, S["mbar"]), e)
    for wi, ai in zip(S["w"], S["a"]):
        F -= wi * _cs(lambda z: lng(z, ai), e)
    return F


FLUIDS = {
    "propane": (dict(m=[2.0020], sigma=[3.6184], eps=[208.11]), [1.0]),
    "C1C2C3": (dict(m=[1.0, 1.6069, 2.0020], sigma=[3.7039, 3.5206, 3.6184], eps=[150.03, 191.42, 208.11]), [0.5, 0.3, 0.2]),
    # Gross & Sadowski 2001 n-decane; a long, asymmetric mixture with methane
    "C1-nC10": (dict(m=[1.0, 4.6627], sigma=[3.7039, 3.8384], eps=[150.03, 243.87]), [0.7, 0.3]),
}


def check_decomposition():
    print("1. decomposition vs direct model: max |F_dec - eta Z_direct| / |eta Z|")
    for name, (p, x) in FLUIDS.items():
        worst = 0
        for T in (150.0, 300.0, 600.0):
            S = scalars(p, T, x)
            etas = np.linspace(1e-4, 0.72, 200)
            Fd = F_decomposed(S, etas)
            Fz = np.array([e * Z_direct(p, T, e / S["q"], x) for e in etas])
            worst = max(worst, np.max(np.abs(Fd - Fz) / np.abs(Fz)))
        print(f"   {name:8s} {worst:.1e}")


# ------------------------------------------------------------------ singularities -> degree
def singularities(mbar, a):
    """Complex singularities of the eta-pieces: eta = 1 (pole/branch), zeros of 1/C1, zeros of g."""
    e = np.polynomial.Polynomial([0, 1])
    om = 1 - e
    # 1/C1 * (1-e)^4 (2-e)^2 = (1-e)^4(2-e)^2 + m (8e-2e^2)(2-e)^2 + (1-m)(20e-27e^2+12e^3-2e^4)(1-e)^2
    den = om ** 4 * (2 - e) ** 2 + mbar * (8 * e - 2 * e ** 2) * (2 - e) ** 2 + (1 - mbar) * (20 * e - 27 * e ** 2 + 12 * e ** 3 - 2 * e ** 4) * om ** 2
    # g (1-e)^3 = (1-e)^2 + a e (1-e) + (2a^2/9) e^2
    gnum = om ** 2 + a * e * om + 2 * a * a / 9 * e ** 2
    return dict(C1=den.roots(), g=gnum.roots())


def bernstein_rho(z, lo, hi):
    """Bernstein-ellipse parameter of the singularity z for the interval [lo, hi]."""
    u = (2 * z - lo - hi) / (hi - lo)
    r = u + np.sqrt(u * u - 1 + 0j)
    return max(abs(r), 1 / abs(r))


def degree_needed(sings, lo, hi, tol=1e-13):
    rho = min(bernstein_rho(z, lo, hi) for z in sings)
    return int(np.ceil(np.log(1 / tol) / np.log(rho))), rho


def check_singularities():
    print("\n2. nearest complex singularities of the eta-pieces (eta = 1 is always one)")
    for mbar in (1.0, 2.0, 4.66, 10.0, 30.0):
        s = singularities(mbar, 1.5)
        near = sorted(s["C1"], key=lambda z: abs(z - 0.37))[:2]
        print(f"   mbar={mbar:5.2f}: 1/C1 zeros nearest the domain: " + ", ".join(f"{z:.3f}" for z in near))
    for a in (1.0, 1.5, 2.0, 3.0):
        s = singularities(2.0, a)
        print(f"   a={a:4.2f}: g zeros: " + ", ".join(f"{z:.3f}" for z in s["g"]))
    print("\n   degree for 1e-13 on [0, 0.74] split into k equal pieces (worst piece), mbar in {1..30}, a in {1..3}:")
    all_s = [1.0 + 0j]
    for mbar in np.linspace(1, 30, 30):
        all_s += list(singularities(mbar, 1.5)["C1"])
    for a in np.linspace(1, 3, 21):
        all_s += list(singularities(2.0, a)["g"])
    for k in (1, 2, 4, 8, 16):
        edges = np.linspace(0, 0.74, k + 1)
        worst = max(degree_needed(all_s, edges[i], edges[i + 1])[0] for i in range(k))
        print(f"   k={k:2d} pieces: n = {worst}")


# ------------------------------------------------------------------ rootfinding
def cheb_pieces(f, edges, n):
    return [(lo, hi, C.chebinterpolate(lambda u, lo=lo, hi=hi: f(lo + (hi - lo) * (u + 1) / 2), n)) for lo, hi in zip(edges[:-1], edges[1:])]


def cheb_roots(pieces, target):
    roots = []
    for lo, hi, c in pieces:
        c = c.copy()
        c[0] -= target
        for r in C.chebroots(c):
            if abs(r.imag) < 1e-8 and -1 - 1e-10 <= r.real <= 1 + 1e-10:
                roots.append(lo + (hi - lo) * (r.real + 1) / 2)
    roots = sorted(roots)
    out = []
    for r in roots:  # merge duplicates at piece joins
        if not out or abs(r - out[-1]) > 1e-9:
            out.append(r)
    return out


def scan_roots(f, target, lo=1e-8, hi=0.74, n=200001):
    e = np.linspace(lo, hi, n)
    g = f(e) - target
    idx = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
    out = []
    for i in idx:  # bisect
        a, b = e[i], e[i + 1]
        for _ in range(60):
            mid = 0.5 * (a + b)
            if np.sign(f(np.array([mid]))[0] - target) == np.sign(f(np.array([a]))[0] - target):
                a = mid
            else:
                b = mid
        out.append(0.5 * (a + b))
    return out


def check_roots(k=8, n=12):
    print(f"\n3. all roots of F(eta) = p q/(RT) on [0, 0.74], {k} pieces x degree {n}, vs brute-force scan")
    edges = np.linspace(0, 0.74, k + 1)
    worst_fit = 0
    for name, (p, x) in FLUIDS.items():
        for T in (150.0, 250.0, 350.0, 500.0):
            S = scalars(p, T, x)
            f = lambda e: F_decomposed(S, np.atleast_1d(e))
            pieces = cheb_pieces(f, edges, n)
            ee = np.linspace(0, 0.74, 4001)[1:]
            fit = np.concatenate([C.chebval(2 * (ee[(ee >= lo) & (ee <= hi)] - lo) / (hi - lo) - 1, c) for lo, hi, c in pieces])
            ref = np.concatenate([f(ee[(ee >= lo) & (ee <= hi)]) for lo, hi, _ in pieces])
            worst_fit = max(worst_fit, np.max(np.abs(fit - ref) / np.maximum(np.abs(ref), 1e-300)))
            for pMPa in (0.01, 0.1, 1.0, 10.0, 100.0):
                target = pMPa * 1e6 * S["q"] / (R * T)
                rc, rs = cheb_roots(pieces, target), scan_roots(f, target)
                ok = len(rc) == len(rs) and all(abs(a - b) < 1e-10 * max(b, 1e-6) + 1e-14 for a, b in zip(rc, rs))
                tag = "" if ok else "   <-- MISMATCH"
                if len(rs) > 1 or not ok:
                    print(f"   {name:8s} T={T:5.0f} p={pMPa:6.2f} MPa  cheb {['%.6f' % r for r in rc]}  scan {['%.6f' % r for r in rs]}{tag}")
    print(f"   (single-root cases matched silently)  worst relative fit error of F on [0,0.74]: {worst_fit:.1e}")


# ------------------------------------------------------------------ pole removal + adaptive pieces
def P_C1(e, mbar):
    """Polynomial with 1/C1 = P / ((1-e)^4 (2-e)^2); its zeros are the C1 poles."""
    om = 1 - e
    return om ** 4 * (2 - e) ** 2 + mbar * (8 * e - 2 * e ** 2) * (2 - e) ** 2 + (1 - mbar) * (20 * e - 27 * e ** 2 + 12 * e ** 3 - 2 * e ** 4) * om ** 2


def build_edges(sings, n, tol=1e-13, hi_end=0.74):
    """Greedy: widest pieces from 0 upward whose Bernstein bound gives degree <= n."""
    edges, lo = [0.0], 0.0
    while lo < hi_end:
        a, b = lo, hi_end
        if degree_needed(sings, lo, b, tol)[0] <= n:
            edges.append(b)
            break
        for _ in range(60):
            mid = 0.5 * (a + b)
            if degree_needed(sings, lo, mid, tol)[0] <= n:
                a = mid
            else:
                b = mid
        edges.append(a)
        lo = a
    return np.array(edges)


def check_pole_removal(n=12):
    print(f"\n4. pieces needed for degree n={n}, tol 1e-13 (greedy edges from the singularity set)")
    print("   mbar   a     F (C1 poles present)   G = P^2 F (C1 poles removed)")
    for mbar in (1.0, 2.0, 4.66, 10.0, 30.0):
        for a in (1.0, 1.5, 3.0):
            sg = [1.0 + 0j] + list(singularities(mbar, a)["g"])
            sF = sg + list(singularities(mbar, a)["C1"])
            print(f"   {mbar:5.2f} {a:4.1f}   {len(build_edges(sF, n)) - 1:6d}                 {len(build_edges(sg, n)) - 1:6d}")


def g_pieces(S, edges, n):
    """Per piece: Chebyshev coefficients (in eta) of P^2 F = eta * [P^2 Z] and of P^2.
    Q = P^2 Z is O(1) and is what gets fitted; the factor eta is applied exactly
    (Chebyshev multiply-by-x in the piece variable), so the absolute error near eta = 0
    scales with eta and the gas root keeps full relative precision."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        emap = lambda u: lo + (hi - lo) * (u + 1) / 2
        Q = C.chebinterpolate(lambda u: P_C1(emap(u), S["mbar"]) ** 2 * F_decomposed(S, emap(u)) / emap(u), n)
        # eta = (hi+lo)/2 + (hi-lo)/2 * u  ->  eta*Q exactly
        PF = C.chebadd(C.chebmulx(Q) * (hi - lo) / 2, Q * (hi + lo) / 2)
        P2 = C.chebinterpolate(lambda u: P_C1(emap(u), S["mbar"]) ** 2, 12)  # degree-12 polynomial: exact
        out.append((lo, hi, PF, P2))
    return out


def g_roots(pieces, target):
    roots = []
    for lo, hi, cPF, cP2 in pieces:
        c = C.chebsub(cPF, target * cP2)
        for r in C.chebroots(c):
            if abs(r.imag) < 1e-8 and -1 - 1e-10 <= r.real <= 1 + 1e-10:
                roots.append(lo + (hi - lo) * (r.real + 1) / 2)
    out = []
    for r in sorted(roots):
        if not out or abs(r - out[-1]) > 1e-9 * max(r, 1e-6):
            out.append(r)
    return out


def check_roots_G(n=12):
    print(f"\n5. rootfinding on G = P^2 (F - target), greedy pieces at n={n}; relative root error vs scan+bisection")
    cases = dict(FLUIDS)
    cases["nC10"] = (dict(m=[4.6627], sigma=[3.8384], eps=[243.87]), [1.0])
    worst_err, worst_fit, nmulti, maxroots, total, npieces = 0.0, 0.0, 0, 0, 0, []
    for name, (p, x) in cases.items():
        for T in (40.0, 60.0, 80.0, 100.0, 150.0, 250.0, 350.0, 500.0, 800.0):
            S = scalars(p, T, x)
            sg = [1.0 + 0j] + [z for ai in S["a"] for z in singularities(S["mbar"], ai)["g"]]
            edges = build_edges(sg, n)
            npieces.append(len(edges) - 1)
            pieces = g_pieces(S, edges, n)
            f = lambda e: F_decomposed(S, np.atleast_1d(e))
            ee = np.linspace(1e-6, 0.74, 3001)
            G = P_C1(ee, S["mbar"]) ** 2 * f(ee)
            Gfit = np.empty_like(ee)
            for lo, hi, cPF, _ in pieces:
                k = (ee >= lo) & (ee <= hi)
                Gfit[k] = C.chebval(2 * (ee[k] - lo) / (hi - lo) - 1, cPF)
            worst_fit = max(worst_fit, np.max(np.abs(Gfit - G) / np.abs(G)))
            for pMPa in (1e-4, 1e-3, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0):
                target = pMPa * 1e6 * S["q"] / (R * T)
                rc, rs = g_roots(pieces, target), scan_roots(f, target)
                total += 1
                maxroots = max(maxroots, len(rs))
                if len(rc) != len(rs):
                    print(f"   COUNT MISMATCH {name} T={T} p={pMPa}: cheb {rc} scan {rs}")
                    continue
                for a, b in zip(rc, rs):
                    worst_err = max(worst_err, abs(a - b) / b)
                if len(rs) > 3:
                    print(f"   {len(rs)} roots: {name} T={T:.0f} K p={pMPa:g} MPa  eta = " + ", ".join(f"{r:.5f}" for r in rs))
                nmulti += len(rs) > 1
    print(f"   {total} (fluid, T, p) cases; {nmulti} with >1 root; max roots seen {maxroots}")
    print(f"   worst relative root error {worst_err:.1e}; worst pointwise relative fit error of P^2 F {worst_fit:.1e}; pieces per (T,x): {min(npieces)}-{max(npieces)}")


if __name__ == "__main__":
    check_decomposition()
    check_singularities()
    check_pole_removal()
    check_roots_G()
