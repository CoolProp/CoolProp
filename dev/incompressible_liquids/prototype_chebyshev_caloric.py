"""Prototype: Chebyshev basis for the incompressible caloric fits (rho, cp).

Answers the question "if the data fits switch to Chebyshev polynomials, can
we still get entropy?" -- empirically, on real shipped fluids, before any
C++ or schema work.

The current backend integrates the centered-monomial cp fit analytically:
h = Int cp dT (Polynomial2DFrac::integral) and s = Int cp/T dT
(fracIntCentral: a binomial expansion with one ln T term). Both are exact
but the entropy path is numerically fragile at high order and historically
bred the T == Tbase singularity bugs (#1578). In a Chebyshev basis:

  h        the antiderivative of a Chebyshev series is another Chebyshev
           series (stable recurrence, no centering, no singularity);
  s        two independent routes, cross-checked against each other and
           against adaptive quadrature:
           route A: re-expand cp(T)/T by sampling at Chebyshev-Lobatto
                    nodes (cp/T is analytic on [Tmin,Tmax], T > 0, so the
                    expansion converges geometrically), then integrate;
           route B: exact closed form. With T = alpha*u + beta mapping
                    [-1,1] -> [Tmin,Tmax]: cp(T)/T dT = cp(u)/(u-u0) du,
                    u0 = -beta/alpha (< -1, outside the domain). Chebyshev
                    division (numpy chebdiv) by (u - u0) gives quotient q
                    and scalar remainder r with cp(u) = (u-u0) q(u) + r, so
                    Int cp/T dT = Int q du + r * ln(T) + const
                    -- the stable Chebyshev analogue of the current
                    fracIntCentralDvector, with exactly one log term.

The composition direction is unchanged: a low-order ordinary polynomial in
(x - xbase) per Chebyshev T-coefficient, collapsed to an effective 1D
T-expansion at the queried composition -- the same row-combination the
current Polynomial2DFrac does, so the eventual C++ structure carries over.

Run as a script for the per-fluid ENTROPY table (the headline result) plus
fit-quality/conditioning tables, or under pytest for the pass/fail checks.
Needs numpy + scipy only, like the rest of the fitting pipeline.
"""

import math
import sys

import numpy as np

try:
    from scipy.integrate import quad
except ImportError:  # pragma: no cover
    quad = None

from numpy.polynomial import chebyshev as ncheb

from CPIncomp import getPureFluids, getSolutionFluids, getSecCoolFluids

# ----------------------------------------------------------------------
# Chebyshev caloric model
# ----------------------------------------------------------------------

# Route-A auxiliary expansion of cp(T)/T: cp/T is analytic on [Tmin, Tmax]
# (the 1/T pole sits at u0 = -(Tmax+Tmin)/(Tmax-Tmin) < -1, outside the
# mapped domain), so its Chebyshev coefficients decay geometrically at rate
# ~|u0 + sqrt(u0^2 - 1)|^-n. The rate depends on how close u0 is to -1 --
# wide relative T-ranges (LiqNa: 400..2500 K, u0 = -1.38) converge slower
# than narrow ones -- so the degree is chosen adaptively until the
# coefficient tail is at rounding level, rather than fixed.
AUX_TAIL_TARGET = 1e-13
AUX_MAX_DEGREE = 96


class ChebyshevCaloric2D:
    """rho/cp as Chebyshev-in-T (x) monomial-in-(x - xbase) fits, plus the
    derived caloric properties h and s, mirroring the formulas of
    IncompressibleBackend::raw_calc_hmass / raw_calc_smass."""

    def __init__(self, Tvec, xvec, rho_grid, cp_grid, deg_T, deg_x):
        Tvec = np.asarray(Tvec, dtype=float)
        xvec = np.asarray(xvec, dtype=float)
        self.Tmin, self.Tmax = float(np.min(Tvec)), float(np.max(Tvec))
        self.xbase = 0.5 * (float(np.min(xvec)) + float(np.max(xvec)))
        self.deg_T, self.deg_x = deg_T, deg_x
        self.rho_C, self.rho_rms, self.rho_cond = self._fit2d(Tvec, xvec, rho_grid, deg_T, deg_x)
        self.cp_C, self.cp_rms, self.cp_cond = self._fit2d(Tvec, xvec, cp_grid, deg_T, deg_x)

    # -- fitting -------------------------------------------------------

    def _design(self, T, x, deg_T, deg_x):
        u = (2.0 * T - (self.Tmax + self.Tmin)) / (self.Tmax - self.Tmin)
        cheb_cols = ncheb.chebvander(u, deg_T)  # (npts, deg_T+1)
        xpow = np.vander(x - self.xbase, deg_x + 1, increasing=True)  # (npts, deg_x+1)
        # column (k, j) = T_k(u) * (x - xbase)^j
        return (cheb_cols[:, :, None] * xpow[:, None, :]).reshape(len(T), -1)

    def _fit2d(self, Tvec, xvec, grid, deg_T, deg_x):
        TT, XX = np.meshgrid(Tvec, xvec, indexing="ij")
        T, x, z = TT.ravel(), XX.ravel(), np.asarray(grid, dtype=float).ravel()
        mask = np.isfinite(z)
        A = self._design(T[mask], x[mask], deg_T, deg_x)
        coeffs, *_ = np.linalg.lstsq(A, z[mask], rcond=None)
        rms = float(np.sqrt(np.mean(np.square(A @ coeffs - z[mask]))))
        cond = float(np.linalg.cond(A))
        return coeffs.reshape(deg_T + 1, deg_x + 1), rms, cond

    # -- 2D -> 1D collapse at fixed composition (Polynomial2DFrac-style) --

    def _collapse(self, C, x):
        dx = np.power(x - self.xbase, np.arange(C.shape[1]))
        return ncheb.Chebyshev(C @ dx, domain=[self.Tmin, self.Tmax])

    def rho(self, T, x=0.0):
        return self._collapse(self.rho_C, x)(T)

    def cp(self, T, x=0.0):
        return self._collapse(self.cp_C, x)(T)

    def drhodT(self, T, x=0.0):
        return self._collapse(self.rho_C, x).deriv()(T)

    # -- caloric building blocks ----------------------------------------

    def h_T_part(self, x=0.0):
        """Indefinite Int cp dT as a domain-mapped Chebyshev (exact)."""
        return self._collapse(self.cp_C, x).integ()

    def s_T_part_routeA(self, x=0.0):
        """Indefinite Int cp/T dT: Lobatto re-expansion of cp/T, integrated.

        The auxiliary degree grows until the relative coefficient tail hits
        rounding level, so the geometric-decay claim is enforced, not assumed.
        """
        cp1d = self._collapse(self.cp_C, x)
        degree = self.deg_T + 8
        while True:
            aux = ncheb.Chebyshev.interpolate(lambda T: cp1d(T) / T, degree, domain=[self.Tmin, self.Tmax])
            tail = np.max(np.abs(aux.coef[-3:])) / np.max(np.abs(aux.coef))
            if tail < AUX_TAIL_TARGET or degree >= AUX_MAX_DEGREE:
                return aux.integ(), tail
            degree += 16

    def s_T_part_routeB(self, x=0.0):
        """Exact closed form via Chebyshev division: returns callable S(T).

        cp(u)/(u - u0) = q(u) + r/(u - u0), so
        Int cp/T dT = Int q du + r*ln(u - u0) = Int q du + r*ln(T/alpha).

        Exact in exact arithmetic, but the remainder r is the cp series
        evaluated at u0, OUTSIDE [-1,1], where T_n(u0) grows like
        |u0 + sqrt(u0^2-1)|^n -- so for narrow T-ranges far from 0 K (|u0|
        large) the log term and the quotient integral are individually large
        and cancel, costing digits. Empirically ~1e-10 relative on IceEA
        (u0 ~ -26) vs ~1e-16 for route A. This is the same cancellation
        family that makes the current monomial fracIntCentral fragile, which
        is exactly why route A (sampling, no extrapolated evaluation) is the
        recommended production path and route B only the cross-check.
        """
        cp1d = self._collapse(self.cp_C, x)
        alpha = 0.5 * (self.Tmax - self.Tmin)
        beta = 0.5 * (self.Tmax + self.Tmin)
        u0 = -beta / alpha  # < -1 for any Tmin > 0
        quotient, remainder = ncheb.chebdiv(cp1d.coef, np.array([-u0, 1.0]))
        r = remainder[0]  # == cp-series evaluated at u0
        q_int = ncheb.Chebyshev(ncheb.chebint(quotient), domain=[-1.0, 1.0])

        def S(T):
            u = (np.asarray(T, dtype=float) - beta) / alpha
            return q_int(u) + r * np.log(u - u0)  # u - u0 == T/alpha > 0

        return S

    # -- full properties (raw_calc_hmass / raw_calc_smass mirrors) -------

    def hmass_raw(self, T, p, x=0.0):
        rho = self.rho(T, x)
        return self.h_T_part(x)(T) + p * (1.0 / rho) * (1.0 + (T / rho) * self.drhodT(T, x))

    def smass_raw(self, T, p, x=0.0, route="A"):
        s_T = self.s_T_part_routeA(x)[0] if route == "A" else self.s_T_part_routeB(x)
        return s_T(T) + p * self.drhodT(T, x) / self.rho(T, x) ** 2

    def smass(self, T, p, x=0.0, T_ref=293.15, p_ref=101325.0, s_ref=0.0, route="A"):
        """Reference-pinned entropy VALUE, exactly like IncompressibleBackend:
        one raw evaluation at the state, one at the reference, subtract."""
        return s_ref + self.smass_raw(T, p, x, route) - self.smass_raw(T_ref, p_ref, x, route)


# ----------------------------------------------------------------------
# Fluid data extraction (reuses the CPIncomp classes unchanged)
# ----------------------------------------------------------------------


_FLUID_CACHE = {}


def _get_fluid(name):
    if not _FLUID_CACHE:
        for getter in (getPureFluids, getSolutionFluids, getSecCoolFluids):
            for obj in getter():
                _FLUID_CACHE[obj.name] = obj
    return _FLUID_CACHE[name]


def _grids(fluid, prop):
    """(Tvec, xvec, grid) for one property, from data or from the committed fit."""
    data = getattr(fluid, prop)
    Tvec = np.asarray(fluid.temperature.data, dtype=float).ravel()
    xcand = fluid.concentration.data
    xvec = np.asarray(xcand, dtype=float).ravel() if xcand is not None and np.size(xcand) else np.array([0.0])
    if data.data is not None:
        grid = np.asarray(data.data, dtype=float).reshape(len(Tvec), len(xvec))
        return Tvec, xvec, grid
    raise ValueError("{0}.{1} has no tabular data".format(fluid.name, prop))


def _seccool_grids(fluid, prop):
    """SecCool fluids load their grids during fitFluid(); reuse those arrays."""
    data = getattr(fluid, prop)
    if data.data is None:
        fluid.fitFluid()
    data = getattr(fluid, prop)
    Tvec = np.asarray(data.xData, dtype=float).ravel()
    xvec = np.asarray(data.yData, dtype=float).ravel()
    grid = np.asarray(data.data, dtype=float).reshape(len(Tvec), len(xvec))
    return Tvec, xvec, grid


def build_model(name, deg_T=8, deg_x=None):
    fluid = _get_fluid(name)
    loader = _seccool_grids if hasattr(fluid, "fitFluid") else _grids
    T_rho, x_rho, rho = loader(fluid, "density")
    T_cp, x_cp, cp = loader(fluid, "specific_heat")
    if not (np.array_equal(T_rho, T_cp) and np.array_equal(x_rho, x_cp)):
        raise ValueError("{0}: rho and cp grids differ".format(name))
    if deg_x is None:
        deg_x = 0 if len(x_rho) == 1 else min(5, len(x_rho) - 1)
    deg_T_eff = min(deg_T, len(T_rho) - 1)
    return ChebyshevCaloric2D(T_rho, x_rho, rho, cp, deg_T_eff, deg_x), fluid


FLUIDS = ["Water", "LiqNa", "LiBr", "IceEA"]  # smooth pure / huge T-range / 2D solution / pathological cp


def _x_probe(model):
    # a representative off-base composition (pure fluids: 0)
    return model.xbase * 1.2 if model.xbase != 0.0 else 0.0


# ----------------------------------------------------------------------
# The validation battery (pytest-collectable; also drives the CLI tables)
# ----------------------------------------------------------------------


def _entropy_report(name, deg_T=8):
    model, _fluid = build_model(name, deg_T=deg_T)
    x = _x_probe(model)
    T0, T1 = model.Tmin + 0.05 * (model.Tmax - model.Tmin), model.Tmax - 0.05 * (model.Tmax - model.Tmin)
    p = 1.0e6

    sA_fun, aux_tail = model.s_T_part_routeA(x)
    sB_fun = model.s_T_part_routeB(x)
    dA = sA_fun(T1) - sA_fun(T0)
    dB = sB_fun(T1) - sB_fun(T0)
    cp1d = model._collapse(model.cp_C, x)
    dQ, _err = quad(lambda T: cp1d(T) / T, T0, T1, epsabs=1e-13, epsrel=1e-13, limit=200)

    scale = abs(dQ)
    routeAB = abs(dA - dB) / scale
    routeAQ = abs(dA - dQ) / scale

    # ds/dT == cp/T at the fit level (route A derivative closes the loop)
    Tm = 0.5 * (T0 + T1)
    dsdT = sA_fun.deriv()(Tm)
    dsdT_rel = abs(dsdT - cp1d(Tm) / Tm) / abs(cp1d(Tm) / Tm)

    # dh/dT == cp
    h_fun = model.h_T_part(x)
    dhdT_rel = abs(h_fun.deriv()(Tm) - cp1d(Tm)) / abs(cp1d(Tm))

    # Maxwell: (ds/dp)_T == -(d(1/rho)/dT)_p ; s is linear in p, so the
    # p-slope of smass_raw must equal drhodT/rho^2 exactly
    slope = (model.smass_raw(Tm, p, x) - model.smass_raw(Tm, 0.0, x)) / p
    maxwell_rel = abs(slope - model.drhodT(Tm, x) / model.rho(Tm, x) ** 2) / abs(slope)

    # reference pinning: s(T_ref, p_ref) == s_ref bit-consistently
    pin = model.smass(293.15 if model.Tmin < 293.15 < model.Tmax else Tm, 101325.0, x,
                      T_ref=293.15 if model.Tmin < 293.15 < model.Tmax else Tm)

    return {
        "fluid": name, "x": x, "T0": T0, "T1": T1,
        "routeA_vs_routeB": routeAB, "routeA_vs_quad": routeAQ,
        "dsdT_vs_cpT": dsdT_rel, "dhdT_vs_cp": dhdT_rel,
        "maxwell": maxwell_rel, "ref_pin": abs(pin), "aux_tail": aux_tail,
    }


def _fit_report(name):
    rows, seen = [], set()
    for deg_T in (4, 8, 12):
        model, _ = build_model(name, deg_T=deg_T)
        if model.deg_T in seen:  # short data vectors cap the usable degree
            continue
        seen.add(model.deg_T)
        rows.append((model.deg_T, model.cp_rms, model.rho_rms, model.cp_cond, ""))
    return rows


# ---- pytest entry points ---------------------------------------------


def test_entropy_routes_agree():
    for name in FLUIDS:
        rep = _entropy_report(name)
        # Route A vs adaptive quadrature of the same fit is the hard gate.
        assert rep["routeA_vs_quad"] < 1e-12, (name, rep)
        # Route B is the independent closed-form cross-check; its documented
        # cancellation (see s_T_part_routeB) caps its own accuracy near 1e-10
        # on narrow-range fluids, so gate the agreement at 1e-9.
        assert rep["routeA_vs_routeB"] < 1e-9, (name, rep)


def test_derivative_consistency():
    for name in FLUIDS:
        rep = _entropy_report(name)
        assert rep["dhdT_vs_cp"] < 1e-12, (name, rep)
        assert rep["dsdT_vs_cpT"] < 1e-10, (name, rep)
        assert rep["maxwell"] < 1e-12, (name, rep)


def test_reference_state_pins_to_zero():
    for name in FLUIDS:
        rep = _entropy_report(name)
        assert rep["ref_pin"] == 0.0, (name, rep)  # bit-exact cancellation


def test_aux_expansion_tail_decays():
    for name in FLUIDS:
        rep = _entropy_report(name)
        assert rep["aux_tail"] < 1e-13, (name, rep)


def test_ice_slurry_cp_improves_with_order():
    # The point of the basis switch for pathological cp(T): raising the
    # T-order keeps improving the fit while conditioning stays sane. (The
    # IceEA grid has few T-points, so the usable degree is data-capped;
    # compare the lowest available degree against the highest.)
    rows = _fit_report("IceEA")
    rms = {deg: cp_rms for deg, cp_rms, *_rest in rows}
    lo, hi = min(rms), max(rms)
    assert hi > lo, rows
    assert rms[hi] < 0.5 * rms[lo], rows
    assert max(r[3] for r in rows) < 1e9, rows


# ---- CLI --------------------------------------------------------------

if __name__ == "__main__":
    if quad is None:
        sys.exit("scipy is required")

    print("== ENTROPY (headline): route A (cp/T re-expansion) vs route B (chebdiv closed form) vs quadrature ==")
    hdr = "{0:8s} {1:>9s} {2:>12s} {3:>12s} {4:>12s} {5:>12s} {6:>12s} {7:>10s} {8:>10s}"
    print(hdr.format("fluid", "x", "|A-B|/s", "|A-quad|/s", "ds/dT~cp/T", "dh/dT~cp", "Maxwell", "ref-pin", "aux tail"))
    failures = 0
    for name in FLUIDS:
        r = _entropy_report(name)
        ok = r["routeA_vs_quad"] < 1e-12 and r["routeA_vs_routeB"] < 1e-9 and r["ref_pin"] == 0.0
        failures += 0 if ok else 1
        print("{0:8s} {1:9.4f} {2:12.2e} {3:12.2e} {4:12.2e} {5:12.2e} {6:12.2e} {7:10.1e} {8:10.1e}{9}".format(
            r["fluid"], r["x"], r["routeA_vs_routeB"], r["routeA_vs_quad"], r["dsdT_vs_cpT"],
            r["dhdT_vs_cp"], r["maxwell"], r["ref_pin"], r["aux_tail"], "" if ok else "  <-- FAIL"))

    print("\n== cp fit quality vs Chebyshev T-order (RMS in J/kg/K; design-matrix condition number) ==")
    for name in FLUIDS:
        print(name)
        for deg, cp_rms, rho_rms, cond, note in _fit_report(name):
            print("   deg_T={0:2d}  cp RMS={1:10.4g}  rho RMS={2:10.4g}  cond={3:10.3g}  {4}".format(
                deg, cp_rms, rho_rms, cond, note))

    sys.exit(1 if failures else 0)
