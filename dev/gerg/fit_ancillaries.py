#!/usr/bin/env python3
"""Fit saturation ancillaries to the GERG-2004/GERG-2008 PURE equations of state.

WHY THIS EXISTS AT ALL
----------------------
CoolProp already ships saturation ancillaries for every fluid GERG names --
but those belong to the *reference* equation of state for that fluid
(Setzmann-Wagner for methane, Span-Wagner for carbon dioxide, IAPWS-95 for
water, ...), not to GERG's shortened technical form.  Borrowing them would
work numerically -- an ancillary is only a starting guess -- but the whole
point of the GERG backend is that nothing in it traces back to a
non-GERG equation.  So they are refitted here against the GERG pure EOS
itself, using teqp as the reference implementation.

An ancillary is a SEED, NOT AN ANSWER.  Nothing in the backend returns an
ancillary value to the caller: `FlashRoutines::QT_flash` hands them to
`SaturationSolvers::saturation_T_pure_Maxwell`, which iterates to the true
GERG saturation state.  (This is also why no superancillary may ever be
attached to a GERG fluid: `FlashRoutines::sat_superanc_path_applies` routes
straight to the Chebyshev expansion and RETURNS THAT as the answer, which for
a GERG fluid would be a reference-EOS saturation density labelled GERG-2008.
See src/Backends/GERG/GERGBackend.cpp and task-10.)

WHAT IT PRODUCES
----------------
src/Backends/GERG/GERGAncillaries.h -- a generated C++ table consumed by
GERG::get_ancillary() (declared in GERGData.h, defined in GERGBackend.cpp).
The schema mirrors dev/fluids/*.json's ANCILLARIES block exactly, so
make_gerg_fluid can build CoolProp's own SaturationAncillaryFunction from it
and the same evaluation path as every other fluid is reused:

    rhoL   type "rhoLnoexp"  using_tau_r = false
           rho' = reducing_value * (1 + sum_i n_i theta^t_i)
    rhoV   type "rhoV"       using_tau_r = true
           rho'' = reducing_value * exp(T_r/T * sum_i n_i theta^t_i)
    pV     type "pV"         using_tau_r = true
           p_sat = reducing_value * exp(T_r/T * sum_i n_i theta^t_i)

    theta = 1 - T/T_r

REDUCING POINT != CRITICAL POINT (load-bearing, and the reason a t = 0 term
appears in every fit).  For most CoolProp fluids the ancillary reducing point
IS the EOS critical point, so rhoL(T_c) = rhoV(T_c) = rho_c and p(T_c) = p_c
fall out of the functional form for free (theta = 0 kills every t > 0 term).
That is NOT true for GERG: Table A3.5's reducing parameters are fitted
quantities, and the true critical point of GERG's short 12-term form can sit
far away from them -- measured with teqp's own critical solver, up to
+1.10 K in T_c (n-heptane) and -5.43% in rho_c (n-butane); see the
"reducing vs. true critical" column printed by --report.  Forcing the
ancillary through the tabulated reducing point would therefore put a
several-percent error into the seed exactly where the VLE iteration is
worst conditioned.  A theta^0 term makes the value at T = T_r a free fitted
parameter instead, so the fit follows the EOS rather than the table.

T_r = max(true T_c, tabulated T_reduce), and that maximum is load-bearing in
both directions:

  * It can never be BELOW the tabulated reducing temperature, because
    SaturationAncillaryFunction::evaluate returns NaN for T > T_r
    (Ancillaries.cpp, guarding pow(negative, fractional)) while QT_flash's
    upper saturation bound is the backend's T_critical(), i.e. the tabulated
    value.  A lower T_r would hand NaN seeds to the VLE solver in the sliver
    between the two.
  * Where the true critical temperature is the HIGHER of the pair (propane,
    n-butane, oxygen, argon, n-heptane, n-octane, and GERG-2004 isopentane --
    up to +1.10 K), using it puts theta = 0 exactly at the end of the
    saturation curve, which is what makes the classical theta^(1/2) scaling
    representable.  Forcing theta = 0 to land 1.1 K short of the critical
    point instead leaves a square-root branch point INSIDE the fit domain and
    costs a factor of five in accuracy (propane: 7.3e-3 that way versus
    1.4e-3 with T_r = Tc(true); run --report to reproduce).

reducing_value stays the TABULATED reducing density / the pressure at the
tabulated reducing state, so that the ancillary's reducing value and the
fluid's rhomolar_critical()/p_critical() are the same number.  The theta^0
term absorbs the resulting offset.

RUN (teqp 0.23.2; the system python3's teqp 0.15.3 predates the GERG factory
registration and fails with "Unknown kind:GERG2008resid"):

    python3 -m venv /tmp/gergenv
    /tmp/gergenv/bin/pip install "teqp>=0.18" numpy
    /tmp/gergenv/bin/python dev/gerg/fit_ancillaries.py \
        > src/Backends/GERG/GERGAncillaries.h
    uvx clang-format@18.1.8 -i src/Backends/GERG/GERGAncillaries.h

    # human-readable fit report (deviations per fluid) on stderr/stdout:
    /tmp/gergenv/bin/python dev/gerg/fit_ancillaries.py --report

REFIT WHEN THE COEFFICIENT TABLES CHANGE.  These fits are tied to the
specific GERG coefficient tables in GERGData.h/GERGBackend.cpp.  Editing a
pure-fluid n/t/d/l row or a Table A3.5 reducing parameter invalidates the
corresponding ancillary; rerun this script.

METHOD
------
1.  VLE trace.  For each of the 23 distinct pure EOS (21 GERG-2008 fluids,
    plus the GERG-2004 variants of carbon monoxide and isopentane, whose
    reducing parameters GERG-2008 changed), the saturation curve is traced
    from just below the true critical temperature downwards.  Each step is
    seeded from the previous one -- the technique ~/Code/fastchebpure's
    fastcheb.cpp uses for the same job -- with a from-scratch fallback:
    scan for the two spinodal densities, bisect on p between them with
    equal-chemical-potential as the residual, then polish with a damped
    2-D Newton on (rhoL, rhoV).  The from-scratch path needs no guess at
    all, so a step that fails never strands the trace.

2.  Regression.  Each ancillary is linear in n once the exponents t are
    fixed, so t is chosen by forward selection from a fixed candidate
    ladder (greedy, one term at a time, keeping whichever candidate most
    reduces the weighted sum of squares) and n falls out of a weighted
    linear least squares.  Weights make the objective RELATIVE error in
    rho/p rather than absolute error in the fitted combination.

3.  Reporting.  The maximum relative deviation over the traced points is
    printed per fluid and embedded in the generated header.  A fit worse
    than MAX_ACCEPTABLE_DEV aborts rather than shipping.
"""

import argparse
import math
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")  # teqp 0.23 emits FutureWarnings for top-level fns

import numpy as np  # noqa: E402
import teqp  # noqa: E402

# The GERG gas constant, teqp GERG.hpp:470 -- NOT CODATA R.  Matches
# GERGData.h's R_GERG; the pressures written into sat_min_liquid/sat_min_vapor
# would be wrong in the 7th digit with CODATA.
R_GERG = 8.314472

Z1 = np.array([1.0])

# GERG2004::component_names / GERG2008::component_names order, teqp
# GERG.hpp:431 and :973 -- the same order GERGData.h::component_names() uses.
GERG2004_NAMES = [
    "methane", "nitrogen", "carbondioxide", "ethane", "propane", "n-butane",
    "isobutane", "n-pentane", "isopentane", "n-hexane", "n-heptane", "n-octane",
    "hydrogen", "oxygen", "carbonmonoxide", "water", "helium", "argon",
]  # fmt: skip
GERG2008_NAMES = GERG2004_NAMES + ["hydrogensulfide", "n-nonane", "n-decane"]

# The only two components whose PURE EOS differs between the two models: both
# are reducing-parameter changes (GERGData.h pure_info_2008_overrides), and
# both move the saturation curve, so each needs its own fit.  Mirrors the
# base-plus-overrides shape of pure_info_2004/pure_info_2008_overrides, so the
# generated table is a 2004 base plus a 2008 override set.
CHANGED_IN_2008 = ["carbonmonoxide", "isopentane"]
NEW_IN_2008 = ["hydrogensulfide", "n-nonane", "n-decane"]

# Any fluid whose fit is worse than this over the traced range is a bug, not a
# fit -- Task 11's own test asserts the ancillary is within 2% of the
# converged saturation density, and the brief's ship gate is 1%.
MAX_ACCEPTABLE_DEV = 0.005  # 0.5% relative

# Candidate exponents for forward selection.  t = 0 is deliberately included
# (see the module docstring): it makes the value at T = T_r free, which is
# what lets the fit follow GERG's true critical point instead of the
# tabulated reducing point.  The half-integer ladder covers the classical
# (analytic-EOS) critical exponent 1/2 and its harmonics; the large exponents
# pick up the low-temperature tail.
CANDIDATE_T = [
    0.0, 0.25, 0.35, 0.5, 0.7, 0.9, 1.2, 1.5, 2.0, 2.6, 3.4, 4.5, 6.0, 8.0,
    11.0, 15.0, 20.0, 27.0,
]  # fmt: skip

# Two exponents that are nearly equal span nearly the same function, so the
# least squares can only tell them apart by giving them huge, opposite-signed
# coefficients -- a fit that is accurate on paper but evaluates as the
# difference of two ~1e3 terms whose sum is ~1.  Forward selection therefore
# refuses a candidate that sits too close to an exponent already chosen.  The
# ladder above is pre-spaced; this is the guard that keeps a future edit to it
# from silently reintroducing the problem.
MIN_EXPONENT_GAP = 0.15
MIN_EXPONENT_RATIO = 1.3

MAX_TERMS = 10

# Below this saturation pressure the vapour density underflows the useful
# range of a double (rho'' < 1e-9 mol/m^3 for every fluid here) and the
# equal-chemical-potential residual stops being resolvable, so the trace
# stops.  This is a numerical floor, NOT a statement about where the GERG EOS
# stops being valid.
P_FLOOR_PA = 1e-3

# Never trace below this fraction of the true critical temperature even if the
# solver would keep going: the fit has to spend its terms where the backend
# can actually be asked for a saturation state.
T_FLOOR_FRACTION = 0.30

N_TRACE_POINTS = 250


# ---------------------------------------------------------------------------
# EOS wrapper
# ---------------------------------------------------------------------------
class Pure:
    """One pure GERG EOS, with the handful of derivatives VLE needs.

    teqp's accessors are get_Ar<TAU order><DELTA order> and are already
    multiplied by the corresponding power of delta, i.e. get_Ar01 is
    delta*(dalpha^r/ddelta)_tau and get_Ar02 is delta^2*(d2alpha^r/ddelta2)_tau.
    """

    def __init__(self, year, name):
        self.year = year
        self.name = name
        self.m = teqp.make_model({"kind": f"GERG{year}resid", "model": {"names": [name]}})
        # The reducing state as tabulated (Table A3.5) -- what GERGData.h
        # stores and what the backend reports as its critical point.
        self.T_reduce = self.m.get_Tcvec()[0]
        self.rho_reduce = 1.0 / self.m.get_vcvec()[0]
        # The TRUE critical point of this short EOS, which is generally not
        # the same thing.  Bounds the trace, and sets T_r when it is the
        # larger of the two (see the module docstring).
        self.Tc_true, self.rhoc_true = self.m.solve_pure_critical(self.T_reduce, self.rho_reduce)
        self.T_r = max(self.Tc_true, self.T_reduce)

    def p(self, T, rho):
        return rho * R_GERG * T * (1.0 + self.m.get_Ar01(T, rho, Z1))

    def dpdrho(self, T, rho):
        return R_GERG * T * (1.0 + 2 * self.m.get_Ar01(T, rho, Z1) + self.m.get_Ar02(T, rho, Z1))

    def mu_over_RT(self, T, rho):
        """Chemical potential in RT units, dropping the terms that depend on T
        alone -- they are identical in both phases at a given T and cancel out
        of the equilibrium condition."""
        return math.log(rho) + self.m.get_Ar00(T, rho, Z1) + self.m.get_Ar01(T, rho, Z1)


# ---------------------------------------------------------------------------
# VLE
# ---------------------------------------------------------------------------
def _bisect(f, a, b, rtol=1e-14, maxit=300):
    fa, fb = f(a), f(b)
    if not (np.isfinite(fa) and np.isfinite(fb)) or fa * fb > 0:
        raise RuntimeError("no sign change in bracket")
    for _ in range(maxit):
        c = 0.5 * (a + b)
        fc = f(c)
        if not np.isfinite(fc):
            raise RuntimeError("non-finite residual")
        if fa * fc <= 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
        if abs(b - a) <= rtol * max(abs(c), 1e-300):
            break
    return 0.5 * (a + b)


def spinodals(P, T):
    """Densities bounding the mechanically unstable region (dp/drho < 0).

    Returned as (vapour-side, liquid-side).  None when the isotherm is
    monotone, i.e. there is no VLE at this temperature.
    """
    grid = np.geomspace(1e-9 * P.rho_reduce, 6.0 * P.rho_reduce, 900)
    d = np.array([P.dpdrho(T, r) for r in grid])
    neg = np.where(d < 0)[0]
    if len(neg) == 0 or neg[0] == 0 or neg[-1] == len(grid) - 1:
        return None
    i0, i1 = neg[0], neg[-1]
    return (
        _bisect(lambda r: P.dpdrho(T, r), grid[i0 - 1], grid[i0]),
        _bisect(lambda r: P.dpdrho(T, r), grid[i1], grid[i1 + 1]),
    )


def sat_from_scratch(P, T):
    """Guess-free saturation solve: bisect on p between the two spinodal
    pressures with equal chemical potential as the residual.  Slow (a few
    hundred EOS evaluations) but it cannot be stranded by a bad seed, which
    is what makes the stepping loop below safe to run unattended."""
    sp = spinodals(P, T)
    if sp is None:
        return None
    rho_sp_v, rho_sp_l = sp
    p_sp_v, p_sp_l = P.p(T, rho_sp_v), P.p(T, rho_sp_l)

    def rhoV_of(p):
        lo = min(1e-14 * P.rho_reduce, 0.5 * p / (R_GERG * T))
        for _ in range(400):
            if P.p(T, lo) < p:
                break
            lo *= 0.1
        return _bisect(lambda r: P.p(T, r) - p, lo, rho_sp_v)

    def rhoL_of(p):
        hi = rho_sp_l
        for _ in range(400):
            hi *= 1.2
            if P.p(T, hi) > p:
                break
        return _bisect(lambda r: P.p(T, r) - p, rho_sp_l, hi)

    def g(p):
        return P.mu_over_RT(T, rhoL_of(p)) - P.mu_over_RT(T, rhoV_of(p))

    # The saturation pressure lies between the two spinodal pressures.  The
    # liquid spinodal pressure is routinely negative at low temperature, in
    # which case the physical lower bound is p -> 0+.
    a = p_sp_l * (1 + 1e-12) if p_sp_l > 0 else p_sp_v * 1e-13
    b = p_sp_v * (1 - 1e-12)
    try:
        ga, gb = g(a), g(b)
    except RuntimeError:
        return None
    if ga * gb > 0:
        return None
    for _ in range(300):
        c = math.sqrt(a * b)  # geometric bisection: p spans many decades
        try:
            gc = g(c)
        except RuntimeError:
            return None
        if ga * gc <= 0:
            b, gb = c, gc
        else:
            a, ga = c, gc
        if b - a <= 1e-15 * b:
            break
    p = 0.5 * (a + b)
    try:
        return rhoL_of(p), rhoV_of(p)
    except RuntimeError:
        return None


def newton_polish(P, T, rhoL, rhoV, tol=1e-13, maxit=80):
    """Damped 2-D Newton on (rhoL, rhoV) with equal pressure and equal
    chemical potential.  Also the stepping path: seeded from the previous
    temperature's answer it converges in a handful of iterations."""
    for _ in range(maxit):
        f1 = (P.p(T, rhoL) - P.p(T, rhoV)) / (R_GERG * T * P.rho_reduce)
        f2 = P.mu_over_RT(T, rhoL) - P.mu_over_RT(T, rhoV)
        if not (np.isfinite(f1) and np.isfinite(f2)):
            return None
        if abs(f1) < tol and abs(f2) < tol:
            return (rhoL, rhoV) if rhoL > rhoV > 0 else None
        dpL, dpV = P.dpdrho(T, rhoL), P.dpdrho(T, rhoV)
        j11 = dpL / (R_GERG * T * P.rho_reduce)
        j12 = -dpV / (R_GERG * T * P.rho_reduce)
        j21 = dpL / (R_GERG * T * rhoL)
        j22 = -dpV / (R_GERG * T * rhoV)
        det = j11 * j22 - j12 * j21
        if det == 0 or not np.isfinite(det):
            return None
        dL = (-f1 * j22 + f2 * j12) / det
        dV = (-j11 * f2 + j21 * f1) / det
        s = 1.0
        while rhoL + s * dL <= 0 or rhoV + s * dV <= 0:
            s *= 0.5
            if s < 1e-10:
                return None
        rhoL += s * dL
        rhoV += s * dV
    return None


def trace(P):
    """Trace the saturation curve downwards from just below the true critical
    temperature.  Returns a list of (T, rhoL, rhoV, p), highest T first."""
    theta_top = 1e-4  # closer than this and the two roots are not separable in double precision
    T_floor = T_FLOOR_FRACTION * P.Tc_true
    theta_bot = 1.0 - T_floor / P.Tc_true
    rows = []
    prev = None
    for theta in np.geomspace(theta_top, theta_bot, N_TRACE_POINTS):
        T = P.Tc_true * (1.0 - theta)
        sol = newton_polish(P, T, *prev) if prev is not None else None
        if sol is None:
            sol = sat_from_scratch(P, T)
            if sol is not None:
                sol = newton_polish(P, T, *sol) or sol
        if sol is None:
            break
        rhoL, rhoV = sol
        if not (rhoL > rhoV > 0) or (rhoL - rhoV) / P.rho_reduce < 1e-7:
            break
        p = P.p(T, rhoV)
        if not np.isfinite(p) or p < P_FLOOR_PA:
            break
        rows.append((T, rhoL, rhoV, p))
        prev = sol
    return rows


# ---------------------------------------------------------------------------
# Regression
# ---------------------------------------------------------------------------
def _fit_linear(theta, y, w, exponents):
    A = np.column_stack([w * np.power(theta, t) for t in exponents])
    n, *_ = np.linalg.lstsq(A, w * y, rcond=None)
    resid = A @ n - w * y
    return n, float(resid @ resid)


def _well_separated(t, chosen):
    """True when exponent t is far enough from every already-chosen exponent
    to be numerically distinguishable (see MIN_EXPONENT_GAP)."""
    for t0 in chosen:
        if abs(t - t0) < MIN_EXPONENT_GAP:
            return False
        lo, hi = min(t, t0), max(t, t0)
        if lo > 0 and hi / lo < MIN_EXPONENT_RATIO:
            return False
    return True


def forward_select(theta, y, w, max_terms=MAX_TERMS):
    """Greedy forward selection of exponents, then a final least squares.

    The exponents are the nonlinear part of the fit; once they are fixed the
    coefficients are a linear least squares.  Greedy selection is enough here
    because the shipped gate is the measured deviation, not optimality.
    """
    chosen = []
    best_n = None
    best_sse = np.inf
    for _ in range(max_terms):
        cand_best = None
        for t in CANDIDATE_T:
            if not _well_separated(t, chosen):
                continue
            try:
                n, sse = _fit_linear(theta, y, w, chosen + [t])
            except np.linalg.LinAlgError:
                continue
            if not np.isfinite(sse):
                continue
            if cand_best is None or sse < cand_best[2]:
                cand_best = (t, n, sse)
        if cand_best is None or cand_best[2] >= best_sse * (1 - 1e-10):
            break
        chosen.append(cand_best[0])
        best_n, best_sse = cand_best[1], cand_best[2]
    order = np.argsort(chosen)
    return [chosen[i] for i in order], [float(best_n[i]) for i in order]


def fit_rhoL(rows, T_r, reducing_value):
    T = np.array([r[0] for r in rows])
    rho = np.array([r[1] for r in rows])
    theta = 1.0 - T / T_r
    # rho' = reducing_value*(1 + S), so S = rho'/reducing_value - 1; weighting
    # by reducing_value/rho' turns the residual in S into relative error in rho'.
    y = rho / reducing_value - 1.0
    w = reducing_value / rho
    t, n = forward_select(theta, y, w)
    pred = reducing_value * (1.0 + sum(ni * np.power(theta, ti) for ni, ti in zip(n, t)))
    return t, n, float(np.max(np.abs(pred / rho - 1.0)))


def _fit_exponential(T, value, T_r, reducing_value):
    theta = 1.0 - T / T_r
    tau = T_r / T
    # value = reducing_value*exp(tau*S) => S = ln(value/reducing_value)/tau, and
    # relative error in value is tau*|dS|, hence the weight.
    y = np.log(value / reducing_value) / tau
    w = tau
    t, n = forward_select(theta, y, w)
    pred = reducing_value * np.exp(tau * sum(ni * np.power(theta, ti) for ni, ti in zip(n, t)))
    return t, n, float(np.max(np.abs(pred / value - 1.0)))


def fit_rhoV(rows, T_r, reducing_value):
    return _fit_exponential(
        np.array([r[0] for r in rows]), np.array([r[2] for r in rows]), T_r, reducing_value
    )


def fit_pV(rows, T_r, reducing_value):
    return _fit_exponential(
        np.array([r[0] for r in rows]), np.array([r[3] for r in rows]), T_r, reducing_value
    )


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
class FluidFit:
    def __init__(self, year, name, strict=True):
        self.year = year
        self.name = name
        P = Pure(year, name)
        self.P = P
        rows = trace(P)
        if len(rows) < 30:
            raise RuntimeError(f"{year}/{name}: only {len(rows)} usable VLE points")
        rows = sorted(rows, key=lambda r: r[0])  # ascending T
        self.rows = rows
        self.Tmin = rows[0][0]
        # Tmax is T_r: SaturationAncillaryFunction only uses Tmax to bracket
        # invert(), and T_r is the largest temperature at which evaluate()
        # returns a number at all.
        self.Tmax = P.T_r
        # p at the tabulated reducing state, computed from the EOS -- exactly
        # what make_gerg_fluid stores in EOS.reduce.p and what the backend
        # reports as p_critical().  Reusing it keeps the ancillary's
        # reducing_value and the fluid's p_critical() the same number.
        self.p_reduce = P.p(P.T_reduce, P.rho_reduce)

        # Only points strictly below T_r can enter the fit: theta <= 0 there
        # and pow(theta, non-integer) is undefined for negative theta.  When
        # T_r is the true critical temperature this drops nothing but the
        # single point the trace already placed at Tc*(1 - 1e-4).
        fit_rows = [r for r in rows if r[0] < P.T_r * (1 - 1e-9)]
        if len(fit_rows) < 30:
            raise RuntimeError(f"{year}/{name}: only {len(fit_rows)} points below T_r")
        self.fit_rows = fit_rows

        self.t_rhoL, self.n_rhoL, self.dev_rhoL = fit_rhoL(fit_rows, P.T_r, P.rho_reduce)
        self.t_rhoV, self.n_rhoV, self.dev_rhoV = fit_rhoV(fit_rows, P.T_r, P.rho_reduce)
        self.t_pV, self.n_pV, self.dev_pV = fit_pV(fit_rows, P.T_r, self.p_reduce)

        self.worst = max(self.dev_rhoL, self.dev_rhoV, self.dev_pV)
        if strict and self.worst > MAX_ACCEPTABLE_DEV:
            raise RuntimeError(
                f"{year}/{name}: worst relative deviation {self.worst:.3e} exceeds "
                f"{MAX_ACCEPTABLE_DEV:.3e}; add terms or raise Tmin"
            )

    @property
    def sat_min(self):
        """(T, p, rhoL, rhoV) at the lowest traced temperature."""
        T, rhoL, rhoV, p = self.rows[0]
        return T, p, rhoL, rhoV


def _fmt(x):
    # float() first: repr(np.float64(...)) is "np.float64(57.169)" on numpy 2,
    # which is not a C++ literal.
    return repr(float(x))


def _vec(values):
    return "{" + ", ".join(_fmt(v) for v in values) + "}"


#: The aggregate initialisers emitted by _anc() and emit_row() are POSITIONAL,
#: so reordering either struct in GERGData.h would still compile and would
#: silently pair every value with the wrong field.  These are the orders this
#: generator emits; check_struct_field_order() below refuses to generate if
#: the header no longer agrees.
EMITTED_ANCILLARY_FIELDS = ["n", "t", "reducing_value", "T_r", "Tmin", "Tmax", "max_rel_dev", "type"]
EMITTED_SAT_END_FIELDS = ["T_K", "p_Pa", "rhoL_molm3", "rhoV_molm3"]

#: `struct NAME\n{ ... };`, non-greedy so it stops at the first closing brace.
_STRUCT_RE = r"struct\s+%s\s*\{(.*?)\}\s*;"
#: One leading type token, optionally templated (std::vector<double>), then the
#: comma-separated declarators.
_DECL_TYPE_RE = re.compile(r"^[A-Za-z_][\w:]*(?:<[^>]*>)?\s+")


def declared_struct_fields(blob, struct_name):
    """Field names of `struct struct_name`, in declaration order."""
    m = re.search(_STRUCT_RE % struct_name, blob, re.S)
    if m is None:
        raise SystemExit(f"fit_ancillaries.py: struct {struct_name} not found in GERGData.h")
    body = re.sub(r"//[^\n]*", "", m.group(1))          # line comments
    body = re.sub(r"/\*.*?\*/", "", body, flags=re.S)    # block comments
    fields = []
    for stmt in body.split(";"):
        stmt = " ".join(stmt.split())
        if not stmt:
            continue
        rest = _DECL_TYPE_RE.sub("", stmt, count=1)
        if rest == stmt:                                 # no type token -> not a declaration
            continue
        fields.extend(name.strip().lstrip("*&") for name in rest.split(","))
    return [f for f in fields if f]


def check_struct_field_order():
    """Refuse to generate if GERGData.h no longer matches what we emit.

    A reordered struct is invisible to the compiler here -- both aggregates are
    all-double apart from the two vectors and the trailing string -- so without
    this the corruption would only surface as wrong saturation numbers.
    """
    header = Path(__file__).resolve().parents[2] / "src" / "Backends" / "GERG" / "GERGData.h"
    blob = header.read_text(encoding="utf-8")
    for struct_name, emitted in (
        ("AncillaryCoeffs", EMITTED_ANCILLARY_FIELDS),
        ("SatEndState", EMITTED_SAT_END_FIELDS),
    ):
        declared = declared_struct_fields(blob, struct_name)
        if declared != emitted:
            raise SystemExit(
                f"fit_ancillaries.py: {struct_name} field order changed in {header}.\n"
                f"  declared: {declared}\n"
                f"  emitted:  {emitted}\n"
                "The brace-initialisers in _anc()/emit_row() are positional, so generating "
                "now would silently pair each value with the wrong field. Update them together."
            )


def _anc(n, t, reducing_value, T_r, Tmin, Tmax, dev):
    """One AncillaryCoeffs brace-initialiser.  Field order MUST match the
    struct declaration in GERGData.h: n, t, reducing_value, T_r, Tmin, Tmax,
    max_rel_dev.  `type` is the trailing member and is deliberately left
    default-initialised here -- get_ancillary() fills it in from `which`, so
    the table cannot disagree with the accessor about which form it is."""
    scalars = ", ".join(_fmt(v) for v in (reducing_value, T_r, Tmin, Tmax, dev))
    return "{" + _vec(n) + ", " + _vec(t) + ", " + scalars + "}"


def emit_row(f):
    T, p, rhoL, rhoV = f.sat_min
    body = ", ".join(
        [
            # T_r is P.T_r == max(true Tc, tabulated T_reduce) -- the SAME
            # temperature the fit was regressed against.  Emitting T_reduce
            # here instead compiles and evaluates fine, and silently shifts
            # theta for the seven fluids where the two differ.
            _anc(f.n_rhoL, f.t_rhoL, f.P.rho_reduce, f.P.T_r, f.Tmin, f.Tmax, f.dev_rhoL),
            _anc(f.n_rhoV, f.t_rhoV, f.P.rho_reduce, f.P.T_r, f.Tmin, f.Tmax, f.dev_rhoV),
            _anc(f.n_pV, f.t_pV, f.p_reduce, f.P.T_r, f.Tmin, f.Tmax, f.dev_pV),
            "{" + ", ".join(_fmt(v) for v in (T, p, rhoL, rhoV)) + "}",
        ]
    )
    return '{"' + f.name + '", {' + body + "}}"


HEADER = '''// GENERATED FILE -- do not edit by hand.
//
// Saturation ancillaries FITTED AGAINST THE GERG PURE EQUATIONS OF STATE,
// produced by dev/gerg/fit_ancillaries.py from teqp {teqp_version}
// (https://github.com/usnistgov/teqp).  CoolProp's build has no dependency on
// teqp; this header is committed, exactly like GERGReferenceValues.h.
//
// These are NOT CoolProp's shipped ancillaries for the same fluids.  Those
// belong to each fluid's REFERENCE equation of state (Setzmann-Wagner for
// methane, IAPWS-95 for water, ...), and a backend named GERG-2004/GERG-2008
// must not carry data traceable to a different equation.  They are refitted
// here against GERG's own shortened technical form.  REFIT WHENEVER A GERG
// COEFFICIENT TABLE CHANGES -- see dev/gerg/README.md.
//
// An ancillary is a SEED, not an answer: QT_flash hands rhoL/rhoV/pV to
// SaturationSolvers::saturation_T_pure_Maxwell, which iterates to the true
// GERG saturation state, and no ancillary value is ever returned to a caller.
// (The Task 11 test "GERG saturation is independent of the ancillary seed"
// pins that property by perturbing these coefficients and checking the
// converged answer does not move.)
//
// A t = 0 exponent appears in every fit on purpose.  GERG's reducing
// parameters (Table A3.5) are fitted quantities, not the true critical point
// of the short EOS -- they differ by up to +1.10 K and -5.43% in density --
// so the usual "ancillary is exact at the critical point" property does not
// hold and the value at T = T_r is left as a free fitted parameter.
//
// Worst relative deviation over the fitted range, across all 23 pure EOS:
// {worst_all:.3e} (rhoL {worst_rhoL:.3e}, rhoV {worst_rhoV:.3e}, pV {worst_pV:.3e}).
// Per-fluid values are stored in each AncillaryCoeffs::max_rel_dev.
#ifndef GERGANCILLARIES_H_
#define GERGANCILLARIES_H_

#include <map>
#include <string>

#include "GERGData.h"

namespace CoolProp {{
namespace GERG {{
namespace detail {{

/// The three ancillaries for one pure EOS plus the saturation state at the
/// low-temperature end of the fitted range.  Field order matches the
/// brace-initialised rows below; reordering silently corrupts every fluid.
struct AncillarySet
{{
    AncillaryCoeffs rhoL, rhoV, pV;
    SatEndState sat_min;
}};

/// Keyed on the GERG component name.  GERG-2004 is the base table and
/// GERG-2008 overrides it, mirroring pure_info_2004/pure_info_2008_overrides
/// in GERGData.h.
inline const std::map<std::string, AncillarySet>& ancillaries_2004() {{
    static const std::map<std::string, AncillarySet> data = {{
{rows_2004}
    }};
    return data;
}}

/// Entries GERG-2008 changes or adds relative to GERG-2004.  Carbon monoxide
/// and isopentane appear here because GERG-2008 moved their reducing
/// parameters, which moves the whole saturation curve.
inline const std::map<std::string, AncillarySet>& ancillaries_2008_overrides() {{
    static const std::map<std::string, AncillarySet> data = {{
{rows_2008}
    }};
    return data;
}}

}}  // namespace detail
}}  // namespace GERG
}}  // namespace CoolProp
#endif /* GERGANCILLARIES_H_ */
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", action="store_true", help="print the fit report instead of the header")
    args = ap.parse_args()

    check_struct_field_order()

    fits_2004 = {}
    fits_2008 = {}
    report = []

    for name in GERG2004_NAMES:
        f = FluidFit(2004, name, strict=not args.report)
        fits_2004[name] = f
        report.append(f)
    for name in CHANGED_IN_2008 + NEW_IN_2008:
        f = FluidFit(2008, name, strict=not args.report)
        fits_2008[name] = f
        report.append(f)

    if args.report:
        print(
            f"{'model/fluid':26s} {'Tmin':>9s} {'T_r':>9s} {'npts':>5s} "
            f"{'dev rhoL':>10s} {'dev rhoV':>10s} {'dev pV':>10s}  "
            f"{'dTc(true-red)':>13s} {'drhoc':>9s}"
        )
        for f in report:
            print(
                f"{f.year}/{f.name:20s} {f.Tmin:9.3f} {f.Tmax:9.3f} {len(f.fit_rows):5d} "
                f"{f.dev_rhoL:10.3e} {f.dev_rhoV:10.3e} {f.dev_pV:10.3e}  "
                f"{f.P.Tc_true - f.P.T_reduce:+13.4f} "
                f"{(f.P.rhoc_true / f.P.rho_reduce - 1) * 100:+8.4f}%"
            )
        worst = max(r.worst for r in report)
        print(f"\nworst relative deviation across all {len(report)} pure EOS: {worst:.3e}")
        return

    rows_2004 = ",\n".join("      " + emit_row(fits_2004[n]) for n in GERG2004_NAMES)
    rows_2008 = ",\n".join("      " + emit_row(fits_2008[n]) for n in CHANGED_IN_2008 + NEW_IN_2008)
    sys.stdout.write(
        HEADER.format(
            teqp_version=teqp.__version__,
            worst_all=max(r.worst for r in report),
            worst_rhoL=max(r.dev_rhoL for r in report),
            worst_rhoV=max(r.dev_rhoV for r in report),
            worst_pV=max(r.dev_pV for r in report),
            rows_2004=rows_2004,
            rows_2008=rows_2008,
        )
    )


if __name__ == "__main__":
    main()
