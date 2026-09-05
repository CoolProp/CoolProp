#!/usr/bin/env python3
"""Compute the acentric factor of every GERG pure EOS from the GERG EOS itself.

WHY THIS EXISTS
---------------
GERG-2004/GERG-2008 publish no acentric factor.  ``make_gerg_fluid`` used to
leave ``EquationOfState::acentric`` at the ``_HUGE`` sentinel, and that
sentinel is read DIRECTLY -- not through the overridable
``calc_acentric_factor`` accessor -- by

  * ``SaturationSolvers::x_and_y_from_K`` / the Wilson K-factor seed, which is
    every mixture QT/PQ/DQ flash and ``build_phase_envelope()``;
  * ``FlashRoutines::T_DP_PengRobinson``, i.e. every ``DmolarP`` flash that
    lands in the gas phase;
  * ``HelmholtzEOSMixtureBackend::solver_rho_Tp_SRK``, which is what
    ``solver_rho_Tp`` seeds from -- the pure-fluid ``HmolarP`` / ``PSmolar`` /
    ``PUmolar`` failures.

Borrowing CoolProp's tabulated acentric factor for the same fluid would work
numerically and would break the one property this backend exists to have:
nothing in it traces back to a non-GERG equation.  CoolProp's value belongs to
that fluid's REFERENCE equation of state (Setzmann-Wagner for methane,
IAPWS-95 for water, ...), not to GERG's shortened technical form.

So it is DERIVED from the GERG EOS, using Pitzer's defining relation

    omega = -1 - log10( p_sat(0.7 * T_c) / p_c )

which needs nothing but a saturation solve and a critical pressure -- both of
which the GERG pure EOS supplies on its own.

TWO JUDGEMENT CALLS, both deliberate
------------------------------------
1.  ``p_sat(0.7*T_c)`` comes from a CONVERGED saturation solve on the GERG
    EOS, NOT from the fitted ``pV`` ancillary in GERGAncillaries.h.  The
    ancillary is a SEED, never an answer (see fit_ancillaries.py's docstring
    and the Task 11 test "GERG saturation is independent of the ancillary
    seed"); making a shipped constant depend on it would make omega
    fit-dependent, and a refit with a different term ladder would silently
    move it.  The converged solve is GERG-exact instead.  Measured cost of the
    choice: the ancillary would move omega by at most 9.1e-05 (nitrogen),
    which is small -- the point is provenance, not magnitude.

    Reliability was checked, not assumed: the guess-free solve
    (``sat_from_scratch``: spinodal scan, geometric bisection on p with
    equal-chemical-potential as the residual, then a damped 2-D Newton
    polish) converges at 0.7*T_c for ALL 23 pure EOS, with a final
    |p'-p''|/p'' residual of at most 5.5e-12.  No fluid needs an ancillary
    fallback and none is implemented -- a silent fallback is exactly the kind
    of fail-open this backend refuses elsewhere.  If a future coefficient
    change makes some fluid fail to converge, this script ABORTS.

2.  ``p_c`` is the pressure the EOS produces at GERG's TABULATED reducing
    state (Table A3.5), i.e. exactly the number ``make_gerg_fluid`` stores in
    ``EOS.reduce.p``/``fluid.crit.p`` and exactly what the backend answers
    from ``p_critical()``.  It is NOT the pressure at the true critical point
    of the shortened form, which can sit up to +1.10 K / -5.4% away.

    The reason is CONSISTENCY WITH THE CONSUMERS, not convenience.  The
    Wilson correlation the seed uses is

        ln K_i = ln(p_c,i / p) + 5.373 (1 + omega_i) (1 - T_c,i / T)

    and it reads ``p_c,i`` and ``T_c,i`` from the SAME fluid object that
    carries omega.  Using the tabulated reducing point in both places makes
    the correlation reproduce p_sat/p exactly at T = 0.7 T_c; mixing a
    true-critical p_c into omega while the correlation keeps reading the
    reducing-point p_c would put the inconsistency straight into the seed.
    The same argument applies to ``solver_rho_Tp_SRK``, which builds a and b
    from ``EOS.reduce.T``/``EOS.reduce.p`` and alpha from ``EOS.acentric``.
    ``0.7 * T_c`` likewise uses the tabulated reducing temperature.

    Measured size of the choice: using the true critical pressure instead
    would move omega by up to 6.7e-03 (n-heptane).  ``EOS.reduce.p`` is finite
    and positive for all 23 (0.227 MPa for helium to 22.06 MPa for water);
    this script asserts that.

HELIUM AND HYDROGEN
-------------------
``make_gerg_fluid`` caps ``EOS.limits.Tmin`` at the component's own reducing
temperature, so for helium (T_c = 5.1953 K) and hydrogen (T_c = 33.19 K) the
enforced Tmin EQUALS T_c and 0.7*T_c (3.6367 K / 23.2330 K) is below the range
``check_gerg_range_of_validity`` lets a caller reach at run time.

That is not an obstacle to this quantity and it is not smuggled past one.  The
acentric factor is a DEFINITIONAL, compile-time constant of the equation, not
a state a caller asks for: it is evaluated here, once, offline, and only the
resulting scalar is shipped.  The GERG helium and hydrogen equations are
perfectly well behaved at 0.7*T_c -- the same VLE solver traces both all the
way down to 0.30*T_c when fitting their ancillaries -- so there is nothing
numerically delicate about the evaluation either.  What the run-time Tmin
encodes is the published operating envelope of the MIXTURE model, and
evaluating an equation below its recommended range to obtain a defining
constant of that equation is a different act from answering a user's query
there.  Both values are therefore computed exactly like the other 21, and the
consequence is stated in Web/coolprop/GERG.rst rather than left implicit.

WHAT IT PRODUCES
----------------
src/Backends/GERG/GERGAcentric.h -- a generated C++ table consumed by
``GERG::get_acentric_factor()`` (declared in GERGData.h, defined in
GERGBackend.cpp) and written into ``EOS.acentric`` by ``make_gerg_fluid``.

23 rows, not 21: GERG-2008 moves the reducing parameters of carbon monoxide
and isopentane, which moves their whole saturation curve, so each model needs
its own value.  The table is a GERG-2004 base plus a GERG-2008 override set,
mirroring pure_info_2004()/pure_info_2008_overrides() in GERGData.h.

RUN (teqp 0.23.2; the system python3's teqp 0.15.3 predates the GERG factory
registration and fails with "Unknown kind:GERG2008resid"):

    python3 -m venv /tmp/gergenv
    /tmp/gergenv/bin/pip install "teqp>=0.18" numpy
    /tmp/gergenv/bin/python dev/gerg/compute_acentric.py \
        > src/Backends/GERG/GERGAcentric.h
    uvx clang-format@18.1.8 -i src/Backends/GERG/GERGAcentric.h

    # human-readable report, including the comparison against published
    # acentric factors and against the two rejected alternatives:
    /tmp/gergenv/bin/python dev/gerg/compute_acentric.py --report

RECOMPUTE WHEN THE COEFFICIENT TABLES CHANGE.  These values are tied to the
GERG coefficient tables in GERGData.h/GERGBackend.cpp exactly as the
ancillaries are; editing a pure-fluid n/t/d/l row or a Table A3.5 reducing
parameter invalidates them.  dev/gerg/verify_transcription.py re-derives the
whole table and diffs it against the committed header, so a stale table is
caught mechanically rather than by remembering to rerun this.

REPRODUCIBILITY.  Unlike the ancillary fits (a least-squares solve, whose last
two or three digits move with the numpy version), this is a bisection followed
by a Newton polish to a fixed residual, so it reproduces to ~1e-12 across
environments.  verify_transcription.py compares with a 1e-9 ABSOLUTE tolerance
(|committed - re-derived|, verify_transcription.py:842) to leave room for that
without letting a real change through.
"""

import argparse
import math
import sys
import warnings

warnings.filterwarnings("ignore")  # teqp 0.23 emits FutureWarnings for top-level fns

import teqp  # noqa: E402

from fit_ancillaries import (  # noqa: E402
    CHANGED_IN_2008,
    GERG2004_NAMES,
    NEW_IN_2008,
    Pure,
    newton_polish,
    sat_from_scratch,
)

#: T/T_c at which the saturation pressure is evaluated -- Pitzer's definition.
OMEGA_T_FRACTION = 0.7

#: Largest tolerated |p' - p''| / p'' at the converged saturation state.  The
#: measured worst case across all 23 pure EOS is 5.5e-12; this is a factor of
#: ~200 of headroom, and a solve that misses it ABORTS rather than shipping a
#: half-converged number.
MAX_SAT_RESIDUAL = 1e-9

#: Published acentric factors, for the sanity column of --report ONLY.  These
#: are NEVER written into the generated header and nothing in the backend
#: reads them: they belong to each fluid's reference equation of state, which
#: is precisely the provenance this backend refuses.  They are here so that a
#: gross error (wrong sign, wrong order of magnitude, a fluid mixed up with
#: another) is visible in the report.  Source: Poling, Prausnitz & O'Connell,
#: "The Properties of Gases and Liquids", 5th ed., Appendix A.
PUBLISHED_OMEGA = {
    "methane": 0.011,
    "nitrogen": 0.037,
    "carbondioxide": 0.224,
    "ethane": 0.099,
    "propane": 0.152,
    "n-butane": 0.200,
    "isobutane": 0.186,
    "n-pentane": 0.251,
    "isopentane": 0.227,
    "n-hexane": 0.299,
    "n-heptane": 0.349,
    "n-octane": 0.398,
    "hydrogen": -0.219,
    "oxygen": 0.022,
    "carbonmonoxide": 0.050,
    "water": 0.344,
    "helium": -0.385,
    "argon": -0.002,
    "hydrogensulfide": 0.100,
    "n-nonane": 0.443,
    "n-decane": 0.489,
}


class Acentric:
    """omega for one (model, component), plus everything --report shows."""

    def __init__(self, year, name):
        self.year = year
        self.name = name
        P = Pure(year, name)
        self.P = P

        # The TABULATED reducing state, which is what the backend reports as
        # its critical point -- see judgement call 2 in the module docstring.
        self.Tc = P.T_reduce
        self.pc = P.p(P.T_reduce, P.rho_reduce)
        if not (math.isfinite(self.pc) and self.pc > 0.0):
            raise RuntimeError(f"{year}/{name}: p at the reducing state is {self.pc}, not a usable p_c")

        self.T = OMEGA_T_FRACTION * self.Tc

        # Guess-free solve, then polish.  Deliberately no ancillary fallback:
        # see judgement call 1.
        sol = sat_from_scratch(P, self.T)
        if sol is None:
            raise RuntimeError(f"{year}/{name}: no saturation state found at T = {self.T} K ({OMEGA_T_FRACTION}*T_c)")
        polished = newton_polish(P, self.T, *sol)
        if polished is None:
            raise RuntimeError(f"{year}/{name}: the Newton polish diverged at T = {self.T} K")
        self.rhoL, self.rhoV = polished

        pL, pV = P.p(self.T, self.rhoL), P.p(self.T, self.rhoV)
        self.residual = abs(pL - pV) / abs(pV)
        if not (self.residual <= MAX_SAT_RESIDUAL):
            raise RuntimeError(
                f"{year}/{name}: saturation residual {self.residual:.3e} exceeds {MAX_SAT_RESIDUAL:.3e} at T = {self.T} K"
            )
        # p'' is the vapour-side pressure; the two agree to `residual`.
        self.psat = pV
        if not (math.isfinite(self.psat) and self.psat > 0.0):
            raise RuntimeError(f"{year}/{name}: p_sat is {self.psat}")

        self.omega = -1.0 - math.log10(self.psat / self.pc)
        if not math.isfinite(self.omega):
            raise RuntimeError(f"{year}/{name}: omega is {self.omega}")

    # -- report-only alternatives, never emitted -------------------------
    def omega_true_critical(self):
        """omega with p_c taken at the TRUE critical point of the short form
        instead of the tabulated reducing state.  Rejected: see judgement
        call 2."""
        pc_true = self.P.p(self.P.Tc_true, self.P.rhoc_true)
        return -1.0 - math.log10(self.psat / pc_true)


def compute_all():
    """Every (model, component) whose pure EOS is distinct, in the order the
    generated header lists them."""
    rows = []
    for name in GERG2004_NAMES:
        rows.append(Acentric(2004, name))
    for name in CHANGED_IN_2008 + NEW_IN_2008:
        rows.append(Acentric(2008, name))
    return rows


def _fmt(x):
    # float() first: repr(np.float64(...)) is "np.float64(0.011)" on numpy 2,
    # which is not a C++ literal.
    return repr(float(x))


HEADER = '''// GENERATED FILE -- do not edit by hand.
//
// Acentric factors DERIVED FROM THE GERG PURE EQUATIONS OF STATE, produced by
// dev/gerg/compute_acentric.py from teqp {teqp_version}
// (https://github.com/usnistgov/teqp).  CoolProp's build has no dependency on
// teqp; this header is committed, exactly like GERGAncillaries.h and
// GERGReferenceValues.h.
//
// GERG-2004/GERG-2008 publish no acentric factor.  These are NOT CoolProp's
// tabulated acentric factors for the same fluids -- those belong to each
// fluid's REFERENCE equation of state (Setzmann-Wagner for methane, IAPWS-95
// for water, ...), and a backend named GERG-2004/GERG-2008 must not carry
// data traceable to a different equation.  Each value here is Pitzer's
// defining relation evaluated on GERG's own shortened technical form:
//
//     omega = -1 - log10( p_sat(0.7 * T_c) / p_c )
//
//   * p_sat(0.7*T_c) is a CONVERGED saturation solve on the GERG pure EOS
//     (worst |p' - p''|/p'' across all 23: {worst_residual:.2e}), NOT the fitted pV
//     ancillary -- an ancillary is a seed, never an answer.
//   * p_c is the pressure the EOS produces at GERG's TABULATED reducing state
//     (Table A3.5), i.e. the same number make_gerg_fluid stores in
//     EOS.reduce.p and the backend reports from p_critical(); T_c likewise.
//     That consistency is load-bearing -- the Wilson K-factor seed reads
//     p_c/T_c and omega off the same fluid.
//
// See dev/gerg/compute_acentric.py for both judgement calls in full, and for
// why helium and hydrogen -- whose enforced Tmin equals their T_c, so 0.7*T_c
// is below the range a caller can reach -- are computed the same way as the
// other 21 rather than left out.
//
// RECOMPUTE WHENEVER A GERG COEFFICIENT TABLE CHANGES.
// dev/gerg/verify_transcription.py re-derives this whole table and diffs it
// against these literals, so a stale table is caught mechanically.
//
// 23 rows, not 21: GERG-2008 moves the reducing parameters of carbon monoxide
// and isopentane, which moves their saturation curve, so each model has its
// own value.
#ifndef GERGACENTRIC_H_
#define GERGACENTRIC_H_

#include <map>
#include <string>

#include "GERGData.h"

namespace CoolProp {{
namespace GERG {{
namespace detail {{

/// Keyed on the GERG component name.  GERG-2004 is the base table and
/// GERG-2008 overrides it, mirroring pure_info_2004/pure_info_2008_overrides
/// in GERGData.h.
inline const std::map<std::string, double>& acentric_2004() {{
    static const std::map<std::string, double> data = {{
{rows_2004}
    }};
    return data;
}}

/// Entries GERG-2008 changes or adds relative to GERG-2004.
inline const std::map<std::string, double>& acentric_2008_overrides() {{
    static const std::map<std::string, double> data = {{
{rows_2008}
    }};
    return data;
}}

}}  // namespace detail
}}  // namespace GERG
}}  // namespace CoolProp
#endif /* GERGACENTRIC_H_ */
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", action="store_true", help="print the human-readable report instead of the header")
    args = ap.parse_args()

    rows = compute_all()

    if args.report:
        print(
            f"{'model/fluid':27s} {'T_c/K':>10s} {'p_c/MPa':>9s} {'0.7T_c/K':>9s} "
            f"{'p_sat/Pa':>12s} {'omega':>9s} {'published':>10s} {'diff':>9s} "
            f"{'omega(pc_true)':>15s} {'sat resid':>10s}"
        )
        for r in rows:
            pub = PUBLISHED_OMEGA.get(r.name, float("nan"))
            print(
                f"{r.year}/{r.name:21s} {r.Tc:10.4f} {r.pc / 1e6:9.5f} {r.T:9.4f} "
                f"{r.psat:12.6g} {r.omega:9.5f} {pub:10.3f} {r.omega - pub:+9.4f} "
                f"{r.omega_true_critical():15.5f} {r.residual:10.2e}"
            )
        worst_pub = max(abs(r.omega - PUBLISHED_OMEGA[r.name]) for r in rows)
        worst_true = max(abs(r.omega_true_critical() - r.omega) for r in rows)
        print(f"\n{len(rows)} pure EOS")
        print(f"largest |omega - published|                 : {worst_pub:.3e}")
        print(f"largest |omega(true p_c) - omega(reducing)| : {worst_true:.3e}  (rejected alternative)")
        print(f"worst saturation residual |p'-p''|/p''      : {max(r.residual for r in rows):.3e}")
        return

    by_key = {(r.year, r.name): r for r in rows}
    rows_2004 = ",\n".join(f'      {{"{n}", {_fmt(by_key[(2004, n)].omega)}}}' for n in GERG2004_NAMES)
    rows_2008 = ",\n".join(f'      {{"{n}", {_fmt(by_key[(2008, n)].omega)}}}' for n in CHANGED_IN_2008 + NEW_IN_2008)
    sys.stdout.write(
        HEADER.format(
            teqp_version=teqp.__version__,
            worst_residual=max(r.residual for r in rows),
            rows_2004=rows_2004,
            rows_2008=rows_2008,
        )
    )


if __name__ == "__main__":
    main()
