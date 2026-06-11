"""Superancillary bindings (CoolProp-1tbe.5).

Locks the nanobind ports of the Chebyshev rootfinding building blocks
(``ChebyshevExpansion``, ``ChebyshevApproximation1D``) and the ``SuperAncillary``
saturation evaluator, mirroring the legacy Cython surface.  These are advertised
in the v8 changelog and exercised by the SuperAncillary docs notebook
(``CP.SuperAncillary(json.dumps(...))`` + ``eval_sat_many``), so a missing or
broken binding must not ship.

Runs against an *installed* CoolProp (legacy or nanobind); both expose the same
classes.
"""
import json

import numpy as np
import pytest

from CoolProp import CoolProp as CP

pytestmark = pytest.mark.skipif(
    not hasattr(CP, "SuperAncillary"), reason="superancillary classes not bound in this build"
)


def _superanc(fluid):
    jEOS = json.loads(CP.get_fluid_param_string(fluid, "JSON"))[0]["EOS"][0]
    if "SUPERANCILLARY" not in jEOS:
        pytest.skip(f"{fluid} has no embedded SUPERANCILLARY")
    return CP.SuperAncillary(json.dumps(jEOS["SUPERANCILLARY"]))


@pytest.mark.parametrize("fluid", ["Propane", "Water"])
def test_superancillary_eval_sat_many_matches_eos(fluid):
    # The EOS itself uses superancillaries internally, so the SA evaluator must
    # agree with PropsSI to (near) machine precision across the saturation line.
    sa = _superanc(fluid)
    Tt = CP.PropsSI(fluid, "Ttriple")
    Tc = CP.PropsSI(fluid, "Tcrit")
    T = np.linspace(Tt + 5.0, Tc - 5.0, 12)

    rhoL = np.zeros_like(T)
    rhoV = np.zeros_like(T)
    p = np.zeros_like(T)
    sa.eval_sat_many(T, "D", 0, rhoL)
    sa.eval_sat_many(T, "D", 1, rhoV)
    sa.eval_sat_many(T, "P", 1, p)

    for i, Ti in enumerate(T):
        assert rhoL[i] == pytest.approx(CP.PropsSI("Dmolar", "T", Ti, "Q", 0, "HEOS::" + fluid), rel=1e-6)
        assert rhoV[i] == pytest.approx(CP.PropsSI("Dmolar", "T", Ti, "Q", 1, "HEOS::" + fluid), rel=1e-6)
        assert p[i] == pytest.approx(CP.PropsSI("P", "T", Ti, "Q", 1, "HEOS::" + fluid), rel=1e-6)


def test_superancillary_eval_sat_scalar_and_errors():
    sa = _superanc("Propane")
    T = 0.5 * (CP.PropsSI("Propane", "Ttriple") + CP.PropsSI("Propane", "Tcrit"))
    assert sa.eval_sat(T, "P", 1) == pytest.approx(CP.PropsSI("P", "T", T, "Q", 1, "HEOS::Propane"), rel=1e-6)
    # Empty property string is rejected
    with pytest.raises(ValueError):
        sa.eval_sat(T, "", 1)
    # Mismatched array lengths are rejected
    with pytest.raises(ValueError):
        sa.eval_sat_many(np.array([T, T]), "P", 1, np.zeros(1))


def test_chebyshev_expansion_eval_and_invert():
    from numpy.polynomial import chebyshev as C

    xmin, xmax = 1.0, 3.0

    def f(x):  # monotonic increasing on [1, 3]
        return x ** 3 + 2.0 * x

    xs = np.linspace(xmin, xmax, 200)
    coef = list(C.Chebyshev.fit(xs, f(xs), 12, domain=[xmin, xmax]).coef)
    e = CP.ChebyshevExpansion(xmin, xmax, coef)

    assert e.xmin() == xmin and e.xmax() == xmax
    x = np.linspace(xmin, xmax, 9)
    y = np.zeros_like(x)
    e.eval_many(x, y)
    assert np.allclose(y, f(x), rtol=1e-9)

    # solve_for_x inverts eval (f(2) -> x == 2)
    assert e.solve_for_x(f(2.0), xmin, xmax, 53, 100, 1e-12) == pytest.approx(2.0, abs=1e-8)

    # vectorized inversion + step counts
    yq = f(np.array([1.5, 2.0, 2.5]))
    xout = np.zeros_like(yq)
    counts = np.zeros(yq.size, dtype=np.uintp)
    e.solve_for_x_many(yq, xmin, xmax, 53, 100, 1e-12, xout, counts)
    assert np.allclose(xout, [1.5, 2.0, 2.5], atol=1e-8)
    assert np.all(counts > 0)


def test_chebyshev_approximation1d():
    from numpy.polynomial import chebyshev as C

    xmin, xmax = 1.0, 3.0

    def f(x):
        return x ** 3 + 2.0 * x

    xs = np.linspace(xmin, xmax, 200)
    e = CP.ChebyshevExpansion(xmin, xmax, list(C.Chebyshev.fit(xs, f(xs), 12, domain=[xmin, xmax]).coef))
    a = CP.ChebyshevApproximation1D([e])

    assert a.xmin() == xmin and a.xmax() == xmax
    assert a.is_monotonic() is True

    x = np.linspace(xmin, xmax, 9)
    y = np.zeros_like(x)
    a.eval_many(x, y)
    assert np.allclose(y, f(x), rtol=1e-9)

    sols = a.get_x_for_y(f(2.0), 53, 100, 1e-12)
    assert len(sols) == 1
    assert sols[0][0] == pytest.approx(2.0, abs=1e-8)

    counts = np.zeros(3, dtype=np.uintp)
    a.count_x_for_y_many(f(np.array([1.5, 2.0, 2.5])), 53, 100, 1e-12, counts)
    assert list(counts) == [1, 1, 1]

    intervals = a.monotonic_intervals()
    assert len(intervals) == 1
    assert intervals[0].xmin == pytest.approx(xmin)
    assert intervals[0].xmax == pytest.approx(xmax)
