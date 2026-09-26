"""Coverage for IsoLine.calc_sat_range's near-critical rescue (COO-12).

When a saturation flash fails within 1 K (QT) or 100 Pa (PQ) of the critical
point, calc_sat_range substitutes the critical point for the failed point.
That branch used to end in ``pass`` instead of ``continue``, so the
substitution was immediately overwritten with NaN and never had any effect:
quality lines stopped short of the critical point with a gap.

The flash failure is injected with a thin proxy around a real state, so the
test does not depend on which points a given solver version happens to fail
at.  ``CoolProp.Plots.Common`` imports matplotlib, so everything here skips
without it.
"""
import warnings

import numpy as np
import pytest

import CoolProp


@pytest.fixture(scope="module")
def Common():
    """CoolProp.Plots.Common, or skip if matplotlib is absent"""
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from CoolProp.Plots import Common
    return Common


class _FailingState(object):
    """Delegates to a real state, but its update() raises when told to"""

    def __init__(self, state, should_fail):
        self._wrapped = state
        self._should_fail = should_fail

    def update(self, pair, one, two):
        if self._should_fail(pair, one, two):
            raise ValueError("injected saturation-flash failure")
        return self._wrapped.update(pair, one, two)

    def __getattr__(self, name):
        return getattr(self._wrapped, name)


def _isoline(Common, Q):
    iso = Common.IsoLine(CoolProp.iQ, CoolProp.iHmass, CoolProp.iP, value=Q,
                         state=CoolProp.AbstractState("HEOS", "Propane"))
    crit = iso.critical_state
    return iso, crit.keyed_output(CoolProp.iT), crit.keyed_output(CoolProp.iP)


@pytest.mark.parametrize("Q", [0.0, 1.0])
def test_failure_near_critical_temperature_is_replaced_by_critical_point(Common, Q):
    iso, Tc, _ = _isoline(Common, Q)
    # Last point 0.5 K below Tc: inside the 1 K rescue window.  One point
    # 5 K below Tc also fails, and is outside it, so must still be NaN.
    Trange = np.array([Tc - 20.0, Tc - 5.0, Tc - 0.5])
    iso._state = _FailingState(iso.state, lambda pair, one, two: two > Tc - 6.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        iso.calc_sat_range(Trange=Trange)
    xc = iso.critical_state.keyed_output(CoolProp.iHmass)
    yc = iso.critical_state.keyed_output(CoolProp.iP)
    assert np.isfinite(iso.x[0]) and np.isfinite(iso.y[0])
    assert np.isnan(iso.x[1]) and np.isnan(iso.y[1])
    assert iso.x[2] == xc and iso.y[2] == yc


def test_failure_near_critical_pressure_is_replaced_by_critical_point(Common):
    iso, _, pc = _isoline(Common, 0.0)
    Prange = np.array([0.5 * pc, pc - 50.0])  # 50 Pa below pc: inside the 100 Pa window
    iso._state = _FailingState(iso.state, lambda pair, one, two: one > 0.9 * pc)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        iso.calc_sat_range(Prange=Prange)
    assert np.isfinite(iso.x[0])
    assert iso.x[1] == iso.critical_state.keyed_output(CoolProp.iHmass)
    assert iso.y[1] == iso.critical_state.keyed_output(CoolProp.iP)
