"""Coverage for PropertyPlot.draw_process.

draw_process used to warn on every call that it was "not tested".  This
runs the compression cycle from the docs page (Web/coolprop/python-cycles.rst)
through it and checks that it draws the process line and the state markers
at the cycle's own states, and that it no longer warns.  ``CoolProp.Plots``
imports matplotlib, so everything here skips without it.
"""
import warnings

import numpy as np
import pytest

import CoolProp


@pytest.fixture
def cycle_plot():
    """(PropertyPlot, StateContainer) for the R134a compression cycle in the docs"""
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from CoolProp.Plots import PropertyPlot, SimpleCompressionCycle

    pp = PropertyPlot('HEOS::R134a', 'PH', unit_system='EUR')
    cycle = SimpleCompressionCycle('HEOS::R134a', 'PH', unit_system='EUR')
    T0, T2 = 280, 310
    pp.state.update(CoolProp.QT_INPUTS, 0.0, T0 - 10)
    p0 = pp.state.keyed_output(CoolProp.iP)
    pp.state.update(CoolProp.QT_INPUTS, 1.0, T2 + 15)
    p2 = pp.state.keyed_output(CoolProp.iP)
    cycle.simple_solve(T0, p0, T2, p2, 0.7, SI=True)
    cycle.steps = 50
    yield pp, cycle.get_state_changes()
    plt.close(cycle.figure)
    plt.close(pp.figure)


def test_draw_process_draws_cycle_without_warning(cycle_plot):
    pp, sc = cycle_plot
    before = len(pp.axis.lines)

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        pp.draw_process(sc)

    # One line for the process, one for the markers of the repeated points
    assert len(pp.axis.lines) == before + 2
    process, markers = pp.axis.lines[before:]

    # The process line runs through the states that differ from the one before
    # them (the first compared with the last, as the cycle closes), the markers
    # sit on the repeated ones, all in the axis units
    moved, repeated = [], []
    old = sc[len(sc) - 1]
    for i in sc:
        (repeated if sc[i] == old else moved).append(sc[i])
        old = sc[i]
    assert moved and repeated
    dimx = pp.system[pp.x_index]
    dimy = pp.system[pp.y_index]
    for line, states in ((process, moved), (markers, repeated)):
        x = dimx.from_SI(np.array([s[pp.x_index] for s in states]))
        y = dimy.from_SI(np.array([s[pp.y_index] for s in states]))
        np.testing.assert_allclose(line.get_xdata(), x)
        np.testing.assert_allclose(line.get_ydata(), y)
    assert process.get_linestyle() == '-'
    assert markers.get_marker() == 'o'


def test_draw_process_accepts_line_options(cycle_plot):
    pp, sc = cycle_plot
    before = len(pp.axis.lines)
    pp.draw_process(sc, line_opts={'color': 'blue', 'lw': 2.5})
    process = pp.axis.lines[before]
    assert process.get_color() == 'blue'
    assert process.get_linewidth() == 2.5
