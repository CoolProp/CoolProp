"""Coverage for CoolProp.Plots.Common.IsoLineTracer (bd CoolProp-dgxi).

The tracer walks along an isoline instead of flashing every point cold, which
is what makes a mixture property plot finish in seconds rather than minutes
(GH discussion #3269).  What these tests pin down is that the speed does not
cost correctness, and in particular the three ways it could:

* the traced point has to be the same state the flash would have found;
* the tracer has to stay on the right side of the saturation curve -- which
  for a near-azeotropic blend cannot be decided from the bubble and dew values
  of the swept property, because they agree to within round-off;
* an isoline the tracer cannot bracket has to be flashed, not traced.  A
  phase-imposed Newton iteration converges happily on the metastable extension
  of the EOS inside the two-phase dome, and the answer looks entirely
  self-consistent, so nothing downstream would notice.

``CoolProp.Plots.Common`` imports matplotlib, so everything here skips without
it.
"""
import pytest

import CoolProp

AbstractState = CoolProp.AbstractState


@pytest.fixture(scope="module")
def Common():
    """The worktree's CoolProp.Plots.Common, or skip if matplotlib is absent"""
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    from CoolProp.Plots import Common
    return Common


def _tracer(Common, backend, fluid, index1, index2, iso_index):
    """Build a tracer on a fresh state of this fluid"""
    return Common.IsoLineTracer(AbstractState(backend, fluid), index1, index2, iso_index)


def _reconstruct(backend, fluid, state, mole_fractions=None):
    """Rebuild a traced state from (T, Q) or (T, rho) on a fresh object

    Neither route uses a two-dimensional flash, so this checks the solvers the
    tracer used to get there.  The single-phase route imposes the phase
    because the free ``DmassT`` path on a mixture can return a state at a
    different density than it was given (bd CoolProp-1gth).

    Note what this deliberately does *not* prove: re-evaluating the same
    ``(T, rho)`` reproduces a metastable root just as faithfully as a stable
    one, so a root that should have been two-phase reconstructs perfectly.
    Checking for that needs the saturation densities at the point's own
    temperature -- see TestBranchSelection and TestUnbracketableInputs.
    """
    fresh = AbstractState(backend, fluid)
    if mole_fractions is not None:
        fresh.set_mole_fractions(mole_fractions)
    Q = state.Q()
    if 0.0 <= Q <= 1.0:
        fresh.update(CoolProp.QT_INPUTS, Q, state.T())
    else:
        fresh.specify_phase(CoolProp.iphase_gas)
        fresh.update(CoolProp.DmassT_INPUTS, state.rhomass(), state.T())
    return fresh


def _relerr(a, b):
    """Relative difference, safe when either value is near zero"""
    return abs(a - b) / max(abs(a), abs(b), 1e-8)


def _isoline(Common, fluid, graph, iso_type, value, points, tracing=True):
    """Calculate one isoline the way PropertyPlot.calc_isolines does"""
    from CoolProp.Plots import PropertyPlot
    import numpy as np
    plot = PropertyPlot(fluid, graph, unit_system='SI', tp_limits='ACHP')
    limits = plot._get_axis_limits()
    xvals = plot.generate_ranges(plot._x_index, limits[0], limits[1], points)
    yvals = plot.generate_ranges(plot._y_index, limits[2], limits[3], points)
    line = Common.IsoLine(iso_type, plot._x_index, plot._y_index,
                          value=value, state=plot._state, tracing=tracing)
    line.calc_range(xvals, yvals)
    return np.asarray(line.x), np.asarray(line.y)


class TestSinglePhase:
    def test_isentrope_matches_the_flash(self, Common):
        """A traced superheated isentrope reproduces the PSmass flash"""
        smass = 7000.0  # J/kg/K, superheated steam
        tracer = _tracer(Common, "HEOS", "Water", CoolProp.iP, CoolProp.iSmass, CoolProp.iSmass)
        flash = AbstractState("HEOS", "Water")
        for p in [1e5, 3e5, 1e6, 3e6, 1e7]:
            tracer.update(p, smass)
            flash.update(CoolProp.PSmass_INPUTS, p, smass)
            assert _relerr(tracer.keyed_output(CoolProp.iHmass), flash.hmass()) < 1e-9
            assert _relerr(tracer.keyed_output(CoolProp.iT), flash.T()) < 1e-9

    def test_traced_state_reconstructs(self, Common):
        """Every traced point rebuilds to itself from (T, rho) alone"""
        tracer = _tracer(Common, "HEOS", "R134a", CoolProp.iP, CoolProp.iSmass, CoolProp.iSmass)
        smass = 1900.0
        for p in [5e4, 1e5, 2e5, 5e5, 1e6]:
            tracer.update(p, smass)
            fresh = _reconstruct("HEOS", "R134a", tracer.state)
            assert _relerr(fresh.smass(), smass) < 1e-8
            assert _relerr(fresh.p(), p) < 1e-8


class TestTwoPhase:
    def test_quality_solve_on_a_wide_glide_mixture(self, Common):
        """Inside the dome the tracer solves for the quality, not a density root"""
        fluid = "Methane&Ethane"
        state = AbstractState("HEOS", fluid)
        state.set_mole_fractions([0.6, 0.4])
        sat = AbstractState("HEOS", fluid)
        sat.set_mole_fractions([0.6, 0.4])
        p = 2e6
        sat.update(CoolProp.PQ_INPUTS, p, 0.0)
        s_liq = sat.smass()
        sat.update(CoolProp.PQ_INPUTS, p, 1.0)
        s_vap = sat.smass()
        assert s_vap - s_liq > 1.0, "the test needs a real two-phase span"

        tracer = Common.IsoLineTracer(state, CoolProp.iP, CoolProp.iSmass, CoolProp.iP)
        for frac in [0.2, 0.5, 0.8]:
            smass = s_liq + frac * (s_vap - s_liq)
            tracer.update(p, smass)
            assert 0.0 < tracer.state.Q() < 1.0
            assert _relerr(tracer.keyed_output(CoolProp.iSmass), smass) < 1e-8
            assert _relerr(tracer.keyed_output(CoolProp.iP), p) < 1e-8
            fresh = _reconstruct("HEOS", fluid, tracer.state, [0.6, 0.4])
            assert _relerr(fresh.smass(), smass) < 1e-8
            assert _relerr(fresh.p(), p) < 1e-8


class TestBranchSelection:
    def test_near_azeotrope_isotherm_crosses_to_the_liquid(self, Common):
        """R513A brackets p by psat on both sides; the branch must not be guessed

        Near 305 K the bubble and dew pressures of this blend come within
        1.4e-9 relative of each other, which is inside BRACKET_RTOL, so the
        bracket reads as degenerate and its two values cannot say which side
        of the saturation curve a point is on.  The branch has to come from
        the slope of p along density instead.  Pick the temperature carefully:
        a few degrees away the gap is 1e-6 or wider and the value comparison
        would carry the test on its own.
        """
        T = 305.0
        state = AbstractState("HEOS", "R513A.mix")
        sat = AbstractState("HEOS", "R513A.mix")
        sat.update(CoolProp.QT_INPUTS, 0.0, T)
        psat = sat.p()
        tracer = Common.IsoLineTracer(state, CoolProp.iP, CoolProp.iT, CoolProp.iT)
        flash = AbstractState("HEOS", "R513A.mix")
        for p, expect_liquid in [(0.5 * psat, False), (0.9 * psat, False),
                                 (1.1 * psat, True), (2.0 * psat, True)]:
            tracer.update(p, T)
            rho = tracer.state.rhomass()
            sat.update(CoolProp.QT_INPUTS, 1.0 if not expect_liquid else 0.0, T)
            if expect_liquid:
                assert rho > sat.rhomass() * (1.0 - 1e-6), "expected the liquid branch at p={0}".format(p)
            else:
                assert rho < sat.rhomass() * (1.0 + 1e-6), "expected the vapour branch at p={0}".format(p)
            flash.update(CoolProp.PT_INPUTS, p, T)
            assert _relerr(tracer.keyed_output(CoolProp.iHmass), flash.hmass()) < 1e-8

    def test_the_slope_is_what_decides_the_branch_there(self, Common):
        """Without the derivative probe the same points cannot be placed

        This is what stops the previous test passing for the wrong reason: at
        a temperature where the bracket is not degenerate, the value order
        alone would do the job and the slope would be dead weight.
        """
        T = 305.0
        state = AbstractState("HEOS", "R513A.mix")
        sat = AbstractState("HEOS", "R513A.mix")
        sat.update(CoolProp.QT_INPUTS, 0.0, T)
        psat = sat.p()
        tracer = Common.IsoLineTracer(state, CoolProp.iP, CoolProp.iT, CoolProp.iT)
        tracer._branch_slope = lambda *args: None
        for p in [0.5 * psat, 0.9 * psat, 1.1 * psat, 2.0 * psat]:
            with pytest.raises(ValueError):
                tracer.update(p, T)

    def test_wrong_branch_is_handed_back(self, Common):
        """A root inside the dome is rejected rather than reported"""
        state = AbstractState("HEOS", "Water")
        tracer = Common.IsoLineTracer(state, CoolProp.iP, CoolProp.iT, CoolProp.iT)
        sat = AbstractState("HEOS", "Water")
        sat.update(CoolProp.QT_INPUTS, 0.0, 400.0)
        rho_liq = sat.rhomolar()
        sat.update(CoolProp.QT_INPUTS, 1.0, 400.0)
        rho_vap = sat.rhomolar()
        bracket = {'liq': (sat.p(), 400.0, rho_liq), 'vap': (sat.p(), 400.0, rho_vap)}
        mid = 0.5 * (rho_liq + rho_vap)
        with pytest.raises(ValueError):
            tracer._check_branch(bracket, tracer.LIQUID, 400.0, mid)
        with pytest.raises(ValueError):
            tracer._check_branch(bracket, tracer.VAPOUR, 400.0, mid)
        # the saturation densities themselves are on their own branch
        tracer._check_branch(bracket, tracer.LIQUID, 400.0, rho_liq)
        tracer._check_branch(bracket, tracer.VAPOUR, 400.0, rho_vap)

    def test_an_undecidable_unbracketed_root_is_handed_back(self, Common):
        """A dome check that cannot answer is not the same as "no dome"

        For an isoline indexed by temperature, the check runs the very call at
        the very temperature that building the bracket already failed on, so
        it can add nothing.  Reading that as "beyond the dome" left a band on
        R504 -- above where its saturation solver converges, below its real
        critical point at 335.5 K, which its phase envelope stops 37 K short
        of -- where roots inside the dome were accepted unexamined.
        """
        tracer = Common.IsoLineTracer(AbstractState("HEOS", "R504.mix"),
                                      CoolProp.iP, CoolProp.iT, CoolProp.iT)
        sat = AbstractState("HEOS", "R504.mix")
        undecidable = []
        for T in (325.0, 330.0, 334.0):
            try:
                sat.update(CoolProp.QT_INPUTS, 0.0, T)
            except Exception:
                undecidable.append(T)
        assert undecidable, "this fluid is meant to defeat the saturation solver here"
        assert tracer._dome_T_limit > max(undecidable), \
            "the band has to be below the temperature where a dome can still exist"
        for T in undecidable:
            with pytest.raises(ValueError):
                tracer._check_outside_the_dome(T, 5000.0)
        # and well above the critical point there really is no dome, cheaply
        tracer._check_outside_the_dome(tracer._dome_T_limit * 1.5, 5000.0)

    def test_an_unbracketed_root_is_still_checked(self, Common):
        """Without a bracket the dome is checked at the root's own temperature

        The reported extent of the phase envelope is not enough on its own:
        for R504 the trace stops about 20 K short of the real cricondentherm,
        and points in between used to be accepted unexamined -- inside the
        dome, as it turns out.
        """
        tracer = Common.IsoLineTracer(AbstractState("HEOS", "R504.mix"),
                                      CoolProp.iP, CoolProp.iSmass, CoolProp.iP)
        sat = AbstractState("HEOS", "R504.mix")
        T = 314.54                       # beyond the envelope, inside the dome
        sat.update(CoolProp.QT_INPUTS, 0.0, T); rho_liq = sat.rhomolar()
        sat.update(CoolProp.QT_INPUTS, 1.0, T); rho_vap = sat.rhomolar()
        assert rho_liq > rho_vap
        with pytest.raises(ValueError):
            tracer._check_branch(None, None, T, 0.5 * (rho_liq + rho_vap))
        # and a root that really is outside passes
        tracer._check_branch(None, None, T, rho_liq * 1.05)
        tracer._check_branch(None, None, T, rho_vap * 0.95)


class TestIsoLineIntegration:
    def test_tracing_is_on_by_default_and_switching_it_off_reaches_the_flash(self, Common):
        """The flag has to change what calc_range does, not just its own value"""
        import numpy as np
        line = Common.IsoLine(CoolProp.iT, CoolProp.iHmass, CoolProp.iP, value=300.0,
                              state=AbstractState("HEOS", "Water"))
        assert line.tracing is True
        line.tracing = False
        assert line.tracing is False

        calls = []
        original = Common.IsoLineTracer.update

        def counting(self, v1, v2):
            """Record that the tracer was reached, then behave normally"""
            calls.append(1)
            return original(self, v1, v2)

        Common.IsoLineTracer.update = counting
        try:
            args = (Common, 'HEOS::Water', 'PH', CoolProp.iSmass, 7000.0, 20)
            off_x, off_y = _isoline(*args, tracing=False)
            assert calls == [], "tracing=False still went through the tracer"
            on_x, on_y = _isoline(*args, tracing=True)
            assert calls, "tracing=True never reached the tracer"
        finally:
            Common.IsoLineTracer.update = original
        both = np.isfinite(off_x) & np.isfinite(off_y)
        assert both.any()
        assert np.allclose(on_x[both], off_x[both], rtol=1e-7)
        assert np.allclose(on_y[both], off_y[both], rtol=1e-7)


class TestUnbracketableInputs:
    """An isoline whose inputs are neither T nor p must not be traced

    A phase-imposed Newton iteration has nothing to stop it converging on the
    metastable extension of the EOS inside the two-phase dome, and without T
    or p among the inputs there is no saturation bracket to rule that out.
    Constant-enthalpy lines on a T-s diagram are the case that bites: they are
    ordinary plot furniture and they cross the dome.
    """

    def test_supports_requires_temperature_or_pressure(self, Common):
        """Only pairs that can be bracketed are accepted"""
        assert Common.IsoLineTracer.supports(CoolProp.iP, CoolProp.iSmass)
        assert Common.IsoLineTracer.supports(CoolProp.iHmass, CoolProp.iT)
        assert not Common.IsoLineTracer.supports(CoolProp.iHmass, CoolProp.iSmass)
        assert not Common.IsoLineTracer.supports(CoolProp.iDmass, CoolProp.iSmass)

    def test_constructor_refuses_an_unbracketable_pair(self, Common):
        """supports() is enforced, not merely advisory"""
        with pytest.raises(ValueError):
            Common.IsoLineTracer(AbstractState("HEOS", "Water"),
                                 CoolProp.iHmass, CoolProp.iSmass, CoolProp.iHmass)

    def test_constant_enthalpy_on_a_ts_diagram_matches_the_flash(self, Common):
        """Every point of this line is flashed, so it must equal the flash"""
        import numpy as np
        value = 1.86222e6
        state = AbstractState("HEOS", "Water")
        line = Common.IsoLine(CoolProp.iHmass, CoolProp.iSmass, CoolProp.iT,
                              value=value, state=state)
        _, _, _, pair = line.get_update_pair()
        assert pair == CoolProp.HmassSmass_INPUTS
        svals = np.linspace(4000.0, 7000.0, 40)
        line.calc_range(svals, None)
        flash = AbstractState("HEOS", "Water")
        checked = 0
        for smass, T in zip(svals, line.y):
            if not np.isfinite(T):
                continue
            flash.update(CoolProp.HmassSmass_INPUTS, value, smass)
            assert _relerr(T, flash.T()) < 1e-6, (
                "traced {0} K where the flash gives {1} K at Q={2}".format(T, flash.T(), flash.Q()))
            checked += 1
        assert checked > 20, "the line has to actually cross the region of interest"


class TestNoMetastableRoots:
    """No point the tracer calls single phase may sit inside the dome

    This is the check that a reconstruction from (T, rho) cannot make: the
    metastable root re-evaluates to itself perfectly.  Only the saturation
    densities at the point's own temperature can tell the two apart.
    """

    @pytest.mark.parametrize("fluid,graph,iso_type,value", [
        ("HEOS::Water", 'PH', CoolProp.iSmass, 6000.0),
        ("HEOS::Water", 'TS', CoolProp.iP, 5.0e5),
        ("HEOS::R513A.mix", 'PH', CoolProp.iSmass, 1750.0),
        ("HEOS::R513A.mix", 'PH', CoolProp.iT, 300.0),
    ])
    def test_traced_points_are_outside_the_dome(self, Common, fluid, graph, iso_type, value):
        """Every single-phase root sits outside the saturation densities at its own T"""
        from CoolProp.Plots import PropertyPlot
        import numpy as np
        backend, name = fluid.split("::")
        plot = PropertyPlot(fluid, graph, unit_system='SI', tp_limits='ACHP')
        limits = plot._get_axis_limits()
        xvals = plot.generate_ranges(plot._x_index, limits[0], limits[1], 60)
        yvals = plot.generate_ranges(plot._y_index, limits[2], limits[3], 60)
        line = Common.IsoLine(iso_type, plot._x_index, plot._y_index,
                              value=value, state=plot._state)
        ipos, xpos, ypos, _ = line.get_update_pair()
        keys = [v for (_, v) in sorted(zip([ipos, xpos, ypos],
                                           [iso_type, plot._x_index, plot._y_index]))]
        vals = [v for (_, v) in sorted(zip([ipos, xpos, ypos],
                                           [np.array(value), xvals, yvals]))]
        if vals[0].size < vals[1].size:
            vals[0] = np.resize(vals[0], vals[1].shape)
        elif vals[1].size < vals[0].size:
            vals[1] = np.resize(vals[1], vals[0].shape)

        tracer = Common.IsoLineTracer(plot._state, keys[0], keys[1], iso_type)
        sat = AbstractState(backend, name)
        checked = 0
        for v0, v1 in zip(vals[0], vals[1]):
            try:
                tracer.update(v0, v1)
            except Exception:
                continue
            state = tracer.state
            if 0.0 <= state.Q() <= 1.0:
                continue                       # solved as two phase, which is fine
            T, rho = state.T(), state.rhomass()
            try:
                sat.update(CoolProp.QT_INPUTS, 0.0, T); rho_liq = sat.rhomass()
                sat.update(CoolProp.QT_INPUTS, 1.0, T); rho_vap = sat.rhomass()
            except Exception:
                continue                       # no dome at this T to be inside of
            checked += 1
            assert not rho_vap * (1.0 + 1e-6) < rho < rho_liq * (1.0 - 1e-6), (
                "traced a single-phase root at rho={0} inside the dome "
                "({1} .. {2}) at T={3}".format(rho, rho_vap, rho_liq, T))
        assert checked >= 5, "the sweep never reached a temperature that has a dome"


class TestBracketIsNotSilentlySkipped:
    """A bracket that could not be built must never read as "no bracket needed"

    ``_get_bracket`` returns None for a point beyond the dome, where an
    unchecked single-phase root is safe, and raises everywhere else.  Two ways
    that distinction has been lost, both of which left whole isolines traced
    with the branch check disabled:

    * recording the saturation value before the build, so that a build which
      raised left the *previous* outcome (None) cached against the new value.
      An isoline in T or p has one saturation value for its whole length, so a
      single failure disabled the check for every remaining point;
    * caching the phase envelope between tracers.  Building it is what makes
      the mixture saturation flashes on that state converge, so the first
      isoline of a fluid bracketed and the rest did not.
    """

    def test_a_failed_build_keeps_failing_instead_of_going_unchecked(self, Common):
        """A cached bracket failure must keep raising, not read as no-dome"""
        tracer = Common.IsoLineTracer(AbstractState("HEOS", "R513A.mix"),
                                      CoolProp.iP, CoolProp.iSmass, CoolProp.iP)
        p, smass = 1.0e6, 1750.0
        tracer.update(p, smass)                 # works, and caches a real bracket
        assert tracer._bracket is not None

        def always_fails(sat_value, target_index):
            """Stand in for a saturation solver that will not converge"""
            raise ValueError("saturation solver gave up")

        tracer._build_bracket = always_fails
        tracer._sat_value = None                # force one rebuild, which fails
        for _ in range(3):                      # same saturation value each time
            with pytest.raises(ValueError):
                tracer.update(p, smass)
            assert tracer._bracket is None
            assert tracer._bracket_failure is not None

    @pytest.mark.parametrize("attempt", [1, 2, 3])
    def test_every_tracer_for_a_fluid_brackets_the_same(self, Common, attempt):
        """The Nth tracer in a process must behave like the first

        3.4 MPa is below R513A's cricondenbar but high enough that PQ_INPUTS
        needs the phase envelope to converge, which is exactly where a cached
        envelope used to leave later isolines unbracketed.
        """
        from CoolProp.Plots import PropertyPlot
        import numpy as np
        plot = PropertyPlot('HEOS::R513A.mix', 'TS', unit_system='SI', tp_limits='ACHP')
        limits = plot._get_axis_limits()
        svals = plot.generate_ranges(plot._x_index, limits[0], limits[1], 40)
        pressure = 3.4e6
        flash = AbstractState("HEOS", "R513A.mix")
        for _ in range(attempt):                # burn N-1 tracers first
            tracer = Common.IsoLineTracer(plot._state, CoolProp.iP, CoolProp.iSmass, CoolProp.iP)
        unchecked = compared = 0
        for smass in svals:
            try:
                tracer.update(pressure, smass)
            except Exception:
                continue
            if tracer._bracket is None:
                unchecked += 1
                continue
            try:
                flash.update(CoolProp.PSmass_INPUTS, pressure, smass)
            except Exception:
                continue           # the flash cannot do this point; that is the premise
            compared += 1
            assert _relerr(tracer.state.T(), flash.T()) < 1e-6
        assert unchecked == 0, "{0} points traced without a saturation bracket".format(unchecked)
        assert compared > 20


class TestFluidModelIsCarriedOver:
    def test_custom_binary_interaction_parameters_reach_the_tracer(self, Common):
        """The tracer clones the state; a modified kij must come with it

        Otherwise one isoline mixes two fluid models -- traced points from the
        defaults, flashed points from the caller's parameters.
        """
        state = AbstractState("HEOS", "R1234yf&R134a")
        state.set_mole_fractions([0.5, 0.5])
        state.set_binary_interaction_double(0, 1, "betaT", 1.05)
        default = AbstractState("HEOS", "R1234yf&R134a")
        default.set_mole_fractions([0.5, 0.5])

        p, T = 1.0e6, 320.0
        state.update(CoolProp.PT_INPUTS, p, T)
        default.update(CoolProp.PT_INPUTS, p, T)
        assert _relerr(state.hmass(), default.hmass()) > 1e-9, \
            "the modified parameter has to actually change the answer"

        tracer = Common.IsoLineTracer(state, CoolProp.iP, CoolProp.iT, CoolProp.iT)
        tracer.update(p, T)
        assert _relerr(tracer.keyed_output(CoolProp.iHmass), state.hmass()) < 1e-9


class TestGapFilling:
    def test_mixture_isentrope_has_far_fewer_gaps(self, Common):
        """Flashing an R513A isentrope cold leaves holes in it; tracing mostly does not

        Not zero: a handful of near-critical points cannot be bracketed and
        cannot be flashed either, and those are handed on as non-finite rather
        than guessed at.  The point is the order of magnitude.
        """
        import numpy as np
        args = (Common, 'HEOS::R513A.mix', 'PH', CoolProp.iSmass, 1950.0, 40)
        x, y = _isoline(*args, tracing=True)
        xf, yf = _isoline(*args, tracing=False)
        traced_holes = int(np.sum(~np.isfinite(x) | ~np.isfinite(y)))
        flashed_holes = int(np.sum(~np.isfinite(xf) | ~np.isfinite(yf)))
        assert flashed_holes > 0, "this line is meant to be one the cold flash cannot do"
        assert traced_holes * 3 < flashed_holes, (
            "tracing left {0} holes against the flash's {1}".format(traced_holes, flashed_holes))
        # and where both landed, they agree
        both = (np.isfinite(xf) & np.isfinite(yf) & np.isfinite(x) & np.isfinite(y))
        assert both.sum() > 20
        assert np.allclose(x[both], xf[both], rtol=1e-7)
        assert np.allclose(y[both], yf[both], rtol=1e-7)
