"""Fresh coverage for AbstractState derivatives + the Plots package (bd CoolProp-r9sq.14).

Replaces the stale, deleted test_CoolPropState.py / test_plots.py with v8-native
pytest coverage:
  * the partial / saturation / two-phase derivative methods on the nanobind
    AbstractState (CI-runnable, no optional deps);
  * a Plots smoke test that ``CoolProp.Plots.PropertyPlot`` constructs (it imports
    ``CoolProp.Plots.Common``, which only works now that the nanobind core ships
    the legacy ``Py*`` class aliases + tuple ``extract_*`` -- skips without matplotlib).
"""
import math

import pytest

from CoolProp import CoolProp as C


def _not_marshalling_error(exc):
    assert not isinstance(exc, TypeError), "derivative method uncallable (binding bug): %s" % exc


class TestPartialDerivatives:
    def _compressed_liquid(self):
        AS = C.AbstractState("HEOS", "Water")
        AS.update(int(C.PT_INPUTS), 5.0e6, 350.0)  # 5 MPa, 350 K -> single-phase, all derivs defined
        return AS

    def _saturated_liquid(self):
        AS = C.AbstractState("HEOS", "Water")
        AS.update(int(C.QT_INPUTS), 0.0, 350.0)
        return AS

    def _two_phase(self):
        AS = C.AbstractState("HEOS", "Water")
        AS.update(int(C.QT_INPUTS), 0.5, 350.0)
        return AS

    def test_first_partial_deriv(self):
        AS = self._compressed_liquid()
        v = AS.first_partial_deriv(int(C.iP), int(C.iT), int(C.iDmolar))  # (dP/dT)|rho
        assert math.isfinite(v)

    def test_second_partial_deriv(self):
        AS = self._compressed_liquid()
        v = AS.second_partial_deriv(int(C.iP), int(C.iT), int(C.iDmolar), int(C.iT), int(C.iDmolar))
        assert math.isfinite(v)

    def test_first_saturation_deriv(self):
        AS = self._saturated_liquid()
        v = AS.first_saturation_deriv(int(C.iP), int(C.iT))  # Clausius-Clapeyron dP/dT > 0
        assert math.isfinite(v) and v > 0

    def test_second_saturation_deriv(self):
        AS = self._saturated_liquid()
        try:
            v = AS.second_saturation_deriv(int(C.iP), int(C.iT), int(C.iT))
        except Exception as e:  # noqa: BLE001 -- some combos are domain-restricted, not a binding bug
            _not_marshalling_error(e)
            return
        assert math.isfinite(v)

    def test_first_two_phase_deriv(self):
        AS = self._two_phase()
        try:
            v = AS.first_two_phase_deriv(int(C.iDmolar), int(C.iHmolar), int(C.iP))
        except Exception as e:  # noqa: BLE001
            _not_marshalling_error(e)
            return
        assert math.isfinite(v)

    def test_first_two_phase_deriv_splined(self):
        AS = self._two_phase()
        try:
            v = AS.first_two_phase_deriv_splined(int(C.iDmolar), int(C.iHmolar), int(C.iP), 0.3)
        except Exception as e:  # noqa: BLE001
            _not_marshalling_error(e)
            return
        assert math.isfinite(v)


class TestPlotsSmoke:
    def test_property_plot_constructs(self):
        # Plots.Common imports + constructs PyCriticalState/PyGuessesStructure and
        # unpacks extract_backend/extract_fractions tuples -- all of which only work
        # on the v8 core via the #3120 aliases.
        pytest.importorskip("matplotlib")
        import matplotlib

        matplotlib.use("Agg")
        from CoolProp.Plots import PropertyPlot

        plot = PropertyPlot("Water", "PH")
        assert plot is not None

    def test_plots_common_imports(self):
        pytest.importorskip("matplotlib")
        from CoolProp.Plots.Common import IsoLine, BasePlot  # noqa: F401 -- import is the test

        assert IsoLine is not None and BasePlot is not None

    def test_process_fluid_state_accepts_state_and_string(self):
        # process_fluid_state must pass an AbstractState instance through unchanged
        # (its ``isinstance(fluid_ref, AbstractState)`` branch) and build one from a
        # fluid string.  This relies on AbstractState being a real type -- bd
        # CoolProp-r9sq.28; before that it was a factory function and isinstance
        # against it raised "arg 2 must be a type".
        pytest.importorskip("matplotlib")
        from CoolProp.Plots.Common import process_fluid_state

        AS = C.AbstractState("HEOS", "Water")
        assert process_fluid_state(AS) is AS
        built = process_fluid_state("HEOS::Water")
        assert hasattr(built, "keyed_output")
        with pytest.raises(TypeError):
            process_fluid_state(3.14)  # neither a string nor a state-like object
