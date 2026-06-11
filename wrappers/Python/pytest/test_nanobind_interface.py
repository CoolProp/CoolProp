"""
pytest suite for the v8 (nanobind) CoolProp Python interface.

Runs against an *installed* CoolProp built with ``-DCOOLPROP_NANOBIND=ON`` (see
README.md in this directory).  Most assertions hold for the legacy build too;
the v8-specific cases (the ``Props``/``HAProps`` removals) are skipped when run
against a legacy wheel so the file is safe to point at either.

Covers: the high-level API, the AbstractState low-level surface including the
q2sh parity additions (idealgas/residual decompositions, fluid-constant
accessors, ``neff``, the ``set_reference_state`` dispatcher), the module-level
convenience wrappers, the capsule ``State`` shim unit conventions (the kPa/kJ /
g/mol / kW/m/K split that the original PDSim regression hinged on), and the
approved v8 removals.
"""
import math

import pytest

import CoolProp
from CoolProp import CoolProp as CP

# v8 (nanobind) build is identified by the free ``Props`` having been removed.
is_legacy = hasattr(CP, "Props")
v8_only = pytest.mark.skipif(is_legacy, reason="v8 (nanobind) build only")


# --------------------------------------------------------------------------- #
# High-level API
# --------------------------------------------------------------------------- #
def test_propssi_liquid_water():
    # liquid water density at 300 K, 1 atm
    assert CP.PropsSI("D", "T", 300.0, "P", 101325.0, "Water") == pytest.approx(996.5, rel=1e-3)


def test_propssi_co2():
    assert CP.PropsSI("D", "T", 320.0, "P", 6e6, "CO2") == pytest.approx(139.1, rel=1e-2)


def test_hapropssi_wetbulb():
    assert CP.HAPropsSI("B", "T", 300.0, "P", 101325.0, "W", 0.01) == pytest.approx(291.706, rel=1e-4)


# --------------------------------------------------------------------------- #
# AbstractState low-level surface
# --------------------------------------------------------------------------- #
@pytest.fixture
def water():
    return CP.AbstractState("HEOS", "Water")


def test_abstractstate_update_pt(water):
    water.update(CP.PT_INPUTS, 101325.0, 300.0)
    assert water.rhomass() == pytest.approx(996.5, rel=1e-3)
    assert water.T() == pytest.approx(300.0)


def test_abstractstate_update_failure_raises_valueerror(water):
    # A failing AbstractState.update() must raise ValueError, matching the
    # legacy Cython wrapper.  Without the nanobind exception translator the
    # underlying CoolPropBaseError surfaces as RuntimeError, silently breaking
    # the canonical `try: ... except ValueError` pattern (incl. CoolProp's own
    # shipped Plots code).  CoolProp-1tbe.2.
    with pytest.raises(ValueError):
        water.update(CP.PT_INPUTS, -1.0, -1.0)


def test_get_fluid_constant(water):
    assert water.get_fluid_constant(0, CP.iT_critical) == pytest.approx(647.096, rel=1e-4)


def test_neff_finite(water):
    water.update(CP.DmassT_INPUTS, 2.0, 600.0)
    assert math.isfinite(water.neff())


# --- q2sh parity additions: the idealgas + residual decomposition ---
@pytest.mark.parametrize(
    "D, T, label",
    [(0.5, 400.0, "superheated vapor"), (950.0, 400.0, "compressed liquid"), (2.0, 600.0, "vapor")],
)
def test_h_eq_ideal_plus_residual_singlephase(D, T, label):
    AS = CP.AbstractState("HEOS", "Water")
    AS.update(CP.DmassT_INPUTS, D, T)
    assert AS.hmolar() == pytest.approx(AS.hmolar_idealgas() + AS.hmolar_residual(), rel=1e-9)
    assert AS.smolar() == pytest.approx(AS.smolar_idealgas() + AS.smolar_residual(), rel=1e-9)


def test_h_neq_ideal_plus_residual_twophase():
    # At (D=2, T=400) water is two-phase (p == Psat); hmolar() is the
    # quality-weighted mixture, so the EOS-at-(tau,delta) parts do NOT sum to it.
    AS = CP.AbstractState("HEOS", "Water")
    AS.update(CP.DmassT_INPUTS, 2.0, 400.0)
    assert 0.0 < AS.Q() < 1.0
    assert AS.hmolar() != pytest.approx(AS.hmolar_idealgas() + AS.hmolar_residual(), rel=1e-3)


# --- set_reference_state unified (FluidName, *args) dispatcher ---
def test_set_reference_state_string_dispatch():
    try:
        CP.set_reference_state("Propane", "IIR")  # arity 1 -> set_reference_stateS
    finally:
        CP.set_reference_state("Propane", "DEF")  # reset global state


def test_set_reference_state_bad_arity_raises():
    with pytest.raises(ValueError):
        CP.set_reference_state("Propane", 1.0, 2.0, 3.0)  # arity 3 -> invalid


# --------------------------------------------------------------------------- #
# Module-level convenience wrappers
# --------------------------------------------------------------------------- #
def test_fluidslist():
    fl = CP.FluidsList()
    assert "Water" in fl
    assert "" not in fl  # the empty-input split fix: no stray '' entry


def test_get_aliases():
    assert "R744" in CP.get_aliases("CO2")


def test_get_refpropname():
    assert CP.get_REFPROPname("R134a") == "R134A"


def test_get_bibtexkey():
    assert CP.get_BibTeXKey("Water", "EOS")  # non-empty string


def test_config_int_bound():
    # smoke: the int-config accessors are bound (a value round-trip needs a known
    # integer config key, which the public enum does not expose stably)
    assert callable(CP.get_config_int) and callable(CP.set_config_int)


# --------------------------------------------------------------------------- #
# Capsule State shim -- the legacy unit conventions PDSim depends on
# --------------------------------------------------------------------------- #
def test_state_unit_conventions():
    from CoolProp.State import State

    S = State("Water", {"T": 320.0, "D": 1000.0})
    # kPa / kJ getters vs the SI Props()
    assert S.get_p() == pytest.approx(S.Props(CP.iP) / 1000.0, rel=1e-9)
    assert S.get_h() == pytest.approx(S.Props(CP.iHmass) / 1000.0, rel=1e-9)
    # get_MM is g/mol == 1000x the SI molar mass from the same handle (~18.015)
    assert S.get_MM() == pytest.approx(S.pAS.keyed_output(CP.imolar_mass) * 1000.0, rel=1e-9)
    assert S.get_MM() == pytest.approx(18.015, rel=1e-3)
    # get_rho and .pAS are SI
    assert S.get_rho() == pytest.approx(1000.0, rel=1e-9)
    assert S.pAS.rhomass() == pytest.approx(1000.0, rel=1e-9)


# --------------------------------------------------------------------------- #
# Approved v8 removals
# --------------------------------------------------------------------------- #
@v8_only
def test_free_props_removed():
    assert not hasattr(CP, "Props")  # non-SI Props (deprecated for years) gone


@v8_only
def test_core_haprops_removed():
    assert not hasattr(CP, "HAProps")  # non-SI humid air gone from the core


@v8_only
def test_haprops_submodule_raises():
    import CoolProp.HumidAirProp as HA

    with pytest.raises(NotImplementedError):
        HA.HAProps("B", "T", 300.0, "P", 101325.0, "W", 0.01)


def test_hapropssi_kept_on_submodule():
    import CoolProp.HumidAirProp as HA

    assert HA.HAPropsSI("B", "T", 300.0, "P", 101325.0, "W", 0.01) == pytest.approx(291.706, rel=1e-4)


# --------------------------------------------------------------------------- #
# Import-tree smoke (the r9sq.3 parity, exercised at runtime)
# --------------------------------------------------------------------------- #
def test_toplevel_surface_present():
    for name in ("AbstractState", "CoolProp", "HumidAirProp", "State", "get", "copy_BibTeX_library",
                 "__fluids__", "__version__"):
        assert hasattr(CoolProp, name), name
