"""State-shim parity regression suite (bd CoolProp-r9sq.26).

Locks the legacy ``CoolProp.State.State`` surface that PDSim and other downstream
code depend on but the v8 capsule shim had dropped: the widened constructor
(``params=None`` / ``phase=`` / ``backend=`` / ``BACKEND::Fluid``), ``set_Fluid``,
``Phase``, and the refrigeration quantities ``Tsat`` / ``subcooling`` /
``superheat``.  Units follow the legacy convention (pressures in kPa).
"""
import pytest

from CoolProp.State import State


class TestConstructorWidening:
    def test_none_params_no_update(self):
        S = State("Water", None)  # bare construction (no state update) -- used by get_Tsat/copy
        assert S.Fluid == b"Water"

    def test_phase_kwarg_accepted(self):
        # Legacy accepted a (deprecated) phase= flag; must not TypeError.
        S = State("Water", {"T": 300.0, "P": 101.325}, phase="liquid")
        assert S.T == pytest.approx(300.0)

    def test_backend_colon_syntax(self):
        S = State("HEOS::Water", {"T": 300.0, "D": 1000.0})
        assert S.T == pytest.approx(300.0)

    def test_backend_kwarg(self):
        S = State("Water", {"T": 300.0, "D": 1000.0}, backend="HEOS")
        assert S.rho == pytest.approx(1000.0)


class TestSetFluid:
    def test_set_fluid_recreates(self):
        S = State("Water", {"T": 300.0, "D": 1000.0})
        S.set_Fluid("R134a", "HEOS")
        S.update({"T": 300.0, "Q": 0.0})
        assert S.Fluid == b"R134a"

    def test_bracket_mixture_construction(self):
        # Bracketed mole-fraction mixtures are fully supported: the fractions are
        # parsed out and set on the handle (like the legacy State.set_Fluid).
        import math

        S = State("R32[0.5]&R125[0.5]", {"T": 300.0, "P": 2000.0})  # P in kPa
        assert b"R32" in S.Fluid and b"R125" in S.Fluid
        assert S.T == pytest.approx(300.0)
        assert math.isfinite(S.rho) and S.rho > 0

    def test_set_fluid_to_mixture(self):
        import math

        S = State("Water", {"T": 300.0, "D": 1000.0})
        S.set_Fluid("R32[0.7]&R125[0.3]", "HEOS")
        S.update({"T": 300.0, "P": 2000.0})
        assert math.isfinite(S.rho) and S.rho > 0


class TestPhaseAndSaturation:
    def test_phase_returns_int(self):
        S = State("Water", {"T": 300.0, "D": 1000.0})
        assert isinstance(S.Phase(), int)

    def test_Tsat_two_phase(self):
        S = State("Water", {"P": 101.325, "Q": 0.0})
        assert S.get_Tsat(1.0) == pytest.approx(373.12, abs=1.0)
        assert S.Tsat == pytest.approx(373.12, abs=1.0)

    def test_Tsat_above_critical_is_none(self):
        S = State("Water", {"T": 800.0, "P": 50000.0})  # 50 MPa > 22 MPa critical
        assert S.get_Tsat(1.0) is None

    def test_superheat_positive(self):
        S = State("Water", {"T": 450.0, "P": 101.325})  # superheated steam at 1 atm
        sh = S.get_superheat()
        assert sh is not None and sh > 0
        assert S.superheat == pytest.approx(sh)

    def test_subcooling_positive(self):
        S = State("Water", {"T": 300.0, "P": 101.325})  # subcooled liquid at 1 atm
        sc = S.get_subcooling()
        assert sc is not None and sc > 0
        assert S.subcooling == pytest.approx(sc)

    def test_superheat_none_outside_two_phase(self):
        S = State("Water", {"T": 800.0, "P": 50000.0})
        assert S.get_superheat() is None


class TestMiscParity:
    def test_prandtl_property(self):
        S = State("Water", {"T": 300.0, "P": 101.325})
        assert S.Prandtl == pytest.approx(S.cp * S.visc / S.k)

    def test_props_negative_key_raises(self):
        S = State("Water", {"T": 300.0, "D": 1000.0})
        with pytest.raises(ValueError):
            S.Props(-1)
