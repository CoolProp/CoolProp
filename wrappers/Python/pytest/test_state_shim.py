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
        S = State("Water", {"T": 300.0, "P": 101.325}, phase="liquid")
        assert S.T == pytest.approx(300.0)

    def test_phase_is_imposed_not_swallowed(self):
        # The phase= flag must actually be imposed (like legacy), not silently
        # ignored: at a near-saturation P,T it forces the liquid vs (metastable)
        # gas root, giving very different densities.
        liq = State("Water", {"T": 372.0, "P": 101.325}, phase="liquid")
        gas = State("Water", {"T": 372.0, "P": 101.325}, phase="gas")
        assert liq.rho > 100.0 * gas.rho  # ~960 kg/m^3 vs ~0.6 kg/m^3

    def test_view_specify_phase(self):
        # The pAS view exposes specify_phase (legacy State.pAS.specify_phase).
        from CoolProp import CoolProp as C

        S = State("Water", {"T": 300.0, "D": 1000.0})
        S.pAS.specify_phase(int(C.iphase_liquid))  # must not raise

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

    def test_copy_pure_fluid_independent(self):
        # copy() rebuilds from the full spec and is an independent object.
        S = State("Water", {"T": 350.0, "P": 2000.0})
        c = S.copy()
        assert c.Fluid == b"Water"
        assert c.rho == pytest.approx(S.rho, rel=1e-9)
        c.update({"T": 400.0, "P": 2000.0})
        assert S.T == pytest.approx(350.0)  # the original is untouched

    def test_get_Tsat_mixture_keeps_composition(self):
        # get_Tsat rebuilds a temp State from the full spec, so a mixture keeps its
        # fractions -- without that the temp state would have no composition and
        # fail. (copy() of a mixture is separately limited by CoolProp's DmassT
        # flash not supporting mixtures, exactly as the legacy State was.)
        import math

        S = State("R32[0.5]&R125[0.5]", {"T": 280.0, "P": 800.0})
        Tsat = S.get_Tsat(1.0)
        assert Tsat is None or (math.isfinite(Tsat) and Tsat > 0)

    def test_failed_set_fluid_keeps_state_usable(self):
        # A failed re-fluid (bad name) must leave the prior state + .pAS valid --
        # the new handle is made before the old one is freed (no use-after-free).
        S = State("Water", {"T": 300.0, "D": 1000.0})
        rho_before = S.rho
        with pytest.raises(Exception):  # noqa: PT011 -- any failure, just must not corrupt state
            S.set_Fluid("NotARealFluid", "HEOS")
        assert S.rho == pytest.approx(rho_before, rel=1e-9)
        assert S.pAS.rhomass() == pytest.approx(rho_before, rel=1e-9)

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


class TestGenericInputPairs:
    """update() accepts any pair generate_update_pair resolves, not just the 5
    common ones -- matching the legacy State.update (bd CoolProp-r9sq.26)."""

    def test_PS_pair(self):
        # {P, S}: round-trip a superheated-vapor state through pressure + entropy.
        ref = State("Water", {"T": 400.0, "P": 200.0})  # 200 kPa, vapor
        S = State("Water", {"P": 200.0, "S": ref.s})    # S in legacy kJ/kg/K
        assert S.T == pytest.approx(400.0, rel=1e-3)

    def test_DP_pair(self):
        # {D, P}: density + pressure (note the non-default key order, too).
        ref = State("Water", {"T": 300.0, "D": 1000.0})
        S = State("Water", {"D": 1000.0, "P": ref.p})   # P in kPa
        assert S.T == pytest.approx(300.0, rel=1e-3)

    def test_PH_still_works(self):
        # The previously-special-cased pair must still behave identically.
        ref = State("Water", {"T": 350.0, "P": 101.325})
        S = State("Water", {"P": 101.325, "H": ref.h})  # H in kJ/kg
        assert S.rho == pytest.approx(ref.rho, rel=1e-3)

    def test_unsupported_key_raises(self):
        with pytest.raises(ValueError):
            State("Water", {"T": 300.0, "X": 5.0})

    def test_one_input_raises(self):
        with pytest.raises(ValueError):
            State("Water", {"T": 300.0})


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
