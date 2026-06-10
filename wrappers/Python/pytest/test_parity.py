"""Legacy-vs-nanobind parity regression suite (bd CoolProp-r9sq.16/.21/.22/.23/.24).

These lock the behavioral parity gaps a systematic legacy-vs-nanobind audit
surfaced, so they cannot silently regress:

  r9sq.21  scalar PropsSI: int/numpy-scalar args must return a float (not a
           1-element ndarray); a failed evaluation must raise ValueError (not
           return inf, not RuntimeError).
  r9sq.22  HAPropsSI: int/numpy-scalar args -> float; ndarray in -> ndarray out;
           list in -> list out; failed eval -> ValueError.
  r9sq.23  enum completeness: the mass-basis parameters/input_pairs, fluid_types
           and fast_evaluate_status that the hand-written enums had dropped.
  r9sq.24  out-reference / unregistered-return methods are callable and return
           the legacy tuple/dict/object shapes (not TypeError / cast errors).
  r9sq.16  set_reference_state D-form accepts numpy/array-typed numeric args.

The bar is the legacy *contract* (type/shape/exception class), not exact values
-- the numbers are the C++ core's concern and are covered elsewhere.
"""
import numpy as np
import pytest

from CoolProp import CoolProp as C
from CoolProp.CoolProp import PropsSI, HAPropsSI


# --------------------------------------------------------------------------- #
# r9sq.21 -- scalar PropsSI returns float, raises ValueError on failure
# --------------------------------------------------------------------------- #
class TestPropsSIScalar:
    def test_int_args_return_float(self):
        # The canonical call: int pressure / quality literals. Must be a float,
        # not a 1-element ndarray (the regression that broke ~every example).
        r = PropsSI("T", "P", 101325, "Q", 0, "Water")
        assert isinstance(r, float), f"int-arg PropsSI returned {type(r).__name__}"

    def test_float_args_return_float(self):
        assert isinstance(PropsSI("T", "P", 101325.0, "Q", 0.0, "Water"), float)

    def test_numpy_scalar_args_return_float(self):
        r = PropsSI("T", "P", np.float64(101325), "Q", np.int64(0), "Water")
        assert isinstance(r, float)

    def test_mixed_int_float_return_float(self):
        assert isinstance(PropsSI("D", "T", 300, "P", 101325.0, "Water"), float)

    def test_list_args_return_ndarray(self):
        r = PropsSI("T", "P", [101325.0, 200000.0], "Q", [0.0, 0.0], "Water")
        assert isinstance(r, np.ndarray) and r.shape == (2,)

    def test_ndarray_args_return_ndarray(self):
        P = np.array([1e5, 2e5, 3e5])
        r = PropsSI("T", "P", P, "Q", 0.0, "Water")
        assert isinstance(r, np.ndarray) and r.shape == (3,)

    def test_bad_fluid_raises_valueerror(self):
        with pytest.raises(ValueError):
            PropsSI("T", "P", 101325.0, "Q", 0.0, "Nonexistium")

    def test_bad_output_raises_valueerror(self):
        with pytest.raises(ValueError):
            PropsSI("ZZZ", "P", 101325.0, "Q", 0.0, "Water")

    def test_out_of_range_raises_valueerror(self):
        with pytest.raises(ValueError):
            PropsSI("T", "P", 101325.0, "Q", 5.0, "Water")  # Q out of [0,1]

    def test_bad_input_never_returns_inf(self):
        # Belt-and-suspenders: whatever happens, must not silently return inf.
        try:
            v = PropsSI("T", "P", 101325.0, "Q", 0.0, "Nonexistium")
        except ValueError:
            return
        assert np.isfinite(v), "scalar PropsSI returned a non-finite value instead of raising"

    def test_two_arg_trivial_form(self):
        assert isinstance(PropsSI("Tcrit", "Water"), float)

    def test_two_arg_bad_fluid_raises_valueerror(self):
        with pytest.raises(ValueError):
            PropsSI("Tcrit", "Nonexistium")

    def test_array_bad_eval_raises_valueerror(self):
        with pytest.raises(ValueError):
            PropsSI("ZZZ", "P", [1e5], "Q", [0.0], "Water")


# --------------------------------------------------------------------------- #
# vecPropsSI -- multi-output PropsSI: a sequence first argument (list of output
# names) returns a matrix, np.squeeze'd exactly like the legacy Cython wrapper
# (issue #2736).  Shapes, with n = broadcast input length, m = number of outputs:
#   n>1, m>1 -> (n, m)    n>1, m==1 -> (n,)
#   n==1, m>1 -> (m,)     n==1, m==1 -> (1,)
# --------------------------------------------------------------------------- #
class TestPropsSIMultiOutput:
    def test_issue_2736_2d_shape_and_values(self):
        # The exact reported call: 2 outputs x 3 input points -> (3, 2).
        T = [300.0, 600.0, 900.0]
        r = PropsSI(["Dmass", "viscosity"], "T", T, "P", 101325, "Water")
        assert isinstance(r, np.ndarray) and r.shape == (3, 2)
        # Column j must equal the single-output PropsSI for that output.
        for i, Ti in enumerate(T):
            assert r[i, 0] == pytest.approx(PropsSI("Dmass", "T", Ti, "P", 101325, "Water"))
            assert r[i, 1] == pytest.approx(PropsSI("viscosity", "T", Ti, "P", 101325, "Water"))

    def test_list_outputs_scalar_inputs_1d(self):
        # m>1, n==1 -> squeeze to (m,)
        r = PropsSI(["T", "Dmass"], "P", 101325, "Q", 0, "Water")
        assert isinstance(r, np.ndarray) and r.shape == (2,)
        assert r[0] == pytest.approx(PropsSI("T", "P", 101325, "Q", 0, "Water"))
        assert r[1] == pytest.approx(PropsSI("Dmass", "P", 101325, "Q", 0, "Water"))

    def test_single_element_list_output_vector_inputs_1d(self):
        # m==1, n>1 -> squeeze to (n,)
        r = PropsSI(["T"], "P", [1e5, 2e5], "Q", [0.0, 0.0], "Water")
        assert isinstance(r, np.ndarray) and r.shape == (2,)

    def test_single_element_list_output_scalar_inputs_1d(self):
        # m==1, n==1 -> never collapses below 1-D: (1,)
        r = PropsSI(["T"], "P", 101325, "Q", 0, "Water")
        assert isinstance(r, np.ndarray) and r.shape == (1,)
        assert r[0] == pytest.approx(PropsSI("T", "P", 101325, "Q", 0, "Water"))

    def test_tuple_outputs_like_list(self):
        r = PropsSI(("Dmass", "viscosity"), "T", [300.0, 600.0], "P", 101325, "Water")
        assert isinstance(r, np.ndarray) and r.shape == (2, 2)

    def test_ndarray_str_outputs(self):
        r = PropsSI(np.array(["Dmass", "viscosity"]), "T", [300.0, 600.0], "P", 101325, "Water")
        assert isinstance(r, np.ndarray) and r.shape == (2, 2)

    def test_scalar_broadcast_against_vector(self):
        # One vector input, one scalar input, multiple outputs -> (n, m).
        r = PropsSI(["T", "Dmass"], "P", [1e5, 2e5, 3e5], "Q", 0, "Water")
        assert isinstance(r, np.ndarray) and r.shape == (3, 2)

    def test_empty_input_returns_empty(self):
        r = PropsSI(["T", "Dmass"], "P", [], "Q", [], "Water")
        assert isinstance(r, np.ndarray) and r.shape == (0,)

    def test_failed_eval_raises_valueerror(self):
        # An invalid output name yields an empty matrix -> ValueError (parity with
        # the single-output vectorized path, not RuntimeError).  Use a VALID fluid
        # so the raise is provoked by the bad output, not by backend init.
        with pytest.raises(ValueError):
            PropsSI(["ZZZ"], "P", [1e5], "Q", [0.0], "Water")


# --------------------------------------------------------------------------- #
# r9sq.22 -- HAPropsSI scalar/array typing
# --------------------------------------------------------------------------- #
class TestHAPropsSI:
    def test_int_arg_returns_float(self):
        # Flagship doc example: integer pressure literal -> must be a float so
        # `h/1000` works (the regression returned a 1-element list).
        h = HAPropsSI("H", "T", 298.15, "P", 101325, "R", 0.5)
        assert isinstance(h, float)

    def test_all_float_returns_float(self):
        assert isinstance(HAPropsSI("H", "T", 298.15, "P", 101325.0, "R", 0.5), float)

    def test_numpy_scalar_returns_float(self):
        assert isinstance(HAPropsSI("H", "T", 298.15, "P", np.float64(101325), "R", 0.5), float)

    def test_list_input_returns_list(self):
        out = HAPropsSI("W", "R", 0.5, "P", 101325.0, "T", [290.0, 300.0, 310.0])
        assert isinstance(out, list) and len(out) == 3

    def test_ndarray_input_returns_ndarray(self):
        T = np.linspace(290, 310, 5)
        out = HAPropsSI("W", "R", 0.5, "P", 101325.0, "T", T)
        assert isinstance(out, np.ndarray) and out.shape == (5,)

    def test_bad_input_raises_valueerror(self):
        with pytest.raises(ValueError):
            HAPropsSI("T", "P", 101325.0, "R", 2.0, "B", 290.0)  # R out of range


# --------------------------------------------------------------------------- #
# r9sq.23 -- enum completeness
# --------------------------------------------------------------------------- #
class TestEnumCompleteness:
    @pytest.mark.parametrize(
        "name",
        [
            "iQmass", "iHmolar_idealgas", "iSmolar_idealgas", "iUmolar_idealgas",
            "iHmass_idealgas", "iSmass_idealgas", "iUmass_idealgas",
            "iisentropic_expansion_coefficient", "idalphar_dtau_constdelta",
            "ialpha0", "idalpha0_ddelta_consttau", "id2alpha0_ddelta2_consttau",
            "id3alpha0_ddelta3_consttau",
        ],
    )
    def test_missing_parameters_present(self, name):
        assert hasattr(C, name), f"parameters::{name} not exported"

    @pytest.mark.parametrize(
        "name",
        [
            "QmassT_INPUTS", "PQmass_INPUTS", "QmassSmolar_INPUTS", "QmassSmass_INPUTS",
            "HmolarQmass_INPUTS", "HmassQmass_INPUTS", "DmolarQmass_INPUTS", "DmassQmass_INPUTS",
        ],
    )
    def test_missing_input_pairs_present(self, name):
        assert hasattr(C, name), f"input_pairs::{name} not exported"

    def test_mass_basis_quality_pair_is_usable(self):
        # The whole point of QmassT_INPUTS: drive update() with a mass-basis quality.
        AS = C.AbstractState("HEOS", "Water")
        AS.update(int(C.QmassT_INPUTS), 0.0, 373.0)
        assert AS.Qmass() == pytest.approx(0.0, abs=1e-9)

    def test_fluid_types_enum_present(self):
        for n in ("FLUID_TYPE_PURE", "FLUID_TYPE_PSEUDOPURE", "FLUID_TYPE_UNDEFINED"):
            assert hasattr(C, n), f"fluid_types::{n} not exported"

    def test_fast_evaluate_status_enum_present(self):
        for n in ("fast_evaluate_ok", "fast_evaluate_out_of_range", "fast_evaluate_internal_error"):
            assert hasattr(C, n), f"fast_evaluate_status::{n} not exported"


# --------------------------------------------------------------------------- #
# r9sq.24 -- out-reference / unregistered-return marshalling
# --------------------------------------------------------------------------- #
def _assert_not_marshalling_error(exc):
    """A raised exception is acceptable IFF it is a physics/domain error, NOT the
    binding-marshalling bug (uncallable out-param method / unconvertible return)."""
    msg = str(exc)
    assert not isinstance(exc, TypeError), f"out-param method uncallable (marshalling bug): {msg}"
    assert "Unable to convert" not in msg, f"unregistered return type (marshalling bug): {msg}"


class TestOutParamMarshalling:
    def _water(self):
        AS = C.AbstractState("HEOS", "Water")
        AS.update(int(C.PT_INPUTS), 101325.0, 300.0)
        return AS

    def test_true_critical_point_returns_pair(self):
        AS = C.AbstractState("HEOS", "Water")
        try:
            r = AS.true_critical_point()
        except Exception as e:  # noqa: BLE001
            _assert_not_marshalling_error(e)
            return
        assert len(r) == 2 and all(np.isfinite(x) for x in r)

    def test_viscosity_contributions_returns_dict(self):
        AS = self._water()
        try:
            d = AS.viscosity_contributions()
        except Exception as e:  # noqa: BLE001 -- backend may not implement; that's OK
            _assert_not_marshalling_error(e)
            return
        assert set(d) == {"dilute", "initial_density", "residual", "critical"}

    def test_conductivity_contributions_returns_dict(self):
        AS = self._water()
        try:
            d = AS.conductivity_contributions()
        except Exception as e:  # noqa: BLE001
            _assert_not_marshalling_error(e)
            return
        assert set(d) == {"dilute", "initial_density", "residual", "critical"}

    def test_conformal_state_returns_dict(self):
        AS = self._water()
        try:
            d = AS.conformal_state("Water", 300.0, 1000.0)
        except Exception as e:  # noqa: BLE001 -- ECS may not converge; not a marshalling bug
            _assert_not_marshalling_error(e)
            return
        assert set(d) == {"T", "rhomolar"}

    def test_get_spinodal_data_returns_struct(self):
        assert hasattr(C, "SpinodalData"), "SpinodalData type not registered"
        AS = C.AbstractState("HEOS", "Water")
        try:
            AS.build_spinodal()
            sd = AS.get_spinodal_data()
        except Exception as e:  # noqa: BLE001 -- build may not be supported; not a marshalling bug
            _assert_not_marshalling_error(e)
            return
        for attr in ("tau", "delta", "M1"):
            assert hasattr(sd, attr), f"SpinodalData missing {attr}"

    def test_extract_backend_returns_pair(self):
        backend, fluid = C.extract_backend("REFPROP::Water")
        assert backend == "REFPROP" and fluid == "Water"

    def test_extract_backend_default(self):
        backend, fluid = C.extract_backend("Water")
        assert fluid == "Water"  # backend defaults (e.g. '?')

    def test_extract_fractions_returns_pair(self):
        fluids, fracs = C.extract_fractions("R32[0.7]&R125[0.3]")
        assert fluids == ["R32", "R125"]
        assert list(fracs) == pytest.approx([0.7, 0.3])


# --------------------------------------------------------------------------- #
# r9sq.18/.24 -- legacy Py*-prefixed class aliases (C++ names primary)
# --------------------------------------------------------------------------- #
class TestLegacyClassAliases:
    @pytest.mark.parametrize(
        "py_name,cpp_name",
        [
            ("PyCriticalState", "CriticalState"),
            ("PyGuessesStructure", "GuessesStructure"),
            ("PyPhaseEnvelopeData", "PhaseEnvelopeData"),
            ("PySpinodalData", "SpinodalData"),
        ],
    )
    def test_alias_is_cpp_named_class(self, py_name, cpp_name):
        assert hasattr(C, py_name), f"backwards-compat alias {py_name} missing"
        assert getattr(C, py_name) is getattr(C, cpp_name), f"{py_name} should alias {cpp_name}"

    @pytest.mark.parametrize("py_name", ["PyCriticalState", "PyGuessesStructure"])
    def test_alias_is_constructible(self, py_name):
        # CoolProp.Plots.Common constructs these directly.
        getattr(C, py_name)()


# --------------------------------------------------------------------------- #
# r9sq.16 -- set_reference_state D-form accepts numpy/array numeric args
# --------------------------------------------------------------------------- #
class TestSetReferenceState:
    def test_D_form_accepts_numpy_scalars(self):
        from CoolProp.CoolProp import set_reference_state
        try:
            # numpy scalars must coerce via __float__ (no std::bad_cast).
            set_reference_state("Water", np.float64(300.0), np.float64(1000.0), 100.0, 0.5)
        finally:
            set_reference_state("Water", "DEF")

    def test_D_form_with_propssi_density(self):
        # The exact HighLevelAPI doc pattern: density from PropsSI fed straight in.
        from CoolProp.CoolProp import set_reference_state
        Dmolar = PropsSI("Dmolar", "T", 298.15, "P", 101325, "NH3")
        try:
            set_reference_state("NH3", 298.15, Dmolar, -60000.12, 314.159)
        finally:
            set_reference_state("NH3", "DEF")

    def test_string_form_still_works(self):
        # The 1-arg (string) dispatch must still work (no bad_cast). ASHRAE is only
        # valid for fluids that exist at -40 C, so use R134a for it.
        from CoolProp.CoolProp import set_reference_state
        try:
            set_reference_state("R134a", "ASHRAE")
        finally:
            set_reference_state("R134a", "DEF")
