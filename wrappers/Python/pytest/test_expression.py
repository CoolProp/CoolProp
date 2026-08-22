"""
pytest suite for the runtime transport-property expression DSL surfaced as
``CoolProp.CoolProp.Expression``.

Runs against an *installed* CoolProp built with ``-DCOOLPROP_NANOBIND=ON`` (see
README.md in this directory); the class does not exist on the legacy Cython
wheel, so the whole module skips there.

The point of this surface is authoring: you paste the same ``"type":
"expression"`` block you would put in a fluid JSON file, evaluate it at a state,
and see the number -- no rebuild, no fluid-file edit.  These tests pin that the
block text is accepted verbatim, that the reported inputs match the formula, and
that the numbers agree with what the same algebra gives through ``PropsSI``.
"""
import json
import math

import pytest

from CoolProp import CoolProp as CP

Expression = getattr(CP, "Expression", None)
pytestmark = pytest.mark.skipif(Expression is None, reason="v8 (nanobind) build only")


def block(formula, constants=None, arrays=None):
    """Build the JSON text of a `"type": "expression"` transport block."""
    b = {"type": "expression", "formula": formula}
    if constants:
        b["constants"] = constants
    if arrays:
        b["arrays"] = arrays
    return json.dumps(b)


# --------------------------------------------------------------------------- #
# Compiling and inspecting
# --------------------------------------------------------------------------- #
def test_constant_only_formula_needs_no_fluid():
    e = Expression(block("2 + 3*4"))
    assert e.required_inputs() == []
    assert e.evaluate(300.0, 1e4) == pytest.approx(14.0)


def test_required_inputs_reports_the_names_the_formula_uses():
    e = Expression(block("p*2 + T"))
    # Reported in first-reference order, spelled the way the formula spells them.
    assert e.required_inputs() == ["p", "T"]
    assert Expression(block("rhomass*molar_mass")).required_inputs() == ["rhomass", "molar_mass"]


def test_constants_and_arrays_and_let_and_sum():
    e = Expression(
        block(
            "let delta = rhomolar/rhomolar_reduce\nsum(i: a[i]*delta^d[i])",
            constants={"rhomolar_reduce": 10000.0},
            arrays={"a": [2.0, 3.0], "d": [1.0, 2.0]},
        )
    )
    assert e.required_inputs() == ["rhomolar"]
    assert e.evaluate(300.0, 5000.0) == pytest.approx(2.0 * 0.5 + 3.0 * 0.25)


def test_bad_formula_raises_at_construction():
    with pytest.raises(ValueError):
        Expression(block("2 +"))
    with pytest.raises(ValueError):
        Expression(block("nonexistent_name"))
    with pytest.raises(ValueError):
        Expression('{"type": "expression"}')  # no "formula"


def test_a_transport_output_name_does_not_resolve():
    # `V` is a valid PropsSI output; resolving it inside a transport formula would
    # re-enter the correlation being defined.  It must be a compile error, not a
    # recursion at eval time.
    with pytest.raises(ValueError):
        Expression(block("V"))


# --------------------------------------------------------------------------- #
# Evaluating against a fluid
# --------------------------------------------------------------------------- #
def test_state_variables_come_from_the_fluid_state():
    e = Expression(block("T"))
    assert e.evaluate(321.5, 1e4, "R123") == pytest.approx(321.5)
    e = Expression(block("rhomolar"))
    assert e.evaluate(321.5, 1234.0, "R123") == pytest.approx(1234.0)


def test_molar_mass_and_pressure_come_from_the_EOS():
    T, rhomolar = 400.0, 2000.0
    assert Expression(block("molar_mass")).evaluate(T, rhomolar, "R123") == pytest.approx(
        CP.PropsSI("molar_mass", "R123")
    )
    assert Expression(block("p")).evaluate(T, rhomolar, "R123") == pytest.approx(
        CP.PropsSI("P", "T", T, "Dmolar", rhomolar, "R123"), rel=1e-12
    )
    # rhomass is the EOS's own mass density at that state, not a re-derivation.
    assert Expression(block("rhomass")).evaluate(T, rhomolar, "R123") == pytest.approx(
        CP.PropsSI("Dmass", "T", T, "Dmolar", rhomolar, "R123"), rel=1e-12
    )


def test_a_fluid_is_required_only_when_the_formula_needs_one():
    Expression(block("T*rhomolar")).evaluate(300.0, 1e4)  # fine, no fluid
    with pytest.raises(ValueError):
        Expression(block("p")).evaluate(300.0, 1e4)


def test_backend_qualified_fluid_name():
    e = Expression(block("p"))
    assert e.evaluate(400.0, 2000.0, "HEOS::R123") == pytest.approx(
        e.evaluate(400.0, 2000.0, "R123"), rel=1e-12
    )


def test_reproduces_the_hardcoded_dilute_viscosity_of_R123():
    """R123's dilute viscosity is the powers_of_T form: sum a[i]*T^t[i].

    Lift its coefficients straight out of the fluid's own JSON, express the same
    algebra as a DSL block, and check the DSL agrees -- the authoring workflow
    this class exists for, end to end.
    """
    dilute = json.loads(CP.get_fluid_param_string("R123", "JSON"))[0]["TRANSPORT"]["viscosity"]["dilute"]
    assert dilute["type"] == "powers_of_T"

    e = Expression(block("sum(i: a[i]*T^t[i])", arrays={"a": dilute["a"], "t": dilute["t"]}))
    assert e.required_inputs() == ["T"]
    for T in (250.0, 300.0, 400.0, 500.0):
        expected = sum(a * T ** t for a, t in zip(dilute["a"], dilute["t"]))
        got = e.evaluate(T, 0.1, "R123")
        assert math.isfinite(got)
        assert got == pytest.approx(expected, rel=1e-12)


def test_a_compiled_block_is_reusable_across_states():
    e = Expression(block("T*rhomolar"))
    assert [e.evaluate(T, 100.0) for T in (300.0, 400.0)] == pytest.approx([30000.0, 40000.0])
