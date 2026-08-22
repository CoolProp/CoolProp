"""
pytest suite for the runtime transport-property expression DSL surfaced as
``CoolProp.CoolProp.Expression``.

Runs against an *installed* CoolProp built with ``-DCOOLPROP_NANOBIND=ON`` (see
README.md in this directory); the class does not exist on the legacy Cython
wheel, so the whole module skips there.

The point of this surface is authoring: you paste the same ``"type":
"expression"`` block you would put in a fluid JSON file, evaluate it against an
``AbstractState`` you set up the usual way, and see the number -- no rebuild, no
fluid-file edit.  These tests pin that the block text is accepted verbatim, that
the reported inputs match the formula, and that the numbers agree with what the
same algebra gives through ``PropsSI``.
"""
import json
import math

import pytest

from CoolProp import CoolProp as CP

# Probe the BUILD, not the symbol: skipping on `Expression is None` would turn a
# dropped/renamed binding into a green skip.  `Props` was removed in v8, so its
# presence identifies the legacy Cython wheel -- the only build that legitimately
# lacks Expression.  (Same idiom as test_nanobind_interface.py.)
is_legacy = hasattr(CP, "Props")
pytestmark = pytest.mark.skipif(is_legacy, reason="v8 (nanobind) build only")

Expression = getattr(CP, "Expression", None)


def state(fluid="R123", rhomolar=1e4, T=300.0):
    """An AbstractState sat at (rhomolar, T) -- the caller's job, not the block's."""
    AS = CP.AbstractState("HEOS", fluid)
    AS.update(CP.DmolarT_INPUTS, rhomolar, T)
    return AS


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
def test_constant_only_formula_reads_no_state():
    e = Expression(block("2 + 3*4"))
    assert e.required_inputs() == []
    assert e.evaluate(state()) == pytest.approx(14.0)


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
    assert e.evaluate(state(rhomolar=5000.0)) == pytest.approx(2.0 * 0.5 + 3.0 * 0.25)


def test_the_binding_exists_on_a_v8_build():
    # Guards the skip above: on a nanobind build the class must be there.
    assert Expression is not None


def test_bad_formula_raises_at_construction():
    with pytest.raises(ValueError):
        Expression(block("2 +"))
    with pytest.raises(ValueError):
        Expression(block("nonexistent_name"))
    with pytest.raises(ValueError):
        Expression('{"type": "expression"}')  # no "formula"


def test_a_constant_colliding_with_an_input_is_rejected():
    # Inputs resolve before constants, so a constant named `T` would be dead data
    # and the formula would quietly mean the state temperature.  Hard error.
    for name in ("T", "rhomolar", "rhomass", "molar_mass", "p"):
        with pytest.raises(ValueError):
            Expression(block("1", constants={name: 1.0}))
    Expression(block("T_reduce", constants={"T_reduce": 132.0}))  # near-miss is fine


def test_a_transport_output_name_does_not_resolve():
    # `V` is a valid PropsSI output; resolving it inside a transport formula would
    # re-enter the correlation being defined.  It must be a compile error, not a
    # recursion at eval time.
    with pytest.raises(ValueError):
        Expression(block("V"))


# --------------------------------------------------------------------------- #
# Evaluating against a fluid
# --------------------------------------------------------------------------- #
def test_state_variables_come_from_the_AbstractState():
    assert Expression(block("T")).evaluate(state(T=321.5)) == pytest.approx(321.5)
    assert Expression(block("rhomolar")).evaluate(state(rhomolar=1234.0)) == pytest.approx(1234.0)


def test_the_block_follows_the_state_object():
    # Re-set the state through the ordinary API and the same compiled block tracks
    # it -- including through an input pair the block knows nothing about.
    e = Expression(block("T"))
    AS = CP.AbstractState("HEOS", "R123")
    AS.update(CP.DmolarT_INPUTS, 1e4, 300.0)
    assert e.evaluate(AS) == pytest.approx(300.0)
    AS.update(CP.PT_INPUTS, 1e5, 350.0)
    assert e.evaluate(AS) == pytest.approx(350.0)


def test_molar_mass_and_pressure_come_from_the_EOS():
    T, rhomolar = 400.0, 2000.0
    AS = state(T=T, rhomolar=rhomolar)
    assert Expression(block("molar_mass")).evaluate(AS) == pytest.approx(CP.PropsSI("molar_mass", "R123"))
    assert Expression(block("p")).evaluate(AS) == pytest.approx(
        CP.PropsSI("P", "T", T, "Dmolar", rhomolar, "R123"), rel=1e-12
    )
    # rhomass is the EOS's own mass density at that state, not a re-derivation.
    assert Expression(block("rhomass")).evaluate(AS) == pytest.approx(
        CP.PropsSI("Dmass", "T", T, "Dmolar", rhomolar, "R123"), rel=1e-12
    )


def test_a_mixture_state_is_just_another_state():
    # Nothing in the block is pure-fluid-specific -- the caller owns the state, so
    # a mixture composition is simply what the formula reads.
    AS = CP.AbstractState("HEOS", "R32&R125")
    AS.set_mole_fractions([0.4, 0.6])
    AS.update(CP.PT_INPUTS, 1e5, 350.0)  # HEOS mixtures do not take DmolarT
    # NOT compared against AS.rhomass(): that IS rhomolar*molar_mass internally, so
    # it would hold even if the mixture molar mass were wrong.  Rebuild it from the
    # pure-component values.
    M_mix = 0.4 * CP.PropsSI("molar_mass", "R32") + 0.6 * CP.PropsSI("molar_mass", "R125")
    assert AS.rhomolar() > 0.0
    assert Expression(block("molar_mass*rhomolar")).evaluate(AS) == pytest.approx(M_mix * AS.rhomolar(), rel=1e-12)


def test_evaluating_against_an_unset_state_raises():
    # AbstractState.p() on a fresh state hands back -inf rather than raising, and
    # the formula would propagate it into a plausible-looking answer.
    e = Expression(block("p*2"))
    with pytest.raises(ValueError, match="finite"):
        e.evaluate(CP.AbstractState("HEOS", "R123"))  # never update()d


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
    AS = CP.AbstractState("HEOS", "R123")
    for T in (250.0, 300.0, 400.0, 500.0):
        AS.update(CP.DmolarT_INPUTS, 0.1, T)
        expected = sum(a * T ** t for a, t in zip(dilute["a"], dilute["t"]))
        got = e.evaluate(AS)
        assert math.isfinite(got)
        assert got == pytest.approx(expected, rel=1e-12)


def test_a_compiled_block_is_reusable_across_states():
    e = Expression(block("T*rhomolar"))
    AS = CP.AbstractState("HEOS", "R123")
    out = []
    for T in (300.0, 400.0):
        AS.update(CP.DmolarT_INPUTS, 100.0, T)
        out.append(e.evaluate(AS))
    assert out == pytest.approx([30000.0, 40000.0])
