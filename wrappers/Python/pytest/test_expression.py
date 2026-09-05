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
import re
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


# Every thermodynamic name a block reads must be declared in "state_variables";
# anything undeclared is the author's own.  The suite derives the declaration from
# the formula so each test can say what it is actually about, but the contract
# itself is pinned explicitly in test_state_variables_* below.
_STATE = ("T", "P", "Dmolar", "Dmass", "molar_mass", "Smolar_residual", "Bvirial", "dBvirial_dT")


def block(formula, constants=None, arrays=None, state=None):
    """Build the JSON text of a `"type": "expression"` transport block."""
    b = {"type": "expression", "formula": formula}
    if state is None:
        local = set(re.findall(r"\blet\s+(\w+)", formula)) | set(constants or ()) | set(arrays or ())
        state = []
        for m in re.finditer(r"\b([A-Za-z_]\w*)\b(\s*\[)?", formula):
            tok = m.group(1)
            if m.group(2) or tok in local or tok not in _STATE or tok in state:
                continue
            state.append(tok)
    # `state=[]` must emit an EXPLICIT empty declaration, not omit the key -- those
    # are different inputs to the parser and both need exercising.
    if state is not None:
        b["state_variables"] = list(state)
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
    e = Expression(block("P*2 + T"))
    # Reported in first-reference order (NOT declaration order); only the canonical
    # spelling is accepted, so it round-trips whatever the author wrote.
    assert e.required_inputs() == ["P", "T"]
    assert Expression(block("Dmass*molar_mass")).required_inputs() == ["Dmass", "molar_mass"]


def test_constants_and_arrays_and_let_and_sum():
    e = Expression(
        block(
            "let delta = Dmolar/rhomolar_reduce\nsum(i: a[i]*delta^d[i])",
            constants={"rhomolar_reduce": 10000.0},
            arrays={"a": [2.0, 3.0], "d": [1.0, 2.0]},
        )
    )
    assert e.required_inputs() == ["Dmolar"]
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


def test_a_constant_colliding_with_a_declared_state_variable_is_rejected():
    # A declared name resolves before constants, so a constant of that name would be
    # dead data and the formula would quietly mean the state value.  Hard error --
    # but only when the block actually declared it.
    for name in ("T", "Dmolar", "Dmass", "molar_mass", "P"):
        with pytest.raises(ValueError):
            Expression(block(name, constants={name: 1.0}, state=[name]))
        # Undeclared, the very same constant is simply the author's own number.
        assert Expression(block(name, constants={name: 1.0}, state=[])).evaluate(state()) == 1.0
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
    assert Expression(block("Dmolar")).evaluate(state(rhomolar=1234.0)) == pytest.approx(1234.0)


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
    assert Expression(block("P")).evaluate(AS) == pytest.approx(
        CP.PropsSI("P", "T", T, "Dmolar", rhomolar, "R123"), rel=1e-12
    )
    # rhomass is the EOS's own mass density at that state, not a re-derivation.
    assert Expression(block("Dmass")).evaluate(AS) == pytest.approx(
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
    assert Expression(block("molar_mass*Dmolar")).evaluate(AS) == pytest.approx(M_mix * AS.rhomolar(), rel=1e-12)


def test_evaluating_against_an_unset_state_raises():
    # AbstractState.p() on a fresh state hands back -inf rather than raising, and
    # the formula would propagate it into a plausible-looking answer.
    e = Expression(block("P*2"))
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
    e = Expression(block("T*Dmolar"))
    AS = CP.AbstractState("HEOS", "R123")
    out = []
    for T in (300.0, 400.0):
        AS.update(CP.DmolarT_INPUTS, 100.0, T)
        out.append(e.evaluate(AS))
    assert out == pytest.approx([30000.0, 40000.0])


# --------------------------------------------------------------------------- #
# Declared state variables
# --------------------------------------------------------------------------- #
def test_state_variables_use_coolprops_own_names():
    """The DSL keeps no vocabulary of its own; the names are CoolProp's."""
    assert Expression(block("T", state=["T"])).required_inputs() == ["T"]
    assert Expression(block("Dmolar", state=["Dmolar"])).required_inputs() == ["Dmolar"]
    # The DSL once spelled these `rhomolar`/`rhomass`/`p`.  It no longer invents
    # names -- and the invented lowercase `p` is what used to collide with the
    # exponent array every viscosity paper writes as p_i.
    for legacy in ("rhomolar", "rhomass", "p"):
        with pytest.raises(ValueError):
            Expression(block(legacy, state=[legacy]))


def test_state_variables_are_opt_in_so_undeclared_names_are_yours():
    """A block that never asks for pressure keeps `p` for its own coefficients."""
    e = Expression(block("sum(i: n[i]*p[i])", arrays={"n": [2.0, 3.0], "p": [5.0, 7.0]}, state=[]))
    assert e.required_inputs() == []
    assert e.evaluate(state()) == pytest.approx(2.0 * 5.0 + 3.0 * 7.0)
    # Declare pressure and the same name is a collision, not a silent shadowing.
    with pytest.raises(ValueError):
        Expression(block("sum(i: n[i]*p[i])", arrays={"n": [2.0], "p": [5.0]}, state=["P"]))


def test_reading_an_undeclared_quantity_names_the_fix():
    with pytest.raises(ValueError, match="state_variables"):
        Expression(block("Bvirial", state=[]))


def test_two_classes_stay_refused_even_when_declared():
    """Opt-in is not a licence: these cannot be honoured whatever the author asks."""
    for nm in ("V", "viscosity", "conductivity"):  # re-enters the correlation
        with pytest.raises(ValueError):
            Expression(block(nm, state=[nm]))
    for nm in ("T_critical", "rhomolar_critical", "T_reducing"):  # config-dependent
        with pytest.raises(ValueError):
            Expression(block(nm, state=[nm]))
        # ...but each stays usable as an ordinary frozen constant.
        assert Expression(block(nm, constants={nm: 7.0}, state=[])).evaluate(state()) == 7.0


def test_let_renames_a_state_variable_to_whatever_reads_best():
    """Dropping the invented aliases costs nothing: renaming is in the language."""
    e = Expression(block("let rho = Dmolar\nlet tau = Tc/T\nrho*tau",
                         constants={"Tc": 456.83}, state=["Dmolar", "T"]))
    assert e.required_inputs() == ["Dmolar", "T"]
    assert e.evaluate(state(rhomolar=1e4, T=300.0)) == pytest.approx(1e4 * 456.83 / 300.0)
    # A `let` cannot silently shadow a DECLARED name: the declaration then goes
    # unread, which is already an error.
    with pytest.raises(ValueError):
        Expression(block("let T = 5\nT", state=["T"]))
    # Undeclared, the same name is simply a local -- no state involved at all.
    assert Expression(block("let T = 5\nT", state=[])).evaluate(state()) == 5.0


def test_only_the_exposed_set_is_declarable():
    """Opt-in: the DSL exposes an explicit set and refuses everything else."""
    EXPOSED = ("T", "P", "Dmolar", "Dmass", "molar_mass", "Smolar_residual", "Bvirial", "dBvirial_dT")
    for nm in EXPOSED:
        assert Expression(block(nm, state=[nm])).required_inputs() == [nm]
    # CoolProp's back-compat aliases are a trap in a formula full of coefficient
    # names -- `A` is the speed of sound, `D` is Dmass, `M` is molar_mass -- and they
    # are refused for free, because the exposed set holds canonical names only.
    for alias in ("A", "C", "D", "G", "M", "O", "S", "U", "DMOLAR"):
        with pytest.raises(ValueError, match="not a state variable the DSL exposes"):
            Expression(block(alias, state=[alias]))
    # ...so each stays free for the author, which is the point.
    assert Expression(block("A*M", constants={"A": 3.0, "M": 4.0}, state=[])).evaluate(state()) == 12.0
    # The refusal names the whole set, so it answers "then what do I write?" --
    # including for this DSL's own retired spellings.
    for old, now in (("rhomolar", "Dmolar"), ("rhomass", "Dmass"), ("p", "P")):
        with pytest.raises(ValueError, match=now):
            Expression(block(old, state=[old]))


def test_tau_and_delta_are_the_reducing_state_by_another_name():
    """keyed_output computes them as _reducing.T/_T and _rhomolar/_reducing.rhomolar."""
    for nm in ("Tau", "Delta", "T_reducing", "rhomolar_reducing"):
        with pytest.raises(ValueError):
            Expression(block(nm, state=[nm]))


def test_non_state_quantities_are_refused():
    for nm in ("Phase", "Q", "Qmass"):          # enum ordinal / sentinel-valued
        with pytest.raises(ValueError):
            Expression(block(nm, state=[nm]))
    for nm in ("T_freeze", "GWP100", "ODP"):    # fluid metadata, not state
        with pytest.raises(ValueError):
            Expression(block(nm, state=[nm]))

