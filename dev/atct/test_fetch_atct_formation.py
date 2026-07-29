"""Offline tests for the ATcT regenerator. No network access."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_atct_formation import (  # noqa: E402
    AmbiguousSpecies,
    AtctParseError,
    AtctRow,
    parse_atct_rows,
    select_gas_row,
)

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "atct_rows_sample.html"


@pytest.fixture(scope="module")
def rows():
    return parse_atct_rows(FIXTURE.read_text(encoding="utf-8"))


def by_cas(rows, cas):
    return [r for r in rows if r.cas == cas]


def test_parses_a_normal_row(rows):
    methane = by_cas(rows, "74-82-8")
    assert len(methane) == 1
    assert methane[0].name == "Methane"
    assert methane[0].formula == "CH4"
    assert methane[0].phase == "(g)"
    assert methane[0].dhf298_kJ_per_mol == pytest.approx(-74.513)
    assert methane[0].uncertainty_kJ_per_mol == pytest.approx(0.043)
    assert methane[0].atct_id == "74-82-8*0"


def test_element_rows_are_not_dropped(rows):
    """Elements report uncertainty 'exact' with an empty units cell.

    A '±'-anchored regex silently drops all nine elements.
    """
    argon = [r for r in by_cas(rows, "7440-37-1") if r.phase == "(g)"]
    assert len(argon) == 1
    assert argon[0].dhf298_kJ_per_mol == 0.0
    assert argon[0].uncertainty_kJ_per_mol == 0.0


def test_aqueous_row_is_not_selected_as_gas(rows):
    """CAS is not a unique key: argon has (g) and (aq, undissoc) rows."""
    selected = select_gas_row(by_cas(rows, "7440-37-1"))
    assert selected.phase == "(g)"
    assert selected.atct_id == "7440-37-1*0"


def test_exact_gas_row_wins_over_qualified_rows(rows):
    """Water has (g), (l) and spin-qualified rows; the exact (g) row wins."""
    selected = select_gas_row(by_cas(rows, "7732-18-5"))
    assert selected.phase == "(g)"
    assert selected.dhf298_kJ_per_mol == pytest.approx(-241.808)


def test_sole_qualified_gas_row_is_accepted(rows):
    """cis-2-Butene exists only as (g, cis); an exact-(g)-only filter drops it."""
    selected = select_gas_row(by_cas(rows, "590-18-1"))
    assert selected.phase == "(g, cis)"
    assert selected.dhf298_kJ_per_mol == pytest.approx(-6.85)


def test_spin_qualifier_selects_the_right_row(rows):
    selected = select_gas_row(by_cas(rows, "1333-74-0"), qualifier="para")
    assert selected.phase == "(g, para)"
    assert selected.dhf298_kJ_per_mol == pytest.approx(-0.058)
    assert selected.atct_id == "1333-74-0*2"


def test_missing_qualifier_returns_none(rows):
    assert select_gas_row(by_cas(rows, "1333-74-0"), qualifier="nosuch") is None


def test_ambiguous_gas_rows_raise():
    duplicated = [
        AtctRow("111-11-1", "Fake", "XX", "(g)", 1.0, 0.1, "111-11-1*0"),
        AtctRow("111-11-1", "Fake", "XX", "(g)", 2.0, 0.1, "111-11-1*1"),
    ]
    with pytest.raises(AmbiguousSpecies):
        select_gas_row(duplicated)


# --- Fix round 1 -------------------------------------------------------
#
# Regressions for a review pass on Task 1: a formula that itself starts
# with '(' broke the first-'('-anchored phase split (Finding 1); `name`
# skipped HTML-unescaping that every other text field gets (Finding 2);
# a missing name anchor or ATcT-ID cell degraded silently to an empty
# string instead of failing loud (Finding 3); and select_gas_row() never
# checked that its rows actually share one CAS number, the precondition
# its own docstring claims to defend (Finding 4).


def test_paren_leading_formula_is_not_dropped(rows):
    """'(CH3)2NH  (g)' must not be split on its own leading '('.

    A first-'('-anchored split turns this into formula='' and
    phase='(CH3)2NH  (g)', which then matches neither '(g)' nor the
    '(g,' prefix, so select_gas_row() silently reports the species as
    absent.
    """
    dma = by_cas(rows, "124-40-3")
    assert len(dma) == 1
    assert dma[0].name == "Dimethylamine"
    assert dma[0].formula == "(CH3)2NH"
    assert dma[0].phase == "(g)"
    assert dma[0].dhf298_kJ_per_mol == pytest.approx(-16.97)
    assert dma[0].uncertainty_kJ_per_mol == pytest.approx(0.45)
    assert dma[0].atct_id == "124-40-3*0"

    selected = select_gas_row(dma)
    assert selected is not None
    assert selected.phase == "(g)"
    assert selected.formula == "(CH3)2NH"


def test_multi_paren_formula_with_lin_qualifier_is_not_dropped(rows):
    """'(F)(HF)  (g, lin)' exercises a formula with two parenthesized
    groups of its own, plus a previously unseen '(g, lin)' phase
    qualifier, all resolved through the sole-qualified-gas-row branch.
    """
    fhf = by_cas(rows, "12528-21-1")
    assert len(fhf) == 1
    assert fhf[0].formula == "(F)(HF)"
    assert fhf[0].phase == "(g, lin)"

    selected = select_gas_row(fhf)
    assert selected is not None
    assert selected.phase == "(g, lin)"


def _single_row_html(cas, name_cell, formula, atct_id_cell="900-00-0*0"):
    """Build one minimal <tr> in the same cell shape as the fixture rows,
    for exercising a single malformed cell in isolation.
    """
    return (
        '<tr id="s0n0c0 i0 CAS%s">' % cas
        + '<td class="bkgName"><span class="Name"> %s</span></td>' % name_cell
        + '<td class="bkgFormula> <span class="Formula"> '
        '<button onclick="copyname(\'x\');" type="button">%s </button> </span></td>' % formula
        + '<td class="bkgImage"><span class="Name">'
        '<img class="lazy" data-original="images/0.png" alt="X" style="margin:5px" /></span></td>'
        + '<td class="bkgDHf0"><span class="DHf0">0</span></td>'
        + '<td class="bkgDHf298"><span class="DHf298">-1.0</span></td>'
        + '<td class="bkgUncert"><span class="Uncert">&plusmn; 0.10</span></td>'
        + '<td class="bkgUnits"><span class="Units">kJ/mol</span></td>'
        + '<td class="bkgMass"><span class="Mass">1.0 &plusmn;<br />0.1</span></td>'
        + '<td class="bkgATcTID"><span class="ATcTID">%s</span></td>' % atct_id_cell
        + "</tr>"
    )


def test_name_is_html_unescaped():
    """`formula` is explicitly unescaped and the `_span` fields unescape
    internally; `name` must not be the one field left raw.
    """
    row = _single_row_html(
        "900-00-0",
        '<a href="species/?species_number=1">Bis(2-chloroethyl) &amp; ether</a>',
        "XX (g)",
    )
    parsed = parse_atct_rows(row)
    assert parsed[0].name == "Bis(2-chloroethyl) & ether"


def test_missing_name_anchor_raises():
    """Without the species_number anchor, a row must fail loud rather than
    silently yielding a blank-named row.
    """
    row = _single_row_html("900-00-1", "NoAnchorHere", "XX (g)")
    with pytest.raises(AtctParseError):
        parse_atct_rows(row)


def test_missing_atct_id_raises():
    """An empty ATcT-ID cell must fail loud rather than silently yielding a
    row with a blank atct_id.
    """
    row = _single_row_html(
        "900-00-2",
        '<a href="species/?species_number=1">Placeholder</a>',
        "XX (g)",
        atct_id_cell="",
    )
    with pytest.raises(AtctParseError):
        parse_atct_rows(row)


def test_select_gas_row_rejects_mixed_cas():
    """select_gas_row() documents that its rows share one CAS number but
    never checked it -- exactly the "CAS is not a unique key" trap this
    function otherwise defends against.
    """
    mixed = [
        AtctRow("111-11-1", "Fake", "XX", "(g)", 1.0, 0.1, "111-11-1*0"),
        AtctRow("222-22-2", "Other", "YY", "(g)", 2.0, 0.1, "222-22-2*0"),
    ]
    with pytest.raises(AtctParseError):
        select_gas_row(mixed)
