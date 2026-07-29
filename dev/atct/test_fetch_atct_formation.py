"""Offline tests for the ATcT regenerator. No network access."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_atct_formation import (  # noqa: E402
    AmbiguousSpecies,
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
