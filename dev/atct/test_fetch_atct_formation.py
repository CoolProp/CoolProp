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
from fetch_atct_formation import (  # noqa: E402
    FluidRef,
    bind,
    clear_standard_state,
    load_coolprop_fluids,
    normalize_cas,
)

FLUIDS_DIR = Path(__file__).resolve().parents[2] / "dev" / "fluids"

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
#
# Fix round 2 (found running the real 3442-row page in Task 3): Finding 3's
# "missing name anchor" half was itself wrong.  ATcT's own documentation
# says the Species Name is a link "for species other than elements in
# their standard states" -- so Boron, Sulfur, Br2, I2 and the aqueous
# Hydron ion render their name as plain text with no <a> at all.  That is
# a valid row, not a data-integrity failure; requiring the anchor was too
# narrow an extraction, not a fail-loud guard worth keeping as originally
# written.  The ATcT-ID half of Finding 3 was correctly calibrated (it
# never fired on the full page) and is unchanged.


def test_unlinked_element_row_from_the_real_page_is_not_dropped(rows):
    """Boron's real row from the live ATcT 1.220 page (verbatim, not
    synthesized): a condensed reference-state element with a plain-text
    name and no species-detail link.  Must parse, and must not be
    selectable as a gas row -- it is (cr,l), not (g).
    """
    boron = by_cas(rows, "7440-42-8")
    assert len(boron) == 1
    assert boron[0].name == "Boron"
    assert boron[0].formula == "B"
    assert boron[0].phase == "(cr,l)"
    assert boron[0].dhf298_kJ_per_mol == pytest.approx(0.0)
    assert boron[0].uncertainty_kJ_per_mol == pytest.approx(0.0)
    assert boron[0].atct_id == "7440-42-8*500"

    assert select_gas_row(boron) is None


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


def test_unlinked_name_parses_to_plain_text():
    """ATcT's own page documentation: 'For species other than elements in
    their standard states, the Species Name acts as a link' -- so an
    element in its standard state (e.g. Boron, Sulfur) has NO <a> anchor
    around its name, just plain text.  That is a valid row, not an error;
    an earlier version of this parser required the anchor and raised here,
    which was too narrow an extraction, not a real data-integrity problem.
    """
    row = _single_row_html("900-00-1", "PlainTextName", "XX (g)")
    parsed = parse_atct_rows(row)
    assert len(parsed) == 1
    assert parsed[0].name == "PlainTextName"


def test_empty_name_cell_raises():
    """A row with a genuinely empty name cell must still fail loud."""
    row = _single_row_html("900-00-1", "", "XX (g)")
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


# --- Task 2: CoolProp <-> ATcT binder -----------------------------------


def test_normalize_cas_passes_plain_cas_through():
    assert normalize_cas("74-82-8") == ("74-82-8", None)


def test_normalize_cas_splits_spin_isomer_suffixes():
    """CoolProp writes ParaHydrogen's CAS as '1333-74-0p'."""
    assert normalize_cas("1333-74-0p") == ("1333-74-0", "para")
    assert normalize_cas("1333-74-0o") == ("1333-74-0", "ortho")


def test_load_coolprop_fluids_reads_the_real_tree():
    fluids = load_coolprop_fluids(FLUIDS_DIR)
    assert fluids["Methane"].cas == "74-82-8"
    assert fluids["ParaHydrogen"].cas == "1333-74-0p"
    assert len(fluids) > 100


def test_bind_matches_and_reports_absent(rows):
    fluids = {
        "Methane": FluidRef("Methane", "74-82-8", Path("Methane.json")),
        "ParaHydrogen": FluidRef("ParaHydrogen", "1333-74-0p", Path("ParaHydrogen.json")),
        "R134a": FluidRef("R134a", "811-97-2", Path("R134a.json")),
    }
    result = bind(fluids, rows)
    assert result.errors == []
    assert result.matched["Methane"].atct_id == "74-82-8*0"
    # ATcT publishes a real spin-resolved value; it is NOT zero.
    assert result.matched["ParaHydrogen"].dhf298_kJ_per_mol == pytest.approx(-0.058)
    assert "R134a" in result.absent
    assert "811-97-2" in result.absent["R134a"]


def test_bind_records_a_fluid_with_no_cas_as_absent(rows):
    fluids = {"Air": FluidRef("Air", "", Path("Air.json"))}
    result = bind(fluids, rows)
    assert "Air" in result.absent
    assert result.errors == []


def test_bind_collects_ambiguity_as_an_error():
    fluids = {"Fake": FluidRef("Fake", "111-11-1", Path("Fake.json"))}
    duplicated = [
        AtctRow("111-11-1", "Fake", "XX", "(g)", 1.0, 0.1, "111-11-1*0"),
        AtctRow("111-11-1", "Fake", "XX", "(g)", 2.0, 0.1, "111-11-1*1"),
    ]
    result = bind(fluids, duplicated)
    assert result.matched == {}
    assert len(result.errors) == 1
    assert "Fake" in result.errors[0]


# --- Task 3: writer, coverage ledger, CLI --------------------------------

from fetch_atct_formation import (  # noqa: E402
    compare_ledger,
    coverage_ledger,
    main,
    standard_state_block,
    write_standard_state,
)


def test_standard_state_block_converts_kJ_to_J(rows):
    methane = [r for r in rows if r.cas == "74-82-8"][0]
    block = standard_state_block(methane, "1.220")
    assert block["T"] == 298.15
    assert block["p"] == 100000.0
    assert block["phase"] == "ideal_gas"
    assert block["hmolar_formation"]["value"] == pytest.approx(-74513.0)
    assert block["hmolar_formation"]["uncertainty"] == pytest.approx(43.0)
    assert block["hmolar_formation"]["units"] == "J/mol"
    assert block["hmolar_formation"]["source"] == "ATcT"
    assert block["hmolar_formation"]["version"] == "1.220"
    assert block["hmolar_formation"]["id"] == "74-82-8*0"


def test_write_standard_state_preserves_canonical_formatting(tmp_path, rows):
    """The writer must not reformat the file; only the new block may appear.

    `original` is read from the real dev/fluids/Methane.json, which this
    same tool's --write step regenerates in place -- so once that has run
    (in CI, or in any workspace after this feature merges), the "before"
    snapshot already carries INFO.STANDARD_STATE.  Stripping it from
    `before` too (it is a no-op when the block is not yet present) keeps
    this test asserting its real invariant -- nothing else in the file
    changes -- independent of whether Methane.json has already been
    regenerated once or a hundred times.
    """
    import json as _json

    original = FLUIDS_DIR / "Methane.json"
    scratch = tmp_path / "Methane.json"
    scratch.write_text(original.read_text(encoding="utf-8"), encoding="utf-8")

    methane = [r for r in rows if r.cas == "74-82-8"][0]
    write_standard_state(scratch, methane, "1.220")

    before = _json.loads(original.read_text(encoding="utf-8"))
    before.get("INFO", {}).pop("STANDARD_STATE", None)
    after = _json.loads(scratch.read_text(encoding="utf-8"))
    assert "STANDARD_STATE" in after["INFO"]
    del after["INFO"]["STANDARD_STATE"]
    assert after == before
    # 2-space indent, no trailing newline.
    assert not scratch.read_text(encoding="utf-8").endswith("\n")


def test_insert_then_clear_restores_the_file_byte_for_byte(tmp_path, rows):
    """Insert-then-remove must be a textual no-op, not merely a parse no-op.

    The parsed-equality assertion above cannot see reordering: `after ==
    before` holds however the keys are permuted.  That is exactly how
    json_options' sort_keys=True went unnoticed -- dev/fluids/*.json has
    non-canonical ordering inside SUPERANCILLARY (injected by
    dev/scripts/inject_superancillary.py, which writes indent=2 with no
    sort_keys), so re-serializing sorted rewrote 2608 lines across 76 files
    while every parsed-equality test stayed green.

    Comparing the TEXT closes that hole: with a sorting writer this fails on
    the first fluid whose stored order is not already sorted.
    """
    # R134a is genuinely absent from ATcT, so its committed file never carries
    # a block -- an unambiguous baseline no matter how often the regenerator
    # has run.  Its SUPERANCILLARY ordering is non-canonical like almost every
    # other fluid, which is what gives this test its teeth.
    original = FLUIDS_DIR / "R134a.json"
    scratch = tmp_path / "R134a.json"
    before = original.read_text(encoding="utf-8")
    scratch.write_text(before, encoding="utf-8")
    assert "STANDARD_STATE" not in before

    methane = [r for r in rows if r.cas == "74-82-8"][0]
    write_standard_state(scratch, methane, "1.220")
    assert scratch.read_text(encoding="utf-8") != before  # the block really landed

    assert clear_standard_state(scratch) is True
    assert scratch.read_text(encoding="utf-8") == before


def test_compare_ledger_flags_a_dropped_fluid():
    expected = {"Methane": {"state": "matched", "atct_id": "74-82-8*0"}}
    actual = {"Methane": {"state": "absent", "reason": "no ATcT species for CAS 74-82-8"}}
    problems = compare_ledger(expected, actual)
    assert len(problems) == 1
    assert "Methane" in problems[0]


def test_compare_ledger_is_quiet_when_identical():
    ledger = {"Methane": {"state": "matched", "atct_id": "74-82-8*0"}}
    assert compare_ledger(ledger, dict(ledger)) == []


def test_coverage_ledger_records_both_states(rows):
    fluids = {
        "Methane": FluidRef("Methane", "74-82-8", Path("Methane.json")),
        "R134a": FluidRef("R134a", "811-97-2", Path("R134a.json")),
    }
    ledger = coverage_ledger(bind(fluids, rows))
    assert ledger["Methane"] == {"state": "matched", "atct_id": "74-82-8*0"}
    assert ledger["R134a"]["state"] == "absent"


# --- Fix round 1 (post-merge review) ------------------------------------
#
# A reviewer corrupted a committed ledger entry (flipped Methane's state to
# "absent") and ran --update-ledger: exit 0, the corruption silently
# overwritten, zero console output about what changed.  --write alone
# correctly gates on compare_ledger(), but --update-ledger is exactly the
# brief's own Step 5 bootstrap invocation (paired with --write in one shot)
# and exactly what a future ATcT version bump will reach for -- so a fluid
# quietly flipping matched -> absent between versions would leave no trace
# in the console, the PR diff, or CI logs.  main() now diffs the existing
# ledger against the freshly computed one BEFORE overwriting it and prints
# every difference as labeled informational output; it still exits 0 and
# still writes, since --update-ledger's whole point is to accept the new
# state, but a human (or CI log) reviewing the run now has something to
# actually read.


def test_update_ledger_reports_changes_against_a_corrupted_ledger(tmp_path, monkeypatch, capsys):
    """--update-ledger must surface, not silently swallow, what it overwrites.

    main() resolves its ledger path from its own module __file__, not a
    CLI flag, so isolating it from the real committed
    dev/atct/expected_coverage.json means monkeypatching __file__ to a
    scratch location -- Path(...).resolve() does not require the path to
    exist, so this is a clean way to run main() fully offline against a
    disposable ledger.
    """
    import json as _json

    # Reach the module object via sys.modules rather than a second
    # "import fetch_atct_formation": importing it both ways trips
    # CodeQL py/import-and-import-from, and the module is already
    # loaded by the "from fetch_atct_formation import ..." above.
    module = sys.modules["fetch_atct_formation"]

    monkeypatch.setattr(module, "__file__", str(tmp_path / "fetch_atct_formation.py"))

    fluids_dir = tmp_path / "fluids"
    fluids_dir.mkdir()
    (fluids_dir / "Methane.json").write_text(
        _json.dumps({"INFO": {"NAME": "Methane", "CAS": "74-82-8"}}), encoding="utf-8"
    )

    # Pre-existing ledger, corrupted: Methane recorded as absent when the
    # page (FIXTURE) actually has a clean gas-phase row for it.
    ledger_path = tmp_path / "expected_coverage.json"
    ledger_path.write_text(
        _json.dumps({"Methane": {"state": "absent", "reason": "corrupted by hand"}}),
        encoding="utf-8",
    )

    exit_code = main(
        [
            "--version", "1.220",
            "--cache", str(FIXTURE),
            "--fluids-dir", str(fluids_dir),
            "--update-ledger",
        ]
    )
    out = capsys.readouterr().out

    assert exit_code == 0
    assert "INFO:" in out
    assert "Methane" in out
    # The old (corrupted) and new (correct) states must both be visible in
    # the printed diff, not just a bare "something changed".
    assert "absent" in out
    assert "matched" in out

    new_ledger = _json.loads(ledger_path.read_text(encoding="utf-8"))
    assert new_ledger["Methane"]["state"] == "matched"


def test_update_ledger_with_no_existing_ledger_says_so(tmp_path, monkeypatch, capsys):
    """A bootstrap run (no ledger file yet) must say so explicitly, not just
    silently write one -- otherwise "first ever run" and "corrupted entry
    silently overwritten" look identical in the console output.
    """
    import json as _json

    # Reach the module object via sys.modules rather than a second
    # "import fetch_atct_formation": importing it both ways trips
    # CodeQL py/import-and-import-from, and the module is already
    # loaded by the "from fetch_atct_formation import ..." above.
    module = sys.modules["fetch_atct_formation"]

    monkeypatch.setattr(module, "__file__", str(tmp_path / "fetch_atct_formation.py"))

    fluids_dir = tmp_path / "fluids"
    fluids_dir.mkdir()
    (fluids_dir / "Methane.json").write_text(
        _json.dumps({"INFO": {"NAME": "Methane", "CAS": "74-82-8"}}), encoding="utf-8"
    )
    # No expected_coverage.json written in tmp_path -- genuinely first run.

    exit_code = main(
        [
            "--version", "1.220",
            "--cache", str(FIXTURE),
            "--fluids-dir", str(fluids_dir),
            "--update-ledger",
        ]
    )
    out = capsys.readouterr().out

    assert exit_code == 0
    assert "no existing ledger" in out
    assert (tmp_path / "expected_coverage.json").exists()


def test_measured_row_with_no_units_cell_raises():
    """An empty units cell is legitimate only on the 'exact' element rows.

    `_span` returns "" both for a genuinely empty cell and for a cell whose
    span is missing outright, so accepting "" unconditionally made the units
    check pass for every row on exactly the page change most likely to
    accompany a units change -- a restyle that renames the CSS class.  The
    values would then be ingested as kJ/mol whatever the page actually said,
    and neither the coverage ledger (state and ID unchanged) nor the C++
    plausibility test (|dHf| < 2000 kJ/mol) would notice.
    """
    row = _single_row_html("900-00-0", "<a href=\'#\'>Fakium</a>", "Fk (g)").replace(
        '<span class="Units">kJ/mol</span>', '<span class="UnitsRenamed">kcal/mol</span>'
    )
    with pytest.raises(AtctParseError) as excinfo:
        parse_atct_rows(row)
    assert "missing units" in str(excinfo.value)


def test_exact_row_with_no_units_cell_still_parses():
    """The converse of the check above: the element rows must keep working.

    ATcT publishes elements in their standard state with an 'exact'
    uncertainty and no units at all; rejecting those would silently drop
    every element from the output.
    """
    row = _single_row_html("900-00-0", "<a href=\'#\'>Fakium</a>", "Fk (g)")
    row = row.replace('<span class="Uncert">&plusmn; 0.10</span>', '<span class="Uncert">exact</span>')
    row = row.replace('<span class="Units">kJ/mol</span>', '<span class="Units"></span>')
    parsed = parse_atct_rows(row)
    assert len(parsed) == 1
    assert parsed[0].uncertainty_kJ_per_mol == 0.0


def test_clear_standard_state_removes_a_stale_block(tmp_path):
    """A fluid that drops out of ATcT must lose its block, not keep it.

    Writing is matched-only, so without an explicit removal pass a fluid
    flipping matched -> absent between versions would keep serving its
    previous value under the OLD version stamp while the ledger and the
    report both recorded it as absent.
    """
    import json as _json

    path = tmp_path / "Fakium.json"
    path.write_text(
        _json.dumps(
            {
                "INFO": {
                    "NAME": "Fakium",
                    "CAS": "900-00-0",
                    "STANDARD_STATE": {"hmolar_formation": {"value": -1000.0, "version": "1.220"}},
                }
            }
        ),
        encoding="utf-8",
    )

    assert clear_standard_state(path) is True
    doc = _json.loads(path.read_text(encoding="utf-8"))
    assert "STANDARD_STATE" not in doc["INFO"]
    # Nothing else may be disturbed by the removal.
    assert doc["INFO"]["NAME"] == "Fakium"
    assert doc["INFO"]["CAS"] == "900-00-0"


def test_clear_standard_state_leaves_another_sources_data_alone(tmp_path):
    """Removal must take this tool's quantity, not the whole shared block.

    dev/atct/README.md's Scope section plans an entropy tier from a different
    source (ATcT publishes no standard entropies).  Once that lands,
    STANDARD_STATE is shared, and deleting the object wholesale would destroy
    the other source's data for every fluid ATcT does not cover.
    """
    import json as _json

    path = tmp_path / "Fakium.json"
    path.write_text(
        _json.dumps(
            {
                "INFO": {
                    "NAME": "Fakium",
                    "STANDARD_STATE": {
                        "T": 298.15,
                        "phase": "ideal_gas",
                        "hmolar_formation": {"value": -1000.0, "source": "ATcT"},
                        "smolar_standard": {"value": 186.25, "source": "SomeOtherDatabase"},
                    },
                }
            }
        ),
        encoding="utf-8",
    )

    assert clear_standard_state(path) is True
    block = _json.loads(path.read_text(encoding="utf-8"))["INFO"]["STANDARD_STATE"]
    assert "hmolar_formation" not in block
    assert block["smolar_standard"] == {"value": 186.25, "source": "SomeOtherDatabase"}
    assert block["T"] == 298.15  # scaffolding kept: another quantity still needs it


def test_clear_standard_state_will_not_delete_another_sources_formation_value(tmp_path):
    """If hmolar_formation did not come from ATcT, this tool does not own it."""
    import json as _json

    path = tmp_path / "Fakium.json"
    payload = _json.dumps(
        {"INFO": {"NAME": "Fakium",
                  "STANDARD_STATE": {"hmolar_formation": {"value": -1.0, "source": "CODATA"}}}}
    )
    path.write_text(payload, encoding="utf-8")

    assert clear_standard_state(path) is False
    assert path.read_text(encoding="utf-8") == payload


def test_clear_standard_state_is_a_noop_without_a_block(tmp_path):
    import json as _json

    path = tmp_path / "Fakium.json"
    payload = _json.dumps({"INFO": {"NAME": "Fakium", "CAS": "900-00-0"}})
    path.write_text(payload, encoding="utf-8")

    assert clear_standard_state(path) is False
    # Untouched, byte for byte -- an absent fluid with no block must not be
    # rewritten just to reserialize it.
    assert path.read_text(encoding="utf-8") == payload


def test_write_clears_the_block_of_a_fluid_that_left_the_source(tmp_path, monkeypatch, capsys):
    """End-to-end: --write must remove, not preserve, a now-stale block.

    Fakium carries a block from a previous run but has no row on the page,
    so bind() records it as absent.  Before the removal pass, --write left
    the block in place and the fluid kept reporting a value that the ledger
    and atct_report.csv both denied.
    """
    import json as _json

    module = sys.modules["fetch_atct_formation"]
    monkeypatch.setattr(module, "__file__", str(tmp_path / "fetch_atct_formation.py"))

    fluids_dir = tmp_path / "fluids"
    fluids_dir.mkdir()
    (fluids_dir / "Methane.json").write_text(
        _json.dumps({"INFO": {"NAME": "Methane", "CAS": "74-82-8"}}), encoding="utf-8"
    )
    stale = fluids_dir / "Fakium.json"
    stale.write_text(
        _json.dumps(
            {
                "INFO": {
                    "NAME": "Fakium",
                    "CAS": "900-00-0",
                    "STANDARD_STATE": {"hmolar_formation": {"value": -1000.0, "version": "1.219"}},
                }
            }
        ),
        encoding="utf-8",
    )

    exit_code = main(
        [
            "--version", "1.220",
            "--cache", str(FIXTURE),
            "--fluids-dir", str(fluids_dir),
            "--update-ledger",
            "--write",
            "--allow-removals",
        ]
    )
    captured = capsys.readouterr()

    assert exit_code == 0
    # A deletion is reported on stderr, not buried in progress output.
    assert "REMOVED the ATcT formation value from Fakium" in captured.err

    assert "STANDARD_STATE" not in _json.loads(stale.read_text(encoding="utf-8"))["INFO"]
    # ...and the matched fluid still gained its block on the same pass.
    assert "STANDARD_STATE" in _json.loads((fluids_dir / "Methane.json").read_text(encoding="utf-8"))["INFO"]


def test_write_refuses_to_delete_blocks_without_allow_removals(tmp_path, monkeypatch, capsys):
    """Deleting a committed value must require an explicit opt-in.

    Every path that lands a fluid in `absent` is quiet: parse_atct_rows()
    skips rows missing a formula button or a DHf298 span, the only floor is
    one surviving row, and fetch_index() caches the page before any sanity
    check.  So a truncated download reports most fluids absent, and paired
    with --update-ledger (which rewrites the ledger AND the report on the same
    pass) the deletion would be self-consistent and invisible.  Verified
    against a degraded page: 76 blocks went to 1, exit 0, nothing on stderr.

    The gate must fire BEFORE anything is written, so a refused run leaves the
    tree untouched rather than half-updated.
    """
    import json as _json

    module = sys.modules["fetch_atct_formation"]
    monkeypatch.setattr(module, "__file__", str(tmp_path / "fetch_atct_formation.py"))

    fluids_dir = tmp_path / "fluids"
    fluids_dir.mkdir()
    methane = fluids_dir / "Methane.json"
    methane.write_text(_json.dumps({"INFO": {"NAME": "Methane", "CAS": "74-82-8"}}), encoding="utf-8")
    stale = fluids_dir / "Fakium.json"
    stale_text = _json.dumps(
        {"INFO": {"NAME": "Fakium", "CAS": "900-00-0",
                  "STANDARD_STATE": {"hmolar_formation": {"value": -1000.0, "version": "1.219"}}}}
    )
    stale.write_text(stale_text, encoding="utf-8")

    exit_code = main(
        [
            "--version", "1.220",
            "--cache", str(FIXTURE),
            "--fluids-dir", str(fluids_dir),
            "--update-ledger",
            "--write",
        ]
    )
    captured = capsys.readouterr()

    assert exit_code == 1
    assert "refusing to delete committed values" in captured.err
    assert "Fakium" in captured.err
    assert "--allow-removals" in captured.err

    # Nothing written: the stale file is byte-identical and the matched fluid
    # did NOT get its block, so a refused run is not a partial run.
    assert stale.read_text(encoding="utf-8") == stale_text
    assert "STANDARD_STATE" not in _json.loads(methane.read_text(encoding="utf-8"))["INFO"]


def test_refused_run_leaves_the_ledger_untouched(tmp_path, monkeypatch, capsys):
    """"nothing was written" must include expected_coverage.json.

    The ledger branch used to run before the gate, so a refused
    --update-ledger --write still overwrote the ledger with the degraded
    coverage.  That silences the loudest tripwire: a later plain --write sails
    past compare_ledger because the ledger now agrees with the bad page,
    leaving only the gate the operator is already being tempted to override.
    """
    import json as _json

    module = sys.modules["fetch_atct_formation"]
    monkeypatch.setattr(module, "__file__", str(tmp_path / "fetch_atct_formation.py"))

    fluids_dir = tmp_path / "fluids"
    fluids_dir.mkdir()
    (fluids_dir / "Methane.json").write_text(
        _json.dumps({"INFO": {"NAME": "Methane", "CAS": "74-82-8"}}), encoding="utf-8"
    )
    (fluids_dir / "Fakium.json").write_text(
        _json.dumps({"INFO": {"NAME": "Fakium", "CAS": "900-00-0",
                              "STANDARD_STATE": {"hmolar_formation": {"value": -1.0, "source": "ATcT"}}}}),
        encoding="utf-8",
    )
    ledger_path = tmp_path / "expected_coverage.json"
    ledger_text = _json.dumps({"Fakium": {"state": "matched", "atct_id": "900-00-0*0"}})
    ledger_path.write_text(ledger_text, encoding="utf-8")

    exit_code = main(["--version", "1.220", "--cache", str(FIXTURE),
                      "--fluids-dir", str(fluids_dir), "--update-ledger", "--write"])
    assert exit_code == 1
    assert ledger_path.read_text(encoding="utf-8") == ledger_text


def test_write_refuses_to_overwrite_another_sources_formation_value(tmp_path, monkeypatch, capsys):
    """A MATCHED fluid's foreign value must be as protected as an absent one's.

    clear_standard_state refused to delete another source's hmolar_formation,
    but write_standard_state replaced the whole STANDARD_STATE object -- so the
    same value was safe on the ~60 absent fluids and destroyed on all 76
    matched ones, with no flag and no message.
    """
    import json as _json

    module = sys.modules["fetch_atct_formation"]
    monkeypatch.setattr(module, "__file__", str(tmp_path / "fetch_atct_formation.py"))

    fluids_dir = tmp_path / "fluids"
    fluids_dir.mkdir()
    methane = fluids_dir / "Methane.json"
    payload = _json.dumps({"INFO": {"NAME": "Methane", "CAS": "74-82-8",
                                    "STANDARD_STATE": {"hmolar_formation": {"value": -1.0, "source": "CODATA"}}}})
    methane.write_text(payload, encoding="utf-8")

    exit_code = main(["--version", "1.220", "--cache", str(FIXTURE),
                      "--fluids-dir", str(fluids_dir), "--update-ledger", "--write"])
    captured = capsys.readouterr()

    assert exit_code == 1
    assert "another source" in captured.err
    assert "CODATA" in captured.err
    assert methane.read_text(encoding="utf-8") == payload


def test_write_merges_into_a_block_shared_with_another_quantity(tmp_path, rows):
    """Writing must add our quantity, not replace the shared object.

    The planned entropy tier will live in the same STANDARD_STATE block; a
    replacing writer destroys it for every matched fluid on the first run.
    """
    import json as _json

    path = tmp_path / "Methane.json"
    path.write_text(
        _json.dumps({"INFO": {"NAME": "Methane", "CAS": "74-82-8",
                              "STANDARD_STATE": {"smolar_standard": {"value": 186.25,
                                                                     "source": "SomeOtherDatabase"}}}}),
        encoding="utf-8",
    )
    methane = [r for r in rows if r.cas == "74-82-8"][0]
    write_standard_state(path, methane, "1.220")

    block = _json.loads(path.read_text(encoding="utf-8"))["INFO"]["STANDARD_STATE"]
    assert block["hmolar_formation"]["value"] == -74513.0
    assert block["smolar_standard"] == {"value": 186.25, "source": "SomeOtherDatabase"}


def test_writers_preserve_a_trailing_newline(tmp_path, rows):
    """6 of 137 fluid files end with a newline; inserting a block must not strip it."""
    import json as _json

    path = tmp_path / "Methane.json"
    before = _json.dumps({"INFO": {"NAME": "Methane", "CAS": "74-82-8"}}, indent=2) + "\n"
    path.write_text(before, encoding="utf-8")

    methane = [r for r in rows if r.cas == "74-82-8"][0]
    write_standard_state(path, methane, "1.220")
    assert path.read_text(encoding="utf-8").endswith("\n")

    assert clear_standard_state(path) is True
    assert path.read_text(encoding="utf-8") == before
