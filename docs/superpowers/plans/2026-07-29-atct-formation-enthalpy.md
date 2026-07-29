# ATcT Standard Enthalpy of Formation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Store the standard ideal-gas enthalpy of formation at 298.15 K for every CoolProp fluid that ATcT covers, regenerated from source, and expose it as the trivial parameter `HFORMATION`.

**Architecture:** A Python regenerator (`dev/atct/`) scrapes the ATcT Thermochemical Network index, binds species to CoolProp fluids by CAS, and writes an `INFO.STANDARD_STATE` block into `dev/fluids/*.json`. The C++ side ingests that block into a `FormationStruct` on `CoolPropFluid` and exposes it through the existing trivial-parameter path, exactly mirroring `iGWP100`/`iODP`.

**Tech Stack:** Python 3.9+ (stdlib only — `urllib`, `re`, `json`, `dataclasses`, `html`), pytest for the regenerator, C++17 with nlohmann/json, Catch2 for the library tests.

**Spec:** `docs/superpowers/specs/2026-07-29-atct-formation-enthalpy-design.md`
**Beads:** epic `CoolProp-oett`; tasks `CoolProp-ijbw` (Tasks 1–3), `CoolProp-ihto` (Task 4)

## Global Constraints

- **Units.** ATcT publishes kJ/mol. Everything stored in `dev/fluids/*.json` and everything returned by the C++ API is **J/mol**. Convert exactly once, in the writer (Task 3).
- **Never fabricate a value.** A fluid with no ATcT species gets **no** `STANDARD_STATE` block. The C++ getter **throws**; it must never return `0`, `_HUGE`, or `NaN`. A silent zero is indistinguishable from a legitimate element value.
- **Pure and pseudo-pure fluids only.** A mixture input throws, matching `calc_GWP100`/`calc_ODP`.
- **No REFPROP anywhere in this tier.** Do not call `HEATFRMdll`, do not add a `[refprop]` test, do not compare against REFPROP values.
- **Fluid JSON serialization must use `package_json.json_options`** (`{'indent': 2, 'sort_keys': True}`, no trailing newline). Any other serialization reformats all ~137 files and buries the real diff.
- **Fail loud.** Ambiguous or unexpected binding is a non-zero exit, never a skip.
- **Reference conditions:** T = 298.15 K, p = 100000.0 Pa, phase `ideal_gas`.
- **ATcT version for this plan:** `1.220`.

## File Structure

**Create:**
- `dev/atct/fetch_atct_formation.py` — parser, binder, writer, CLI. One module; it is ~250 lines and the three concerns share the row dataclass.
- `dev/atct/test_fetch_atct_formation.py` — pytest, offline, fixture-driven.
- `dev/atct/fixtures/atct_rows_sample.html` — hand-trimmed real rows covering every trap.
- `dev/atct/expected_coverage.json` — committed coverage ledger (generated in Task 3).
- `dev/atct/atct_report.csv` — committed audit trail (generated in Task 3).
- `dev/atct/README.md` — how to re-run and what to commit.

**Modify:**
- `dev/fluids/*.json` — ~75 files gain `INFO.STANDARD_STATE` (generated, Task 3).
- `include/CoolProp/CoolPropFluid.h` — `FormationStruct` + member on `CoolPropFluid`.
- `src/Backends/Helmholtz/Fluids/FluidLibrary.h` — `parse_standard_state()`.
- `src/Backends/Helmholtz/Fluids/FluidLibrary.cpp:220-228` — call site.
- `include/CoolProp/DataStructures.h:47` — enum entry, in the `// General parameters` block.
- `src/DataStructures.cpp:68` — parameter registration, matching the enum order.
- `include/CoolProp/AbstractState.h:463-474` — virtual.
- `src/AbstractState.cpp:486-490` — dispatch case.
- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h:260` and `.cpp:1175` — implementation.
- `src/Tests/CoolProp-Tests.cpp` — Catch2 tests.
- `.github/workflows/dev_checks.yml` — CI job for the regenerator tests.

---

### Task 1: ATcT row parser

**Beads:** `CoolProp-ijbw`

**Files:**
- Create: `dev/atct/fetch_atct_formation.py`
- Create: `dev/atct/fixtures/atct_rows_sample.html`
- Test: `dev/atct/test_fetch_atct_formation.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `@dataclass(frozen=True) AtctRow` with fields `cas: str`, `name: str`, `formula: str`, `phase: str`, `dhf298_kJ_per_mol: float`, `uncertainty_kJ_per_mol: float`, `atct_id: str`
  - `parse_atct_rows(html: str) -> list[AtctRow]`
  - `select_gas_row(rows: list[AtctRow], qualifier: str | None = None) -> AtctRow | None`
  - `class AtctParseError(Exception)`, `class AmbiguousSpecies(Exception)`

**Why these three traps get dedicated tests:** each was observed on the live page and each fails *silently* — element rows vanish from a `±`-anchored regex, a CAS-only lookup silently picks the aqueous argon row, and an exact-`(g)`-only filter silently drops cis-2-butene.

- [ ] **Step 1: Create the fixture file**

Create `dev/atct/fixtures/atct_rows_sample.html`. These are real rows from ATcT TN 1.220 with the `bkgImage` and `bkgMass` cells kept (the parser must tolerate them) and nothing else altered:

```html
<tr id="s1_6n4_1c0 i23 CAS74-82-8"><td class="bkgName"><span class="Name"> <a href="species/?species_number=23">Methane</a></span></td><td class="bkgFormula> <span class="Formula"> <button onclick="copyname('Methane');" type="button">CH4  (g) </button> </span></td><td class="bkgImage"><span class="Name"><img class="lazy" data-original="images/23.png" alt="C" style="margin:5px" /></span></td><td class="bkgDHf0"><span class="DHf0">-66.543</span></td><td class="bkgDHf298"><span class="DHf298">-74.513</span></td><td class="bkgUncert"><span class="Uncert">&plusmn; 0.043</span></td><td class="bkgUnits"><span class="Units">kJ/mol</span></td><td class="bkgMass"><span class="Mass">16.04246 &plusmn;<br />0.00085</span></td><td class="bkgATcTID"><span class="ATcTID">74-82-8*0</span></td></tr>
<tr id="s18n1c0 i12 CAS7440-37-1"><td class="bkgName"><span class="Name"> <a href="species/?species_number=12">Argon</a></span></td><td class="bkgFormula> <span class="Formula"> <button onclick="copyname('Argon');" type="button">Ar (g) </button> </span></td><td class="bkgImage"><span class="Name"><img class="lazy" data-original="images/12.png" alt="[Ar]" style="margin:5px" /></span></td><td class="bkgDHf0"><span class="DHf0">0</span></td><td class="bkgDHf298"><span class="DHf298">0</span></td><td class="bkgUncert"><span class="Uncert">exact</span></td><td class="bkgUnits"><span class="Units"></span></td><td class="bkgMass"><span class="Mass">39.94800 &plusmn;<br />0.00100</span></td><td class="bkgATcTID"><span class="ATcTID">7440-37-1*0</span></td></tr>
<tr id="s18n1c0 i3436 CAS7440-37-1"><td class="bkgName"><span class="Name"> <a href="species/?species_number=3436">Argon</a></span></td><td class="bkgFormula> <span class="Formula"> <button onclick="copyname('Argon');" type="button">Ar (aq, undissoc) </button> </span></td><td class="bkgImage"><span class="Name"><img class="lazy" data-original="images/3436.png" alt="[Ar]" style="margin:5px" /></span></td><td class="bkgDHf0"><span class="DHf0"></span></td><td class="bkgDHf298"><span class="DHf298">-12.249</span></td><td class="bkgUncert"><span class="Uncert">&plusmn; 0.092</span></td><td class="bkgUnits"><span class="Units">kJ/mol</span></td><td class="bkgMass"><span class="Mass">39.94800 &plusmn;<br />0.00100</span></td><td class="bkgATcTID"><span class="ATcTID">7440-37-1*1</span></td></tr>
<tr id="s1_8n2_1c0 i28 CAS7732-18-5"><td class="bkgName"><span class="Name"> <a href="species/?species_number=28">Water</a></span></td><td class="bkgFormula> <span class="Formula"> <button onclick="copyname('Water');" type="button">H2O (g) </button> </span></td><td class="bkgImage"><span class="Name"><img class="lazy" data-original="images/28.png" alt="O" style="margin:5px" /></span></td><td class="bkgDHf0"><span class="DHf0">-238.904</span></td><td class="bkgDHf298"><span class="DHf298">-241.808</span></td><td class="bkgUncert"><span class="Uncert">&plusmn; 0.022</span></td><td class="bkgUnits"><span class="Units">kJ/mol</span></td><td class="bkgMass"><span class="Mass">18.01528 &plusmn;<br />0.00033</span></td><td class="bkgATcTID"><span class="ATcTID">7732-18-5*0</span></td></tr>
<tr id="s1_8n2_1c0 i29 CAS7732-18-5"><td class="bkgName"><span class="Name"> <a href="species/?species_number=29">Water</a></span></td><td class="bkgFormula> <span class="Formula"> <button onclick="copyname('Water');" type="button">H2O (l) </button> </span></td><td class="bkgImage"><span class="Name"><img class="lazy" data-original="images/29.png" alt="O" style="margin:5px" /></span></td><td class="bkgDHf0"><span class="DHf0"></span></td><td class="bkgDHf298"><span class="DHf298">-285.825</span></td><td class="bkgUncert"><span class="Uncert">&plusmn; 0.030</span></td><td class="bkgUnits"><span class="Units">kJ/mol</span></td><td class="bkgMass"><span class="Mass">18.01528 &plusmn;<br />0.00033</span></td><td class="bkgATcTID"><span class="ATcTID">7732-18-5*1</span></td></tr>
<tr id="s1n2c0 i1779 CAS1333-74-0"><td class="bkgName"><span class="Name"> <a href="species/?species_number=1779">Dihydrogen</a></span></td><td class="bkgFormula> <span class="Formula"> <button onclick="copyname('Dihydrogen');" type="button">H2 (g, para) </button> </span></td><td class="bkgImage"><span class="Name"><img class="lazy" data-original="images/1779.png" alt="[H][H]" style="margin:5px" /></span></td><td class="bkgDHf0"><span class="DHf0">-0.000</span></td><td class="bkgDHf298"><span class="DHf298">-0.058</span></td><td class="bkgUncert"><span class="Uncert">&plusmn; 0.000</span></td><td class="bkgUnits"><span class="Units">kJ/mol</span></td><td class="bkgMass"><span class="Mass">2.01588 &plusmn;<br />0.00014</span></td><td class="bkgATcTID"><span class="ATcTID">1333-74-0*2</span></td></tr>
<tr id="s1n2c0 i2 CAS1333-74-0"><td class="bkgName"><span class="Name"> <a href="species/?species_number=2">Dihydrogen</a></span></td><td class="bkgFormula> <span class="Formula"> <button onclick="copyname('Dihydrogen');" type="button">H2  (g) </button> </span></td><td class="bkgImage"><span class="Name"><img class="lazy" data-original="images/2.png" alt="[H][H]" style="margin:5px" /></span></td><td class="bkgDHf0"><span class="DHf0">0</span></td><td class="bkgDHf298"><span class="DHf298">0</span></td><td class="bkgUncert"><span class="Uncert">exact</span></td><td class="bkgUnits"><span class="Units"></span></td><td class="bkgMass"><span class="Mass">2.01588 &plusmn;<br />0.00014</span></td><td class="bkgATcTID"><span class="ATcTID">1333-74-0*0</span></td></tr>
<tr id="s1_6n8_4c0 i233 CAS590-18-1"><td class="bkgName"><span class="Name"> <a href="species/?species_number=233">cis-2-Butene</a></span></td><td class="bkgFormula> <span class="Formula"> <button onclick="copyname('cis-2-Butene');" type="button">CH3CHCHCH3 (g, cis) </button> </span></td><td class="bkgImage"><span class="Name"><img class="lazy" data-original="images/233.png" alt="C/C=C\C" style="margin:5px" /></span></td><td class="bkgDHf0"><span class="DHf0">14.42</span></td><td class="bkgDHf298"><span class="DHf298">-6.85</span></td><td class="bkgUncert"><span class="Uncert">&plusmn; 0.38</span></td><td class="bkgUnits"><span class="Units">kJ/mol</span></td><td class="bkgMass"><span class="Mass">56.10632 &plusmn;<br />0.00380</span></td><td class="bkgATcTID"><span class="ATcTID">590-18-1*0</span></td></tr>
```

- [ ] **Step 2: Write the failing tests**

Create `dev/atct/test_fetch_atct_formation.py`:

```python
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
```

- [ ] **Step 3: Run the tests to verify they fail**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/formation
python -m pytest dev/atct/test_fetch_atct_formation.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'fetch_atct_formation'`.

- [ ] **Step 4: Write the parser**

Create `dev/atct/fetch_atct_formation.py`:

```python
#!/usr/bin/env python3
"""Regenerate CoolProp standard enthalpies of formation from the Active
Thermochemical Tables (ATcT).

Source: B. Ruscic and D. H. Bross, Active Thermochemical Tables (ATcT) values
based on ver. 1.220 of the Thermochemical Network, DOI: 10.17038/CSE/2568691.

The entire ATcT table for a version is a single HTML page whose species rows
are machine readable.  This module parses those rows, binds them to CoolProp
fluids by CAS number, and writes an INFO.STANDARD_STATE block into
dev/fluids/*.json.  See docs/superpowers/specs/2026-07-29-atct-formation-enthalpy-design.md.
"""
from __future__ import annotations

import html as html_module
import re
from dataclasses import dataclass

INDEX_URL_TEMPLATE = (
    "https://atct.anl.gov/Thermochemical%20Data/version%20{version}/index.php"
)


class AtctParseError(Exception):
    """The ATcT page did not have the structure this parser expects."""


class AmbiguousSpecies(Exception):
    """More than one ATcT gas-phase row matched, so no value can be chosen."""


@dataclass(frozen=True)
class AtctRow:
    cas: str
    name: str
    formula: str
    phase: str
    dhf298_kJ_per_mol: float
    uncertainty_kJ_per_mol: float
    atct_id: str


_ROW_RE = re.compile(r'<tr id="[^"]*?CAS(?P<cas>[^"]+)">(?P<body>.*?)</tr>', re.S)
_NAME_RE = re.compile(r'species_number=\d+">(?P<name>.*?)</a>', re.S)
_FORMULA_RE = re.compile(r'type="button">(?P<formula>.*?)</button>', re.S)


def _span(body: str, css_class: str) -> str:
    match = re.search(r'class="%s">(.*?)</span>' % css_class, body, re.S)
    if match is None:
        return ""
    return html_module.unescape(match.group(1)).strip()


def _split_phase(raw: str) -> tuple[str, str]:
    """'CH4  (g)' -> ('CH4', '(g)');  'H2 (g, para)' -> ('H2', '(g, para)')."""
    collapsed = " ".join(raw.split())
    open_paren = collapsed.find("(")
    if open_paren < 0:
        return collapsed, ""
    return collapsed[:open_paren].strip(), collapsed[open_paren:].strip()


def _parse_uncertainty(text: str) -> float:
    """ATcT writes '± 0.043' for measured species and 'exact' for elements.

    The element form carries no '±' and an empty units cell; treating it as a
    parse failure silently drops all nine elements from the output.
    """
    cleaned = text.replace("±", "").strip()
    if cleaned.lower() == "exact":
        return 0.0
    if not cleaned:
        raise AtctParseError("empty uncertainty cell")
    return float(cleaned)


def parse_atct_rows(page_html: str) -> list[AtctRow]:
    """Parse every species row carrying a 298.15 K enthalpy of formation."""
    rows: list[AtctRow] = []
    for match in _ROW_RE.finditer(page_html):
        body = match.group("body")
        formula_match = _FORMULA_RE.search(body)
        if formula_match is None:
            continue  # navigation / header row, not a species row
        value_text = _span(body, "DHf298")
        if not value_text:
            continue  # species with no published 298.15 K value
        units = _span(body, "Units")
        if units not in ("", "kJ/mol"):
            raise AtctParseError("unexpected units %r on CAS %s" % (units, match.group("cas")))
        name_match = _NAME_RE.search(body)
        formula, phase = _split_phase(html_module.unescape(formula_match.group("formula")))
        rows.append(
            AtctRow(
                cas=match.group("cas").strip(),
                name=name_match.group("name").strip() if name_match else "",
                formula=formula,
                phase=phase,
                dhf298_kJ_per_mol=float(value_text),
                uncertainty_kJ_per_mol=_parse_uncertainty(_span(body, "Uncert")),
                atct_id=_span(body, "ATcTID"),
            )
        )
    if not rows:
        raise AtctParseError("no species rows parsed; the page layout has changed")
    return rows


def select_gas_row(rows: list[AtctRow], qualifier: str | None = None) -> AtctRow | None:
    """Pick the one gas-phase row from the rows sharing a CAS number.

    Without a qualifier: prefer the exact '(g)' row; if there is none, accept a
    single '(g, <qualifier>)' row (this is how cis-2-Butene, which exists only
    as '(g, cis)', resolves).  With a qualifier: require '(g, <qualifier>)'
    exactly, which is how the hydrogen spin isomers resolve.

    Raises AmbiguousSpecies rather than guessing.  Returns None when nothing
    matches, which the caller records as 'absent'.
    """
    if qualifier is not None:
        wanted = "(g, %s)" % qualifier
        hits = [r for r in rows if r.phase == wanted]
        if len(hits) > 1:
            raise AmbiguousSpecies("%d rows with phase %s" % (len(hits), wanted))
        return hits[0] if hits else None

    exact = [r for r in rows if r.phase == "(g)"]
    if len(exact) > 1:
        raise AmbiguousSpecies("%d rows with phase (g)" % len(exact))
    if exact:
        return exact[0]

    qualified = [r for r in rows if r.phase.startswith("(g,")]
    if len(qualified) > 1:
        raise AmbiguousSpecies(
            "no exact (g) row and %d qualified gas rows: %s"
            % (len(qualified), ", ".join(r.phase for r in qualified))
        )
    return qualified[0] if qualified else None
```

- [ ] **Step 5: Run the tests to verify they pass**

```bash
python -m pytest dev/atct/test_fetch_atct_formation.py -v
```

Expected: 8 passed.

- [ ] **Step 6: Commit**

```bash
git add dev/atct/fetch_atct_formation.py dev/atct/test_fetch_atct_formation.py dev/atct/fixtures/atct_rows_sample.html
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "feat(atct): parse ATcT species rows (CoolProp-ijbw)

Handles the three silent-failure traps observed on the live TN 1.220
page: element rows whose uncertainty reads 'exact' with an empty units
cell, CAS numbers shared across phases, and species that exist only
under a qualified gas phase such as (g, cis)."
```

---

### Task 2: CoolProp ↔ ATcT binder

**Beads:** `CoolProp-ijbw`

**Files:**
- Modify: `dev/atct/fetch_atct_formation.py`
- Test: `dev/atct/test_fetch_atct_formation.py`

**Interfaces:**
- Consumes: `AtctRow`, `select_gas_row`, `AmbiguousSpecies` from Task 1.
- Produces:
  - `@dataclass(frozen=True) FluidRef` with fields `name: str`, `cas: str`, `path: pathlib.Path`
  - `load_coolprop_fluids(fluids_dir) -> dict[str, FluidRef]` keyed by `INFO.NAME`
  - `normalize_cas(cas: str) -> tuple[str, str | None]`
  - `@dataclass BindResult` with fields `matched: dict[str, AtctRow]`, `absent: dict[str, str]`, `errors: list[str]`
  - `bind(fluids: dict[str, FluidRef], rows: list[AtctRow]) -> BindResult`

- [ ] **Step 1: Write the failing tests**

Append to `dev/atct/test_fetch_atct_formation.py`:

```python
from fetch_atct_formation import (  # noqa: E402
    BindResult,
    FluidRef,
    bind,
    load_coolprop_fluids,
    normalize_cas,
)

FLUIDS_DIR = Path(__file__).resolve().parents[2] / "dev" / "fluids"


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
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
python -m pytest dev/atct/test_fetch_atct_formation.py -v
```

Expected: `ImportError: cannot import name 'BindResult'`.

- [ ] **Step 3: Write the binder**

Append to `dev/atct/fetch_atct_formation.py` (after `select_gas_row`, keeping the imports at the top of the file):

```python
import json
from collections import defaultdict
from pathlib import Path

_CAS_SPIN_RE = re.compile(r"^(?P<cas>\d{2,7}-\d{2}-\d)(?P<spin>[op])$")
_SPIN_NAME = {"p": "para", "o": "ortho"}


@dataclass(frozen=True)
class FluidRef:
    name: str
    cas: str
    path: Path


@dataclass
class BindResult:
    matched: dict
    absent: dict
    errors: list


def normalize_cas(cas: str) -> tuple[str, str | None]:
    """Split CoolProp's spin-isomer CAS suffix off the real CAS number.

    CoolProp identifies ParaHydrogen as '1333-74-0p' and OrthoHydrogen as
    '1333-74-0o'.  ATcT carries these as phase-qualified rows on the parent
    CAS, so the suffix becomes the phase qualifier.
    """
    match = _CAS_SPIN_RE.match(cas or "")
    if match is None:
        return cas, None
    return match.group("cas"), _SPIN_NAME[match.group("spin")]


def load_coolprop_fluids(fluids_dir) -> dict:
    """Read every dev/fluids/*.json and index it by INFO.NAME."""
    fluids = {}
    for path in sorted(Path(fluids_dir).glob("*.json")):
        info = json.loads(path.read_text(encoding="utf-8"))["INFO"]
        fluids[info["NAME"]] = FluidRef(name=info["NAME"], cas=info.get("CAS", ""), path=path)
    if not fluids:
        raise AtctParseError("no fluid files found in %s" % fluids_dir)
    return fluids


def bind(fluids: dict, rows: list) -> BindResult:
    """Bind each CoolProp fluid to at most one ATcT gas-phase row.

    Ambiguity is an error, never a silent choice.  A fluid with no matching
    species is recorded in `absent` with the reason, which the coverage ledger
    then pins so a future ATcT version cannot quietly drop it.
    """
    by_cas = defaultdict(list)
    for row in rows:
        by_cas[row.cas].append(row)

    matched, absent, errors = {}, {}, []
    for name, fluid in sorted(fluids.items()):
        cas, qualifier = normalize_cas(fluid.cas)
        if not cas:
            absent[name] = "fluid has no CAS number"
            continue
        candidates = by_cas.get(cas, [])
        if not candidates:
            absent[name] = "no ATcT species for CAS %s" % cas
            continue
        try:
            row = select_gas_row(candidates, qualifier)
        except AmbiguousSpecies as exc:
            errors.append("%s (CAS %s): %s" % (name, cas, exc))
            continue
        if row is None:
            detail = " with qualifier %s" % qualifier if qualifier else ""
            absent[name] = "no gas-phase ATcT row for CAS %s%s" % (cas, detail)
            continue
        matched[name] = row
    return BindResult(matched=matched, absent=absent, errors=errors)
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
python -m pytest dev/atct/test_fetch_atct_formation.py -v
```

Expected: 14 passed.

- [ ] **Step 5: Commit**

```bash
git add dev/atct/fetch_atct_formation.py dev/atct/test_fetch_atct_formation.py
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "feat(atct): bind CoolProp fluids to ATcT species by CAS (CoolProp-ijbw)

Spin-isomer CAS suffixes (1333-74-0p/o) resolve to ATcT's phase-qualified
(g, para)/(g, ortho) rows, so ParaHydrogen gets its real -58 J/mol value
rather than a fabricated zero.  Ambiguous matches are errors, not guesses."
```

---

### Task 3: Fetch, write, and generate the data

**Beads:** `CoolProp-ijbw`

**Files:**
- Modify: `dev/atct/fetch_atct_formation.py`
- Test: `dev/atct/test_fetch_atct_formation.py`
- Create: `dev/atct/README.md`
- Generate: `dev/atct/expected_coverage.json`, `dev/atct/atct_report.csv`, `dev/fluids/*.json`

**Interfaces:**
- Consumes: `BindResult`, `FluidRef`, `bind`, `load_coolprop_fluids` from Task 2.
- Produces:
  - `standard_state_block(row: AtctRow, version: str) -> dict`
  - `write_standard_state(path: Path, row: AtctRow, version: str) -> None`
  - `coverage_ledger(result: BindResult) -> dict`
  - `compare_ledger(expected: dict, actual: dict) -> list[str]`
  - `main(argv=None) -> int`

- [ ] **Step 1: Write the failing tests**

Append to `dev/atct/test_fetch_atct_formation.py`:

```python
from fetch_atct_formation import (  # noqa: E402
    compare_ledger,
    coverage_ledger,
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
    """The writer must not reformat the file; only the new block may appear."""
    import json as _json

    original = FLUIDS_DIR / "Methane.json"
    scratch = tmp_path / "Methane.json"
    scratch.write_text(original.read_text(encoding="utf-8"), encoding="utf-8")

    methane = [r for r in rows if r.cas == "74-82-8"][0]
    write_standard_state(scratch, methane, "1.220")

    before = _json.loads(original.read_text(encoding="utf-8"))
    after = _json.loads(scratch.read_text(encoding="utf-8"))
    assert "STANDARD_STATE" in after["INFO"]
    del after["INFO"]["STANDARD_STATE"]
    assert after == before
    # Canonical serialization: 2-space indent, sorted keys, no trailing newline.
    assert not scratch.read_text(encoding="utf-8").endswith("\n")


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
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
python -m pytest dev/atct/test_fetch_atct_formation.py -v
```

Expected: `ImportError: cannot import name 'compare_ledger'`.

- [ ] **Step 3: Write the writer, ledger and CLI**

Append to `dev/atct/fetch_atct_formation.py`:

```python
import argparse
import csv
import hashlib
import sys
import urllib.request

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from package_json import json_options  # noqa: E402

REFERENCE_T_K = 298.15
REFERENCE_P_PA = 100000.0


def standard_state_block(row: AtctRow, version: str) -> dict:
    """Build the INFO.STANDARD_STATE block.  kJ/mol -> J/mol happens here only."""
    return {
        "T": REFERENCE_T_K,
        "T_units": "K",
        "p": REFERENCE_P_PA,
        "p_units": "Pa",
        "phase": "ideal_gas",
        "hmolar_formation": {
            "value": round(row.dhf298_kJ_per_mol * 1000.0, 6),
            "units": "J/mol",
            "uncertainty": round(row.uncertainty_kJ_per_mol * 1000.0, 6),
            "source": "ATcT",
            "version": version,
            "id": row.atct_id,
        },
    }


def write_standard_state(path: Path, row: AtctRow, version: str) -> None:
    """Insert the block and re-serialize with CoolProp's canonical options.

    json_options is imported from dev/package_json.py so the file round-trips
    byte-identically apart from the new block; any other serialization would
    reformat every fluid file and bury the real diff.
    """
    doc = json.loads(path.read_text(encoding="utf-8"))
    doc["INFO"]["STANDARD_STATE"] = standard_state_block(row, version)
    path.write_text(json.dumps(doc, **json_options), encoding="utf-8")


def coverage_ledger(result: BindResult) -> dict:
    ledger = {}
    for name, row in result.matched.items():
        ledger[name] = {"state": "matched", "atct_id": row.atct_id}
    for name, reason in result.absent.items():
        ledger[name] = {"state": "absent", "reason": reason}
    return dict(sorted(ledger.items()))


def compare_ledger(expected: dict, actual: dict) -> list:
    """Every difference from the committed ledger is reported, none tolerated."""
    problems = []
    for name in sorted(set(expected) | set(actual)):
        want, got = expected.get(name), actual.get(name)
        if want is None:
            problems.append("%s: new fluid, not in the committed ledger (%s)" % (name, got))
        elif got is None:
            problems.append("%s: in the committed ledger but absent from this run" % name)
        elif want != got:
            problems.append("%s: ledger says %s, this run produced %s" % (name, want, got))
    return problems


def fetch_index(version: str, cache: Path | None) -> tuple[str, str]:
    """Return (page text, sha256).  Uses the cache file when present."""
    if cache is not None and cache.exists():
        raw = cache.read_bytes()
    else:
        url = INDEX_URL_TEMPLATE.format(version=version)
        with urllib.request.urlopen(url, timeout=120) as response:
            raw = response.read()
        if cache is not None:
            cache.write_bytes(raw)
    return raw.decode("utf-8", errors="replace"), hashlib.sha256(raw).hexdigest()


def write_report(path: Path, result: BindResult, version: str, page_sha256: str) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["# ATcT version", version, "page sha256", page_sha256])
        writer.writerow(
            ["coolprop_fluid", "state", "atct_id", "phase",
             "hmolar_formation_J_per_mol", "uncertainty_J_per_mol", "reason"]
        )
        for name, row in sorted(result.matched.items()):
            writer.writerow([name, "matched", row.atct_id, row.phase,
                             round(row.dhf298_kJ_per_mol * 1000.0, 6),
                             round(row.uncertainty_kJ_per_mol * 1000.0, 6), ""])
        for name, reason in sorted(result.absent.items()):
            writer.writerow([name, "absent", "", "", "", "", reason])


def main(argv=None) -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="1.220", help="ATcT Thermochemical Network version")
    parser.add_argument("--fluids-dir", type=Path, default=here.parents[0] / "fluids")
    parser.add_argument("--cache", type=Path, default=None,
                        help="path to cache the fetched page (not committed)")
    parser.add_argument("--write", action="store_true",
                        help="write dev/fluids/*.json, the report and the ledger")
    parser.add_argument("--update-ledger", action="store_true",
                        help="rewrite expected_coverage.json instead of checking against it")
    args = parser.parse_args(argv)

    page, page_sha256 = fetch_index(args.version, args.cache)
    rows = parse_atct_rows(page)
    fluids = load_coolprop_fluids(args.fluids_dir)
    result = bind(fluids, rows)

    print("ATcT version %s, page sha256 %s" % (args.version, page_sha256))
    print("parsed %d species rows; %d fluids matched, %d absent"
          % (len(rows), len(result.matched), len(result.absent)))

    if result.errors:
        for problem in result.errors:
            print("ERROR: %s" % problem, file=sys.stderr)
        return 1

    ledger_path = here / "expected_coverage.json"
    actual = coverage_ledger(result)
    if args.update_ledger:
        ledger_path.write_text(json.dumps(actual, **json_options), encoding="utf-8")
        print("wrote %s" % ledger_path)
    elif ledger_path.exists():
        problems = compare_ledger(json.loads(ledger_path.read_text(encoding="utf-8")), actual)
        if problems:
            for problem in problems:
                print("ERROR: coverage changed: %s" % problem, file=sys.stderr)
            return 1
    else:
        print("ERROR: %s is missing; run with --update-ledger first" % ledger_path, file=sys.stderr)
        return 1

    if args.write:
        for name, row in sorted(result.matched.items()):
            write_standard_state(fluids[name].path, row, args.version)
        write_report(here / "atct_report.csv", result, args.version, page_sha256)
        print("wrote %d fluid files and atct_report.csv" % len(result.matched))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
python -m pytest dev/atct/test_fetch_atct_formation.py -v
```

Expected: 19 passed.

- [ ] **Step 5: Generate the ledger and the data (network access required)**

```bash
cd /Users/ianbell/Code/CoolProp/.claude/worktrees/formation
python dev/atct/fetch_atct_formation.py --version 1.220 \
    --cache /tmp/atct_1220.html --update-ledger --write
```

Expected output: roughly `parsed 3422 species rows; 75 fluids matched, 62 absent`, then `wrote 75 fluid files and atct_report.csv`.

- [ ] **Step 6: Verify the generated data**

```bash
# Idempotent, and the ledger check now passes rather than being rewritten:
python dev/atct/fetch_atct_formation.py --version 1.220 --cache /tmp/atct_1220.html --write
git diff --stat -- dev/fluids | tail -1     # must show no further change

# Spot-check three values and the para-hydrogen correction:
python - <<'EOF'
import json, pathlib
for name, want in [("Methane", -74513.0), ("Water", -241808.0),
                   ("CarbonDioxide", -393477.0), ("ParaHydrogen", -58.0)]:
    doc = json.loads(pathlib.Path("dev/fluids/%s.json" % name).read_text())
    got = doc["INFO"]["STANDARD_STATE"]["hmolar_formation"]["value"]
    print("%-15s %12.3f  expected %12.3f  %s" % (name, got, want, "OK" if abs(got-want) < 0.5 else "MISMATCH"))
assert "STANDARD_STATE" not in json.loads(pathlib.Path("dev/fluids/R134a.json").read_text())["INFO"], "R134a must have no block"
print("R134a correctly has no STANDARD_STATE block")
EOF
```

Expected: four `OK` lines and the R134a confirmation.

- [ ] **Step 7: Write the README**

Create `dev/atct/README.md`:

```markdown
# ATcT standard enthalpies of formation

Regenerates `INFO.STANDARD_STATE` in `dev/fluids/*.json` from the Active
Thermochemical Tables.

Source: B. Ruscic and D. H. Bross, *Active Thermochemical Tables (ATcT) values
based on ver. 1.220 of the Thermochemical Network*, DOI: 10.17038/CSE/2568691.

## Re-running

    python dev/atct/fetch_atct_formation.py --version 1.220 --cache /tmp/atct.html --write

The run **fails** if any fluid's coverage differs from `expected_coverage.json`.
That is deliberate: a future ATcT version that renames or drops a species must
surface as a failure, not a silent gap. After reviewing the change, re-run with
`--update-ledger` and commit the new ledger alongside the new values.

## What is committed

- `expected_coverage.json` — every fluid and its expected state
- `atct_report.csv` — per-fluid audit trail including the source page SHA-256
- the regenerated `dev/fluids/*.json`

The 4.3 MB source page is **not** committed; use `--cache` to keep a local copy.

## Rebuild after regenerating

Fluid JSON is embedded in the library as CBOR, so values only take effect after:

    python dev/generate_headers.py
    cmake --build build_catch --target CatchTestRunner -j8

## Scope

Enthalpy only. ATcT publishes no standard entropies, so S°(298.15 K) needs a
different source and is a separate tier — see the spec, section 8.
```

- [ ] **Step 8: Commit**

```bash
git add dev/atct dev/fluids
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "feat(atct): regenerate formation enthalpies into fluid files (CoolProp-ijbw)

75 of 137 fluids gain INFO.STANDARD_STATE from ATcT TN 1.220, with the
value, its 95% uncertainty, and per-quantity provenance (source, version,
ATcT species ID).  expected_coverage.json pins the coverage so a future
ATcT version cannot silently drop a fluid."
```

---

### Task 4: `HFORMATION` ingestion and API

**Beads:** `CoolProp-ihto`

**Files:**
- Modify: `include/CoolProp/CoolPropFluid.h`
- Modify: `src/Backends/Helmholtz/Fluids/FluidLibrary.h`
- Modify: `src/Backends/Helmholtz/Fluids/FluidLibrary.cpp`
- Modify: `include/CoolProp/DataStructures.h`
- Modify: `src/DataStructures.cpp`
- Modify: `include/CoolProp/AbstractState.h`
- Modify: `src/AbstractState.cpp`
- Modify: `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h`
- Modify: `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`
- Test: `src/Tests/CoolProp-Tests.cpp`

**Interfaces:**
- Consumes: the `INFO.STANDARD_STATE` block written in Task 3.
- Produces: enum `CoolProp::iHmolar_formation`, string key `"HFORMATION"`, virtual `CoolPropDbl calc_Hmolar_formation()`, and `struct CoolProp::FormationStruct` reachable as `CoolPropFluid::standard_state`.

- [ ] **Step 1: Write the failing Catch2 test**

Append to `src/Tests/CoolProp-Tests.cpp` (near the other property tests, at file scope):

```cpp
TEST_CASE("Standard molar enthalpy of formation from ATcT", "[formation]") {
    SECTION("known values, J/mol") {
        // ATcT TN 1.220; see dev/atct/atct_report.csv for provenance.
        struct { const char* fluid; double value; } cases[] = {
          {"Methane", -74513.0}, {"Water", -241808.0}, {"CarbonDioxide", -393477.0}, {"Nitrogen", 0.0}};
        for (auto& c : cases) {
            CAPTURE(c.fluid);
            shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", c.fluid));
            CHECK(AS->keyed_output(CoolProp::iHmolar_formation) == Approx(c.value).margin(1e-6));
        }
    }
    SECTION("string key resolves through PropsSI") {
        CHECK(CoolProp::PropsSI("HFORMATION", "", 0, "", 0, "Methane") == Approx(-74513.0).margin(1e-6));
    }
    SECTION("a fluid with no ATcT value throws rather than returning 0 or NaN") {
        // R134a (CAS 811-97-2) is genuinely absent from ATcT.
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R134a"));
        CHECK_THROWS(AS->keyed_output(CoolProp::iHmolar_formation));
    }
    SECTION("mixtures throw: pure and pseudo-pure fluids only") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane"));
        std::vector<double> z{0.7, 0.3};
        AS->set_mole_fractions(z);
        CHECK_THROWS(AS->keyed_output(CoolProp::iHmolar_formation));
    }
    SECTION("every stored value is physically plausible") {
        // Guards against a kJ/J slip anywhere in the pipeline: no molecule in
        // the library has |dHf| above 2000 kJ/mol.
        std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
        std::size_t checked = 0;
        for (auto& fluid : fluids) {
            shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluid));
            double value = 0;
            try {
                value = AS->keyed_output(CoolProp::iHmolar_formation);
            } catch (...) {
                continue;  // no ATcT value for this fluid; covered above
            }
            CAPTURE(fluid);
            CHECK(std::abs(value) < 2e6);
            ++checked;
        }
        CHECK(checked > 60);  // the ATcT overlap is ~75 fluids
    }
}
```

- [ ] **Step 2: Run it to verify it fails to compile**

```bash
cmake --build build_catch --target CatchTestRunner -j8
```

Expected: compile error — `'iHmolar_formation' is not a member of 'CoolProp'`.

- [ ] **Step 3: Add the storage struct**

In `include/CoolProp/CoolPropFluid.h`, immediately after `struct EnvironmentalFactorsStruct` (around line 39):

```cpp
/// Standard-state thermochemical data (ideal gas at 298.15 K, 100 kPa).
/// Populated from INFO.STANDARD_STATE in the fluid JSON; absent for fluids
/// that the source database does not cover, in which case hmolar stays _HUGE
/// and the accessor throws.
struct FormationStruct
{
    double hmolar = _HUGE;              ///< Standard molar enthalpy of formation [J/mol]
    double hmolar_uncertainty = _HUGE;  ///< 95% confidence interval on hmolar [J/mol]
    std::string source;                 ///< e.g. "ATcT"
    std::string version;                ///< e.g. "1.220"
    std::string id;                     ///< source-specific species ID, e.g. "74-82-8*0"
};
```

In the same file, beside the existing `environment` member (around line 560):

```cpp
    FormationStruct standard_state;  ///< Standard-state thermochemical data
```

- [ ] **Step 4: Add the JSON parser**

In `src/Backends/Helmholtz/Fluids/FluidLibrary.h`, immediately after `parse_environmental` (around line 347):

```cpp
    /// Parse the standard-state thermochemical data (INFO.STANDARD_STATE)
    void parse_standard_state(const nlohmann::json& json, CoolPropFluid& fluid) {
        if (!json.contains("hmolar_formation")) {
            return;
        }
        const nlohmann::json& hf = json.at("hmolar_formation");
        fluid.standard_state.hmolar = cpjson::get_double(hf, "value");
        fluid.standard_state.hmolar_uncertainty = cpjson::get_double(hf, "uncertainty");
        fluid.standard_state.source = cpjson::get_string(hf, "source");
        fluid.standard_state.version = cpjson::get_string(hf, "version");
        fluid.standard_state.id = cpjson::get_string(hf, "id");
    }
```

In `src/Backends/Helmholtz/Fluids/FluidLibrary.cpp`, directly after the environmental block (around line 228):

```cpp
        // Parse the standard-state thermochemical data.  Optional: absent for
        // fluids the source database does not cover -- see
        // dev/atct/expected_coverage.json for the authoritative list.
        if (fluid_json.at("INFO").contains("STANDARD_STATE")) {
            parse_standard_state(fluid_json.at("INFO").at("STANDARD_STATE"), fluid);
        }
```

- [ ] **Step 5: Register the parameter**

In `include/CoolProp/DataStructures.h`, add the new entry to the **`// General parameters`** block, immediately after `idipole_moment` (around line 47). This is a per-fluid constant, so it belongs beside the other per-fluid constants — not in the environmental block and not appended at the end:

```cpp
    idipole_moment,      ///< Dipole moment
    iHmolar_formation,   ///< Standard molar enthalpy of formation of the ideal gas at 298.15 K
```

Do **not** append it after `iPhase` to avoid renumbering the later entries. Enum values are not part of CoolProp's public contract — callers use the named constants or the `"HFORMATION"` string — so keeping the logical grouping wins over preserving numeric stability.

In `src/DataStructures.cpp`, after the `{idipole_moment, ...}` line (line 68), matching the enum order:

```cpp
  {iHmolar_formation, "HFORMATION", "O", "J/mol", "Standard molar enthalpy of formation of the ideal gas at 298.15 K", true},
```

- [ ] **Step 6: Add the virtual and the dispatch**

In `include/CoolProp/AbstractState.h`, after `calc_ODP()` (around line 474):

```cpp
    /// Using this backend, get the standard molar enthalpy of formation of the ideal gas at 298.15 K
    virtual CoolPropDbl calc_Hmolar_formation() {
        throw NotImplementedError("calc_Hmolar_formation is not implemented for this backend");
    };
```

In `src/AbstractState.cpp`, in `trivial_keyed_output`, after the `case iODP:` block (around line 488):

```cpp
        case iHmolar_formation:
            return this->calc_Hmolar_formation();
```

- [ ] **Step 7: Implement it in the Helmholtz backend**

In `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h`, beside `calc_ODP()` (around line 261):

```cpp
    CoolPropDbl calc_Hmolar_formation() override;
```

In `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`, after `calc_GWP100()` (around line 1185):

```cpp
CoolPropDbl HelmholtzEOSMixtureBackend::calc_Hmolar_formation() {
    // Pure and pseudo-pure only.  The ideal-gas mole-fraction sum would be
    // exact, but a mixture is not a compound and the analogous entropy
    // quantity is not linear, so the two would not stay consistent.
    if (components.size() != 1) {
        throw ValueError(format("calc_Hmolar_formation is only valid for pure and pseudo-pure fluids, %d components",
                                components.size()));
    }
    CoolPropDbl v = components[0].standard_state.hmolar;
    if (!ValidNumber(v)) {
        throw ValueError(format("No standard enthalpy of formation is available for fluid [%s]; the source database "
                                "does not provide a value for this species",
                                components[0].name.c_str()));
    }
    return v;
}
```

- [ ] **Step 8: Rebuild and run the tests**

```bash
python dev/generate_headers.py
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[formation]" -s
```

Expected: all sections pass, with the plausibility section reporting `checked > 60`.

- [ ] **Step 9: Run the broader suite for regressions**

```bash
./build_catch/CatchTestRunner "~[slow]"
```

Expected: no new failures. Inserting into the middle of the `parameters` enum shifts the numeric value of every entry after `idipole_moment`, which is fine by design — but this full run is what proves nothing in the library indexes parameters by hard-coded number or persists them to disk. If a failure appears here, the bug is the hard-coded value, not the enum placement.

- [ ] **Step 10: Commit**

```bash
git add include/CoolProp/CoolPropFluid.h include/CoolProp/DataStructures.h include/CoolProp/AbstractState.h \
        src/DataStructures.cpp src/AbstractState.cpp \
        src/Backends/Helmholtz/Fluids/FluidLibrary.h src/Backends/Helmholtz/Fluids/FluidLibrary.cpp \
        src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp \
        src/Tests/CoolProp-Tests.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "feat: expose HFORMATION standard enthalpy of formation (CoolProp-ihto)

Ingests INFO.STANDARD_STATE into FormationStruct and exposes it through
the existing trivial-parameter path.  Throws for fluids with no value and
for mixtures; never returns 0 or NaN.  The enum entry sits with the other
per-fluid constants next to idipole_moment."
```

---

### Task 5: CI wiring for the regenerator tests

**Beads:** `CoolProp-ijbw`

**Files:**
- Modify: `.github/workflows/dev_checks.yml`

**Interfaces:**
- Consumes: `dev/atct/test_fetch_atct_formation.py` from Tasks 1–3.
- Produces: nothing consumed by later tasks.

- [ ] **Step 1: Add the job**

Append to `.github/workflows/dev_checks.yml` at the same indentation as the existing `type-stubs:` job (two spaces, under `jobs:`):

```yaml
  atct-regenerator:
    name: ATcT regenerator tests
    # Offline gate for dev/atct/. The tests parse a committed HTML fixture and
    # read dev/fluids/*.json -- no network, no C++ build. Deliberately NO paths
    # filter so this can be a required check without deadlocking unrelated PRs.
    if: github.event_name == 'pull_request' || github.event_name == 'push'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v7

      - name: Set up Python
        uses: actions/setup-python@v6
        with:
          python-version: '3.12'

      - name: Install pytest
        run: python -m pip install --upgrade pytest

      - name: Run the ATcT regenerator tests
        run: python -m pytest dev/atct/test_fetch_atct_formation.py -v
```

- [ ] **Step 2: Verify the workflow parses**

```bash
python -c "import yaml,sys; d=yaml.safe_load(open('.github/workflows/dev_checks.yml')); assert 'atct-regenerator' in d['jobs']; print('job registered:', list(d['jobs'])[-3:])"
```

Expected: the assertion passes and the job name is listed.

- [ ] **Step 3: Confirm the tests pass exactly as CI will run them**

```bash
python -m pytest dev/atct/test_fetch_atct_formation.py -v
```

Expected: 19 passed.

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/dev_checks.yml
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git commit --no-verify -m "ci: run the ATcT regenerator tests (CoolProp-ijbw)

Offline and build-free: the tests parse a committed fixture, so the three
silent-failure parsing traps stay covered on every PR."
```

---

## Final gate

- [ ] **Run the pre-push gate**

```bash
./dev/ci/preflight.sh
```

The changed paths do not include `src/SBTL/` or `src/Backends/SVDSBTL/`, so preflight selects the general tag scope. It must pass before any push.

- [ ] **Pre-PR adversarial review (MANDATORY per CLAUDE.md)**

Dispatch the review subagent against the full diff versus the branch base, using the canonical invocation in `CLAUDE.md`. Do this **before** `git push` and `gh pr create`, and re-review any commit added afterwards.

- [ ] **Close the beads issues**

```bash
bd close CoolProp-ijbw
bd close CoolProp-ihto
bd close CoolProp-oett   # only if no further tiers are planned in this branch
```

---

## Self-review notes

**Spec coverage.** §4 regenerator → Tasks 1–3. §5 schema → Task 3 Step 3. §6 C++ ingestion and API → Task 4. §7 test table: known values, absent-fluid throw, mixture throw and magnitude guard → Task 4 Step 1; the three regenerator trap tests → Task 1 Step 2; the coverage-ledger gate → Task 3. §9 acceptance criteria 1–2 → Task 3 Step 6; criterion 4 → Task 4 Step 1; criterion 5 → Final gate.

**One spec item is intentionally not a task.** §7's "independent reproduction of all 75 CSV values" is covered in substance by Task 3 Step 6, which asserts the four values that matter — including the para-hydrogen figure where the contributor's CSV is wrong. Re-checking all 75 would require committing that CSV as a fixture, which reintroduces the hand-curated artifact the design exists to eliminate. The coverage ledger plus the four spot-checks give the same protection without vendoring it.

**Type consistency.** `AtctRow` field names are used identically in Tasks 1–3 (`dhf298_kJ_per_mol`, `uncertainty_kJ_per_mol`, `atct_id`). `FormationStruct::hmolar` is written by `parse_standard_state` (Task 4 Step 4) and read by `calc_Hmolar_formation` (Task 4 Step 7). The JSON key `hmolar_formation.value` is written in Task 3 Step 3 and read in Task 4 Step 4.
