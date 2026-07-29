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
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

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
_PHASE_RE = re.compile(r"^(?P<formula>.+?)\s+(?P<phase>\([^()]*\))$")


def _span(body: str, css_class: str) -> str:
    match = re.search(r'class="%s">(.*?)</span>' % css_class, body, re.S)
    if match is None:
        return ""
    return html_module.unescape(match.group(1)).strip()


def _split_phase(raw: str) -> tuple[str, str]:
    """'CH4  (g)' -> ('CH4', '(g)');  'H2 (g, para)' -> ('H2', '(g, para)');
    '(CH3)2NH (g)' -> ('(CH3)2NH', '(g)').

    Anchors on the trailing whitespace-separated parenthesized group rather
    than the first '(' in the string.  Some real ATcT formulas themselves
    begin with '(' -- '(CH3)2NH (g)', '(F)(HF) (g, lin)' -- and splitting on
    the first '(' mis-parses those into an empty formula and a bogus phase
    string, which then fails both the exact-'(g)' and '(g,'-prefix checks in
    select_gas_row() and makes the species silently vanish as 'absent'.
    """
    collapsed = " ".join(raw.split())
    match = _PHASE_RE.match(collapsed)
    if match is None:
        return collapsed, ""
    return match.group("formula").strip(), match.group("phase").strip()


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
        cas = match.group("cas").strip()
        if units not in ("", "kJ/mol"):
            raise AtctParseError("unexpected units %r on CAS %s" % (units, cas))
        name_match = _NAME_RE.search(body)
        if name_match is None:
            # A blank name would otherwise degrade quietly into a row that
            # still parses and still gets selected -- fail loud instead.
            raise AtctParseError("missing species-name anchor on CAS %s" % cas)
        name = html_module.unescape(name_match.group("name")).strip()
        formula, phase = _split_phase(html_module.unescape(formula_match.group("formula")))
        atct_id = _span(body, "ATcTID")
        if not atct_id:
            raise AtctParseError("missing ATcT ID on CAS %s" % cas)
        rows.append(
            AtctRow(
                cas=cas,
                name=name,
                formula=formula,
                phase=phase,
                dhf298_kJ_per_mol=float(value_text),
                uncertainty_kJ_per_mol=_parse_uncertainty(_span(body, "Uncert")),
                atct_id=atct_id,
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

    Requires that every row shares one CAS number.  "CAS is not a unique
    key" is the exact trap this function exists to defend against (argon's
    (g) and (aq, undissoc) rows share a CAS, for instance), so a caller
    accidentally handing it rows for two different species is a violated
    precondition, not a legitimate multi-row ambiguity -- it is raised as
    AtctParseError (an invalid-input-shape error) rather than
    AmbiguousSpecies (which models a real choice among candidate gas rows
    for one species that cannot be resolved).
    """
    cas_values = {r.cas for r in rows}
    if len(cas_values) > 1:
        raise AtctParseError(
            "select_gas_row() requires rows for a single CAS number; got %s"
            % ", ".join(sorted(cas_values))
        )

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
