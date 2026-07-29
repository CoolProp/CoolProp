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
