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

import argparse
import csv
import hashlib
import html as html_module
import json
import re
import sys
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

# dev/package_json.py is repo-local, not an installed package; put dev/ on
# sys.path before importing it so json_options (the canonical serialization
# options) is available without duplicating them here.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from package_json import json_options  # noqa: E402

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
_NAME_CELL_RE = re.compile(r'class="bkgName">(?P<cell>.*?)</td>', re.S)
_TAG_RE = re.compile(r"<[^>]+>")
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


def _extract_name(body: str, cas: str) -> str:
    """Extract the species name from the bkgName cell's text content.

    ATcT's own page documentation: "For species other than elements in
    their standard states, the Species Name acts as a link" -- so elements
    in their standard state (Boron, Sulfur, Br2, I2, ...) render as plain
    text with no <a> anchor, e.g. '<span class="Name">Boron</span>', while
    every other species wraps the name in a species-detail link.  Extracting
    text from the whole bkgName cell -- tolerating the anchor when present
    rather than requiring it -- covers both shapes; anchoring the extraction
    on the <a> tag itself (an earlier version of this parser did that)
    mishandles the unlinked element rows, which the small fixture sample
    never happened to exercise.

    Matches on the 'bkgName' CSS class specifically (not just any
    class="Name" span) because the bkgImage cell later in the same row also
    carries a nested <span class="Name">...</span> around its <img> alt
    text; matching the first "Name" span anywhere in the row would pick
    that up instead on rows this regex is not anchored to a td boundary.
    """
    cell_match = _NAME_CELL_RE.search(body)
    if cell_match is None:
        raise AtctParseError("missing name cell on CAS %s" % cas)
    text = html_module.unescape(_TAG_RE.sub("", cell_match.group("cell"))).strip()
    if not text:
        raise AtctParseError("empty name cell on CAS %s" % cas)
    return text


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
        name = _extract_name(body, cas)
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
        # atct.anl.gov 403s urllib's default "Python-urllib/x.y" User-Agent
        # (ordinary Cloudflare bot-defense on an unidentified client); a
        # plain, honest identifying UA is standard practice for a script
        # fetching a public data page, not a bypass of any access control.
        request = urllib.request.Request(
            url, headers={"User-Agent": "CoolProp-atct-fetcher/1.0 (+https://github.com/CoolProp/CoolProp)"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
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
        # --update-ledger overwrites the pinned coverage, so it is exactly
        # the invocation someone reaches for on a version bump -- the brief's
        # own bootstrap command pairs it with --write in one shot.  Diffing
        # against whatever ledger is about to be overwritten and printing
        # every difference (informational: this still exits 0 and still
        # writes) is what gives the two-step review discipline described in
        # the README something a reviewer can actually check, instead of a
        # silent overwrite that leaves no trace in the console, the PR diff,
        # or CI logs if a fluid quietly flips matched -> absent.
        if ledger_path.exists():
            changes = compare_ledger(json.loads(ledger_path.read_text(encoding="utf-8")), actual)
            if changes:
                print("INFO: --update-ledger is overwriting %d changed entr%s:"
                      % (len(changes), "y" if len(changes) == 1 else "ies"))
                for change in changes:
                    print("INFO:   %s" % change)
            else:
                print("INFO: --update-ledger: no changes from the existing ledger")
        else:
            print("INFO: no existing ledger at %s; writing %d entries" % (ledger_path, len(actual)))
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
