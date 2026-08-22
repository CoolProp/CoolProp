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
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path

# dev/package_json.py is repo-local, not an installed package; put dev/ on
# sys.path before importing it so json_options (the canonical serialization
# options) is available without duplicating them here.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from package_json import json_options  # noqa: E402

# dev/fluids/*.json is not key-sorted below the top level, so writing a file
# with package_json's json_options (which sets sort_keys=True) reorders every
# pre-existing key as a side effect -- 2608 lines of pure churn across 76 files,
# burying the one block this tool actually adds.  Preserve the file's existing
# key order instead and append the new block, so the diff is the addition and
# nothing else.  indent is still taken from json_options so the two agree on
# formatting; only the ordering differs.
#
# This is local to this tool, NOT a settled repo convention: dev/package_json.py,
# dev/scripts/inject_InChI.py, dev/scripts/fit_rational_functions.py and
# dev/scripts/convert_FLD.py all still write dev/fluids with sort_keys=True,
# while dev/scripts/inject_superancillary.py writes indent=2 with no sort_keys
# and ran last over ~127 files -- which is why the tree is unsorted in the
# first place.  Any of those reintroduces the churn.  Reconciling them is a
# separate change; re-sorting the tree without doing so would just regress on
# the next superancillary injection.
_WRITE_OPTIONS = {**json_options, "sort_keys": False}

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


def _is_exact(text: str) -> bool:
    """ATcT's marker for elements in their standard state.

    One definition shared by the units gate and the uncertainty parser.  They
    previously disagreed about whether a stray '±' counted as part of the
    token; harmless only because the divergence happened to be fail-closed.
    """
    return text.replace("±", "").strip().lower() == "exact"


def _parse_uncertainty(text: str) -> float:
    """ATcT writes '± 0.043' for measured species and 'exact' for elements.

    The element form carries no '±' and an empty units cell; treating it as a
    parse failure silently drops all nine elements from the output.
    """
    if _is_exact(text):
        return 0.0
    cleaned = text.replace("±", "").strip()
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
        uncertainty_text = _span(body, "Uncert")
        # An empty units cell is legitimate ONLY on the "exact" element rows,
        # which carry no units at all.  _span() also returns "" when the span
        # is missing outright -- a renamed CSS class, a restyled table -- so
        # accepting "" unconditionally would make this check, whose entire
        # purpose is to catch a units change, pass for every row on exactly
        # the page change most likely to accompany one.  Tie the exemption to
        # the uncertainty text that justifies it.
        if units == "":
            if not _is_exact(uncertainty_text):
                raise AtctParseError(
                    "missing units on CAS %s (uncertainty %r); the units cell is only "
                    "allowed to be empty on 'exact' element rows" % (cas, uncertainty_text)
                )
        elif units != "kJ/mol":
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
                uncertainty_kJ_per_mol=_parse_uncertainty(uncertainty_text),
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


# The reference-state scaffolding standard_state_block() writes around the
# quantity itself; shared with clear_standard_state() so the two agree on what
# belongs to this tool.
_REFERENCE_STATE_KEYS = frozenset({"T", "T_units", "p", "p_units", "phase"})

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


def _load_ordered(path: Path):
    """Read a fluid file preserving key order, and report its line terminator.

    Returns (doc, trailing_newline).  The terminator comes back with the
    document deliberately: capturing it separately is how R1130(E).json lost
    its newline, and a helper that returns only the doc is exactly what the
    next writer would reach for.
    """
    text = path.read_text(encoding="utf-8")
    return json.loads(text, object_pairs_hook=OrderedDict), text.endswith("\n")


def write_standard_state(path: Path, row: AtctRow, version: str) -> None:
    """Insert the block, preserving the file's existing key order.

    Appends STANDARD_STATE as the last INFO key and rewrites with
    sort_keys=False, so every other line of the file is untouched and the diff
    is exactly the added block.  Two files are not perfectly round-trippable
    for reasons that predate this tool -- Chlorine.json carries a duplicate
    "source_eos_hash" key, and Chlorine/R1130(E) hold float literals with more
    digits than Python's shortest round-trip repr -- so those pick up a few
    incidental normalization lines.
    """
    doc, trailing_newline = _load_ordered(path)
    block = doc["INFO"].get("STANDARD_STATE")
    fresh = standard_state_block(row, version)
    if block is None:
        doc["INFO"]["STANDARD_STATE"] = fresh
    elif not isinstance(block, dict):
        # Refuse hand-edited garbage rather than overwriting it, matching
        # clear_standard_state.  main() gates this, so reaching here is a bug.
        raise AtctParseError("%s has a non-object INFO.STANDARD_STATE; refusing to overwrite it" % path.name)
    else:
        # MERGE, do not replace: assigning the whole object deletes every other
        # quantity sharing STANDARD_STATE, and the design commits to sharing it
        # (spec section "Shared state at the top" -- a later entropy tier drops
        # a smolar sibling in).
        #
        # But merging `fresh` wholesale rewrote the SHARED REFERENCE STATE too.
        # A block holding another source's quantity at p=101325 silently became
        # p=100000, leaving THEIR value annotated with OUR reference state --
        # worse than deletion, because it still reads as valid data.  The gate
        # could not see it: it only ever inspected hmolar_formation.
        #
        # Write our quantity, and only fill scaffolding that is absent.  A
        # conflicting value is left exactly as found; main() refuses the run so
        # a human resolves it, because there is no correct automatic answer
        # when two sources disagree about the reference state.
        block["hmolar_formation"] = fresh["hmolar_formation"]
        for key in _REFERENCE_STATE_KEYS:
            block.setdefault(key, fresh[key])
    _write_doc(path, doc, trailing_newline)


def _formation_source(path: Path):
    """Read-only: the `source` of this file's hmolar_formation, if it has one.

    Returns (has_formation, source).  `source` is None both when the key is
    absent and when there is no block at all, so callers must check the flag.
    """
    block = json.loads(path.read_text(encoding="utf-8"))["INFO"].get("STANDARD_STATE")
    if not isinstance(block, dict):
        return False, None
    hf = block.get("hmolar_formation")
    if not isinstance(hf, dict):
        # A scalar or string here is hand-edited garbage; report it as present
        # but foreign so the gates refuse to touch it rather than crash on
        # .get().
        return (hf is not None), "<malformed>"
    return True, hf.get("source")


def _is_ours(source) -> bool:
    """Did this tool write the value?

    A block with no `source` is treated as ours: every block this tool has
    ever written carries source="ATcT", so an unsourced one predates the field
    or was hand-copied from our output.  A block naming a different source is
    someone else's and is never touched without --allow-removals.
    """
    return source in (None, "ATcT")


def owns_formation(path: Path) -> bool:
    """Read-only: does this file carry an hmolar_formation THIS tool may remove?

    Mirrors clear_standard_state()'s ownership test exactly, so main() can
    decide about removals before writing anything.  A broader predicate here
    would make --write demand --allow-removals for deletions that would never
    happen -- which teaches the operator to pass the flag by reflex, and that
    is how a guard fails open in practice.
    """
    has, source = _formation_source(path)
    return has and _is_ours(source)


def reference_state_conflicts(path: Path, row: AtctRow, version: str) -> dict:
    """Read-only: scaffolding keys whose existing value disagrees with ours.

    The reference state is shared by every quantity in the block, so writing
    ours over a differing one silently re-annotates somebody else's data.
    Returns {key: (existing, ours)}; empty when there is nothing to argue about.
    """
    block = json.loads(path.read_text(encoding="utf-8"))["INFO"].get("STANDARD_STATE")
    if not isinstance(block, dict):
        return {}
    fresh = standard_state_block(row, version)
    return {k: (block[k], fresh[k]) for k in _REFERENCE_STATE_KEYS if k in block and block[k] != fresh[k]}


def has_unwritable_block(path: Path) -> bool:
    """Read-only: is INFO.STANDARD_STATE present but not an object?"""
    info = json.loads(path.read_text(encoding="utf-8"))["INFO"]
    return "STANDARD_STATE" in info and not isinstance(info["STANDARD_STATE"], dict)


def foreign_formation(path: Path):
    """Read-only: the foreign `source` owning this file's hmolar_formation, else None."""
    has, source = _formation_source(path)
    return source if (has and not _is_ours(source)) else None


def _write_doc(path: Path, doc, trailing_newline: bool) -> None:
    """Serialize, restoring the file's original trailing-newline state.

    7 of the 137 fluid files end with a newline and the rest do not; dropping
    it turns "insert a block" into "insert a block and strip the terminator",
    which is exactly the kind of incidental change this tool exists to avoid.
    """
    text = json.dumps(doc, **_WRITE_OPTIONS)
    path.write_text(text + "\n" if trailing_newline else text, encoding="utf-8")


def clear_standard_state(path: Path) -> bool:
    """Drop INFO.STANDARD_STATE from a fluid the current version does not cover.

    Returns True if a block was actually removed.

    Writing is matched-only, so without this a fluid that flips matched ->
    absent between ATcT versions keeps serving its previous value: the ledger
    and atct_report.csv would say "absent" while dev/fluids/<Fluid>.json still
    carried a block stamped with the OLD version, and HFORMATION would keep
    returning it forever.  The ledger gate cannot catch that -- by the time
    --write runs, the ledger has already been brought into agreement with the
    new bind result.  Removal has to happen on the same pass that writes.
    """
    doc, trailing_newline = _load_ordered(path)
    block = doc["INFO"].get("STANDARD_STATE")
    if not isinstance(block, dict) or "hmolar_formation" not in block:
        return False
    if not isinstance(block["hmolar_formation"], dict):
        return False  # hand-edited garbage; leave it for a human
    # Remove only THIS tool's quantity, and only if ATcT wrote it.  The README's
    # Scope section plans an entropy tier from a different source; deleting the
    # whole STANDARD_STATE object would take that source's data with it for
    # every fluid ATcT happens not to cover.
    if not _is_ours(block["hmolar_formation"].get("source")):
        return False
    del block["hmolar_formation"]
    # What is left is either nothing, or the reference-state scaffolding that
    # standard_state_block() writes alongside the quantity.  If no OTHER
    # quantity is using it, the wrapper is ours too and goes with it -- so
    # insert-then-clear restores the file byte for byte.  If a future tier has
    # added its own key, keep the wrapper and leave that data alone.
    if not set(block) - _REFERENCE_STATE_KEYS:
        del doc["INFO"]["STANDARD_STATE"]
    _write_doc(path, doc, trailing_newline)
    return True


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
    parser.add_argument("--allow-removals", action="store_true",
                        help="authorize destructive changes, all refused by default: deleting an "
                             "ATcT formation value from a fluid this run reports absent, "
                             "overwriting one owned by another source, or overwriting a "
                             "non-object STANDARD_STATE. Does NOT authorize rewriting a "
                             "conflicting reference state, which is never automatic.")
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

    # ---- destructive-change gate -------------------------------------------
    # This runs BEFORE the ledger branch, not just before the fluid writes.
    # --update-ledger rewrites expected_coverage.json, and that file is the one
    # artifact whose job is to make the NEXT run fail loudly.  Overwriting it
    # during a run that then refuses to proceed would silence the loudest
    # tripwire: a later plain --write would sail past compare_ledger because
    # the ledger now agrees with the degraded page.  So a refused run must
    # leave the ledger untouched too -- "nothing was written" has to be true.
    ledger_path = here / "expected_coverage.json"
    # Gate BOTH destructive modes.  --update-ledger alone rewrites the pinned
    # coverage, and that file is the tripwire everything else leans on, so it
    # is not a read-only invocation.
    if args.write or args.update_ledger:
        # Removals: absent fluids whose file still holds an ATcT value.
        to_clear = [name for name in sorted(result.absent) if owns_formation(fluids[name].path)]
        # Overwrites: matched fluids holding SOMEONE ELSE'S formation value.
        foreign = [(n, foreign_formation(fluids[n].path)) for n in sorted(result.matched)]
        foreign = [(n, s) for n, s in foreign if s is not None]
        # Unwritable: a hand-edited non-object block the writer must not clobber.
        garbage = [n for n in sorted(result.matched) if has_unwritable_block(fluids[n].path)]
        # Reference-state disagreements: the block is SHARED, so overwriting
        # T/p/phase re-annotates another source's quantity.  There is no correct
        # automatic resolution, so this is refused even with --allow-removals.
        conflicts = []
        if args.write:
            for name in sorted(result.matched):
                clash = reference_state_conflicts(fluids[name].path, result.matched[name], args.version)
                if clash:
                    conflicts.append((name, clash))

        blocked = (to_clear or foreign or garbage) and not args.allow_removals
        if blocked or conflicts:
            if to_clear and not args.allow_removals:
                print("ERROR: %d fluid(s) carry an ATcT formation value but are absent from this"
                      " run; refusing to delete committed values." % len(to_clear), file=sys.stderr)
                for name in to_clear:
                    print("ERROR:   delete %s (%s)" % (name, result.absent[name]), file=sys.stderr)
            if foreign and not args.allow_removals:
                print("ERROR: %d matched fluid(s) already carry an hmolar_formation from another"
                      " source; refusing to overwrite it." % len(foreign), file=sys.stderr)
                for name, src in foreign:
                    print("ERROR:   overwrite %s (source %r)" % (name, src), file=sys.stderr)
            if garbage and not args.allow_removals:
                print("ERROR: %d matched fluid(s) have a non-object INFO.STANDARD_STATE;"
                      " refusing to overwrite it." % len(garbage), file=sys.stderr)
                for name in garbage:
                    print("ERROR:   overwrite %s" % name, file=sys.stderr)
            if conflicts:
                print("ERROR: %d matched fluid(s) declare a reference state that disagrees with"
                      " ATcT's; another quantity in the same block is annotated by it, so this"
                      " cannot be resolved automatically -- fix the file by hand."
                      % len(conflicts), file=sys.stderr)
                for name, clash in conflicts:
                    for key, (existing, ours) in sorted(clash.items()):
                        print("ERROR:   %s: %s is %r, ATcT uses %r" % (name, key, existing, ours), file=sys.stderr)
            hint = ("" if conflicts else
                    "  If these changes are genuinely correct, re-run with --allow-removals.")
            print("ERROR: nothing was written -- not the fluid files, not the report, and not %s.%s"
                  % (ledger_path, hint), file=sys.stderr)
            return 1
    else:
        to_clear = []

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
        cleared = []
        for name in to_clear:
            if clear_standard_state(fluids[name].path):
                # stderr, and only when something actually went: printing
                # unconditionally produced "REMOVED ... X" next to "cleared 0
                # stale blocks" in the same run.
                cleared.append(name)
                print("REMOVED the ATcT formation value from %s (%s)"
                      % (name, result.absent[name]), file=sys.stderr)
        write_report(here / "atct_report.csv", result, args.version, page_sha256)
        print("wrote %d fluid files, cleared %d stale block%s, and atct_report.csv"
              % (len(result.matched), len(cleared), "" if len(cleared) == 1 else "s"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
