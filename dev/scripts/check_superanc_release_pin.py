#!/usr/bin/env python3
"""
Gate the fastchebpure release pinned by the docs superancillary plots against
CoolProp's current EOS.

The deviation plots in `Web/scripts/fluid_properties.Superancillary.py`
download a *pinned* fastchebpure release and divide CoolProp's superancillary
by that release's extended-precision reference data.  fastchebpure generates
that reference straight from CoolProp's `dev/fluids/*.json` at the commit its
`externals/CoolProp` submodule was pinned to when the release was tagged.  If
an EOS is updated in CoolProp *after* that pin -- the chicken-and-egg case
where the EOS (and its freshly-fit SA) lands here before fastchebpure can be
regenerated and re-tagged -- the pinned release describes a *different*
equation of state for that fluid, and the deviation plot becomes meaningless
(this is exactly issue #3063 / R1233zd(E): a ~4 % "deviation" that was really
the old Mondejar EOS vs the new international-standard EOS).

This is the EOS<->docs-reference twin of `check_superanc_freshness.py`, which
covers EOS<->embedded-SA.  Both speak the same `source_eos_hash` contract --
the FNV-1a structural hash of `EOS[0]` with the `SUPERANCILLARY` subtree
removed, defined in `inject_superanc_check_points.py::eos_fnv1a_hex`.
fastchebpure stamps that hash into every `output/{fluid}_exps.json` from
release `2026.04.23` onward (fastchebpure commit dcdd88c); we recompute it from
the current `dev/fluids` EOS and require equality.

`PENDING_UPSTREAM` lists fluids knowingly awaiting an upstream fastchebpure
regeneration.  A hash mismatch for an allowlisted fluid is reported but does
NOT fail the gate -- otherwise the transition window (EOS lands here, upstream
not yet re-tagged) would deadlock every PR.  Each entry MUST cite a tracking
issue and should be deleted once the bumped release is pinned; the gate warns
loudly if an allowlisted fluid has since started matching (stale exemption).

Fluids the pinned release does not ship at all, or ships without a stamped
hash (a pre-2026.04.23 release format), are reported as warnings, not
failures: the docs script already degrades those to an annotated placeholder
rather than a wrong curve, so they cannot mislead.

The gate also requires every OTHER `outputversion = '...'` pin under Web/
(currently the SuperAncillary.ipynb docs notebook, which downloads the same
release archive independently) to equal the docs-script pin.  The notebook is
executed at docs-build time with `--allow-errors`, so a stale pin there fails
soft: the build succeeds and the rendered plots silently compare against an
old release (the same #3063 failure mode this gate exists to prevent).

Exit status 0 if every covered fluid matches (or is allowlisted) and all
secondary pins agree; 1 if any non-allowlisted fluid's pinned reference is
stale or a secondary pin has drifted.
"""
import os
import re
import sys
import json
import zipfile
import tempfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from inject_superanc_check_points import eos_fnv1a_hex  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
FLUIDS = REPO / 'dev' / 'fluids'
DOCS_SCRIPT = REPO / 'Web' / 'scripts' / 'fluid_properties.Superancillary.py'

# Fluids knowingly awaiting an upstream fastchebpure regeneration.  Map the
# dev/fluids stem to a human-readable reason that cites the tracking issue.
# A mismatch here is reported but does not fail the gate.  Keep this empty
# whenever the pinned release is fully current.
PENDING_UPSTREAM = {
    # 'R1243zf': 'awaiting fastchebpure regen, CoolProp/fastchebpure#6',
}


def pinned_version():
    """Parse `outputversion = '...'` out of the docs plot script so this gate
    always tracks exactly the release the docs actually download."""
    text = DOCS_SCRIPT.read_text()
    m = re.search(r"""^outputversion\s*=\s*['"]([^'"]+)['"]""", text, re.M)
    if not m:
        raise RuntimeError(f'could not find outputversion in {DOCS_SCRIPT}')
    return m.group(1)


# Files under Web/ that are KNOWN to carry their own `outputversion` pin.
# Each must yield at least one parseable pin equal to the docs-script pin;
# a file listed here in which the scan finds NO pin at all fails the gate
# too (a rename/f-string refactor/deleted cell must not silently drop the
# file out of coverage -- that is the fail-open this check exists to close).
# If a pin is legitimately removed (e.g. derived from the docs script
# instead), delete the entry here in the same commit.
EXPECTED_SECONDARY_PIN_FILES = [
    Path('Web') / 'coolprop' / 'SuperAncillary.ipynb',
]


def secondary_pin_drift(version):
    """Return [(path, problem), ...] for every `outputversion = ...`
    assignment found in notebooks or scripts under Web/ (outside the
    canonical docs script) whose value does not equal the docs-script pin,
    plus an entry for every EXPECTED_SECONDARY_PIN_FILES member in which no
    pin could be parsed at all.

    Notebooks embed their source as JSON strings, so the assignment is
    matched unanchored against the raw file text (an `^` anchor would never
    match mid-string), with optional escaped quotes (a double-quoted pin
    appears as \\" in raw notebook text) and an optional 1-2 letter string
    prefix (f/r/rb...).  The captured value may be empty or non-numeric --
    a malformed pin must FAIL the equality check, not evade the match."""
    drift = []
    found_in = set()
    pat = re.compile(r"""outputversion\s*=\s*[A-Za-z]{0,2}(?:\\?['"])([^'"\\]*)(?:\\?['"])""")
    for path in sorted((REPO / 'Web').rglob('*')):
        if path.suffix not in ('.ipynb', '.py') or '.ipynb_checkpoints' in path.parts:
            continue
        if path == DOCS_SCRIPT or path.resolve() == DOCS_SCRIPT.resolve():
            continue
        rel = path.relative_to(REPO)
        for found in pat.findall(path.read_text(encoding='utf-8')):
            found_in.add(rel)
            if found != version:
                drift.append((rel, repr(found)))
    for rel in EXPECTED_SECONDARY_PIN_FILES:
        if rel not in found_in:
            drift.append((rel, 'no parseable outputversion pin found '
                               '(expected one; see EXPECTED_SECONDARY_PIN_FILES)'))
    return drift


def release_data(version):
    """Return (check_stems, exps_hash_ci) for the pinned fastchebpure release.

    check_stems: the exact set of stems shipped as outputcheck/{stem}_check.json.
        This is the set the docs plot actually consumes -- it renders a real
        deviation curve iff `{fluid}_check.json` exists (case-sensitively, as the
        docs script does `os.path.exists`), otherwise an annotated placeholder.
        Keying coverage on this set (not on the expansion set) is what makes the
        gate faithful: it verifies exactly the fluids the docs actually plot.

    exps_hash_ci: {stem.lower(): source_eos_hash_or_None} from the expansion
        files output/{stem}_exps.json (the hash is stamped here, not in the
        check file).  Keyed case-folded because fastchebpure's _exps and _check
        casing can diverge for the same fluid (e.g. output/R1234YF_exps.json
        alongside outputcheck/R1234yf_check.json); an exact-case lookup would
        miss the hash and silently drop a docs-covered fluid into the "absent"
        bucket -- the very false-PASS this gate must avoid.

    Downloads the release zip once, atomically (cached in the system temp dir so
    repeated local runs do not re-fetch ~40 MB, and an interrupted download
    cannot leave a truncated zip that later runs would trust)."""
    cache = Path(tempfile.gettempdir()) / f'fastchebpure-{version}.zip'
    if not cache.exists():
        url = f'https://github.com/CoolProp/fastchebpure/archive/refs/tags/{version}.zip'
        print(f'Downloading pinned fastchebpure release {version} ...')
        part = cache.with_suffix('.zip.part')
        urllib.request.urlretrieve(url, part)
        os.replace(part, cache)
    check_stems = set()
    exps_hash_ci = {}
    with zipfile.ZipFile(cache) as z:
        for name in z.namelist():
            if '/outputcheck/' in name and name.endswith('_check.json'):
                check_stems.add(name.rsplit('/', 1)[1][:-len('_check.json')])
            elif '/output/' in name and name.endswith('_exps.json'):
                stem = name.rsplit('/', 1)[1][:-len('_exps.json')]
                exps_hash_ci[stem.lower()] = json.loads(z.read(name)).get('source_eos_hash')
    return check_stems, exps_hash_ci


def main():
    version = pinned_version()
    check_stems, exps_hash_ci = release_data(version)

    matched = []
    stale = []        # (fluid, current_hash, release_hash) -- fails the gate
    allowlisted = []  # (fluid, current_hash, release_hash) -- reported only
    absent = []       # fluid has SA here but pinned release does not ship it
    unhashed = []     # release ships fluid but with no source_eos_hash
    stale_exempt = [] # in PENDING_UPSTREAM yet now matching -> drop the entry

    for path in sorted(FLUIDS.glob('*.json')):
        d = json.loads(path.read_text())
        if 'EOS' not in d or not d['EOS']:
            continue
        eos = d['EOS'][0]
        if 'SUPERANCILLARY' not in eos:
            continue
        fluid = path.stem
        current = eos_fnv1a_hex(eos)
        # Coverage is decided by the check-file set the docs actually consume,
        # not the expansion set: a fluid the docs plot for real must be verified.
        if fluid not in check_stems:
            absent.append(fluid)
            continue
        # Hash lives in the expansion file; look it up case-folded so a casing
        # divergence (R1234yf_check.json vs R1234YF_exps.json) still verifies.
        stamped = exps_hash_ci.get(fluid.lower())
        if stamped is None:
            unhashed.append(fluid)
            continue
        if stamped == current:
            matched.append(fluid)
            if fluid in PENDING_UPSTREAM:
                stale_exempt.append(fluid)
        elif fluid in PENDING_UPSTREAM:
            allowlisted.append((fluid, current, stamped))
        else:
            stale.append((fluid, current, stamped))

    print(f'Docs-pinned fastchebpure release: {version}')
    print(f'Checked {len(matched) + len(stale) + len(allowlisted)} hashed '
          f'SA-bearing fluids against the pinned release.')

    if absent:
        print(f'\nNOTE: {len(absent)} SA-bearing fluid(s) are not shipped by the '
              'pinned release (docs render the placeholder for these):')
        print('  ' + ', '.join(absent))
    if unhashed:
        print(f'\nNOTE: {len(unhashed)} fluid(s) shipped without a source_eos_hash '
              '(pre-2026.04.23 release format); cannot verify:')
        print('  ' + ', '.join(unhashed))
    if allowlisted:
        print(f'\nALLOWLISTED ({len(allowlisted)} stale but pending upstream, '
              'not failing the gate):')
        for fluid, cur, rl in allowlisted:
            print(f'  {fluid}: current={cur} release={rl}  [{PENDING_UPSTREAM[fluid]}]')
    if stale_exempt:
        print(f'\nWARNING: {len(stale_exempt)} fluid(s) are in PENDING_UPSTREAM but '
              'now MATCH the pinned release -- remove their stale exemption:')
        print('  ' + ', '.join(stale_exempt))

    drift = secondary_pin_drift(version)
    if drift:
        print(f'\nDRIFTED SECONDARY PIN(S): {len(drift)} other outputversion '
              f'pin(s) under Web/ disagree with the docs-script pin ({version}):')
        for path, found in drift:
            print(f'  {path}: {found}')
        print('These download the fastchebpure release independently of the docs '
              'script; a stale pin renders plots against an old release without '
              'failing the docs build (the notebook executes with --allow-errors). '
              'Bump them to match.')

    if stale:
        print(f'\nSTALE: {len(stale)} fluid(s) whose pinned docs reference was '
              'built from a different EOS than CoolProp currently ships:')
        for fluid, cur, rl in stale:
            print(f'  {fluid}: current={cur}  release={rl}')
        print()
        print('The deviation plots for these fluids would compare CoolProp\'s '
              'superancillary against a different equation of state (cf. #3063).')
        print('Fix: bump the pin in Web/scripts/fluid_properties.Superancillary.py '
              'to a fastchebpure release regenerated against the current EOS, or -- '
              'if the regeneration is still pending upstream -- add the fluid(s) to '
              'PENDING_UPSTREAM with a tracking-issue reference.')

    if stale or drift:
        return 1

    print('\nAll hashed, shipped SA-bearing fluids match the pinned release EOS, '
          'and all secondary pins agree.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
