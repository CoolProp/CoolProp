#!/usr/bin/env python3
"""
Verify that each fluid's stored SUPERANCILLARY.source_eos_hash still matches
the hash of its current EOS[0] (with the SUPERANCILLARY subtree excluded).

`source_eos_hash` is stamped at SA-generation time by fastchebpure (emitted
under the same name in its output/*_exps.json) and copied verbatim by
`fitcheb inject` and dev/scripts/inject_superanc_check_points.py. If anyone
has edited the fluid's EOS since (gas constant, alpha0/alphar coefficients,
reducing state, ...), the current hash no longer matches, and the SA is
stale — bump the fastchebpure pin to regenerate.

This check is a thin Python twin of the C++ test "Superancillary
source_eos_hash matches current EOS at bit level"; both are wired to the
same byte-stream contract documented in
inject_superanc_check_points.py::eos_fnv1a_hex. Having one of each lets
the dev-time check run without a CatchTestRunner build and makes CI-side
pre-merge gating trivial.

Exits with status 0 if every fluid's hash matches (or no hash is stored);
status 1 if any fluid is stale.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from inject_superanc_check_points import eos_fnv1a_hex  # noqa: E402

FLUIDS = Path(__file__).resolve().parents[2] / 'dev' / 'fluids'


def main():
    stale = []
    checked = skipped = 0
    for path in sorted(FLUIDS.glob('*.json')):
        d = json.loads(path.read_text())
        if 'EOS' not in d or not d['EOS']:
            skipped += 1
            continue
        eos = d['EOS'][0]
        sa = eos.get('SUPERANCILLARY') or {}
        stored = sa.get('source_eos_hash')
        if stored is None:
            skipped += 1
            continue
        current = eos_fnv1a_hex(eos)
        if current != stored:
            stale.append((path.stem, stored, current))
        checked += 1

    print(f'Checked {checked} fluids with stored source_eos_hash; skipped {skipped}.')
    if stale:
        print(f'STALE: {len(stale)} fluid(s) whose EOS differs from the hash '
              'stamped at SA-generation time:')
        for name, stored, current in stale:
            print(f'  {name}: stored={stored}  current={current}')
        print()
        print('Regenerate the superancillary for these fluids against the '
              'current EOS by bumping the fastchebpure pin.')
        return 1
    print('All superancillary source_eos_hash stamps match the current EOS.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
