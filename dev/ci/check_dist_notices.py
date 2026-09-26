#!/usr/bin/env python3
"""Assert every built Python distribution ships THIRD_PARTY_NOTICES.md.

scikit-build-core silently drops a ``wheel.license-files`` entry that matches
nothing, so a renamed/excluded notices file would otherwise ship wheels
without the third-party notices their bundled components require.

Usage: check_dist_notices.py DIST_DIR REFERENCE_NOTICES
Fails if DIST_DIR contains no wheels or no sdist, or if any wheel's
``*.dist-info/licenses/THIRD_PARTY_NOTICES.md`` or any sdist's root
``THIRD_PARTY_NOTICES.md`` is missing or differs from REFERENCE_NOTICES
(the repository copy), so an empty or stale file cannot pass.
"""

import sys
import tarfile
import zipfile
from pathlib import Path

NOTICES = "THIRD_PARTY_NOTICES.md"


def main(dist_dir: Path, reference: Path) -> int:
    expected = reference.read_bytes()
    wheels = sorted(dist_dir.rglob("*.whl"))
    sdists = sorted(dist_dir.rglob("*.tar.gz"))
    if not wheels or not sdists:
        print(f"ERROR: expected wheels and an sdist under {dist_dir}; found "
              f"{len(wheels)} wheel(s), {len(sdists)} sdist(s)", file=sys.stderr)
        return 1

    problems = []
    for whl in wheels:
        with zipfile.ZipFile(whl) as zf:
            hits = [n for n in zf.namelist() if n.endswith(f".dist-info/licenses/{NOTICES}")]
            if len(hits) != 1:
                problems.append(f"{whl.name}: {len(hits)} copies of {NOTICES} in dist-info/licenses")
            elif zf.read(hits[0]) != expected:
                problems.append(f"{whl.name}: {NOTICES} differs from {reference}")
    for sdist in sdists:
        with tarfile.open(sdist) as tf:
            # sdist layout is <name>-<version>/<file>
            hits = [m for m in tf.getmembers() if m.name.count("/") == 1 and m.name.endswith(f"/{NOTICES}")]
            if len(hits) != 1:
                problems.append(f"{sdist.name}: {len(hits)} copies of {NOTICES} at the root")
            elif tf.extractfile(hits[0]).read() != expected:
                problems.append(f"{sdist.name}: {NOTICES} differs from {reference}")

    for msg in problems:
        print(f"ERROR: {msg}", file=sys.stderr)
    print(f"Checked {len(wheels)} wheel(s) and {len(sdists)} sdist(s); {len(problems)} problem(s).")
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
