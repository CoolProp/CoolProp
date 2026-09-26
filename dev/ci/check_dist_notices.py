#!/usr/bin/env python3
"""Assert every built Python distribution ships THIRD_PARTY_NOTICES.md.

scikit-build-core silently drops a ``wheel.license-files`` entry that matches
nothing, so a renamed/excluded notices file would otherwise ship wheels
without the third-party notices their bundled components require.

Usage: check_dist_notices.py DIST_DIR
Fails if DIST_DIR contains no wheels, if any wheel lacks
``*.dist-info/licenses/THIRD_PARTY_NOTICES.md``, or if any sdist lacks
``THIRD_PARTY_NOTICES.md`` at its root.
"""

import sys
import tarfile
import zipfile
from pathlib import Path

NOTICES = "THIRD_PARTY_NOTICES.md"


def main(dist_dir: Path) -> int:
    wheels = sorted(dist_dir.rglob("*.whl"))
    sdists = sorted(dist_dir.rglob("*.tar.gz"))
    if not wheels:
        print(f"ERROR: no wheels found under {dist_dir}", file=sys.stderr)
        return 1

    bad = []
    for whl in wheels:
        with zipfile.ZipFile(whl) as zf:
            names = zf.namelist()
        if not any(n.endswith(f".dist-info/licenses/{NOTICES}") for n in names):
            bad.append(whl)
    for sdist in sdists:
        with tarfile.open(sdist) as tf:
            # sdist layout is <name>-<version>/<file>
            names = tf.getnames()
        if not any(n.count("/") == 1 and n.endswith(f"/{NOTICES}") for n in names):
            bad.append(sdist)

    for path in bad:
        print(f"ERROR: {path.name} is missing {NOTICES}", file=sys.stderr)
    print(f"Checked {len(wheels)} wheel(s) and {len(sdists)} sdist(s); {len(bad)} missing {NOTICES}.")
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(Path(sys.argv[1])))
