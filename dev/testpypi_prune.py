"""Compute which TestPyPI (or PyPI) CoolProp releases to prune.

TestPyPI enforces a per-project storage quota, and every non-tag CI run
uploads the full wheel matrix (~43 files per ``.postNNN`` release).  To stay
under the ceiling we keep only the newest ``--keep`` *pre/post* releases and
delete the rest.

This script is **read-only**.  It queries the public JSON API, decides the
delete set, and emits an anchored alternation regex on stdout for
``pypi-cleanup -r`` to consume.  Final (non-dev/post/pre) releases are NEVER
candidates for deletion, so a real tagged release can't be pruned by accident.

Matching note: ``pypi-cleanup`` deletes a version ``k`` when ``regex.match(k)``
is truthy (see its source) -- ``re.match`` anchors the *start* only, so we
emit ``^(?:...)$`` (with a trailing ``$``) to avoid a prefix collision between
sibling timestamps.  The alternatives are the PEP 440-*normalized* version
strings, because ``pypi-cleanup`` matches against the normalized versions from
the ``/simple/`` index rather than the raw ``/pypi/.../json`` keys.

Usage::

    python dev/testpypi_prune.py --keep 10            # regex -> stdout, plan -> stderr
    python dev/testpypi_prune.py --keep 10 --pypi     # target real PyPI instead
"""
import sys
import argparse

import requests
from packaging import version


def fetch_release_keys(package, pypi=False):
    """Return the raw release-version strings from the JSON API."""
    host = "https://pypi.org" if pypi else "https://test.pypi.org"
    response = requests.get(f"{host}/pypi/{package}/json", timeout=30)
    response.raise_for_status()
    return list(response.json()["releases"].keys())


def select_delete_set(release_keys, keep):
    """Split release keys into (keep_keys, delete_keys).

    Only dev/pre/post releases are candidates.  Candidates are ordered by
    parsed version (for ``.postNNN`` timestamps that is chronological), the
    newest ``keep`` are retained, and the remainder are returned as the
    delete set.  Final releases are returned in neither list -- they are
    untouchable.
    """
    if keep < 1:
        raise ValueError("--keep must be >= 1 (refusing to delete everything)")

    # Map each release to its PEP 440-normalized string; skip anything
    # unparseable rather than risk a mis-targeted delete.  We key the regex off
    # the *normalized* form (str(version.parse(...))) because pypi-cleanup
    # matches against the normalized versions from the /simple/ index, not the
    # raw /pypi/.../json keys -- a non-normalized key (e.g. "6.4.1-1") would
    # otherwise silently fail to match and the prune would under-delete.
    candidates = []
    for key in release_keys:
        try:
            v = version.parse(key)
        except version.InvalidVersion:
            continue
        if v.is_prerelease or v.is_postrelease:
            candidates.append((v, str(v)))

    candidates.sort(key=lambda pair: pair[0])
    keep_pairs = candidates[-keep:] if keep < len(candidates) else candidates
    delete_pairs = candidates[:-keep] if keep < len(candidates) else []
    return [s for _, s in keep_pairs], [s for _, s in delete_pairs]


def build_regex(delete_keys):
    """Anchored alternation matching exactly the delete-set keys."""
    import re
    alts = "|".join(re.escape(k) for k in delete_keys)
    return f"^(?:{alts})$"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute the TestPyPI/PyPI prune regex (keep newest N pre/post releases)")
    parser.add_argument("--package", default="CoolProp", help="package name")
    parser.add_argument("--keep", type=int, default=10,
                        help="number of newest pre/post releases to retain")
    parser.add_argument("--pypi", action="store_true",
                        help="target real PyPI instead of TestPyPI")
    args = parser.parse_args()

    keys = fetch_release_keys(args.package, pypi=args.pypi)
    keep_keys, delete_keys = select_delete_set(keys, args.keep)

    remote = "PyPI" if args.pypi else "TestPyPI"
    print(f"{remote} {args.package}: {len(keys)} total releases; "
          f"keeping {len(keep_keys)} newest pre/post, deleting {len(delete_keys)}.",
          file=sys.stderr)
    for k in delete_keys:
        print(f"  DELETE {k}", file=sys.stderr)
    for k in keep_keys:
        print(f"  keep   {k}", file=sys.stderr)

    # stdout = machine-readable regex (empty when nothing to prune).
    if delete_keys:
        print(build_regex(delete_keys))
