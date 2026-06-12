#!/usr/bin/env python3
"""Assert a ``v*`` tag matches the CMakeLists version and is publishable to PyPI.

On a tagged release the wheels publish to **real PyPI** with a version derived
from ``CMakeLists.txt`` (via ``dev/extract_version.py``), NOT from the tag name.
So nothing currently couples the tag to the published version: tagging e.g.
``v8.0.0b1`` while ``CMakeLists.txt`` still says ``7.2.1dev`` would silently
publish ``7.2.1.dev0`` to real PyPI (the #2988 class of bug -- CoolProp-1tbe.3).

This gate, run on the tag path before the ``--pypi`` version is computed,
enforces two things:

  1. the tag's PEP 440 version equals the CMakeLists version (so ``v8.0.0b1``
     requires ``COOLPROP_VERSION_* = 8.0.0b1``); and
  2. that version is NOT a developmental (``.dev``) release -- a ``dev`` revision
     must never reach real PyPI.  Final, beta (``b1``) and rc (``rc1``) markers
     are all allowed, so a ``v8.0.0b1`` beta tag passes.

Usage:
    check_tag_version.py <git-ref-or-tag>     # e.g. refs/tags/v8.0.0b1

Exit 0 if the tag is consistent and publishable; non-zero with a GitHub-Actions
``::error::`` annotation otherwise.
"""
import re
import sys
from pathlib import Path

from packaging import version

CMAKELISTS = Path(__file__).resolve().parents[2] / "CMakeLists.txt"


def cmake_version() -> str:
    """Build the version string from CMakeLists.txt the same way the build does.

    Mirrors dev/extract_version.py / dev/coolprop_version_provider.py: REVISION
    may be unquoted (``dev`` / ``b1``) or a quoted string (``""`` on a final
    release); both are accepted, and it is concatenated directly onto the patch
    (so ``8.0.0`` + ``b1`` -> ``8.0.0b1``, which packaging normalizes to a beta).
    """
    text = CMAKELISTS.read_text(encoding="utf-8")

    def grab(key: str, pat: str) -> str:
        m = re.search(rf'set\(COOLPROP_VERSION_{key}\s+"?{pat}"?\)', text)
        if m is None:
            raise RuntimeError(f"could not parse COOLPROP_VERSION_{key} from {CMAKELISTS}")
        return m.group(1)

    major = grab("MAJOR", r"(\d+)")
    minor = grab("MINOR", r"(\d+)")
    patch = grab("PATCH", r"(\d+)")
    rev = grab("REVISION", r"(\w*)")
    return f"{major}.{minor}.{patch}{rev}"


def main(ref: str) -> int:
    tag = ref.rsplit("/", 1)[-1]  # refs/tags/v8.0.0b1 -> v8.0.0b1
    if tag[:1].lower() == "v":
        tag = tag[1:]

    try:
        tag_v = version.Version(tag)
    except version.InvalidVersion as e:
        print(f"::error::tag '{tag}' is not a PEP 440 version: {e}", file=sys.stderr)
        return 1

    cmake_str = cmake_version()
    try:
        cmake_v = version.Version(cmake_str)
    except version.InvalidVersion as e:
        print(f"::error::CMakeLists version '{cmake_str}' is not a PEP 440 version: {e}", file=sys.stderr)
        return 1

    if tag_v != cmake_v:
        print(
            f"::error::tag version {tag_v} does not match CMakeLists.txt version {cmake_v}. "
            f"Bump COOLPROP_VERSION_MAJOR/MINOR/PATCH/REVISION to match the tag before tagging "
            f"(the published PyPI version comes from CMakeLists.txt, not the tag).",
            file=sys.stderr,
        )
        return 1

    if cmake_v.is_devrelease:
        print(
            f"::error::refusing to publish a developmental version ({cmake_v}) to real PyPI. "
            f"Set COOLPROP_VERSION_REVISION to a release marker (empty), or a beta/rc (e.g. b1, rc1), not 'dev'.",
            file=sys.stderr,
        )
        return 1

    kind = "pre-release" if cmake_v.is_prerelease else "final release"
    print(f"OK: tag {tag_v} matches CMakeLists.txt {cmake_v} ({kind}); publishable to PyPI.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: check_tag_version.py <git-ref-or-tag>", file=sys.stderr)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
