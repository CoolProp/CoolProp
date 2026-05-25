"""scikit-build-core dynamic-metadata version provider for CoolProp.

Why this exists
---------------
The Python package version is *dynamic* (``dynamic = ["version"]`` in
``pyproject.toml``).  scikit-build-core resolves it during the PEP 517
metadata hook (``prepare_metadata_for_build_wheel``), which runs **before
and independently of** the CMake build — so the version that names the
wheel and lands in the package metadata cannot be produced by CMake at
compile time.  This provider is the supported place for the build to
*call Python* to compute that version.

Resolution order
----------------
1. ``$COOLPROP_VERSION_OVERRIDE`` — if set and non-empty, used verbatim.
   CI sets this from ``dev/extract_version.py`` so every TestPyPI/PyPI
   upload gets a unique version (the previous mechanism wrote a ``.version``
   file that nothing read, so every build was ``7.2.1.dev0`` and TestPyPI
   rejected all but the first with ``400 File already exists``).
2. Otherwise the base version parsed from ``CMakeLists.txt``
   (``{MAJOR}.{MINOR}.{PATCH}{REVISION}``).  This keeps local / offline
   ``pip install .`` builds working with no env var and no tracked version
   file — identical to the historical default (e.g. ``7.2.1.dev0``).

CMakeLists.txt remains the single source of truth for the *base* version
(and the C++ library version); this provider only layers on the optional
CI-supplied unique suffix.

Wired up in ``pyproject.toml``::

    [tool.scikit-build.metadata.version]
    provider = "coolprop_version_provider"
    provider-path = "dev"
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Mapping, Optional

_ROOT = Path(__file__).resolve().parent.parent
_CMAKELISTS = _ROOT / "CMakeLists.txt"

# Mirrors the regex that pyproject.toml previously applied directly to
# CMakeLists.txt.  REVISION is optional/empty for tagged releases (giving a
# clean "7.2.1") and "dev" on the development line (giving "7.2.1.dev0").
_VERSION_RE = re.compile(
    r"set\(COOLPROP_VERSION_MAJOR\s+(?P<major>\d+)\).*?"
    r"set\(COOLPROP_VERSION_MINOR\s+(?P<minor>\d+)\).*?"
    r"set\(COOLPROP_VERSION_PATCH\s+(?P<patch>\d+)\).*?"
    r"set\(COOLPROP_VERSION_REVISION\s+(?P<dev>\w*)\)",
    re.DOTALL,
)


def _base_version_from_cmake() -> str:
    text = _CMAKELISTS.read_text(encoding="utf-8")
    match = _VERSION_RE.search(text)
    if match is None:
        msg = f"Could not parse COOLPROP_VERSION_* from {_CMAKELISTS}"
        raise RuntimeError(msg)
    return "{major}.{minor}.{patch}{dev}".format(**match.groupdict())


def dynamic_metadata(
    field: str,
    settings: Optional[Mapping[str, Any]] = None,
) -> str:
    if field != "version":
        msg = f"This provider only supplies 'version', not {field!r}"
        raise RuntimeError(msg)

    override = os.environ.get("COOLPROP_VERSION_OVERRIDE", "").strip()
    if override:
        return override

    return _base_version_from_cmake()
