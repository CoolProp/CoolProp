#!/usr/bin/env python3
"""Post-process a stubgen-pyx-generated .pyi into valid, self-contained typing.

stubgen-pyx faithfully transcribes the Cython declarations of CoolProp.pyx
(which ``include``s HumidAirProp.pyx and AbstractState.pyx, so the whole
CoolProp.CoolProp module lands in one stub).  Several of the emitted types are
Cython/C++ constructs that are NOT valid Python typing and must be rewritten:

  * Non-importable imports: ``import cython``; ``from libcpp import bool``;
    ``from libcpp.string import string``; ``from libcpp.vector import vector``;
    ``from . import constants_header`` / ``superancillary`` / ``from .typedefs
    import *`` (Cython .pxd modules with no shipped Python/.pyi counterpart).
  * ``CoolPropDbl`` (a ``ctypedef double``) leaks unresolved; the derivative
    methods that return it come out as ``-> ...``.
  * ``string`` / ``string_or_size_t`` / ``vector`` C++ types in annotations.
  * the enum types (parameters/phases/input_pairs/configuration_keys) are plain
    ``int`` at the Python boundary (callers pass module-level int constants such
    as PT_INPUTS / iT / iphase_liquid).
  * ``supanc.SuperAncillaryHolder`` — an unshipped helper type -> ``object``.
  * Cython-only constructs meaningless in a stub: ``@cython.*`` decorators,
    ``__cinit__`` (-> ``__init__``), ``__dealloc__`` (dropped).

Substitutions are applied ONLY in annotation position (after ``:`` in a
parameter / variable, or after ``->``) so an identifier that merely *happens* to
be named ``string`` (e.g. ``def set_departure_functions(string): ...``) is left
untouched.

The transform is deterministic, which is what lets the CI drift gate compare a
fresh regeneration against the committed .pyi via ``git diff --exit-code``.

Usage:  postprocess.py <input.pyi> <output.pyi>
"""
from __future__ import annotations

import re
import sys

HEADER = '''\
# AUTO-GENERATED — DO NOT EDIT BY HAND.
# Regenerate with: dev/stubs/gen_stubs.sh
# Source: wrappers/Python/CoolProp/CoolProp.pyx (which includes HumidAirProp.pyx
# and AbstractState.pyx) via stubgen-pyx + dev/stubs/postprocess.py.
#
# Enum-like arguments (parameters / phases / input_pairs / configuration_keys)
# are plain ``int`` at the Python boundary: pass the module-level constants from
# CoolProp.constants (e.g. PT_INPUTS, iT, iphase_liquid), which are integers.
from __future__ import annotations

from typing import Any, overload

from numpy import ndarray
'''

# Map a leaked type *name* (as it appears in an annotation) to valid typing.
# Applied only in annotation position (see _sub_annotations).
_TYPE_MAP = {
    "string_or_size_t": "str | int",
    "CoolPropDbl": "float",
    "string": "str",
    "vector": "list[float]",
    "constants_header.parameters": "int",
    "constants_header.phases": "int",
    "constants_header.input_pairs": "int",
    "constants_header.configuration_keys": "int",
    "parameters": "int",
    "phases": "int",
    "input_pairs": "int",
    "configuration_keys": "int",
    "supanc.SuperAncillaryHolder": "object",
    "bool": "bool",  # identity, but normalises any libcpp-bool reference
}

# Longest keys first so qualified names win over their bare suffixes.
_TYPE_KEYS = sorted(_TYPE_MAP, key=len, reverse=True)
_TYPE_ALT = "|".join(re.escape(k) for k in _TYPE_KEYS)

# An annotation token: ``: <type>`` (param/var) or ``-> <type>`` (return).  The
# type may carry a ``[...]`` subscript and a ``= default`` / ``,`` / ``)`` tail
# which we leave alone.  We only rewrite the leading type *name*.
_ANNOT = re.compile(r"(:\s*|->\s*)(" + _TYPE_ALT + r")\b")
_ELLIPSIS_RET = re.compile(r"->\s*\.\.\.(\s*:)")


def _sub_annotations(text: str) -> str:
    """Rewrite leaked types in annotation position, on code lines only.

    Substitution is applied per line and skips docstring bodies, so prose like
    ``:returns: phases of the mixture`` or ``mode: string id`` inside a
    docstring is never mangled — only real signatures / attribute annotations
    are touched.
    """
    out: list[str] = []
    in_doc = False
    delim = ""
    for line in text.splitlines():
        if in_doc:
            out.append(line)
            if delim in line:  # docstring closes
                in_doc, delim = False, ""
            continue
        stripped = line.lstrip()
        opener = '"""' if stripped.startswith('"""') else ("'''" if stripped.startswith("'''") else "")
        if opener:
            # Docstring line: leave verbatim; enter doc-mode unless it also
            # closes on the same line (one-line docstring).
            if opener not in stripped[len(opener):]:
                in_doc, delim = True, opener
            out.append(line)
            continue
        line = _ELLIPSIS_RET.sub(r"-> float\1", line)  # unresolved CoolPropDbl returns
        out.append(_ANNOT.sub(lambda m: m.group(1) + _TYPE_MAP[m.group(2)], line))
    return "\n".join(out)


def _strip_header(text: str) -> str:
    """Drop everything before the first top-level ``class``/``def``."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if re.match(r"^(class |def )", line):
            return "\n".join(lines[i:])
    return text


def _drop_methods(text: str, names: tuple[str, ...]) -> str:
    """Remove a def (and its indented body) for the given method names."""
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        m = re.match(r"^(\s*)def (\w+)\(", lines[i])
        if m and m.group(2) in names:
            indent = m.group(1)
            i += 1
            while i < len(lines) and (
                lines[i].strip() == ""
                or len(lines[i]) - len(lines[i].lstrip()) > len(indent)
            ):
                i += 1
            continue
        out.append(lines[i])
        i += 1
    return "\n".join(out)


def _drop_toplevel_defs(text: str, names: set[str]) -> str:
    """Drop module-level (unindented) ``def NAME`` blocks for the given names."""
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        m = re.match(r"^def (\w+)\(", lines[i])
        if m and m.group(1) in names:
            i += 1
            while i < len(lines) and (lines[i].strip() == "" or lines[i][:1] in " \t"):
                i += 1
            continue
        out.append(lines[i])
        i += 1
    return "\n".join(out)


def postprocess(text: str, overloads: str | None = None) -> str:
    text = _strip_header(text)
    text = _drop_methods(text, ("__dealloc__",))
    text = re.sub(r"\bdef __cinit__\b", "def __init__", text)
    # drop Cython-only decorator lines and any leftover alias assignments
    text = "\n".join(
        l
        for l in text.splitlines()
        if not re.match(r"^\s*@cython\.", l)
        and not re.match(r"^\w+\s*:\s*_TypeAlias", l)
    )
    text = _sub_annotations(text)
    body = text.strip("\n")
    if overloads:
        # Replace the generated (untyped) defs that the fragment redefines.
        replaced = set(re.findall(r"^def (\w+)\(", overloads, re.M))
        body = _drop_toplevel_defs(body, replaced).strip("\n")
        body += "\n\n\n" + overloads.strip("\n")
    return HEADER + "\n\n" + body + "\n"


def main(argv: list[str]) -> int:
    if len(argv) not in (3, 4):
        print(__doc__)
        print(f"error: expected 2-3 args, got {len(argv) - 1}", file=sys.stderr)
        return 2
    with open(argv[1], encoding="utf-8") as fh:
        raw = fh.read()
    overloads = None
    if len(argv) == 4:
        with open(argv[3], encoding="utf-8") as fh:
            overloads = fh.read()
    with open(argv[2], "w", encoding="utf-8") as fh:
        fh.write(postprocess(raw, overloads))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
