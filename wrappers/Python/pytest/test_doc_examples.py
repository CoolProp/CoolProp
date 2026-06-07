"""
Doc-example compat suite for the v8 (nanobind) CoolProp Python interface
(bd CoolProp-r9sq.13).

Rather than building the whole Sphinx site, this harvests the executable code
from the docs' ``.. ipython::`` directives and replays it under pytest -- a fast,
CI-cheap, comprehensive integration test of the high-level API as it is actually
used in the documentation. Any example that fails to run against the installed
CoolProp surfaces a compat gap (the original motivation: the array-``PropsSI`` /
2-arg-``PropsSI`` regressions would have shown up here as failures).

Bar: run-without-exception (the values are the C++ core's concern). One test per
doc page, replaying that page's ``.ipython`` blocks in order in a shared
namespace (matching the Sphinx directive's per-document state). Pages whose
examples need an absent optional dependency (matplotlib, pybtex, scipy) or
REFPROP are skipped, not failed.
"""
import os
import re
import glob

import pytest

# Locate the repo's Web/ docs relative to this file (dev-only test; the .rst
# sources are not shipped in the wheel).
_HERE = os.path.dirname(os.path.abspath(__file__))
_WEB = os.path.normpath(os.path.join(_HERE, "..", "..", "..", "Web"))

# Prompt lines inside an .ipython block: "In [1]: code" or "   ...: code".
_PROMPT = re.compile(r"^\s*(?:In \[\d+\]|\.\.\.)\s*:\s?(.*)$")



def _extract_ipython_code(rst_text):
    """Return the concatenated executable code from all non-verbatim
    ``.. ipython::`` blocks in a single .rst document, in order."""
    lines = rst_text.splitlines()
    code_chunks = []
    i, n = 0, len(lines)
    while i < n:
        if lines[i].lstrip().startswith(".. ipython::"):
            i += 1
            verbatim = False
            # directive option lines (":verbatim:", ":okexcept:", ...)
            while i < n and lines[i].strip().startswith(":"):
                if "verbatim" in lines[i]:
                    verbatim = True
                i += 1
            # the indented directive body (until a non-indented, non-blank line)
            body = []
            while i < n:
                ln = lines[i]
                if ln.strip() == "":
                    body.append("")
                    i += 1
                    continue
                if not ln.startswith((" ", "\t")):
                    break
                body.append(ln)
                i += 1
            if not verbatim:
                chunk = _parse_block(body)
                if chunk.strip():
                    code_chunks.append(chunk)
        else:
            i += 1
    return "\n".join(code_chunks)


def _parse_block(body):
    """Pull the Python out of one block's body: keep only prompt lines, strip
    the ``In [n]:``/``...:`` prefix (preserving the code's own indentation),
    and drop ``@verbatim`` sub-blocks."""
    out = []
    skip = False
    for ln in body:
        s = ln.strip()
        if s.startswith("@verbatim"):
            skip = True
            continue
        if s.startswith("@"):  # @suppress/@doctest/@savefig markers -> drop the marker
            continue
        m = _PROMPT.match(ln)
        if m is None:
            # blank line ends a verbatim run; non-prompt text is output/comment
            if s == "":
                skip = False
            continue
        if skip:
            continue
        code = m.group(1)
        # Drop IPython magics / shell escapes (e.g. ``%run``, ``%timeit``, ``!cmd``);
        # they are not Python and are not part of the API surface under test.
        if code.lstrip().startswith(("%", "!")):
            continue
        out.append(code)
    return "\n".join(out)


def _doc_pages():
    pages = []
    for path in sorted(glob.glob(os.path.join(_WEB, "**", "*.rst"), recursive=True)):
        try:
            text = open(path, encoding="utf-8").read()
        except (OSError, UnicodeDecodeError):
            continue
        if ".. ipython::" in text:
            pages.append(path)
    return pages


_PAGES = _doc_pages()

# Compat gaps the harvester surfaced; xfail (strict) until fixed -- when the gap
# is closed the page xpasses, turning the suite red so the xfail entry is removed.
_XFAIL = {
    "coolprop/LowLevelAPI.rst": "bd CoolProp-r9sq.15 (generate_update_pair binding broken)",
    "coolprop/HighLevelAPI.rst": "bd CoolProp-r9sq.16 (set_reference_state D-form -> std::bad_cast)",
}


def _params():
    params = []
    for p in _PAGES:
        rel = os.path.relpath(p, _WEB).replace(os.sep, "/")  # forward slashes on Windows too
        marks = [pytest.mark.xfail(reason=_XFAIL[rel], strict=True)] if rel in _XFAIL else []
        params.append(pytest.param(p, id=rel, marks=marks))
    return params


@pytest.mark.skipif(not _PAGES, reason="Web/ docs not found next to the repo")
@pytest.mark.parametrize("rst_path", _params())
def test_doc_page_examples(rst_path):
    code = _extract_ipython_code(open(rst_path, encoding="utf-8").read())
    if not code.strip():
        pytest.skip("no executable .ipython code")
    # Headless plotting; never pop a window.
    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        pass
    ns = {"__name__": "__doc_example__"}
    try:
        exec(compile(code, rst_path, "exec"), ns)
    except Exception as e:  # noqa: BLE001 -- this is the compat signal
        # Skip ONLY for a genuinely-optional dependency or REFPROP. A broken
        # CoolProp import (or a removed API imported by an example) is a real
        # compat failure and must NOT be masked as a skip.
        optional_roots = {"matplotlib", "pybtex", "scipy", "pandas"}
        missing_root = (getattr(e, "name", None) or "").split(".")[0]
        msg = str(e).lower()
        is_optional = (
            (isinstance(e, ModuleNotFoundError) and missing_root in optional_roots)
            or (isinstance(e, ImportError) and any(r in msg for r in optional_roots))
            or "refprop" in msg
        )
        if is_optional:
            pytest.skip("requires an absent optional dependency / REFPROP: %s" % str(e)[:80])
        raise
