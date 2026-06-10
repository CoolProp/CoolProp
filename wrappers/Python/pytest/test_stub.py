"""
Confirm the v8 nanobind type stub is valid and matches the runtime (bd
CoolProp-r9sq.10).

Note on tooling: nanobind functions are opaque to ``inspect.signature`` (they
present as ``(*args, **kwargs)``), so ``mypy.stubtest``'s signature comparison
does NOT apply to a nanobind stub. Instead we confirm the properties that
actually matter:
  1. ``CoolProp.pyi`` + ``py.typed`` ship and the stub parses;
  2. every public runtime symbol is present in the stub (symbol parity);
  3. a type checker resolves the key polymorphic overloads -- ``PropsSI`` scalar
     -> ``float`` and array -> ``ndarray`` -- which is the point of the stub.
"""
import ast
import os

import pytest

import CoolProp.CoolProp as CP

_PKG = os.path.dirname(CP.__file__)
_PYI = os.path.join(_PKG, "CoolProp.pyi")


def _stub_toplevel_names():
    names = set()
    with open(_PYI, encoding="utf-8") as f:
        tree = ast.parse(f.read())
    for n in tree.body:
        if isinstance(n, (ast.FunctionDef, ast.ClassDef)):
            names.add(n.name)
        elif isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name):
            names.add(n.target.id)
        elif isinstance(n, ast.Assign):
            names.update(t.id for t in n.targets if isinstance(t, ast.Name))
    return names


def test_stub_and_marker_ship_and_parse():
    assert os.path.exists(_PYI), "CoolProp.pyi not shipped in the wheel"
    assert os.path.exists(os.path.join(_PKG, "py.typed")), "py.typed marker not shipped"
    with open(_PYI, encoding="utf-8") as f:
        ast.parse(f.read())  # valid Python (raises on failure)


def test_stub_symbol_parity_with_runtime():
    names = _stub_toplevel_names()
    runtime = {n for n in dir(CP) if not n.startswith("__")}
    missing = sorted(runtime - names)
    assert not missing, "runtime symbols missing from the stub: %s" % missing


def test_stub_overloads_typecheck():
    """mypy resolves PropsSI scalar->float and array->ndarray against the stub."""
    mypy_api = pytest.importorskip("mypy.api")
    snippet = (
        "from CoolProp.CoolProp import PropsSI\n"
        "reveal_type(PropsSI('Tcrit', 'Water'))\n"
        "reveal_type(PropsSI('D', 'T', 300.0, 'P', 101325.0, 'Water'))\n"
        "reveal_type(PropsSI('D', 'T', [300.0, 310.0], 'P', 101325.0, 'Water'))\n"
    )
    out, _err, _code = mypy_api.run(["-c", snippet])
    revealed = [ln for ln in out.splitlines() if "Revealed type" in ln]
    assert len(revealed) == 3, out
    assert "float" in revealed[0]                       # 2-arg trivial
    assert "float" in revealed[1]                       # scalar inputs
    assert "ndarray" in revealed[2]                     # array input
