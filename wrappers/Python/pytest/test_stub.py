"""
Confirm the v8 nanobind type stub is valid and matches the runtime (bd
CoolProp-r9sq.10).

Note on tooling: nanobind functions are opaque to ``inspect.signature`` (they
present as ``(*args, **kwargs)``), so ``mypy.stubtest``'s signature comparison
does NOT apply to a nanobind stub. Instead we confirm the properties that
actually matter:
  1. the stubs + ``py.typed`` ship and the stubs parse;
  2. every public runtime symbol is present in the stub (symbol parity);
  3. a type checker resolves the key polymorphic overloads -- ``PropsSI`` scalar
     -> ``float`` and array -> ``ndarray`` -- which is the point of the stub;
  4. a type checker finds the ``Plot`` submodule where the runtime has it,
     ``CoolProp.CoolProp.Plot``.

The module has a submodule, so its stubs are a stub package next to the
extension module: ``CoolProp/CoolProp/__init__.pyi`` and
``CoolProp/CoolProp/Plot.pyi`` (see ``wrappers/Python/CMakeLists.txt``).
"""
import ast
import os

import pytest

import CoolProp.CoolProp as CP

_PKG = os.path.dirname(CP.__file__)
_STUB_PACKAGE = os.path.join(_PKG, "CoolProp")              # stubs of CoolProp.CoolProp
_PYI = os.path.join(_STUB_PACKAGE, "__init__.pyi")
_PLOT_PYI = os.path.join(_STUB_PACKAGE, "Plot.pyi")         # the CoolProp.CoolProp.Plot submodule


def _stub_toplevel_names(path=_PYI):
    names = set()
    with open(path, encoding="utf-8") as f:
        tree = ast.parse(f.read())
    for n in tree.body:
        if isinstance(n, (ast.FunctionDef, ast.ClassDef)):
            names.add(n.name)
        elif isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name):
            names.add(n.target.id)
        elif isinstance(n, ast.Assign):
            names.update(t.id for t in n.targets if isinstance(t, ast.Name))
        elif isinstance(n, ast.ImportFrom) and n.level == 1:
            # A submodule re-exported as ``from . import Plot as Plot`` is a
            # symbol of the module too; stubgen writes it that way.
            names.update(a.asname or a.name for a in n.names)
    return names


def test_stub_and_marker_ship_and_parse():
    assert os.path.exists(_PYI), "CoolProp/CoolProp/__init__.pyi not shipped in the wheel"
    assert os.path.exists(os.path.join(_PKG, "py.typed")), "py.typed marker not shipped"
    with open(_PYI, encoding="utf-8") as f:
        ast.parse(f.read())  # valid Python (raises on failure)


def test_stubs_form_a_stub_package():
    """A stub module CoolProp/CoolProp.pyi next to the stub package would hide it
    from a type checker, and CoolProp/Plot.pyi would describe CoolProp.Plot, which
    does not exist at runtime."""
    for stale in ("CoolProp.pyi", "Plot.pyi"):
        assert not os.path.exists(os.path.join(_PKG, stale)), "stale stub %s shipped" % stale
    # The stub folder holds no Python file, so the extension module stays what
    # "import CoolProp.CoolProp" loads.
    assert not [f for f in os.listdir(_STUB_PACKAGE) if f.endswith(".py")]
    assert CP.__file__.endswith((".so", ".pyd")), CP.__file__


def test_stub_symbol_parity_with_runtime():
    names = _stub_toplevel_names()
    runtime = {n for n in dir(CP) if not n.startswith("__")}
    missing = sorted(runtime - names)
    assert not missing, "runtime symbols missing from the stub: %s" % missing


def test_plot_stub_ships_and_matches_runtime():
    """The submodule stub is a committed copy, so check it against the runtime by hand."""
    assert os.path.exists(_PLOT_PYI), "Plot.pyi not shipped in the wheel"
    names = _stub_toplevel_names(_PLOT_PYI)
    runtime = {n for n in dir(CP.Plot) if not n.startswith("_")}
    missing = sorted(runtime - names)
    assert not missing, "runtime symbols of CoolProp.CoolProp.Plot missing from Plot.pyi: %s" % missing
    stale = sorted(names - runtime)
    assert not stale, "Plot.pyi names symbols the runtime does not have: %s" % stale


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


def test_plot_stub_is_found_at_the_runtime_path():
    """mypy reads the Plot stub for CoolProp.CoolProp.Plot and knows no CoolProp.Plot."""
    mypy_api = pytest.importorskip("mypy.api")
    snippet = (
        "from CoolProp.CoolProp import iP, iHmass, iT\n"
        "from CoolProp.CoolProp.Plot import PropertyPlot\n"
        "import CoolProp.CoolProp as CP\n"
        "plot = PropertyPlot('HEOS::R134a', iP, iHmass)\n"
        "reveal_type(plot.calc_isolines(iT, [300.0], 50))\n"
        "reveal_type(CP.Plot.CalcSettings().points)\n"
    )
    out, _err, _code = mypy_api.run(["--no-incremental", "-c", snippet])
    revealed = [ln for ln in out.splitlines() if "Revealed type" in ln]
    assert len(revealed) == 2, out
    assert "list[CoolProp.CoolProp.Plot.Isoline]" in revealed[0], out
    assert revealed[1].endswith(('"int"', '"builtins.int"')), out   # mypy 2.x and 1.x spellings
    assert "error" not in out, out
    out, _err, code = mypy_api.run(["--no-incremental", "-c", "import CoolProp.Plot\n"])
    assert code != 0 and "CoolProp.Plot" in out, out
