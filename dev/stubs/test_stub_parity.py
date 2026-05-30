"""Keep CoolProp.pyi in sync with the compiled CoolProp.CoolProp module.

This is the cheap, dependency-light half of the stub sync strategy (it needs
only a built CoolProp, not the Cython-3.x generator).  It asserts that the stub
and the runtime agree on the public surface of the module:

  * the AbstractState class (the bulk of the typed API), and
  * the module-level functions (PropsSI, PhaseSI, ...).

Drift in either direction is caught: a stub that names a method the runtime no
longer has, or a runtime method/function never regenerated into the stub.

The exhaustive signature/type check is the CI "drift gate" that regenerates the
stub and diffs it; see dev/stubs/gen_stubs.sh.

This file lives in dev/stubs/ (NOT inside the CoolProp package) on purpose: a
test placed under wrappers/Python/CoolProp/tests/ is imported by pytest as
``CoolProp.tests.<name>``, which resolves against the *installed* CoolProp
package, not the source tree — so it would fail to import.  Here it imports
standalone and picks up CoolProp from the environment.
"""
from __future__ import annotations

import ast
import pathlib

import pytest

# AbstractState et al. live in the CoolProp.CoolProp extension module
# (AbstractState.pyx is ``include``d into CoolProp.pyx); there is no separate
# CoolProp.AbstractState module.
cpcp = pytest.importorskip("CoolProp.CoolProp")

# dev/stubs/ -> repo root -> wrappers/Python/CoolProp/CoolProp.pyi
STUB_PATH = (
    pathlib.Path(__file__).resolve().parents[2]
    / "wrappers" / "Python" / "CoolProp" / "CoolProp.pyi"
)


def _stub_tree() -> ast.Module:
    return ast.parse(STUB_PATH.read_text(encoding="utf-8"))


def _stub_class_methods(tree: ast.Module, cls_name: str) -> set[str]:
    cls = next(
        (n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == cls_name),
        None,
    )
    assert cls is not None, f"{cls_name} class missing from CoolProp.pyi"
    return {
        n.name
        for n in cls.body
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
        and not n.name.startswith("_")
    }


def _stub_module_funcs(tree: ast.Module) -> set[str]:
    return {
        n.name
        for n in tree.body
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
        and not n.name.startswith("_")
    }


def _runtime_methods(obj) -> set[str]:
    return {
        name
        for name in dir(obj)
        if not name.startswith("_") and callable(getattr(obj, name, None))
    }


def test_stub_file_exists_and_parses():
    assert STUB_PATH.is_file(), f"missing stub: {STUB_PATH}"
    ast.parse(STUB_PATH.read_text(encoding="utf-8"))  # raises on garbage


def test_abstractstate_methods_match_runtime():
    stub = _stub_class_methods(_stub_tree(), "AbstractState")
    runtime = _runtime_methods(cpcp.AbstractState)
    stale = stub - runtime
    missing = runtime - stub
    assert not stale, (
        f"CoolProp.pyi:AbstractState declares methods absent at runtime: "
        f"{sorted(stale)} — stale stub; run dev/stubs/gen_stubs.sh"
    )
    assert not missing, (
        f"runtime AbstractState has public methods missing from the stub: "
        f"{sorted(missing)} — regenerate with dev/stubs/gen_stubs.sh"
    )


def test_module_functions_present_at_runtime():
    stub_funcs = _stub_module_funcs(_stub_tree())
    runtime = {
        n for n in dir(cpcp)
        if not n.startswith("_") and callable(getattr(cpcp, n, None))
    }
    stale = stub_funcs - runtime
    assert not stale, (
        f"CoolProp.pyi declares module functions absent at runtime: "
        f"{sorted(stale)} — stale stub; run dev/stubs/gen_stubs.sh"
    )
