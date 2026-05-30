"""Static type-correctness smoke for the shipped CoolProp.pyi stub.

This file is NOT executed.  A static type checker (pyright or mypy) checks it to
prove the stub expresses the intended types.  Run, e.g.:

    pyright dev/stubs/typecheck_smoke.py
    # or
    mypy dev/stubs/typecheck_smoke.py

`assert_type` is a no-op at runtime but makes the type checker fail if the
inferred type differs, so a regression in CoolProp.pyi breaks CI here.
"""
from __future__ import annotations

from typing import assert_type

from numpy import ndarray

# AbstractState and the module functions all live in CoolProp.CoolProp.
from CoolProp.CoolProp import AbstractState, PropsSI

# --- AbstractState: scalar property getters are float -------------------------
state = AbstractState("HEOS", "Water")
assert_type(state.T(), float)
assert_type(state.p(), float)
assert_type(state.rhomolar(), float)
assert_type(state.keyed_output(0), float)
assert_type(state.first_partial_deriv(0, 0, 0), float)

# --- PropsSI overload resolution: trivial 2-arg, scalar 6-arg, array 6-arg ----
assert_type(PropsSI("Tcrit", "Water"), float)
assert_type(PropsSI("T", "P", 101325.0, "Q", 0.0, "Water"), float)
assert_type(PropsSI("T", "P", [101325.0, 202650.0], "Q", 0.0, "Water"), ndarray)
