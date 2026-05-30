# --- hand-written overloads (do not regenerate; spliced by postprocess.py) ---
#
# The high-level property functions are typed by hand because stubgen cannot
# express their scalar/array polymorphism (overloads are a stub-only construct,
# and Cython annotation_typing would coerce the args at runtime).
#
# NOTE on overload overlap: PropsSI/HAPropsSI/HAProps return ndarray iff ANY
# numeric input is array-like — which cannot be expressed without the scalar
# overload (-> float) and the array overload (whose `_Num` includes float, to
# admit mixed scalar+array calls -> ndarray) overlapping on an all-float call.
# Ordering (scalar first) makes the all-scalar case resolve to float under both
# pyright and mypy; the array overload carries `# type: ignore[overload-overlap]`
# so mypy does not flag the unavoidable overlap.  Per a .pyi stub, NO
# implementation line follows the overloads (that is itself a mypy error).

from typing import Sequence

_Scalar = float
_Array = Sequence[float] | ndarray
_Num = _Scalar | _Array


# PropsSI: two-argument trivial-property lookup -> float; six-argument form is
# scalar->float or array->ndarray.  (in1/Output and in3/in5 are required for the
# six-arg form; the two-arg form is _Props1SI(Output, FluidName).)
@overload
def PropsSI(Output: str, FluidName: str) -> float: ...
@overload
def PropsSI(Output: str, Name1: str, Prop1: _Scalar, Name2: str, Prop2: _Scalar, FluidName: str) -> float: ...  # type: ignore[overload-overlap]
@overload
def PropsSI(Output: str, Name1: str, Prop1: _Num, Name2: str, Prop2: _Num, FluidName: str) -> ndarray: ...


# PhaseSI does NOT vectorize (per its docstring) — scalar inputs, returns a
# single phase string.
def PhaseSI(Name1: str, Prop1: _Scalar, Name2: str, Prop2: _Scalar, FluidName: str) -> str: ...


# Props is the deprecated scalar-only predecessor of PropsSI (no array path).
def Props(Output: str, Name1: str, Prop1: _Scalar, Name2: str, Prop2: _Scalar, FluidName: str) -> float: ...


# HAPropsSI / HAProps DO vectorize (scalar->float, array->ndarray).
@overload
def HAPropsSI(Output: str, Name1: str, Value1: _Scalar, Name2: str, Value2: _Scalar, Name3: str, Value3: _Scalar) -> float: ...  # type: ignore[overload-overlap]
@overload
def HAPropsSI(Output: str, Name1: str, Value1: _Num, Name2: str, Value2: _Num, Name3: str, Value3: _Num) -> ndarray: ...


@overload
def HAProps(Output: str, Name1: str, Value1: _Scalar, Name2: str, Value2: _Scalar, Name3: str, Value3: _Scalar) -> float: ...  # type: ignore[overload-overlap]
@overload
def HAProps(Output: str, Name1: str, Value1: _Num, Name2: str, Value2: _Num, Name3: str, Value3: _Num) -> ndarray: ...
