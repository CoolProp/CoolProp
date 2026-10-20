"""
Numerics of property diagrams: isolines, the dome and the limits of a chart
"""

from collections.abc import Sequence
import enum
from typing import overload

import CoolProp.CoolProp


class Scale(enum.Enum):
    Lin = 0

    Log = 1

class TPLimits(enum.Enum):
    NONE = 0

    DEF = 1

    ACHP = 2

    ORC = 3

class Range:
    def __init__(self, min: float, max: float) -> None: ...

    @property
    def min(self) -> float: ...

    @min.setter
    def min(self, arg: float, /) -> None: ...

    @property
    def max(self) -> float: ...

    @max.setter
    def max(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

class CalcSettings:
    """
    How the isolines are computed: points per line, the solver (auto, tracer or flash), the margin to the critical point and the spacing of the quality lines (critical or linear).
    """

    def __init__(self) -> None: ...

    @property
    def points(self) -> int: ...

    @points.setter
    def points(self, arg: int, /) -> None: ...

    @property
    def saturation_points(self) -> int: ...

    @saturation_points.setter
    def saturation_points(self, arg: int, /) -> None: ...

    @property
    def solver(self) -> str: ...

    @solver.setter
    def solver(self, arg: str, /) -> None: ...

    @property
    def critical_margin(self) -> float: ...

    @critical_margin.setter
    def critical_margin(self, arg: float, /) -> None: ...

    @property
    def saturation_spacing(self) -> str: ...

    @saturation_spacing.setter
    def saturation_spacing(self, arg: str, /) -> None: ...

class Isoline:
    """
    One computed line: x and y in SI units with NaN for a gap, and how it was computed.
    """

    @property
    def x(self) -> list[float]: ...

    @property
    def y(self) -> list[float]: ...

    @property
    def value(self) -> float: ...

    @property
    def n_failed(self) -> int: ...

    @property
    def n_fallback(self) -> int: ...

    @property
    def n_inserted(self) -> int: ...

    @property
    def n_evaluated(self) -> int: ...

    @property
    def n_outside(self) -> int: ...

    @property
    def method(self) -> str: ...

    @property
    def reason(self) -> str: ...

    def size(self) -> int: ...

    def n_valid(self) -> int: ...

    def n_gaps(self) -> int: ...

    def __len__(self) -> int: ...

class Axis:
    @property
    def scale(self) -> Scale: ...

    @property
    def min(self) -> float: ...

    @property
    def max(self) -> float: ...

class PropertyPlot:
    """
    A property diagram of one fluid: its T-p box, its axes and the isolines it can compute.
    The property keys are CoolProp.parameters (iT, iP, iHmass, iSmass, iDmass, iUmass, iQ).
    Not safe to use from several threads at once (calc_isolines updates the plot's state and
    releases the GIL); give each thread its own PropertyPlot.
    """

    @overload
    def __init__(self, fluid_name: str, ykey: CoolProp.CoolProp.parameters, xkey: CoolProp.CoolProp.parameters, tp_limits: TPLimits = TPLimits.DEF) -> None: ...

    @overload
    def __init__(self, fluid_name: str, ykey: CoolProp.CoolProp.parameters, xkey: CoolProp.CoolProp.parameters, T_min_K: float, T_max_K: float, p_min_Pa: float, p_max_Pa: float) -> None: ...

    @property
    def xaxis(self) -> Axis: ...

    @property
    def yaxis(self) -> Axis: ...

    def isoline_range(self, key: CoolProp.CoolProp.parameters) -> Range: ...

    @overload
    def calc_isolines(self, key: CoolProp.CoolProp.parameters, values: Sequence[float], points: int = 250) -> list[Isoline]: ...

    @overload
    def calc_isolines(self, key: CoolProp.CoolProp.parameters, values: Sequence[float], settings: CalcSettings) -> list[Isoline]: ...

    def supported_isoline_keys(self) -> list[CoolProp.CoolProp.parameters]: ...

    def sweep_key(self, key: CoolProp.CoolProp.parameters) -> CoolProp.CoolProp.parameters: ...

    def value_at(self, key: CoolProp.CoolProp.parameters, xvalue: float, yvalue: float, phase: CoolProp.CoolProp.phases = CoolProp.CoolProp.phases.iphase_not_imposed) -> float: ...

    def tp_box(self) -> list[float]: ...

    def is_mixture(self) -> bool: ...

@overload
def generate_values_in_range(scale: Scale, range: Range, count: int) -> list[float]: ...

@overload
def generate_values_in_range(type: CoolProp.CoolProp.parameters, range: Range, count: int) -> list[float]: ...
