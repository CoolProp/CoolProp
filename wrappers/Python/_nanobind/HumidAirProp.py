"""
``CoolProp.HumidAirProp`` -- humid-air properties (v8 nanobind package).

The SI entry points come straight from the nanobind core (the same C++
``HumidAir::`` functions the legacy Cython module wrapped), so this keeps the
legacy import path ``CoolProp.HumidAirProp.HAPropsSI`` working unchanged.

The legacy non-SI ``HAProps`` was removed in v8 (CoolProp is SI-only); calling
it raises a clear error pointing at ``HAPropsSI`` rather than silently 404-ing.
"""
from .CoolProp import HAPropsSI, HAProps_Aux, cair_sat

__all__ = ['HAPropsSI', 'HAProps_Aux', 'cair_sat', 'HAProps']


def HAProps(*args, **kwargs):
    """Removed in v8 -- use :func:`HAPropsSI` (SI units) instead."""
    raise NotImplementedError(
        "HAProps (non-SI humid-air properties) was removed in CoolProp v8; "
        "use HAPropsSI with SI units instead."
    )
