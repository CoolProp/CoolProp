from __future__ import absolute_import
"""
Package init for the nanobind-based CoolProp build (COOLPROP_NANOBIND=ON).

The hand-written Cython AbstractState layer is gone; `CoolProp.CoolProp` is the
nanobind core (which also exports the `_capi` PyCapsule), and `CoolProp.State` is
the link-free capsule shim.  This intentionally stays minimal -- the full set of
legacy conveniences (AbstractState/HumidAirProp re-exports, etc.) is the nanobind
parity-completion work (bd CoolProp-q2sh), not State packaging.
"""
# Order matters: import the core first so `_capi` is present on the package
# before the State shim (which does `import CoolProp`) is imported below.
from .CoolProp import *      # nanobind core: PropsSI, the bound API, enums, _capi
from . import CoolProp       # submodule handle for the helpers below
from . import State          # the capsule State shim (CoolProp.State.State)

try:
    from . import constants
    from .constants import *
except Exception:
    pass

__version__ = CoolProp.get_global_param_string('version')
__gitrevision__ = CoolProp.get_global_param_string('gitrevision')


def get_include_directory():
    """Path to the bundled CoolProp headers (for compiling extensions that
    cimport CoolProp.State or include the C-ABI header)."""
    import os
    head, _ = os.path.split(__file__)
    return os.path.join(head, 'include')
