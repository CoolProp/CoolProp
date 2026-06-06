"""
Package init for the nanobind-based CoolProp build (COOLPROP_NANOBIND=ON).

The hand-written Cython AbstractState layer is gone; `CoolProp.CoolProp` is the
nanobind core (which also exports the `_capi` PyCapsule), and `CoolProp.State` is
the link-free capsule shim.  This intentionally stays minimal -- the full set of
legacy conveniences (AbstractState/HumidAirProp re-exports, etc.) is the nanobind
parity-completion work (bd CoolProp-q2sh), not State packaging.
"""
from __future__ import absolute_import

# Order matters: import the core first so its `_capi` capsule exists before the
# State shim (which does `import CoolProp`) is imported below.  `import *` also
# brings get_global_param_string into the package namespace.
from .CoolProp import *      # nanobind core: PropsSI, the bound API, enums, _capi
from . import State          # the capsule State shim (CoolProp.State.State)
from . import constants      # generated runtime enum values
from .constants import *

__version__ = get_global_param_string('version')
__gitrevision__ = get_global_param_string('gitrevision')


def get_include_directory():
    """Path to the bundled CoolProp headers (for compiling extensions that
    cimport CoolProp.State or include the C-ABI header)."""
    import os
    head, _ = os.path.split(__file__)
    return os.path.join(head, 'include')
