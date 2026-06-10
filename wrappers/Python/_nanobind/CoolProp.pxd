# Re-export the cdef `State` class on the `CoolProp.CoolProp` cimport path.
#
# PDSim's core/state_flooded does `from CoolProp.CoolProp cimport State` (see
# dev/pdsim_cimport_contract/SURFACE.md).  In the legacy Cython package the
# `State` cdef class lived in the CoolProp.CoolProp module itself; in v8 the
# nanobind core IS CoolProp.CoolProp (a compiled extension with no cimportable
# cdef classes) and `State` is the frozen capsule shim in CoolProp.State.  This
# .pxd bridges the two: re-exporting `State` here lets the documented
# `from CoolProp.CoolProp cimport State` resolve to the same cdef class as
# `from CoolProp.State cimport State`, at both Cython compile time and runtime
# (the type object is sourced from the CoolProp.State module).  CoolProp-1tbe.6.
from CoolProp.State cimport State
