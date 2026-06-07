# distutils: language = c++
# cython: language_level=3
"""
Proves the frozen `State.pxd` is usable the way PDSim uses it: `cimport` the
cdef class, hold a `cdef State`, and call its methods (and `.pAS`'s) at the
C/vtable level -- not through the Python wrapper.  If this compiles, PDSim's
`from CoolProp.State cimport State; cdef State S; S.get_p()` pattern compiles
against the shim.
"""
from State cimport State, _AbstractStateView


cpdef dict check(object fluid, double T, double rho):
    cdef State S = State(fluid, {'T': T, 'D': rho})   # cdef-typed, cimported
    S.update_Trho(T, rho)
    cdef _AbstractStateView pas = S.pAS               # cdef attr access
    return {
        'p': S.get_p(),            # kPa (cdef call)
        'h': S.get_h(),            # kJ/kg
        'T_': S.T_,                # cached cdef attr
        'rhomass': pas.rhomass(),  # SI, through cdef pAS
        'dpdT': S.get_dpdT(),
    }
