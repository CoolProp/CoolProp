# Frozen cimport contract for the compat `State` shim (full surface).
# Mirrors the members PDSim cimports, per dev/pdsim_cimport_contract/SURFACE.md.
# Method args are typed `long` so PDSim's `constants_header` enum values pass
# through by implicit enum->int conversion (same as the old State class).

cdef class _AbstractStateView:
    cdef void* handle
    cdef object fluids
    cpdef double keyed_output(self, long key) except *
    cpdef double rhomass(self) except *
    cpdef double cpmass(self) except *
    cpdef double cvmass(self) except *
    cpdef double T(self) except *
    cpdef double p(self) except *
    cpdef update(self, long input_pair, double value1, double value2)
    cpdef specify_phase(self, long phase)
    cpdef double first_partial_deriv(self, long Of, long Wrt, long Constant) except *
    cpdef fluid_names(self)

cdef class State:
    cdef void* handle
    cdef object _fluids
    cdef object _fractions   # parsed mixture mole fractions (None for pure), for copy()
    cdef readonly _AbstractStateView pAS
    cdef readonly double T_, p_, rho_
    cdef readonly bytes Fluid, phase
    cdef _refresh(self)
    cpdef set_Fluid(self, Fluid, backend)
    cpdef long Phase(self) except *
    cpdef get_Tsat(self, double Q=*)
    cpdef get_subcooling(self)
    cpdef get_superheat(self)
    cpdef update(self, dict params)
    cpdef update_Trho(self, double T, double rho)
    cpdef update_ph(self, double p, double h)
    cpdef State copy(self)
    cpdef double Props(self, long key) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_h(self) except *
    cpdef double get_rho(self) except *
    cpdef double get_s(self) except *
    cpdef double get_u(self) except *
    cpdef double get_cp(self) except *
    cpdef double get_cp0(self) except *
    cpdef double get_cv(self) except *
    cpdef double get_MM(self) except *
    cpdef double get_dpdT(self) except *
    cpdef double get_visc(self) except *
    cpdef double get_cond(self) except *
    cpdef double get_speed_sound(self) except *
    cpdef double get_Q(self) except *
