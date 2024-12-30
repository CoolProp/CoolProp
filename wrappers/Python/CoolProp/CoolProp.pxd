import cython
cimport cython

from libcpp.vector cimport vector

from .typedefs cimport *

include "AbstractState.pxd"
       
cdef class State:
    cdef public AbstractState pAS
    cdef readonly bytes Fluid, phase
    cdef int iFluid,iParam1,iParam2,iOutput
    cdef double T_, rho_, p_
    
    cpdef set_Fluid(self, Fluid, backend)
    cpdef speed_test(self, int N)
    cpdef update(self, dict params)
    cpdef update_ph(self, double p, double h)
    cpdef update_Trho(self, double T, double rho)
    cpdef State copy(self)
    cpdef double Props(self, constants_header.parameters iOutput) except *
    cpdef long Phase(self) except *
    cpdef double get_Q(self) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_h(self) except *
    cpdef double get_rho(self) except *
    cpdef double get_s(self) except *
    cpdef double get_u(self) except *
    cpdef double get_visc(self) except *
    cpdef double get_cond(self) except *
    cpdef double get_cp(self) except *
    cpdef double get_cp0(self) except *
    cpdef double get_cv(self) except *
    cpdef double get_MM(self) except *
    cpdef double get_dpdT(self) except *
    cpdef double get_speed_sound(self) except *
    cpdef get_Tsat(self, double Q = *)
    cpdef get_subcooling(self)
    cpdef get_superheat(self)
