from libcpp cimport bool 
from libcpp.string cimport string

# A header defining the AbstractState class
cimport cAbstractState

cdef class AbstractState:
    cdef cAbstractState.AbstractState *thisptr     # hold a C++ instance which we're wrapping
    cpdef update(self, long iInput1, double Value1, double Value2)
    
    ## ---------------------------------------- 
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double T(self)
    cpdef double p(self)
    cpdef double rhomolar(self)
    cpdef double hmolar(self)
    cpdef double smolar(self)
    cpdef double cpmolar(self)
    cpdef double cvmolar(self)
    cpdef double speed_sound(self)
    
    cpdef double keyed_output(self, long)