#cython: embedsignature = True, c_string_type=str, c_string_encoding=ascii
from __future__ import division
        
cdef class AbstractState:
    
    def __cinit__(self, string backend, string fluid):
        self.thisptr = cAbstractState.factory(backend, fluid)
        
    def __dealloc__(self):
        del self.thisptr
    
    cpdef update(self, long ipair, double Value1, double Value2):
        self.thisptr.update(ipair, Value1, Value2)
    
    ## ----------------------------------------	
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double keyed_output(self, long iOutput) except *: 
        return self.thisptr.keyed_output(iOutput)
    
    cpdef double T(self) except *: 
        return self.thisptr.T()
    cpdef double p(self) except *: 
        return self.thisptr.p()
    cpdef double rhomolar(self) except *: 
        return self.thisptr.rhomolar()
#     cpdef double rhomass(self) except *: 
#         return self.thisptr.rhomass()
    cpdef double hmolar(self) except *: 
        return self.thisptr.hmolar()
    cpdef double smolar(self) except *: 
        return self.thisptr.smolar()
    cpdef double cpmolar(self) except *: 
        return self.thisptr.cpmolar()
    cpdef double cvmolar(self) except *: 
        return self.thisptr.cvmolar()
#     cpdef double hmass(self) except *: 
#         return self.thisptr.hmass()
#     cpdef double smass(self) except *: 
#         return self.thisptr.smass()
#     cpdef double cpmass(self) except *: 
#         return self.thisptr.cpmass()
#     cpdef double cvmass(self) except *: 
#         return self.thisptr.cvmass()
    cpdef double speed_sound(self) except *: 
        return self.thisptr.speed_sound()
    
        