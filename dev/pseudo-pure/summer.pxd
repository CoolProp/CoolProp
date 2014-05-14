cimport numpy as np
cimport cython

cdef class PPF_summer:
    cdef public np.ndarray n,t,d,L,tau,delta
    cdef int IL0
    
    cpdef set_constants(self, np.ndarray t, np.ndarray d, np.ndarray l, int IL0)
    cpdef dphir_dDelta(self)