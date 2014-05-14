
cimport cython

import numpy as np
from libc.math cimport exp, log

@cython.boundscheck(False)
cdef class PPF_summer:
    
    def __init__(self, tau, delta):
        self.tau = tau
        self.delta = delta
        
    cpdef set_constants(self, np.ndarray t, np.ndarray d, np.ndarray l, int IL0):
        self.t = t
        self.d = d
        self.L = l
        self.IL0 = IL0

    cpdef dphir_dDelta(self):
        cdef int i, j
        
        dphir_dDelta = np.zeros_like(self.tau)
        cdef double [:] dphir_dDelta_view = dphir_dDelta
        cdef long [:] d = self.d
        cdef double [:] t = self.t
        cdef long [:] L = self.L
        cdef double [:] n = self.n
        
        cdef double [:] tau = self.tau
        cdef double [:] delta = self.delta
        
        for j in xrange(len(self.tau)):
            log_tau = log(tau[j])
            log_delta = log(delta[j])
            for i in xrange(len(self.t)):
                if i < self.IL0:
                    dphir_dDelta_view[j] += d[i]*n[i]*exp((d[i]-1)*log_delta+t[i]*log_tau)
                else:
                    dphir_dDelta_view[j] += n[i]*exp((d[i]-1)*log_delta+t[i]*log_tau-delta[j]**L[i])*(d[i]-L[i]*delta[j]**L[i])
                
        return dphir_dDelta
        