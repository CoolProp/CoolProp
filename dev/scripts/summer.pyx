
cimport cython

import numpy as np


@cython.boundscheck(False)
cpdef  sum_function( double [:] B, double [:] x, double [:] n, double [:] out):

    
    cdef int i,j
    cdef double aa
    cdef double do
    cdef int Nx = len(x)
    cdef int Nn = len(n)
    
    out[:] = 0
    
    for i in xrange(Nx):
        for j in range(Nn):
            out[i] += B[j]*x[i]**n[j]
        