# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import the TACS module
from tacs.TACS cimport *
from tacs.constitutive cimport *

cdef extern from "mpi-compat.h":
    pass

cdef extern from "PSTopo.h":        
    cdef cppclass PSTopo(PlaneStressStiffness):
        PSTopo(TacsScalar, TacsScalar, TacsScalar, TacsScalar,
               double, double, int*, double*, int)
    
cdef class ptopo(PlaneStress):
    cdef PSTopo *self_ptr
    def __cinit__(self, double rho, double E, double nu, double ys,
                  double q, double eps,
                  np.ndarray[int, ndim=1, mode='c'] nodes,
                  np.ndarray[double, ndim=1, mode='c'] weights):
        '''Multimaterial topology optimization'''
        assert(len(nodes) == len(weights))
        self.self_ptr = new PSTopo(rho, E, nu, ys, q, eps, <int*>nodes.data,
                                        <double*>weights.data, len(nodes))
        self.ptr = self.self_ptr
        self.ptr.incref()
        return
