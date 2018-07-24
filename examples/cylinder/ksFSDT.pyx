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

cdef extern from "ksFSDT.h":        
    cdef cppclass ksFSDTStiffness(FSDTStiffness):
        ksFSDTStiffness(TacsScalar, TacsScalar, TacsScalar, TacsScalar, 
                        TacsScalar, TacsScalar, TacsScalar, 
                        int, TacsScalar, TacsScalar)
    
cdef class ksFSDT(FSDT):
    cdef ksFSDTStiffness *self_ptr
    def __cinit__(self, ksweight, rho, E, nu, kcorr, ys, t, tNum, minT, maxT):
        '''Multimaterial topology optimization'''
        self.self_ptr = new ksFSDTStiffness(ksweight, rho, E, nu, kcorr, ys, 
                                            t, tNum, minT, maxT)
        self.ptr = self.self_ptr
        self.ptr.incref()
        return
