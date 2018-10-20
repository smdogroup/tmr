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
from tacs.functions cimport *

cdef extern from "mpi-compat.h":
    pass

cdef extern from "PoissonGradFunc.h":
    enum PoissonGradFunctionType"TACSPoissonGradFunction::PoissonGradFunctionType":
        POISSON_CONTINUOUS"TACSPoissonGradFunction::CONTINUOUS"
        POISSON_DISCRETE"TACSPoissonGradFunction::DISCRETE"
        POISSON_PNORM_CONTINUOUS"TACSPoissonGradFunction::PNORM_CONTINUOUS"
        POISSON_PNORM_DISCRETE"TACSPoissonGradFunction::PNORM_DISCRETE"
    
    cdef cppclass TACSPoissonGradFunction(TACSFunction):
        TACSPoissonGradFunction(TACSAssembler*, TacsScalar, const TacsScalar*)
        void setGradFunctionType(PoissonGradFunctionType)
            
cdef class KSPoissonGrad(Function):
    cdef TACSPoissonGradFunction *self_ptr
    def __cinit__(self, Assembler assembler, double ksweight, direction):
        '''Multimaterial topology optimization'''
        cdef TacsScalar a[2]
        a[0] = direction[0]
        a[1] = direction[1]
        self.self_ptr = new TACSPoissonGradFunction(assembler.ptr, ksweight, a)
        self.ptr = self.self_ptr
        self.ptr.incref()
        return

    def setGradFunctionType(self, ftype='continuous'):
        if ftype == 'discrete':
            self.self_ptr.setGradFunctionType(POISSON_DISCRETE)
        elif ftype == 'continuous':
            self.self_ptr.setGradFunctionType(POISSON_CONTINUOUS)
        elif ftype == 'pnorm-discrete':
            self.self_ptr.setGradFunctionType(POISSON_PNORM_DISCRETE)
        elif ftype == 'pnorm-continuous':
            self.self_ptr.setGradFunctionType(POISSON_PNORM_CONTINUOUS)
        return
