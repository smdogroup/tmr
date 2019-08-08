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
from tacs.elements cimport *

cdef extern from "mpi-compat.h":
    pass

cdef extern from "PSTopo.h":
    cdef cppclass LobattoQuad2(TACSElement):
        LobattoQuad2(PlaneStressStiffness*)

    cdef cppclass LobattoQuad3(TACSElement):
        LobattoQuad3(PlaneStressStiffness*)

    cdef cppclass LobattoQuad4(TACSElement):
        LobattoQuad4(PlaneStressStiffness*)

    cdef cppclass LobattoQuad5(TACSElement):
        LobattoQuad5(PlaneStressStiffness*)
        
    cdef cppclass PSTopo(PlaneStressStiffness):
        PSTopo(TacsScalar, TacsScalar, TacsScalar, TacsScalar,
               double, double, int*, double*, int)
        TacsScalar getDensity()

    cdef cppclass PSTopo4(PlaneStressStiffness):
        PSTopo4(TacsScalar, TacsScalar, TacsScalar, TacsScalar,
                double, double, int**, double**, int*)

cdef class pstopo(PlaneStress):
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

    def getDensity(self):
        return self.self_ptr.getDensity()

cdef class LobattoQuad(Element):
    def __cinit__(self, int order, PlaneStress stiff):
        cdef PlaneStressStiffness *con = _dynamicPlaneStress(stiff.ptr)
        if order == 2:
            self.ptr = new LobattoQuad2(con)
            self.ptr.incref()
        elif order == 3:
            self.ptr = new LobattoQuad3(con)
            self.ptr.incref()
        elif order == 4:
            self.ptr = new LobattoQuad4(con)
            self.ptr.incref()
        elif order == 5:
            self.ptr = new LobattoQuad5(con)
            self.ptr.incref()
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

cdef class pstopo4(PlaneStress):
    cdef PSTopo4 *self_ptr
    def __cinit__(self, double rho, double E, double nu, double ys,
                  double q, double eps,
                  np.ndarray[int, ndim=1, mode='c'] nodes1,
                  np.ndarray[double, ndim=1, mode='c'] weights1,
                  np.ndarray[int, ndim=1, mode='c'] nodes2,
                  np.ndarray[double, ndim=1, mode='c'] weights2,
                  np.ndarray[int, ndim=1, mode='c'] nodes3,
                  np.ndarray[double, ndim=1, mode='c'] weights3,
                  np.ndarray[int, ndim=1, mode='c'] nodes4,
                  np.ndarray[double, ndim=1, mode='c'] weights4):
        '''Multimaterial topology optimization'''
        assert(len(nodes1) == len(weights1))
        assert(len(nodes2) == len(weights2))
        assert(len(nodes3) == len(weights3))
        assert(len(nodes4) == len(weights4))
        cdef int *nodes[4]
        cdef double *weights[4]
        cdef int nweights[4]

        # Set the weights
        nweights[0] = len(weights1)
        nweights[1] = len(weights2)
        nweights[2] = len(weights3)
        nweights[3] = len(weights4)

        # Set the weights
        nodes[0] = <int*>nodes1.data
        nodes[1] = <int*>nodes2.data
        nodes[2] = <int*>nodes3.data
        nodes[3] = <int*>nodes4.data

        # Set the weights
        weights[0] = <double*>weights1.data
        weights[1] = <double*>weights2.data
        weights[2] = <double*>weights3.data
        weights[3] = <double*>weights4.data

        self.self_ptr = new PSTopo4(rho, E, nu, ys, q, eps, nodes, weights, nweights)
        self.ptr = self.self_ptr
        self.ptr.incref()
        return
