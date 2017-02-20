# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy 
cimport numpy as np
import numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import the definitions
from TMR cimport *

# Include the mpi4py header
cdef extern from "mpi-compat.h":
   pass

cdef class Vertex:
    cdef TMRVertex *ptr
    def __cinit__(self):
        self.ptr = NULL
    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

cdef class Curve:
    cdef TMRCurve *ptr
    def __cinit__(self):
        self.ptr = NULL
    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
    def writeToVTK(self, char* filename):
        self.ptr.writeToVTK(filename)

cdef class Pcurve:
    cdef TMRPcurve *ptr
    def __cinit__(self):
        self.ptr = NULL
    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

cdef class Surface:
    cdef TMRSurface *ptr
    def __cinit__(self):
        self.ptr = NULL
    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
    def writeToVTK(self, char* filename):
        self.ptr.writeToVTK(filename)

cdef _init_Vertex(TMRVertex *ptr):
    vertex = Vertex()
    vertex.ptr = ptr
    vertex.ptr.incref()
    return vertex

cdef _init_Curve(TMRCurve *ptr):
    curve = Curve()
    curve.ptr = ptr
    curve.ptr.incref()
    return curve

cdef _init_Surface(TMRSurface *ptr):
    surface = Surface()
    surface.ptr = ptr
    surface.ptr.incref()
    return surface

cdef class BsplineCurve(Curve):
    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] pts, int k=4):
        cdef int nctl = pts.shape[0]
        cdef ku = k
        if ku > nctl:
            ku = nctl
        cdef TMRPoint* p = <TMRPoint*>malloc(nctl*sizeof(TMRPoint))
        for i in range(nctl):
            p[i].x = pts[i,0]
            p[i].y = pts[i,1]
            p[i].z = pts[i,2]
        self.ptr = new TMRBsplineCurve(nctl, ku, p)
        self.ptr.incref()
        free(p)

cdef class BsplinePcurve(Pcurve):
    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] pts, int k=4):
        cdef int nctl = pts.shape[0]
        cdef ku = k
        if ku > nctl:
            ku = nctl
        self.ptr = new TMRBsplinePcurve(nctl, ku, <double*>pts.data)
        self.ptr.incref()

cdef class BsplineSurface(Surface):
    def __cinit__(self, np.ndarray[double, ndim=3, mode='c'] pts, int ku=4, int kv=4):
        cdef int nx = pts.shape[0]
        cdef int ny = pts.shape[1]
        cdef kx = ku
        cdef ky = kv
        if kx > nx:
            kx = nx
        if ky > ny:
            ky = ny
        cdef TMRPoint* p = <TMRPoint*>malloc(nx*ny*sizeof(TMRPoint))
        for j in range(ny):
            for i in range(nx):
                p[i + j*nx].x = pts[i,j,0]
                p[i + j*nx].y = pts[i,j,1]
                p[i].z = pts[i,2]
        self.ptr = new TMRBsplineSurface(nx, ny, kx, ky, p)
        self.ptr.incref()
        free(p)

cdef class CurveInterpolation:
    cdef TMRCurveInterpolation *ptr
    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] pts):
        cdef int nctl = pts.shape[0]
        cdef TMRPoint* p = <TMRPoint*>malloc(nctl*sizeof(TMRPoint))
        for i in range(nctl):
            p[i].x = pts[i,0]
            p[i].y = pts[i,1]
            p[i].z = pts[i,2]
        self.ptr = new TMRCurveInterpolation(p, nctl)
        self.ptr.incref()
        free(p)

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setNumControlPoints(self, int nctl):
        self.ptr.setNumControlPoints(nctl)
        return

    def createCurve(self, int ku):
        cdef TMRBsplineCurve *curve = self.ptr.createCurve(ku)
        return _init_Curve(curve)

cdef class CurveLofter:
    cdef TMRCurveLofter *ptr
    def __cinit__(self, curves):
        cdef int ncurves = len(curves)
        cdef TMRBsplineCurve **crvs = NULL
        cdef TMRBsplineCurve *bspline = NULL
        crvs = <TMRBsplineCurve**>malloc(ncurves*sizeof(TMRBsplineCurve*))
        for i in range(ncurves):
            bspline =  _dynamicBsplineCurve((<Curve>curves[i]).ptr)
            if bspline != NULL:
               crvs[i] = bspline
            else:
               raise ValueError('All lofting curves must be of type BsplineCurve')
        self.ptr = new TMRCurveLofter(crvs, ncurves)
        self.ptr.incref()
        free(crvs)

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def createSurface(self, int kv):
        cdef TMRSurface *surf = self.ptr.createSurface(kv)
        return _init_Surface(surf)
