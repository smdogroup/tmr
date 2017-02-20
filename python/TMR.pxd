# For MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import string stuff
from libc.string cimport const_char

# Import numpy 
cimport numpy as np
import numpy as np

cdef extern from "TMRBase.h":
    cdef cppclass TMREntity:
        void incref()
        void decref()
        void setAttribute(char *)

    cdef cppclass TMRPoint:
        double x
        double y
        double z

    void TMRInitialize()
    void TMRFinalize()

cdef extern from "TMRGeometry.h":
    cdef cppclass TMRVertex(TMREntity):
        pass

    cdef cppclass TMRCurve(TMREntity):
        void setVertices(TMRVertex*, TMRVertex*)
        void getVertices(TMRVertex**, TMRVertex**)
        void writeToVTK(char*)

    cdef cppclass TMRSurface(TMREntity):
        void addCurveSegment(int, TMRCurve**, int**)
        void writeToVTK(char*)

    cdef cppclass TMRPcurve(TMREntity):
        pass

    cdef cppclass TMRVertexFromPoint(TMRVertex):
        TMRVertexFromPoint(TMRPoint)

    cdef cppclass TMRVertexFromCurve(TMRVertex):
        TMRVertexFromCurve(TMRCurve*, double)
        TMRVertexFromCurve(TMRCurve*, TMRPoint)

    cdef cppclass TMRVertexFromSurface(TMRVertex):
        TMRVertexFromSurface(TMRSurface*, double, double)
        TMRVertexFromSurface(TMRSurface*, TMRPoint)

    cdef cppclass TMRCurveFromSurface(TMRCurve):
        TMRCurveFromSurface(TMRSurface*, TMRPcurve*)

    cdef cppclass TMRSplitCurve(TMRCurve):
        TMRSplitCurve(TMRCurve*, double, double)
        TMRSplitCurve(TMRCurve*, TMRVertex*, TMRVertex*)

cdef extern from "TMRBspline.h":
    cdef cppclass TMRBsplineCurve(TMRCurve):
        TMRBsplineCurve(int, int, TMRPoint*)
        TMRBsplineCurve(int, int, double*, TMRPoint*)
        TMRBsplineCurve(int, int, double*, double*, TMRPoint*)

    cdef cppclass TMRBsplineSurface(TMRSurface):
        TMRBsplineSurface(int, int, int, int, TMRPoint*)
        TMRBsplineSurface(int, int, int, int, double*, double*, TMRPoint*)
        TMRBsplineSurface(int, int, int, int, double*, double*, double*, TMRPoint*)

    cdef cppclass TMRBsplinePcurve(TMRPcurve):
        TMRBsplinePcurve(int, int, double*)
        TMRBsplinePcurve(int, int, double*, double*)
        TMRBsplinePcurve(int, int, double*, double*, double*)

    cdef cppclass TMRCurveInterpolation(TMREntity):
        TMRCurveInterpolation(TMRPoint*, int)
        void setNumControlPoints(int)
        TMRBsplineCurve *createCurve(int)
        
    cdef cppclass TMRCurveLofter(TMREntity):
        TMRCurveLofter(TMRBsplineCurve**, int)
        TMRBsplineSurface* createSurface(int)

cdef extern from "":
    TMRBsplineCurve* _dynamicBsplineCurve "dynamic_cast<TMRBsplineCurve*>"(TMRCurve*)
