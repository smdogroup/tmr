# For MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import string stuff
from libc.string cimport const_char

# Import numpy 
cimport numpy as np
import numpy as np

# Import TACS
from tacs.TACS cimport *

cdef extern from "<stdint.h>":
    ctypedef signed int int32_t

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
        # void setVertices(TMRVertex*, TMRVertex*)
        # void getVertices(TMRVertex**, TMRVertex**)
        void writeToVTK(char*)

    cdef cppclass TMRSurface(TMREntity):
        void addCurveSegment(int, TMRCurve**, int*)
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

    cdef cppclass TMRGeometry(TMREntity):
        TMRGeometry(int, TMRVertex**, int, TMRCurve**, 
                    int, TMRSurface**)

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

cdef extern from "TMRTopology.h":
    cdef cppclass TMRTopology(TMREntity):
        pass
    cdef cppclass TMRVertex(TMREntity):
        pass
    cdef cppclass TMREdge(TMREntity):
        TMREdge()
        void setVertices(TMRVertex*,TMRVertex*)
        void getVertices(TMRVertex**, TMRVertex**)
        void writeToVTK(char*)
        
    cdef cppclass TMRFace(TMREntity):
        TMRFace(int)
        int getNumEdgeLoops()
        void setMaster(TMRFace*)
    cdef cppclass TMRVolume(TMREntity):
        TMRVolume(int, TMRFace**,const int*)
        void getFaces(int*,TMRFace***,const int**)
        void writeToVTK(char*)
        void updateOrientation()
    cdef cppclass TMRModel(TMREntity):
        TMRModel(int, TMRVertex**,
                 int, TMREdge**,
                 int, TMRFace**,
                 int, TMRVolume**)
        void getVertices(int*,TMRVertex***)
        void getEdges(int*,TMREdge***)
        void getFaces(int*,TMRFace***)
        void getVolumes(int*,TMRVolume***)
 
cdef extern from "TMRMesh.h":
    cdef cppclass TMRMesh(TMREntity):
        TMRMesh(MPI_Comm,TMRModel*)
        void mesh(double)
        int getMeshPoints(TMRPoint**)
        int getMeshConnectivity(int*,const int**,
                                int*,const int**)
        TMRModel *createModelFromMesh()
    cdef cppclass TMRMeshOptions:
        TMRMeshOptions()

cdef extern from "TMRQuadrant.h":
    cdef cppclass TMRQuadrant:
        int childId()
        void getSibling(int, TMRQuadrant*)
        void parent(TMRQuadrant*)
        void edgeNeighbor(int, TMRQuadrant*)
        void cornerNeighbor(int, TMRQuadrant*)
        int contains(TMRQuadrant*)
        int32_t x
        int32_t y
        int32_t level
        int32_t face
        int32_t tag

    cdef cppclass TMRQuadrantArray:
        TMRQuadrantArray(TMRQuadrant*, int)
        TMRQuadrantArray* duplicate()
        void getArray(TMRQuadrant**, int*)
        void sort()
        TMRQuadrant* contains(TMRQuadrant *q, int)

cdef extern from "TMRQuadForest.h":
    cdef cppclass TMRQuadForest(TMREntity):
        TMRQuadForest(MPI_Comm)
        MPI_Comm getMPIComm()
        void setTopology(TMRTopology*)
        void setConnectivity(int, const int*, int)
        void setFullConnectivity(int, int, int, const int*, const int*)
        void repartition()
        void createTrees(int)
        void createRandomTrees(int)
        void refine(int*)
        TMRQuadForest *duplicate()
        TMRQuadForest *coarsen()
        void balance(int)
        void createNodes(int)
        TMRQuadrantArray* getQuadsWithAttribute(const char*)
        TMRQuadrantArray* getNodesWithAttribute(const char*)
        void createMeshConn(const int**, const int*)
        int getDepNodeConn(const int**, const int**, const double**)
        void createInterpolation(TMRQuadForest*, TACSBVecInterp*)


cdef extern from "TMROctant.h":
    cdef cppclass TMROctant:
        int childId()
        void getSibling(int, TMROctant*)
        void parent(TMROctant*)
        void edgeNeighbor(int, TMROctant*)
        void cornerNeighbor(int, TMROctant*)
        void faceNeighbor(int, TMROctant)
        int contains(TMROctant*)
        int32_t x
        int32_t y
        int32_t z
        int32_t level
        int32_t block
        int32_t tag

    cdef cppclass TMROctantArray:
        TMROctantArray(TMROctant*, int)
        TMROctantArray* duplicate()
        void getArray(TMROctant**, int*)
        void sort()
        TMROctant* contains(TMROctant *q, int)
        
cdef extern from "TMROctForest.h":
    cdef cppclass TMROctForest(TMREntity):
        TMROctForest(MPI_Comm)
        MPI_Comm getMPIComm()
        void setTopology(TMRTopology*)
        void setConnectivity(int, const int*, int)
        void setFullConnectivity(int, int, int, const int*, const int*)
        void repartition()
        void createTrees(int)
        void createRandomTrees(int)
        void refine(int*)
        TMROctForest *duplicate()
        TMROctForest *coarsen()
        void balance(int)
        void createNodes(int)
        TMROctantArray* getOctsWithAttribute(const char*)
        TMROctantArray* getNodesWithAttribute(const char*)
        void createMeshConn(const int**, const int*)
        int getDepNodeConn(const int**, const int**, const double**)
        void createInterpolation(TMROctForest*, TACSBVecInterp*)

cdef extern from "TMR_TACSCreator.h":
    cdef cppclass TMRBoundaryConditions(TMREntity):
        TMRBoundaryConditions()
        void addBoundaryCondition(const char*,int,const int*,
                                  const TacsScalar*)
        int getNumBoundaryConditions()

    cdef cppclass TMRQuadTACSCreator(TMREntity):
        TMRQuadTACSCreator(TMRBoundaryConditions*)

    cdef cppclass TMROctTACSCreator(TMREntity):
        TMROctTACSCreator(TMRBoundaryConditions*)

cdef extern from "TMROctStiffness.h":
    cdef cppclass TMRStiffnessProperties:
        TMRStiffnessProperties()
    
# cdef extern from "TMR_TACSTopoCreator.c":
#     cdef cppclass TMROctTACSTopoCreator(TMROctTACSCreator):
#         TMROctTACSTopoCreator(TMRBoundaryConditions*,
#                               TMRStiffnessProperties,
#                               TMROctForest*,const char*,
#                               SolidShellWrapper*)

cdef extern from "SolidShellWrapper.h":
    cdef cppclass SolidShellWrapper(TACSElement):
        pass
    
cdef extern from "TMROpenCascade.h":
    cdef TMRModel* TMR_LoadModelFromSTEPFile(const char*) 
