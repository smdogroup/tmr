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

cdef extern from "TMRTopology.h":
    cdef cppclass TMRTopology(TMREntity):
        TMRTopology(MPI_Comm, TMRModel*)
            
    cdef cppclass TMRVertex(TMREntity):
        pass
    
    cdef cppclass TMREdge(TMREntity):
        TMREdge()
        void setSource(TMREdge*)
        void getSource(TMREdge**)
        void setVertices(TMRVertex*, TMRVertex*)
        void getVertices(TMRVertex**, TMRVertex**)
        void writeToVTK(char*)
        
    cdef cppclass TMRFace(TMREntity):
        TMRFace()
        TMRFace(int)
        void setSource(TMRVolume*, TMRFace*)
        void getSource(int*, TMRVolume**, TMRFace**)
        int getNumEdgeLoops()
        void addEdgeLoop(TMREdgeLoop*)
        void writeToVTK(char*)

    cdef cppclass TMREdgeLoop(TMREntity):
        TMREdgeLoop(int, TMREdge**, const int*)
        
    cdef cppclass TMRVolume(TMREntity):
        TMRVolume(int, TMRFace**, const int*)
        void getFaces(int*, TMRFace***, const int**)
        void writeToVTK(char*)
        
    cdef cppclass TMRModel(TMREntity):
        TMRModel(int, TMRVertex**,
                 int, TMREdge**,
                 int, TMRFace**,
                 int, TMRVolume**)
        void getVertices(int*,TMRVertex***)
        void getEdges(int*,TMREdge***)
        void getFaces(int*,TMRFace***)
        void getVolumes(int*,TMRVolume***)
        
cdef extern from "TMRGeometry.h":
    cdef cppclass TMRCurve(TMREntity):
        void writeToVTK(char*)

    cdef cppclass TMRSurface(TMREntity):
        void addCurveSegment(int, TMRCurve**, int*)
        void writeToVTK(char*)

    cdef cppclass TMRPcurve(TMREntity):
        pass

cdef extern from "TMRNativeTopology.h":
    cdef cppclass TMRVertexFromPoint(TMRVertex):
        TMRVertexFromPoint(TMRPoint)

    cdef cppclass TMRVertexFromEdge(TMRVertex):
        TMRVertexFromEdge(TMREdge*, double)
        TMRVertexFromEdge(TMREdge*, TMRPoint)

    cdef cppclass TMRVertexFromFace(TMRVertex):
        TMRVertexFromFace(TMRFace*, double, double)
        TMRVertexFromFace(TMRFace*, TMRPoint)

    cdef cppclass TMREdgeFromFace(TMREdge):
        TMREdgeFromFace(TMRFace*, TMRPcurve*)

    cdef cppclass TMRSplitEdge(TMREdge):
        TMRSplitEdge(TMREdge*, double, double)
        TMRSplitEdge(TMREdge*, TMRVertex*, TMRVertex*)

    cdef cppclass TMREdgeFromCurve(TMREdge):
        TMREdgeFromCurve(TMRCurve*)

    cdef cppclass TMRFaceFromSurface(TMRFace):
        TMRFaceFromSurface(TMRSurface*)

cdef extern from "TMRBspline.h":
    cdef cppclass TMRBsplineCurve(TMRCurve):
        TMRBsplineCurve(int, int, TMRPoint*)
        TMRBsplineCurve(int, int, double*, TMRPoint*)
        TMRBsplineCurve(int, int, double*, double*, TMRPoint*)

    cdef cppclass TMRBsplineSurface(TMRSurface):
        TMRBsplineSurface(int, int, int, int, TMRPoint*)
        TMRBsplineSurface(int, int, int, int,
                          double*, double*, TMRPoint*)
        TMRBsplineSurface(int, int, int, int,
                          double*, double*, double*, TMRPoint*)

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
 
cdef extern from "TMRMesh.h":    
    cdef cppclass TMRMesh(TMREntity):     
        TMRMesh(MPI_Comm, TMRModel*)
        void mesh(double)
        void mesh(TMRMeshOptions, double)
        int getMeshPoints(TMRPoint**)
        int getMeshConnectivity(int*, const int**,
                                int*, const int**)
        TMRModel *createModelFromMesh()
        void writeToVTK(const char*, int)
        void writeToBDF(const char*, int)

    cdef cppclass TMREdgeMesh(TMREntity):
        TMREdgeMesh(MPI_Comm, TMREdge*)
        void mesh(TMRMeshOptions, double)

    cdef cppclass TMRFaceMesh(TMREntity):
        TMRFaceMesh(MPI_Comm, TMRFace*)
        void mesh(TMRMeshOptions, double)

    cdef cppclass TMRMeshOptions:
        TMRMeshOptions()
        int triangularize_print_level
        int triangularize_print_iter
        int write_mesh_quality_histogram
        int num_smoothing_steps
        double frontal_quality_factor

        # Write intermediate surface meshes to file
        int write_init_domain_triangle
        int write_triangularize_intermediate
        int write_pre_smooth_triangle
        int write_post_smooth_triangle
        int write_dual_recombine
        int write_pre_smooth_quad
        int write_post_smooth_quad
        int write_quad_dual

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
        void createRandomTrees(int, int, int)
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
        void writeToVTK(const char*)
        void writeForestToVTK(const char*)

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
        void createRandomTrees(int, int, int)
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
        void writeToVTK(const char*)
        void writeForestToVTK(const char*)

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
        TacsScalar rho, E, nu, q
        
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
    cdef TMRModel* TMR_LoadModelFromSTEPFile(const char*, int)
