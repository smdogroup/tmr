# For MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import string stuff
from libc.string cimport const_char

# Import numpy 
cimport numpy as np
import numpy as np

# Import TACS
from paropt.ParOpt cimport *
from tacs.TACS cimport *
from tacs.constitutive cimport *

cdef extern from "<stdint.h>":
    ctypedef signed int int32_t

cdef extern from "TMRBase.h":
    enum:
        TMR_MAX_LEVEL"TMR_MAX_LEVEL"
        
    cdef cppclass TMREntity:
        void incref()
        void decref()
        void setAttribute(char*)
        int getEntityId()
        
    cdef cppclass TMRPoint:
        double x
        double y
        double z

    cdef cppclass TMRIndexWeight:
        int index
        double weight

    void TMRInitialize()
    void TMRFinalize()

    enum TMRInterpolationType:
        TMR_UNIFORM_POINTS
        TMR_GAUSS_LOBATTO_POINTS

cdef extern from "TMRTopology.h":
    cdef cppclass TMRTopology(TMREntity):
        TMRTopology(MPI_Comm, TMRModel*)
        void getVolume(int, TMRVolume**)
        void getFace(int, TMRFace**)
        void getEdge(int, TMREdge**)
        void getVertex(int, TMRVertex**)

    cdef cppclass TMRVertex(TMREntity):
        TMRVertex()
        int evalPoint(TMRPoint*)
        int setNodeNum(int*)
        int getNodeNum(int*)
    
    cdef cppclass TMREdge(TMREntity):
        TMREdge()
        void setSource(TMREdge*)
        void getSource(TMREdge**)
        void setVertices(TMRVertex*, TMRVertex*)
        void getVertices(TMRVertex**, TMRVertex**)
        void setMesh(TMREdgeMesh*)
        void getMesh(TMREdgeMesh**)
        void writeToVTK(char*)
        
    cdef cppclass TMRFace(TMREntity):
        TMRFace()
        TMRFace(int)
        void setSource(TMRVolume*, TMRFace*)
        void getSource(int*, TMRVolume**, TMRFace**)
        int getNumEdgeLoops()
        void addEdgeLoop(TMREdgeLoop*)
        void getEdgeLoop(int, TMREdgeLoop**)
        void setMesh(TMRFaceMesh*)
        void getMesh(TMRFaceMesh**)
        void writeToVTK(char*)

    cdef cppclass TMREdgeLoop(TMREntity):
        TMREdgeLoop(int, TMREdge**, const int*)
        void getEdgeLoop(int*, TMREdge***, const int**)

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
        TMREdgeFromFace(TMRFace*, TMRPcurve*, int)
        void addEdgeFromFace(TMRFace*, TMRPcurve*)

    cdef cppclass TMRSplitEdge(TMREdge):
        TMRSplitEdge(TMREdge*, double, double)
        TMRSplitEdge(TMREdge*, TMRVertex*, TMRVertex*)

    cdef cppclass TMREdgeFromCurve(TMREdge):
        TMREdgeFromCurve(TMRCurve*)

    cdef cppclass TMRFaceFromSurface(TMRFace):
        TMRFaceFromSurface(TMRSurface*)

    cdef cppclass TMRTFIFace(TMRFace):
        TMRTFIFace(TMREdge**, const int*, TMRVertex**)

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
    void _deleteMe "delete [] "(int *array)
    TMRBsplineCurve* _dynamicBsplineCurve "dynamic_cast<TMRBsplineCurve*>"(TMRCurve*)
    TMREdgeFromFace* _dynamicEdgeFromFace "dynamic_cast<TMREdgeFromFace*>"(TMREdge*)
    TMRTopoProblem* _dynamicTopoProblem "dynamic_cast<TMRTopoProblem*>"(ParOptProblem*)
    ParOptBVecWrap* _dynamicParOptBVecWrap "dynamic_cast<ParOptBVecWrap*>"(ParOptVec*)

cdef extern from "TMRMesh.h":
    enum TMRFaceMeshType:
        TMR_NO_MESH
        TMR_STRUCTURED
        TMR_UNSTRUCTURED
        TMR_TRIANGLE

    cdef cppclass TMRElementFeatureSize(TMREntity):
        TMRElementFeatureSize()
        TMRElementFeatureSize(double)

    cdef cppclass TMRLinearElementSize(TMRElementFeatureSize):
        TMRLinearElementSize(double, double,
                             double, double, double, double)

    cdef cppclass TMRMesh(TMREntity):     
        TMRMesh(MPI_Comm, TMRModel*)
        void mesh(TMRMeshOptions, double)
        void mesh(TMRMeshOptions, TMRElementFeatureSize*)
        int getMeshPoints(TMRPoint**)
        int getQuadConnectivity(int*, const int**)
        int getTriConnectivity(int*, const int**)
        int getHexConnectivity(int*, const int**)
        int getTetConnectivity(int*, const int**)

        TMRModel *createModelFromMesh()
        void writeToVTK(const char*, int)
        void writeToBDF(const char*, int)

    cdef cppclass TMREdgeMesh(TMREntity):
        TMREdgeMesh(MPI_Comm, TMREdge*)
        void mesh(TMRMeshOptions, TMRElementFeatureSize*)

    cdef cppclass TMRFaceMesh(TMREntity):
        TMRFaceMesh(MPI_Comm, TMRFace*)
        void mesh(TMRMeshOptions, TMRElementFeatureSize*)

    cdef cppclass TMRMeshOptions:
        TMRMeshOptions()
        TMRFaceMeshType mesh_type_default
        int triangularize_print_level
        int triangularize_print_iter
        int write_mesh_quality_histogram
        int num_smoothing_steps
        double frontal_quality_factor
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
        TMRQuadForest(MPI_Comm, int, TMRInterpolationType)
        MPI_Comm getMPIComm()
        void setTopology(TMRTopology*)
        void setConnectivity(int, const int*, int)
        void setFullConnectivity(int, int, int, const int*, const int*)
        void repartition()
        void createTrees(int)
        void createRandomTrees(int, int, int)
        void refine(int*, int, int)
        TMRQuadForest *duplicate()
        TMRQuadForest *coarsen()
        void balance(int)
        void createNodes()
        int getMeshOrder()
        void setMeshOrder(int, TMRInterpolationType)
        void getNodeConn(const int**, int*)
        int getDepNodeConn(const int**, const int**, const double**)
        TMRQuadrantArray* getQuadsWithAttribute(const char*)
        int getNodesWithAttribute(const char*, int**)
        void createInterpolation(TMRQuadForest*, TACSBVecInterp*)
        int getOwnedNodeRange(const int**)
        void getQuadrants(TMRQuadrantArray**)
        int getPoints(TMRPoint**)
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
        TMROctant* contains(TMROctant*, int)
        
cdef extern from "TMROctForest.h":
    cdef cppclass TMROctForest(TMREntity):
        TMROctForest(MPI_Comm, int, TMRInterpolationType)
        MPI_Comm getMPIComm()
        void setTopology(TMRTopology*)
        void setConnectivity(int, const int*, int)
        void setFullConnectivity(int, int, int, const int*, const int*)
        void repartition()
        void createTrees(int)
        void createRandomTrees(int, int, int)
        void refine(int*, int, int)
        TMROctForest *duplicate()
        TMROctForest *coarsen()
        void balance(int)
        void createNodes()
        int getMeshOrder()
        void setMeshOrder(int, TMRInterpolationType)
        void getNodeConn(const int**, int*)
        int getDepNodeConn(const int**, const int**, const double**)
        TMROctantArray* getOctsWithAttribute(const char*)
        int getNodesWithAttribute(const char*, int**)
        void createInterpolation(TMROctForest*, TACSBVecInterp*)
        int getOwnedNodeRange(const int**)
        void getOctants(TMROctantArray**)
        int getPoints(TMRPoint**)
        void writeToVTK(const char*)
        void writeForestToVTK(const char*)

cdef extern from "TMR_TACSCreator.h":
    cdef cppclass TMRBoundaryConditions(TMREntity):
        TMRBoundaryConditions()
        void addBoundaryCondition(const char*, int, const int*,
                                  const TacsScalar*)
        int getNumBoundaryConditions()

    cdef cppclass TMRQuadTACSCreator(TMREntity):
        TMRQuadTACSCreator(TMRBoundaryConditions*)

    cdef cppclass TMROctTACSCreator(TMREntity):
        TMROctTACSCreator(TMRBoundaryConditions*)

cdef extern from "TMROctStiffness.h":
    cdef cppclass TMRStiffnessProperties(TMREntity):
        TMRStiffnessProperties(int, TacsScalar*, TacsScalar*, TacsScalar*)
        int nmats
        
    cdef cppclass TMROctStiffness(SolidStiffness):
        TMROctStiffness(TMRIndexWeight*, int, TMRStiffnessProperties*,
                        double, double)
        
cdef extern from "TMRQuadStiffness.h":
    cdef cppclass TMRQuadStiffness(PlaneStressStiffness):
        TMRQuadStiffness(TMRIndexWeight*, int, TacsScalar, TacsScalar,
                         TacsScalar, double)

cdef extern from "SolidShellWrapper.h":
    cdef cppclass SolidShellWrapper(TACSElement):
        pass
    
cdef extern from "TMROpenCascade.h":
    cdef TMRModel* TMR_LoadModelFromSTEPFile(const char*, int)

cdef extern from "TMR_RefinementTools.h":
    void TMR_CreateTACSMg(int, TACSAssembler**,
                          TMROctForest**, TACSMg**, int, int)
    void TMR_CreateTACSMg(int, TACSAssembler**,
                          TMRQuadForest**, TACSMg**, int, int)
    void TMR_ComputeReconSolution(TACSAssembler*, TMRQuadForest*,
                                  TACSAssembler*, TACSBVec*, TACSBVec*)
    void TMR_ComputeReconSolution(TACSAssembler*, TMROctForest*,
                                  TACSAssembler*, TACSBVec*, TACSBVec*)
    double TMR_StrainEnergyErrorEst(TMRQuadForest*, TACSAssembler*, double*)
    double TMR_StrainEnergyErrorEst(TMROctForest*, TACSAssembler*, double*)
    double TMR_AdjointErrorEst(TACSAssembler*, TACSAssembler*, TACSBVec*,
                               TMRQuadForest*, double*, double*)
    double TMR_AdjointErrorEst(TACSAssembler*, TACSAssembler*, TACSBVec*,
                               TMROctForest*, double*, double*)

cdef extern from "TMRCyCreator.h":
    ctypedef TACSElement* (*createquadelements)(void*, int, TMRQuadrant*)
    ctypedef TACSElement* (*createoctelements)(void*, int, TMROctant*)
    ctypedef TACSElement* (*createquadtopoelements)(
        void*, int, TMRQuadrant*, TMRIndexWeight*, int)
    ctypedef TACSElement* (*createocttopoelements)(
        void*, int, TMROctant*, TMRIndexWeight*, int)

    cdef cppclass TMRCyQuadCreator(TMREntity):
        TMRCyQuadCreator(TMRBoundaryConditions*)
        void setSelfPointer(void*)
        void setCreateQuadElement( 
            TACSElement* (*createquadelements)(void*, int, TMRQuadrant*))
        TACSAssembler *createTACS(TMRQuadForest*, OrderingType)

    cdef cppclass TMRCyOctCreator(TMREntity):
        TMRCyOctCreator(TMRBoundaryConditions*)
        void setSelfPointer(void*)
        void setCreateOctElement( 
            TACSElement* (*createoctelements)(void*, int, TMROctant*))
        TACSAssembler *createTACS(TMROctForest*, OrderingType)

    cdef cppclass TMRCyTopoQuadCreator(TMREntity):
        TMRCyTopoQuadCreator(TMRBoundaryConditions*, TMRQuadForest*)
        void setSelfPointer(void*)
        void setCreateQuadTopoElement( 
            TACSElement* (*createquadtopoelements)(
                void*, int, TMRQuadrant*, TMRIndexWeight*, int))
        TACSAssembler *createTACS(int, TMRQuadForest*, OrderingType)
        void getFilter(TMRQuadForest**)
        void getMap(TACSVarMap**)
        void getIndices(TACSBVecIndices**)

    cdef cppclass TMRCyTopoOctCreator(TMREntity):
        TMRCyTopoOctCreator(TMRBoundaryConditions*, TMROctForest*)
        void setSelfPointer(void*)
        void setCreateOctTopoElement( 
            TACSElement* (*createocttopoelements)(
                void*, int, TMROctant*, TMRIndexWeight*, int))
        TACSAssembler *createTACS(int, TMROctForest*, OrderingType, double)
        void getFilter(TMROctForest**)
        void getMap(TACSVarMap**)
        void getIndices(TACSBVecIndices**)

cdef extern from "TMRTopoProblem.h":
    cdef cppclass TMRTopoProblem(ParOptProblem):
        TMRTopoProblem(int, TACSAssembler**, TMROctForest**, 
                       TACSVarMap**, TACSBVecIndices**, TACSMg*, int)
        TMRTopoProblem(int, TACSAssembler**, TMRQuadForest**, 
                       TACSVarMap**, TACSBVecIndices**, TACSMg*, int)
        void setLoadCases(TACSBVec**, int)
        int getNumLoadCases()
        void addConstraints(int, TACSFunction**,
                            const TacsScalar*, const TacsScalar*, int)
        void addLinearConstraints(ParOptVec**, TacsScalar*, int)
        void addFrequencyConstraint(double, int, TacsScalar,
                                    TacsScalar, TacsScalar, int, double)
        void setObjective(const TacsScalar*)
        void setObjective(const TacsScalar*, TACSFunction**)
        void initialize()
        void setPrefix(const char*)
        void setInitDesignVars(ParOptVec*,ParOptVec*,ParOptVec*)
        void setIterationCounter(int)
        TACSBVec* createVolumeVec()
        TACSBVec* createAreaVec()
        ParOptVec* createDesignVec()

    cdef cppclass ParOptBVecWrap(ParOptVec):
        ParOptBVecWrap(TACSBVec*)
        TACSBVec *vec

cdef extern from "TMR_STLTools.h":
    int TMR_GenerateBinFile(const char*, TMROctForest*,
                            TACSBVec*, int, double)
