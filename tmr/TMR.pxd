# distutils: language=c++

# This file is part of the package TMR for adaptive mesh refinement.

# Copyright (C) 2015 Georgia Tech Research Corporation.
# Additional copyright (C) 2015 Graeme Kennedy.
# All rights reserved.

# TMR is licensed under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# For MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import string stuff
from libc.string cimport const_char
from libc.stdint cimport int32_t, int16_t

# Import the python version information
from cpython.version cimport PY_MAJOR_VERSION

# Import numpy
cimport numpy as np
import numpy as np

# Import TACS
from paropt.ParOpt cimport *
from tacs.TACS cimport *
from tacs.constitutive cimport *
from egads4py.egads cimport *

cdef inline char* tmr_convert_str_to_chars(s):
   if isinstance(s, unicode):
      s = (<unicode>s).encode('utf-8')
   return s

cdef inline str tmr_convert_char_to_str(const char* s):
    if s == NULL:
        return None
    elif PY_MAJOR_VERSION >= 3:
        return s.decode('utf-8')
    else:
        return str(s)

cdef extern from "TMRBase.h":
    enum:
        TMR_MAX_LEVEL"TMR_MAX_LEVEL"

    cdef cppclass TMREntity:
        void incref()
        void decref()
        void setName(const char*)
        const char* getName()
        int getEntityId()

    cdef cppclass TMRPoint:
        double x
        double y
        double z

    cdef cppclass TMRIndexWeight:
        int index
        double weight

    cdef cppclass TMR_STLTriangle:
        TMRPoint p[3]

    void TMRInitialize()
    int TMRIsInitialized()
    void TMRFinalize()

    enum TMRInterpolationType:
        TMR_UNIFORM_POINTS
        TMR_GAUSS_LOBATTO_POINTS
        TMR_BERNSTEIN_POINTS

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
        void setCopySource(TMRVertex*)
        void getCopySource(TMRVertex**)
        int isSame(TMRVertex*)

    cdef cppclass TMREdge(TMREntity):
        TMREdge()
        void setSource(TMREdge*)
        void getSource(TMREdge**)
        void getRange(double*, double*)
        int evalPoint(double, TMRPoint*)
        int evalDeriv(double, TMRPoint*, TMRPoint*)
        int invEvalPoint(TMRPoint, double*)
        void setVertices(TMRVertex*, TMRVertex*)
        void getVertices(TMRVertex**, TMRVertex**)
        int isDegenerate()
        void setMesh(TMREdgeMesh*)
        void getMesh(TMREdgeMesh**)
        void setCopySource(TMREdge*)
        void getCopySource(TMREdge**)
        void writeToVTK(char*)
        int isSame(TMREdge*)

    cdef cppclass TMRFace(TMREntity):
        TMRFace()
        TMRFace(int)
        void getRange(double*, double*, double*, double*)
        int evalPoint(double, double, TMRPoint*)
        void setSource(TMRVolume*, TMRFace*)
        void getSource(TMRVolume**, TMRFace**)
        int getNumEdgeLoops()
        void addEdgeLoop(int, TMREdgeLoop*)
        void getEdgeLoop(int, TMREdgeLoop**)
        void setMesh(TMRFaceMesh*)
        void getMesh(TMRFaceMesh**)
        void setCopySource(int, TMRFace*)
        void getCopySource(int*, TMRFace**)
        void writeToVTK(char*)
        int isSame(TMRFace*)

    cdef cppclass TMREdgeLoop(TMREntity):
        TMREdgeLoop(int, TMREdge**, const int*)
        void getEdgeLoop(int*, TMREdge***, const int**)

    cdef cppclass TMRVolume(TMREntity):
        TMRVolume(int, TMRFace**)
        void getFaces(int*, TMRFace***)
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

    cdef cppclass TMRTFIEdge(TMREdge):
        TMRTFIEdge(TMRVertex*, TMRVertex*)

    cdef cppclass TMRTFIFace(TMRFace):
        TMRTFIFace(TMREdge**, const int*, TMRVertex**)

cdef extern from "TMRBspline.h":
    cdef cppclass TMRBsplineCurve(TMRCurve):
        TMRBsplineCurve(int, int, TMRPoint*)
        TMRBsplineCurve(int, int, double*, TMRPoint*)
        TMRBsplineCurve(int, int, double*, double*, TMRPoint*)
        void getData(int*, int*, const double**,
                     const double**, const TMRPoint**)
        void split(double, TMRBsplineCurve**, TMRBsplineCurve**)

    cdef cppclass TMRBsplineSurface(TMRSurface):
        TMRBsplineSurface(int, int, int, int, TMRPoint*)
        TMRBsplineSurface(int, int, int, int,
                          double*, double*, TMRPoint*)
        TMRBsplineSurface(int, int, int, int,
                          double*, double*, double*, TMRPoint*)
        void getData(int*, int*, int*, int*, const double**, const double**,
                     const double**, const TMRPoint**)

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
    TMRBsplineSurface* _dynamicBsplineSurface "dynamic_cast<TMRBsplineSurface*>"(TMRSurface*)
    TMREdgeFromFace* _dynamicEdgeFromFace "dynamic_cast<TMREdgeFromFace*>"(TMREdge*)
    TMRTopoProblem* _dynamicTopoProblem "dynamic_cast<TMRTopoProblem*>"(ParOptProblem*)
    ParOptBVecWrap* _dynamicParOptBVecWrap "dynamic_cast<ParOptBVecWrap*>"(ParOptVec*)

cdef extern from "TMRFeatureSize.h":
    cdef cppclass TMRElementFeatureSize(TMREntity):
        TMRElementFeatureSize()
        TMRElementFeatureSize(double)
        double getFeatureSize(TMRPoint)

    cdef cppclass TMRLinearElementSize(TMRElementFeatureSize):
        TMRLinearElementSize(double, double,
                             double, double, double, double)

    cdef cppclass TMRBoxFeatureSize(TMRElementFeatureSize):
        TMRBoxFeatureSize(TMRPoint, TMRPoint, double, double)
        void addBox(TMRPoint, TMRPoint, double)

    cdef cppclass TMRPointFeatureSize(TMRElementFeatureSize):
        TMRPointFeatureSize(int, TMRPoint*, double*, double, double, int)

    cdef cppclass TMRPointLocator(TMREntity):
        TMRPointLocator(int, TMRPoint*)
        void locateClosest(int, TMRPoint, int*, int*, double*)

cdef class PointLocator:
    cdef TMRPointLocator *ptr

cdef extern from "TMREdgeMesh.h":
    cdef cppclass TMREdgeMesh(TMREntity):
        TMREdgeMesh(MPI_Comm, TMREdge*, TMRPoint*, int)
        void mesh(TMRMeshOptions, TMRElementFeatureSize*)
        void writeToVTK(const char*)

cdef extern from "TMRFaceMesh.h":
    cdef cppclass TMRFaceMesh(TMREntity):
        TMRFaceMesh(MPI_Comm, TMRFace*, TMRPoint*, int, int*, int)
        void mesh(TMRMeshOptions, TMRElementFeatureSize*)
        void writeToVTK(const char*)

cdef extern from "TMRVolumeMesh.h":
    cdef cppclass TMRVolumeMesh(TMREntity):
        TMRVolumeMesh(MPI_Comm, TMRVolume*)
        void mesh(TMRMeshOptions)
        void writeToVTK(const char*)

cdef extern from "TMRMesh.h":
    enum TMRFaceMeshType:
        TMR_NO_MESH
        TMR_STRUCTURED
        TMR_UNSTRUCTURED
        TMR_TRIANGLE

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
        void writeToBDF(const char*, int, TMRBoundaryConditions*)

    cdef cppclass TMRMeshOptions:
        TMRMeshOptions()
        TMRFaceMeshType mesh_type_default
        int triangularize_print_level
        int triangularize_print_iter
        int write_mesh_quality_histogram
        int num_smoothing_steps
        double frontal_quality_factor
        int reset_mesh_objects
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
        int32_t face
        int32_t x
        int32_t y
        int32_t tag
        int16_t level
        int16_t info

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
        TMRTopology* getTopology()
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
        TMRInterpolationType getInterpType()
        void setMeshOrder(int, TMRInterpolationType)
        void getNodeConn(const int**, int*)
        int getDepNodeConn(const int**, const int**, const double**)
        TMRQuadrantArray* getQuadsWithName(const char*)
        int getNodesWithName(const char*, int**)
        void createInterpolation(TMRQuadForest*, TACSBVecInterp*)
        int getOwnedNodeRange(const int**)
        void getQuadrants(TMRQuadrantArray**)
        int getPoints(TMRPoint**)
        int getLocalNodeNumber(int);
        int getExtPreOffset()
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
        int32_t block
        int32_t x
        int32_t y
        int32_t z
        int32_t tag
        int16_t level
        int16_t info

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
        TMRTopology* getTopology()
        void setConnectivity(int, const int*, int)
        void setFullConnectivity(int, int, int, const int*, const int*)
        void repartition(int)
        void createTrees(int)
        void createRandomTrees(int, int, int)
        void refine(int*, int, int)
        TMROctForest *duplicate()
        TMROctForest *coarsen()
        void balance(int)
        void createNodes()
        int getMeshOrder()
        TMRInterpolationType getInterpType()
        void setMeshOrder(int, TMRInterpolationType)
        void getNodeConn(const int**, int*)
        int getDepNodeConn(const int**, const int**, const double**)
        TMROctantArray* getOctsWithName(const char*)
        int getNodesWithName(const char*, int**)
        void createInterpolation(TMROctForest*, TACSBVecInterp*)
        int getOwnedNodeRange(const int**)
        void getOctants(TMROctantArray**)
        int getPoints(TMRPoint**)
        int getExtPreOffset()
        void writeToVTK(const char*)
        void writeForestToVTK(const char*)

cdef extern from "TMRBoundaryConditions.h":
    cdef cppclass TMRBoundaryConditions(TMREntity):
        TMRBoundaryConditions()
        void addBoundaryCondition(const char*, int, const int*,
                                  const TacsScalar*)
        int getNumBoundaryConditions()

cdef extern from "TMR_TACSCreator.h":
    cdef cppclass TMRQuadTACSCreator(TMREntity):
        TMRQuadTACSCreator(TMRBoundaryConditions*, int, TMRQuadForest*)
        TMRQuadForest* getFilter()

    cdef cppclass TMROctTACSCreator(TMREntity):
        TMROctTACSCreator(TMRBoundaryConditions*, int, TMROctForest*)
        TMROctForest* getFilter()

cdef extern from "TMROpenCascade.h":
    cdef TMRModel* TMR_LoadModelFromIGESFile(const char*, int)
    cdef TMRModel* TMR_LoadModelFromSTEPFile(const char*, int)

cdef extern from "TMREgads.h" namespace "TMR_EgadsInterface":
    cdef TMRModel* TMR_ConvertEGADSModel"TMR_EgadsInterface::TMR_ConvertEGADSModel"(ego, int)
    cdef TMRModel* TMR_LoadModelFromEGADSFile"TMR_EgadsInterface::TMR_LoadModelFromEGADSFile"(const char*, int)

cdef extern from "TMR_RefinementTools.h":
    void TMR_CreateTACSMg(int, TACSAssembler**,
                          TMRQuadForest**, TACSMg**, double, int, int, int)
    void TMR_ComputeInterpSolution(TMRQuadForest*, TACSAssembler*,
                                   TMRQuadForest*, TACSAssembler*,
                                   TACSBVec*, TACSBVec*)
    void TMR_ComputeReconSolution(TMRQuadForest*, TACSAssembler*,
                                  TMRQuadForest*, TACSAssembler*,
                                  TACSBVec*, TACSBVec*, int)
    double TMR_StrainEnergyErrorEst(TMRQuadForest*, TACSAssembler*,
                                    TMRQuadForest *, TACSAssembler*,
                                    double*)
    double TMR_AdjointErrorEst(TMRQuadForest*, TACSAssembler*,
                               TMRQuadForest*, TACSAssembler*,
                               TACSBVec*, TACSBVec*, double*, double*)

    void TMR_CreateTACSMg(int, TACSAssembler**,
                          TMROctForest**, TACSMg**, double, int, int, int)
    void TMR_ComputeInterpSolution(TMROctForest*, TACSAssembler*,
                                   TMROctForest*, TACSAssembler*,
                                   TACSBVec*, TACSBVec*)
    void TMR_ComputeReconSolution(TMROctForest*, TACSAssembler*,
                                  TMROctForest*, TACSAssembler*,
                                  TACSBVec*, TACSBVec*, int)
    double TMR_StrainEnergyErrorEst(TMROctForest*, TACSAssembler*,
                                    TMROctForest *, TACSAssembler*,
                                    double*)
    double TMR_AdjointErrorEst(TMROctForest*, TACSAssembler*,
                               TMROctForest*, TACSAssembler*,
                               TACSBVec*, TACSBVec*, double*, double*)

cdef extern from "TMRCyCreator.h":
    ctypedef TACSElement* (*createquadelements)(void*, int, TMRQuadrant*)
    ctypedef TACSElement* (*createoctelements)(void*, int, TMROctant*)
    ctypedef TACSElement* (*createquadtopoelements)(
        void*, int, TMRQuadrant*, int, TMRIndexWeight*)
    ctypedef TACSElement* (*createocttopoelements)(
        void*, int, TMROctant*, int, TMRIndexWeight*)

    cdef cppclass TMRCyQuadCreator(TMRQuadTACSCreator):
        TMRCyQuadCreator(TMRBoundaryConditions*, int, TMRQuadForest*)
        void setSelfPointer(void*)
        void setCreateQuadElement(
            TACSElement* (*createquadelements)(void*, int, TMRQuadrant*))
        TACSAssembler *createTACS(TMRQuadForest*, OrderingType)

    cdef cppclass TMRCyOctCreator(TMROctTACSCreator):
        TMRCyOctCreator(TMRBoundaryConditions*, int, TMROctForest*)
        void setSelfPointer(void*)
        void setCreateOctElement(
            TACSElement* (*createoctelements)(void*, int, TMROctant*))
        TACSAssembler *createTACS(TMROctForest*, OrderingType)

    cdef cppclass TMRCyTopoQuadCreator(TMRQuadTACSCreator):
        TMRCyTopoQuadCreator(TMRBoundaryConditions*, int, TMRQuadForest*)
        void setSelfPointer(void*)
        void setCreateQuadTopoElement(
            TACSElement* (*createquadtopoelements)(
                void*, int, TMRQuadrant*, int, TMRIndexWeight*))
        TACSAssembler *createTACS(TMRQuadForest*, OrderingType)

    cdef cppclass TMRCyTopoOctCreator(TMROctTACSCreator):
        TMRCyTopoOctCreator(TMRBoundaryConditions*, int, TMROctForest*)
        void setSelfPointer(void*)
        void setCreateOctTopoElement(
            TACSElement* (*createocttopoelements)(
                void*, int, TMROctant*, int, TMRIndexWeight*))
        TACSAssembler *createTACS(TMROctForest*, OrderingType)

    cdef cppclass TMRCyTopoQuadConformCreator(TMRQuadTACSCreator):
       TMRCyTopoQuadConformCreator(TMRBoundaryConditions*, int, TMRQuadForest*,
                                   int, TMRInterpolationType)
       void setSelfPointer(void*)
       void setCreateQuadTopoElement(
          TACSElement* (*createquadtopoelements)(
             void*, int, TMRQuadrant*, int, const int*, TMRQuadForest*))
       TACSAssembler *createTACS(TMRQuadForest*, OrderingType)

    cdef cppclass TMRCyTopoOctConformCreator(TMROctTACSCreator):
        TMRCyTopoOctConformCreator(TMRBoundaryConditions*, int, TMROctForest*,
                                   int, TMRInterpolationType)
        void setSelfPointer(void*)
        void setCreateOctTopoElement(
            TACSElement* (*createocttopoelements)(
                void*, int, TMROctant*, int, const int*, TMROctForest*))
        TACSAssembler *createTACS(TMROctForest*, OrderingType)

cdef extern from "TMRTopoFilter.h":
    cdef cppclass TMRTopoFilter(TMREntity):
        TACSAssembler* getAssembler()
        TMRQuadForest* getFilterQuadForest()
        TMROctForest* getFilterOctForest()
        void setDesignVars(TACSBVec*)
        void addValues(TACSBVec*)
        void applyFilter(TACSBVec*, TACSBVec*)
        void applyTranspose(TACSBVec*, TACSBVec*)

cdef class TopoFilter:
    cdef TMRTopoFilter *ptr

cdef inline _init_TopoFilter(TMRTopoFilter *ptr):
   fltr = TopoFilter()
   fltr.ptr = ptr
   fltr.ptr.incref()
   return fltr

cdef extern from "TMRLagrangeFilter.h":
    cdef cppclass TMRLagrangeFilter(TMRTopoFilter):
        TMRLagrangeFilter(int, TACSAssembler**, TMROctForest**)
        TMRLagrangeFilter(int, TACSAssembler**, TMRQuadForest**)

cdef extern from "TMRConformFilter.h":
    cdef cppclass TMRConformFilter(TMRTopoFilter):
        TMRConformFilter(int, TACSAssembler**, TMROctForest**)
        TMRConformFilter(int, TACSAssembler**, TMRQuadForest**)

cdef extern from "TMRHelmholtzFilter.h":
    cdef cppclass TMRHelmholtzFilter(TMRTopoFilter):
        TMRHelmholtzFilter(double, int, TACSAssembler**, TMROctForest**)
        TMRHelmholtzFilter(double, int, TACSAssembler**, TMRQuadForest**)

cdef extern from "TMRMatrixFilter.h":
    cdef cppclass TMRMatrixFilter(TMRTopoFilter):
        TMRMatrixFilter(double, int, int, TACSAssembler**, TMROctForest**)
        TMRMatrixFilter(double, int, int, TACSAssembler**, TMRQuadForest**)

cdef extern from "TMROctConstitutive.h":
    cdef enum TMRTopoPenaltyType:
        TMR_SIMP_PENALTY
        TMR_RAMP_PENALTY

    cdef cppclass TMRStiffnessProperties(TMREntity):
        TMRStiffnessProperties(int, TACSMaterialProperties**,
                               double, double, double,
                               TMRTopoPenaltyType, double, double,
                               double, double, double, double, int)
        int nmats
        TACSMaterialProperties **props
        TMRTopoPenaltyType penalty_type
        double stiffness_penalty_value
        double stiffness_offset
        double mass_penalty_value
        double conduction_penalty_value
        double temperature_penalty_value
        double stress_relax_value
        double ks_penalty
        double beta
        double xoffset
        int use_project

    cdef cppclass TMROctConstitutive(TACSSolidConstitutive):
        TMROctConstitutive(TMRStiffnessProperties*, TMROctForest*)

cdef extern from "TMRQuadConstitutive.h":
    cdef cppclass TMRQuadConstitutive(TACSPlaneStressConstitutive):
        TMRQuadConstitutive(TMRStiffnessProperties*, TMRQuadForest*)

cdef extern from "TMRHelmholtzPUFilter.h":
    ctypedef int (*getinteriorstencil)( void*, int, int,
                                        TacsScalar*, double* )
    ctypedef int (*getboundarystencil)( void*, int, TacsScalar*, int,
                                        TacsScalar*, double* )

    cdef cppclass TMRCallbackHelmholtzPUFilter(TMRTopoFilter):
        TMRCallbackHelmholtzPUFilter(int, int, TACSAssembler**,
                                     TMROctForest**)
        TMRCallbackHelmholtzPUFilter(int, int, TACSAssembler**,
                                     TMRQuadForest**)
        void initialize()
        void setSelfPointer(void*)
        void setGetInteriorStencil(getinteriorstencil)
        void setGetBoundaryStencil(getboundarystencil)

cdef extern from "TMRTopoProblem.h":
    ctypedef void (*writeoutputcallback)(void*, const char*, int,
                                         TMROctForest*, TMRQuadForest*,
                                         TACSBVec*)
    ctypedef void (*constraintcallback)(void*, TMRTopoFilter*, TACSMg*,
                                        int, TacsScalar*)
    ctypedef void (*constraintgradientcallback)(void*, TMRTopoFilter*, TACSMg*,
                                                int, TACSBVec**)
    ctypedef void (*objectivecallback)(void*, TMRTopoFilter*, TACSMg*,
                                       TacsScalar*)
    ctypedef void (*objectivegradientcallback)(void*, TMRTopoFilter*, TACSMg*,
                                               TACSBVec*)
    ctypedef void (*qncorrectioncallback)(int, void*, ParOptVec*, ParOptScalar*,
                                          ParOptVec*, ParOptVec*, ParOptVec*)

    cdef cppclass TMRTopoProblem(ParOptProblem):
        TMRTopoProblem(TMRTopoFilter*, TACSMg*, int, double)
        TACSAssembler *getAssembler()
        TMRQuadForest* getFilterQuadForest()
        TMROctForest* getFilterOctForest()
        TMRTopoFilter* getTopoFilter()
        TACSMg* getMg()
        void setLoadCases(TACSBVec**, int)
        int getNumLoadCases()
        void addConstraints(int, TACSFunction**,
                            const TacsScalar*, const TacsScalar*, int)
        void addLinearConstraints(ParOptVec**, TacsScalar*, int)
        void addFrequencyConstraint(double, int, TacsScalar,
                                    TacsScalar, TacsScalar, int,
                                    double, int, int, double, double,
                                    int, JDRecycleType)
        void addBucklingConstraint(double, int, TacsScalar,
                                   TacsScalar, TacsScalar, int, double)
        void addConstraintCallback(int, int, void*,
                                   void (*constraintcallback)(void*, TMRTopoFilter*, TACSMg*,
                                                              int, TacsScalar*),
                                   void*,
                                   void (*constraintgradientcallback)(void*, TMRTopoFilter*, TACSMg*,
                                                                      int, TACSBVec**))
        void addQnCorrectionCallback(int, void*,
                                     void (*qncorrectioncallback)(int, void*, ParOptVec*, ParOptScalar*,
                                           ParOptVec*, ParOptVec*, ParOptVec*))
        void setObjective(const TacsScalar*)
        void setObjective(const TacsScalar*, TACSFunction**)
        void addObjectiveCallback(void*,
                                  void (*objectivecallback)(void*, TMRTopoFilter*, TACSMg*,
                                                            TacsScalar*),
                                  void*,
                                  void (*objectivegradientcallback)(void*, TMRTopoFilter*, TACSMg*,
                                                                    TACSBVec*))
        void initialize()
        void setPrefix(const char*)
        void setInitDesignVars(ParOptVec*,ParOptVec*,ParOptVec*)
        void setIterationCounter(int)
        ParOptVec* createDesignVec()
        void setF5OutputFlags(int, ElementType, int)
        void setF5EigenOutputFlags(int, ElementType, int)
        void setOutputCallback(void*,
            void (*writeoutputcallback)(void*, const char*, int,
                                        TMROctForest*, TMRQuadForest*,
                                        TACSBVec*))
        void useQnCorrectionComplianceObj()
        void addNonDesignMass(ParOptVec*)
        int evalObjCon(ParOptVec*, ParOptScalar*, ParOptScalar*)
        int evalObjConGradient(ParOptVec*, ParOptVec*, ParOptVec**)

    cdef cppclass ParOptBVecWrap(ParOptVec):
        ParOptBVecWrap(TACSBVec*)
        TACSBVec *vec

cdef extern from "TMR_STLTools.h":
    int TMR_GenerateBinFile(const char*, TMROctForest*,
                            TACSBVec*, int, double)
    int TMR_GenerateSTLTriangles(int, TMROctForest*, TACSBVec*,
                                 int, double, int*, TMR_STLTriangle**)

cdef extern from "TMRApproximateDistance.h":
    void TMRApproximateDistance(TMRQuadForest*, int, double, double,
                                TACSBVec*, const char*, double*)
    void TMRApproximateDistance(TMROctForest*, int, double, double,
                                TACSBVec*, const char*, double*)
