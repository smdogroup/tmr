/*
  This file is part of the package TMR for adaptive mesh refinement.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#ifndef TMR_FACE_MESH_H
#define TMR_FACE_MESH_H

#include <map>
#include "TMRMesh.h"
#include "TMREdgeMesh.h"

/*
  The TMRFaceMesh class: This is used to generate quadrilateral and
  triangular face meshes
*/
class TMRFaceMesh : public TMREntity {
 public:
  TMRFaceMesh( MPI_Comm _comm, TMRFace *face,
               TMRPoint *_X=NULL, int _npts=0,
               const int *_quads=NULL, int _nquads=0 );
  ~TMRFaceMesh();

  // Retrieve the underlying surface
  void getFace( TMRFace **_surface );

  // Mesh the underlying geometric object
  void mesh( TMRMeshOptions options,
             TMRElementFeatureSize *fs );

  // Return the type of the underlying mesh
  TMRFaceMeshType getMeshType(){
    return mesh_type;
  }

  // Retrieve the mesh points
  void getMeshPoints( int *_npts, const double **_pts, TMRPoint **X );

  // Order the mesh points uniquely
  int setNodeNums( int *num );
  int getNodeNums( const int **_vars );
  int getNumFixedPoints();

  // Get points indexed via index on an edge or structured face
  void getSourceToTargetMapping( const int **_source_to_target );
  int getFaceIndexFromEdge( TMREdge *e, int idx );
  int getStructuredFaceIndex( TMREdge *e1, int idx1,
                              TMREdge *e2, int idx2 );

  // Retrieve the local connectivity from this surface mesh
  int getQuadConnectivity( const int **_quads );
  int getTriConnectivity( const int **_tris );

  // Write the quadrilateral mesh to a VTK format
  void writeToVTK( const char *filename );

  // Print the quadrilateral quality
  void addMeshQuality( int nbins, int count[] );
  void printMeshQuality();

 private:
  // Set the prescribed mesh
  void setPrescribedMesh( const TMRPoint *_X, int _npts,
                          const int *_quads, int _nquads );

  // Write the segments to a VTK file in parameter space
  void writeSegmentsToVTK( const char *filename,
                           int npts, const double *params,
                           int nsegs, const int segs[] );

  // Print the triangle quality
  void printTriQuality( int ntris, const int tris[] );
  void writeTrisToVTK( const char *filename,
                       int ntris, const int tris[] );

  // Write the dual mesh - used for recombination - to a file
  void writeDualToVTK( const char *filename, int nodes_per_elem,
                       int nelems, const int elems[],
                       int num_dual_edges, const int dual_edges[],
                       const TMRPoint *p );

  // Recombine the mesh to a quadrilateral mesh based on the
  // quad-Blossom algorithm
  void recombine( int ntris, const int tris[], const int tri_neighbors[],
                  const int node_to_tri_ptr[], const int node_to_tris[],
                  int num_edges, const int dual_edges[],
                  int *_num_quads, int **_new_quads,
                  TMRMeshOptions options );

  // Simplify the quadrilateral mesh to remove points that make
  // for a poor quadrilateral mesh
  void simplifyQuads( int flag );

  // Compute recombined triangle information
  int getRecombinedQuad( const int tris[], const int trineighbors[],
                         int t1, int t2, int quad[] );
  double computeRecombinedQuality( const int tris[],
                                   const int trineighbors[],
                                   int t1, int t2, const TMRPoint *p );

  // Compute the quadrilateral quality
  double computeQuadQuality( const int *quad, const TMRPoint *p );
  double computeTriQuality( const int *tri, const TMRPoint *p );

  // Compute the locations of the hole points
  void computeHolePts( const int nloops, int hole_pt,
                       const int *loop_pt_offset,
                       const int *segments, double *params );

  // Map the source face to the target face
  void mapSourceToTarget( TMRMeshOptions options, const double *params );

  // Map the copy source face to the target face
  int mapCopyToTarget( TMRMeshOptions options, const double *params );

  // Create a structured mesh
  void createStructuredMesh( TMRMeshOptions options, const double *params );

  // Create an unstructured mesh
  void createUnstructuredMesh( TMRMeshOptions options,
                               TMRElementFeatureSize *fs,
                               TMRFaceMeshType mesh_type,
                               const int total_num_pts, const int nholes,
                               const double *params,
                               const int nsegs, const int *segments,
                               const int num_degen, const int *degen,
                               int *npts, double **param_pts, TMRPoint **Xpts,
                               int *nquads, int **mesh_quads,
                               int *ntris, int **mesh_tris );

  // The underlying surface
  MPI_Comm comm;
  TMRFace *face;

  // The actual type of mesh used to mesh the structure
  TMRFaceMeshType mesh_type;

  // Points in the mesh
  int num_fixed_pts; // number of fixed points on the boundary
  int num_points; // The number of point locations
  double *pts; // The parametric node locations
  TMRPoint *X; // The physical node locations
  int *vars; // The global variable numbers

  // Source to target mapping - exists only when source/target faces
  // are set in the TMRFace class
  int *source_to_target; // Source to target mapping

  // Copy to target mapping - exists only when source/target indices
  // are set in the TMRFace class
  int *copy_to_target;

  // Record whether this mesh is prescribed
  int prescribed_mesh;

  // Quadrilateral surface mesh information
  int num_quads;
  int *quads;

  // Triangle mesh surface information
  int num_tris;
  int *tris;
};

#endif // TMR_FACE_MESH_H
