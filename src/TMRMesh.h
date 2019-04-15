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

#ifndef TMR_MESH_H
#define TMR_MESH_H

#include <map>
#include "TMRBase.h"
#include "TMRGeometry.h"
#include "TMRTopology.h"
#include "TMRFeatureSize.h"

/*
  The type of face mesh algorithm to apply
*/
enum TMRFaceMeshType { TMR_NO_MESH,
                       TMR_STRUCTURED,
                       TMR_UNSTRUCTURED,
                       TMR_TRIANGLE };

/*
  Methods for computing the mesh connectivity and dual connectivity
  information from other data.
*/
void TMR_ComputeNodeToElems( int nnodes, int nelems,
                             int num_elem_nodes,
                             const int conn[],
                             int **_ptr,
                             int **_node_to_elems );
void TMR_ComputeTriEdges( int nnodes, int ntris,
                          const int tris[],
                          int *num_tri_edges,
                          int **_tri_edges,
                          int **_tri_neighbors,
                          int **_dual_edges,
                          int **_node_to_tri_ptr=NULL,
                          int **_node_to_tris=NULL,
                          int **_tri_edge_nums=NULL );
void TMR_ComputeQuadEdges( int nnodes, int nquads,
                           const int quads[],
                           int *_num_quad_edges,
                           int **_quad_edges,
                           int **_quad_neighbors=NULL,
                           int **_dual_edges=NULL,
                           int **_quad_edge_nums=NULL );
void TMR_ComputeHexEdgesAndFaces( int nnodes, int nhex,
                                  const int hex[],
                                  int *_num_hex_edges,
                                  int **_hex_edges,
                                  int **_hex_edge_nums,
                                  int *_num_hex_faces,
                                  int **_hex_faces,
                                  int **_hex_face_nums );
/*
  Global options for meshing
*/
class TMRMeshOptions {
 public:
  enum TriangleSmoothingType { TMR_LAPLACIAN, TMR_SPRING };

  /*
    Create the mesh options objects with the default settings
  */
  TMRMeshOptions(){
    // Set the default print level
    triangularize_print_level = 0;
    triangularize_print_iter = 1000;
    write_mesh_quality_histogram = 0;

    // Set the default meshing options
    mesh_type_default = TMR_STRUCTURED;
    num_smoothing_steps = 10;
    tri_smoothing_type = TMR_LAPLACIAN;
    frontal_quality_factor = 1.5;

    // By default, reset the mesh objects
    reset_mesh_objects = 1;

    // By default, write nothing to any files
    write_init_domain_triangle = 0;
    write_triangularize_intermediate = 0;
    write_pre_smooth_triangle = 0;
    write_post_smooth_triangle = 0;
    write_dual_recombine = 0;
    write_pre_smooth_quad = 0;
    write_post_smooth_quad = 0;
    write_quad_dual = 0;
  }

  // Set the print level for the triangularize code
  int triangularize_print_level;
  int triangularize_print_iter;

  // Set the write level for the quality histogram
  int write_mesh_quality_histogram;

  // Options to control the meshing algorithm
  TMRFaceMeshType mesh_type_default;
  int num_smoothing_steps;
  TriangleSmoothingType tri_smoothing_type;
  double frontal_quality_factor;

  // Reset the mesh objects in each geometry object
  int reset_mesh_objects;

  // Write intermediate surface meshes to file
  int write_init_domain_triangle;
  int write_triangularize_intermediate;
  int write_pre_smooth_triangle;
  int write_post_smooth_triangle;
  int write_dual_recombine;
  int write_pre_smooth_quad;
  int write_post_smooth_quad;
  int write_quad_dual;
};

/*
  The mesh for a geometric curve
*/
class TMREdgeMesh : public TMREntity {
 public:
  TMREdgeMesh( MPI_Comm _comm, TMREdge *edge,
               TMRPoint *_X=NULL, int _npts=0 );
  ~TMREdgeMesh();

  // Is this edge mesh degenerate
  int isDegenerate(){ return edge->isDegenerate(); }

  // Retrieve the underlying curve
  void getEdge( TMREdge **_edge );

  // Mesh the geometric object
  void mesh( TMRMeshOptions options,
             TMRElementFeatureSize *fs );

  // Order the mesh points uniquely
  int setNodeNums( int *num );
  int getNodeNums( const int **_vars );

  // Retrieve the mesh points
  void getMeshPoints( int *_npts, const double **_pts, TMRPoint **X );

 private:
  MPI_Comm comm;
  TMREdge *edge;

  // Keep track if this mesh is prescribed
  int prescribed_mesh;

  // The parametric locations of the points that are obtained from
  // meshing the curve
  int npts; // number of points along the curve
  double *pts; // Parametric node locations
  TMRPoint *X; // Physical node locations
  int *vars; // Global node variable numbers
};

/*
  Triangle info class

  This class is used to generate a triangular mesh
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

  // Set the mapping of the source/copy to the target (this) mesh
  void setMeshFromMapping( TMRMeshOptions options, const double *params,
                           TMRFace *src, const int rel_orient,
                           std::map<TMREdge*, TMREdge*> &src_to_target_edge,
                           std::map<TMREdge*, int> &src_to_target_orient,
                           int **_src_to_target );

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

/*
  TMRVolumeMesh class

  This is the class that contains the volume mesh. This mesh
  takes in the arguments needed to form a volume mesh
*/
class TMRVolumeMesh : public TMREntity {
 public:
  TMRVolumeMesh( MPI_Comm _comm, TMRVolume *volume );
  ~TMRVolumeMesh();

  // Create the volume mesh
  int mesh( TMRMeshOptions options );

  // Retrieve the mesh points
  void getMeshPoints( int *_npts, TMRPoint **X );

  // Retrieve the local connectivity from this volume mesh
  int getHexConnectivity( const int **hex );
  int getTetConnectivity( const int **tets );

  // Order the mesh points uniquely
  int setNodeNums( int *num );
  int getNodeNums( const int **_vars );

  // Write the volume mesh to a VTK file
  void writeToVTK( const char *filename );

 private:
  // Create a tetrahedral mesh (if possible)
  int tetMesh( TMRMeshOptions options );

  // The underlying volume
  MPI_Comm comm;
  TMRVolume *volume;

  // Set the number of face loops
  int num_swept_faces;
  TMREdge **swept_edges;
  TMRFace **swept_faces;

  // Number of points swept through-thickness
  int num_swept_pts;

  // Keep the bottom/top surfaces (master/target) in the mesh for
  // future reference
  TMRFace *target, *source;
  int target_orient, source_orient;

  int num_points; // The number of points
  TMRPoint *X; // The physical node locations
  int *vars; // The global variable numbers

  // Hexahedral mesh information
  int num_hex;
  int *hex;

  // Tetrahedral mesh information
  int num_tet;
  int *tet;
};

/*
  Mesh the geometry model.

  This class handles the meshing for surface objects without any
  additional information. For hexahedral meshes, the model must
*/
class TMRMesh : public TMREntity {
 public:
  static const int TMR_QUAD = 1;
  static const int TMR_HEX = 2;

  TMRMesh( MPI_Comm _comm, TMRModel *_geo );
  ~TMRMesh();

  // Mesh the underlying geometry
  void mesh( TMRMeshOptions options, double htarget );
  void mesh( TMRMeshOptions options,
             TMRElementFeatureSize *fs );

  // Write the mesh to a VTK file
  void writeToVTK( const char *filename,
                   int flag=(TMRMesh::TMR_QUAD | TMRMesh::TMR_HEX) );

  // Write the mesh to a BDF file
  void writeToBDF( const char *filename,
                   int flag=(TMRMesh::TMR_QUAD | TMRMesh::TMR_HEX) );

  // Retrieve the mesh components
  int getMeshPoints( TMRPoint **_X );
  void getQuadConnectivity( int *_nquads, const int **_quads );
  void getTriConnectivity( int *_ntris, const int **_tris );
  void getHexConnectivity( int *_nhex, const int **_hex );

  // Create a topology object (with underlying mesh geometry)
  TMRModel* createModelFromMesh();

 private:
  // Allocate and initialize the underlying mesh
  void initMesh( int count_nodes=0 );

  // Reset the mesh
  void resetMesh();

  // The underlying geometry object
  MPI_Comm comm;
  TMRModel *geo;

  // The number of nodes/positions in the mesh
  int num_nodes;
  TMRPoint *X;

  // The number of quads
  int num_quads;
  int *quads;

  // The number of triangles
  int num_tris;
  int *tris;

  // The number of hexahedral elements in the mesh
  int num_hex;
  int *hex;

  // The number of tetrahedral elements
  int num_tet;
  int *tet;
};

#endif // TMR_MESH_H
