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

#include "TMRTopology.h"
#include "TMRFeatureSize.h"
#include "TMRBoundaryConditions.h"

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
void TMR_ComputePlanarTriEdges( int nnodes, int ntris,
                                const int tris[],
                                int *num_tri_edges,
                                int **_tri_edges,
                                int **_tri_neighbors,
                                int **_dual_edges,
                                int **_node_to_tri_ptr=NULL,
                                int **_node_to_tris=NULL,
                                int **_tri_edge_nums=NULL );
void TMR_ComputePlanarQuadEdges( int nnodes, int nquads,
                                 const int quads[],
                                 int *_num_quad_edges,
                                 int **_quad_edges,
                                 int **_quad_neighbors=NULL,
                                 int **_dual_edges=NULL,
                                 int **_quad_edge_nums=NULL );
void TMR_ComputeQuadEdges( int nnodes, int nquads,
                           const int quads[],
                           int *_num_quad_edges,
                           int **_quad_edges );
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
                   int flag=(TMRMesh::TMR_QUAD | TMRMesh::TMR_HEX),
                   TMRBoundaryConditions *bcs=NULL );

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
