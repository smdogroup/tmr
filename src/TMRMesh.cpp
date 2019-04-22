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

#include <math.h>
#include <stdio.h>

#include "TMRMesh.h"
#include "TMREdgeMesh.h"
#include "TMRFaceMesh.h"
#include "TMRVolumeMesh.h"
#include "TMRNativeTopology.h"
#include "TMRTriangularize.h"
#include "TMRMeshSmoothing.h"
#include "TMRBspline.h"
#include "tmrlapack.h"

/*
  The triangle nodes and edges are ordered locally as follows. Note
  that the edges are ordered based on the node across the triangle.

       2
      / \
   e1/   \e0
    /     \
   /       \
  0 ------- 1
       e2

  The edge nodes for a given edge index in a triangle
*/
const int tri_edge_nodes[][2] = {{1, 2}, {2, 0}, {0, 1}};

/*
  The edges that are connected to a given node
*/
const int tri_node_edges[][2] = {{1, 2}, {0, 2}, {0, 1}};


/*
  The nodes in a quadrilateral element are ordered locally as follows:

  3 --------- 2
  |           |
  |           |
  |           |
  |           |
  0 --------- 1
*/

/*
  The local node numbers for each edge in the quadrilateral
*/
const int quad_edge_nodes[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

/*
  The local node numbers for each edge in a flipped quadrilateral
*/
const int flipped_quad_edge_nodes[][2] = {{0, 3}, {3, 2}, {2, 1}, {1, 0}};

/*
  The nodes in a hexahedral element are ordered as follows:

       7 ---------- 6
      /|           /|
     / |          / |
    /  |         /  |
   4 ---------- 5   |
   |   |        |   |
   |   3 -------|-- 2
   |  /         |  /
   | /          | /
   |/           |/
   0 ---------- 1
*/

/*
  The start/end nodes corresponding to each edge in a
  hexahedral elements
*/
const int hex_edge_nodes[12][2] =
  {{0,1}, {3,2}, {4,5}, {7,6},  // x-aligned edges
   {0,3}, {1,2}, {4,7}, {5,6},  // y-aligned edges
   {0,4}, {1,5}, {3,7}, {2,6}}; // z-aligned edges

/*
  The nodes corresponding to each face in a hexahedral
  element
*/
const int hex_face_nodes[6][4] =
  {{0,3,7,4}, {1,2,6,5},  // x-faces
   {0,4,5,1}, {3,7,6,2},  // y-faces
   {0,1,2,3}, {4,5,6,7}}; // z-faces

/*
  Possible orientations for two connecting faces
*/
const int face_orientations[8][4] =
  {{0,1,2,3}, {1,2,3,0},
   {2,3,0,1}, {3,0,1,2},
   {3,2,1,0}, {0,3,2,1},
   {1,0,3,2}, {2,1,0,3}};

/*
  Hexahedral local node to local coordinate ordering transformation
*/
const int hex_coordinate_order[] =
  {0, 1, 3, 2, 4, 5, 7, 6};

/*
  Compare coordinate pairs of points. This uses a Morton ordering
  comparison to facilitate sorting/searching the list of edges.

  This function can be used by the stdlib functions qsort and
  bsearch to sort/search the edge pairs.
*/
static int compare_edges( const void *avoid, const void *bvoid ){
  const int *a = static_cast<const int*>(avoid);
  const int *b = static_cast<const int*>(bvoid);

  // Extract the x/y locations for the a and b points
  int ax = a[1], ay = a[2];
  int bx = b[1], by = b[2];

  int xxor = ax ^ bx;
  int yxor = ay ^ by;
  int sor = xxor | yxor;

  // Note that here we do not distinguish between levels
  // Check for the most-significant bit
  int discrim = 0;
  if (xxor > (sor ^ xxor)){
    discrim = ax - bx;
  }
  else {
    discrim = ay - by;
  }

  return discrim;
}

/*
  Compare faces. The faces themselves must be sorted such that
  the indices are increasing.

  This is used to determine the resulting face number given the
  corresponding list of nodes.
*/
static int compare_faces( const void *avoid, const void *bvoid ){
  const int *a = static_cast<const int*>(avoid);
  const int *b = static_cast<const int*>(bvoid);

  if (a[1] == b[1]){
    if (a[2] == b[2]){
      if (a[3] == b[3]){
        return a[4] - b[4];
      }
      else {
        return a[3] - b[3];
      }
    }
    else {
      return a[2] - b[2];
    }
  }
  else {
    return a[1] - b[1];
  }

  return 0;
}

/*
  Sort the nodes witin a face array
*/
static inline void sort_face_nodes( int *node ){
  int tmp[4];
  if (node[0] < node[1]){
    tmp[0] = node[0];
    tmp[1] = node[1];
  }
  else {
    tmp[0] = node[1];
    tmp[1] = node[0];
  }

  if (node[2] < node[3]){
    tmp[2] = node[2];
    tmp[3] = node[3];
  }
  else {
    tmp[2] = node[3];
    tmp[3] = node[2];
  }

  // Merge the two arrays
  int i = 0, j = 2;
  for ( int k = 0; k < 4; k++ ){
    if (i >= 2){
      node[k] = tmp[j];
      j++;
    }
    else if (j >= 4){
      node[k] = tmp[i];
      i++;
    }
    else if (tmp[i] <= tmp[j]){
      node[k] = tmp[i];
      i++;
    }
    else {
      node[k] = tmp[j];
      j++;
    }
  }
}
/*
  Compute a node to triangle or node to quad data structure
*/
void TMR_ComputeNodeToElems( int nnodes, int nelems,
                             int num_elem_nodes,
                             const int conn[],
                             int **_ptr,
                             int **_node_to_elems ){
  // Set the pointer
  int *ptr = new int[ nnodes+1 ];
  memset(ptr, 0, (nnodes+1)*sizeof(int));

  // Count up the references
  const int conn_size = nelems*num_elem_nodes;
  for ( int i = 0; i < conn_size; i++ ){
    if (conn[i] >= 0){
      ptr[conn[i]+1]++;
    }
  }

  // Set the pointer into the quad array
  for ( int i = 0; i < nnodes; i++ ){
    ptr[i+1] += ptr[i];
  }

  // Compute the node to quads
  int *node_to_elems = new int[ ptr[nnodes] ];
  const int *conn_ptr = conn;
  for ( int i = 0; i < nelems; i++ ){
    for ( int j = 0; j < num_elem_nodes; j++ ){
      int node = conn_ptr[0];
      if (node >= 0){
        node_to_elems[ptr[node]] = i;
        ptr[node]++;
        conn_ptr++;
      }
    }
  }

  // Reset the pointer array
  for ( int i = nnodes-1; i >= 0; i-- ){
    ptr[i+1] = ptr[i];
  }
  ptr[0] = 0;

  // Set the output points
  *_ptr = ptr;
  *_node_to_elems = node_to_elems;
}

/*
  Compute all of the edges within the triangular mesh
*/
void TMR_ComputeTriEdges( int nnodes, int ntris,
                          const int tris[],
                          int *num_tri_edges,
                          int **_tri_edges,
                          int **_tri_neighbors,
                          int **_dual_edges,
                          int **_node_to_tri_ptr,
                          int **_node_to_tris,
                          int **_tri_edge_nums ){
  // Compute the edges in the triangular mesh
  int *ptr;
  int *node_to_tris;
  TMR_ComputeNodeToElems(nnodes, ntris, 3, tris, &ptr, &node_to_tris);

  // Now compute the neighbors for each triangle
  int *tri_edge_nums = new int[ 3*ntris ];
  for ( int i = 0; i < 3*ntris; i++ ){
    tri_edge_nums[i] = -1;
  }

  // Allocate the array for the triangle neighbors
  int *tri_neighbors = new int[ 3*ntris ];

  int ne = 0;
  for ( int i = 0; i < ntris; i++ ){
    // Search through each edge of the each triangle
    for ( int j = 0; j < 3; j++ ){
      if (tri_edge_nums[3*i+j] < 0){
        tri_edge_nums[3*i+j] = ne;

        // Triangle edges that have no neighbors are labeled with a -1
        tri_neighbors[3*i+j] = -1;

        int e0[2];
        e0[0] = tris[3*i + tri_edge_nodes[j][0]];
        e0[1] = tris[3*i + tri_edge_nodes[j][1]];

        // Search for the neighboring quad that shares this edge
        int kp = ptr[e0[0]];
        int kpend = ptr[e0[0]+1];
        for ( ; kp < kpend; kp++ ){
          // Find the potential quad neighbor
          int n = node_to_tris[kp];

          // Don't count the same edge twice
          if (n == i){ continue; }

          // Flag to indicate that we have found the other edge (there
          // will only be at most one other match since this is
          // planar in parameter space)
          int quit = 0;

          // Search over all the edges on this quad, and see
          // if they match
          for ( int e = 0; e < 3; e++ ){
            int e1[2];
            e1[0] = tris[3*n + tri_edge_nodes[e][0]];
            e1[1] = tris[3*n + tri_edge_nodes[e][1]];

            // Check if the adjacent edge matches in either direction
            if ((e0[0] == e1[0] && e0[1] == e1[1]) ||
                (e0[0] == e1[1] && e0[1] == e1[0])){
              // Label the other edge that shares this same node
              tri_edge_nums[3*n+e] = ne;

              // Set the triangle neighbors
              tri_neighbors[3*n+e] = i;
              tri_neighbors[3*i+j] = n;

              quit = 1;
            }
          }
          if (quit){ break; }
        }

        // Increment the edge number
        ne++;
      }
    }
  }

  // Free the data that is no longer required
  if (_node_to_tri_ptr){
    *_node_to_tri_ptr = ptr;
  }
  else {
    delete [] ptr;
  }
  if (_node_to_tris){
    *_node_to_tris = node_to_tris;
  }
  else {
    delete [] node_to_tris;
  }

  // Now we have a unique list of edge numbers and the total number of
  // edges, we can construct the unique edge list
  int *tri_edges = new int[ 2*ne ];
  int *dual_edges = new int[ 2*ne ];
  for ( int i = 0; i < ntris; i++ ){
    for ( int j = 0; j < 3; j++ ){
      // Get the unique edge number
      int n = tri_edge_nums[3*i+j];
      tri_edges[2*n] = tris[3*i + tri_edge_nodes[j][0]];
      tri_edges[2*n+1] = tris[3*i + tri_edge_nodes[j][1]];

      // Set the dual edge numbers - connecting triangles to other
      // triangles. Note that some elements of this array will be -1.
      dual_edges[2*n] = i;
      dual_edges[2*n+1] = tri_neighbors[3*i+j];
    }
  }

  // Free the edge numbers or delete them
  if (_tri_edge_nums){
    *_tri_edge_nums = tri_edge_nums;
  }
  else {
    delete [] tri_edge_nums;
  }

  // Set the number of triangle edges and the triangle edges themselves
  *num_tri_edges = ne;
  *_tri_edges = tri_edges;
  *_tri_neighbors = tri_neighbors;
  *_dual_edges = dual_edges;
}

/*
  Compute the connectivity between the edges
*/
void TMR_ComputeQuadEdges( int nnodes, int nquads,
                           const int quads[],
                           int *_num_quad_edges,
                           int **_quad_edges,
                           int **_quad_neighbors,
                           int **_dual_edges,
                           int **_quad_edge_nums ){
  // Compute the connectivity from nodes to quads
  int *ptr;
  int *node_to_quads;
  TMR_ComputeNodeToElems(nnodes, nquads, 4, quads, &ptr, &node_to_quads);

  // Now compute the neighbors for each quad
  int *quad_edge_nums = new int[ 4*nquads ];
  for ( int i = 0; i < 4*nquads; i++ ){
    quad_edge_nums[i] = -1;
  }

  // Allocate an array for the quad neighbors
  int *quad_neighbors = NULL;
  if (_quad_neighbors){
    quad_neighbors = new int[ 4*nquads ];
  }

  int ne = 0;
  for ( int i = 0; i < nquads; i++ ){
    // Search through each edge of the each quadrilateral
    for ( int j = 0; j < 4; j++ ){
      if (quad_edge_nums[4*i+j] < 0){
        quad_edge_nums[4*i+j] = ne;

        // Set the quad neighbor index to -1
        if (quad_neighbors){
          quad_neighbors[4*i+j] = -1;
        }

        // Extract the edge nodes for the j-th edge
        int e0[2];
        e0[0] = quads[4*i + quad_edge_nodes[j][0]];
        e0[1] = quads[4*i + quad_edge_nodes[j][1]];

        // Search for the neighboring quad that shares this edge
        int kp = ptr[e0[0]];
        int kpend = ptr[e0[0]+1];

        for ( ; kp < kpend; kp++ ){
          // Find the potential quad neighbor
          int n = node_to_quads[kp];

          // Don't count the same edge twice
          if (n == i){ continue; }

          // Flag to indicate that we have found the other edge (there
          // will only be at most one other match since this is planar)
          int quit = 0;

          // Search over all the edges on this quad, and see
          // if they match
          for ( int e = 0; e < 4; e++ ){
            // Extract the edge nodes for the e-th edge
            int e1[2];
            e1[0] = quads[4*n + quad_edge_nodes[e][0]];
            e1[1] = quads[4*n + quad_edge_nodes[e][1]];

            // Check if the adjacent edge matches in either direction
            if ((e0[0] == e1[0] && e0[1] == e1[1]) ||
                (e0[0] == e1[1] && e0[1] == e1[0])){
              // Label the other edge that shares this same node
              quad_edge_nums[4*n+e] = ne;

              // Set the quad neighbors
              if (quad_neighbors){
                quad_neighbors[4*n+e] = i;
                quad_neighbors[4*i+j] = n;
              }

              quit = 1;
            }
          }
          if (quit){ break; }
        }

        // Increment the edge number
        ne++;
      }
    }
  }

  // Free the pointers
  delete [] ptr;
  delete [] node_to_quads;

  // Now we have a unique list of edge numbers and the total number of
  // edges, we can construct the unique edge list
  int *quad_edges = new int[ 2*ne ];
  int *dual_edges = NULL;
  if (_dual_edges){
    dual_edges = new int[ 2*ne ];
  }
  for ( int i = 0; i < nquads; i++ ){
    for ( int j = 0; j < 4; j++ ){
      // Get the unique edge number
      int n = quad_edge_nums[4*i+j];
      quad_edges[2*n] = quads[4*i + quad_edge_nodes[j][0]];
      quad_edges[2*n+1] = quads[4*i + quad_edge_nodes[j][1]];

      // Set the dual edges
      if (dual_edges){
        dual_edges[2*n] = i;
        dual_edges[2*n+1] = quad_neighbors[4*i+j];
      }
    }
  }

  // Free the data
  if (_quad_edge_nums){
    *_quad_edge_nums = quad_edge_nums;
  }
  else {
    delete [] quad_edge_nums;
  }

  *_num_quad_edges = ne;
  *_quad_edges = quad_edges;
  if (_quad_neighbors){
    *_quad_neighbors = quad_neighbors;
  }
  if (_dual_edges){
    *_dual_edges = dual_edges;
  }
}

/*
  Compute the connectivity between edges within a hexahedral mesh
*/
void TMR_ComputeHexEdgesAndFaces( int nnodes, int nhex,
                                  const int hex[],
                                  int *_num_hex_edges,
                                  int **_hex_edges,
                                  int **_hex_edge_nums,
                                  int *_num_hex_faces,
                                  int **_hex_faces,
                                  int **_hex_face_nums ){
  // Compute the connectivity from nodes to hex
  int *ptr;
  int *node_to_hex;
  TMR_ComputeNodeToElems(nnodes, nhex, 8, hex, &ptr, &node_to_hex);

  // Comput the neighbors for each hex
  int *hex_edge_nums = new int[ 12*nhex ];
  int *hex_face_nums = new int[ 6*nhex ];
  for ( int i = 0; i < 12*nhex; i++ ){
    hex_edge_nums[i] = -1;
  }
  for ( int i = 0; i < 6*nhex; i++ ){
    hex_face_nums[i] = -1;
  }

  int edge_num = 0, face_num = 0;
  for ( int i = 0; i < nhex; i++ ){
    // Search through each hexahedral element for the edge
    for ( int j = 0; j < 12; j++ ){
      if (hex_edge_nums[12*i+j] < 0){
        hex_edge_nums[12*i+j] = edge_num;

        // Find the edge nodes for the new edge
        int e[2];
        e[0] = hex[8*i + hex_edge_nodes[j][0]];
        e[1] = hex[8*i + hex_edge_nodes[j][1]];

        // Find the node that shares the edge
        int kp = ptr[e[0]];
        int kpend =ptr[e[0]+1];

        for ( ; kp < kpend; kp++ ){
          // Find the potential hex neighbor
          int ihex = node_to_hex[kp];
          if (i == ihex){
            continue;
          }

          // Search over all the edges on this hex, and see if any of
          // the edges match
          for ( int k = 0; k < 12; k++ ){
            int e2[2];
            e2[0] = hex[8*ihex + hex_edge_nodes[k][0]];
            e2[1] = hex[8*ihex + hex_edge_nodes[k][1]];

            // Check if the adjacent edge matches in either direction
            if ((e[0] == e2[0] && e[1] == e2[1]) ||
                (e[0] == e2[1] && e[1] == e2[0])){
              hex_edge_nums[12*ihex + k] = edge_num;
            }
          }
        }

        // Increment the edge number
        edge_num++;
      }
    }

    // Search through adjacent hexahedral elements for the faces
    for ( int j = 0; j < 6; j++ ){
      if (hex_face_nums[6*i+j] < 0){
        hex_face_nums[6*i+j] = face_num;

        // Find the nodes associated with face j
        int f[4];
        f[0] = hex[8*i + hex_face_nodes[j][0]];
        f[1] = hex[8*i + hex_face_nodes[j][1]];
        f[2] = hex[8*i + hex_face_nodes[j][2]];
        f[3] = hex[8*i + hex_face_nodes[j][3]];

        // Find a node that shares the face
        int kp = ptr[f[0]];
        int kpend = ptr[f[0]+1];

        for ( ; kp < kpend; kp++ ){
          // Find the potential hex neighbor
          int ihex = node_to_hex[kp];
          if (i == ihex){
            continue;
          }

          // Search over all the faces on this hex, and see if any of
          // the edges match
          for ( int k = 0; k < 6; k++ ){
            int f2[4];
            f2[0] = hex[8*ihex + hex_face_nodes[k][0]];
            f2[1] = hex[8*ihex + hex_face_nodes[k][1]];
            f2[2] = hex[8*ihex + hex_face_nodes[k][2]];
            f2[3] = hex[8*ihex + hex_face_nodes[k][3]];

            // Check if the adjacent edge matches in either direction
            for ( int ort = 0; ort < 8; ort++ ){
              if (f[0] == f2[face_orientations[ort][0]] &&
                  f[1] == f2[face_orientations[ort][1]] &&
                  f[2] == f2[face_orientations[ort][2]] &&
                  f[3] == f2[face_orientations[ort][3]]){
                hex_face_nums[6*ihex + k] = face_num;
                break;
              }
            }
          }
        }

        // Increment the face number
        face_num++;
      }
    }
  }

  // Free the node->hex data
  delete [] ptr;
  delete [] node_to_hex;

  if (_num_hex_edges){
    *_num_hex_edges = edge_num;
  }
  if (_num_hex_faces){
    *_num_hex_faces = face_num;
  }

  // Create the hex edges/faces
  int *hex_edges = new int[ 2*edge_num ];
  int *hex_faces = new int[ 4*face_num ];
  for ( int i = 0; i < nhex; i++ ){
    for ( int j = 0; j < 12; j++ ){
      int e = hex_edge_nums[12*i+j];
      hex_edges[2*e] = hex[8*i + hex_edge_nodes[j][0]];
      hex_edges[2*e+1] = hex[8*i + hex_edge_nodes[j][1]];
    }

    for ( int j = 0; j < 6; j++ ){
      int f = hex_face_nums[6*i+j];
      hex_faces[4*f] = hex[8*i + hex_face_nodes[j][0]];
      hex_faces[4*f+1] = hex[8*i + hex_face_nodes[j][1]];
      hex_faces[4*f+2] = hex[8*i + hex_face_nodes[j][2]];
      hex_faces[4*f+3] = hex[8*i + hex_face_nodes[j][3]];
    }
  }

  // Set the pointers to the data that may (or may not)
  // be returned from this function call.
  if (_hex_edges){
    *_hex_edges = hex_edges;
  }
  else {
    delete [] hex_edges;
  }
  if (_hex_faces){
    *_hex_faces = hex_faces;
  }
  else {
    delete [] hex_faces;
  }

  // Free the allocated arrays
  if (_hex_edge_nums){
    *_hex_edge_nums = hex_edge_nums;
  }
  else {
    delete [] hex_edge_nums;
  }
  if (_hex_face_nums){
    *_hex_face_nums = hex_face_nums;
  }
  else {
    delete [] hex_face_nums;
  }
}

/*
  Mesh the given geometry and retrieve either a regular mesh
*/
TMRMesh::TMRMesh( MPI_Comm _comm, TMRModel *_geo ){
  // Initialize the TMR-specific MPI data types
  if (!TMRIsInitialized()){
    TMRInitialize();
  }

  // Copy the communicator
  comm = _comm;
  geo = _geo;
  geo->incref();

  // Set the mesh properties
  num_nodes = 0;
  num_quads = 0;
  num_tris = 0;
  num_hex = 0;
  num_tet = 0;
  quads = NULL;
  tris = NULL;
  hex = NULL;
  tet = NULL;
  X = NULL;
}

TMRMesh::~TMRMesh(){
  // Free the meshes created here...
  int num_edges;
  TMREdge **edges;
  geo->getEdges(&num_edges, &edges);
  for ( int i = 0; i < num_edges; i++ ){
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);
    if (mesh){ delete mesh; }
    edges[i]->setMesh(NULL);
  }

  int num_faces;
  TMRFace **faces;
  geo->getFaces(&num_faces, &faces);
  for ( int i = 0; i < num_faces; i++ ){
    TMRFaceMesh *mesh = NULL;
    faces[i]->getMesh(&mesh);
    if (mesh){ delete mesh; }
    faces[i]->setMesh(NULL);
  }

  int num_volumes;
  TMRVolume **volumes;
  geo->getVolumes(&num_volumes, &volumes);
  for ( int i = 0; i < num_volumes; i++ ){
    TMRVolumeMesh *mesh = NULL;
    volumes[i]->getMesh(&mesh);
    if (mesh){ delete mesh; }
    volumes[i]->setMesh(NULL);
  }

  geo->decref();
  if (quads){ delete [] quads; }
  if (tris){ delete [] tris; }
  if (hex){ delete [] hex; }
  if (tet){ delete [] tet; }
  if (X){ delete [] X; }
}

/*
  Call the underlying mesher with the default options
*/
void TMRMesh::mesh( TMRMeshOptions options, double htarget ){
  TMRElementFeatureSize *fs = new TMRElementFeatureSize(htarget);
  fs->incref();
  mesh(options, fs);
  fs->decref();
}

/*
  Clear the old mesh - if any exists
*/
void TMRMesh::resetMesh(){
  // Reset the node numbers for the vertices
  int num_vertices = 0;
  TMRVertex **vertices;
  geo->getVertices(&num_vertices, &vertices);
  for ( int i = 0; i < num_vertices; i++ ){
    vertices[i]->resetNodeNum();
  }

  // Reset the edge meshes
  int num_edges;
  TMREdge **edges;
  geo->getEdges(&num_edges, &edges);
  for ( int i = 0; i < num_edges; i++ ){
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);
    if (mesh){
      mesh->decref();
      edges[i]->setMesh(NULL);
    }
  }

  // Reset the surface meshes
  int num_faces;
  TMRFace **faces;
  geo->getFaces(&num_faces, &faces);
  for ( int i = 0; i < num_faces; i++ ){
    TMRFaceMesh *mesh = NULL;
    faces[i]->getMesh(&mesh);
    if (mesh){
      mesh->decref();
      faces[i]->setMesh(NULL);
    }
  }

  // Reset the volume meshes (if any)
  int num_volumes;
  TMRVolume **volumes;
  geo->getVolumes(&num_volumes, &volumes);
  for ( int i = 0; i < num_volumes; i++ ){
    TMRVolumeMesh *mesh = NULL;
    volumes[i]->getMesh(&mesh);
    if (mesh){
      mesh->decref();
      volumes[i]->setMesh(NULL);
    }
  }
}

/*
  Mesh the underlying geometry
*/
void TMRMesh::mesh( TMRMeshOptions options,
                    TMRElementFeatureSize *fs ){
  // Reset the meshes within the mesh
  if (options.reset_mesh_objects){
    resetMesh();
  }

  // Mesh the curves
  int num_edges;
  TMREdge **edges;
  geo->getEdges(&num_edges, &edges);
  for ( int i = 0; i < num_edges; i++ ){
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);
    if (!mesh){
      mesh = new TMREdgeMesh(comm, edges[i]);
      mesh->mesh(options, fs);
      edges[i]->setMesh(mesh);
    }
  }

  // Mesh the surface
  int num_faces;
  TMRFace **faces;
  geo->getFaces(&num_faces, &faces);
  for ( int i = 0; i < num_faces; i++ ){
    TMRFaceMesh *mesh = NULL;
    faces[i]->getMesh(&mesh);
    if (!mesh){
      mesh = new TMRFaceMesh(comm, faces[i]);
      mesh->mesh(options, fs);
      faces[i]->setMesh(mesh);
    }
  }

  // Update target/source relationships
  int num_volumes;
  TMRVolume **volumes;
  geo->getVolumes(&num_volumes, &volumes);
  for ( int i = 0; i < num_volumes; i++ ){
    TMRVolumeMesh *mesh = NULL;
    volumes[i]->getMesh(&mesh);
    if (!mesh){
      mesh = new TMRVolumeMesh(comm, volumes[i]);
      int fail = mesh->mesh(options);
      if (fail){
        const char *name = volumes[i]->getName();
        if (name){
          fprintf(stderr,
                  "TMRMesh: Volume meshing failed for object %s\n",
                  name);
        }
        else {
          fprintf(stderr,
                  "TMRMesh: Volume meshing failed for volume %d\n", i);
        }
      }
      else {
        volumes[i]->setMesh(mesh);
      }
    }
  }

  // Now that we're done meshing, go ahead and uniquely order
  // the nodes in the mesh
  int num = 0;
  int num_vertices = 0;
  TMRVertex **vertices;
  geo->getVertices(&num_vertices, &vertices);
  for ( int i = 0; i < num_vertices; i++ ){
    vertices[i]->setNodeNum(&num);
  }

  // Order the edges
  for ( int i = 0; i < num_edges; i++ ){
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);
    mesh->setNodeNums(&num);
  }

  // Order the faces
  for ( int i = 0; i < num_faces; i++ ){
    TMRFaceMesh *mesh = NULL;
    faces[i]->getMesh(&mesh);
    mesh->setNodeNums(&num);
  }

  // Order the volumes
  for ( int i = 0; i < num_volumes; i++ ){
    TMRVolumeMesh *mesh = NULL;
    volumes[i]->getMesh(&mesh);
    mesh->setNodeNums(&num);
  }

  // Set the number of nodes in the mesh
  num_nodes = num;

  // Count up the number of quadrilaterals or hex in the mesh
  num_quads = 0;
  num_tris = 0;
  num_hex = 0;
  num_tet = 0;

  if (num_volumes > 0){
    for ( int i = 0; i < num_volumes; i++ ){
      TMRVolumeMesh *mesh = NULL;
      volumes[i]->getMesh(&mesh);
      num_hex += mesh->getHexConnectivity(NULL);
      num_tet += mesh->getTetConnectivity(NULL);
    }
  }
  if (num_faces > 0){
    for ( int i = 0; i < num_faces; i++ ){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);
      TMRFace *copy_face = NULL;
      faces[i]->getCopySource(NULL, &copy_face);
      if (!copy_face){
        num_quads += mesh->getQuadConnectivity(NULL);
        num_tris += mesh->getTriConnectivity(NULL);
      }
    }

    // Add the mesh quality to the bins from each mesh
    const int nbins = 20;
    int bins[nbins];
    memset(bins, 0, nbins*sizeof(int));
    for ( int i = 0; i < num_faces; i++ ){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);
      TMRFace *copy_face = NULL;
      faces[i]->getCopySource(NULL, &copy_face);
      if (!copy_face){
        mesh->addMeshQuality(nbins, bins);
      }
    }

    // Sum up the total number of elements
    int total = 0;
    for ( int i = 0; i < nbins; i++ ){
      total += bins[i];
    }

    // Get the MPI rank and print out the quality on the root processor
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    if (options.write_mesh_quality_histogram && mpi_rank == 0){
      printf("Quality   # elements   percentage\n");
      for ( int k = 0; k < nbins; k++ ){
        printf("< %.2f    %10d   %10.3f\n",
               1.0*(k+1)/nbins, bins[k], 100.0*bins[k]/total);
      }
      printf("          %10d\n", total);
    }
  }
}

/*
  Allocate and initialize the global mesh using the global ordering
*/
void TMRMesh::initMesh( int count_nodes ){
  // Allocate the global arrays
  X = new TMRPoint[ num_nodes ];

  if (num_hex > 0){
    hex = new int[ 8*num_hex ];

    int *count = NULL;
    if (count_nodes){
      count = new int[ num_nodes ];
      memset(count, 0, num_nodes*sizeof(int));
    }

    // Retrieve the surface information
    int num_volumes;
    TMRVolume **volumes;
    geo->getVolumes(&num_volumes, &volumes);

    // Set the values into the global arrays
    int *h = hex;
    for ( int i = 0; i < num_volumes; i++ ){
      // Get the mesh
      TMRVolumeMesh *mesh = NULL;
      volumes[i]->getMesh(&mesh);

      // Get the local mesh points
      int npts;
      TMRPoint *Xpts;
      mesh->getMeshPoints(&npts, &Xpts);

      // Get the local quadrilateral connectivity
      const int *hex_local;
      int nlocal = mesh->getHexConnectivity(&hex_local);

      // Get the local to global variable numbering
      const int *vars;
      mesh->getNodeNums(&vars);

      // Set the quadrilateral connectivity
      for ( int j = 0; j < 8*nlocal; j++, h++ ){
        h[0] = vars[hex_local[j]];
      }

      // Set the node locations
      for ( int j = 0; j < npts; j++ ){
        X[vars[j]] = Xpts[j];
      }

      if (count_nodes){
        for ( int j = 0; j < npts; j++ ){
          count[vars[j]]++;
        }
      }
    }

    if (count_nodes){
      // Set the count to the number of variables
      for ( int j = 0; j < num_nodes; j++ ){
        if (count[j] == 0){
          printf("TMRMesh error: Node %d not referenced\n", j);
        }
        if (count[j] >= 2){
          printf("TMRMesh error: Node %d referenced more than once %d\n",
                 j, count[j]);
        }
      }
      delete [] count;
    }
  }
  if (num_quads > 0 || num_tris > 0){
    if (num_tris > 0){
      tris = new int[ 3*num_tris ];
    }
    if (num_quads > 0){
      quads = new int[ 4*num_quads ];
    }

    // Retrieve the surface information
    int num_faces;
    TMRFace **faces;
    geo->getFaces(&num_faces, &faces);

    // Set the values into the global arrays
    int *q = quads;
    int *t = tris;
    for ( int i = 0; i < num_faces; i++ ){
      TMRFace *copy_face = NULL;
      faces[i]->getCopySource(NULL, &copy_face);
      if (!copy_face){
        // Get the mesh
        TMRFaceMesh *mesh = NULL;
        faces[i]->getMesh(&mesh);

        // Get the local mesh points
        int npts;
        TMRPoint *Xpts;
        mesh->getMeshPoints(&npts, NULL, &Xpts);

        // Get the local quadrilateral connectivity
        const int *quad_local, *tri_local;
        int nquad_local = mesh->getQuadConnectivity(&quad_local);
        int ntri_local = mesh->getTriConnectivity(&tri_local);

        // Get the local to global variable numbering
        const int *vars;
        mesh->getNodeNums(&vars);

        // Set the quadrilateral connectivity
        if (faces[i]->getOrientation() > 0){
          for ( int j = 0; j < 4*nquad_local; j++, q++ ){
            q[0] = vars[quad_local[j]];
          }
          for ( int j = 0; j < 3*ntri_local; j++, t++ ){
            t[0] = vars[tri_local[j]];
          }
        }
        else {
          for ( int j = 0; j < nquad_local; j++ ){
            for ( int k = 0; k < 4; k++ ){
              q[k] = vars[quad_local[4*j + k]];
            }

            // Flip the orientation of the quad
            int tmp = q[1];
            q[1] = q[3];
            q[3] = tmp;
            q += 4;
          }
          for ( int j = 0; j < ntri_local; j++ ){
            for ( int k = 0; k < 3; k++ ){
              t[k] = vars[tri_local[3*j + k]];
            }

            // Flip the orientation of the quad
            int tmp = t[1];
            t[1] = t[2];
            t[2] = tmp;
            t += 3;
          }
        }

        // Set the node locations
        for ( int j = 0; j < npts; j++ ){
          X[vars[j]] = Xpts[j];
        }
      }
    }
  }
}

/*
  Retrieve the mesh points (allocate them if they do not exist)
*/
int TMRMesh::getMeshPoints( TMRPoint **_X ){
  if (_X){
    if (!X){ initMesh(); }
    *_X = X;
  }
  return num_nodes;
}

/*
  Retrieve the underlying mesh connectivity (allocate it if it does
  not already exist)
*/
void TMRMesh::getQuadConnectivity( int *_nquads, const int **_quads ){
  if (!X){ initMesh(); }
  if (_nquads){ *_nquads = num_quads; }
  if (_quads){ *_quads = quads; }
}

/*
  Get the triangluar mesh connectivity
*/
void TMRMesh::getTriConnectivity( int *_ntris, const int **_tris ){
  if (!X){ initMesh(); }
  if (_ntris){ *_ntris = num_tris; }
  if (_tris){ *_tris = tris; }
}

/*
  Get the hexahedral connectivity
*/
void TMRMesh::getHexConnectivity( int *_nhex, const int **_hex ){
  if (!X){ initMesh(); }
  if (_nhex){ *_nhex = num_hex; }
  if (_hex){ *_hex = hex; }
}

/*
  Print out the mesh to a VTK file
*/
void TMRMesh::writeToVTK( const char *filename, int flag ){
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0 && num_nodes > 0){
    // Check whether to print out just quads or hex or both
    int nquad = num_quads;
    int nhex = num_hex;
    int ntris = num_tris;
    if (!(flag & TMR_QUAD)){
      nquad = 0;
    }
    if (!(flag & TMR_HEX)){
      nhex = 0;
    }

    if (!X){ initMesh(); }
    FILE *fp = fopen(filename, "w");
    if (fp){
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

      // Write out the points
      fprintf(fp, "POINTS %d float\n", num_nodes);
      for ( int k = 0; k < num_nodes; k++ ){
        fprintf(fp, "%e %e %e\n", X[k].x, X[k].y, X[k].z);
      }

      fprintf(fp, "\nCELLS %d %d\n", nquad + ntris + nhex,
              5*nquad + 4*ntris + 9*nhex);

      // Write out the cell connectivities
      for ( int k = 0; k < nquad; k++ ){
        fprintf(fp, "4 %d %d %d %d\n",
                quads[4*k], quads[4*k+1],
                quads[4*k+2], quads[4*k+3]);
      }
      for ( int k = 0; k < ntris; k++ ){
        fprintf(fp, "3 %d %d %d\n",
                tris[3*k], tris[3*k+1], tris[3*k+2]);
      }
      for ( int k = 0; k < nhex; k++ ){
        fprintf(fp, "8 %d %d %d %d %d %d %d %d\n",
                hex[8*k], hex[8*k+1], hex[8*k+2], hex[8*k+3],
                hex[8*k+4], hex[8*k+5], hex[8*k+6], hex[8*k+7]);
      }

      // All quadrilaterals
      fprintf(fp, "\nCELL_TYPES %d\n", nquad + ntris + nhex);
      for ( int k = 0; k < nquad; k++ ){
        fprintf(fp, "%d\n", 9);
      }
      for ( int k = 0; k < ntris; k++ ){
        fprintf(fp, "%d\n", 5);
      }
      for ( int k = 0; k < nhex; k++ ){
        fprintf(fp, "%d\n", 12);
      }

      fclose(fp);
    }
  }
}

/*
  Write the bulk data file with material properties
*/
void TMRMesh::writeToBDF( const char *filename, int flag ){
  // Static string for the beginning of the file
  const char nastran_file_header[] =
    "$ Generated by TMR\n$ NASTRAN input deck\nSOL 101\nBEGIN BULK\n";

  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0 && num_nodes > 0){
    if (!X){ initMesh(); }

    FILE *fp = fopen(filename, "w");
    if (fp){
      fprintf(fp, nastran_file_header);
      fprintf(fp, "$ Grid data\n");

      // Write out the coordinates to the BDF file
      int coord_disp = 0, coord_id = 0, seid = 0;
      for ( int i = 0; i < num_nodes; i++ ){
        fprintf(fp, "%-8s%16d%16d%16.9f%16.9f*%7d\n",
                "GRID*", i+1, coord_id,
                X[i].x, X[i].y, i+1);
        fprintf(fp, "*%7d%16.9f%16d%16s%16d        \n",
  	      i+1, X[i].z, coord_disp, " ", seid);
      }

      if (num_quads > 0 && (flag & TMR_QUAD)){
        int num_faces;
        TMRFace **faces;
        geo->getFaces(&num_faces, &faces);

        // Write out the element data for each component
        for ( int i = 0, j = 0; i < num_faces; i++ ){
          TMRFaceMesh *mesh = NULL;
          faces[i]->getMesh(&mesh);

          // Loop over all possible edges in the surface mesh
          const int *quad_local;
          int nlocal = mesh->getQuadConnectivity(&quad_local);

          // Get the local to global variable numbering
          const int *vars;
          mesh->getNodeNums(&vars);

          // Set the face id == name
          char descript[128];
          snprintf(descript, sizeof(descript), "FACE%d",
                   faces[i]->getEntityId());
          if (faces[i]->getName() != NULL){
            strncpy(descript, faces[i]->getName(), sizeof(descript));
          }

          // Print a local description of the face - use the entity
          // data if it exists, otherwise use the id value
          fprintf(fp, "%-41s","$       Shell element data");
          fprintf(fp, "%s\n", descript);

          if (faces[i]->getOrientation() > 0){
            for ( int k = 0; k < nlocal; k++, j++ ){
              int part = i+1;
              fprintf(fp, "%-8s%8d%8d%8d%8d%8d%8d%8d\n", "CQUADR",
                      j+1, part, vars[quad_local[4*k]]+1,
                      vars[quad_local[4*k+1]]+1,
                      vars[quad_local[4*k+2]]+1,
                      vars[quad_local[4*k+3]]+1, part);
            }
          }
          else {
            // Print out the nodes in the reversed orientation
            for ( int k = 0; k < nlocal; k++, j++ ){
              int part = i+1;
              fprintf(fp, "%-8s%8d%8d%8d%8d%8d%8d%8d\n", "CQUADR",
                      j+1, part, vars[quad_local[4*k]]+1,
                      vars[quad_local[4*k+3]]+1,
                      vars[quad_local[4*k+2]]+1,
                      vars[quad_local[4*k+1]]+1, part);
            }
          }
        }
      }
      if (num_hex > 0 && (flag & TMR_HEX)){
        int num_volumes;
        TMRVolume **volumes;
        geo->getVolumes(&num_volumes, &volumes);

        // Write out the element data for each component
        for ( int i = 0, j = 0; i < num_volumes; i++ ){
          TMRVolumeMesh *mesh = NULL;
          volumes[i]->getMesh(&mesh);

          // Loop over all possible edges in the surface mesh
          const int *hex_local;
          int nlocal = mesh->getHexConnectivity(&hex_local);

          // Get the local to global variable numbering
          const int *vars;
          mesh->getNodeNums(&vars);

          // Set the face id == name
          char descript[128];
          snprintf(descript, sizeof(descript), "VOLUME%d",
                   volumes[i]->getEntityId());
          if (volumes[i]->getName() != NULL){
            strncpy(descript, volumes[i]->getName(), sizeof(descript));
          }

          // Print a local description of the face - use the entity
          // data if it exists, otherwise use the id value
          int part = i+1;
          fprintf(fp, "%-41s","$       Volume element data");
          fprintf(fp, "%s\n", descript);
          // fprintf(fp, "%-8s%8d%8d%8d\n", "PSOLID", part, part, 0);

          for ( int k = 0; k < nlocal; k++, j++ ){
            fprintf(fp, "%-8s%8d%8d%8d%8d%8d%8d%8d%8d\n", "CHEXA",
                    j+1, part,
                    vars[hex_local[8*k]]+1,
                    vars[hex_local[8*k+1]]+1,
                    vars[hex_local[8*k+2]]+1,
                    vars[hex_local[8*k+3]]+1,
                    vars[hex_local[8*k+4]]+1,
                    vars[hex_local[8*k+5]]+1);
            fprintf(fp, "%-8s%8d%8d\n", " ",
                    vars[hex_local[8*k+6]]+1,
                    vars[hex_local[8*k+7]]+1);
          }
        }
      }

      // Signal end of bulk data section and close file handle
      fprintf(fp, "ENDDATA\n");
      fclose(fp);
    }
  }
}

/*
  Create the topology object, generating the vertices, edges and faces
  for each element in the underlying mesh.
*/
TMRModel* TMRMesh::createModelFromMesh(){
  // Initialize the mesh
  initMesh();

  // Create vertices
  TMRVertex **new_verts = new TMRVertex*[ num_nodes ];
  memset(new_verts, 0, num_nodes*sizeof(TMRVertex*));

  // Copy over the vertices in the original geometry
  int num_vertices;
  TMRVertex **vertices;
  geo->getVertices(&num_vertices, &vertices);
  for ( int i = 0; i < num_vertices; i++ ){
    TMRVertex *copy_vert;
    vertices[i]->getCopySource(&copy_vert);
    if (!copy_vert){
      int vnum;
      vertices[i]->getNodeNum(&vnum);
      new_verts[vnum] = vertices[i];
    }
  }

  // Create the vertices on the edges
  int num_edges;
  TMREdge **edges;
  geo->getEdges(&num_edges, &edges);
  for ( int i = 0; i < num_edges; i++ ){
    TMREdge *copy_edge;
    edges[i]->getCopySource(&copy_edge);
    if (!copy_edge){
      TMREdgeMesh *mesh = NULL;
      edges[i]->getMesh(&mesh);

      // Get the variable numbers
      const int *vars;
      mesh->getNodeNums(&vars);

      // Get the parametric points associated with the mesh
      int npts;
      const double *tpts;
      mesh->getMeshPoints(&npts, &tpts, NULL);
      for ( int j = 1; j < npts-1; j++ ){
        new_verts[vars[j]] = new TMRVertexFromEdge(edges[i], tpts[j]);
      }
    }
  }

  // Create the vertices from the faces
  int num_faces;
  TMRFace **faces;
  geo->getFaces(&num_faces, &faces);
  for ( int i = 0; i < num_faces; i++ ){
    TMRFace *copy_face;
    faces[i]->getCopySource(NULL, &copy_face);
    if (!copy_face){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);

      // Get the variable numbers
      const int *vars;
      mesh->getNodeNums(&vars);

      // Get the mesh points
      int npts;
      const double *pts;
      mesh->getMeshPoints(&npts, &pts, NULL);

      // Get the point-offset for this surface
      int offset = mesh->getNumFixedPoints();
      for ( int j = offset; j < npts; j++ ){
        new_verts[vars[j]] = new TMRVertexFromFace(faces[i],
                                                   pts[2*j], pts[2*j+1]);
      }
    }
  }

  // Create all the remaining nodes. These are associated with the
  // TMRVolumeMesh since they do not lie on an edge or face
  if (num_hex > 0){
    for ( int i = 0; i < num_nodes; i++ ){
      if (!new_verts[i]){
        new_verts[i] = new TMRVertexFromPoint(X[i]);
      }
    }
  }

  // The edges within the quadrilateral mesh
  int num_quad_edges = 0;
  int *quad_edges = NULL;

  // The edges within the hexahedral mesh
  int num_hex_edges = 0, num_hex_faces = 0;
  int *hex_edges = NULL, *hex_faces = NULL;
  int *hex_edge_nums = NULL, *hex_face_nums = NULL;

  if (num_hex > 0){
    // Compute the hexahedral edges and surfaces
    TMR_ComputeHexEdgesAndFaces(num_nodes, num_hex, hex,
                                &num_hex_edges, &hex_edges, &hex_edge_nums,
                                &num_hex_faces, &hex_faces, &hex_face_nums);
  }
  else {
    // Compute the quadrilateral edges
    TMR_ComputeQuadEdges(num_nodes, num_quads, quads,
                         &num_quad_edges, &quad_edges);
  }

  // Set a pointer to all of the edges
  int num_mesh_faces = num_quads;
  int num_mesh_edges = num_quad_edges;

  // Set the pointer to all of the edges within the mesh
  // This points to either the quad edges (if they were
  // computed) or the hex edges
  const int *mesh_edges = quad_edges;
  const int *mesh_faces = quads;
  if (hex_edges){
    num_mesh_edges = num_hex_edges;
    num_mesh_faces = num_hex_faces;

    mesh_edges = hex_edges;
    mesh_faces = hex_faces;
  }

  // Create a searchable array that stores the edge number and the two
  // connecting node numbers. The node numbers are sorted such that
  // the lowest one comes first. This enables fast sorting/searching
  int *sorted_edges = new int[ 3*num_mesh_edges ];

  // Populate the sorted_edges array, sorting each edge as we go...
  for ( int i = 0; i < num_mesh_edges; i++ ){
    sorted_edges[3*i] = i;

    // Set the quad edge, setting the lower edge number first
    if (mesh_edges[2*i] < mesh_edges[2*i+1]){
      sorted_edges[3*i+1] = mesh_edges[2*i];
      sorted_edges[3*i+2] = mesh_edges[2*i+1];
    }
    else {
      sorted_edges[3*i+1] = mesh_edges[2*i+1];
      sorted_edges[3*i+2] = mesh_edges[2*i];
    }
  }

  // Sort all of the edges/nodes so that they can be easily searched
  qsort(sorted_edges, num_mesh_edges, 3*sizeof(int), compare_edges);

  // Keep track of the orientation of the mesh edges.  edge_dir[i]
  // > 0 means that the ordering of the vertices in the sorted_edges
  // array is consistent with the orientation in the mesh. edge_dir[i]
  // < 0 means the edge is flipped relative to the sorted_edges array.
  int *edge_dir = new int[ num_mesh_edges ];
  memset(edge_dir, 0, num_mesh_edges*sizeof(int));

  // Create the edges for the faces/hexahedral elements for the
  // mesh from the underlying edges
  TMREdge **new_edges = new TMREdge*[ num_mesh_edges ];
  memset(new_edges, 0, num_mesh_edges*sizeof(TMREdge*));

  // Loop over the geometric edges within the model
  for ( int i = 0; i < num_edges; i++ ){
    // See if there is an underlying copy edge
    TMREdge *copy_edge = NULL;
    edges[i]->getCopySource(&copy_edge);

    TMREdgeMesh *mesh = NULL;
    if (copy_edge){
      copy_edge->getMesh(&mesh);
    }
    else {
      edges[i]->getMesh(&mesh);
    }

    // Get the global variables associated with the edge
    const int *vars;
    mesh->getNodeNums(&vars);

    // Get the parametric points associated with the mesh
    int npts;
    const double *tpts;
    mesh->getMeshPoints(&npts, &tpts, NULL);
    for ( int j = 0; j < npts-1; j++ ){
      // Find the edge associated with this curve
      int edge[3];
      edge[0] = 0;
      if (vars[j] < vars[j+1]){
        edge[1] = vars[j];
        edge[2] = vars[j+1];
      }
      else {
        edge[1] = vars[j+1];
        edge[2] = vars[j];
      }

      // Find the associated edge number
      int *res = (int*)bsearch(edge, sorted_edges, num_mesh_edges,
                               3*sizeof(int), compare_edges);

      if (res){
        int edge_num = res[0];
        if (!new_edges[edge_num]){
          // Check whether the node ordering is consistent with the
          // edge orientation. If not, tag this edge as reversed.
          if (vars[j] < vars[j+1]){
            edge_dir[edge_num] = 1;
          }
          else {
            edge_dir[edge_num] = -1;
          }

          new_edges[edge_num] =
            new TMRSplitEdge(edges[i], tpts[j], tpts[j+1]);
        }
      }
      else {
        fprintf(stderr,
                "TMRMesh error: Could not find edge (%d, %d) to split\n",
                edge[1], edge[2]);
      }
    }
  }

  // Create the edges on the surface using Pcurve/CurveFromSurface
  for ( int i = 0; i < num_faces; i++ ){
    TMRFace *copy_face = NULL;
    faces[i]->getCopySource(NULL, &copy_face);

    if (!copy_face){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);

      // Get the parametric points associated with the surface mesh
      int npts;
      const double *pts;
      mesh->getMeshPoints(&npts, &pts, NULL);

      // Loop over all possible edges in the surface mesh
      const int *quad_local;
      int nlocal = mesh->getQuadConnectivity(&quad_local);

      // Get the global variables associated with the local mesh
      const int *vars;
      mesh->getNodeNums(&vars);

      for ( int j = 0; j < nlocal; j++ ){
        for ( int k = 0; k < 4; k++ ){
          // The local variable numbers
          int l1 = 0, l2 = 0;
          if (faces[i]->getOrientation() > 0){
            l1 = quad_local[4*j + quad_edge_nodes[k][0]];
            l2 = quad_local[4*j + quad_edge_nodes[k][1]];
          }
          else {
            l1 = quad_local[4*j + flipped_quad_edge_nodes[k][0]];
            l2 = quad_local[4*j + flipped_quad_edge_nodes[k][1]];
          }

          // Get the global edge number
          int edge[3];
          edge[0] = 0;
          if (vars[l1] < vars[l2]){
            edge[1] = vars[l1];
            edge[2] = vars[l2];
          }
          else {
            edge[1] = vars[l2];
            edge[2] = vars[l1];
          }

          // Find the associated edge number
          int *res = (int*)bsearch(edge, sorted_edges, num_mesh_edges,
                                   3*sizeof(int), compare_edges);

          if (res){
            // Get the edge number
            int edge_num = res[0];
            if (!new_edges[edge_num]){
              // These edges are constructed such that they are
              // always in the forward orientation.
              edge_dir[edge_num] = 1;

              // Create the TMRBsplinePcurve on this edge
              double cpts[4];
              if (vars[l1] < vars[l2]){
                cpts[0] = pts[2*l1];
                cpts[1] = pts[2*l1+1];
                cpts[2] = pts[2*l2];
                cpts[3] = pts[2*l2+1];
              }
              else {
                cpts[0] = pts[2*l2];
                cpts[1] = pts[2*l2+1];
                cpts[2] = pts[2*l1];
                cpts[3] = pts[2*l1+1];
              }
              TMRBsplinePcurve *pcurve = new TMRBsplinePcurve(2, 2, cpts);
              new_edges[edge_num] = new TMREdgeFromFace(faces[i], pcurve);
            }
          }
          else {
            fprintf(stderr,
                    "TMRMesh error: Could not find edge (%d, %d) for Pcurve\n",
                    edge[1], edge[2]);
          }
        }
      }
    }
  }

  // Create the hexahedral elements
  if (num_hex > 0){
    for ( int i = 0; i < num_mesh_edges; i++ ){
      if (!new_edges[i]){
        // Get the edges
        int v1 = mesh_edges[2*i];
        int v2 = mesh_edges[2*i+1];

        // Set the edge direction and create the edge. Note that the
        // edge is always created in the positive orientation so that
        // edge_dir[i] = 1.
        edge_dir[i] = 1;
        if (v1 < v2){
          new_edges[i] = new TMRTFIEdge(new_verts[v1], new_verts[v2]);
        }
        else {
          new_edges[i] = new TMRTFIEdge(new_verts[v2], new_verts[v1]);
        }
      }
    }
  }

  int edge_count = 0;
  for ( int i = 0; i < num_mesh_edges; i++ ){
    if (!new_edges[i]){
      printf("edge %d (%d, %d) not created\n", i, mesh_edges[2*i], mesh_edges[2*i+1]);
      edge_count++;
    }
  }
  printf("edge_count = %d\n", edge_count);


  FILE *fp = fopen("unused_edges.vtk", "w");
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "vtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d float\n", 2*edge_count);
  int p3 = -1;
  for ( int i = 0; i < num_mesh_edges; i++ ){
    if (!new_edges[i]){
      int p1 = mesh_edges[2*i];
      int p2 = mesh_edges[2*i+1];
      if (i == 2){
        p3 = p1;
      }

      fprintf(fp, "%e %e %e\n", X[p1].x, X[p1].y, X[p1].z);
      fprintf(fp, "%e %e %e\n", X[p2].x, X[p2].y, X[p2].z);
    }
  }
  fprintf(fp, "\nCELLS %d %d\n", edge_count, 3*edge_count);
  for ( int i = 0; i < edge_count; i++ ){
    fprintf(fp, "%d %d %d\n", 2, 2*i, 2*i+1);
  }
  fprintf(fp, "\nCELL_TYPES %d\n", edge_count);
  for ( int i = 0; i < edge_count; i++ ){
    fprintf(fp, "%d\n", 3);
  }
  fclose(fp);

  // Find the closest node to p3
  double dist = 1e20;
  int index = -1;
  for ( int i = 0; i < num_nodes; i++ ){
    if (i != p3){
      double h = ((X[i].x - X[p3].x)*(X[i].x - X[p3].x) +
                  (X[i].y - X[p3].y)*(X[i].y - X[p3].y) +
                  (X[i].z - X[p3].z)*(X[i].z - X[p3].z));
      if (h < dist){
        dist = h;
        index = i;
      }
    }
  }

  printf("min dist point = %d dist = %15.8e\n", index, sqrt(dist));



  // Create the TMRFace objects
  TMRFace **new_faces = new TMRFace*[ num_mesh_faces ];
  memset(new_faces, 0, num_mesh_faces*sizeof(TMRFace*));

  // Create the face array
  int *sorted_faces = NULL;

  if (num_hex > 0){
    // Allocate an array to store the searchable array that
    // stores the face nodes
    sorted_faces = new int[ 5*num_hex_faces ];

    for ( int i = 0; i < num_hex_faces; i++ ){
      sorted_faces[5*i] = i;
      for ( int k = 0; k < 4; k++ ){
        sorted_faces[5*i+1+k] = hex_faces[4*i+k];
      }
      sort_face_nodes(&sorted_faces[5*i+1]);
    }

    // Sort all of the faces so that they can be easily searched
    qsort(sorted_faces, num_mesh_faces,
          5*sizeof(int), compare_faces);
  }

  // Allocate the faces
  int face_num = 0;
  for ( int i = 0; i < num_faces; i++ ){
    TMRFace *copy_face = NULL;
    faces[i]->getCopySource(NULL, &copy_face);

    if (!copy_face){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);

      // Get the parametric points associated with the surface mesh
      int npts;
      const double *pts;
      mesh->getMeshPoints(&npts, &pts, NULL);

      // Loop over all possible edges in the surface mesh
      const int *quad_local;
      int nlocal = mesh->getQuadConnectivity(&quad_local);

      // Get the global variables associated with the local mesh
      const int *vars;
      mesh->getNodeNums(&vars);

      // Iterate over all of the edges, creating the appropriate faces
      for ( int j = 0; j < nlocal; j++ ){
        TMREdge *c[4];
        int dir[4];
        TMRVertex *v[4];

        // Set the nodes for this quad
        if (faces[i]->getOrientation() > 0){
          v[0] = new_verts[vars[quad_local[4*j]]];
          v[1] = new_verts[vars[quad_local[4*j+1]]];
          v[2] = new_verts[vars[quad_local[4*j+2]]];
          v[3] = new_verts[vars[quad_local[4*j+3]]];
        }
        else {
          v[0] = new_verts[vars[quad_local[4*j]]];
          v[1] = new_verts[vars[quad_local[4*j+3]]];
          v[2] = new_verts[vars[quad_local[4*j+2]]];
          v[3] = new_verts[vars[quad_local[4*j+1]]];
        }

        // Loop over and search for the edges associated with this quad
        for ( int k = 0; k < 4; k++ ){
          // The edge variable numbers
          int l1 = 0, l2 = 0;
          if (faces[i]->getOrientation() > 0){
            l1 = quad_local[4*j + quad_edge_nodes[k][0]];
            l2 = quad_local[4*j + quad_edge_nodes[k][1]];
          }
          else {
            l1 = quad_local[4*j + flipped_quad_edge_nodes[k][0]];
            l2 = quad_local[4*j + flipped_quad_edge_nodes[k][1]];
          }

          // Get the global edge number
          int edge[3];
          edge[0] = 0;
          if (vars[l1] < vars[l2]){
            edge[1] = vars[l1];
            edge[2] = vars[l2];
          }
          else {
            edge[1] = vars[l2];
            edge[2] = vars[l1];
          }

          // Find the associated edge number
          int *res = (int*)bsearch(edge, sorted_edges, num_mesh_edges,
                                   3*sizeof(int), compare_edges);

          if (res){
            // Get the global edge number
            int edge_num = res[0];
            c[k] = new_edges[edge_num];

            if (vars[l1] < vars[l2]){
              dir[k] = 1;
            }
            else {
              dir[k] = -1;
            }
            dir[k] *= edge_dir[edge_num];
          }
          else {
            fprintf(stderr,
                    "TMRMesh error: Could not find edge (%d, %d) for surface\n",
                    edge[1], edge[2]);
          }
        }

        // If this is a hexahedral mesh, then we need to be consistent
        // with how the faces are ordered. This code searches for the
        // face number within the sorted_faces array to obtain the
        // required face number
        if (sorted_faces){
          // Set the nodes associated with this face and sort them
          int face[5];
          face[0] = 0;
          for ( int k = 0; k < 4; k++ ){
            face[k+1] = vars[quad_local[4*j + k]];
          }
          sort_face_nodes(&face[1]);

          // Search for the face
          int *res = (int*)bsearch(face, sorted_faces, num_mesh_faces,
                                   5*sizeof(int), compare_faces);

          // Set the face number
          if (res){
            face_num = res[0];
          }
        }

        // Create the parametric TFI surface
        new_faces[face_num] = new TMRParametricTFIFace(faces[i], c, dir, v);
        face_num++;
      }
    }
  }

  // Create the remaining faces
  if (num_hex > 0){
    for ( int i = 0; i < num_mesh_faces; i++ ){
      if (!new_faces[i]){
        // The edge, direction and vertex information
        TMREdge *c[4];
        int dir[4];
        TMRVertex *v[4];

        // Loop over all of the edges/nodes associated with
        // this hexahedral face
        for ( int k = 0; k < 4; k++ ){
          v[k] = new_verts[mesh_faces[4*i + k]];

          // Get the edge numbers associated with this edge
          int e0 = mesh_faces[4*i + quad_edge_nodes[k][0]];
          int e1 = mesh_faces[4*i + quad_edge_nodes[k][1]];

          // Get the global edge number
          int edge[3];
          edge[0] = 0;
          dir[k] = 1;
          if (e0 < e1){
            edge[1] = e0;
            edge[2] = e1;
          }
          else {
            dir[k] = -1;
            edge[1] = e1;
            edge[2] = e0;
          }

          // Find the associated edge number
          int *res = (int*)bsearch(edge, sorted_edges, num_mesh_edges,
                                   3*sizeof(int), compare_edges);

          if (res){
            // Get the global edge number
            int edge_num = res[0];
            c[k] = new_edges[edge_num];
            dir[k] *= edge_dir[edge_num];
          }
          else {
            fprintf(stderr,
                    "TMRMesh error: Could not find edge (%d, %d) for hex\n",
                    edge[1], edge[2]);
          }
        }

        new_faces[i] = new TMRTFIFace(c, dir, v);
      }
    }
  }

  // Free all of the edge search information
  delete [] sorted_edges;
  if (sorted_faces){
    delete [] sorted_faces;
  }

  TMRVolume **new_volumes = NULL;

  if (num_hex > 0){
    // Create the new volume array
    new_volumes = new TMRVolume*[ num_hex ];

    for ( int i = 0; i < num_hex; i++ ){
      // Get the edges
      TMRVertex *v[8];
      TMREdge *e[12];
      TMRFace *f[6];
      int edir[12], fdir[6];

      // Get the vertices
      for ( int j = 0; j < 8; j++ ){
        v[j] = new_verts[hex[8*i + hex_coordinate_order[j]]];
      }

      // Get the edges and their directions
      for ( int j = 0; j < 12; j++ ){
        int edge_num = hex_edge_nums[12*i + j];
        edir[j] = 1;
        if (hex[8*i + hex_edge_nodes[j][0]] >
            hex[8*i + hex_edge_nodes[j][1]]){
          edir[j] = -1;
        }
        edir[j] *= edge_dir[edge_num];
        e[j] = new_edges[edge_num];
      }

      // Get the faces and their orientation
      for ( int j = 0; j < 6; j++ ){
        int face_num = hex_face_nums[6*i + j];
        f[j] = new_faces[face_num];
      }

      // Allocate the new transfinite interpolation face
      new_volumes[i] = new TMRTFIVolume(f, e, edir, v);
    }
  }

  // Free the edge directions
  delete [] edge_dir;

  // Free all the auxiliary connectivity data
  if (quad_edges){
    delete [] quad_edges;
  }
  if (hex_edges){
    delete [] hex_edges;
    delete [] hex_faces;
    delete [] hex_edge_nums;
    delete [] hex_face_nums;
  }

  // Create the geometry object
  TMRModel *geo = NULL;

  if (num_hex > 0){
    geo = new TMRModel(num_nodes, new_verts,
                       num_mesh_edges, new_edges,
                       num_mesh_faces, new_faces,
                       num_hex, new_volumes);
  }
  else {
    geo = new TMRModel(num_nodes, new_verts,
                       num_mesh_edges, new_edges,
                       num_mesh_faces, new_faces);
  }

  // Free the data that was locally allocated
  delete [] new_verts;
  delete [] new_edges;
  delete [] new_faces;
  delete [] new_volumes;

  return geo;
}
