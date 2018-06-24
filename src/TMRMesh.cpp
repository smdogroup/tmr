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

#include "TMRMesh.h"
#include "TMRNativeTopology.h"
#include "TMRTriangularize.h"
#include "TMRMeshSmoothing.h"
#include "TMRPerfectMatchInterface.h"
#include "TMRBspline.h"
#include "tmrlapack.h"
#include <math.h>
#include <stdio.h>
#include <map>

#ifdef TMR_USE_NETGEN
// The namespace is required because of the way nglib is compiled by
// default. This is a funny way to do it, but who am I to complain.
namespace nglib {
#include "nglib.h"
}
#endif // TMR_USE_NETGEN

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
const int face_orient[8][4] =
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
static void computeNodeToElems( int nnodes, int nelems,
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
static void computeTriEdges( int nnodes, int ntris,
                             const int tris[],
                             int *num_tri_edges,
                             int **_tri_edges,
                             int **_tri_neighbors,
                             int **_dual_edges,
                             int **_node_to_tri_ptr=NULL,
                             int **_node_to_tris=NULL,
                             int **_tri_edge_nums=NULL ){
  // Compute the edges in the triangular mesh
  int *ptr;
  int *node_to_tris;
  computeNodeToElems(nnodes, ntris, 3, tris, &ptr, &node_to_tris);

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
static void computeQuadEdges( int nnodes, int nquads,
                              const int quads[],
                              int *_num_quad_edges,
                              int **_quad_edges,
                              int **_quad_neighbors=NULL,
                              int **_dual_edges=NULL,
                              int **_quad_edge_nums=NULL ){
  // Compute the connectivity from nodes to quads
  int *ptr;
  int *node_to_quads;
  computeNodeToElems(nnodes, nquads, 4, quads, &ptr, &node_to_quads);

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
static void computeHexEdgesAndFaces( int nnodes, int nhex,
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
  computeNodeToElems(nnodes, nhex, 8, hex, &ptr, &node_to_hex);

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
              if (f[0] == f2[face_orient[ort][0]] &&
                  f[1] == f2[face_orient[ort][1]] &&
                  f[2] == f2[face_orient[ort][2]] &&
                  f[3] == f2[face_orient[ort][3]]){
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
  An integral entry for the linked list
*/
class IntegralPt {
 public:
  double t;
  double dist;
  IntegralPt *next;
};

/*
  Evaluate the distance between two points
*/
double pointDist( TMRPoint *a, TMRPoint *b ){
  return sqrt((a->x - b->x)*(a->x - b->x) +
              (a->y - b->y)*(a->y - b->y) +
              (a->z - b->z)*(a->z - b->z));
}

/*
  Recursive integration on an edge with an adaptive error control to
  ensure that the integral is computed with sufficient accuracy.

  input:
  t1, t2:  the limits of integration for this interval
  tol:     the absolute error measure
  ncalls:  the recursion depth
  pt:      pointer into the linked list
*/
void integrateEdge( TMREdge *edge, TMRElementFeatureSize *fs,
                    double t1, double h1, TMRPoint p1,
                    double t2, double tol,
                    int ncalls, IntegralPt **_pt ){
  // Dereference the pointer to the integral point
  IntegralPt *pt = *_pt;

  // Find the mid point of the interval
  TMRPoint pmid;
  double tmid = 0.5*(t1 + t2);
  edge->evalPoint(tmid, &pmid);
  double hmid = fs->getFeatureSize(pmid);

  // Evaluate the point at the end of the interval
  TMRPoint p2;
  edge->evalPoint(t2, &p2);
  double h2 = fs->getFeatureSize(p2);

  // Evaluate the approximate integral contributions
  double int1 = 2.0*pointDist(&p1, &pmid)/(h1 + hmid);
  double int2 = 4.0*pointDist(&pmid, &p2)/(h1 + 2.0*hmid + h2);
  double int3 = 2.0*pointDist(&p1, &p2)/(hmid + h2);

  // Compute the integration error
  double error = fabs(int3 - int1 - int2);

  if (((ncalls > 6) && (error < tol)) || (ncalls > 20)){
    // Add the mid point
    pt->next = new IntegralPt;
    pt->next->dist = pt->dist + int1;
    pt->next->t = tmid;
    pt->next->next = NULL;
    pt = pt->next;

    // Add the final point p2
    pt->next = new IntegralPt;
    pt->next->dist = pt->dist + int2;
    pt->next->t = t2;
    pt->next->next = NULL;
    pt = pt->next;

    // Set the pointer to the end of the linked list
    *_pt = pt;
  }
  else {
    // Continue the recursive integration
    integrateEdge(edge, fs, t1, h1, p1, tmid, tol, ncalls+1, _pt);
    integrateEdge(edge, fs, tmid, hmid, pmid, t2, tol, ncalls+1, _pt);
  }
}

/*
  Integrate along the edge adaptively, creating a list
*/
double integrateEdge( TMREdge *edge, TMRElementFeatureSize *fs,
                      double t1, double t2, double tol,
                      double **_tvals, double **_dist,
                      int *_nvals ){
  *_tvals = NULL;
  *_dist = NULL;
  *_nvals = 0;

  // Allocate the entry in the linked list
  IntegralPt *root = new IntegralPt;
  root->next = NULL;
  root->dist = 0.0;
  root->t = t1;

  // Evaluate the first point
  TMRPoint p1;
  edge->evalPoint(t1, &p1);
  double h1 = fs->getFeatureSize(p1);

  // Integrate over the edge
  IntegralPt *pt = root;
  integrateEdge(edge, fs, t1, h1, p1, t2, tol, 0, &pt);

  // Count up and allocate the num
  int count = 1;
  IntegralPt *curr = root;
  while (curr->next){
    curr = curr->next;
    count++;
  }

  // Allocate arrays to store the parametric location/distance data
  double *tvals = new double[ count ];
  double *dist = new double[ count ];

  // Scan through the linked list, read out the values of the
  // parameter and its integral and delete the entries as we go...
  count = 0;
  curr = root;
  tvals[count] = curr->t;
  dist[count] = curr->dist;
  count++;

  while (curr->next){
    IntegralPt *tmp = curr;
    curr = curr->next;
    tvals[count] = curr->t;
    dist[count] = curr->dist;
    count++;
    delete tmp;
  }

  double len = curr->dist;
  delete curr;

  // Set the pointers for the output
  *_nvals = count;
  *_tvals = tvals;
  *_dist = dist;

  return len;
}

/*
  Create the element feature size.

  This corresponds to a constant feature size across the entire
  domain. Element grading can be implemented by overriding this class
  to modify the feature size as a function of position.
*/
TMRElementFeatureSize::TMRElementFeatureSize( double _hmin ){
  hmin = _hmin;
}

TMRElementFeatureSize::~TMRElementFeatureSize(){}

/*
  Return the feature size - hmin
*/
double TMRElementFeatureSize::getFeatureSize( TMRPoint pt ){
  return hmin;
}

/*
  Create a feature size dependency that is linear but does not
  exceed hmin or hmax anywhere in the domain
*/
TMRLinearElementSize::TMRLinearElementSize( double _hmin, double _hmax,
                                            double _c,
                                            double _ax,
                                            double _ay,
                                            double _az ):
TMRElementFeatureSize(_hmin){
  hmax = _hmax;
  c = _c;
  ax = _ax;
  ay = _ay;
  az = _az;
}

TMRLinearElementSize::~TMRLinearElementSize(){}

/*
  Get the feature size based on the input coefficients
*/
double TMRLinearElementSize::getFeatureSize( TMRPoint pt ){
  double h = c + ax*pt.x + ay*pt.y + az*pt.z;
  if (h < hmin){ h = hmin; }
  if (h > hmax){ h = hmax; }
  return h;
}

/*
  Create the feature size within a box
*/
TMRBoxFeatureSize::TMRBoxFeatureSize( TMRPoint p1, TMRPoint p2,
                                      double _hmin, double _hmax ):
TMRElementFeatureSize(_hmin){
  hmax = _hmax;

  TMRPoint m1, d1;
  m1.x = 0.5*(p1.x + p2.x);
  m1.y = 0.5*(p1.y + p2.y);
  m1.z = 0.5*(p1.z + p2.z);
  d1.x = 0.5*fabs(p1.x - p2.x);
  d1.y = 0.5*fabs(p1.y - p2.y);
  d1.z = 0.5*fabs(p1.z - p2.z);

  // Allocate the root for the boxes
  list_current = new BoxList();
  list_root = list_current;

  // Set the initial number of boxes
  num_boxes = 0;

  // Allocate a box
  list_current->boxes[num_boxes].m = m1;
  list_current->boxes[num_boxes].d = d1;
  list_current->boxes[num_boxes].h = hmax;

  // Add a box to the data structure
  root = new BoxNode(&(list_current->boxes[num_boxes]), m1, d1);
  num_boxes++;
}

/*
  Free the data allocated by the feature-size object
*/
TMRBoxFeatureSize::~TMRBoxFeatureSize(){
  delete root;
  while (list_root){
    BoxList *tmp = list_root;
    list_root = list_root->next;
    delete tmp;
  }
}

/*
  Add a box that is designed to constrain the feature size
*/
void TMRBoxFeatureSize::addBox( TMRPoint p1, TMRPoint p2, double hval ){
  TMRPoint m1, d1;
  m1.x = 0.5*(p1.x + p2.x);
  m1.y = 0.5*(p1.y + p2.y);
  m1.z = 0.5*(p1.z + p2.z);
  d1.x = 0.5*fabs(p1.x - p2.x);
  d1.y = 0.5*fabs(p1.y - p2.y);
  d1.z = 0.5*fabs(p1.z - p2.z);

  // Check if we need to expand the array of boxes
  if (num_boxes >= MAX_LIST_BOXES){
    list_current->next = new BoxList();
    list_current = list_current->next;
    num_boxes = 0;
  }

  // Allocate a box
  list_current->boxes[num_boxes].m = m1;
  list_current->boxes[num_boxes].d = d1;
  list_current->boxes[num_boxes].h = hval;

  // Add a box to the data structure
  root->addBox(&(list_current->boxes[num_boxes]));
  num_boxes++;
}

/*
  Retrieve the feature size
*/
double TMRBoxFeatureSize::getFeatureSize( TMRPoint pt ){
  double h = hmax;
  root->getSize(pt, &h);
  if (h < hmin){ h = hmin; }
  if (h > hmax){ h = hmax; }
  return h;
}

/*
  Check if the box contains the point
*/
int TMRBoxFeatureSize::BoxSize::contains( TMRPoint p ){
  // Compute the lower/upper bounds
  double xl = m.x - d.x;
  double xu = m.x + d.x;
  double yl = m.y - d.y;
  double yu = m.y + d.y;
  double zl = m.z - d.z;
  double zu = m.z + d.z;

  // Return true if the box contains the point
  if ((p.x >= xl && p.x <= xu) &&
      (p.y >= yl && p.y <= yu) &&
      (p.z >= zl && p.z <= zu)){
    return 1;
  }
  return 0;
}

/*
  Create the box node and allocate a single box (for now)
*/
TMRBoxFeatureSize::BoxNode::BoxNode( BoxSize *cover,
                                     TMRPoint m1,
                                     TMRPoint d1 ){
  m = m1;
  d = d1;
  num_boxes = 1;
  boxes = new BoxSize*[ MAX_NUM_BOXES ];
  boxes[0] = cover;
  memset(c, 0, sizeof(c));
}

/*
  Deallocate all of the data
*/
TMRBoxFeatureSize::BoxNode::~BoxNode(){
  if (boxes){
    delete [] boxes;
  }
  for ( int i = 0; i < 8; i++ ){
    if (c[i]){
      delete c[i];
    }
  }
}

/*
  Add a box to the octree data structure
*/
void TMRBoxFeatureSize::BoxNode::addBox( BoxSize *ptr ){
  // Get the dimensions for the box
  double xl = ptr->m.x - ptr->d.x;
  double xu = ptr->m.x + ptr->d.x;
  double yl = ptr->m.y - ptr->d.y;
  double yu = ptr->m.y + ptr->d.y;
  double zl = ptr->m.z - ptr->d.z;
  double zu = ptr->m.z + ptr->d.z;

  // Get the dimensions for this node
  double nxl = m.x - d.x;
  double nxu = m.x + d.x;
  double nyl = m.y - d.y;
  double nyu = m.y + d.y;
  double nzl = m.z - d.z;
  double nzu = m.z + d.z;

  // Check if this box coverss any part of this node, if not we're
  // done
  if (!((xu >= nxl && xl <= nxu) &&
        (yu >= nyl && yl <= nyu) &&
        (zu >= nzl && zl <= nzu))){
    return;
  }

  // Check if the boxes at this level are currently in use
  if (boxes){
    // Check if this box covers the entire node or not
    if ((xl <= nxl && xu >= nxu) &&
        (yl <= nyl && yu >= nyu) &&
        (zl <= nzl && zu >= nzu)){
      // This box fully covers the node. Check if the h-dimension is
      // more strict than the current covering box
      if (boxes[0]->h > ptr->h){
        boxes[0] = ptr;
      }
      return;
    }

    // Otherwise, the box only covers part of the node and we have
    // to add it to the list of local boxes
    if (num_boxes < MAX_NUM_BOXES){
      boxes[num_boxes] = ptr;
      num_boxes++;
      return;
    }
    else {
      // Allocate new children
      for ( int k = 0; k < 2; k++ ){
        double zl = m.z + (k-1)*d.z;
        double zu = m.z + k*d.z;

        for ( int j = 0; j < 2; j++ ){
          double yl = m.y + (j-1)*d.y;
          double yu = m.y + j*d.y;

          for ( int i = 0; i < 2; i++ ){
            double xl = m.x + (i-1)*d.x;
            double xu = m.x + i*d.x;

            // Set the mid-points and distances for each child
            TMRPoint m1, d1;
            m1.x = 0.5*(xu + xl);
            m1.y = 0.5*(yu + yl);
            m1.z = 0.5*(zu + zl);
            d1.x = 0.5*(xu - xl);
            d1.y = 0.5*(yu - yl);
            d1.z = 0.5*(zu - zl);

            // Create the child nodes
            c[i + 2*j + 4*k] = new BoxNode(boxes[0], m1, d1);
          }
        }
      }

      // Add the boxes to the rest of the tree
      while (num_boxes > 0){
        for ( int k = 0; k < 2; k++ ){
          for ( int j = 0; j < 2; j++ ){
            for ( int i = 0; i < 2; i++ ){
              c[i + 2*j + 4*k]->addBox(boxes[num_boxes-1]);
            }
          }
        }
        num_boxes--;
      }

      // Free the boxes
      delete [] boxes;
      boxes = NULL;
    }
  }

  // Add the boxes to the child - if applicable
  for ( int k = 0; k < 2; k++ ){
    for ( int j = 0; j < 2; j++ ){
      for ( int i = 0; i < 2; i++ ){
        c[i + 2*j + 4*k]->addBox(ptr);
      }
    }
  }
}

/*
  Get the element size based on the extent of the box
*/
void TMRBoxFeatureSize::BoxNode::getSize( TMRPoint p, double *hval ){
  // Scan through the tree and find the most-constraining box size
  for ( int k = 0; k < 2; k++ ){
    double zl = m.z + (k-1)*d.z;
    double zu = m.z + k*d.z;

    for ( int j = 0; j < 2; j++ ){
      double yl = m.y + (j-1)*d.y;
      double yu = m.y + j*d.y;

      for ( int i = 0; i < 2; i++ ){
        double xl = m.x + (i-1)*d.x;
        double xu = m.x + i*d.x;

        // Check if the point lies in the box and the box exists
        if ((p.x >= xl && p.x <= xu) &&
            (p.y >= yl && p.y <= yu) &&
            (p.z >= zl && p.z <= zu)){
          if (c[i + 2*j + 4*k]){
            return c[i + 2*j + 4*k]->getSize(p, hval);
          }
        }
      }
    }
  }

  for ( int i = 0; i < num_boxes; i++ ){
    if (boxes[i]->contains(p) && *hval > boxes[i]->h){
      *hval = boxes[i]->h;
    }
  }
}

/*
  Create a mesh along curve
*/
TMREdgeMesh::TMREdgeMesh( MPI_Comm _comm, TMREdge *_edge ){
  comm = _comm;
  edge = _edge;
  edge->incref();

  npts = 0;
  pts = NULL;
  X = NULL;
  vars = NULL;
}

/*
  Destroy the mesh for this curve, and free the underlying data
*/
TMREdgeMesh::~TMREdgeMesh(){
  edge->decref();
  if (pts){ delete [] pts; }
  if (X){ delete [] X; }
  if (vars){ delete [] vars; }
}

/*
  Find the points along the edge for the mesh
*/
void TMREdgeMesh::mesh( TMRMeshOptions options,
                        TMRElementFeatureSize *fs ){
  int mpi_rank, mpi_size;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Get the source edge
  TMREdge *source;
  edge->getSource(&source);

  // Figure out if there is a source edge and whether or not it has
  // been meshed.
  npts = -1;
  if (source && source != edge){
    TMREdgeMesh *mesh;
    source->getMesh(&mesh);
    if (!mesh){
      mesh = new TMREdgeMesh(comm, source);
      mesh->mesh(options, fs);
      source->setMesh(mesh);
    }

    // Retrieve the number of points along the source edge
    mesh->getMeshPoints(&npts, NULL, NULL);
  }

  if (mpi_rank == 0){
    // Get the limits of integration that will be used
    double tmin, tmax;
    edge->getRange(&tmin, &tmax);

    // Get the associated vertices
    TMRVertex *v1, *v2;
    edge->getVertices(&v1, &v2);

    if (!edge->isDegenerate()){
      // Set the integration error tolerance
      double integration_eps = 1e-8;

      // Integrate along the curve to obtain the distance function such
      // that dist(tvals[i]) = int_{tmin}^{tvals[i]} ||d{C(t)}dt||_{2} dt
      int nvals;
      double *dist, *tvals;
      integrateEdge(edge, fs, tmin, tmax, integration_eps,
                    &tvals, &dist, &nvals);

      // Only compute the number of points if there is no source edge
      if (npts < 0){
        // Compute the number of points along this curve
        npts = (int)(ceil(dist[nvals-1]));
        if (npts < 2){ npts = 2; }

        // If we have an even number of points, increment by one to ensure
        // that we have an even number of segments along the boundary
        if (npts % 2 != 1){ npts++; }

        // If the start/end vertex are the same, then the minimum number
        // of points is 5
        if ((v1 == v2) && npts < 5){
          npts = 5;
        }
      }

      // The average non-dimensional distance between points
      double d = dist[nvals-1]/(npts-1);

      // Allocate the parametric points that will be used
      pts = new double[ npts ];

      // Set the starting/end location of the points
      pts[0] = tmin;
      pts[npts-1] = tmax;

      // Perform the integration so that the points are evenly spaced
      // along the curve
      for ( int j = 1, k = 1; (j < nvals && k < npts-1); j++ ){
        while ((k < npts-1) &&
               (dist[j-1] <= d*k && d*k < dist[j])){
          double u = 0.0;
          if (dist[j] > dist[j-1]){
            u = (d*k - dist[j-1])/(dist[j] - dist[j-1]);
          }
          pts[k] = tvals[j-1] + (tvals[j] - tvals[j-1])*u;
          k++;
        }
      }

      // Free the integration result
      delete [] tvals;
      delete [] dist;
    }
    else {
      // This is a degenerate edge
      npts = 2;
      pts = new double[ npts ];
      pts[0] = tmin;
      pts[1] = tmax;
    }

    // Allocate the points
    X = new TMRPoint[ npts ];
    for ( int i = 0; i < npts; i++ ){
      edge->evalPoint(pts[i], &X[i]);
    }
  }

  if (mpi_size > 1){
    // Broadcast the number of points to all the processors
    MPI_Bcast(&npts, 1, MPI_INT, 0, comm);

    if (mpi_rank != 0){
      pts = new double[ npts ];
      X = new TMRPoint[ npts ];
    }

    // Broadcast the parametric locations and points
    MPI_Bcast(pts, npts, MPI_DOUBLE, 0, comm);
    MPI_Bcast(X, npts, TMRPoint_MPI_type, 0, comm);
  }
}

/*
  Order the internal mesh points and return the number of owned
  points that were ordered.
*/
int TMREdgeMesh::setNodeNums( int *num ){
  if (!vars && pts){
    // Retrieve the vertices
    TMRVertex *v1, *v2;
    edge->getVertices(&v1, &v2);

    // Allocate/set the node numbers
    vars = new int[ npts ];

    // Get the variable numbers
    v1->getNodeNum(&vars[0]);
    v2->getNodeNum(&vars[npts-1]);

    // Set the internal node numbers
    for ( int i = 1; i < npts-1; i++ ){
      vars[i] = *num;
      (*num)++;
    }

    return npts-2;
  }

  return 0;
}

/*
  Retrieve the internal node numbers, associated with the same
  set of points/same order as the getMeshPoints code
*/
int TMREdgeMesh::getNodeNums( const int **_vars ){
  if (_vars){
    *_vars = vars;
  }
  return npts;
}

/*
  Get the mesh points
*/
void TMREdgeMesh::getMeshPoints( int *_npts,
                                 const double **_pts,
                                 TMRPoint **_X ){
  if (_npts){ *_npts = npts; }
  if (_pts){ *_pts = pts; }
  if (_X){ *_X = X; }
}

/*
  Create the surface mesh object

  This call does not create the underlying surface. You can set
  meshing options into the class after it is created, and then create
  the mesh.

  Note that the curve/edge meshes must be meshed before calling this
  object.
*/
TMRFaceMesh::TMRFaceMesh( MPI_Comm _comm, TMRFace *_face ){
  comm = _comm;
  face = _face;
  face->incref();
  mesh_type = TMR_NO_MESH;

  num_fixed_pts = 0;
  num_points = 0;
  num_quads = 0;
  num_tris = 0;

  // NULL things that will be used later
  pts = NULL;
  X = NULL;
  vars = NULL;
  quads = NULL;
  tris = NULL;
}

/*
  Free the data associated with the surface mesh
*/
TMRFaceMesh::~TMRFaceMesh(){
  face->decref();
  if (pts){ delete [] pts; }
  if (X){ delete [] X; }
  if (vars){ delete [] vars; }
  if (quads){ delete [] quads; }
  if (tris){ delete [] tris; }
}

/*
  Create the surface mesh
*/
void TMRFaceMesh::mesh( TMRMeshOptions options,
                        TMRElementFeatureSize *fs ){
  int mpi_rank, mpi_size;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Set the default mesh type
  TMRFaceMeshType _mesh_type = options.mesh_type_default;
  if (_mesh_type == TMR_NO_MESH){
    _mesh_type = TMR_STRUCTURED;
  }

  // Get the source face and its orientation relative to this
  // face. Note that the source face may be NULL in which case the
  // source orientation is meaningless.
  int source_dir;
  TMRVolume *source_volume;
  TMRFace *source;
  face->getSource(&source_dir, &source_volume, &source);

  if (source){
    // If the face mesh for the source does not yet exist, create it...
    TMRFaceMesh *face_mesh;
    source->getMesh(&face_mesh);
    if (!face_mesh){
      face_mesh = new TMRFaceMesh(comm, source);
      face_mesh->mesh(options, fs);
      source->setMesh(face_mesh);
    }
  }

  // First check if the conditions for a structured mesh are satisfied
  if (_mesh_type == TMR_STRUCTURED){
    int nloops = face->getNumEdgeLoops();
    if (nloops != 1){
      _mesh_type = TMR_UNSTRUCTURED;
    }

    // Get the first edge loop and the edges in the loop
    TMREdgeLoop *loop;
    face->getEdgeLoop(0, &loop);
    int nedges;
    TMREdge **edges;
    loop->getEdgeLoop(&nedges, &edges, NULL);
    if (nedges != 4){
      _mesh_type = TMR_UNSTRUCTURED;
    }
    for ( int k = 0; k < nedges; k++ ){
      if (edges[k]->isDegenerate()){
        _mesh_type = TMR_UNSTRUCTURED;
      }
    }

    if (nedges == 4){
      // Check that parallel edges have the same number of nodes
      int nx1 = 0, nx2 = 0, ny1 = 0, ny2 = 0;
      TMREdgeMesh *mesh;
      for ( int k = 0; k < 4; k++ ){
        edges[k]->getMesh(&mesh);
        if (!mesh){
          fprintf(stderr,
                  "TMRFaceMesh error: Edge mesh does not exist\n");
        }
      }

      edges[0]->getMesh(&mesh);
      mesh->getMeshPoints(&nx1, NULL, NULL);
      edges[2]->getMesh(&mesh);
      mesh->getMeshPoints(&nx2, NULL, NULL);
      if (nx1 != nx2){
        _mesh_type = TMR_UNSTRUCTURED;
      }

      edges[1]->getMesh(&mesh);
      mesh->getMeshPoints(&ny1, NULL, NULL);
      edges[3]->getMesh(&mesh);
      mesh->getMeshPoints(&ny2, NULL, NULL);
      if (ny1 != ny2){
        _mesh_type = TMR_UNSTRUCTURED;
      }
    }
  }

  // Record the mesh type
  mesh_type = _mesh_type;

  if (mpi_rank == 0){
    // Count up the number of points and segments from the curves that
    // bound the surface. Keep track of the number of points = the
    // number of segments.
    int total_num_pts = 0;

    // Keep track of the number of closed loop cycles in the domain
    int nloops = face->getNumEdgeLoops();

    // The number of degenerate edges
    int num_degen = 0;

    // The total number of edges in the model
    int total_nedges = 0;

    // Get all of the edges and count up the mesh points
    for ( int k = 0; k < nloops; k++ ){
      TMREdgeLoop *loop;
      face->getEdgeLoop(k, &loop);
      int nedges;
      TMREdge **edges;
      loop->getEdgeLoop(&nedges, &edges, NULL);

      // Keep track of the total number of edges attached to this
      // surface object
      total_nedges += nedges;

      for ( int i = 0; i < nedges; i++ ){
        // Count whether this edge is degenerate
        if (edges[i]->isDegenerate()){
          num_degen++;
        }

        // Check whether the edge mesh exists - it has to!
        TMREdgeMesh *mesh = NULL;
        edges[i]->getMesh(&mesh);
        if (!mesh){
          fprintf(stderr,
                  "TMRFaceMesh error: Edge mesh does not exist\n");
        }

        // Get the number of points associated with the curve
        int npts;
        mesh->getMeshPoints(&npts, NULL, NULL);

        // Update the total number of points
        total_num_pts += npts-1;
      }
    }

    // The number of holes is equal to the number of loops-1. One loop
    // bounds the domain, the other loops cut out holes in the domain.
    // Note that the domain must be contiguous.
    int nholes = nloops-1;

    // All the boundary loops are closed, therefore, the total number
    // of segments is equal to the total number of points
    int nsegs = total_num_pts;

    // Allocate the points and the number of segments based on the
    // number of holes
    double *params = new double[ 2*(total_num_pts + nholes) ];
    int *segments = new int[ 2*nsegs ];

    // Start entering the points from the end of the last hole entry in
    // the parameter points array.
    int pt = 0;

    int init_loop_pt = 0; // What point value did this loop start on?
    int hole_pt = total_num_pts; // What hole are we on?

    // Set up the degenerate edges
    int *degen = NULL;
    if (num_degen > 0){
      degen = new int[ 2*num_degen ];
    }
    num_degen = 0;

    for ( int k = 0; k < nloops; k++ ){
      // Set the offset to the initial point/segment on this loop
      init_loop_pt = pt;

      // Get the curve information for this loop segment
      TMREdgeLoop *loop;
      face->getEdgeLoop(k, &loop);
      int nedges;
      TMREdge **edges;
      const int *dir;
      loop->getEdgeLoop(&nedges, &edges, &dir);

      for ( int i = 0; i < nedges; i++ ){
        // Retrieve the underlying curve mesh
        TMREdge *edge = edges[i];
        TMREdgeMesh *mesh = NULL;
        edge->getMesh(&mesh);

        // Get the mesh points corresponding to this curve
        int npts;
        const double *tpts;
        mesh->getMeshPoints(&npts, &tpts, NULL);

        // Find the point on the curve
        if (dir[i] > 0){
          for ( int j = 0; j < npts-1; j++ ){
            edge->getParamsOnFace(face, tpts[j], dir[i],
                                  &params[2*pt], &params[2*pt+1]);
            segments[2*pt] = pt;
            segments[2*pt+1] = pt+1;
            if (edge->isDegenerate()){
              degen[2*num_degen] = pt;
              degen[2*num_degen+1] = pt+1;
              num_degen++;
            }
            pt++;
          }
        }
        else {
          // Reverse the parameter values on the edge
          for ( int j = npts-1; j >= 1; j-- ){
            edge->getParamsOnFace(face, tpts[j], dir[i],
                                  &params[2*pt], &params[2*pt+1]);
            segments[2*pt] = pt;
            segments[2*pt+1] = pt+1;
            if (edge->isDegenerate()){
              degen[2*num_degen] = pt;
              degen[2*num_degen+1] = pt+1;
              num_degen++;
            }
            pt++;
          }
        }
      }

      // Close off the loop by connecting the segment back to the
      // initial loop point
      segments[2*(pt-1)+1] = init_loop_pt;

      // Compute the area enclosed by the loop. If the area is
      // positive, it is the domain boundary. If the area is negative,
      // we have a hole!  Note that this assumes that the polygon
      // creating the hole is not self-intersecting. (In reality we
      // compute twice the area since we omit the 1/2 factor.)
      double Area = 0.0;
      for ( int i = init_loop_pt; i < pt; i++ ){
        int s1 = segments[2*i];
        int s2 = segments[2*i+1];
        const double x1 = params[2*s1];
        const double y1 = params[2*s1+1];
        const double x2 = params[2*s2];
        const double y2 = params[2*s2+1];
        Area += (x1*y2 - x2*y1);
      }

      // Check the area constraint
      if (Area < 0.0){
        // This is a hole! Compute an approximate position for the hole.
        // Note that this may not work in all cases so beware.
        int s1 = segments[2*init_loop_pt];
        int s2 = segments[2*init_loop_pt+1];
        const double x1 = params[2*s1];
        const double y1 = params[2*s1+1];
        const double x2 = params[2*s2];
        const double y2 = params[2*s2+1];
        const double dx = x2 - x1;
        const double dy = y2 - y1;

        // This is arbitrary and won't work in general if we have a very
        // thin sliver for a hole...
        double frac = 0.01;

        // Set the average location for the hole
        params[2*hole_pt] = 0.5*(x1 + x2) + frac*dy;
        params[2*hole_pt+1] = 0.5*(y1 + y2) - frac*dx;

        // Increment the hole pointer
        hole_pt++;
      }
    }

    // Set the total number of fixed points. These are the points that
    // will not be smoothed and constitute the boundary nodes. Note
    // that the Triangularize class removes the holes from the domain
    // automatically.  The boundary points are guaranteed to be
    // ordered first.
    num_fixed_pts = total_num_pts - num_degen;

    if (source){
      // Create the source map of edges and keep track of their local
      // directions relative to the source surface
      std::map<TMREdge*, int> source_edges;
      for ( int k = 0; k < source->getNumEdgeLoops(); k++ ){
        TMREdgeLoop *loop;
        source->getEdgeLoop(k, &loop);

        // Get the number of edges/edges from the source loop
        int nedges;
        TMREdge **edges;
        const int *dir;
        loop->getEdgeLoop(&nedges, &edges, &dir);
        for ( int j = 0; j < nedges; j++ ){
          source_edges[edges[j]] = dir[j];
        }
      }

      // Create the target map of edges and keep track of their local
      // directions relative to the target surface
      std::map<TMREdge*, int> target_edges;
      for ( int k = 0; k < face->getNumEdgeLoops(); k++ ){
        TMREdgeLoop *loop;
        face->getEdgeLoop(k, &loop);

        // Get the number of edges/edges from the source loop
        int nedges;
        TMREdge **edges;
        const int *dir;
        loop->getEdgeLoop(&nedges, &edges, &dir);
        for ( int j = 0; j < nedges; j++ ){
          target_edges[edges[j]] = dir[j];
        }
      }

      // Keep track of the source-to-target edge and target-to-source
      // edge mappings as well as their relative orientations
      std::map<TMREdge*, int> target_edge_dir;
      std::map<TMREdge*, TMREdge*> source_to_target_edge;

      // Loop over the faces that are within the source volume
      int num_faces;
      TMRFace **faces;
      source_volume->getFaces(&num_faces, &faces, NULL);

      for ( int i = 0; i < num_faces; i++ ){
        // Check that this is not a target or source face
        if (faces[i] != source && faces[i] != face){
          // Find the source and target edge shared by the
          TMREdge *sedge = NULL, *tedge = NULL;
          int sdir = 0, tdir = 0;
          for ( int k = 0; k < faces[i]->getNumEdgeLoops(); k++ ){
            sedge = tedge = NULL;
            sdir = tdir = 0;

            // Get the edge loop
            TMREdgeLoop *loop;
            faces[i]->getEdgeLoop(k, &loop);

            // Get the number of edges/edges from the source loop
            int nedges;
            TMREdge **edges;
            const int *dir;
            loop->getEdgeLoop(&nedges, &edges, &dir);

            // Determine which edge is shared
            for ( int j = 0; j < nedges; j++ ){
              if (target_edges.count(edges[j]) > 0){
                tedge = edges[j];
                tdir = dir[j];
              }
              if (source_edges.count(edges[j]) > 0){
                sedge = edges[j];
                sdir = dir[j];
              }
            }

            if (sedge && tedge){
              break;
            }
          }

          if (sedge && tedge){
            // Compute the relative source-to-target directions
            int tmp = source_edges[sedge]*target_edges[tedge];
            target_edge_dir[tedge] = -sdir*tdir*tmp;

            // Source to target and target to source edges
            source_to_target_edge[sedge] = tedge;
          }
        }
      }

      // Now, count up the number of nodes that the target index must
      // be offset
      int target_offset = 0;
      std::map<TMREdge*, int> target_edge_offset;
      for ( int k = 0; k < face->getNumEdgeLoops(); k++ ){
        TMREdgeLoop *loop;
        face->getEdgeLoop(k, &loop);

        // Get the edges within this loop
        int nedges;
        TMREdge **edges;
        const int *dir;
        loop->getEdgeLoop(&nedges, &edges, &dir);

        for ( int i = 0; i < nedges; i++ ){
          // Retrieve the underlying curve mesh
          TMREdge *edge = edges[i];
          TMREdgeMesh *mesh = NULL;
          edge->getMesh(&mesh);

          // Get the mesh points corresponding to this curve
          int npts;
          const double *tpts;
          mesh->getMeshPoints(&npts, &tpts, NULL);

          target_edge_offset[edge] = target_offset;
          if (!edge->isDegenerate()){
            target_offset += npts-1;
          }
        }
      }

      // March through the sources loop, and compute the source to
      // target ordering
      int *source_to_target = new int[ num_fixed_pts ];
      int source_offset = 0;
      for ( int k = 0; k < source->getNumEdgeLoops(); k++ ){
        TMREdgeLoop *loop;
        source->getEdgeLoop(k, &loop);

        // Get the edges within this loop
        int nedges;
        TMREdge **edges;
        const int *dir;
        loop->getEdgeLoop(&nedges, &edges, &dir);

        for ( int i = 0; i < nedges; i++ ){
          // Retrieve the underlying mesh
          TMREdge *edge = edges[i];
          TMREdge *tedge = source_to_target_edge[edge];

          // Get the offset for the target edge
          int offset = target_edge_offset[tedge];

          // Retrieve the source mesh
          TMREdgeMesh *mesh = NULL;
          edge->getMesh(&mesh);

          // Get the mesh points corresponding to this curve
          int npts;
          mesh->getMeshPoints(&npts, NULL, NULL);

          // source:       target:
          // 0 -- 1 -> 2   6 <- 5 -- 4
          // |         |   |         |
          // 7         3   7         3
          // |         |   |         |
          // 6 <- 5 -- 4   0 -- 1 -> 2

          if (!edge->isDegenerate()){
            if (target_edge_dir[tedge] > 0){
              for ( int j = 0; j < npts-1; j++ ){
                source_to_target[source_offset + j] = offset + j;
              }
            }
            else {
              for ( int j = 0; j < npts-1; j++ ){
                source_to_target[source_offset + j] = offset + npts-1 - j;
              }

              // Get the previous target edge in the loop. This will give
              // the first number from the last edge loop.
              TMREdge *init_edge = NULL;
              if (i == 0){
                init_edge = source_to_target_edge[edges[nedges-1]];
              }
              else {
                init_edge = source_to_target_edge[edges[i-1]];
              }
              int init_offset = target_edge_offset[init_edge];
              source_to_target[source_offset] = init_offset;
            }
            // Increment the offset to the source
            source_offset += npts-1;
          }
        }
      }

      // Create the face mesh
      TMRFaceMesh *face_mesh;
      source->getMesh(&face_mesh);

      // Compute the total number of points
      mesh_type = face_mesh->mesh_type;
      num_points = face_mesh->num_points;
      num_fixed_pts = face_mesh->num_fixed_pts;
      num_quads = face_mesh->num_quads;
      num_tris = face_mesh->num_tris;

      // Allocate the array for the parametric locations
      pts = new double[ 2*num_points ];

      // Copy the points from around the boundaries
      for ( int i = 0; i < num_fixed_pts; i++ ){
        pts[2*i] = params[2*i];
        pts[2*i+1] = params[2*i+1];
      }

      // Compute a least squares transformation between the two
      // surfaces
      double N[16], A[4];
      double sc[2], tc[3];
      memset(N, 0, 16*sizeof(double));
      memset(A, 0, 4*sizeof(double));
      sc[0] = sc[1] = 0.0;
      tc[0] = tc[1] = 0.0;

      for ( int k = 0; k < num_fixed_pts; k++ ){
        sc[0] += face_mesh->pts[2*k];
        sc[1] += face_mesh->pts[2*k+1];
        tc[0] += pts[2*k];
        tc[1] += pts[2*k+1];
      }
      sc[0] = sc[0]/num_fixed_pts;
      sc[1] = sc[1]/num_fixed_pts;
      tc[0] = tc[0]/num_fixed_pts;
      tc[1] = tc[1]/num_fixed_pts;

      for ( int k = 0; k < num_fixed_pts; k++ ){
        double uS[2], uT[2];
        uS[0] = face_mesh->pts[2*k] - sc[0];
        uS[1] = face_mesh->pts[2*k+1] - sc[1];

        // Compute the source->target index number
        int kt = source_to_target[k];
        uT[0] = pts[2*kt] - tc[0];
        uT[1] = pts[2*kt+1] - tc[1];

        // Add the terms to the matrix/right-hand-side
        for ( int i = 0; i < 4; i++ ){
          for ( int j = 0; j < 4; j++ ){
            double B = uS[i % 2]*uS[j % 2];
            if ((i/2) == (j/2)){
              N[i + 4*j] += B;
            }
          }
          A[i] += uS[i % 2]*uT[i/2];
        }
      }

      // Factor the least-squares matrix and perform the transformation
      int ipiv[4];
      int n = 4, one = 1, info;
      TmrLAPACKdgetrf(&n, &n, N, &n, ipiv, &info);
      TmrLAPACKdgetrs("N", &n, &one, N, &n, ipiv, A, &n, &info);

      // Set the interior points based on the linear transformation
      for ( int k = num_fixed_pts; k < num_points; k++ ){
        double uS = face_mesh->pts[2*k] - sc[0];
        double vS = face_mesh->pts[2*k+1] - sc[1];
        pts[2*k] = A[0]*uS + A[1]*vS + tc[0];
        pts[2*k+1] = A[2]*uS + A[3]*vS + tc[1];
      }

      // Copy the quadrilateral mesh (if any)
      if (num_quads > 0){
        quads = new int[ 4*num_quads ];
        memcpy(quads, face_mesh->quads, 4*num_quads*sizeof(int));

        // Adjust the quadrilateral ordering at the boundary
        for ( int i = 0; i < 4*num_quads; i++ ){
          if (quads[i] < num_fixed_pts){
            quads[i] = source_to_target[quads[i]];
          }
        }

        // Flip the orientation of the quads to match the orientation
        // of the face
        if (source_dir < 0){
          for ( int i = 0; i < num_quads; i++ ){
            int tmp = quads[4*i+1];
            quads[4*i+1] = quads[4*i+3];
            quads[4*i+3] = tmp;
          }
        }
      }

      // Copy the triangular mesh (if any)
      if (num_tris > 0){
        tris = new int[ 3*num_tris ];
        memcpy(tris, face_mesh->tris, 3*num_tris*sizeof(int));

        // Adjust the triangle ordering at the boundary
        for ( int i = 0; i < 3*num_tris; i++ ){
          if (tris[i] < num_fixed_pts){
            tris[i] = source_to_target[tris[i]];
          }
        }

        // Flip the orientation of the triangles to match the
        // orientation of the face
        if (source_dir < 0){
          for ( int i = 0; i < num_tris; i++ ){
            int tmp = tris[3*i+1];
            tris[4*i+1] = tris[4*i+2];
            tris[4*i+2] = tmp;
          }
        }
      }

      // Free the data
      delete [] source_to_target;

      // Evaluate the points
      X = new TMRPoint[ num_points ];
      for ( int i = 0; i < num_points; i++ ){
        face->evalPoint(pts[2*i], pts[2*i+1], &X[i]);
      }

      if (num_quads > 0){
        // Smooth the copied mesh on the new surface
        int *pts_to_quad_ptr;
        int *pts_to_quads;
        computeNodeToElems(num_points, num_quads, 4, quads,
                           &pts_to_quad_ptr, &pts_to_quads);

        // Smooth the mesh using a local optimization of node locations
        quadSmoothing(options.num_smoothing_steps, num_fixed_pts,
                      num_points, pts_to_quad_ptr, pts_to_quads,
                      num_quads, quads, pts, X, face);

        // Free the connectivity information
        delete [] pts_to_quad_ptr;
        delete [] pts_to_quads;
      }
      else if (num_tris > 0){
        // Compute the triangle edges and neighbors in the dual mesh
        int num_tri_edges;
        int *tri_edges, *tri_neighbors, *dual_edges;
        computeTriEdges(num_points, num_tris, tris,
                        &num_tri_edges, &tri_edges,
                        &tri_neighbors, &dual_edges);

        // Smooth the resulting triangular mesh
        if (options.tri_smoothing_type == TMRMeshOptions::TMR_LAPLACIAN){
          laplacianSmoothing(options.num_smoothing_steps, num_fixed_pts,
                             num_tri_edges, tri_edges,
                             num_points, pts, X, face);
        }
        else {
          double alpha = 0.1;
          springSmoothing(options.num_smoothing_steps, alpha,
                          num_fixed_pts, num_tri_edges, tri_edges,
                          num_points, pts, X, face);
        }

        delete [] tri_edges;
        delete [] tri_neighbors;
        delete [] dual_edges;
      }
    }
    else if (mesh_type == TMR_STRUCTURED){
      // Use a straightforward interpolation technique to obtain the
      // structured parametric locations in terms of the boundary
      // point parametric locations. We do not perform checks here
      // since we already know that the surface has four edges and the
      // nodes on those edges can be used for a structured mesh

      // Get the first edge loop and the edges in the loop
      TMREdgeLoop *loop;
      face->getEdgeLoop(0, &loop);

      // Get the edges associated with the edge loop
      TMREdge **edges;
      loop->getEdgeLoop(NULL, &edges, NULL);

      // Get the number of nodes for the x/y edges
      int nx = 0, ny = 0;
      TMREdgeMesh *mesh;
      edges[0]->getMesh(&mesh);
      mesh->getMeshPoints(&nx, NULL, NULL);
      edges[1]->getMesh(&mesh);
      mesh->getMeshPoints(&ny, NULL, NULL);

      // Compute the total number of points
      num_points = nx*ny;
      num_quads = (nx-1)*(ny-1);

      // Create the connectivity information
      quads = new int[ 4*num_quads ];

      int *q = quads;
      for ( int j = 0; j < ny-1; j++ ){
        for ( int i = 0; i < nx-1; i++ ){
          // Compute the connectivity as if the element is on the
          // interior of the mesh
          q[0] = num_fixed_pts + (i-1) + (j-1)*(nx-2);
          q[1] = num_fixed_pts + i + (j-1)*(nx-2);
          q[2] = num_fixed_pts + i + j*(nx-2);
          q[3] = num_fixed_pts + (i-1) + j*(nx-2);

          // Adjust the ordering for the nodes on the boundary
          if (j == ny-2){
            q[2] = 2*nx + ny - 4 - i;
            q[3] = 2*nx + ny - 3 - i;
          }
          if (i == 0){
            q[0] = 2*nx + 2*ny - 4 - j;
            q[3] = 2*nx + 2*ny - 5 - j;
          }
          if (i == nx-2){
            q[1] = nx - 1 + j;
            q[2] = nx - 1 + j+1;
          }
          if (j == 0){
            q[0] = i;
            q[1] = i+1;
          }
          q += 4;
        }
      }

      // Now set the parametric locations on the interior
      pts = new double[ 2*num_points ];

      // Copy the points from around the boundaries
      for ( int i = 0; i < num_fixed_pts; i++ ){
        pts[2*i] = params[2*i];
        pts[2*i+1] = params[2*i+1];
      }

      // Use a transfinite interpolation to determine the parametric
      // points where the interior nodes should be placed.
      for ( int j = 1; j < ny-1; j++ ){
        for ( int i = 1; i < nx-1; i++ ){
          double u = 1.0*i/(nx-1);
          double v = 1.0*j/(ny-1);

          // Compute the weights on the corners
          double c1 = (1.0 - u)*(1.0 - v);
          double c2 = u*(1.0 - v);
          double c3 = u*v;
          double c4 = (1.0 - u)*v;

          // Compute the weights on the curves
          double w1 = (1.0 - v);
          double w2 = u;
          double w3 = v;
          double w4 = (1.0 - u);

          // New parametric point
          int p = num_fixed_pts + i-1 + (j-1)*(nx-2);

          // Boundary points that we're interpolating from
          int p1 = i;
          int p2 = nx-1 + j;
          int p3 = 2*nx + ny - 3 - i;
          int p4 = 2*nx + 2*ny - 4 - j;

          // Evaluate the parametric points based on the transfinite
          // interpolation
          pts[2*p] =
            ((w1*pts[2*p1] + w2*pts[2*p2] +
              w3*pts[2*p3] + w4*pts[2*p4]) -
             (c1*pts[0] + c2*pts[2*(nx-1)] +
              c3*pts[2*(nx+ny-2)] + c4*pts[2*(2*nx+ny-3)]));

          pts[2*p+1] =
            ((w1*pts[2*p1+1] + w2*pts[2*p2+1] +
              w3*pts[2*p3+1] + w4*pts[2*p4+1]) -
             (c1*pts[1] + c2*pts[2*(nx-1)+1] +
              c3*pts[2*(nx+ny-2)+1] + c4*pts[2*(2*nx+ny-3)+1]));
        }
      }

      // Allocate and evaluate the new physical point locations
      X = new TMRPoint[ num_points ];
      for ( int i = 0; i < num_points; i++ ){
        face->evalPoint(pts[2*i], pts[2*i+1], &X[i]);
      }
    }
    else {
      // Here mesh_type == TMR_TRIANGLE or TMR_UNSTRUCTURED

      // Create the triangularization class
      TMRTriangularize *tri =
        new TMRTriangularize(total_num_pts + nholes, params, nholes,
                             nsegs, segments, face);
      tri->incref();

      // Set the frontal quality factor
      tri->setFrontalQualityFactor(options.frontal_quality_factor);

      if (options.write_init_domain_triangle){
        char filename[256];
        sprintf(filename, "init_domain_triangle%d.vtk",
                face->getEntityId());
        tri->writeToVTK(filename);
      }

      // Create the mesh using the frontal algorithm
      tri->frontal(options, fs);

      // Free the degenerate triangles and reorder the mesh
      if (num_degen > 0){
        tri->removeDegenerateEdges(num_degen, degen);
        delete [] degen;
      }

      if (options.write_pre_smooth_triangle){
        char filename[256];
        sprintf(filename, "pre_smooth_triangle%d.vtk",
                face->getEntityId());
        tri->writeToVTK(filename);
      }

      // Extract the triangularization
      int ntris, *mesh_tris;
      tri->getMesh(&num_points, &ntris, &mesh_tris, &pts, &X);
      tri->decref();

      if (ntris == 0){
        fprintf(stderr,
                "TMRTriangularize warning: No triangles for mesh id %d\n",
                face->getEntityId());
      }

      if (ntris > 0){
        // Compute the triangle edges and neighbors in the dual mesh
        int num_tri_edges;
        int *tri_edges, *tri_neighbors, *dual_edges;
        int *node_to_tri_ptr, *node_to_tris;
        computeTriEdges(num_points, ntris, mesh_tris,
                        &num_tri_edges, &tri_edges,
                        &tri_neighbors, &dual_edges,
                        &node_to_tri_ptr, &node_to_tris);

        // Smooth the resulting triangular mesh
        if (options.tri_smoothing_type == TMRMeshOptions::TMR_LAPLACIAN){
          laplacianSmoothing(options.num_smoothing_steps, num_fixed_pts,
                             num_tri_edges, tri_edges,
                             num_points, pts, X, face);
        }
        else {
          double alpha = 0.1;
          springSmoothing(options.num_smoothing_steps, alpha,
                          num_fixed_pts, num_tri_edges, tri_edges,
                          num_points, pts, X, face);
        }

        if (options.write_post_smooth_triangle){
          char filename[256];
          sprintf(filename, "post_smooth_triangle%d.vtk",
                  face->getEntityId());
          writeTrisToVTK(filename, ntris, mesh_tris);
        }

        if (mesh_type == TMR_TRIANGLE){
          num_tris = ntris;
          tris = mesh_tris;

          // Free the allocated data
          delete [] tri_edges;
          delete [] tri_neighbors;
          delete [] dual_edges;
          delete [] node_to_tri_ptr;
          delete [] node_to_tris;
        }
        else { // mesh_type == TMR_UNSTRUCTURED
          // Recombine the mesh into a quadrilateral mesh
          if (ntris % 2 == 0){
            recombine(ntris, mesh_tris, tri_neighbors,
                      node_to_tri_ptr, node_to_tris,
                      num_tri_edges, dual_edges, &num_quads, &quads, options);
          }
          else {
            fprintf(stderr, "TMRFaceMesh error: Odd number of triangles, \
cannot perform recombination\n");
          }

          // Free the triangular mesh data
          delete [] mesh_tris;
          delete [] tri_edges;
          delete [] tri_neighbors;
          delete [] dual_edges;
          delete [] node_to_tri_ptr;
          delete [] node_to_tris;

          // Simplify the new quadrilateral mesh by removing points/quads
          // with poor quality/connectivity
          simplifyQuads();

          // Simplify a second time (for good measure)
          simplifyQuads();
        }

        if (options.write_pre_smooth_quad){
          char filename[256];
          sprintf(filename, "pre_smooth_quad%d.vtk",
                  face->getEntityId());
          writeToVTK(filename);
        }

        int *pts_to_quad_ptr;
        int *pts_to_quads;
        computeNodeToElems(num_points, num_quads, 4, quads,
                           &pts_to_quad_ptr, &pts_to_quads);

        // Smooth the mesh using a local optimization of node locations
        quadSmoothing(options.num_smoothing_steps, num_fixed_pts,
                      num_points, pts_to_quad_ptr, pts_to_quads,
                      num_quads, quads, pts, X, face);

        // Free the connectivity information
        delete [] pts_to_quad_ptr;
        delete [] pts_to_quads;

        if (options.write_post_smooth_quad){
          char filename[256];
          sprintf(filename, "post_smooth_quad%d.vtk",
                  face->getEntityId());
          writeToVTK(filename);
        }

        // Write out the dual of the final quadrilateral mesh
        if (options.write_quad_dual){
          int num_quad_edges;
          int *quad_edges;
          int *quad_neighbors, *quad_dual;
          computeQuadEdges(num_points, num_quads, quads,
                           &num_quad_edges, &quad_edges,
                           &quad_neighbors, &quad_dual);

          char filename[256];
          sprintf(filename, "quad_dual%d.vtk",
                  face->getEntityId());
          writeDualToVTK(filename, 4, num_quads, quads,
                         num_quad_edges, quad_dual, X);

          delete [] quad_edges;
          delete [] quad_neighbors;
          delete [] quad_dual;
        }
      }
    }

    // Free the parameter/segment information
    delete [] params;
    delete [] segments;
  }

  if (mpi_size > 1){
    // Broadcast the number of points to all the processors
    int temp[3];
    temp[0] = num_points;
    temp[1] = num_quads;
    temp[2] = num_fixed_pts;
    MPI_Bcast(temp, 3, MPI_INT, 0, comm);
    num_points = temp[0];
    num_quads = temp[1];
    num_fixed_pts = temp[2];

    if (mpi_rank != 0){
      pts = new double[ 2*num_points ];
      X = new TMRPoint[ num_points ];
      quads = new int[ 4*num_quads ];
    }

    // Broadcast the parametric locations and points
    MPI_Bcast(pts, 2*num_points, MPI_DOUBLE, 0, comm);
    MPI_Bcast(X, num_points, TMRPoint_MPI_type, 0, comm);
    MPI_Bcast(quads, 4*num_quads, MPI_INT, 0, comm);
  }
}

/*
  Retrieve the mesh points and parametric locations
*/
void TMRFaceMesh::getMeshPoints( int *_npts, const double **_pts,
                                 TMRPoint **_X ){
  if (_npts){ *_npts = num_points; }
  if (_pts){ *_pts = pts; }
  if (_X){ *_X = X; }
}

/*
  Set the node numbers internally
*/
int TMRFaceMesh::setNodeNums( int *num ){
  if (!vars){
    vars = new int[ num_points ];

    // Retrieve the boundary node numbers from the surface loops
    int pt = 0;
    for ( int k = 0; k < face->getNumEdgeLoops(); k++ ){
      // Get the curve information for this loop segment
      TMREdgeLoop *loop;
      face->getEdgeLoop(k, &loop);

      int nedges;
      TMREdge **edges;
      const int *dir;
      loop->getEdgeLoop(&nedges, &edges, &dir);

      for ( int i = 0; i < nedges; i++ ){
        // Retrieve the underlying curve mesh
        TMREdgeMesh *mesh = NULL;
        edges[i]->getMesh(&mesh);

        // Retrieve the variable numbers for this loop
        const int *edge_vars;
        int npts = mesh->getNodeNums(&edge_vars);

        if (edges[i]->isDegenerate()){
          vars[pt] = edge_vars[0];
        }
        else {
          // Find the point on the curve
          if (dir[i] > 0){
            for ( int j = 0; j < npts-1; j++, pt++ ){
              vars[pt] = edge_vars[j];
            }
          }
          else {
            for ( int j = npts-1; j >= 1; j--, pt++ ){
              vars[pt] = edge_vars[j];
            }
          }
        }
      }
    }

    // Now order the variables as they arrive
    for ( ; pt < num_points; pt++ ){
      vars[pt] = *num;
      (*num)++;
    }

    // Return the number of points that have been allocated
    return num_points - num_fixed_pts;
  }

  return 0;
}

/*
  Retrieve the mapping between the local connectivity and the global
  node numbers
*/
int TMRFaceMesh::getNodeNums( const int **_vars ){
  if (_vars){ *_vars = vars; }
  return num_points;
}

/*
  Get the number of fixed points that are not ordered by this surface
  mesh
*/
int TMRFaceMesh::getNumFixedPoints(){
  return num_fixed_pts;
}

/*
  Get the local quad connectivity
*/
int TMRFaceMesh::getQuadConnectivity( const int **_quads ){
  if (_quads){ *_quads = quads; }
  return num_quads;
}

/*
  Get the local triangle connectivity
*/
int TMRFaceMesh::getTriConnectivity( const int **_tris ){
  if (_tris){ *_tris = tris; }
  return num_tris;
}

/*
  Get the quadrilateral elemnet obtained by combining the triangles
  t1 and t2 together

        . --- .
      / |   /
    /   | /
  . --- .
*/
int TMRFaceMesh::getRecombinedQuad( const int triangles[],
                                    const int trineighbors[],
                                    int t1, int t2,
                                    int quad[] ){
  int fail = 0;

  // Find the common edge between the two tirangles
  const int shift[3] = {2, 0, 1};
  int e1 = 0, e2 = 0;
  for ( ; e1 < 3; e1++ ){
    if (trineighbors[3*t1+shift[e1]] == t2){ break; }
  }
  for ( ; e2 < 3; e2++ ){
    if (trineighbors[3*t2+shift[e2]] == t1){ break; }
  }

  if (e1 >= 3 || e2 >= 3){
    fail = 1;
    return fail;
  }

  // Order the triangle, omitting the common edge
  for ( int j = 0, i = 0; i < 3; i++, j++ ){
    quad[j] = triangles[3*t1+i];
    if (i == e1){
      j++;
    }
  }

  // Set the node contributed by the second triangle
  e2 += 2;
  if (e2 >= 3){ e2 -= 3; }
  quad[e1+1] = triangles[3*t2+e2];

  return fail;
}

/*
  Compute the quality of a quadrilateral element
*/
double TMRFaceMesh::computeQuadQuality( const int *quad,
                                        const TMRPoint *p ){
  // Compute the maximum of fabs(0.5*M_PI - alpha)
  double max_val = 0.0;

  for ( int k = 0; k < 4; k++ ){
    int prev = k-1;
    if (prev < 0){ prev = 3; }
    int next = k+1;
    if (next > 3){ next = 0; }

    TMRPoint a;
    a.x = p[quad[k]].x - p[quad[prev]].x;
    a.y = p[quad[k]].y - p[quad[prev]].y;
    a.z = p[quad[k]].z - p[quad[prev]].z;

    TMRPoint b;
    b.x = p[quad[next]].x - p[quad[k]].x;
    b.y = p[quad[next]].y - p[quad[k]].y;
    b.z = p[quad[next]].z - p[quad[k]].z;

    // Compute the internal angle between the
    double alpha = M_PI - acos(a.dot(b)/sqrt(a.dot(a)*b.dot(b)));
    double val = fabs(0.5*M_PI - alpha);
    if (val > max_val){
      max_val = val;
    }
  }

  // Compute the quality
  double eta = 1.0 - (2.0/M_PI)*max_val;
  if (eta < 0.0){
    eta = 0.0;
  }

  return eta;
}

/*
  Compute the recombined quality
*/
double TMRFaceMesh::computeRecombinedQuality( const int triangles[],
                                              const int tri_neighbors[],
                                              int t1, int t2,
                                              const TMRPoint *p ){
  // Find the combined quadrilateral from the two given triangles
  int quad[4];
  int fail = getRecombinedQuad(triangles, tri_neighbors, t1, t2, quad);
  if (fail){
    return 0.0;
  }

  return computeQuadQuality(quad, p);
}

/*
  Compute the quality of a quadrilateral element
*/
double TMRFaceMesh::computeTriQuality( const int *tri,
                                       const TMRPoint *p ){
  // Compute the maximum of fabs(M_PI/3 - alpha)
  double max_val = 0.0;

  for ( int k = 0; k < 3; k++ ){
    int prev = k-1;
    if (prev < 0){ prev = 2; }
    int next = k+1;
    if (next > 2){ next = 0; }

    TMRPoint a;
    a.x = p[tri[k]].x - p[tri[prev]].x;
    a.y = p[tri[k]].y - p[tri[prev]].y;
    a.z = p[tri[k]].z - p[tri[prev]].z;

    TMRPoint b;
    b.x = p[tri[next]].x - p[tri[k]].x;
    b.y = p[tri[next]].y - p[tri[k]].y;
    b.z = p[tri[next]].z - p[tri[k]].z;

    // Compute the internal angle
    double alpha = M_PI - acos(a.dot(b)/sqrt(a.dot(a)*b.dot(b)));
    double val = fabs(M_PI/3.0 - alpha);
    if (val > max_val){
      max_val = val;
    }
  }

  // Compute the quality
  double eta = 1.0 - (3.0/M_PI)*max_val;
  if (eta < 0.0){
    eta = 0.0;
  }

  return eta;
}

/*
  Recombine the triangulation into a quadrilateral mesh
*/
void TMRFaceMesh::recombine( int ntris, const int triangles[],
                             const int tri_neighbors[],
                             const int node_to_tri_ptr[],
                             const int node_to_tris[],
                             int num_edges, const int dual_edges[],
                             int *_num_quads, int **_new_quads,
                             TMRMeshOptions options ){
  // Allocate the reduced graph weights
  double *weights = new double[ num_edges ];
  int *graph_edges = new int[ 2*num_edges ];

  // Compute the weight associated with each edge by combputing the
  // recombined quality
  const double eps = 0.1;

  int edge_num = 0;
  for ( int i = 0; i < num_edges; i++ ){
    int t1 = dual_edges[2*i];
    int t2 = dual_edges[2*i+1];

    if (t1 >= 0 && t2 >= 0){
      // Compute the weight for this recombination
      double quality =
        computeRecombinedQuality(triangles, tri_neighbors,
                                 t1, t2, X);

      double weight = (1.0 - quality)*(1.0 + 1.0/(quality + eps));
      graph_edges[2*edge_num] = t1;
      graph_edges[2*edge_num+1] = t2;
      weights[edge_num] = weight;
      edge_num++;
    }
  }

  // Keep track of the number of real edges in the graph
  int num_real_edges = edge_num;

  // Determine the extra edges that connect along the boundaries.
  // between unique triangles that are not already adjacent to one
  // another.
  for ( int i = 0; i < ntris; i++ ){
    if (tri_neighbors[3*i] < 0 ||
        tri_neighbors[3*i+1] < 0 ||
        tri_neighbors[3*i+2] < 0){
      // Loop over the edges of the triangle
      for ( int j = 0; j < 3; j++ ){
        // We have a boundary triangle
        if (tri_neighbors[3*i+j] < 0){
          // The leading node that we are searching for
          int ij = tri_edge_nodes[j][1];
          int node = triangles[3*i + ij];

          // Search the triangles adjacent to node
          const int kpend = node_to_tri_ptr[node+1];
          for ( int kp = node_to_tri_ptr[node]; kp < kpend; kp++ ){
            int k = node_to_tris[kp];

            // If this is the same triangle, continue
            if (i == k){
              continue;
            }

            // If the triangles are actually also neighbors, continue
            if (tri_neighbors[3*k] == i ||
                tri_neighbors[3*k+1] == i ||
                tri_neighbors[3*k+2] == i){
              continue;
            }

            // Find the local node number shared between triangles k
            // and triangle i
            int kj = 0;
            if (triangles[3*k+1] == node){ kj = 1; }
            else if (triangles[3*k+2] == node){ kj = 2; }

            // Track from the node to the possible edges
            if (tri_neighbors[3*k + tri_node_edges[kj][0]] < 0 ||
                tri_neighbors[3*k + tri_node_edges[kj][1]] < 0){
              graph_edges[2*edge_num] = i;
              graph_edges[2*edge_num+1] = k;
              weights[edge_num] = 1.0;
              edge_num++;
            }
          }
        }
      }
    }
  }

  // Set the number of edges within the modified dual mesh
  int num_dual_edges = edge_num;

  // Write the dual mesh to a file
  if (options.write_dual_recombine){
    char filename[256];
    sprintf(filename, "dual_recombine%d.vtk",
            face->getEntityId());
    writeDualToVTK(filename, 3, ntris, triangles, num_dual_edges,
                   graph_edges, X);
  }

  // Perform the perfect matching
  int *match = new int[ ntris/2 ];
  int num_match = TMR_PerfectMatchGraph(ntris, num_dual_edges,
                                        graph_edges, weights, match);
  delete [] weights;

  // The quads formed from the original triangles
  int num_quads_from_tris = 0;

  // The new number of quadrilaterals - after every new quad is added
  int num_new_quads = 0;

  // New points array
  int num_new_points = 0;
  double *new_pts = NULL;
  TMRPoint *new_X = NULL;

  for ( int i = 0; i < num_match; i++ ){
    if (match[i] < num_real_edges){
      num_quads_from_tris++;
      num_new_quads++;
    }
    else {
      // We'll add two extra quads for each extra edge
      num_new_quads += 2;
      num_new_points++;
    }
  }

  // We'll be adding new points so allocate new arrays to handle the
  // new points that will be computed...
  if (num_new_points > 0){
    num_new_points += num_points;
    new_pts = new double[ 2*num_new_points ];
    new_X = new TMRPoint[ num_new_points ];
    memcpy(new_pts, pts, 2*num_points*sizeof(double));
    memcpy(new_X, X, num_points*sizeof(TMRPoint));
  }

  // Set the number of quadrilateral elements created in the mesh and
  // record the new quadrilateral element connectivity
  int *new_quads = new int[ 4*num_new_quads ];

  // Recombine the triangles into quadrilateral elements
  num_new_quads = 0;
  for ( int i = 0; i < num_match; i++ ){
    if (match[i] < num_real_edges){
      int t1 = graph_edges[2*match[i]];
      int t2 = graph_edges[2*match[i]+1];

      int fail = getRecombinedQuad(triangles, tri_neighbors, t1, t2,
                                   &new_quads[4*num_new_quads]);
      num_new_quads++;
      if (fail){
        fprintf(stderr,
                "TMRFaceMesh error: \
Quad %d from triangles %d and %d failed\n",
                i, t1, t2);
      }
    }
  }

  // Set the number of new points
  num_new_points = num_points;

  std::map<int, int> new_point_nums;

  // Add the triangles from edges along the boundary. There should
  // only be a handful of these guys...
  for ( int i = 0; i < num_match; i++ ){
    if (match[i] >= num_real_edges){
      // These triangles can only share one common node - the node
      // that we'll now duplicate and move to the interior of the
      // domain. Find the local node index for this shared node in
      // each triangle.
      int t1 = graph_edges[2*match[i]];
      int t2 = graph_edges[2*match[i]+1];

      // j1, j2 are the local nodal indices of the shared node between
      // triangles t1 and t2
      int j1 = 0, j2 = 0;
      for ( int flag = 0; !flag && j1 < 3; j1++ ){
        for ( j2 = 0; j2 < 3; j2++ ){
          if (triangles[3*t1 + j1] == triangles[3*t2 + j2]){
            flag = 1;
            break;
          }
        }
        if (flag){
          break;
        }
      }

      int boundary_pt = triangles[3*t1 + j1];

      // Go through all previous quadrilaterals and adjust the
      // ordering to reflect the duplicated node location
      for ( int k = 0; k < 4*num_new_quads; k++ ){
        if (new_quads[k] == boundary_pt){
          new_quads[k] = num_new_points;
        }
      }

      // Add the first triangle t1 - this triangle is gauranteed to
      // come first when circling the boundary in the CCW direction.
      int n1 = triangles[3*t1+((j1+1) % 3)];
      int n2 = triangles[3*t1+((j1+2) % 3)];
      if (new_point_nums.count(n1)){
        n1 = new_point_nums[n1];
      }
      if (new_point_nums.count(n2)){
        n2 = new_point_nums[n2];
      }
      new_quads[4*num_new_quads] = n1;
      new_quads[4*num_new_quads+1] = n2;
      new_quads[4*num_new_quads+2] = boundary_pt;
      new_quads[4*num_new_quads+3] = num_new_points;
      num_new_quads++;

      // Add the connectivity from the second triangle t2. This
      // triangle will always come second when circling the boundary.
      n1 = triangles[3*t2+((j2+1) % 3)];
      n2 = triangles[3*t2+((j2+2) % 3)];
      if (new_point_nums.count(n1)){
        n1 = new_point_nums[n1];
      }
      if (new_point_nums.count(n2)){
        n2 = new_point_nums[n2];
      }
      new_quads[4*num_new_quads] = n1;
      new_quads[4*num_new_quads+1] = n2;
      new_quads[4*num_new_quads+2] = num_new_points;
      new_quads[4*num_new_quads+3] = boundary_pt;
      num_new_quads++;

      // Add the new boundary points to the map
      new_point_nums[boundary_pt] = num_new_points;

      // Compute the new parameter location by taking the average of
      // the centroid locations for each triangle
      int count = 1;
      new_pts[2*num_new_points] = pts[2*boundary_pt];
      new_pts[2*num_new_points+1] = pts[2*boundary_pt+1];

      int kpend = node_to_tri_ptr[boundary_pt+1];
      for ( int kp = node_to_tri_ptr[boundary_pt]; kp < kpend; kp++ ){
        int k = node_to_tris[kp];
        if (k != t1 && k != t2){
          new_pts[2*num_new_points] += (pts[2*triangles[3*k]] +
                                        pts[2*triangles[3*k+1]] +
                                        pts[2*triangles[3*k+2]])/3.0;
          new_pts[2*num_new_points+1] += (pts[2*triangles[3*k]+1] +
                                          pts[2*triangles[3*k+1]+1] +
                                          pts[2*triangles[3*k+2]+1])/3.0;
          count++;
        }
      }
      if (count > 1){
        new_pts[2*num_new_points] = new_pts[2*num_new_points]/count;
        new_pts[2*num_new_points+1] = new_pts[2*num_new_points+1]/count;
      }
      face->evalPoint(new_pts[2*num_new_points],
                      new_pts[2*num_new_points+1], &new_X[num_new_points]);

      // Increment the number of new points
      num_new_points++;
    }
  }

  delete [] graph_edges;
  delete [] match;

  if (new_pts){
    num_points = num_new_points;
    delete [] pts;
    delete [] X;
    pts = new_pts;
    X = new_X;
  }

  // Set the quads/output
  *_num_quads = num_new_quads;
  *_new_quads = new_quads;
}

/*
  This code performs several topological improvements to try and avoid
  poor mesh quality.

  First, the code attempts to remove topological triangles on
  boundaries where adjacent boundary edges are included in the same
  quad. Adjacent boundary edges in the same quadrilaterla are not
  always bad (for instance when a quad is at a corner), so we check
  whether to adjust the connectivity based on the angle between the
  two edges. Note that this will only work if the second quadrilateral
  does not have an edge on the boundary and requires that the mesh be
  smoothed afterwards.

           x                            x
         /   \                        / / \
        /     \                     /  /   \
       x       x                   x  /     x
     /   \     /                  /   /     /
    /     \   /         ===>     /   /     /
   /       \ /                  /   /     /
  x -- x -- x                  x -- x -- x

  Next, the code identifies and removes adjacent quadrilaterals that
  look like this:

  x --- x             x ---- x
  | \   |             |      |
  |  x  |     ===>    |      |
  |    \|             |      |
  x --- x             x ---- x

  Finally, the code finds and removes quads that look like this:

  x ---- x ---- x         x ---- x ---- x
  |     / \     |         |      |      |
  |    /   \    |         |      |      |
  x - x     x - x   ===>  x ---- x ---- x
  |    \   /    |         |      |      |
  |     \ /     |         |      |      |
  x ---- x ---- x         x ---- x ---- x
*/
void TMRFaceMesh::simplifyQuads(){
  // Compute the node -> quad information
  int *ptr, *pts_to_quads;
  computeNodeToElems(num_points, num_quads, 4, quads,
                     &ptr, &pts_to_quads);

  // First, remove the nodes that are only referred to twice,
  // not on the boundary
  int *new_pt_nums = new int[ num_points ];
  memset(new_pt_nums, 0, num_points*sizeof(int));

  for ( int i = num_fixed_pts; i < num_points; i++ ){
    if (ptr[i+1]-ptr[i] == 2){
      // Retrieve the pointers to the latest quadrilatral
      int q1 = pts_to_quads[ptr[i]];
      int q2 = pts_to_quads[ptr[i]+1];
      int *quad1 = &quads[4*q1];
      int *quad2 = &quads[4*q2];

      // If any of the other points have been deleted, then
      // skip combining this quadrilateral. This call could
      // probably be repeated to remove this element, but this
      // algorithm cannot handle it.
      int skip_me = 0;
      for ( int k = 0; k < 4; k++ ){
        if (quad1[k] < 0 ||
            quad2[k] < 0 ||
            new_pt_nums[quad1[k]] < 0 ||
            new_pt_nums[quad2[k]] < 0){
          skip_me = 1;
          break;
        }
      }

      if (!skip_me){
        // Find the common node between quads
        int k1 = 0, k2 = 0;
        while (quad1[k1] != i) k1++;
        while (quad2[k2] != i) k2++;

        // Adjust k2 so that it points to the node that
        // will be inserted into the connectivity for the
        // first quadrilateral
        k2 += 2;
        if (k2 >= 4){ k2 -= 4; }
        quad1[k1] = quad2[k2];

        // Set the point
        quad2[0] = quad2[1] = quad2[2] = quad2[3] = -1;

        // Set the points that we're going to eliminate
        new_pt_nums[i] = -1;
      }
    }
  }

  // Remove the triangle quadrilaterals
  for ( int i = 0; i < num_quads; i++ ){
    for ( int j1 = 0; j1 < 2; j1++ ){
      // Find the local node numbers that are on opposite
      // sides of the quadrilateral
      int j2 = j1 + 2;

      // Compute the points that are on opposite sides
      // of the quadrilateral
      int p1 = quads[4*i+j1];
      int p2 = quads[4*i+j2];
      if (p1 >= 0 && p2 >= 0 &&
          p1 >= num_fixed_pts &&
          p2 >= num_fixed_pts &&
          ptr[p1+1] - ptr[p1] == 3 &&
          ptr[p2+1] - ptr[p2] == 3){
        // Check whether any of the quadrilaterals which
        // touch either p1 or p2 have negative indices
        int flag = 0;
        for ( int kp = ptr[p1]; kp < ptr[p1+1]; kp++ ){
          int q = pts_to_quads[kp];
          if (quads[4*q] < 0 ||
              new_pt_nums[quads[4*q]] < 0 ||
              new_pt_nums[quads[4*q+1]] < 0 ||
              new_pt_nums[quads[4*q+2]] < 0 ||
              new_pt_nums[quads[4*q+3]] < 0){
            flag = 1;
            break;
          }
        }
        if (flag){ break; }
        for ( int kp = ptr[p2]; kp < ptr[p2+1]; kp++ ){
          int q = pts_to_quads[kp];
          if (quads[4*q] < 0 ||
              new_pt_nums[quads[4*q]] < 0 ||
              new_pt_nums[quads[4*q+1]] < 0 ||
              new_pt_nums[quads[4*q+2]] < 0 ||
              new_pt_nums[quads[4*q+3]] < 0){
            flag = 1;
            break;
          }
        }
        if (flag){ break; }

        // Figure out which node to eliminate
        new_pt_nums[p1] = -1;

        // Set the quadrilaterals that refer to the node p1
        // now instead to refer to the node p2
        for ( int kp = ptr[p1]; kp < ptr[p1+1]; kp++ ){
          int q = pts_to_quads[kp];
          for ( int k = 0; k < 4; k++ ){
            if (quads[4*q+k] == p1){
              quads[4*q+k] = p2;
            }
          }
        }

        // Now eliminate the quad
        quads[4*i] = quads[4*i+1] = quads[4*i+2] = quads[4*i+3] = -1;

        // Quit this quadrilateral
        break;
      }
    }
  }

  // Free the pointer/quad pointer data
  delete [] ptr;
  delete [] pts_to_quads;

  // Remove the points/parameters that have been eliminated
  int pt_num = 0;
  for ( ; pt_num < num_fixed_pts; pt_num++ ){
    new_pt_nums[pt_num] = pt_num;
  }
  for ( int i = num_fixed_pts; i < num_points; i++ ){
    if (new_pt_nums[i] >= 0){
      if (i != pt_num){
        // Copy over the points
        pts[2*pt_num] = pts[2*i];
        pts[2*pt_num+1] = pts[2*i+1];
        X[pt_num] = X[i];
      }
      new_pt_nums[i] = pt_num;
      pt_num++;
    }
  }

  // Set the new number of points
  num_points = pt_num;

  // Set the new connectivity, overwriting the old connectivity
  int quad_num = 0;
  for ( int i = 0; i < num_quads; i++ ){
    if (quads[4*i] >= 0){
      quads[4*quad_num] = new_pt_nums[quads[4*i]];
      quads[4*quad_num+1] = new_pt_nums[quads[4*i+1]];
      quads[4*quad_num+2] = new_pt_nums[quads[4*i+2]];
      quads[4*quad_num+3] = new_pt_nums[quads[4*i+3]];
      quad_num++;
    }
  }

  // Set the new number of quadrilaterals
  num_quads = quad_num;

  // Free the new point numbers
  delete [] new_pt_nums;
}

/*
  Add the quad quality
*/
void TMRFaceMesh::addMeshQuality( int nbins, int bins[] ){
  for ( int i = 0; i < num_quads; i++ ){
    double quality = computeQuadQuality(&quads[4*i], X);

    int k = 0;
    for ( ; k < nbins; k++ ){
      if (quality < 1.0*(k+1)/nbins){
        break;
      }
    }
    if (k == nbins){
      k = nbins-1;
    }
    bins[k]++;
  }

  for ( int i = 0; i < num_tris; i++ ){
    double quality = computeTriQuality(&tris[3*i], X);

    int k = 0;
    for ( ; k < nbins; k++ ){
      if (quality < 1.0*(k+1)/nbins){
        break;
      }
    }
    if (k == nbins){
      k = nbins-1;
    }
    bins[k]++;
  }
}

/*
  Print the quadrilateral quality
*/
void TMRFaceMesh::printMeshQuality(){
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins*sizeof(int));
  addMeshQuality(nbins, bins);

  for ( int i = 0; i < nbins; i++ ){
    total += bins[i];
  }

  printf("Quality   # elements   percentage\n");
  for ( int k = 0; k < nbins; k++ ){
    printf("< %.2f    %10d   %10.3f\n",
           1.0*(k+1)/nbins, bins[k], 100.0*bins[k]/total);
  }
}

/*
  Print the triangle quality
*/
void TMRFaceMesh::printTriQuality( int ntris,
                                   const int triangles[] ){
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins*sizeof(int));
  for ( int i = 0; i < ntris; i++ ){
    double quality = computeTriQuality(&triangles[3*i], X);

    int k = 0;
    for ( ; k < nbins; k++ ){
      if (quality < 1.0*(k+1)/nbins){
        break;
      }
    }
    if (k == nbins){
      k = nbins-1;
    }
    bins[k]++;
  }

  for ( int i = 0; i < nbins; i++ ){
    total += bins[i];
  }

  printf("Quality   # elements   percentage\n");
  for ( int k = 0; k < nbins; k++ ){
    printf("< %.2f    %10d   %10.3f\n",
           1.0*(k+1)/nbins, bins[k], 100.0*bins[k]/total);
  }
}

/*
  Write the quadrilateral mesh to a VTK file
*/
void TMRFaceMesh::writeToVTK( const char *filename ){
  if (num_quads > 0){
    FILE *fp = fopen(filename, "w");
    if (fp){
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

      // Write out the points
      fprintf(fp, "POINTS %d float\n", num_points);
      for ( int k = 0; k < num_points; k++ ){
        fprintf(fp, "%e %e %e\n", X[k].x, X[k].y, X[k].z);
      }

      // Write out the cell values
      fprintf(fp, "\nCELLS %d %d\n", num_quads, 5*num_quads);
      for ( int k = 0; k < num_quads; k++ ){
        fprintf(fp, "4 %d %d %d %d\n", quads[4*k], quads[4*k+1],
                quads[4*k+2], quads[4*k+3]);
      }

      // All quadrilaterals
      fprintf(fp, "\nCELL_TYPES %d\n", num_quads);
      for ( int k = 0; k < num_quads; k++ ){
        fprintf(fp, "%d\n", 9);
      }

      // Print out the rest as fields one-by-one
      fprintf(fp, "CELL_DATA %d\n", num_quads);
      fprintf(fp, "SCALARS quality float 1\n");
      fprintf(fp, "LOOKUP_TABLE default\n");
      for ( int k = 0; k < num_quads; k++ ){
        fprintf(fp, "%e\n", computeQuadQuality(&quads[4*k], X));
      }

      fclose(fp);
    }
  }
  else if (num_tris > 0){
    writeTrisToVTK(filename, num_tris, tris);
  }
}

/*
  Write out the segments to a VTK file
*/
void TMRFaceMesh::writeSegmentsToVTK( const char *filename,
                                      int npts, const double *params,
                                      int nsegs, const int segs[] ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts);
    for ( int k = 0; k < npts; k++ ){
      fprintf(fp, "%e %e 0.0\n", params[2*k], params[2*k+1]);
    }

    // Write out the cell connectivity
    fprintf(fp, "\nCELLS %d %d\n", nsegs, 3*nsegs);
    for ( int k = 0; k < nsegs; k++ ){
      fprintf(fp, "2 %d %d\n", segs[2*k], segs[2*k+1]);
    }

    // Write out the cell types
    fprintf(fp, "\nCELL_TYPES %d\n", nsegs);
    for ( int k = 0; k < nsegs; k++ ){
      fprintf(fp, "%d\n", 3);
    }
    fclose(fp);
  }
}

/*
  Write the output to a VTK file
*/
void TMRFaceMesh::writeTrisToVTK( const char *filename,
                                  int ntris, const int triangles[] ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write out the points
    fprintf(fp, "POINTS %d float\n", num_points);
    for ( int k = 0; k < num_points; k++ ){
      fprintf(fp, "%e %e %e\n", X[k].x, X[k].y, X[k].z);
    }

    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", ntris, 4*ntris);
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "3 %d %d %d\n",
              triangles[3*k], triangles[3*k+1], triangles[3*k+2]);
    }

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", ntris);
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "%d\n", 5);
    }

        // Print out the rest as fields one-by-one
    fprintf(fp, "CELL_DATA %d\n", ntris);
    fprintf(fp, "SCALARS quality float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "%e\n", computeTriQuality(&triangles[3*k], X));
    }

    fclose(fp);
  }
}

/*
  Write out the dual mesh for visualization
*/
void TMRFaceMesh::writeDualToVTK( const char *filename,
                                  int nodes_per_elem,
                                  int nelems, const int elems[],
                                  int num_dual_edges,
                                  const int dual_edges[],
                                  const TMRPoint *p ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write out the points
    fprintf(fp, "POINTS %d float\n", nelems);
    if (nodes_per_elem == 3){
      for ( int k = 0; k < nelems; k++ ){
        fprintf(fp, "%e %e %e\n",
                1.0/3.0*(p[elems[3*k]].x + p[elems[3*k+1]].x +
                         p[elems[3*k+2]].x),
                1.0/3.0*(p[elems[3*k]].y + p[elems[3*k+1]].y +
                         p[elems[3*k+2]].y),
                1.0/3.0*(p[elems[3*k]].z + p[elems[3*k+1]].z +
                         p[elems[3*k+2]].z));
      }
    }
    else { // nodes_per_elem == 4
      for ( int k = 0; k < nelems; k++ ){
        fprintf(fp, "%e %e %e\n",
                0.25*(p[elems[4*k]].x + p[elems[4*k+1]].x +
                      p[elems[4*k+2]].x + p[elems[4*k+3]].x),
                0.25*(p[elems[4*k]].y + p[elems[4*k+1]].y +
                      p[elems[4*k+2]].y + p[elems[4*k+3]].y),
                0.25*(p[elems[4*k]].z + p[elems[4*k+1]].z +
                      p[elems[4*k+2]].z + p[elems[4*k+3]].z));
      }
    }

    // Count up the number of non-degenerate dual edges
    int n = 0;
    for ( int k = 0; k < num_dual_edges; k++ ){
      if (dual_edges[2*k] >= 0 && dual_edges[2*k+1] >= 0){
        n++;
      }
    }

    // Write out the cell connectivity
    fprintf(fp, "\nCELLS %d %d\n", n, 3*n);
    for ( int k = 0; k < num_dual_edges; k++ ){
      if (dual_edges[2*k] >= 0 && dual_edges[2*k+1] >= 0){
        fprintf(fp, "2 %d %d\n", dual_edges[2*k], dual_edges[2*k+1]);
      }
    }

    // Write out the cell types
    fprintf(fp, "\nCELL_TYPES %d\n", n);
    for ( int k = 0; k < n; k++ ){
      fprintf(fp, "%d\n", 3);
    }
    fclose(fp);
  }
}

/*
  Try to create a volume mesh based on the structured/unstructured
  surface meshes
*/
TMRVolumeMesh::TMRVolumeMesh( MPI_Comm _comm,
                              TMRVolume *_volume ){
  comm = _comm;
  volume = _volume;
  volume->incref();

  // Set the number of loops on the source face
  num_face_loops = 0;
  face_loops = NULL;
  face_loop_ptr = NULL;
  face_loop_dir = NULL;
  face_loop_edge_count = NULL;

  // Set the number of points through-thickness
  num_depth_pts = -1;

  // Set the points in the mesh
  num_points = 0;
  num_hex = 0;
  num_tet = 0;

  // Set all the pointer to null
  X = NULL;
  hex = NULL;
  tet = NULL;
  vars = NULL;

  // Set the additional connectivity information that is required for
  // the volume mesh.
  source = target = NULL;
  source_dir = target_dir = 0;
}

TMRVolumeMesh::~TMRVolumeMesh(){
  volume->decref();
  if (X){ delete [] X; }
  if (hex){ delete [] hex; }
  if (tet){ delete [] tet; }
  if (vars){ delete [] vars; }
  if (face_loops){
    for ( int k = 0; k < face_loop_ptr[num_face_loops]; k++ ){
      face_loops[k]->decref();
    }
    delete [] face_loop_ptr;
    delete [] face_loops;
    delete [] face_loop_dir;
    delete [] face_loop_edge_count;
  }
}

/*
  Retrieve the mesh points
*/
void TMRVolumeMesh::getMeshPoints( int *_npts, TMRPoint **_X ){
  if (_npts){ *_npts = num_points; }
  if (_X){ *_X = X; }
}

/*
  Retrieve the local connectivity for this volume mesh
*/
int TMRVolumeMesh::getHexConnectivity( const int **_hex ){
  if (_hex){ *_hex = hex; }
  return num_hex;
}

/*
  Retrieve the local connectivity for this volume mesh
*/
int TMRVolumeMesh::getTetConnectivity( const int **_tet ){
  if (_tet){ *_tet = tet; }
  return num_tet;
}

/*
  Write the volume mesh to a VTK file
*/
void TMRVolumeMesh::writeToVTK( const char *filename ){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  if (mpi_rank == 0 && hex){
    // Write out the connectivity of the mesh to a temp vtk file
    FILE *fp = fopen("volume-mesh.vtk", "w");
    if (fp){
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

      // Write out the points
      fprintf(fp, "POINTS %d float\n", num_points);
      for ( int k = 0; k < num_points; k++ ){
        fprintf(fp, "%e %e %e\n", X[k].x, X[k].y, X[k].z);
      }

      // Write out the cell values
      fprintf(fp, "\nCELLS %d %d\n", num_hex, 9*num_hex);
      for ( int k = 0; k < num_hex; k++ ){
        fprintf(fp, "8 %d %d %d %d %d %d %d %d\n",
                hex[8*k], hex[8*k+1], hex[8*k+2], hex[8*k+3],
                hex[8*k+4], hex[8*k+5], hex[8*k+6], hex[8*k+7]);
      }

      // All hex
      fprintf(fp, "\nCELL_TYPES %d\n", num_hex);
      for ( int k = 0; k < num_hex; k++ ){
        fprintf(fp, "%d\n", 12);
      }
      fclose(fp);
    }
  }
}

/*
  Mesh the volume
*/
int TMRVolumeMesh::mesh( TMRMeshOptions options ){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  if (options.mesh_type_default == TMR_TRIANGLE){
    return tetMesh(options);
  }

  // Keep track of whether the mesh has failed at any time. Try and
  // print out a helpful message.
  int mesh_fail = 0;

  // Get the faces associated with the volume
  int num_faces;
  TMRFace **faces;
  const int *face_dir;
  volume->getFaces(&num_faces, &faces, &face_dir);

  // Set integers for each face to determine whether it is a source or
  // a target face or it is connected to one of them or whether it is
  // structured.
  int *f = new int[ num_faces ];
  memset(f, 0, num_faces*sizeof(int));

  // Set the pointers for the target/source faces and their directions
  // relative to the interior of the hexahedral volume
  target = NULL;
  target_dir = 1;
  source = NULL;
  source_dir = 1;

  // Loop over all the faces to find the sources/targets
  for ( int i = 0; i < num_faces; i++ ){
    TMRFace *src;
    faces[i]->getSource(NULL, NULL, &src);
    if (src){
      // The natural source orientation should point in to the volume,
      // otherwise its orientation must be flipped.
      target = faces[i];
      target_dir = -face_dir[i];

      // Find the source face
      source = src;
      f[i] = 1;
      for ( int j = 0; j < num_faces; j++ ){
        if (faces[j] == source){
          f[j] = 1;

          // The natural orientation for the target face should point out
          // of the volume, otherwise its orientation must be flipped.
          source_dir = face_dir[j];
        }
      }

      // Break here. The meshing algorithm only works with one
      // source/target face pair per volume.
      break;
    }
  }

  // Check if the remaining surface meshes are structured
  for ( int i = 0; i < num_faces; i++ ){
    if (!f[i]){
      TMRFaceMesh *mesh;
      faces[i]->getMesh(&mesh);
      if (!mesh){
        fprintf(stderr,
                "TMRVolumeMesh error: No mesh associated with face %d\n",
                faces[i]->getEntityId());
        mesh_fail = 1;
      }
      else if (mesh->getMeshType() != TMR_STRUCTURED){
        fprintf(stderr,
                "TMRVolumeMesh error: \
Through-thickness meshes must be structured\n");
        fprintf(stderr,
                "Try setting source-target relations on edges and surfaces\n");
        mesh_fail = 1;
      }
    }
  }

  // Check that the source mesh is not a triangular mesh
  TMRFaceMesh *source_mesh;
  source->getMesh(&source_mesh);
  if (source_mesh->getMeshType() == TMR_TRIANGLE){
    fprintf(stderr,
            "TMRVolumeMesh error: Cannot extrude a triangluar mesh\n");
    mesh_fail = 1;
  }

  // Free the f pointer
  delete [] f;

  if (mesh_fail){
    return mesh_fail;
  }

  // Each face that is not either the target or the source, must have
  // only four edges and must be structured. Two of the parallel edges
  // associated with these faces should touch the target and source
  // face, while remaining edges must be parallel and have the same
  // number of nodes. Furthermore, the one parallel edge will be
  // shared by the next face. We loop over the source edge loops and
  // find the connecting faces.

  // Get the number of edge loops
  num_face_loops = source->getNumEdgeLoops();

  // Count up the total number of faces
  face_loop_ptr = new int[ num_face_loops+1 ];

  int count = 0;
  for ( int k = 0; k < num_face_loops; k++ ){
    TMREdgeLoop *source_loop;
    source->getEdgeLoop(k, &source_loop);

    // Get the number of edges for this loop
    int nedges;
    source_loop->getEdgeLoop(&nedges, NULL, NULL);
    count += nedges;
  }

  face_loops = new TMRFace*[ count ];
  face_loop_dir = new int[ count ];
  face_loop_edge_count = new int[ count ];

  // Set the number of points through the depth
  num_depth_pts = -1;

  // Set the first entry in the loop pointer
  face_loop_ptr[0] = 0;

  for ( int k = 0; k < num_face_loops; k++ ){
    // Set the next entry in the loop pointer
    face_loop_ptr[k+1] = face_loop_ptr[k];

    // The target and source edge loops
    TMREdgeLoop *source_loop;
    source->getEdgeLoop(k, &source_loop);

    // Get the edges associated with the source loop
    int nedges;
    TMREdge **edges;
    source_loop->getEdgeLoop(&nedges, &edges, NULL);

    for ( int j = 0; j < nedges; j++ ){
      // Search for the other face object that shares the edge
      // object edges[j] with the source face object
      TMRFace *face = NULL;
      int fdir = 0;
      for ( int i = 0; i < num_faces; i++ ){
        if (!(faces[i] == target || faces[i] == source)){
          // Get the edge loop associated with the face
          TMREdgeLoop *loop;
          faces[i]->getEdgeLoop(0, &loop);

          // Get the edge loops associated with face i
          int n;
          TMREdge **e;
          loop->getEdgeLoop(&n, &e, NULL);

          // Does this edge loop contain edges[j]
          for ( int ii = 0; ii < n; ii++ ){
            if (e[ii] == edges[j]){
              face = faces[i];
              fdir = face_dir[i];
              break;
            }
          }

          if (face){
            break;
          }
        }
      }

      // This face is not the target or source face and therefore must
      // be structured.
      TMREdgeLoop *loop;
      face->getEdgeLoop(0, &loop);

      // Get the edge loop for this face
      TMREdge **face_edges;
      loop->getEdgeLoop(NULL, &face_edges, NULL);

      TMREdgeMesh *mesh = NULL;
      if (face_edges[0] == edges[j] ||
          face_edges[2] == edges[j]){
        face_edges[1]->getMesh(&mesh);
      }
      else if (face_edges[1] == edges[j] ||
               face_edges[3] == edges[j]){
        face_edges[0]->getMesh(&mesh);
      }

      // Get the number of mesh points
      int npts = -1;
      mesh->getMeshPoints(&npts, NULL, NULL);

      // If this is the first time finding the number of points along
      // the depth of this volume, record the number of points.
      // Otherwise, verify that the number of points through the depth
      // of the volume is consistent.
      if (num_depth_pts < 0){
        num_depth_pts = npts;
      }
      else if (num_depth_pts != npts){
        fprintf(stderr,
                "TMRVolumeMesh error: \
Inconsistent number of edge points through-thickness %d != %d\n",
                num_depth_pts, npts);
        mesh_fail = 1;
      }

      // Find the number of nodes on the edges that are parallel to
      // the edge loops surrounding the face
      if (face_edges[0] == edges[j] ||
          face_edges[2] == edges[j]){
        face_edges[0]->getMesh(&mesh);
      }
      else if (face_edges[1] == edges[j] ||
               face_edges[3] == edges[j]){
        face_edges[1]->getMesh(&mesh);
      }

      // Get the number of mesh points from the parallel edge
      mesh->getMeshPoints(&npts, NULL, NULL);

      // Set the number of points
      face_loop_edge_count[face_loop_ptr[k+1]] = npts;

      // Record the face object for later use (during ordering) and
      // incref the face object
      face_loops[face_loop_ptr[k+1]] = face;
      face->incref();
      face_loop_dir[face_loop_ptr[k+1]] = fdir;
      face_loop_ptr[k+1]++;
    }
  }

  // Get the information for the target surface
  TMRFaceMesh *mesh;
  target->getMesh(&mesh);

  // Number of points in the quadrilateral mesh on the surface
  int num_quad_pts = 0;
  TMRPoint *Xtarget;
  mesh->getMeshPoints(&num_quad_pts, NULL, &Xtarget);

  // Get information for the source surface
  source->getMesh(&mesh);

  // Points on the target surface
  TMRPoint *Xbot;
  mesh->getMeshPoints(NULL, NULL, &Xbot);

  // Get the local connectivity on the source surface
  const int *quads;
  int num_quads = mesh->getQuadConnectivity(&quads);

  // The number of hexahedral elements in the mesh
  num_hex = (num_depth_pts-1)*num_quads;
  num_points = num_depth_pts*num_quad_pts;
  hex = new int[ 8*num_hex ];
  X = new TMRPoint[ num_points ];

  // Flip the ordering in the hexahedral elements
  const int flip[] = {0, 3, 2, 1};

  int *h = hex;
  for ( int j = 0; j < num_depth_pts-1; j++ ){
    for ( int i = 0; i < num_quads; i++ ){
      // Set the quadrilateral points in the base layer
      if (source_dir > 0){
        for ( int k = 0; k < 4; k++ ){
          h[k] = j*num_quad_pts + quads[4*i+flip[k]];
        }
        for ( int k = 0; k < 4; k++ ){
          h[4+k] = (j+1)*num_quad_pts + quads[4*i+flip[k]];
        }
      }
      else {
        for ( int k = 0; k < 4; k++ ){
          h[k] = j*num_quad_pts + quads[4*i+k];
        }
        for ( int k = 0; k < 4; k++ ){
          h[4+k] = (j+1)*num_quad_pts + quads[4*i+k];
        }
      }
      h += 8;
    }
  }

  // Set the new coordinates within the hexahedral mesh
  TMRPoint *x = X;
  for ( int j = 0; j < num_depth_pts; j++ ){
    double u = 1.0*j/(num_depth_pts-1);
    for ( int i = 0; i < num_quad_pts; i++ ){
      x[0].x = (1.0 - u)*Xbot[i].x + u*Xtarget[i].x;
      x[0].y = (1.0 - u)*Xbot[i].y + u*Xtarget[i].y;
      x[0].z = (1.0 - u)*Xbot[i].z + u*Xtarget[i].z;
      x += 1;
    }
  }

  return mesh_fail;
}

/*
  Create a tetrahedral mesh
*/
int TMRVolumeMesh::tetMesh( TMRMeshOptions options ){
#ifdef TMR_USE_NETGEN
  nglib::Ng_Init();
  nglib::Ng_Mesh *m = nglib::Ng_NewMesh();

  // First count up the total number of unique points that bound this
  // volume that is about to be meshed
  int num_boundary_pts = 0;

  // Get the faces associated with the volume
  int num_faces;
  TMRFace **faces;
  const int *face_dir;
  volume->getFaces(&num_faces, &faces, &face_dir);

  // Count up the number of boundary nodes for this volume
  for ( int i = 0; i < ; i++ ){
    double pt[3] = {0.0, 0.0, 0.0};
    nglib::Ng_AddPoint(m, pt);
  }

  // Add the surface elements
  for ( int i = 0; i < ntris; i++ ){
    int tri[3];
    tri[0] = tris[3*i]+1;
    tri[1] = tris[3*i+1]+1;
    tri[2] = tris[3*i+2]+1;
    nglib::Ng_AddSurfaceElement(m, NG_TRIG, tri);
  }

  // Set the mesh parameters
  nglib::Ng_Meshing_Parameters mp;
  mp.maxh = htarget;
  mp.fineness = 1;
  mp.second_order = 0;

  // Generate the volume mesh
  nglib::Ng_GenerateVolumeMesh(m, &mp);

  // Get the total number of points and tets in the volume
  num_points = nglib::Ng_GetNP(m);
  num_tet = nglib::Ng_GetNE(m);

  // Allocate space to store everything
  X = new TMRPoint[ num_points ];
  tet = new int[ 4*num_tet ];

  // Retrieve the points
  for ( int i = 0; i < num_points; i++ ){
    double x[3];
    nglib::Ng_GetPoint(m, i+1, x);
    X[i].x = x[0];
    X[i].y = x[1];
    X[i].z = x[2];
  }

  // Retrieve the tets
  for ( int i = 0; i < num_tet; i++ ){
    nglib::Ng_GetVolumeElement(m, i+1, &tet[4*i]);
    for ( int k = 0; k < 4; k++ ){
      tet[4*i+k] -= 1;
    }
  }

  // Free the memory and exit from netgen
  nglib::Ng_DeleteMesh(m);
  nglib::Ng_Exit();
#endif // TMR_USE_NETGEN
  return 0;
}

/*
  Returns the index number for the (i, j) node location along
  the structured edge.

  This the structured face nodes are ordered by first counting around
  the edges of the face counter-clockwise. After these nodes are
  counted, the interior face nodes are ordered in a cartesian order.
  For nx = 5, and ny = 4, this produces the following face numbering:

  11-- 10-- 9 -- 8 -- 7
  |    |    |    |    |
  12-- 17-- 18-- 19-- 6
  |    |    |    |    |
  13-- 14-- 15-- 16-- 5
  |    |    |    |    |
  0 -- 1 -- 2 -- 3 -- 4
*/
static int get_structured_index( const int nx, const int ny,
                                 const int i, const int j ){
  if (j == 0){
    return i;
  }
  else if (i == nx-1){
    return nx - 1 + j;
  }
  else if (j == ny-1){
    return 2*nx + ny - 3 - i;
  }
  else if (i == 0){
    return 2*nx + 2*ny - 4 - j;
  }

  return 2*nx + 2*ny - 4 + (i-1) + (j-1)*(nx-2);
}

/*
  Order the mesh points uniquely

  This code first copies the mesh locations from the target and source
  surfaces and then the surrounding surfaces that are structured. The
  code takes into account the orientations of the surrounding surfaces
  by computing the absolute orientations of the surface meshes. The
  orientation of the source/target surfaces is already accounted for
  within the connectivity.
*/
int TMRVolumeMesh::setNodeNums( int *num ){
  if (!vars){
    // Initially, set all of the nodes to zero
    vars = new int[ num_points ];
    for ( int k = 0; k < num_points; k++ ){
      vars[k] = -1;
    }

    // Get the target surface mesh and target surface mesh points
    TMRFaceMesh *mesh;
    target->getMesh(&mesh);

    // Number of points in the quadrilateral mesh on the surface
    int num_quad_pts = 0;
    mesh->getMeshPoints(&num_quad_pts, NULL, NULL);

    // Set the nodes on the source face
    const int *face_vars;
    mesh->getNodeNums(&face_vars);

    // Set the target surface variable numbers
    for ( int i = 0; i < num_quad_pts; i++ ){
      vars[i + num_quad_pts*(num_depth_pts-1)] = face_vars[i];
    }

    // Get information from the source surface mesh
    source->getMesh(&mesh);
    mesh->getNodeNums(&face_vars);

    // Set the source face variable numbers
    for ( int i = 0; i < num_quad_pts; i++ ){
      vars[i] = face_vars[i];
    }

    // Now the target and source surfaces of the volume have the correct
    // ordering, but the sides are not ordered correctly. Scan through
    // the structured sides (with the same ordering as defined by the
    // edge loops on the source surface) and determine the node
    // ordering.

    // Keep track of the total number of nodes encountered around the
    // edge of the source face.
    int ioffset = 0;

    // Loop over the nodes that are on surfaces...
    for ( int k = 0; k < num_face_loops; k++ ){
      // Get the source edge loop and source edges
      TMREdgeLoop *source_loop;
      source->getEdgeLoop(k, &source_loop);

      // Get the number of source edges
      int nedges;
      TMREdge **edges;
      source_loop->getEdgeLoop(&nedges, &edges, NULL);

      // Set the original ioffset number
      int ioffset_start = ioffset;

      // If the face loop
      int end = face_loop_ptr[k+1];
      for ( int ii = 0, ptr = face_loop_ptr[k]; ptr < end; ii++, ptr++ ){
        TMRFace *face = face_loops[ptr];

        // Get the face variables
        const int *face_vars;
        face->getMesh(&mesh);
        mesh->getNodeNums(&face_vars);

        // Get the edges associated with the face
        TMREdgeLoop *face_loop;
        face->getEdgeLoop(0, &face_loop);

        int num_face_edges;
        TMREdge **face_edges;
        face_loop->getEdgeLoop(&num_face_edges, &face_edges, NULL);

        // Determine where the source edge appears in the face_edge
        // list. This determines the origin location for the local
        // face ordering.
        int orient = 0;
        for ( ; orient < 4; orient++ ){
          if (edges[ii] == face_edges[orient]){
            break;
          }
        }

        // Determine the number of nodes along the orient edge (nx)
        // and the number of nodes in the transverse direction (ny).
        int nx = 0, ny = 0;
        if (orient == 0 || orient == 2){
          nx = face_loop_edge_count[ptr];
          ny = num_depth_pts;
        }
        else {
          nx = num_depth_pts;
          ny = face_loop_edge_count[ptr];
        }

        // Determine whether the directions are consistent
        int abs_face_dir = 0;
        if (face_loop_dir[ptr] < 0){
          abs_face_dir = -1;
        }
        else {
          abs_face_dir = 1;
        }

        // Scan through the depth points
        for ( int j = 0; j < num_depth_pts; j++ ){
          for ( int i = 0; i < face_loop_edge_count[ptr]; i++ ){
            int local = 0;

            // Determine the local number for the node on the boundary
            // of the face connecting the target and source surfaces
            // based on the orientation of the source face mesh.
            if (source_dir < 0){
              // When this is the last node in the loop, then the
              // ordering must loop back on itself.
              if (i == face_loop_edge_count[ptr]-1 &&
                  ptr == face_loop_ptr[k+1]-1){
                local = j*num_quad_pts + ioffset_start;
              }
              else {
                local = j*num_quad_pts + i + ioffset;
              }
            }
            else {
              // The source direction is reversed, so we order each
              // face in reverse as well
              if (i == 0 && ptr == face_loop_ptr[k+1]-1){
                local = j*num_quad_pts + ioffset_start;
              }
              else {
                local = j*num_quad_pts + ioffset +
                  face_loop_edge_count[ptr] - 1 - i;
              }
            }

            int x = 0, y = 0;
            if (orient == 0){
              if (abs_face_dir < 0){
                x = nx - 1 - i;
                y = j;
              }
              else {
                x = i;
                y = j;
              }
            }
            else if (orient == 1){
              if (abs_face_dir < 0){
                x = nx - 1 - j;
                y = ny - 1 - i;
              }
              else {
                x = nx - 1 - j;
                y = i;
              }
            }
            else if (orient == 2){
              if (abs_face_dir < 0){
                x = i;
                y = ny - 1 - j;
              }
              else {
                x = nx - 1 - i;
                y = ny - 1 - j;
              }
            }
            else if (orient == 3){
              if (abs_face_dir < 0){
                x = j;
                y = i;
              }
              else {
                x = j;
                y = ny - 1 - i;
              }
            }

            int index = get_structured_index(nx, ny, x, y);
            vars[local] = face_vars[index];
          }
        }

        // Update the pointer to the offset
        ioffset += face_loop_edge_count[ptr]-1;
      }
    }

    // Now order the variables as they arrive
    int start = *num;
    for ( int k = 0; k < num_points; k++ ){
      if (vars[k] < 0){
        vars[k] = *num;
        (*num)++;
      }
    }

    // Return the number of points that have been allocated
    return *num - start;
  }

  return 0;
}

/*
  Retrieve the global node numbers for the mesh
*/
int TMRVolumeMesh::getNodeNums( const int **_vars ){
  if (_vars){
    *_vars = vars;
  }
  return num_points;
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
  resetMesh();

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
        const char *attr = volumes[i]->getAttribute();
        if (attr){
          fprintf(stderr,
                  "TMRMesh: Volume meshing failed for object %s\n",
                  attr);
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
      num_hex += mesh->getTetConnectivity(NULL);
    }
  }
  if (num_faces > 0){
    for ( int i = 0; i < num_faces; i++ ){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);
      num_quads += mesh->getQuadConnectivity(NULL);
      num_tris += mesh->getTriConnectivity(NULL);
    }

    // Add the mesh quality to the bins from each mesh
    const int nbins = 20;
    int bins[nbins];
    memset(bins, 0, nbins*sizeof(int));
    for ( int i = 0; i < num_faces; i++ ){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);
      mesh->addMeshQuality(nbins, bins);
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

          // Set the face id == attribute
          char descript[128];
          snprintf(descript, sizeof(descript), "FACE%d",
                   faces[i]->getEntityId());
          if (faces[i]->getAttribute() != NULL){
            strncpy(descript, faces[i]->getAttribute(), sizeof(descript));
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

          // Set the face id == attribute
          char descript[128];
          snprintf(descript, sizeof(descript), "VOLUME%d",
                   volumes[i]->getEntityId());
          if (volumes[i]->getAttribute() != NULL){
            strncpy(descript, volumes[i]->getAttribute(),
                    sizeof(descript));
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
  int vnum = 0;
  TMRVertex **new_verts = new TMRVertex*[ num_nodes ];
  memset(new_verts, 0, num_nodes*sizeof(TMRVertex*));

  // Copy over the vertices in the original geometry
  int num_vertices;
  TMRVertex **vertices;
  geo->getVertices(&num_vertices, &vertices);
  for ( int i = 0; i < num_vertices; i++, vnum++ ){
    new_verts[vnum] = vertices[i];
  }

  // Create the vertices on the edges
  int num_edges;
  TMREdge **edges;
  geo->getEdges(&num_edges, &edges);
  for ( int i = 0; i < num_edges; i++ ){
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);

    // Get the parametric points associated with the mesh
    int npts;
    const double *tpts;
    mesh->getMeshPoints(&npts, &tpts, NULL);
    for ( int j = 1; j < npts-1; j++, vnum++ ){
      new_verts[vnum] = new TMRVertexFromEdge(edges[i], tpts[j]);
    }
  }

  // Create the vertices from the faces
  int num_faces;
  TMRFace **faces;
  geo->getFaces(&num_faces, &faces);
  for ( int i = 0; i < num_faces; i++ ){
    TMRFaceMesh *mesh = NULL;
    faces[i]->getMesh(&mesh);

    // Get the mesh points
    int npts;
    const double *pts;
    mesh->getMeshPoints(&npts, &pts, NULL);

    // Get the point-offset for this surface
    int offset = mesh->getNumFixedPoints();
    for ( int j = offset; j < npts; j++, vnum++ ){
      new_verts[vnum] = new TMRVertexFromFace(faces[i],
                                              pts[2*j], pts[2*j+1]);
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
    computeHexEdgesAndFaces(num_nodes, num_hex, hex,
                            &num_hex_edges, &hex_edges, &hex_edge_nums,
                            &num_hex_faces, &hex_faces, &hex_face_nums);
  }
  else {
    // Compute the quadrilateral edges
    computeQuadEdges(num_nodes, num_quads, quads,
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
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);

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
      else {
        fprintf(stderr,
                "TMRMesh error: Could not find edge (%d, %d) to split\n",
                edge[1], edge[2]);
      }
    }
  }

  // Create the edges on the surface using Pcurve/CurveFromSurface
  for ( int i = 0; i < num_faces; i++ ){
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
            // always in the 'positive' orientation.
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
      // face number within the sorted_faces array to obtain the required
      // face number
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
        fdir[j] = 1;
        f[j] = new_faces[face_num];
      }

      // Allocate the new transfinite interpolation face
      new_volumes[i] = new TMRTFIVolume(f, fdir, e, edir, v);
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
