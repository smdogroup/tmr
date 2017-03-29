#include "TMRMesh.h"
#include "TMRNativeTopology.h"
#include "TMRTriangularize.h"
#include "TMRMeshSmoothing.h"
#include "TMRPerfectMatchInterface.h"
#include "TMRBspline.h"
#include <math.h>
#include <stdio.h>



/*
  The nodes in this hexahedral element are ordered as follows:
  
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
   {0,4}, {1,5}, {2,6}, {3,7}}; // z-aligned edges

/*
  The nodes corresponding to each face in a hexahedral
  element
*/
const int hex_face_nodes[6][4] =
  {{0,3,4,7}, {1,2,5,6},  // x-faces
   {0,1,4,5}, {3,2,7,6},  // y-faces
   {0,1,3,2}, {4,5,7,8}}; // z-faces

/*
  Possible orientations for two connecting faces
*/
const int face_orient[8][4] = 
  {{0,1,2,3}, {2,0,3,1},
   {3,2,1,0}, {1,3,0,2},
   {0,2,1,3}, {2,3,0,1},
   {3,1,2,0}, {1,0,3,2}};

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
    else if (j >= 2){
      node[k] = tmp[i];
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
                             int **_node_to_tris=NULL ){
  // Compute the edges in the triangular mesh
  int *ptr;
  int *node_to_tris;
  computeNodeToElems(nnodes, ntris, 3, tris, &ptr, &node_to_tris);
  
  // Now compute the neighbors for each triangle
  int *tri_edge_nums = new int[ 3*ntris ];
  for ( int i = 0; i < 3*ntris; i++ ){
    tri_edge_nums[i] = -1;
  }

  // Quck reference from the edge index to the local node numbering
  const int enodes[3][2] = {{1, 2}, {2, 0}, {0, 1}};

  // Allocate the array for the triangle neighbors
  int *tri_neighbors = new int[ 3*ntris ];
  
  int count = 0;
  int ne = 0;
  for ( int i = 0; i < ntris; i++ ){
    // Search through each edge of the each triangle
    for ( int j = 0; j < 3; j++ ){
      if (tri_edge_nums[3*i+j] < 0){
        tri_edge_nums[3*i+j] = ne;

        // Triangle edges that have no neighbors are labeled with a -1
        tri_neighbors[3*i+j] = -1;

        // Search for the neighboring quad that shares this edge
        int kp = ptr[tris[3*i + enodes[j][0]]];
        int kpend = ptr[tris[3*i + enodes[j][0]]+1];
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
            // Check if the adjacent edge matches in either direction
            if ((tris[3*i+enodes[j][0]] == tris[3*n+enodes[e][0]] &&
                 tris[3*i+enodes[j][1]] == tris[3*n+enodes[e][1]]) ||
                (tris[3*i+enodes[j][0]] == tris[3*n+enodes[e][1]] &&
                 tris[3*i+enodes[j][1]] == tris[3*n+enodes[e][0]])){
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
      tri_edges[2*n] = tris[3*i+enodes[j][0]];
      tri_edges[2*n+1] = tris[3*i+enodes[j][1]];

      // Set the dual edge numbers - connecting triangles to other
      // triangles. Note that some elements of this array will be -1.
      dual_edges[2*n] = i;
      dual_edges[2*n+1] = tri_neighbors[3*i+j];
    }
  } 

  delete [] tri_edge_nums;

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

  // Quck reference from the quad index to the edge
  const int enodes[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

  int count = 0;
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

        // Search for the neighboring quad that shares this edge
        int kp = ptr[quads[4*i + enodes[j][0]]];
        int kpend = ptr[quads[4*i + enodes[j][0]]+1];
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
            // Check if the adjacent edge matches in either direction
            if ((quads[4*i+enodes[j][0]] == quads[4*n+enodes[e][0]] &&
                 quads[4*i+enodes[j][1]] == quads[4*n+enodes[e][1]]) ||
                (quads[4*i+enodes[j][0]] == quads[4*n+enodes[e][1]] &&
                 quads[4*i+enodes[j][1]] == quads[4*n+enodes[e][0]])){
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
      quad_edges[2*n] = quads[4*i+enodes[j][0]];
      quad_edges[2*n+1] = quads[4*i+enodes[j][1]];

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
static void computeHexEdgesAndFaces( int nnodes, int nhexes,
                                     const int hexes[],
                                     int *_num_hex_edges,
                                     int **_hex_edges,
                                     int **_hex_edge_nums,
                                     int *_num_hex_faces,
                                     int **_hex_faces,
                                     int **_hex_face_nums ){
  // Compute the connectivity from nodes to hexes
  int *ptr;
  int *node_to_hexes;
  computeNodeToElems(nnodes, nhexes, 8, hexes, &ptr, &node_to_hexes);

  // Comput the neighbors for each hex
  int *hex_edge_nums = new int[ 12*nhexes ];
  int *hex_face_nums = new int[ 6*nhexes ];
  for ( int i = 0; i < 12*nhexes; i++ ){
    hex_edge_nums[i] = -1;
  }
  for ( int i = 0; i < 6*nhexes; i++ ){
    hex_face_nums[i] = -1;
  }

  int edge_num = 0, face_num = 0;
  for ( int i = 0; i < nhexes; i++ ){
    // Search through each hexahedral element for the edge
    for ( int j = 0; j < 12; j++ ){
      if (hex_edge_nums[12*i+j] < 0){
        hex_edge_nums[12*i+j] = edge_num;

        // Find the edge nodes for the new edge
        int e[2];
        e[0] = hexes[8*i + hex_edge_nodes[j][0]];
        e[1] = hexes[8*i + hex_edge_nodes[j][1]];

        // Find the node that shares the edge
        int kp = ptr[e[0]];
        int kpend =ptr[e[0]+1];

        for ( ; kp < kpend; kp++ ){
          // Find the potential hex neighbor
          int hex = node_to_hexes[kp];
          if (i == hex){
            continue;
          }

          // Search over all the edges on this hex, and see if any of
          // the edges match
          for ( int k = 0; k < 12; k++ ){
            int e2[2];
            e2[0] = hexes[8*hex + hex_edge_nodes[k][0]];
            e2[1] = hexes[8*hex + hex_edge_nodes[k][1]];

            // Check if the adjacent edge matches in either direction
            if ((e[0] == e2[0] && e[1] == e2[1]) ||
                (e[0] == e2[1] && e[1] == e2[0])){
              hex_edge_nums[12*hex + k] = edge_num;
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
        f[0] = hexes[8*i + hex_face_nodes[j][0]];
        f[1] = hexes[8*i + hex_face_nodes[j][1]];
        f[2] = hexes[8*i + hex_face_nodes[j][2]];
        f[3] = hexes[8*i + hex_face_nodes[j][3]];
        
        // Find a node that shares the face
        int kp = ptr[f[0]];
        int kpend =ptr[f[0]+1];

        for ( ; kp < kpend; kp++ ){
          // Find the potential hex neighbor
          int hex = node_to_hexes[kp];
          if (i == hex){
            continue;
          }

          // Search over all the faces on this hex, and see if any of
          // the edges match
          for ( int k = 0; k < 6; k++ ){
            int f2[4];
            f2[0] = hexes[8*hex + hex_face_nodes[k][0]];
            f2[1] = hexes[8*hex + hex_face_nodes[k][1]];
            f2[2] = hexes[8*hex + hex_face_nodes[k][2]];
            f2[3] = hexes[8*hex + hex_face_nodes[k][3]];

            // Check if the adjacent edge matches in either direction
            for ( int ort = 0; ort < 8; ort++ ){
              if (f[0] == f2[face_orient[ort][0]] &&
                  f[1] == f2[face_orient[ort][1]] &&
                  f[2] == f2[face_orient[ort][2]] &&
                  f[3] == f2[face_orient[ort][3]]){
                hex_face_nums[6*hex + k] = face_num;
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

  if (_num_hex_edges){
    *_num_hex_edges = edge_num;
  }
  if (_num_hex_faces){
    *_num_hex_faces = face_num;
  }

  // Create the hex edges/faces
  int *hex_edges = new int[ 2*edge_num ];
  int *hex_faces = new int[ 4*face_num ];
  for ( int i = 0; i < nhexes; i++ ){
    for ( int j = 0; j < 12; j++ ){
      int e = hex_edge_nums[12*i+j];
      hex_edges[2*e] = hexes[8*i + hex_edge_nodes[j][0]];
      hex_edges[2*e+1] = hexes[8*i + hex_edge_nodes[j][1]];
    }

    for ( int j = 0; j < 6; j++ ){
      int f = hex_face_nums[6*i+j];
      hex_faces[4*f] = hexes[8*i + hex_face_nodes[j][0]];
      hex_faces[4*f+1] = hexes[8*i + hex_face_nodes[j][1]];
      hex_faces[4*f+2] = hexes[8*i + hex_face_nodes[j][2]];
      hex_faces[4*f+3] = hexes[8*i + hex_face_nodes[j][3]];
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
  Mesh a curve
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
void TMREdgeMesh::mesh( TMRMeshOptions options, double htarget ){
  int mpi_rank, mpi_size;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Get the master edge
  TMREdge *master;
  edge->getMaster(&master);

  // Figure out if there is a master edge and whether or not it has
  // been meshed.
  npts = -1;
  if (master && master != edge){
    TMREdgeMesh *mesh;
    master->getMesh(&mesh);
    if (!mesh){
      mesh = new TMREdgeMesh(comm, master);
      mesh->mesh(options, htarget);
      master->setMesh(mesh);
    }

    // Retrieve the number of points along the master edge
    mesh->getMeshPoints(&npts, NULL, NULL);
  }

  if (mpi_rank == 0){
    // Get the limits of integration that will be used
    double tmin, tmax;
    edge->getRange(&tmin, &tmax);

    if (!edge->isDegenerate()){
      // Set the integration error tolerance
      double integration_eps = 1e-8;

      // Integrate along the curve to obtain the distance function such
      // that dist(tvals[i]) = int_{tmin}^{tvals[i]} ||d{C(t)}dt||_{2} dt
      int nvals;
      double *dist, *tvals;
      edge->integrate(tmin, tmax, integration_eps, 
                      &tvals, &dist, &nvals);
      
      // Only compute the number of points if there is no master edge
      if (npts < 0){
        // Compute the number of points along this curve
        npts = 1 + (int)(dist[nvals-1]/htarget);
        if (npts < 2){ npts = 2; }

        // If we have an even number of points, increment by one to ensure
        // that we have an even number of segments along the boundary
        if (npts % 2 != 1){ npts++; }
      }

      // The average distance between points
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

    return npts-1;
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
void TMREdgeMesh::getMeshPoints( int *_npts, const double **_pts, 
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
  mesh_type = NO_MESH;

  num_fixed_pts = 0;
  num_points = 0;
  num_quads = 0;

  // NULL things that will be used later
  pts = NULL;
  X = NULL;
  quads = NULL;
  vars = NULL;
}

/*
  Free the data associated with the surface mesh
*/
TMRFaceMesh::~TMRFaceMesh(){
  face->decref();
  if (pts){ delete [] pts; }
  if (X){ delete [] X; }
  if (quads){ delete [] quads; }
  if (vars){ delete [] vars; }
}

/*
  Create the surface mesh
*/
void TMRFaceMesh::mesh( TMRMeshOptions options,
                        double htarget, 
                        TMRFaceMeshType _mesh_type ){
  int mpi_rank, mpi_size;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  if (_mesh_type == NO_MESH){
    _mesh_type = STRUCTURED;
  }
  
  // Get the master face and its orientation relative to this
  // face. Note that the master face may be NULL in which case the
  // master orientation is meaningless.
  int master_dir;
  TMRFace *master;
  face->getMaster(&master_dir, &master);

  if (master){
    // If the face mesh for the master does not yet exist, create it...    
    TMRFaceMesh *face_mesh;
    master->getMesh(&face_mesh);
    if (!face_mesh){
      face_mesh = new TMRFaceMesh(comm, master);
      face_mesh->mesh(options, htarget);
      master->setMesh(face_mesh);
    }
  }

  // First check if the conditions for a structured mesh are satisfied
  if (_mesh_type == STRUCTURED){
    int nloops = face->getNumEdgeLoops();
    if (nloops != 1){
      _mesh_type = UNSTRUCTURED;
    }

    // Get the first edge loop and the edges in the loop
    TMREdgeLoop *loop;
    face->getEdgeLoop(0, &loop);
    int nedges;
    TMREdge **edges;
    loop->getEdgeLoop(&nedges, &edges, NULL);
    if (nedges != 4){
      _mesh_type = UNSTRUCTURED;
    }
    for ( int k = 0; k < nedges; k++ ){
      if (edges[k]->isDegenerate()){
        _mesh_type = UNSTRUCTURED;
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
        _mesh_type = UNSTRUCTURED;
      }

      edges[1]->getMesh(&mesh);
      mesh->getMeshPoints(&ny1, NULL, NULL);
      edges[3]->getMesh(&mesh);
      mesh->getMeshPoints(&ny2, NULL, NULL);
      if (ny1 != ny2){
        _mesh_type = UNSTRUCTURED;
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
    
    // Get all of the meshes
    for ( int k = 0; k < nloops; k++ ){
      TMREdgeLoop *loop;
      face->getEdgeLoop(k, &loop);
      int nedges;
      TMREdge **edges;
      loop->getEdgeLoop(&nedges, &edges, NULL);
      
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
    int init_loop_seg = 0; // What segment value did this loop start on?
    int hole_pt = total_num_pts; // What hole are we on?

    // Set up the degenerate edges
    int *degen = new int[ 2*num_degen ];
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

      // Check if the area constraint
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

    if (master){
      TMRFaceMesh *face_mesh;
      master->getMesh(&face_mesh);

      // Compute the total number of points
      mesh_type = face_mesh->mesh_type;
      num_points = face_mesh->num_points;
      num_fixed_pts = face_mesh->num_fixed_pts;
      num_quads = face_mesh->num_quads;

      quads = new int[ 4*num_quads ];
      memcpy(quads, face_mesh->quads, 4*num_quads*sizeof(int));

      // Allocate the array for the parametric locations
      pts = new double[ 2*num_points ];

      // Copy the points from around the boundaries
      for ( int i = 0; i < num_fixed_pts; i++ ){
        pts[2*i] = params[2*i];
        pts[2*i+1] = params[2*i+1];
      }

      // If the sense of the master face is opposite to the
      // orientation of this face, then reverse the quads
      /*
      if (master_dir < 0){
        for ( int i = 0; i < num_quads; i++ ){
          int tmp = quads[4*i+1];
          quads[4*i+1] = quads[4*i+3];
          quads[4*i+3] = tmp;
        }
      }
      */

      // Copy the interior point locations
      memcpy(&pts[2*num_fixed_pts], &(face_mesh->pts[2*num_fixed_pts]),
             2*(num_points - num_fixed_pts)*sizeof(double));

      // Evaluate the points
      X = new TMRPoint[ num_points ];
      for ( int i = 0; i < num_points; i++ ){
        face->evalPoint(pts[2*i], pts[2*i+1], &X[i]);
      }

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
    else if (mesh_type == STRUCTURED){
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
      tri->frontal(htarget);

      // Free the degenerate triangles and reorder the mesh
      if (num_degen > 0){
        tri->removeDegenerateEdges(num_degen, degen);
      }

      if (options.write_pre_smooth_triangle){
        char filename[256];
        sprintf(filename, "pre_smooth_triangle%d.vtk", 
                face->getEntityId());
        tri->writeToVTK(filename);
      }

      // Extract the triangularization
      int ntris, *tris;
      tri->getMesh(&num_points, &ntris, &tris, &pts, &X);
      tri->decref();

      if (ntris > 0){
        // Compute the triangle edges and neighbors in the dual mesh
        int num_tri_edges;
        int *tri_edges, *tri_neighbors, *dual_edges;
        int *node_to_tri_ptr, *node_to_tris;
        computeTriEdges(num_points, ntris, tris, 
                        &num_tri_edges, &tri_edges,
                        &tri_neighbors, &dual_edges,
                        &node_to_tri_ptr, &node_to_tris);
      
        // Smooth the resulting triangular mesh
        if (options.tri_smoothing_type == TMRMeshOptions::LAPLACIAN){
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
          writeTrisToVTK(filename, ntris, tris);
        }

        // Recombine the mesh into a quadrilateral mesh
        recombine(ntris, tris, tri_neighbors, 
                  node_to_tri_ptr, node_to_tris,
                  num_tri_edges, dual_edges, &num_quads, &quads, options);

        // Free the triangular mesh data
        delete [] tris;
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

    // Free the parameter/segment information
    delete [] params;
    delete [] segments;

    // Flip the orientation of the face so that the normal directions
    // are consistent. If a master face was used to define this face,
    // then the normal direction is already consistent and no
    // quad-flip is performed.
    if (!master && face->getNormalDirection() < 0){
      for ( int i = 0; i < num_quads; i++ ){
        int tmp = quads[4*i+1];
        quads[4*i+1] = quads[4*i+3];
        quads[4*i+3] = tmp;
      }
    }
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
  Get the local connectivity
*/
int TMRFaceMesh::getLocalConnectivity( const int **_quads ){
  if (_quads){ *_quads = quads; }
  return num_quads;
}

/*
  Get the quadrilateral elemnet obtained by combining the triangles
  t1 and t2 together

        . --- .
      / |   /
    /   | /
  . --- . 
*/
int TMRFaceMesh::getRecombinedQuad( const int tris[],
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
    quad[j] = tris[3*t1+i];
    if (i == e1){
      j++;
    }
  }

  // Set the node contributed by the second triangle
  e2 += 2;
  if (e2 >= 3){ e2 -= 3; }
  quad[e1+1] = tris[3*t2+e2];

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
double TMRFaceMesh::computeRecombinedQuality( const int tris[], 
                                              const int tri_neighbors[],
                                              int t1, int t2,
                                              const TMRPoint *p ){
  // Find the combined quadrilateral from the two given triangles
  int quad[4];
  int fail = getRecombinedQuad(tris, tri_neighbors, t1, t2, quad);
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
   
    // Compute the internal angle between the 
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
void TMRFaceMesh::recombine( int ntris, const int tris[],
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
  const double frac = eps/(1.0 + eps);

  int edge_num = 0;
  for ( int i = 0; i < num_edges; i++ ){
    int t1 = dual_edges[2*i];
    int t2 = dual_edges[2*i+1];

    if (t1 >= 0 && t2 >= 0){
      // Compute the weight for this recombination
      double quality = 
        computeRecombinedQuality(tris, tri_neighbors,
                                 t1, t2, X);

      // double weight = frac*(1.0 - quality)*(1.0 + 1.0/(quality + eps));
      double weight = 1.0 - quality;
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
      // Quck reference from the edge index to the local node numbering
      const int enodes[3][2] = {{1, 2}, {2, 0}, {0, 1}};
      const int node_edges[3][2] = {{1, 2}, {0, 2}, {0, 1}};

      // Loop over the edges of the triangle
      for ( int j = 0; j < 3; j++ ){
        // We have a boundary triangle
        if (tri_neighbors[3*i+j] < 0){
          // The leading node that we are searching for
          int ij = enodes[j][1]; 
          int node = tris[3*i + ij];

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
            if (tris[3*k+1] == node){ kj = 1; }
            else if (tris[3*k+2] == node){ kj = 2; }

            // Track from the node to the possible edges
            if (tri_neighbors[3*k + node_edges[kj][0]] < 0 ||
                tri_neighbors[3*k + node_edges[kj][1]] < 0){
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
    writeDualToVTK(filename, 3, ntris, tris, num_dual_edges,
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

      int fail = getRecombinedQuad(tris, tri_neighbors, t1, t2,
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
          if (tris[3*t1 + j1] == tris[3*t2 + j2]){
            flag = 1;
            break;
          }
        }
        if (flag){
          break; 
        }
      }

      int boundary_pt = tris[3*t1 + j1];
      
      // Go through all previous quadrilaterals and adjust the
      // ordering to reflect the duplicated node location
      for ( int k = 0; k < 4*num_quads_from_tris; k++ ){
        if (new_quads[k] == boundary_pt){
          new_quads[k] = num_new_points;
        }
      }

      // Add the first triangle t1 - this triangle is gauranteed to
      // come first when circling the boundary in the CCW direction.
      new_quads[4*num_new_quads] = tris[3*t1+((j1+1) % 3)];
      new_quads[4*num_new_quads+1] = tris[3*t1+((j1+2) % 3)];
      new_quads[4*num_new_quads+2] = boundary_pt;
      new_quads[4*num_new_quads+3] = num_new_points;
      num_new_quads++;
      
      // Add the connectivity from the second triangle t2. This
      // triangle will always come second when circling the boundary.
      new_quads[4*num_new_quads] = tris[3*t2+((j2+1) % 3)];
      new_quads[4*num_new_quads+1] = tris[3*t2+((j2+2) % 3)];
      new_quads[4*num_new_quads+2] = num_new_points;
      new_quads[4*num_new_quads+3] = boundary_pt;
      num_new_quads++;

      // Compute the new parameter location by taking the average of
      // the centroid locations for each triangle
      int count = 1;
      new_pts[2*num_new_points] = pts[2*boundary_pt];
      new_pts[2*num_new_points+1] = pts[2*boundary_pt+1];

      int kpend = node_to_tri_ptr[boundary_pt+1];
      for ( int kp = node_to_tri_ptr[boundary_pt]; kp < kpend; kp++ ){
        int k = node_to_tris[kp];
        if (k != t1 && k != t2){
          new_pts[2*num_new_points] += (pts[2*tris[3*k]] +
                                        pts[2*tris[3*k+1]] +
                                        pts[2*tris[3*k+2]])/3.0;
          new_pts[2*num_new_points+1] += (pts[2*tris[3*k]+1] +
                                          pts[2*tris[3*k+1]+1] +
                                          pts[2*tris[3*k+2]+1])/3.0;
          count++;
        }
      }
      if (count > 1){
        new_pts[2*num_new_points] = new_pts[2*num_new_points]/count;
        new_pts[2*num_new_points+1] = new_pts[2*num_new_points+1]/count;
      }
      face->evalPoint(new_pts[2*num_new_points],
                      new_pts[2*num_new_points+1], &new_X[num_new_points]);
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
  // Compute the pointer -> quad information
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
      for ( int k = 0; k < 4; k++ ){
        if (new_pt_nums[quad1[k]] < 0 || 
            new_pt_nums[quad2[k]] < 0){
          break;
        }
      }

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
  int j = 0;
  for ( ; j < num_fixed_pts; j++ ){
    new_pt_nums[j] = j;
  }
  for ( int i = num_fixed_pts; i < num_points; i++ ){
    if (new_pt_nums[i] >= 0){
      if (i != j){
        // Copy over the points
        pts[2*j] = pts[2*i];
        pts[2*j+1] = pts[2*i+1];
        X[j] = X[i];
      }
      new_pt_nums[i] = j;
      j++;
    }
  }

  // Set the new number of points
  num_points = j;

  // Set the new connectivity, overwriting the old connectivity 
  j = 0;
  for ( int i = 0; i < num_quads; i++ ){
    if (quads[4*i] >= 0){
      quads[4*j] = new_pt_nums[quads[4*i]];
      quads[4*j+1] = new_pt_nums[quads[4*i+1]];
      quads[4*j+2] = new_pt_nums[quads[4*i+2]];
      quads[4*j+3] = new_pt_nums[quads[4*i+3]];
      j++;
    }
  }

  // Set the new number of quadrilaterals
  num_quads = j;

  // Free the new point numbers
  delete [] new_pt_nums;
}

/*
  Add the quad quality
*/
void TMRFaceMesh::addQuadQuality( int nbins, int bins[] ){
  for ( int i = 0; i < num_quads; i++ ){
    double quality = computeQuadQuality(&quads[4*i], X);

    int k = 0;
    for ( ; k < nbins; k++ ){
      if (quality < 1.0*(k+1)/nbins){
        break;
      }
    }
    bins[k]++;
  }
}

/*
  Print the quadrilateral quality
*/
void TMRFaceMesh::printQuadQuality(){
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins*sizeof(int));
  addQuadQuality(nbins, bins);

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
                                   const int tris[] ){
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins*sizeof(int));
  for ( int i = 0; i < ntris; i++ ){
    double quality = computeTriQuality(&tris[3*i], X);

    int k = 0;
    for ( ; k < nbins; k++ ){
      if (quality < 1.0*(k+1)/nbins){
        break;
      }
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
                                  int ntris, const int tris[] ){
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
      fprintf(fp, "3 %d %d %d\n", tris[3*k], tris[3*k+1], tris[3*k+2]);
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
      fprintf(fp, "%e\n", computeTriQuality(&tris[3*k], X));
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
                1.0/3.0*(p[elems[3*k]].x + p[elems[3*k+1]].x + p[elems[3*k+2]].x),
                1.0/3.0*(p[elems[3*k]].y + p[elems[3*k+1]].y + p[elems[3*k+2]].y),
                1.0/3.0*(p[elems[3*k]].z + p[elems[3*k+1]].z + p[elems[3*k+2]].z));
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

  // Set the number of loops on the bottom face
  num_face_loops = 0;
  face_loops = NULL;
  face_loop_ptr = NULL;
  face_loop_edge_count = NULL;

  // Set the number of points through-thickness
  num_depth_pts = -1;

  // Set the points in the mesh
  num_points = 0;
  num_hexes = 0;

  // Set all the pointer to null
  X = NULL;
  hexes = NULL;
  vars = NULL;

  // Set the additional connectivity information that is required for
  // the volume mesh.
  bottom = top = NULL;
}

TMRVolumeMesh::~TMRVolumeMesh(){
  volume->decref();
  if (X){ delete [] X; }
  if (hexes){ delete [] hexes; }
  if (vars){ delete [] vars; }
  if (face_loops){
    for ( int k = 0; k < face_loop_ptr[num_face_loops]; k++ ){
      face_loops[k]->decref();
    }
    delete [] face_loop_ptr;
    delete [] face_loops;
    delete [] face_loop_edge_count;
  }
}

/*
  Mesh the volume
*/
int TMRVolumeMesh::mesh( TMRMeshOptions options ){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Keep track of whether the mesh has failed at any time.
  // Try and print out a helpful message.
  int mesh_fail = 0;

  // Get the faces associated with the volume
  int num_faces;
  TMRFace **faces;
  const int *face_dir;
  volume->getFaces(&num_faces, &faces, &face_dir);

  // Set integers for each face to determine whether it is
  // a master or a slave face or it is connected to one of them
  // or whether it is structured
  int *f = new int[ num_faces ];
  memset(f, 0, num_faces*sizeof(int));

  // Determine if the faces are all referenced
  top = NULL;
  top_dir = 1;
  bottom = NULL;
  bottom_dir = 1;

  for ( int i = 0; i < num_faces; i++ ){
    TMRFace *master;
    faces[i]->getMaster(NULL, &master);
    if (master){
      top = master;
      bottom = faces[i];
      bottom_dir = face_dir[i];

      // Set the remaining flags
      f[i] = 1;
      for ( int j = 0; j < num_faces; j++ ){
        if (faces[j] == master){
          f[j] = 1;
          top_dir = face_dir[j];
        }
      }

      // Break here. The meshing algorithm only works with one
      // master/slave face pair per volume.
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
      else if (mesh->getMeshType() != TMRFaceMesh::STRUCTURED){
        fprintf(stderr,
                "TMRVolumeMesh error: \
Through-thickness meshes must be structured\n");
        fprintf(stderr,
                "Try setting master-slave relations on edges and surfaces\n");
        mesh_fail = 1;
      }
    }
  }

  if (mesh_fail){
    return mesh_fail;
  }

  // Free the f pointer
  delete [] f;

  // Each face that is not either the top or the bottom, must have
  // only four edges and must be structured. Two of the parallel edges
  // associated with these faces should touch the top and bottom face,
  // while remaining edges must be parallel and have the same number
  // of nodes. Furthermore, the one parallel edge will be shared by
  // the next face. We loop over the bottom edge loops and find the
  // connecting faces.

  // Get the number of edge loops
  num_face_loops = bottom->getNumEdgeLoops();

  // Count up the total number of faces that 
  face_loop_ptr = new int[ num_face_loops+1 ];

  int count = 0;
  for ( int k = 0; k < num_face_loops; k++ ){
    TMREdgeLoop *bottom_loop;
    bottom->getEdgeLoop(k, &bottom_loop);

    // Get the number of edges for this loop
    int nedges;
    bottom_loop->getEdgeLoop(&nedges, NULL, NULL);
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

    // The top and bottom edge loops
    TMREdgeLoop *bottom_loop;
    bottom->getEdgeLoop(k, &bottom_loop);

    // Get the edges associated with the bottom loop
    int nedges;
    TMREdge **edges;
    bottom_loop->getEdgeLoop(&nedges, &edges, NULL);

    for ( int j = 0; j < nedges; j++ ){
      // Search for the other face object that shares the edge
      // object edges[j] with the bottom face object
      TMRFace *face = NULL;
      int fdir = 0;
      for ( int i = 0; i < num_faces; i++ ){
        if (!(faces[i] == top || faces[i] == bottom)){
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

      // This face is not the top or bottom face and therefore must
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

  // Get the information for the top surface
  TMRFaceMesh *mesh;
  top->getMesh(&mesh);

  // Number of points in the quadrilateral mesh on the surface
  int num_quad_pts = 0;
  TMRPoint *Xtop;
  mesh->getMeshPoints(&num_quad_pts, NULL, &Xtop);

  // Get information for the bottom surface
  bottom->getMesh(&mesh);
    
  // Points on the top surface
  TMRPoint *Xbot;
  mesh->getMeshPoints(NULL, NULL, &Xbot);

  // Get the local connectivity on the bottom surface
  const int *quads;
  int num_quads = mesh->getLocalConnectivity(&quads);

  // The number of hexahedral elements in the mesh
  num_hexes = (num_depth_pts-1)*num_quads;
  num_points = num_depth_pts*num_quad_pts;
  hexes = new int[ 8*num_hexes ];
  X = new TMRPoint[ num_points ];

  // Flip the  to the hexahedral elements
  const int flip[] = {0, 3, 2, 1};

  int *hex = hexes;
  for ( int j = 0; j < num_depth_pts-1; j++ ){
    for ( int i = 0; i < num_quads; i++ ){
      // Set the quadrilateral points in the base layer
      if (bottom_dir > 0){
        for ( int k = 0; k < 4; k++ ){
          hex[k] = j*num_quad_pts + quads[4*i+flip[k]];
        }
        for ( int k = 0; k < 4; k++ ){
          hex[4+k] = (j+1)*num_quad_pts + quads[4*i+flip[k]];
        }
      }
      else {
        for ( int k = 0; k < 4; k++ ){
          hex[k] = j*num_quad_pts + quads[4*i+k];
        }
        for ( int k = 0; k < 4; k++ ){
          hex[4+k] = (j+1)*num_quad_pts + quads[4*i+k];
        }
      }
      hex += 8;
    }
  }

  // Set the new coordinates within the hexahedral mesh
  TMRPoint *x = X;
  for ( int j = 0; j < num_depth_pts; j++ ){
    double u = 1.0*j/(num_depth_pts-1);
    for ( int i = 0; i < num_quad_pts; i++ ){
      x[0].x = (1.0 - u)*Xbot[i].x + u*Xtop[i].x;
      x[0].y = (1.0 - u)*Xbot[i].y + u*Xtop[i].y;
      x[0].z = (1.0 - u)*Xbot[i].z + u*Xtop[i].z;
      x += 1;
    }
  }

  return mesh_fail;
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
int TMRVolumeMesh::getLocalConnectivity( const int **_hexes ){
  if (_hexes){ *_hexes = hexes; }
  return num_hexes;
}

/*
  Write the volume mesh to a VTK file
*/
void TMRVolumeMesh::writeToVTK( const char *filename ){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  if (mpi_rank == 0 && hexes){
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
      fprintf(fp, "\nCELLS %d %d\n", num_hexes, 9*num_hexes);
      for ( int k = 0; k < num_hexes; k++ ){
        fprintf(fp, "8 %d %d %d %d %d %d %d %d\n", 
                hexes[8*k], hexes[8*k+1], hexes[8*k+2], hexes[8*k+3],
                hexes[8*k+4], hexes[8*k+5], hexes[8*k+6], hexes[8*k+7]);
      }

      // All hexes
      fprintf(fp, "\nCELL_TYPES %d\n", num_hexes);
      for ( int k = 0; k < num_hexes; k++ ){
        fprintf(fp, "%d\n", 12);
      }
      fclose(fp);
    }
  }
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

  This code first copies the mesh locations from the top and bottom
  surfaces and then the surrounding surfaces that are structured. The
  code takes into account the orientations of the surrounding surfaces
  by computing the absolute orientations of the surface meshes. The
  orientation of the bottom/top surfaces is already accounted for
  within the connectivity.
*/
int TMRVolumeMesh::setNodeNums( int *num ){
  if (!vars){
    // Initially, set all of the nodes to zero
    vars = new int[ num_points ];
    for ( int k = 0; k < num_points; k++ ){
      vars[k] = -1;
    }

    // Get the top surface mesh and top surface mesh points
    TMRFaceMesh *mesh;
    top->getMesh(&mesh);

    // Number of points in the quadrilateral mesh on the surface
    int num_quad_pts = 0;
    mesh->getMeshPoints(&num_quad_pts, NULL, NULL);

    // Set the nodes on the bottom face
    const int *face_vars;
    mesh->getNodeNums(&face_vars);

    // Set the top surface variable numbers
    for ( int i = 0; i < num_quad_pts; i++ ){
      vars[i + num_quad_pts*(num_depth_pts-1)] = face_vars[i];
    }

    // Get information from the bottom surface mesh
    bottom->getMesh(&mesh);
    mesh->getNodeNums(&face_vars);

    // Set the bottom face variable numbers
    for ( int i = 0; i < num_quad_pts; i++ ){
      vars[i] = face_vars[i];
    }

    // Now the top and bottom surfaces of the volume have the correct
    // ordering, but the sides are not ordered correctly. Scan through
    // the structured sides (with the same ordering as defined by the
    // edge loops on the bottom surface) and determine 

    // The internal ordering is relative to the bottom face with its
    // normal direction pointing in to the volume. If the bottom face
    // is oriented with the normal facing outwards, then we have to
    // flip the internal ordering of the boundary nodes to make them
    // consistent with the faces.
    int abs_bottom_dir = 0;
    if (bottom_dir*bottom->getNormalDirection() < 0){
      abs_bottom_dir = -1;      
    }
    else {
      abs_bottom_dir = 1;
    }

    // Keep track of the total number of nodes encountered around the
    // edge of the bottom face.
    int ioffset = 0;

    // Loop over the nodes that are on surfaces...
    for ( int k = 0; k < num_face_loops; k++ ){
      // Get the bottom edge loop and bottom edges
      TMREdgeLoop *bottom_loop;
      bottom->getEdgeLoop(k, &bottom_loop);

      // Get the number of bottom edges
      int nedges;
      TMREdge **edges;
      bottom_loop->getEdgeLoop(&nedges, &edges, NULL);

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

        // Determine where the bottom edge appears in the face_edge
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
        if (face_loop_dir[ptr]*face->getNormalDirection() < 0){
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
            // of the face connecting the top and bottom surfaces
            // based on the orientation of the bottom face mesh.
            if (abs_bottom_dir < 0){
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
              // The bottom direction is reversed, so we order each
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
}

/*
  Mesh the given geometry and retrieve either a regular mesh
*/
TMRMesh::TMRMesh( MPI_Comm _comm, TMRModel *_geo ){
  comm = _comm;
  geo = _geo;
  geo->incref();

  // Set the mesh properties
  num_nodes = 0;
  num_quads = 0;
  quads = NULL;
  X = NULL;
}

TMRMesh::~TMRMesh(){
  geo->decref();
  if (quads){ delete [] quads; }
  if (X){ delete [] X; }
}

/*
  Call the underlying mesher with the default options
*/
void TMRMesh::mesh( double htarget ){
  TMRMeshOptions options;
  mesh(options, htarget);
}

/*
  Mesh the underlying geometry
*/
void TMRMesh::mesh( TMRMeshOptions options, double htarget ){
  // Mesh the curves
  int num_edges;
  TMREdge **edges;
  geo->getEdges(&num_edges, &edges);
  for ( int i = 0; i < num_edges; i++ ){
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);
    if (!mesh){
      mesh = new TMREdgeMesh(comm, edges[i]);
      mesh->mesh(options, htarget);
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
      mesh->mesh(options, htarget);
      faces[i]->setMesh(mesh);
    }
  }

  // Check if we can mesh the volume - if any
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
  int num_vertices;
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

  // Count up the number of quadrilaterals or hexes in the mesh
  num_quads = 0;
  num_hexes = 0;

  if (num_volumes > 0){
    for ( int i = 0; i < num_volumes; i++ ){
      TMRVolumeMesh *mesh = NULL;
      volumes[i]->getMesh(&mesh);
      num_hexes += mesh->getLocalConnectivity(NULL);
    }
  }
  if (num_faces > 0){
    for ( int i = 0; i < num_faces; i++ ){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);
      num_quads += mesh->getLocalConnectivity(NULL);
    }

    // Add the quadrilateral quality to the bins from each mesh
    const int nbins = 20;
    int bins[nbins];
    memset(bins, 0, nbins*sizeof(int));
    for ( int i = 0; i < num_faces; i++ ){
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);
      mesh->addQuadQuality(nbins, bins);
    }

    // Sum up the total number of elements
    int total = 0;
    for ( int i = 0; i < nbins; i++ ){
      total += bins[i];
    }

    // Get the MPI rank and print out the quality on the root processor
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    if (mpi_rank == 0){
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
void TMRMesh::initMesh(){
  // Allocate the global arrays
  X = new TMRPoint[ num_nodes ];

  if (num_hexes > 0){
    hexes = new int[ 8*num_hexes ];

    int *count = new int[ num_nodes ];
    memset(count, 0, num_nodes*sizeof(int));

    // Retrieve the surface information
    int num_volumes;
    TMRVolume **volumes;
    geo->getVolumes(&num_volumes, &volumes);

    // Set the values into the global arrays
    int *h = hexes;
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
      int nlocal = mesh->getLocalConnectivity(&hex_local);

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
        count[vars[j]]++;
      }
    }

    // Set the count to the number of variables
    for ( int j = 0; j < num_nodes; j++ ){
      if (count[j] == 0){ 
        printf("Node %d not referenced\n", j); 
      }
      if (count[j] >= 2){ 
        printf("Node %d referenced more than once %d\n", j, count[j]); 
      }
    }
    delete [] count;
  }
  if (num_quads > 0){
    quads = new int[ 4*num_quads ];

    // Retrieve the surface information
    int num_faces;
    TMRFace **faces;
    geo->getFaces(&num_faces, &faces);

    // Set the values into the global arrays
    int *q = quads;
    for ( int i = 0; i < num_faces; i++ ){
      // Get the mesh
      TMRFaceMesh *mesh = NULL;
      faces[i]->getMesh(&mesh);

      // Get the local mesh points
      int npts;
      TMRPoint *Xpts;
      mesh->getMeshPoints(&npts, NULL, &Xpts);

      // Get the local quadrilateral connectivity
      const int *quad_local;
      int nlocal = mesh->getLocalConnectivity(&quad_local);

      // Get the local to global variable numbering
      const int *vars;
      mesh->getNodeNums(&vars);

      // Set the quadrilateral connectivity
      for ( int j = 0; j < 4*nlocal; j++, q++ ){
        q[0] = vars[quad_local[j]];
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
void TMRMesh::getMeshConnectivity( int *_nquads, const int **_quads,
                                   int *_nhexes, const int **_hexes ){
  if (_quads || _hexes){ 
    if (!X){ initMesh(); }
    if (_nquads){ *_nquads = num_quads; }
    if (_quads){ *_quads = quads; }
    if (_nhexes){ *_nhexes = num_hexes; }
    if (_hexes){ *_hexes = hexes; }    
  }
}

/*
  Print out the mesh to a VTK file
*/
void TMRMesh::writeToVTK( const char *filename, int flag ){
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0 && num_nodes > 0){
    // Check whether to print out just quads or hexes or both
    int nquad = num_quads;
    int nhex = num_hexes;
    if (!(flag & TMR_QUADS)){
      nquad = 0;
    }
    if (!(flag & TMR_HEXES)){
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

      fprintf(fp, "\nCELLS %d %d\n", nquad + nhex, 5*nquad + 9*nhex);
    
      // Write out the cell connectivities
      for ( int k = 0; k < nquad; k++ ){
        fprintf(fp, "4 %d %d %d %d\n", 
                quads[4*k], quads[4*k+1], 
                quads[4*k+2], quads[4*k+3]);
      }
      for ( int k = 0; k < nhex; k++ ){
        fprintf(fp, "8 %d %d %d %d %d %d %d %d\n", 
                hexes[8*k], hexes[8*k+1], hexes[8*k+2], hexes[8*k+3],
                hexes[8*k+4], hexes[8*k+5], hexes[8*k+6], hexes[8*k+7]);
      }

      // All quadrilaterals
      fprintf(fp, "\nCELL_TYPES %d\n", nquad + nhex);
      for ( int k = 0; k < nquad; k++ ){
        fprintf(fp, "%d\n", 9);
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

      if (num_quads > 0 && (flag & TMR_QUADS)){
        int num_faces;
        TMRFace **faces;
        geo->getFaces(&num_faces, &faces);

        // Write out the element data for each component
        for ( int i = 0, j = 0; i < num_faces; i++ ){
          TMRFaceMesh *mesh = NULL;
          faces[i]->getMesh(&mesh);

          // Loop over all possible edges in the surface mesh
          const int *quad_local;
          int nlocal = mesh->getLocalConnectivity(&quad_local);

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

          for ( int k = 0; k < nlocal; k++, j++ ){
            int part = i+1;
            fprintf(fp, "%-8s%8d%8d%8d%8d%8d%8d%8d\n", "CQUADR", 
                    j+1, part, vars[quad_local[4*k]]+1, 
                    vars[quad_local[4*k+1]]+1, 
                    vars[quad_local[4*k+2]]+1, 
                    vars[quad_local[4*k+3]]+1, part);
          }
        }
      }
      if (num_hexes > 0 && (flag & TMR_HEXES)){
        int num_volumes;
        TMRVolume **volumes;
        geo->getVolumes(&num_volumes, &volumes);

        // Write out the element data for each component
        for ( int i = 0, j = 0; i < num_volumes; i++ ){
          TMRVolumeMesh *mesh = NULL;
          volumes[i]->getMesh(&mesh);

          // Loop over all possible edges in the surface mesh
          const int *hex_local;
          int nlocal = mesh->getLocalConnectivity(&hex_local);

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
            fprintf(fp, "%-8s%8d%8d%8d%8d%8d%8d%8d%8d*\n", "CHEXA", 
                    j+1, part, 
                    vars[hex_local[8*k]]+1, 
                    vars[hex_local[8*k+1]]+1, 
                    vars[hex_local[8*k+2]]+1, 
                    vars[hex_local[8*k+3]]+1,
                    vars[hex_local[8*k+4]]+1,
                    vars[hex_local[8*k+5]]+1);
            fprintf(fp, "%-8s%8d%8d\n", "*", 
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
      new_verts[vnum] = new TMRVertexFromFace(faces[i], pts[2*j], pts[2*j+1]);
    }
  }

  // Create all the remaining nodes. These are associated with the
  // TMRVolumeMesh since they do not lie on an edge or face
  if (num_hexes > 0){
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

  if (num_hexes > 0){
    // Compute the hexahedral edges and surfaces
    computeHexEdgesAndFaces(num_nodes, num_hexes, hexes,
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

  // Create a searchable array that stores the edge number
  // and the two connecting node numbers. The node numbers
  // are sorted such that the lowest one comes first. This
  // enables fast sorting/searching (based on a Morton ordering).
  int *all_edges = new int[ 3*num_mesh_edges ];

  // Populate the all_edges array, sorting each edge as we go...
  for ( int i = 0; i < num_mesh_edges; i++ ){
    all_edges[3*i] = i;

    // Set the quad edge, setting the lower edge number first
    if (mesh_edges[2*i] < mesh_edges[2*i+1]){
      all_edges[3*i+1] = mesh_edges[2*i];
      all_edges[3*i+2] = mesh_edges[2*i+1];
    }
    else {
      all_edges[3*i+1] = mesh_edges[2*i+1];
      all_edges[3*i+2] = mesh_edges[2*i];      
    }
  }

  // Sort all of the edges/nodes so that they can be easily searched
  qsort(all_edges, num_mesh_edges, 3*sizeof(int), compare_edges);

  // Keep track of the orientation of the mesh edges.
  // all_edge_dir[i] > 0 means that the ordering of the vertices 
  // in the all_edges array is consistent with the orientation
  // in the mesh. edge_dir[i] < 0 means the edge is flipped 
  // relative to the all_edges array.
  int *all_edge_dir = new int[ num_mesh_edges ];
  memset(all_edge_dir, 0, num_mesh_edges*sizeof(int));

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
      int *res = (int*)bsearch(edge, all_edges, num_mesh_edges,
                               3*sizeof(int), compare_edges);

      if (res){
        int edge_num = res[0];

        // Check whether the node ordering is consistent with the
        // edge orientation. If not, tag this edge as reversed.
        if (vars[j] < vars[j+1]){
          all_edge_dir[edge_num] = 1;
        }
        else {
          all_edge_dir[edge_num] = -1;
        }

        new_edges[edge_num] = 
          new TMRSplitEdge(edges[i], tpts[j], tpts[j+1]);
      }
      else {
        fprintf(stderr, "TMRMesh error: Could not find edge (%d, %d)\n",
                edge[1], edge[2]);
      }
    }
  }

  // Set the local quad node to edge information
  const int face_edge_to_nodes[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

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
    int nlocal = mesh->getLocalConnectivity(&quad_local);
    
    // Get the global variables associated with the local mesh
    const int *vars;
    mesh->getNodeNums(&vars);

    for ( int j = 0; j < nlocal; j++ ){
      for ( int k = 0; k < 4; k++ ){
        // The local variable numbers
        int l1 = quad_local[4*j + face_edge_to_nodes[k][0]];
        int l2 = quad_local[4*j + face_edge_to_nodes[k][1]];

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
        int *res = (int*)bsearch(edge, all_edges, num_mesh_edges,
                                 3*sizeof(int), compare_edges);

        if (res){
          // Get the edge number
          int edge_num = res[0];
          if (!new_edges[edge_num]){
            // These edges are constructed such that they are 
            // always in the 'positive' orientation.
            all_edge_dir[edge_num] = 1;

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
          fprintf(stderr, "TMRMesh error: Could not find edge (%d, %d)\n",
                  edge[1], edge[2]);
        }
      }
    } 
  }

  // Create the hexahedral elements
  if (num_hexes > 0){
    for ( int i = 0; i < num_mesh_edges; i++ ){
      if (!new_edges[i]){
        // Get the edges
        int v1 = all_edges[3*i+1];
        int v2 = all_edges[3*i+2];

        // Set the edge direction
        all_edge_dir[i] = 1;
        new_edges[i] = NULL; // new TMRLinearEdge(new_verts[v1], new_verts[v2]);
      }
    }
  }

  // Create the TMRFace objects
  TMRFace **new_faces = new TMRFace*[ num_mesh_faces ];
  memset(new_faces, 0, num_mesh_faces*sizeof(TMRFace*));

  // Create the face array
  int *all_faces = NULL;

  if (num_hexes > 0){
    // Allocate an array to store the searchable array that
    // stores the face nodes
    all_faces = new int[ 5*num_hex_faces ];

    for ( int i = 0; i < num_hex_faces; i++ ){
      all_faces[5*i] = i;
      for ( int k = 0; k < 4; k++ ){
        all_faces[5*i+1+k] = hex_faces[4*i+k];
      }
      sort_face_nodes(&all_faces[5*i+1]);
    }

    // Sort all of the faces so that they can be easily searched
    qsort(all_faces, num_mesh_faces, 
          5*sizeof(int), compare_faces);
  }

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
    int nlocal = mesh->getLocalConnectivity(&quad_local);
    
    // Get the global variables associated with the local mesh
    const int *vars;
    mesh->getNodeNums(&vars);

    // Iterate over all of the edges, creating the appropriate faces
    for ( int j = 0; j < nlocal; j++ ){
      TMREdge *c[4];
      int dir[4];
      TMRVertex *v[4];

      // Loop over all of the edges/nodes associated with this quad
      for ( int k = 0; k < 4; k++ ){
        // Set the vertex number
        int l0 = quad_local[4*j + k];
        v[k] = new_verts[vars[l0]];

        // The edge variable numbers
        int l1 = quad_local[4*j + face_edge_to_nodes[k][0]];
        int l2 = quad_local[4*j + face_edge_to_nodes[k][1]];

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
        int *res = (int*)bsearch(edge, all_edges, num_mesh_edges,
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
          dir[k] *= all_edge_dir[edge_num];
        } 
        else {
          fprintf(stderr, "TMRMesh error: Could not find edge (%d, %d)\n",
                  edge[1], edge[2]);
        }
      }

      // If this is a hexahedral mesh, then we need to be consistent
      // with how the faces are ordered. This code searches for the
      // face number within the all_faces array to obtain the required
      // face number
      if (all_faces){
        // Set the nodes associated with this face and sort them
        int face[5];
        face[0] = 0;
        for ( int k = 0; k < 4; k++ ){
          face[k+1] = vars[quad_local[4*j + k]];
        }
        sort_face_nodes(&face[1]);

        // Search for the face
        int *res = (int*)bsearch(face, all_faces, num_mesh_faces,
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

  // Array to store the new volume elements
  TMRVolume **new_volumes = NULL;

  if (num_hexes > 0){
    // Create the remaining faces from the mesh
    for ( int i = 0; i < num_mesh_faces; i++ ){
      if (!new_faces[i]){
        // new_faces[i] = new TMRTFIFace();
      }
    }
  }

  // Free all of the edge search information
  delete [] all_edges;
  if (all_faces){
    delete [] all_faces;
  }

  if (num_hexes > 0){
    // Create the new volume array
    new_volumes = new TMRVolume*[ num_hexes ];

    for ( int i = 0; i < num_hexes; i++ ){
      // Get the edges
      TMRVertex *v[8];
      TMREdge *e[12];
      TMRFace *f[6];
      int edir[12], fdir[6];

      // Get the vertices
      for ( int j = 0; j < 8; j++ ){
        v[j] = new_verts[hexes[8*i + j]];
      }

      // Get the edges and their directions
      for ( int j = 0; j < 12; j++ ){
        int edge_num = hex_edge_nums[12*i + j];
        edir[j] = 1;
        if (hexes[8*i + hex_edge_nodes[j][0]] >
            hexes[8*i + hex_edge_nodes[j][1]]){
          edir[j] = -1;
        }
        edir[j] *= all_edge_dir[edge_num];
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

  if (num_hexes > 0){
    geo = new TMRModel(num_nodes, new_verts,
                       num_mesh_edges, new_edges,
                       num_mesh_faces, new_faces,
                       num_hexes, new_volumes);
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

  return geo;
}
