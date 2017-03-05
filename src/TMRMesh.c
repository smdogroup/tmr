#include "TMRMesh.h"
#include "TMRNativeTopology.h"
#include "TMRTriangularize.h"
#include "TMRMeshSmoothing.h"
#include "TMRPerfectMatchInterface.h"
#include "TMRBspline.h"
#include <math.h>
#include <stdio.h>

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
  Compute a node to triangle or node to quad data structure
*/
static void computeNodeToElems( int nnodes, int nelems, 
                                int num_elem_nodes, const int conn[], 
                                int **_ptr, int **_node_to_elems ){
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
static void computeTriEdges( int nnodes, int ntris, const int tris[],
                             int *num_tri_edges, int **_tri_edges,
                             int **_tri_neighbors, int **_dual_edges ){
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
  delete [] ptr;
  delete [] node_to_tris;

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
static void computeQuadEdges( int nnodes, int nquads, const int quads[], 
                              int *_num_quad_edges, int **_quad_edges,
                              int **_quad_edge_nums=NULL ){
  // Compute the edges in the quadrilateral mesh
  int *ptr;
  int *node_to_quads;
  computeNodeToElems(nnodes, nquads, 4, quads, &ptr, &node_to_quads);
  
  // Now compute the neighbors for each quad
  int *quad_edge_nums = new int[ 4*nquads ];
  for ( int i = 0; i < 4*nquads; i++ ){
    quad_edge_nums[i] = -1;
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
  for ( int i = 0; i < nquads; i++ ){
    for ( int j = 0; j < 4; j++ ){
      // Get the unique edge number
      int n = quad_edge_nums[4*i+j];
      quad_edges[2*n] = quads[4*i+enodes[j][0]];
      quad_edges[2*n+1] = quads[4*i+enodes[j][1]];
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
}

/*
  Mesh a curve
*/
TMREdgeMesh::TMREdgeMesh( TMREdge *_edge ){
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
  Create a mesh
*/
void TMREdgeMesh::mesh( double htarget ){
  double tmin, tmax;
  edge->getRange(&tmin, &tmax);

  // Set the integration error tolerance
  double integration_eps = 1e-8;

  // Integrate along the curve to obtain the distance function such
  // that dist(tvals[i]) = int_{tmin}^{tvals[i]} ||d{C(t)}dt||_{2} dt
  int nvals;
  double *dist, *tvals;
  edge->integrate(tmin, tmax, integration_eps, &tvals, &dist, &nvals);
  
  // Compute the number of points along this curve
  npts = 1 + (int)(dist[nvals-1]/htarget);
  if (npts < 2){ npts = 2; }

  // If we have an even number of points, increment by one to ensure
  // that we have an even number of segments along the boundary
  if (npts % 2 != 1){ npts++; }

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

  // Allocate the points
  X = new TMRPoint[ npts ];
  for ( int i = 0; i < npts; i++ ){
    edge->evalPoint(pts[i], &X[i]);
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
TMRFaceMesh::TMRFaceMesh( TMRFace *_face ){
  face = _face;
  face->incref();

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

static int tri_mesh_count = 0;

/*
  Create the surface mesh
*/
void TMRFaceMesh::mesh( double htarget ){
  // Count up the number of points and segments from the 
  // curves that bound the surface. Keep track of the number 
  // of points = the number of segments.
  int total_num_pts = 0;

  // Keep track of the number of closed loop cycles in the domain
  int nloops = face->getNumEdgeLoops();

  // Get all of the meshes
  for ( int k = 0; k < nloops; k++ ){
    TMREdgeLoop *loop;
    face->getEdgeLoop(k, &loop);
    int nedges;
    TMREdge **edges;
    loop->getEdgeLoop(&nedges, &edges, NULL);

    for ( int i = 0; i < nedges; i++ ){
      TMREdgeMesh *mesh = NULL;
      edges[i]->getMesh(&mesh);
      if (!mesh){
        mesh = new TMREdgeMesh(edges[i]);
        mesh->mesh(htarget);
        edges[i]->setMesh(mesh);
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
  int seg = 0;
  int pt = 0;

  int init_loop_pt = 0; // What point value did this loop start on? 
  int init_loop_seg = 0; // What segment value did this loop start on?
  int hole_pt = total_num_pts; // What hole are we on?

  for ( int k = 0; k < nloops; k++ ){
    // Set the offset to the initial point/segment on this loop
    init_loop_pt = pt;
    init_loop_seg = seg;

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
          segments[2*seg] = pt;
          segments[2*seg+1] = pt+1;
          seg++;
          pt++;
        }
      }
      else {
        // Reverse the parameter values on the edge
        for ( int j = npts-1; j >= 1; j-- ){
          edge->getParamsOnFace(face, tpts[j], dir[i],
                                &params[2*pt], &params[2*pt+1]);
          segments[2*seg] = pt;
          segments[2*seg+1] = pt+1;
          seg++;
          pt++;
        }
      }
    }

    // Close off the loop by connecting the segment back to the
    // initial loop point
    segments[2*(seg-1)+1] = init_loop_pt;

    /***************************************************************
    // Compute the area enclosed by the loop. If the area is
    // positive, it is the domain boundary. If the area is negative,
    // we have a hole!  Note that this assumes that the polygon
    // creating the hole is not self-intersecting. (In reality we
    // compute twice the area since we omit the 1/2 factor.)
    double Area = 0.0;
    for ( int i = init_loop_seg; i < seg; i++ ){
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
      int s1 = segments[2*init_loop_seg];
      int s2 = segments[2*init_loop_seg+1];
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
    ****************************************************************/
  }

  // Set the total number of fixed points. These are the points that
  // will not be smoothed and constitute the boundary nodes. Note that
  // the Triangularize class removes the holes from the domain
  // automatically.  The boundary points are guaranteed to be ordered
  // first.
  num_fixed_pts = total_num_pts;
  int num_smoothing_steps = 10;

  // Create the triangularization class
  TMRTriangularize *tri = 
    new TMRTriangularize(total_num_pts + nholes, params, nholes,
                         nsegs, segments, face);
  tri->incref();

  // Create the mesh using the frontal algorithm
  tri->frontal(htarget);

  // Extract the triangularization of the domain
  int ntris, *tris;
  tri->getMesh(&num_points, &ntris, &tris, &pts, &X);

  // Free the triangle mesh
  tri->decref();

  if (ntris > 0){
    // Compute the triangle edges and neighbors in the dual mesh
    int num_tri_edges;
    int *tri_edges, *tri_neighbors, *dual_edges;
    computeTriEdges(num_points, ntris, tris, 
                    &num_tri_edges, &tri_edges,
                    &tri_neighbors, &dual_edges);
    
    // Smooth the resulting triangular mesh
    laplacianSmoothing(10*num_smoothing_steps, num_fixed_pts,
                       num_tri_edges, tri_edges,
                       num_points, pts, X, face);

    // Recombine the mesh into a quadrilateral mesh
    recombine(ntris, tris, tri_neighbors,
              num_tri_edges, dual_edges, X, &num_quads, &quads);
    
    // Free the triangular mesh data
    delete [] tris;
    delete [] tri_neighbors;
    delete [] dual_edges;

    // Simplify the new quadrilateral mesh
    simplifyQuads();
    
    // Print the quadrilateral mesh quality
    // printQuadQuality();
    
    // Compute the quadrilateral mesh
    // int num_quad_edges;
    // int *quad_edges;
    // computeQuadEdges(num_points, num_quads, quads,
    //                  &num_quad_edges, &quad_edges);
    //
    // Smooth the quadrilateral mesh using Laplacian smoothing
    // laplacianSmoothing(10*num_smoothing_steps, num_fixed_pts,
    //                    num_quad_edges, quad_edges,
    //                    num_points, pts, X, surface);
    //
    // Smooth the quad mesh using a spring analogy
    // double alpha = 0.1;
    // springQuadSmoothing(10*num_smoothing_steps, alpha, num_fixed_pts,
    //                     num_quads, quads, num_quad_edges, quad_edges,
    //                     num_points, pts, X, surface);
    // delete [] quad_edges;

    int *ptr;
    int *pts_to_quads;
    computeNodeToElems(num_points, num_quads, 4, quads, &ptr, &pts_to_quads);

    // Smooth the mesh using a local optimization of node locations
    num_smoothing_steps = 5;
    quadSmoothing(10*num_smoothing_steps, num_fixed_pts,
                  num_points, ptr, pts_to_quads, num_quads, quads, 
                  pts, X, face);

    // Free the connectivity information
    delete [] ptr;
    delete [] pts_to_quads;

    // Print the quadrilateral mesh quality
    printQuadQuality();
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
                             int num_edges, const int dual_edges[],
                             const TMRPoint *p, int *_num_quads, 
                             int **_new_quads ){
  // Count up the number of dual edges
  int num_dual_edges = 0;
  for ( int i = 0; i < num_edges; i++ ){
    if (dual_edges[2*i] >= 0 && dual_edges[2*i+1] >= 0){
      num_dual_edges++;
    }
  }

  // Allocate the weights
  double *weights = new double[ num_dual_edges ];
  int *graph_edges = new int[ 2*num_dual_edges ];

  // Compute the weight associated with each edge by combputing the
  // recombined quality
  for ( int i = 0, e = 0; i < num_edges; i++ ){
    int t1 = dual_edges[2*i];
    int t2 = dual_edges[2*i+1];

    if (t1 >= 0 && t2 >= 0){
      // Compute the weight for this recombination
      double quality = 
        computeRecombinedQuality(tris, tri_neighbors,
                                 t1, t2, p);

      double weight = (1.0 - quality); // *(1.0 + 1.0/(quality + 0.1));
      graph_edges[2*e] = t1;
      graph_edges[2*e+1] = t2;
      weights[e] = weight;
      e++;
    }
  }

  writeDualToVTK("dual_mesh.vtk", ntris, tris,
                 num_dual_edges, graph_edges, p);

  int *match = new int[ 2*ntris ];
  TMR_PerfectMatchGraph(ntris, num_dual_edges, 
                        graph_edges, weights, match);

  delete [] weights;
  delete [] graph_edges;

  // Set the number of quadrilateral elements created in the mesh and record
  // the new quadrilateral element connectivity
  int num_new_quads = ntris/2;
  int *new_quads = new int[ 4*num_new_quads ];

  // Recombine the quads
  for ( int i = 0; i < num_new_quads; i++ ){
    int fail = getRecombinedQuad(tris, tri_neighbors,
                                 match[2*i], match[2*i+1],
                                 &new_quads[4*i]);
    if (fail){
      printf("Recombined quadrilateral %d between triangles %d and %d failed\n", 
             i, match[2*i], match[2*i+1]);
    }
  }

  delete [] match;

  // Set the quads/output
  *_num_quads = num_new_quads;
  *_new_quads = new_quads;
}

/*
  Identify and remove adjacent quadrilaterals that look like this:

  x --- x             x ---- x
  | \   |             |      |
  |  x  |     ===>    |      |
  |    \|             |      |
  x --- x             x ---- x
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

      // Reset the quadrilateral by removing the
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
  Print the quad quality
*/
void TMRFaceMesh::printQuadQuality(){
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins*sizeof(int));
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
                                  int ntris, const int tris[],
                                  int num_dual_edges, const int dual_edges[],
                                  const TMRPoint *p ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    
    // Write out the points
    fprintf(fp, "POINTS %d float\n", ntris);
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "%e %e %e\n", 
        1.0/3.0*(p[tris[3*k]].x + p[tris[3*k+1]].x + p[tris[3*k+2]].x),
        1.0/3.0*(p[tris[3*k]].y + p[tris[3*k+1]].y + p[tris[3*k+2]].y),
        1.0/3.0*(p[tris[3*k]].z + p[tris[3*k+1]].z + p[tris[3*k+2]].z));
    }

    // Write out the cell connectivity
    fprintf(fp, "\nCELLS %d %d\n", num_dual_edges, 3*num_dual_edges);
    for ( int k = 0; k < num_dual_edges; k++ ){
      fprintf(fp, "2 %d %d\n", dual_edges[2*k], dual_edges[2*k+1]);
    }

    // Write out the cell types
    fprintf(fp, "\nCELL_TYPES %d\n", num_dual_edges);
    for ( int k = 0; k < num_dual_edges; k++ ){
      fprintf(fp, "%d\n", 3);
    }
    fclose(fp);
  }
}

/*
  Mesh the given geometry and retrieve either a regular mesh
*/
TMRMesh::TMRMesh( TMRModel *_geo ){
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
  Mesh the underlying geometry
*/
void TMRMesh::mesh( double htarget ){
  // Mesh the curves
  int num_edges;
  TMREdge **edges;
  geo->getEdges(&num_edges, &edges);
  for ( int i = 0; i < num_edges; i++ ){
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);
    if (!mesh){
      mesh = new TMREdgeMesh(edges[i]);
      mesh->mesh(htarget);
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
      mesh = new TMRFaceMesh(faces[i]);
      mesh->mesh(htarget);
      faces[i]->setMesh(mesh);
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

  // Set the number of nodes in the mesh
  num_nodes = num; 

  // Count up the number of quadrilaterals in the mesh
  num_quads = 0;
  for ( int i = 0; i < num_faces; i++ ){
    TMRFaceMesh *mesh = NULL;
    faces[i]->getMesh(&mesh);
    num_quads += mesh->getLocalConnectivity(NULL);
  }
}

/*
  Allocate and initialize the global mesh using the global ordering
*/
void TMRMesh::initMesh(){
  // Allocate the global arrays
  X = new TMRPoint[ num_nodes ];
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
    for ( int j = 0; j < 4*nlocal; j++ ){
      q[0] = vars[quad_local[j]];
      q++;
    }

    // Set the node locations
    for ( int j = 0; j < npts; j++ ){
      X[vars[j]] = Xpts[j];
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
int TMRMesh::getMeshConnectivity( const int **_quads ){
  if (_quads){ 
    if (!quads){ initMesh(); }
    *_quads = quads; 
  }
  return num_quads;
}

/*
  Create the topology object, generating the vertices, edges and
  faces for each element in the underlying mesh.
*/
TMRModel* TMRMesh::createMeshGeometry(){
  // Initialize the mesh
  initMesh();

  // Compute the quadrilateral edges
  int num_quad_edges;
  int *quad_edges;
  computeQuadEdges(num_nodes, num_quads, quads,
                   &num_quad_edges, &quad_edges);

  int *all_edges = new int[ 3*num_quad_edges ];
  for ( int i = 0; i < num_quad_edges; i++ ){
    all_edges[3*i] = i;

    // Set the quad edge, setting the lower edge number first
    if (quad_edges[2*i] < quad_edges[2*i+1]){
      all_edges[3*i+1] = quad_edges[2*i];
      all_edges[3*i+2] = quad_edges[2*i+1];
    }
    else {
      all_edges[3*i+1] = quad_edges[2*i+1];
      all_edges[3*i+2] = quad_edges[2*i];      
    }
  }

  // Sort all of the edges/nodes so that they can be easily searched
  qsort(all_edges, num_quad_edges, 3*sizeof(int), compare_edges);

  // Create vertices
  int vnum = 0;
  TMRVertex **verts = new TMRVertex*[ num_nodes ];

  // Copy over the vertices in the original geometry
  int num_vertices;
  TMRVertex **vertices;
  geo->getVertices(&num_vertices, &vertices);
  for ( int i = 0; i < num_vertices; i++, vnum++ ){
    verts[vnum] = vertices[i];
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
      verts[vnum] = new TMRVertexFromEdge(edges[i], tpts[j]);
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
      verts[vnum] = new TMRVertexFromFace(faces[i], pts[2*j], pts[2*j+1]);
    }
  }

  // Create the edges on the surface and on the edge
  TMREdge **new_edges = new TMREdge*[ num_quad_edges ];
  memset(new_edges, 0, num_quad_edges*sizeof(TMREdge*));

  // Create the edges for the mesh from the underlying edges
  for ( int i = 0; i < num_edges; i++ ){
    TMREdgeMesh *mesh = NULL;
    edges[i]->getMesh(&mesh);

    // Get the global variables associated with the local curve mesh
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
      int *res = (int*)bsearch(edge, all_edges, num_quad_edges,
                               3*sizeof(int), compare_edges);

      if (res){
        int edge_num = res[0];
        // Tag the edge numbers as negative
        if (quad_edges[2*edge_num] != vars[j] ||
            quad_edges[2*edge_num+1] != vars[j+1]){
          quad_edges[2*edge_num] *= -1;
          quad_edges[2*edge_num+1] *= -1;
        }
        new_edges[edge_num] = new TMRSplitEdge(edges[i], tpts[j], tpts[j+1]);
      }
    }
  }

  // Set the local quad node to edge information
  const int edge_to_nodes[][2] = {{0, 3}, {1, 2}, {0, 1}, {3, 2}};

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
        int l1 = quad_local[4*j + edge_to_nodes[k][0]];
        int l2 = quad_local[4*j + edge_to_nodes[k][1]];

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
        int *res = (int*)bsearch(edge, all_edges, num_quad_edges,
                                 3*sizeof(int), compare_edges);

        if (res){
          // Get the edge number
          int edge_num = res[0];
          if (!new_edges[edge_num]){
            // Create the TMRBsplinePcurve on this edge
            double cpts[4];
            if (vars[l1] == quad_edges[2*edge_num] && 
                vars[l2] == quad_edges[2*edge_num+1]){
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
      }
    } 
  }  

  // Create the surface 
  TMRFace **surfs = new TMRFace*[ num_quads ];

  // Switch the ordering from right-handed to coordinate ordering of the
  // vertices in a quadrilateral
  const int quad_to_coordinate[] = {0, 1, 3, 2};

  // Create the TMRFace objects
  for ( int i = 0, q = 0; i < num_faces; i++ ){
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
    for ( int j = 0; j < nlocal; j++, q++ ){
      TMREdge *c[4];
      int dir[4];
      TMRVertex *v[4];

      // Loop over all of the edges/nodes associated with this vertex
      for ( int k = 0; k < 4; k++ ){
        // Set the vertex number
        int l0 = quad_local[4*j + quad_to_coordinate[k]];
        v[k] = verts[vars[l0]];

        // The local variable numbers
        int l1 = quad_local[4*j + edge_to_nodes[k][0]];
        int l2 = quad_local[4*j + edge_to_nodes[k][1]];

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
        int *res = (int*)bsearch(edge, all_edges, num_quad_edges,
                                 3*sizeof(int), compare_edges);

        if (res){
          // Get the global edge number
          int edge_num = res[0];
          c[k] = new_edges[edge_num];

          // Figure out the edge direction
          int dir_fact = 1;
          int e1 = quad_edges[2*edge_num];
          int e2 = quad_edges[2*edge_num+1];
          if (e1 < 0 || e2 < 0){
            e1 *= -1;
            e2 *= -1;
            dir_fact = -1;
          }
          if (vars[l1] == e1 && vars[l2] == e2){
            dir[k] = dir_fact;
          }
          else {
            dir[k] = -dir_fact;
          }
        }
      }

      // Create the parametric TFI surface
      surfs[q] = new TMRParametricTFIFace(faces[i], c, dir, v);
    }
  } 

  // Create the geometry object
  TMRModel *geo = new TMRModel(num_nodes, verts,
                               num_quad_edges, new_edges,
                               num_quads, surfs);

  // Free the data that was locally allocated
  delete [] verts;
  delete [] new_edges;
  delete [] surfs;
  delete [] all_edges;
  delete [] quad_edges;

  return geo;
}

