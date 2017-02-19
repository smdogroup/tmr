#include "TMRMesh.h"
#include "TMRTriangularize.h"
#include "TMRPerfectMatchInterface.h"
#include <math.h>
#include <stdio.h>

/*
  Add the motion of the point in parameter space
*/
inline void addParamMovement( double alpha,
                              TMRPoint *Xu, TMRPoint *Xv,
                              TMRPoint *d, double *delta ){
  // Compute the sum - in the parameter space - of the motion 
  // contributed by each node
  double g11 = Xu->dot(Xu);
  double g12 = Xu->dot(Xv);
  double g22 = Xv->dot(Xv);

  // Compute the right hand sides
  double b1 = d->dot(Xu);
  double b2 = d->dot(Xv);

  // Compute the inverse of the determinant
  double invdet = alpha/(g11*g22 - g12*g12);

  // Add the contribution
  delta[0] += invdet*(g22*b1 - g12*b2);
  delta[1] += invdet*(g11*b2 - g12*b1);
}

/*
  Smoothing methods for the different points
*/
static void laplacianSmoothing( int nsmooth, int num_fixed_pts,
                                int num_edges, const int *edge_list,
                                int num_pts, double *prm, TMRPoint *p,
                                TMRSurface *surface );
static void springSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                             int num_edges, const int *edge_list,
                             int num_pts, double *prm, TMRPoint *p,
                             TMRSurface *surface );
static void springQuadSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                                 int num_quads, const int *quad_list,
                                 int num_edges, const int *edge_list,
                                 int num_pts, double *prm, TMRPoint *p,
                                 TMRSurface *surface );
static void quadSmoothing( int nsmooth, int num_fixed_pts,
                           int num_pts, const int *ptr, const int *pts_to_quads,
                           int num_quads, const int *quad_list,
                           double *prm, TMRPoint *p,
                           TMRSurface *surface );

/*
  Mesh a curve
*/
TMRCurveMesh::TMRCurveMesh( TMRCurve *_curve ){
  curve = _curve;
  curve->incref();

  npts = 0;
  pts = NULL;
  X = NULL;
}

/*
  Destroy the mesh for this curve, and free the underlying data
*/
TMRCurveMesh::~TMRCurveMesh(){
  curve->decref();
  if (pts){ delete [] pts; }
  if (X){ delete [] X; }
}

/*
  Create a mesh
*/
void TMRCurveMesh::mesh( double htarget ){
  double tmin, tmax;
  curve->getRange(&tmin, &tmax);

  // Set the integration error tolerance
  double integration_eps = 1e-8;

  // Integrate along the curve to obtain the distance function such
  // that dist(tvals[i]) = int_{tmin}^{tvals[i]} ||d{C(t)}dt||_{2} dt
  int nvals;
  double *dist, *tvals;
  curve->integrate(tmin, tmax, integration_eps, &tvals, &dist, &nvals);
  
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
    curve->evalPoint(pts[i], &X[i]);
  }
}

/*
  Get the mesh points
*/
void TMRCurveMesh::getMesh( int *_npts, const double **_pts, 
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
TMRSurfaceMesh::TMRSurfaceMesh( TMRSurface *surf ){
  surface = surf;
  surface->incref();
}

/*
  Free the data associated with the surface mesh
*/
TMRSurfaceMesh::~TMRSurfaceMesh(){
  surface->decref();
}

/*
  Create the surface mesh
*/
void TMRSurfaceMesh::mesh( double htarget ){
  // Count up the number of points and segments from the curves that
  // bound the surface

  // Keep track of the number of points = the number of segments.
  int total_num_pts = 0;

  // Keep track of the number of closed loop cycles in the domain
  int ncycles = surface->getNumSegments();

  // Get all of the meshes
  for ( int k = 0; k < surface->getNumSegments(); k++ ){
    int ncurves;
    TMRCurve **curves;
    surface->getCurveSegment(k, &ncurves, &curves, NULL);

    for ( int i = 0; i < ncurves; i++ ){
      TMRCurveMesh *mesh = NULL;
      curves[i]->getMesh(&mesh);
      if (!mesh){
        mesh = new TMRCurveMesh(curves[i]);
        mesh->mesh(htarget);
        curves[i]->setMesh(mesh);
      }
      
      // Get the number of points associated with the curve
      int npts;
      mesh->getMesh(&npts, NULL, NULL);
      
      // Update the total number of points
      total_num_pts += npts-1;
    }
  }

  // The number of holes is equal to the number of cycles-1. One loop
  // bounds the domain, the other loops cut out holes in the domain. 
  // Note that the domain must be contiguous.
  int nholes = ncycles-1;

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

  for ( int k = 0; k < surface->getNumSegments(); k++ ){
    // Set the offset to the initial point/segment on this loop
    init_loop_pt = pt;
    init_loop_seg = seg;

    // Get the curve information for this loop segment
    int ncurves;
    TMRCurve **curves;
    const int *dir;
    surface->getCurveSegment(k, &ncurves, &curves, &dir);

    for ( int i = 0; i < ncurves; i++ ){
      // Retrieve the underlying curve mesh
      TMRCurveMesh *mesh = NULL;
      curves[i]->getMesh(&mesh);
      
      // Get the mesh points corresponding to this curve
      int npts;
      const double *tpts;
      mesh->getMesh(&npts, &tpts, NULL);
      
      // Find the point on the curve
      if (dir[i] > 0){
        for ( int j = 0; j < npts-1; j++ ){
          curves[i]->getParamsOnSurface(surface, tpts[j], dir[i],
                                        &params[2*pt], &params[2*pt+1]);
          segments[2*seg] = pt;
          segments[2*seg+1] = pt+1;
          seg++;
          pt++;
        }
      }
      else {
        for ( int j = npts-1; j >= 1; j-- ){
          curves[i]->getParamsOnSurface(surface, tpts[j], dir[i],
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
                         nsegs, segments, surface);
  tri->incref();

  // Create the mesh using the frontal algorithm
  tri->frontal(htarget);

  // Extract the triangularization of the domain
  int ntris, *tris;
  tri->getMesh(&num_points, &ntris, &tris, &pts, &X);

  // Free the triangle mesh
  tri->decref();

  // Compute the triangle edges and neighbors in the dual mesh
  int num_tri_edges;
  int *tri_edges, *tri_neighbors, *dual_edges;
  computeTriEdges(num_points, ntris, tris, 
                  &num_tri_edges, &tri_edges,
                  &tri_neighbors, &dual_edges);

  // Smooth the resulting triangular mesh
  laplacianSmoothing(10*num_smoothing_steps, num_fixed_pts,
                     num_tri_edges, tri_edges,
                     num_points, pts, X, surface);

  writeTrisToVTK("bspline_tris.vtk", ntris, tris);

  printTriQuality(ntris, tris);

  // Recombine the mesh into a quadrilateral mesh
  recombine(ntris, tris, tri_neighbors,
            num_tri_edges, dual_edges, X, &num_quads, &quads);

  // Free the triangular mesh data
  delete [] tris;
  delete [] tri_neighbors;
  delete [] dual_edges;

  // Print the quadrilateral mesh quality
  printQuadQuality();

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
  quadSmoothing(10*num_smoothing_steps, num_fixed_pts,
                num_points, ptr, pts_to_quads, num_quads, quads, 
                pts, X, surface);

  // Free the connectivity information
  delete [] ptr;
  delete [] pts_to_quads;

  // Print the quadrilateral mesh quality
  printQuadQuality();
}

/*
  Compute a node to triangle or node to quad data structure
*/
void TMRSurfaceMesh::computeNodeToElems( int nnodes, int nelems, 
                                         int num_elem_nodes,
                                         const int conn[], 
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
void TMRSurfaceMesh::computeTriEdges( int nnodes, int ntris, const int tris[],
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
void TMRSurfaceMesh::computeQuadEdges( int nnodes, int nquads, 
                                       const int quads[],
                                       int *_num_quad_edges,
                                       int **_quad_edges ){
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
  delete [] quad_edge_nums;

  *_num_quad_edges = ne;
  *_quad_edges = quad_edges;
}

/*
  Get the quadrilateral elemnet obtained by combining the triangles
  t1 and t2 together

        . --- .
      / |   /
    /   | /
  . --- . 
*/
int TMRSurfaceMesh::getRecombinedQuad( const int tris[],
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
double TMRSurfaceMesh::computeQuadQuality( const int *quad,
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
double TMRSurfaceMesh::computeRecombinedQuality( const int tris[], 
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
double TMRSurfaceMesh::computeTriQuality( const int *tri,
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
void TMRSurfaceMesh::recombine( int ntris, const int tris[],
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
      double quality = computeRecombinedQuality(tris, tri_neighbors,
                                                t1, t2, p);

      double weight = (1.0 - quality)*(1.0 + 1.0/(quality + 0.1));
      graph_edges[2*e] = t1;
      graph_edges[2*e+1] = t2;
      weights[e] = weight;
      e++;
    }
  }

  int *match = new int[ 2*ntris ];
  TMR_PerfectMatchGraph(ntris, num_dual_edges, graph_edges, weights, match);

  delete [] weights;
  delete [] graph_edges;

  // Set the number of quadrilateral elements created in the mesh and record
  // the new quadrilateral element connectivity
  int num_new_quads = ntris/2;
  int *new_quads = new int[ 4*num_new_quads ];

  // Recombine the quads
  for ( int i = 0; i < num_new_quads; i++ ){
    int fail = getRecombinedQuad(tris, tri_neighbors,
                                 match[2*i], match[2*i+1], &new_quads[4*i]);
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
  Print the quad quality
*/
void TMRSurfaceMesh::printQuadQuality(){
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
void TMRSurfaceMesh::printTriQuality( int ntris,
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
void TMRSurfaceMesh::writeToVTK( const char *filename ){
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
  Write the output to a VTK file
*/
void TMRSurfaceMesh::writeTrisToVTK( const char *filename,
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
  Apply Laplacian smoothing
*/
void laplacianSmoothing( int nsmooth, int num_fixed_pts,
                         int num_edges, const int *edge_list,
                         int num_pts, double *prm, TMRPoint *p,
                         TMRSurface *surface ){
  int *count = new int[ num_pts ];
  double *new_params = new double[ 2*num_pts ];
  TMRPoint *Xu = new TMRPoint[ num_pts ];
  TMRPoint *Xv = new TMRPoint[ num_pts ];

  for ( int iter = 0; iter < nsmooth; iter++ ){
    memset(count, 0, num_pts*sizeof(int));
    memset(new_params, 0, 2*num_pts*sizeof(double));

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
    }    

    // Loop over all the edges
    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      // Add the movement of the node in parameter space
      if (n1 >= num_fixed_pts){
        addParamMovement(1.0, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
        count[n1]++;
      }

      // Add the movement of the second node in parameter space
      if (n2 >= num_fixed_pts){
        addParamMovement(-1.0, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
        count[n2]++;
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      if (count[i] > 0){
        prm[2*i] += new_params[2*i]/count[i];
        prm[2*i+1] += new_params[2*i+1]/count[i];
        surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
      }
    }
  }

  delete [] count;
  delete [] new_params;
  delete [] Xu;
  delete [] Xv;
}

/*
  Apply the spring smoothing analogy
*/
void springSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                      int num_edges, const int *edge_list,
                      int num_pts, double *prm, TMRPoint *p,
                      TMRSurface *surface ){
  double *len = new double[ num_edges ];
  double *new_params = new double[ 2*num_pts ];
  TMRPoint *Xu = new TMRPoint[ num_pts ];
  TMRPoint *Xv = new TMRPoint[ num_pts ];

  for ( int iter = 0; iter < nsmooth; iter++ ){
    double sum = 0.0;
    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      len[i] = sqrt(d.dot(d));
      sum += len[i];
    }
    double len0 = 0.9*sum/num_edges;

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
    }

    memset(new_params, 0, 2*num_pts*sizeof(double));

    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      double scale = (len0 - len[i])/len[i];

      // Add the movement of the node in parameter space
      if (n1 >= num_fixed_pts){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }

      // Add the movement of the second node in parameter space
      if (n2 >= num_fixed_pts){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      prm[2*i] += alpha*new_params[2*i];
      prm[2*i+1] += alpha*new_params[2*i+1];
      surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
    }
  }

  delete [] new_params;
  delete [] len;
  delete [] Xu;
  delete [] Xv;
}

/*
  Perform spring smoothing specifically for a quadrilateral mesh.

  This code also connects springs across the faces of quadrilateral
  elements in an attempt to achieve higher mesh quality. This is 
  not always successful.
*/
void springQuadSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                          int num_quads, const int *quad_list,
                          int num_edges, const int *edge_list,
                          int num_pts, double *prm, TMRPoint *p,
                          TMRSurface *surface ){
  double *len = new double[ num_edges ];
  double *new_params = new double[ 2*num_pts ];
  TMRPoint *Xu = new TMRPoint[ num_pts ];
  TMRPoint *Xv = new TMRPoint[ num_pts ];

  // Compute the squrare root of 2
  const double sqrt2 = 1.1*sqrt(2.0);

  for ( int iter = 0; iter < nsmooth; iter++ ){
    double sum = 0.0;
    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      len[i] = sqrt(d.dot(d));
      sum += len[i];
    }
    double len0 = sum/num_edges;

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
    }

    memset(new_params, 0, 2*num_pts*sizeof(double));

    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      double scale = (len0 - len[i])/len[i];

      // Add the movement of the node in parameter space
      if (n1 >= num_fixed_pts){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }

      // Add the movement of the second node in parameter space
      if (n2 >= num_fixed_pts){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    for ( int i = 0; i < num_quads; i++ ){
      // Add the contribution from the first cross-quad member
      int n1 = quad_list[4*i];
      int n2 = quad_list[4*i+2];

      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;
      double ld = sqrt(d.dot(d));
      double scale = (sqrt2*len0 - ld)/ld;
      if (n1 >= num_fixed_pts){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }
      if (n2 >= num_fixed_pts){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }

      // Add the contribution from the second cross-quad member
      n1 = quad_list[4*i+1];
      n2 = quad_list[4*i+3];

      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;
      ld = sqrt(d.dot(d));
      scale = (sqrt2*len0 - ld)/ld;
      if (n1 >= num_fixed_pts){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }
      if (n2 >= num_fixed_pts){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      prm[2*i] += alpha*new_params[2*i];
      prm[2*i+1] += alpha*new_params[2*i+1];
      surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
    }
  }

  delete [] new_params;
  delete [] len;
  delete [] Xu;
  delete [] Xv;
}

/*
  Smooth the mesh
*/
void quadSmoothing( int nsmooth, int num_fixed_pts,
                    int num_pts, const int *ptr, const int *pts_to_quads,
                    int num_quads, const int *quads,
                    double *prm, TMRPoint *p,
                    TMRSurface *surface ){
  for ( int iter = 0; iter < nsmooth; iter++ ){
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      // Evaluate the derivatives w.r.t. the parameter locations so
      // that we can take movement in the physical plane and convert
      // it to movement in the parametric coordinates
      TMRPoint Xu, Xv;
      surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu, &Xv);

      // Normalize the directions Xu, Xv to form a locally-orthonormal 
      // coordinate frame aligned with the surface
      TMRPoint xdir, ydir;

      // Normalize the x-direction
      double xnorm = sqrt(Xu.dot(Xu));
      xdir.x = Xu.x/xnorm;
      xdir.y = Xu.y/xnorm;
      xdir.z = Xu.z/xnorm;

      // Remove the component of the x-direction from Xv
      double dot = xdir.dot(Xv);
      ydir.x = Xv.x - dot*xdir.x;
      ydir.y = Xv.y - dot*xdir.y;
      ydir.z = Xv.z - dot*xdir.z;

      double ynorm = sqrt(ydir.dot(ydir));
      ydir.x = ydir.x/ynorm;
      ydir.y = ydir.y/ynorm;
      ydir.z = ydir.z/ynorm;

      int N = ptr[i+1] - ptr[i];
      if (N > 0){
        // Loop over the quadrilaterals that reference this point 
        double A = 0.0, B = 0.0;  
        for ( int qp = ptr[i]; qp < ptr[i+1]; qp++ ){
          const int *quad = &quads[4*pts_to_quads[qp]];

          // Pick out the influence triangle points from the quadrilateral 
          // This consists of the base point i and the following two
          int ijk[3];
          if (quad[0] == i){
            ijk[0] = quad[0];  ijk[1] = quad[1];  ijk[2] = quad[3];
          }
          else if (quad[1] == i){
            ijk[0] = quad[1];  ijk[1] = quad[2];  ijk[2] = quad[0];
          }
          else if (quad[2] == i){
            ijk[0] = quad[2];  ijk[1] = quad[3];  ijk[2] = quad[1];
          }
          else {
            ijk[0] = quad[3];  ijk[1] = quad[0];  ijk[2] = quad[2];
          }

          // Now compute the geometric quantities
          // p = yj - yk, q = xk - xj
          double xi = xdir.dot(p[ijk[0]]);
          double yi = ydir.dot(p[ijk[0]]);
          double xj = xdir.dot(p[ijk[1]]);
          double yj = ydir.dot(p[ijk[1]]);
          double xk = xdir.dot(p[ijk[2]]);
          double yk = ydir.dot(p[ijk[2]]);
          double p = yj - yk;
          double q = xk - xj;
          double r = xj*yk - xk*yj;
          double a = 0.5*(p*xi + q*yi + r);
          double b = sqrt(p*p + q*q); 
          A += a;
          B += b;
        }

        double hbar = 2.0*A/B;
        double bbar = B/N;

        // Set the weights
        double w1 = 1.0/(hbar*hbar);
        double w2 = 4.0/(bbar*bbar);

        // The parameters for the Jacobian/right-hand-side
        double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, s5 = 0.0;

        for ( int qp = ptr[i]; qp < ptr[i+1]; qp++ ){
          const int *quad = &quads[4*pts_to_quads[qp]];

          // Pick out the influence triangle points from the quadrilateral 
          // This consists of the base point i and the following two
          int ijk[3];
          if (quad[0] == i){
            ijk[0] = quad[0];  ijk[1] = quad[1];  ijk[2] = quad[3];
          }
          else if (quad[1] == i){
            ijk[0] = quad[1];  ijk[1] = quad[2];  ijk[2] = quad[0];
          }
          else if (quad[2] == i){
            ijk[0] = quad[2];  ijk[1] = quad[3];  ijk[2] = quad[1];
          }
          else {
            ijk[0] = quad[3];  ijk[1] = quad[0];  ijk[2] = quad[2];
          }

          // Now compute the geometric quantities
          // p = yj - yk, q = xk - xj
          double xi = xdir.dot(p[ijk[0]]);
          double yi = ydir.dot(p[ijk[0]]);
          double xj = xdir.dot(p[ijk[1]]);
          double yj = ydir.dot(p[ijk[1]]);
          double xk = xdir.dot(p[ijk[2]]);
          double yk = ydir.dot(p[ijk[2]]);
          double p = yj - yk;
          double q = xk - xj;
          double r = xj*yk - xk*yj;
          double a = 0.5*(p*xi + q*yi + r);
          double b = sqrt(p*p + q*q); 

          // Other quantities derived from the in-plane triangle data
          double xm = 0.5*(xj + xk);
          double ym = 0.5*(yj + yk);
          double binv2 = 1.0/(b*b);

          // Sum up the contributions to the s terms
          s1 += binv2*(w1*p*p + w2*q*q);
          s2 += binv2*p*q*(w1 - w2);
          s3 += binv2*(w1*p*(hbar*b - 2*a) - w2*q*((xi - xm)*q - (yi - ym)*p));
          s4 += binv2*(w1*q*q + w2*p*p);
          s5 += binv2*(w1*q*(hbar*b - 2*a) - w2*p*((yi - ym)*p - (xi - xm)*q));
        }

        // Compute the updates in the physical plane
        double det = s1*s4 - s2*s2;
        double lx = 0.0, ly = 0.0;
        if (det != 0.0){
          det = 1.0/det;
          lx = det*(s3*s4 - s2*s5);
          ly = det*(s1*s5 - s2*s3);
        }

        // Add up the displacements along the local coordinate directions
        TMRPoint dir;
        dir.x = lx*xdir.x + ly*ydir.x;
        dir.y = lx*xdir.y + ly*ydir.y;
        dir.z = lx*xdir.z + ly*ydir.z;

        // Add the parameter movement along the specified direction
        // and compute the update
        addParamMovement(1.0, &Xu, &Xv, &dir, &prm[2*i]);
        surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
      }
    }
  }
}
