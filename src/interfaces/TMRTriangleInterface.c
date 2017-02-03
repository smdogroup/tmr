#include "TMRTriangleInterface.h"
#include "TMRPerfectMatchInterface.h"
#include <math.h>
#include <stdio.h>

// Always use the double precision version of the code
#define REAL double
#define VOID void
#define ANSI_DECLARATORS

// Include Triangle interface
TMR_EXTERN_C_BEGIN
#include "triangle.h"
TMR_EXTERN_C_END

/*
  Free the triangleio class
*/
static void freetriangleio( struct triangulateio *tri ){
  if (tri->pointlist){ free(tri->pointlist); }
  if (tri->pointattributelist){ free(tri->pointattributelist); }
  if (tri->pointmarkerlist){ free(tri->pointmarkerlist); }
  if (tri->trianglelist){ free(tri->trianglelist); }
  if (tri->triangleattributelist){ free(tri->triangleattributelist); }
  if (tri->trianglearealist){ free(tri->trianglearealist); }
  if (tri->neighborlist){ free(tri->neighborlist); }
  if (tri->segmentlist){ free(tri->segmentlist); }
  if (tri->segmentmarkerlist){ free(tri->segmentmarkerlist); }
  if (tri->holelist){ free(tri->holelist); }
  // if (tri->regionlist){ free(tri->regionlist); }
  if (tri->edgelist){ free(tri->edgelist); }
  if (tri->edgemarkerlist){ free(tri->edgemarkerlist); }
  if (tri->normlist){ free(tri->normlist); }
}

/*
  The triangulation class
*/
TMRTriangulation::TMRTriangulation( int _npts,
                                    const double *_params,
                                    const TMRPoint *_pts,
                                    TMRSurface *_surface ){
  // Copy the surface pointer
  surface = _surface;
  surface->incref();

  // Copy the parametric point locations
  npts = _npts;
  params = (double*)malloc(2*npts*sizeof(double));
  memcpy(params, _params, 2*npts*sizeof(double));

  // Copy over the physical point locations
  pts = (TMRPoint*)malloc(npts*sizeof(TMRPoint));
  if (_pts){
    memcpy(pts, _pts, npts*sizeof(TMRPoint));
  }
  else {
    for ( int i = 0; i < npts; i++ ){
      surface->evalPoint(params[2*i], params[2*i+1], &pts[i]);
    }
  }

  // Allocate and set the markers for the points
  ptmarkers = (int*)malloc(npts*sizeof(int));
  for ( int i = 0; i < npts; i++ ){
    ptmarkers[i] = i+1;
  }

  // Set defaults for everything else
  nsegments = 0;
  segments = NULL;

  // Zero the holes
  nholes = 0;
  holes = NULL;

  // Triangles/output
  ntris = 0;
  tris = NULL;
  trineighbors = NULL;
  
  // Zero the number of corners in the mesh
  ncorners = 0;

  // Set the edges in the triangle mesh
  nedges = 0;
  edges = NULL;

  // The dual edges (from the Voronoi diagram)
  dualedges = NULL;

  // Null out the quadrilateral mesh
  nquads = 0;
  quads = NULL;

  // The quadrilateral edges
  nquadedges = 0;
  quadedges = NULL;
}

/*
  Free any allocated memory and delete the object
*/
TMRTriangulation::~TMRTriangulation(){
  if (params){ free(params); }
  if (ptmarkers){ free(ptmarkers); }
  if (segments){ free(segments); }
  if (holes){ free(holes); }
  if (tris){ free(tris); }
  if (trineighbors){ free(trineighbors); }
  if (edges){ free(edges); }
  if (dualedges){ free(dualedges); }
  if (quads){ free(quads); }
  if (quadedges){ free(quadedges); }
}

/*
  Set the segments within the domain. 

  A segment is an impenetrable line within the x-y plane
  that can be used to create boundaries/holes within the mesh.
  Markers can be set so that the edges generated along the
  segments inherit the marker value.
*/
void TMRTriangulation::setSegments( int _nsegments,
                                    const int *_segments ){
  nsegments = _nsegments;
  segments = (int*)malloc(2*nsegments*sizeof(int));
  memcpy(segments, _segments, 2*nsegments*sizeof(int));
}

/*
  Set points within a hole to indicate the presence of a hole
  within the domain
*/
void TMRTriangulation::setHoles( int _nholes,
                                 const double *_holes ){
  nholes = _nholes;
  holes = (double*)malloc(2*nholes*sizeof(double));;
  memcpy(holes, _holes, 2*nholes*sizeof(double));
}

/*
  Retrieve the points used in the triangularization
*/
void TMRTriangulation::getPoints( int *_npts, 
                                  const double **_params,
                                  const TMRPoint **_pts ){
  if (_npts){ *_npts = npts; }
  if (_params){ *_params = params; }
  if (_pts){ *_pts = pts; }
}

/*
  Get the triangulation
*/
void TMRTriangulation::getTriangulation( int *_ntris,
                                         const int **_tris ){
  if (_ntris){ *_ntris = ntris; }
  if (_tris){ *_tris = tris; }
}

/*
  Get the edges
*/
void TMRTriangulation::getEdges( int *_nedges,
                                 const int **_edges ){
  if (_nedges){ *_nedges = nedges; }
  if (_edges){ *_edges = edges; }
}

/*
  Get the dual of the edges
*/
void TMRTriangulation::getDualEdges( int *_nedges,
                                     const int **edgetotris ){
  if (_nedges){ *_nedges = nedges; }
  if (edgetotris){ *edgetotris = dualedges; }
}

/*
  Triangulate the points within the mesh

  This call does not modify the points/markers
*/
void TMRTriangulation::create(){
  // Set up the Triangle data for input/output
  struct triangulateio in, out, vorout;
  memset(&in, 0, sizeof(in));
  memset(&out, 0, sizeof(out));
  memset(&vorout, 0, sizeof(vorout));

  // Set the points for input
  in.numberofpoints = npts;
  in.pointlist = params;
  in.pointmarkerlist = ptmarkers;
    
  // The number of segments
  in.numberofsegments = nsegments;
  in.segmentlist = segments;
  
  // Copy over the holes if any
  in.numberofholes = nholes;
  in.holelist = holes;

  // Set the regions
  in.numberofregions = 1;
  in.regionlist = (double*)malloc(4*sizeof(double));
  in.regionlist[0] = 0.0;
  in.regionlist[1] = 0.0;
  in.regionlist[2] = 1.0; // Regional attribute for the whole mesh
  in.regionlist[3] = 1.0; // Area constraint - unused

  // Perform an initial trianguarlization of the points
  char opts[] = "pczevnYY";
  triangulate(opts, &in, &out, &vorout);

  // Copy out the point information
  npts = out.numberofpoints;
  params = out.pointlist;           out.pointlist = NULL;
  ptmarkers = out.pointmarkerlist;  out.pointmarkerlist = NULL;

  // Copy out the segment information
  nsegments = out.numberofsegments;
  segments = out.segmentlist;       out.segmentlist = NULL;
  
  // Copy out the triangulation information
  ntris = out.numberoftriangles;
  ncorners = out.numberofcorners;
  tris = out.trianglelist;          out.trianglelist = NULL;
  trineighbors = out.neighborlist;  out.neighborlist = NULL;

  // Copy out the edge and edge marker information
  nedges = out.numberofedges;  
  edges = out.edgelist;             out.edgelist = NULL;

  // Record the vorout data
  dualedges = vorout.edgelist;      vorout.edgelist = NULL;

  // Free all of the remaining data
  freetriangleio(&in);
  freetriangleio(&out);
  freetriangleio(&vorout);

  // Allocate a new array of points for the new connectivity
  free(pts);
  pts = (TMRPoint*)malloc(npts*sizeof(TMRPoint));
  for ( int i = 0; i < npts; i++ ){
    surface->evalPoint(params[2*i], params[2*i+1], &pts[i]);
  }
}

/*
  Refine the mesh
*/
void TMRTriangulation::refine( const double areas[] ){
  // Set up the Triangle data for input/output
  struct triangulateio in, out, vorout;
  memset(&in, 0, sizeof(in));
  memset(&out, 0, sizeof(out));
  memset(&vorout, 0, sizeof(vorout));

  // Set the points for input
  in.numberofpoints = npts;
  in.pointlist = params;
  in.pointmarkerlist = ptmarkers;

  // Set the triangles
  in.numberoftriangles = ntris;
  in.numberofcorners = ncorners;
  in.neighborlist = trineighbors;
  in.trianglelist = tris;

  // Set the areas of the triangles as a constraint
  in.trianglearealist = (double*)malloc(ntris*sizeof(double));
  memcpy(in.trianglearealist, areas, ntris*sizeof(double));

  // The number of segments
  in.numberofsegments = nsegments;
  in.segmentlist = segments;
  
  // Copy over the holes if any
  in.numberofholes = nholes;
  in.holelist = holes;

  // Set the regions
  in.numberofregions = 1;
  in.regionlist = (double*)malloc(4*sizeof(double));
  in.regionlist[0] = 0.0;
  in.regionlist[1] = 0.0;
  in.regionlist[2] = 1.0; // Regional attribute for the whole mesh
  in.regionlist[3] = 1.0; // Area constraint - unused

  // Set the edges
  in.numberofedges = nedges;
  in.edgelist = edges;
  
  // Perform an initial trianguarlization of the points
  char opts[] = "prazenvYY";
  triangulate(opts, &in, &out, &vorout);

  // Copy out the point information
  npts = out.numberofpoints;
  params = out.pointlist;           out.pointlist = NULL;
  ptmarkers = out.pointmarkerlist;  out.pointmarkerlist = NULL;

  // Copy out the segment information
  nsegments = out.numberofsegments;
  segments = out.segmentlist;       out.segmentlist = NULL;
  
  // Copy out the triangulation information
  ntris = out.numberoftriangles;
  ncorners = out.numberofcorners;
  tris = out.trianglelist;          out.trianglelist = NULL;
  trineighbors = out.neighborlist;  out.neighborlist = NULL;

  // Copy out the edge and edge marker information
  nedges = out.numberofedges;
  edges = out.edgelist;             out.edgelist = NULL;

  // Record the vorout data
  dualedges = vorout.edgelist;      vorout.edgelist = NULL;

  // Free all of the remaining data
  freetriangleio(&in);
  freetriangleio(&out);
  freetriangleio(&vorout);

  // Allocate a new array of points for the new connectivity
  free(pts);
  pts = (TMRPoint*)malloc(npts*sizeof(TMRPoint));
  for ( int i = 0; i < npts; i++ ){
    surface->evalPoint(params[2*i], params[2*i+1], &pts[i]);
  }
}

/*
  Add the motion of the point in parameter space
*/
inline void addParamMovement( double alpha,
                              TMRPoint *Xu, TMRPoint *Xv,
                              TMRPoint *d, double *delta ){
  // Compute the sum - in the parameter space - of the motion 
  // contributed by each node
  double a11 = Xu->dot(Xu);
  double a12 = Xu->dot(Xv);
  double a22 = Xv->dot(Xv);

  // Compute the right hand sides
  double b1 = d->dot(Xu);
  double b2 = d->dot(Xv);

  // Compute the inverse of the determinant
  double invdet = alpha/(a11*a22 - a12*a12);

  // Add the contribution
  delta[0] += invdet*(a22*b1 - a12*b2);
  delta[1] += invdet*(a11*b2 - a12*b1);
}

/*
  Apply Laplacian smoothing
*/
void TMRTriangulation::laplacianSmoothing( int nsmooth,
                                           int num_edges,
                                           const int *edge_list,
                                           int num_pts,
                                           double *prm,
                                           TMRPoint *p ){
  int *count = (int*)malloc(num_pts*sizeof(int));
  double *new_params = (double*)malloc(2*num_pts*sizeof(double));
  TMRPoint *Xu = (TMRPoint*)malloc(num_pts*sizeof(TMRPoint));
  TMRPoint *Xv = (TMRPoint*)malloc(num_pts*sizeof(TMRPoint));

  for ( int iter = 0; iter < nsmooth; iter++ ){
    memset(count, 0, num_pts*sizeof(int));
    memset(new_params, 0, 2*num_pts*sizeof(double));

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = 0; i < num_pts; i++ ){
      if (ptmarkers[i] == 0){
        surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
      }
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
      if (ptmarkers[n1] == 0){
        addParamMovement(1.0, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
        count[n1]++;
      }

      // Add the movement of the second node in parameter space
      if (ptmarkers[n2] == 0){
        addParamMovement(-1.0, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
        count[n2]++;
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = 0; i < npts; i++ ){
      if (ptmarkers[i] == 0 && count[i] > 0){
        prm[2*i] += new_params[2*i]/count[i];
        prm[2*i+1] += new_params[2*i+1]/count[i];
        surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
      }
    }
  }

  free(count);
  free(new_params);
  free(Xu);
  free(Xv);
}

/*
  Apply the spring smoothing analogy
*/
void TMRTriangulation::springSmoothing( int nsmooth,
                                        double alpha,
                                        int num_edges,
                                        const int *edge_list,
                                        int num_pts,
                                        double *prm,
                                        TMRPoint *p ){

  double *len = (double*)malloc(num_edges*sizeof(double));
  double *new_params = (double*)malloc(2*num_pts*sizeof(double));
  TMRPoint *Xu = (TMRPoint*)malloc(num_pts*sizeof(TMRPoint));
  TMRPoint *Xv = (TMRPoint*)malloc(num_pts*sizeof(TMRPoint));

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
    double len0 = 0.9*sum/nedges;

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = 0; i < num_pts; i++ ){
      if (ptmarkers[i] == 0){
        surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
      }
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
      if (ptmarkers[n1] == 0){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }

      // Add the movement of the second node in parameter space
      if (ptmarkers[n2] == 0){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = 0; i < npts; i++ ){
      if (ptmarkers[i] == 0){
        prm[2*i] += alpha*new_params[2*i];
        prm[2*i+1] += alpha*new_params[2*i+1];
        surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
      }
    }
  }

  free(new_params);
  free(len);
  free(Xu);
  free(Xv);
}

/*
  Apply smoothing to the triangular mesh
*/
void TMRTriangulation::laplacianSmoothing( int nsmooth ){
  laplacianSmoothing(nsmooth, nedges, edges, npts, params, pts);
}

/*
  Apply the spring smoothing to the triangular mesh
*/
void TMRTriangulation::springSmoothing( int nsmooth ){
  double alpha = 0.1;
  springSmoothing(nsmooth, alpha, nedges, edges, npts, params, pts);
}

/*
  Apply smoothing to the triangular mesh
*/
void TMRTriangulation::laplacianQuadSmoothing( int nsmooth ){
  laplacianSmoothing(nsmooth, nquadedges, quadedges, npts, params, pts);
}

/*
  Apply the spring smoothing to the triangular mesh
*/
void TMRTriangulation::springQuadSmoothing( int nsmooth ){
  double alpha = 0.1;
  springSmoothing(nsmooth, alpha, nquadedges, quadedges, npts, params, pts);
}

/*
  Smooth the mesh using laplacian smoothing and then remesh
*/
void TMRTriangulation::remesh(){
  // Set up the Triangle data for input/output
  struct triangulateio in, out, vorout;
  memset(&in, 0, sizeof(in));
  memset(&out, 0, sizeof(out));
  memset(&vorout, 0, sizeof(vorout));

  // Set the points for input
  in.numberofpoints = npts;
  in.pointlist = params;
  in.pointmarkerlist = ptmarkers;

  // The number of segments
  in.numberofsegments = nsegments;
  in.segmentlist = segments;

  // Set the triangles
  in.numberoftriangles = ntris;
  in.numberofcorners = ncorners;
  in.neighborlist = trineighbors;
  in.trianglelist = tris;
  
  // Copy over the holes if any
  in.numberofholes = nholes;
  in.holelist = holes;

  // Set the regions
  in.numberofregions = 1;
  in.regionlist = (double*)malloc(4*sizeof(double));
  in.regionlist[0] = 0.0;
  in.regionlist[1] = 0.0;
  in.regionlist[2] = 1.0; // Regional attribute for the whole mesh
  in.regionlist[3] = 1.0; // Area constraint - unused

  // Set the edges
  in.numberofedges = nedges;
  in.edgelist = edges;

  // Perform an initial trianguarlization of the points
  char opts[] = "zenv";
  triangulate(opts, &in, &out, &vorout);

  // Copy out the point information
  npts = out.numberofpoints;
  params = out.pointlist;           out.pointlist = NULL;
  ptmarkers = out.pointmarkerlist;  out.pointmarkerlist = NULL;

  // Copy out the segment information
  nsegments = out.numberofsegments;
  segments = out.segmentlist;       out.segmentlist = NULL;
  
  // Copy out the triangulation information
  ntris = out.numberoftriangles;
  ncorners = out.numberofcorners;
  tris = out.trianglelist;          out.trianglelist = NULL;
  trineighbors = out.neighborlist;  out.neighborlist = NULL;

  // Copy out the edge and edge marker information
  nedges = out.numberofedges;
  edges = out.edgelist;             out.edgelist = NULL;

  // Record the vorout data
  dualedges = vorout.edgelist;      vorout.edgelist = NULL;

  // Free all of the remaining data
  freetriangleio(&in);
  freetriangleio(&out);
  freetriangleio(&vorout);

  // Allocate a new array of points for the new connectivity
  free(pts);
  pts = (TMRPoint*)malloc(npts*sizeof(TMRPoint));
  for ( int i = 0; i < npts; i++ ){
    surface->evalPoint(params[2*i], params[2*i+1], &pts[i]);
  }
}

/*
  Get the quadrilateral elemnet obtained by combining the triangles
  t1 and t2 together

        . --- .
      / |   /
    /   | /
  . --- . 
*/
int TMRTriangulation::getRecombinedQuad( int t1, int t2,
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
  Compute a node to quad data structure
*/
void TMRTriangulation::computeNodeToQuads( int nnodes, int **_ptr, 
                                           int **_nodetoquads ){
  // Set the pointer
  int *ptr = (int*)malloc((nnodes+1)*sizeof(int));
  memset(ptr, 0, (nnodes+1)*sizeof(int));

  // Count up the references
  for ( int i = 0; i < 4*nquads; i++ ){
    ptr[quads[i]+1]++;
  }

  // Set the pointer into the quad array
  for ( int i = 0; i < nnodes; i++ ){
    ptr[i+1] += ptr[i];
  }

  // Compute the node to quads
  int *nodetoquads = (int*)malloc(ptr[nnodes]*sizeof(int));
  for ( int i = 0; i < nquads; i++ ){
    for ( int j = 0; j < 4; j++ ){
      int node = quads[4*i+j];
      nodetoquads[ptr[node]] = i;
      ptr[node]++;
    }
  }

  // Reset the pointer array
  for ( int i = nnodes-1; i >= 0; i-- ){
    ptr[i+1] = ptr[i];
  }
  ptr[0] = 0;

  // Set the output points
  *_ptr = ptr;
  *_nodetoquads = nodetoquads;
}

/*
  Get the neighbors of the recombined quad mesh
*/
void TMRTriangulation::computeQuadEdges( int nnodes, int nquads, 
                                         const int quads[],
                                         int *_nquadedges,
                                         int **_quadedges ){
  int *ptr, *nodetoquads;
  computeNodeToQuads(nnodes, &ptr, &nodetoquads);

  // Now compute the neighbors for each quad
  int *quadedgenums = (int*)malloc(4*nquads*sizeof(int));
  for ( int i = 0; i < 4*nquads; i++ ){
    quadedgenums[i] = -1;
  }

  // Quck reference from the quad index to the edge
  const int enodes[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

  int count = 0;
  int ne = 0;
  for ( int i = 0; i < nquads; i++ ){
    // Search through each edge of the each quadrilateral
    for ( int j = 0; j < 4; j++ ){
      if (quadedgenums[4*i+j] < 0){
        quadedgenums[4*i+j] = ne;

        // Search for the neighboring quad that shares this edge
        int kp = ptr[quads[4*i + enodes[j][0]]];
        int kpend = ptr[quads[4*i + enodes[j][0]]+1];
        for ( ; kp < kpend; kp++ ){
          // Find the potential quad neighbor
          int n = nodetoquads[kp];

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
              quadedgenums[4*n+e] = ne;
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
  free(ptr);
  free(nodetoquads);

  // Now we have a unique list of edge numbers and the total number of
  // edges, we can construct the unique edge list
  int *quadedges = (int*)malloc(2*ne*sizeof(int));
  for ( int i = 0; i < nquads; i++ ){
    for ( int j = 0; j < 4; j++ ){
      // Get the unique edge number
      int n = quadedgenums[4*i+j];
      quadedges[2*n] = quads[4*i+enodes[j][0]];
      quadedges[2*n+1] = quads[4*i+enodes[j][1]];
    }
  } 

  // Free the data
  free(quadedgenums);

  *_nquadedges = ne;
  *_quadedges = quadedges;
}

/*
  Compute the recombined quality
*/
double TMRTriangulation::computeRecombinedQuality( int t1, int t2,
                                                   const TMRPoint *p ){
  // Find the combined quadrilateral from the two given triangles
  int quad[4];
  int fail = getRecombinedQuad(t1, t2, quad);
  if (fail){
    return 0.0;
  }

  return computeQuadQuality(quad, p);
}

/*
  Compute the quality of a quadrilateral element
*/
double TMRTriangulation::computeQuadQuality( const int *quad,
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
  Recombine the triangulation into a quadrilateral mesh
*/
void TMRTriangulation::recombine(){
  // Count up the number of dual edges
  int ndualedges = 0;
  for ( int i = 0; i < nedges; i++ ){
    if (dualedges[2*i] >= 0 && dualedges[2*i+1] >= 0){
      ndualedges++;
    }
  }

  // Allocate the weights
  double *weights = (double*)malloc(ndualedges*sizeof(double));
  int *graphedges = (int*)malloc(2*ndualedges*sizeof(double));

  // Compute the weight associated with each edge by combputing the
  // recombined quality
  for ( int i = 0, e = 0; i < nedges; i++ ){
    int t1 = dualedges[2*i];
    int t2 = dualedges[2*i+1];

    if (t1 >= 0 && t2 >= 0){
      // Compute the weight for this recombination
      double quality = computeRecombinedQuality(t1, t2, pts);
      double weight = 1000.0*(1.0 - quality)*(1.0 + 1.0/(quality + 0.1));
      graphedges[2*e] = t1;
      graphedges[2*e+1] = t2;
      weights[e] = weight;
      e++;
    }
  }

  int *match = (int*)malloc(2*ntris*sizeof(int));
  TMR_PerfectMatchGraph(ntris, ndualedges, graphedges, weights, match);

  free(weights);
  free(graphedges);

  // Get the required options
  nquads = ntris/2;
  quads = (int*)malloc(4*nquads*sizeof(int));

  // Recombine the quads
  for ( int i = 0; i < nquads; i++ ){
    getRecombinedQuad(match[2*i], match[2*i+1], &quads[4*i]);
  }
  free(match);

  // Compute the quad edges 
  computeQuadEdges(npts, nquads, quads, &nquadedges, &quadedges);
}

/*
  Write the quadrilateral mesh to a VTK file
*/
void TMRTriangulation::writeQuadToVTK( const char *filename ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts);
    for ( int k = 0; k < npts; k++ ){
      fprintf(fp, "%e %e %e\n", pts[k].x, pts[k].y, pts[k].z);
    }
    
    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", nquads, 5*nquads);
    for ( int k = 0; k < nquads; k++ ){
      fprintf(fp, "4 %d %d %d %d\n", quads[4*k], quads[4*k+1], 
              quads[4*k+2], quads[4*k+3]);
    }

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", nquads);
    for ( int k = 0; k < nquads; k++ ){
      fprintf(fp, "%d\n", 9);
    }

    // Print out the rest as fields one-by-one
    fprintf(fp, "CELL_DATA %d\n", nquads);
    fprintf(fp, "SCALARS quality float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for ( int k = 0; k < nquads; k++ ){
      fprintf(fp, "%e\n", computeQuadQuality(&quads[4*k], pts));
    }

    fclose(fp);
  }
}

/*
  Print the quad quality
*/
void TMRTriangulation::printQuadQuality(){
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins*sizeof(int));
  for ( int i = 0; i < nquads; i++ ){
    double quality = computeQuadQuality(&quads[4*i], pts);

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
  Write the output to a VTK file
*/
void TMRTriangulation::writeToVTK( const char *filename ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts);
    for ( int k = 0; k < npts; k++ ){
      fprintf(fp, "%e %e %e\n", pts[k].x, pts[k].y, pts[k].z);
    }
    
    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", ntris, 4*ntris);
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "3 ");
      for ( int j = 0; j < 3; j++ ){
        int node = tris[3*k+j];
        fprintf(fp, "%d ", node);
      }
      fprintf(fp, "\n");
    }

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", ntris);
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "%d\n", 5);
    }

    fclose(fp);
  }
}
