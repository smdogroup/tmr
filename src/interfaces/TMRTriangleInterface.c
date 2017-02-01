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
  The triangulation class
*/
TMRTriangulation::TMRTriangulation( int _npts,
                                    const double *_pts,
                                    const int *_ptmarkers ){
  npts = _npts;
  pts = (double*)malloc(2*npts*sizeof(double));
  memcpy(pts, _pts, 2*npts*sizeof(double));

  // Allocate and set the markers for the points
  ptmarkers = (int*)malloc(npts*sizeof(int));
  if (_ptmarkers){
    memcpy(ptmarkers, _ptmarkers, npts*sizeof(int));
  }
  else {
    for ( int i = 0; i < npts; i++ ){
      ptmarkers[i] = 1;
    }
  }

  // Set defaults for everything else
  nsegments = 0;
  segments = NULL;
  segmarkers = NULL;

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
  edgemarkers = NULL;
  
  // The dual edges (from the Voronoi diagram)
  dualedges = NULL;
}

/*
  Free any allocated memory and delete the object
*/
TMRTriangulation::~TMRTriangulation(){
  if (pts){ free(pts); }
  if (ptmarkers){ free(ptmarkers); }
  if (segments){ free(segments); }
  if (segmarkers){ free(segmarkers); }
  if (holes){ free(holes); }
  if (tris){ free(tris); }
  if (trineighbors){ free(trineighbors); }
  if (edges){ free(edges); }
  if (edgemarkers){ free(edgemarkers); }
  if (dualedges){ free(dualedges); }
}

/*
  Set the segments within the domain. 

  A segment is an impenetrable line within the x-y plane
  that can be used to create boundaries/holes within the mesh.
  Markers can be set so that the edges generated along the
  segments inherit the marker value.
*/
void TMRTriangulation::setSegments( int _nsegments,
                                    const int *_segments,
                                    const int *_segmarkers ){
  nsegments = _nsegments;
  segments = (int*)malloc(2*nsegments*sizeof(int));
  memcpy(segments, _segments, 2*nsegments*sizeof(int));

  if (_segmarkers){
    segmarkers = (int*)malloc(nsegments*sizeof(int));
    memcpy(segmarkers, _segmarkers, nsegments*sizeof(int));
  }
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
void TMRTriangulation::getPoints( int *_npts, const double **_pts ){
  if (_npts){ *_npts = npts; }
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

  // Set the points for input
  in.numberofpoints = npts;
  in.numberofpointattributes = 0;
  in.pointlist = pts;
  in.pointattributelist = NULL;
  in.pointmarkerlist = ptmarkers;
    
  // The number of segments
  in.numberofsegments = nsegments;
  in.segmentlist = segments;
  in.segmentmarkerlist = segmarkers;
  
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

  // Set output data that will be created
  out.pointlist = NULL;
  out.pointattributelist = NULL;
  out.pointmarkerlist = NULL;

  // Triangle information
  out.trianglelist = NULL;
  out.triangleattributelist = NULL;
  out.neighborlist = NULL;

  // Edge information
  out.edgelist = NULL;
  out.edgemarkerlist = NULL;
  out.segmentlist = NULL;
  out.segmentmarkerlist = NULL;

  // NULL out things for the voronoi output
  vorout.pointlist = NULL;
  vorout.pointattributelist = NULL;
  vorout.edgelist = NULL;
  vorout.normlist = NULL;

  // Perform an initial trianguarlization of the points
  char opts[] = "pczevn";
  triangulate(opts, &in, &out, &vorout);

  // Copy out things that have been created
  // Copy out the point information
  npts = out.numberofpoints;
  pts = out.pointlist;
  ptmarkers = out.pointmarkerlist;

  // Copy out the segment information
  nsegments = out.numberofsegments;
  segments = out.segmentlist;
  segmarkers = out.segmentmarkerlist;

  // Copy out the triangulation information
  ntris = out.numberoftriangles;
  ncorners = out.numberofcorners;
  tris = out.trianglelist;
  trineighbors = out.neighborlist;

  // Copy out the edge and edge marker information
  edges = out.edgelist;
  edgemarkers = out.edgemarkerlist;

  // Record the vorout data
  dualedges = vorout.edgelist;

  // Free all of the remaining data
  free(in.pointlist);
  free(in.pointattributelist);
  free(in.pointmarkerlist);
  free(in.regionlist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  free(vorout.pointlist);
  free(vorout.pointattributelist);
  free(vorout.normlist);
}

/*
  Refine the mesh
*/
void TMRTriangulation::refine( const double areas[] ){
  // Set up the Triangle data for input/output
  struct triangulateio in, out, vorout;

    // Set the points for input
  in.numberofpoints = npts;
  in.numberofpointattributes = 0;
  in.pointlist = pts;
  in.pointattributelist = NULL;
  in.pointmarkerlist = ptmarkers;

  // Set the triangles
  in.numberoftriangles = ntris;
  in.numberofcorners = ncorners;
  in.numberoftriangleattributes = 0;
  in.triangleattributelist = NULL;
  in.neighborlist = trineighbors;
  in.trianglelist = tris;

  // Set the areas of the triangles as a constraint
  in.trianglearealist = (double*)malloc(ntris*sizeof(double));
  memcpy(in.trianglearealist, areas, ntris*sizeof(double));

  // The number of segments
  in.numberofsegments = nsegments;
  in.segmentlist = segments;
  in.segmentmarkerlist = segmarkers;
  
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
  in.edgemarkerlist = edgemarkers;
  
  // Set output data that will be created
  out.pointlist = NULL;
  out.pointattributelist = NULL;
  out.pointmarkerlist = NULL;
  out.trianglelist = NULL;
  out.triangleattributelist = NULL;
  out.edgelist = NULL;
  out.neighborlist = NULL;
  out.edgemarkerlist = NULL;
  out.segmentlist = NULL;
  out.segmentmarkerlist = NULL;

  // NULL out things for the voronoi output
  vorout.pointlist = NULL;
  vorout.pointattributelist = NULL;
  vorout.edgelist = NULL;
  vorout.normlist = NULL;

  // Perform an initial trianguarlization of the points
  char opts[] = "prazenvYY";
  triangulate(opts, &in, &out, &vorout);

  // Copy out things that have been created
  npts = out.numberofpoints;
  pts = out.pointlist;
  nsegments = out.numberofsegments;
  segments = out.segmentlist;
  segmarkers = out.segmentmarkerlist;
  ntris = out.numberoftriangles;
  tris = out.trianglelist;
  trineighbors = out.neighborlist;
  nedges = out.numberofedges;
  edges = out.edgelist;
  edgemarkers = out.edgemarkerlist;

  // Record the vorout data
  dualedges = vorout.edgelist;

  // Free all of the remaining data
  free(in.trianglelist);
  free(in.trianglearealist);
  free(in.pointlist);
  free(in.pointattributelist);
  free(in.pointmarkerlist);
  free(in.regionlist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  free(vorout.pointlist);
  free(vorout.pointattributelist);
  free(vorout.normlist);
}

/*
  Smooth the mesh using laplacian smoothing and then remesh
*/
void TMRTriangulation::smooth( int nsmooth ){
  // Allocate a new set of points
  int *count = (int*)malloc(npts*sizeof(int));
  double *spts = (double*)malloc(2*npts*sizeof(double));

  for ( int iter = 0; iter < nsmooth; iter++ ){
    memset(count, 0, npts*sizeof(int));
    memset(spts, 0, 2*npts*sizeof(double));

    // Loop over all the points and 
    for ( int i = 0; i < nedges; i++ ){
      int n1 = edges[2*i];
      int n2 = edges[2*i+1];
      spts[2*n2] += pts[2*n1];
      spts[2*n2+1] += pts[2*n1+1];
      spts[2*n1] += pts[2*n2];
      spts[2*n1+1] += pts[2*n2+1];
      count[n1]++;
      count[n2]++;
    }

    // Set the locations for the new points, keep in place
    for ( int i = 0; i < npts; i++ ){
      if (ptmarkers[i] == 0 && count[i] > 0){
        pts[2*i] = spts[2*i]/count[i];
        pts[2*i+1] = spts[2*i+1]/count[i];
      }
    }
  }

  free(spts);
  free(count);

  // Set up the Triangle data for input/output
  struct triangulateio in, out, vorout;

    // Set the points for input
  in.numberofpoints = npts;
  in.numberofpointattributes = 0;
  in.pointlist = pts;
  in.pointattributelist = NULL;
  in.pointmarkerlist = ptmarkers;

  // Set the triangles
  in.numberoftriangles = ntris;
  in.numberofcorners = ncorners;
  in.numberoftriangleattributes = 0;
  in.triangleattributelist = NULL;
  in.neighborlist = trineighbors;
  in.trianglelist = tris;
  in.trianglearealist = NULL;

  // The number of segments
  in.numberofsegments = nsegments;
  in.segmentlist = segments;
  in.segmentmarkerlist = segmarkers;
  
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
  in.edgemarkerlist = edgemarkers;
  
  // Set output data that will be created
  out.pointlist = NULL;
  out.pointattributelist = NULL;
  out.pointmarkerlist = NULL;
  out.trianglelist = NULL;
  out.triangleattributelist = NULL;
  out.edgelist = NULL;
  out.neighborlist = NULL;
  out.edgemarkerlist = NULL;
  out.segmentlist = NULL;
  out.segmentmarkerlist = NULL;

  // NULL out things for the voronoi output
  vorout.pointlist = NULL;
  vorout.pointattributelist = NULL;
  vorout.edgelist = NULL;
  vorout.normlist = NULL;

  // Perform an initial trianguarlization of the points
  char opts[] = "pzenvYY";
  triangulate(opts, &in, &out, &vorout);

  // Copy out things that have been created
  npts = out.numberofpoints;
  pts = out.pointlist;
  nsegments = out.numberofsegments;
  segments = out.segmentlist;
  segmarkers = out.segmentmarkerlist;
  ntris = out.numberoftriangles;
  tris = out.trianglelist;
  trineighbors = out.neighborlist;
  nedges = out.numberofedges;
  edges = out.edgelist;
  edgemarkers = out.edgemarkerlist;

  // Record the vorout data
  dualedges = vorout.edgelist;

  // Free all of the remaining data
  free(in.trianglelist);
  free(in.trianglearealist);
  free(in.pointlist);
  free(in.pointattributelist);
  free(in.pointmarkerlist);
  free(in.regionlist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  free(vorout.pointlist);
  free(vorout.pointattributelist);
  free(vorout.normlist);
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
  if (eta < 0.0){ eta = 0.0; }

  return eta;
}

/*
  Recombine the triangulation into a quadrilateral mesh
*/
void TMRTriangulation::recombine( const TMRPoint *p ){
  // Count up the number of dual edges
  int ndualedges = 0;
  for ( int i = 0; i < nedges; i++ ){
    if (dualedges[2*i] >= 0 && dualedges[2*i+1] >= 0){
      ndualedges++;
    }
  }

  // Allocate the weights
  double *weights = new double[ ndualedges ];
  int *graphedges = new int[ 2*ndualedges ];

  // Compute the weight associated with each edge by combputing the
  // recombined quality
  for ( int i = 0, e = 0; i < nedges; i++ ){
    int t1 = dualedges[2*i];
    int t2 = dualedges[2*i+1];

    if (t1 >= 0 && t2 >= 0){
      // Compute the weight for this recombination
      double quality = computeRecombinedQuality(t1, t2, p);
      double weight = 1000.0*(1.0 - quality)*(1.0 + 0.1/(quality + 0.1));
      graphedges[2*e] = t1;
      graphedges[2*e+1] = t2;
      weights[e] = weight;
      e++;
    }
  }

  int *match;
  TMR_PerfectMatchGraph(ntris, ndualedges, graphedges, weights, &match);

  delete [] weights;
  delete [] graphedges;

  // Get the required options
  int nquads = ntris/2;
  int *quads = new int[ 4*nquads ];

  for ( int i = 0; i < nquads; i++ ){
    getRecombinedQuad(match[2*i], match[2*i+1], &quads[4*i]);
  }
  delete [] match;

  FILE *fp = fopen("match.vtk", "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts);
    for ( int k = 0; k < npts; k++ ){
      fprintf(fp, "%e %e %e\n", pts[2*k], pts[2*k+1], 0.0);
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

    fclose(fp);
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
      fprintf(fp, "%e %e %e\n", pts[2*k], pts[2*k+1], 0.0);
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
