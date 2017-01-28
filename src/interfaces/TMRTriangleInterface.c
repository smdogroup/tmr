#include "TMRTriangleInterface.h"
#include <stdio.h>

// Always use the double precision version of the code
#define REAL double
#define VOID void
#define ANSI_DECLARATORS

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
  ptmarkers = NULL;
  if (_ptmarkers){
    ptmarkers = (int*)malloc(npts*sizeof(int));
    memcpy(ptmarkers, _ptmarkers, npts*sizeof(int));
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
  out.trianglelist = NULL;
  out.triangleattributelist = NULL;
  out.neighborlist = NULL;
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
  npts = out.numberofpoints;
  pts = out.pointlist;
  nsegments = out.numberofsegments;
  segments = out.segmentlist;
  segmarkers = out.segmentmarkerlist;
  ntris = out.numberoftriangles;
  ncorners = out.numberofcorners;
  tris = out.trianglelist;
  trineighbors = out.neighborlist;
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
  // free(out.triangleattributelist);
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
  out.edgemarkerlist = NULL;
  out.segmentlist = NULL;
  out.segmentmarkerlist = NULL;

  // NULL out things for the voronoi output
  vorout.pointlist = NULL;
  vorout.pointattributelist = NULL;
  vorout.edgelist = NULL;
  vorout.normlist = NULL;

  // Perform an initial trianguarlization of the points
  char opts[] = "prazevYY";
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
  Write the output to a VTK file
*/
void TMRTriangulation::writeToVTK( const char *filename ){
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "vtk output\nASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    
  // Write out the points
  fprintf(fp, "POINTS %d float\n", npts);
  for ( int k = 0; k < npts; k++ ){
    fprintf(fp, "%e %e %e\n", pts[2*k], pts[2*k+1], 0.0);
  }
  
  // Write out the cell values
  fprintf(fp, "\nCELLS %d %d\n",ntris, 4*ntris);
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
