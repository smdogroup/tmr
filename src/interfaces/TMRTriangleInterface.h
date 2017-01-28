#ifndef TMR_TRIANGLE_INTERFACE_H
#define TMR_TRIANGLE_INTERFACE_H

#include "TMRBase.h"

/*
  Triangle info class

  This class is used to generate a triangular mesh
*/
class TMRTriangulation {
 public:
  TMRTriangulation( int _npts, const double *_pts,
                    const int *_ptmarkers=NULL );
  ~TMRTriangulation();

  // Set the segments in the domain
  void setSegments( int _nsegments, const int *_segments,
                    const int *_segmarkers=NULL );

  // Set holes within the domain
  void setHoles( int _nholes, const double *_holes );
  
  // Triangulate the region
  void create();
  
  // Refine the triangulation - maximum area constraints
  void refine( const double areas[] ); 

  // Get the triangulation
  void getTriangulation( int *_ntris, const int **_tris );
  void getEdges( int *_nedges, const int **_edges );
  void getDualEdges( int *_nedges, const int **edgetotris );

  // Write the triangulation (if any) to VTK
  void writeToVTK( const char *filename );
  
 private:
  // Points
  int npts;       // number of points
  double *pts;    // x-y point locations, size: 2*npts
  int *ptmarkers; // point markers, size: npts
  
  // Segments
  int nsegments;   // number of segments
  int *segments;   // segment->points, size: 2*nsegments
  int *segmarkers; // segment marker values, size: nsegments

  // Holes
  int nholes;     // number of holes within the domain
  double *holes;  // a point within each hole

  // Triangles
  int ntris;         // number of triangles
  int *tris;         // triangle vertices
  int *trineighbors; // Triangle neighbors

  // Number of corners in the mesh
  int ncorners;

  // Edges in the triangle mesh
  int nedges;       // number of edges
  int *edges;       // edges->vertices, 2*nedges
  int *edgemarkers; // edge marker values

  // Dual edges (edges->triangles)
  int *dualedges;   // edges->tris, size: 2*nedges
};

#endif // TMR_TRIANGLE_INTERFACE_H
