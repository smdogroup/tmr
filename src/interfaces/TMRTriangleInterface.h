#ifndef TMR_TRIANGLE_INTERFACE_H
#define TMR_TRIANGLE_INTERFACE_H

#include "TMRBase.h"
#include "TMRGeometry.h"

/*
  Triangle info class

  This class is used to generate a triangular mesh
*/
class TMRTriangulation : public TMREntity {
 public:
  TMRTriangulation( int _npts, const double *_params,
                    const TMRPoint *_pts, TMRSurface *_surface );
  ~TMRTriangulation();

  // Set the segments in the domain
  void setSegments( int _nsegments, const int *_segments );

  // Set holes within the domain
  void setHoles( int _nholes, const double *_holes );
  
  // Triangulate the region - adding points in the parameter space
  void create();
  
  // Refine the triangulation - maximum area constraints
  void refine( const double areas[] ); 

  // Smooth the triangular mesh
  void laplacianSmoothing( int nsmooth );
  void springSmoothing( int nsmooth );

  // Smooth the rectangular mesh
  void laplacianQuadSmoothing( int nsmooth );
  void springQuadSmoothing( int nsmooth );

  // Remesh after smoothing operation
  void remesh();

  // Recombine the mesh to a quadrilateral mesh based on Blossom
  void recombine();

  // Get the triangulation
  void getPoints( int *_npts, const double **_params, const TMRPoint **_pts );
  void getTriangulation( int *_ntris, const int **_tris );
  void getEdges( int *_nedges, const int **_edges );
  void getDualEdges( int *_nedges, const int **edgetotris );

  // Write the triangulation (if any) to VTK
  void writeToVTK( const char *filename );
  void writeQuadToVTK( const char *filename );

  // Print the quadrilateral quality
  void printQuadQuality();
  
 private:
  // Apply the smoothing algorithm
  void laplacianSmoothing( int nsmooth,
                           int num_edges, const int *edge_list,
                           int num_pts, double *prm, TMRPoint *p );
  void springSmoothing( int nsmooth, double alpha,
                        int num_edges, const int *edge_list,
                        int num_pts, double *prm, TMRPoint *p );
  void springQuadSmoothing( int nsmooth, double alpha,
                            int num_quads, const int *quad_list,
                            int num_edges, const int *edge_list,
                            int num_pts, double *prm, TMRPoint *p );

  // Compute recombined triangle information
  int getRecombinedQuad( int t1, int t2, int quad[] );
  double computeRecombinedQuality( int t1, int t2,
                                   const TMRPoint *p ); 

  // Compute the quadrilateral quality
  double computeQuadQuality( const int *quad, const TMRPoint *p );

  // Compute the node to quad data structures
  void computeNodeToQuads( int nnodes, int **_ptr, int **_nodetoquads );

  // Compute the edges in the quad mesh
  void computeQuadEdges( int nnodes, int nquads, const int quads[],
                         int *_nquadedges, int **_quadedges );

  // The underlying surface
  TMRSurface *surface;

  // Points
  int npts;       // number of points
  double *params; // parameter locations, size: 2*npts
  TMRPoint *pts;  // physical node locations
  int *ptmarkers; // point markers, size: npts
  
  // Segments
  int nsegments;   // number of segments
  int *segments;   // segment->points, size: 2*nsegments

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

  // Dual edges (edges->triangles)
  int *dualedges;   // edges->tris, size: 2*nedges

  // Quadrilateral mesh information
  int nquads;
  int *quads;

  // Edges of the quadrilateral
  int nquadedges;
  int *quadedges;
};

#endif // TMR_TRIANGLE_INTERFACE_H
