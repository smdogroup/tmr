#ifndef TMR_MESH_H
#define TMR_MESH_H

#include "TMRBase.h"
#include "TMRGeometry.h"

/*
  The mesh for a geometric curve
*/
class TMRCurveMesh : public TMREntity {
 public:
  TMRCurveMesh( TMRCurve *curve );
  ~TMRCurveMesh();

  // Retrieve the underlying curve
  void getCurve( TMRCurve **_curve );

  // Mesh the geometric object
  void mesh( double htarget );

  // Retrieve the mesh points
  void getMesh( int *_npts, double **_pts, TMRPoint **X );

 private:
  TMRCurve *curve;

  // The parametric locations of the points that are obtained from
  // meshing the curve
  int npts;
  double *pts;
  TMRPoint *X;
};

/*
  Triangle info class

  This class is used to generate a triangular mesh
*/
class TMRSurfaceMesh : public TMREntity {
 public:
  TMRSurfaceMesh( TMRSurface *surface );
  ~TMRSurfaceMesh();

  // Retrieve the underlying surface
  void getSurface( TMRSurface **_surface );
  
  // Mesh the underlying geometric object
  void mesh( double htarget );

  // Write the quadrilateral mesh to a VTK format
  void writeToVTK( const char *filename );

  // Print the quadrilateral quality
  void printQuadQuality();

 private:
  // Print the triangle quality
  void printTriQuality( int ntris, const int tris[] );
  void writeTrisToVTK( const char *filename,
                       int ntris, const int tris[] );

  // Compute the connectivity between nodes to corresponding elements
  void computeNodeToElems( int nnodes, int nelems, int numelemnodes,
                           const int conn[], 
                           int **_ptr, int **_nodetoelems );

  // Compute the edges in a triangular or quadrilateral mesh
  void computeTriEdges( int nnodes, int ntris, const int conn[],
                        int *_num_tri_edges, int **_tri_edges,
                        int **_tri_neighbors, int **_dual_edges );
  void computeQuadEdges( int nnodes, int nquads, const int quads[],
                         int *_num_quad_edges, int **_quad_edges );

  // Recombine the mesh to a quadrilateral mesh based on the
  // quad-Blossom algorithm
  void recombine( int ntris, const int tris[], const int tri_neighbors[],
                  int num_edges, const int dual_edges[], const TMRPoint *p,
                  int *_num_quads, int **_new_quads );

  // Compute recombined triangle information
  int getRecombinedQuad( const int tris[], const int trineighbors[],
                         int t1, int t2, int quad[] );
  double computeRecombinedQuality( const int tris[], 
                                   const int trineighbors[],
                                   int t1, int t2, const TMRPoint *p ); 

  // Compute the quadrilateral quality
  double computeQuadQuality( const int *quad, const TMRPoint *p );
  double computeTriQuality( const int *tri, const TMRPoint *p );

  // The underlying surface
  TMRSurface *surface;

  // Points
  int num_fixed_pts; // number of fixed points
  int num_points; // The number of point locations
  double *pts; // The parametric node locations
  TMRPoint *X; // The physical node locations

  // Quadrilateral mesh information
  int num_quads;
  int *quads;
};

#endif // TMR_MESH_H
