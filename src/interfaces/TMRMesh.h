#ifndef TMR_MESH_H
#define TMR_MESH_H

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

 private:
  TMRCurve *curve;

  // The parametric locations of the points that are obtained from
  // meshing the curve
  int npts;
  double *pts;
};

/*
  Triangle info class

  This class is used to generate a triangular mesh
*/
class TMRSurfaceMesh : public TMREntity {
 public:
  TMRSurfaceMesh( int ncurves, TMRCurveMesh **curve_mesh,
                  TMRSurface *surface );
  ~TMRSurfaceMesh();

  // Retrieve the underlying surface
  void getSurface( TMRSurface **_surface );
  
  // Mesh the underlying geometric object
  void mesh( double htarget );

  // Write the triangulation (if any) to VTK
  void writeToVTK( const char *filename );
  void writeQuadToVTK( const char *filename );

  // Print the quadrilateral quality
  void printQuadQuality();

 private:
  // Compute the connectivity between nodes to corresponding elements
  void computeNodeToElems( int nnodes, int nelems, int numelemnodes,
                           const int conn[], 
                           int **_ptr, int **_nodetoelems );

  // Compute the edges in a triangular or quadrilateral mesh
  void computeTriEdges( int nnodes, int ntris, const int conn[],
                        int *_numtriedges, int **_triedges,
                        int **_trineighbors, int **_dualedges );
  void computeQuadEdges( int nnodes, int nquads, const int quads[],
                         int *_nquadedges, int **_quadedges );

  // Recombine the mesh to a quadrilateral mesh based on the
  // quad-Blossom algorithm
  void recombine();

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
  int getRecombinedQuad( const int tris[], const int trineighbors[],
                         int t1, int t2, int quad[] );
  double computeRecombinedQuality( int t1, int t2,
                                   const TMRPoint *p ); 

  // Compute the quadrilateral quality
  double computeQuadQuality( const int *quad, const TMRPoint *p );
  double computeTriQuality( const int *tri, const TMRPoint *p );

  // The underlying surface
  TMRSurface *surface;

  // The 1d meshes that surround the surface
  int num_curves;
  TMRCurveMesh *curve_meshes;

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

#endif // TMR_MESH_H
