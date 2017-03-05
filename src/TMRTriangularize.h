#ifndef TMR_TRIANGULARIZE_H
#define TMR_TRIANGULARIZE_H

#include "TMRBase.h"
#include "TMRTopology.h"

/*
  The rectangular domain used to define the upper/lower limits of the
  quadtree. Points cannot be reliably added outside the domain without
  bad stuff happening.
*/
class TMRQuadDomain {
 public:
  double xlow, xhigh;
  double ylow, yhigh; 
};

/*
  The following class implments a simple quadtree data structure for fast
  geometric searching queries.

  This uses a recursive implementation for searching, adding and deleting
  nodes from the quadtree. The leafs of the quadtree are of fixed size
  so that it is not necessary to 
*/
class TMRQuadNode {
 public:
  static const int MAX_DEPTH = 30;
  static const int NODES_PER_LEVEL = 10;

  // Create the root in a quadtree
  TMRQuadNode( TMRQuadDomain *_domain );

  // Create the root node (or other nodes) 
  ~TMRQuadNode();
 
  // Add/delete nodes from the quadtree
  // ----------------------------------
  void addNode( uint32_t num, const double pt[] );
  int deleteNode( uint32_t num, const double pt[] );

  // Find the closest indexed point to the provided (x,y) location
  // -------------------------------------------------------------
  uint32_t findClosest( const double pt[], double *_dist );

 private:
  // This is only for creating children
  TMRQuadNode( TMRQuadDomain *_domain,
               uint32_t _u, uint32_t _v, int _level );

  // Initialize things
  void initialize( TMRQuadDomain *_domain,
                   uint32_t _u, uint32_t _v, int _level );

  // The recursive call to the quadtree data structure
  void findClosest( const double pt[], uint32_t *index, double *dist );

  // The 2D domain
  TMRQuadDomain *domain;

  // Keep track of the children - if any
  TMRQuadNode *low_left, *low_right;
  TMRQuadNode *up_left, *up_right;

  // The level required 
  int level;

  // The u/v location of the lower left-hand corner of the domain
  // in the parametric space
  uint32_t u, v;

  // The location at which to split incoming point
  double x, y;

  // The point numbers and arrays
  uint32_t num_points;
  uint32_t *pt_nums;
  double *pts;
};

/*
  The basic Triangle class that stores the node numbers associated
  with this triangle, tags and a metric for the triangle quality

  This class is not designed to be externally visible to the user of
  the TMRTriangularize class.
*/
class TMRTriangle {
 public:
  TMRTriangle(){
    u = v = w = 0;
    tag = status = 0;
  }
  TMRTriangle( uint32_t _u, uint32_t _v, uint32_t _w ){
    u = _u;  v = _v;  w = _w;
    tag = status = 0;
  }
  // The indices of this triangle
  uint32_t u, v, w;
  
  // Tag/info values (used to helpfully tag/label triangles)
  uint32_t tag;
  uint32_t status;
  
  // Quality metric
  float quality;
};

/*
  Triangularize the domain/surface using a frontal method based on
  Rebay's paper with Bowyer-Watson Delaunay triangularization
  algorithm.
*/
class TMRTriangularize : public TMREntity {
 public:
  TMRTriangularize( int npts, const double inpts[],
                    int nsegs, const int segs[],
                    TMRFace *surf=NULL );
  TMRTriangularize( int npts, const double inpts[], int nholes,
                    int nsegs, const int segs[],
                    TMRFace *surf=NULL );
  ~TMRTriangularize();

  // Create the frontal mesh with the given mesh spacing
  void frontal( double h );

  // Retrieve the mesh connectivity from the object
  void getMesh( int *_num_points, int *_num_triangles, 
                int **_conn, double **_pts, TMRPoint **_X );
  
  // Write the triangulation to an outputfile
  void writeToVTK( const char *filename );

 private:
  // The Bowyer-Watson algorithm is started with 4 points (2 triangles)
  // that cover the entire domain. These are deleted at the end 
  // of the algorithm.
  static const int FIXED_POINT_OFFSET = 4;

  // TAGS for the triangles
  static const uint32_t NO_STATUS = 0;
  static const uint32_t WAITING = 1;
  static const uint32_t ACTIVE = 2;
  static const uint32_t ACCEPTED = 3;
  static const uint32_t DELETE_ME = 4;

  // Initialize the underlying data structures
  void initialize( int npts, const double inpts[], int nholes,
                   int nsegs, const int segs[],
                   TMRFace *surf );

  // Add a point to the list -- this only adds a point to the list and 
  // returns the new point number, it does not add the point to the 
  uint32_t addPoint( const double pt[] );

  // Add a point to the mesh and re-triangularize the mesh
  // to account for the new point.
  void addPointToMesh( const double pt[] );
  void addPointToMesh( const double pt[], TMRTriangle *tri );

  // Get a hash value for the given edge
  inline uint32_t getEdgeHash( uint32_t u, uint32_t v );

  // Add/delete a triangle from the data structure
  int addTriangle( TMRTriangle tri );
  int deleteTriangle( TMRTriangle tri );

  // Get a hash value for the given triangle
  inline uint32_t getTriangleHash( TMRTriangle *tri );              

  // Add/delete triangles from the active set of triangles
  int addActiveTriangle( TMRTriangle *tri );
  int deleteActiveTriangle( TMRTriangle *tri );

  // Mark all the triangles in the list
  void setTriangleTags( uint32_t tag );
  void tagTriangles( TMRTriangle *tri );

  // Given the two ordered nodes, add the triangle to the list
  void completeMe( uint32_t u, uint32_t v, TMRTriangle **tri );
  void nodeToOneTriangle( uint32_t u, TMRTriangle *tri );

  // Dig cavity: Function a la Shewchuk
  void digCavity( uint32_t u, uint32_t v, uint32_t w );

  // Determine whether this edge is in the PSLG edge list
  void setUpPSLGEdges( int nsegs, const int segs[] ); 
  int edgeInPSLG( uint32_t u, uint32_t v );

  // Does this triangle enclose the point
  int enclosed( const double pt[], uint32_t u, uint32_t v, uint32_t w );
  double inCircle( uint32_t u, uint32_t v, uint32_t w, uint32_t x );

  // Find the enclosing triangle
  void findEnclosing( const double pt[], TMRTriangle **tri );

  // Compute the maximum edge length of the triangle
  double computeCircumcircle( TMRTriangle *tri );
  double computeMaxEdgeLength( TMRTriangle *tri );

  // Compute the intersection
  double computeIntersection( const double m[], const double e[], 
                              uint32_t u, uint32_t v, uint32_t w );

  // The underlying surface
  TMRFace *face;

  // Initial number of boundary points
  uint32_t init_boundary_points;

  // Offset to the points that will be removed from the mesh
  // these are the original background mesh and the holes
  uint32_t fixed_point_offset; 

  // Keep track of the points
  uint32_t num_points; // The current number of points
  uint32_t max_num_points; // The maximum number of points

  // Array of the points that have been set
  double *pts;
  TMRPoint *X;
  TMRTriangle **pts_to_tris;

  // The PSLG edges
  uint32_t num_pslg_edges;
  uint32_t *pslg_edges;

  // The root of the quadtree data structure
  TMRQuadDomain domain;
  TMRQuadNode *root;
  uint32_t search_tag;

  // Keep a doubly-linked list to store the added triangles. The list
  // is doubly-linked facilitate deleting triangles from the list.
  class TriListNode {
  public:
    TMRTriangle tri; // The actual triangle we've allocated
    TriListNode *next, *prev; // Next/previous list entries
  };

  // Keep track of the current set of triangles
  TriListNode *list_start, *list_end;
  TriListNode *list_marker;

  // Keep track of the number of triangles
  uint32_t num_triangles;

  // Keep a hash tabled based on the ordered edges of the triangular
  // mesh. The order must match the counter clockwise ordering of the
  // triangle, making the edge to triangle mapping unique. Each
  // triangle is stored three times within the hash table.
  class EdgeHashNode {
  public:
    uint32_t u, v; // The edge indices
    TriListNode *tri_node; // Pointer to the node within the triangle list
    EdgeHashNode *next; // Next node 
  };

  // Create the array of hash buckets
  EdgeHashNode **buckets;

  // The number of buckets
  int num_buckets;
  int num_hash_nodes;

  // Keep a hash table for the active triangles. This hash is based on the
  // nodes in the triangle. 
  class ActiveHashNode {
  public:
    TMRTriangle *tri;
    ActiveHashNode *next;
  };

  // Active hash buckets
  ActiveHashNode **active_buckets;

  // The number of active buckets/triangles
  int num_active_buckets;
  int num_active_triangles;
};

#endif // TRIANGULARIZE_H
