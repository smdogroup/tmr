#ifndef TMR_TRIANGULARIZE_H
#define TMR_TRIANGULARIZE_H

#include "TMRBase.h"




class TMRTriangle {
 public:
  // The indices of this triangle
  int32_t u, v, w;
  
  // Tag for special values
  int32_t tag;
};

/*
  Hash table to a list of triangles, organized by edge
*/
class TMRTriangularize : public TMREntity {
 public:


  
 private:
  void addVertex( const double pt[] );
  void addVertexFrontal( const double pt[], int32_t tri );



  void addTriangle( TMRTriangle tri );
  void deleteTriangle( TMRTriangle tri );

  // Given the two ordered nodes, add the triangle to the list
  void completeMe( int32_t u, int32_t v, TMRTriangle *tri );
  void nodeToOneTriangle( int32_t u, TMRTriangle *tri );


  // Keep a doubly-linked list to store the added triangles. The list
  // is doubly-linked facilitate deleting triangles from the list.
  class TriListNode {
  public:
    TMRTriangle tri; // The actual triangle we've allocated
    TriListNode *next, *prev; // Next/previous list entries
  };

  // Keep track of the current 
  TriListNode *list_root, *list_current;

  // Keep a hash tabled based on the ordered edges of the triangular
  // mesh. The order must match the counter clockwise ordering of the
  // triangle, making the edge to triangle mapping unique. Each
  // triangle is stored three times within the hash table.
  class EdgeHashNode {
  public:
    int32_t u, v; // The edge indices
    TriListNode *tri_node; // Pointer to the node within the triangle list
    TriHashNode *next; // Next node 
  };

  // Create the array of hash buckets
  EdgeHashNode **hash;
};


/*
  Domain information for 
*/
class TMRQuadDomain {
 public:
  double xlow, xhigh;
  double ylow, yhigh; 
};

/*
  The following class implments a simple quadtree data structure for fast
  geometric searching queries.
*/
class TMRQuadNode {
 public:
  static const int MAX_DEPTH = 30;
  static const int NODES_PER_LEVEL = 10;

  TMRQuadNode( QuadNodeDomain *_domain,
               int32_t _u, int32_t _v, int level );
  ~TMRQuadNode();
 
  // Add/delete nodes from the quadtree
  // ----------------------------------
  void addNode( int num, const double pt[] );
  int deleteNode( int num, const double pt[] );

  // Find the closest indexed point to the provided (x,y) location
  int findClosest( const double pt[], double *_dist );

 private:
  
  // The recursive call to the quadtree data structure
  void findClosest( const double pt[], int *index, double *dist );
  // The domain used 
  TMRQuadDomain *domain;

  // The level required 
  int level;

  // The u/v location of the lower left-hand corner of the domain
  // in the parametric space
  int32_t u, v;

  // The location at which to split incoming point
  double x, y;

  // The point numbers and arrays
  int num_points;
  int *pt_nums;
  double *pts;
};

#endif // TRIANGULARIZE_H
