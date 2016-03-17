#ifndef TMR_OCTANT_H
#define TMR_OCTANT_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

/*
  The following constants define the maximum octant depth and maximum
  element order within the code.
*/
static const int TMR_MAX_NODE_LEVEL = 30;
static const int TMR_LOG2_MAX_ELEMENT_ORDER = 3; 
static const int TMR_MAX_LEVEL = 
  TMR_MAX_NODE_LEVEL - TMR_LOG2_MAX_ELEMENT_ORDER;

/*
  The TMR Octant class

  This class defines an octant that is used to order both the elements
  and nodes within the mesh. The methods can be used to compare
  octants, find the parent, child identification number, and find
  neighbours.
*/
class TMROctant {
 public:
  TMROctant(){}
  ~TMROctant(){}

  int childId();
  void getSibling( int id, TMROctant *sib );
  void parent( TMROctant *parent );
  void faceNeighbor( int face, TMROctant *neighbor );
  void edgeNeighbor( int edge, TMROctant *neighbor );
  void cornerNeighbor( int corner, TMROctant *neighbor );
  int compare( const TMROctant *octant ) const;
  int compareEncoding( const TMROctant *octant ) const;

  int32_t x, y, z; // The x,y,z coordinates
  uint32_t level; // The refinement level
  int32_t tag; // A tag to store additional data
};

/*
  A array of octants that may or may not be sorted 
  
  When the array is sorted, the octants are made unique by discarding
  octants with a smaller level (that have larger side lengths).  After
  the array is sorted, it is searchable either based on elements (when
  use_nodes=0) or by node (use_nodes=1). The difference is that the
  node search ignores the mesh level.
*/
class TMROctantArray {
 public:
  TMROctantArray( TMROctant *array, int size );
  ~TMROctantArray();

  void getArray( TMROctant **_array, int *_size );
  void sort();
  TMROctant* contains( TMROctant *q, int use_nodes=0 );
  void merge( TMROctantArray * list );

 private:
  int is_sorted;
  int size, max_size;
  TMROctant *array;
};

/*
  Create a queue of octants

  This class defines a queue of octants that are used for the balance
  and coarsen operations.
*/
class TMROctantQueue {
 public:
  TMROctantQueue();
  ~TMROctantQueue();

  int length();
  void push( TMROctant *oct );
  TMROctant pop();
  TMROctantArray* toArray();
  
 private:
  // Class that defines an element within the queue
  class OctQueueNode {
  public:
    OctQueueNode(){ next = NULL; }
    TMROctant oct;
    OctQueueNode *next;
  };

  // Keep track of the number of elements in the queue
  int num_elems;
  OctQueueNode *root, *tip;
};

/*
  Build a hash table based on the Morton ordering

  This object enables the creation of a unique set of octants such
  that no two have the same position/level combination. This hash
  table can then be made into an array of unique elements or nodes.
*/
class TMROctantHash {
 public:
  TMROctantHash();
  ~TMROctantHash();

  TMROctantArray * toArray();
  int addOctant( TMROctant *oct );

 private:
  // The minimum bucket size
  static const int min_num_buckets = (1 << 12)-1;

  class OctHashNode {
  public:
    OctHashNode(){ next = NULL; }
    TMROctant oct;
    OctHashNode *next;
  };

  // Keep track of the bucket size
  int num_buckets;
  OctHashNode **hash_buckets;
  
  // Keep track of the number of elements
  int num_elems;

  // Get the buckets to place the octant in
  int getBucket( TMROctant *oct );
};

#endif // TMR_OCTANT_H
