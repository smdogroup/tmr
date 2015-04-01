#ifndef TMR_OCTANT_H
#define TMR_OCTANT_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static const int TMR_MAX_LEVEL = 30;

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

  int x, y, z;
  int level;
};

/*
  A array of octants that may, or may not, be sorted and unique
*/
class TMROctantArray {
 public:
  TMROctantArray( TMROctant *array, int size );
  ~TMROctantArray();

  void getArray( TMROctant **_array, int *_size );
  void sort();
  void sortUnique();
  void merge( TMROctantArray * list );

 private:
  int is_sorted, is_unique;
  int size, max_size;
  TMROctant *array;
};

/*
  Create a list of octants
*/
class TMROctantQueue {
 public:
  TMROctantQueue();

  int length();
  void push( TMROctant *oct );
  TMROctant pop();
  TMROctantArray* toArray();
  
 private:
  class OcQueueNode {
  public:
    OcQueueNode(){
      next = NULL;
    }
    
    TMROctant oct;
    OcQueueNode *next;
  };

  int num_elems;
  OcQueueNode *root, *tip;
};

/*
  Build a hash table based on the octant ordering
*/
class TMROctantHash {
 public:
  TMROctantHash();
  ~TMROctantHash();

  TMROctantArray * toArray();
  int addOctant( TMROctant *oct );

 private:
  // The minimum bucket size
  static const int min_num_buckets = (1 << 10)-1;

  class OcHashNode {
  public:
    OcHashNode(){
      next = NULL;
    }
    
    TMROctant oct;
    OcHashNode *next;
  };

  // Keep track of the bucket size
  int num_buckets;
  OcHashNode **hash_buckets;
  
  // Keep track of the number of elements
  int num_elems;

  // Get the buckets to place the octant in
  int getBucket( TMROctant *oct );
};

#endif // TMR_OCTANT_H
