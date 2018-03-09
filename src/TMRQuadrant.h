/*
  This file is part of the package TMR for adaptive mesh refinement.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#ifndef TMR_QUADRANT_H
#define TMR_QUADRANT_H

#include "TMRBase.h"

/*
  The TMR Quadrant class

  This class defines an quadrant that is used to order both the elements
  and nodes within the mesh. The methods can be used to compare
  quadrants, find the parent, child identification number, and find
  neighbours.
*/
class TMRQuadrant {
 public:
  int childId();
  void getSibling( int id, TMRQuadrant *sib );
  void parent( TMRQuadrant *parent );
  void edgeNeighbor( int edge, TMRQuadrant *neighbor );
  void cornerNeighbor( int corner, TMRQuadrant *neighbor );
  int compare( const TMRQuadrant *quadrant ) const;
  int comparePosition( const TMRQuadrant *quadrant ) const;
  int compareNode( const TMRQuadrant *quadrant ) const;
  int contains( TMRQuadrant *quad );

  int32_t face; // The face owner
  int32_t x, y; // The x,y coordinates
  int32_t tag; // A tag to store additional data
  int16_t level; // The refinement level
  int16_t info; // The info about faces
};

/*
  A array of quadrants that may or may not be sorted 
  
  When the array is sorted, the quadrants are made unique by discarding
  quadrants with a smaller level (that have larger side lengths).  After
  the array is sorted, it is searchable either based on elements (when
  use_nodes=0) or by node (use_nodes=1). The difference is that the
  node search ignores the mesh level.
*/
class TMRQuadrantArray {
 public:
  TMRQuadrantArray( TMRQuadrant *array, int size, 
                    int _use_node_index=0 );
  ~TMRQuadrantArray();

  TMRQuadrantArray* duplicate();
  void getArray( TMRQuadrant **_array, int *_size );
  void sort();
  TMRQuadrant* contains( TMRQuadrant *q, const int use_position=0 );
  void merge( TMRQuadrantArray * list );

 private:
  int use_node_index;
  int is_sorted;
  int size, max_size;
  TMRQuadrant *array;
};

/*
  Create a queue of quadrants

  This class defines a queue of quadrants that are used for the balance
  and coarsen operations.
*/
class TMRQuadrantQueue {
 public:
  TMRQuadrantQueue();
  ~TMRQuadrantQueue();

  int length();
  void push( TMRQuadrant *quad );
  TMRQuadrant pop();
  TMRQuadrantArray* toArray();
  
 private:
  // Class that defines an element within the queue
  class QuadQueueNode {
  public:
    QuadQueueNode(){ next = NULL; }
    TMRQuadrant quad;
    QuadQueueNode *next;
  };

  // Keep track of the number of elements in the queue
  int num_elems;
  QuadQueueNode *root, *tip;
};

/*
  Build a hash table based on the Morton ordering

  This object enables the creation of a unique set of quadrants such
  that no two have the same position/level combination. This hash
  table can then be made into an array of unique elements or nodes.
*/
class TMRQuadrantHash {
 public:
  TMRQuadrantHash( int _use_node_index=0 );
  ~TMRQuadrantHash();

  TMRQuadrantArray* toArray();
  int addQuadrant( TMRQuadrant *quad );

 private:
  // The minimum bucket size
  static const int min_num_buckets = (1 << 12)-1;

  class QuadHashNode {
  public:
    QuadHashNode(){ next = NULL; }
    TMRQuadrant quad;
    QuadHashNode *next;
  };

  // Set the element index/node
  int use_node_index;

  // Keep track of the bucket size
  int num_buckets;
  QuadHashNode **hash_buckets;
  
  // Keep track of the number of elements
  int num_elems;

  // Get the buckets to place the quadrant in
  int getBucket( TMRQuadrant *quad );
};

#endif // TMR_QUADRANT_H
