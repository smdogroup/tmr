#include <stdio.h>
#include "TMRQuadtree.h"

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

/*
  Refine the initial quadtree to the specified depth along all
  coordinate directions 
*/
TMRQuadtree::TMRQuadtree( int refine_level ){
  // Check that the refinement lies within legal bounds
  if (refine_level < 0){
    refine_level = 0;
  }
  else if (refine_level-1 > TMR_MAX_LEVEL){
    refine_level = TMR_MAX_LEVEL-1;
  }

  // Compute the refinement level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - refine_level);

  // Compute the number of quadrants along each edge
  int32_t nx = 1 << refine_level;
  
  // Set the number of quadrants created
  int nquads = nx*nx;

  // Create an array of the quadrants that will be stored
  TMRQuadrant *array = new TMRQuadrant[ nquads ];

  // Create all the quadrants required
  int index = 0;
  for ( int32_t y = 0; y < hmax; y += h ){
    for ( int32_t x = 0; x < hmax; x += h ){
      array[index].x = x;
      array[index].y = y;
      array[index].level = refine_level;
      index++;
    }  
  }

  // Sort the array of quadrants
  elements = new TMRQuadrantArray(array, nquads);
  elements->sort();

  // Zero the array of quadrant nodes
  nodes = NULL;

  // Zero everything else
  order = 2;
}

/*
  Generate a random initial quadtree that can be used for testing.
*/
TMRQuadtree::TMRQuadtree( int nrand, int min_level, int max_level ){
  // Create an array of the quadrants that will be stored
  TMRQuadrant *array = new TMRQuadrant[ nrand ];

  // Set the maximum refinement level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Generate a random number of quadrants along random directions
  for ( int i = 0; i < nrand; i++ ){
    int32_t level = min_level + (rand() % (max_level - min_level + 1));

    const int32_t h = 1 << (TMR_MAX_LEVEL - level);
    int32_t x = h*(rand() % (1 << level));
    int32_t y = h*(rand() % (1 << level));

    if (x >= 0 && x < hmax && y >= 0 && y < hmax){ 
      array[i].x = x;
      array[i].y = y;
      array[i].level = level;
    }
  }

  elements = new TMRQuadrantArray(array, nrand);
  elements->sort();

  // Zero the array of quadrant nodes
  nodes = NULL;

  // Set the default order of the mesh
  order = 2;
}

/*
  Create the quadtree tree from the array
*/
TMRQuadtree::TMRQuadtree( TMRQuadrantArray *_elements ){
  elements = _elements;
  elements->sort();

  // Zero the array of quadrant nodes
  nodes = NULL;

  // Set the default mesh order
  order = 2;
}

/*
  Free the quadtree
*/
TMRQuadtree::~TMRQuadtree(){
  // Free the elements and nodes
  if (elements){ delete elements; }
  if (nodes){ delete nodes; }
}

/*
  Refine the quadtree by adding/subtracting elements

  This code takes in an array of integers that are either positive,
  negative or zero. A positive value means that the cell is refined. A
  negative index means that only the parent of the quadrant is
  included. A zero index indicates that the quadrant is retained.

  If a NULL array is passed in, all quadrants within the mesh are
  refined.

  input:
  refinement:     the refinement level for each element (may be null)
  min_level:      the minimum refinement level
  max_level:      the maximum refinement level
*/
void TMRQuadtree::refine( const int refinement[], 
                          int min_level, int max_level ){
  // Adjust the min and max levels to ensure consistency
  if (min_level < 0){ min_level = 0; }
  if (max_level > TMR_MAX_LEVEL){ max_level = TMR_MAX_LEVEL; }

  // This is just a sanity check
  if (min_level > max_level){ min_level = max_level; }

  // If the nodes were defined, delete them
  if (nodes){ 
    delete nodes; 
    nodes = NULL; 
  }

  // Set the maximum side length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Create a hash table for the refined  t
  TMRQuadrantHash *hash = new TMRQuadrantHash();

  // Get the current array of quadrants
  int size;
  TMRQuadrant *array;
  elements->getArray(&array, &size);

  // Add the element
  if (refinement){
    for ( int i = 0; i < size; i++ ){
      if (refinement[i] == 0){
        TMRQuadrant q;
        array[i].getSibling(0, &q);
        hash->addQuadrant(&q);
      }
      else if (refinement[i] < 0){
        if (array[i].level > min_level){
          TMRQuadrant q;
          array[i].getSibling(0, &q);
          q.level = q.level-1;
          hash->addQuadrant(&q);
        }
        else {
          hash->addQuadrant(&array[i]);
        }
      }
      else if (refinement[i] > 0){
        if (array[i].level < max_level){
          TMRQuadrant c = array[i];
          c.level += 1;
          hash->addQuadrant(&c);
        }
        else {
          hash->addQuadrant(&array[i]);
        }
      }
    }
  }
  else {
    for ( int i = 0; i < size; i++ ){
      if (array[i].level < max_level){
        TMRQuadrant c = array[i];
        c.level += 1;
        hash->addQuadrant(&c);
      }
      else {
        hash->addQuadrant(&array[i]);
      }
    }
  }

  // Now convert the new hash
  TMRQuadrantArray *child0_elems = hash->toArray();
  child0_elems->getArray(&array, &size);

  // Loop over all elements and add their siblings
  for ( int i = 0; i < size; i++ ){
    for ( int j = 0; j < 4; j++ ){
      TMRQuadrant q;
      array[i].getSibling(j, &q);
      if ((q.x >= 0 && q.x < hmax) &&
          (q.y >= 0 && q.y < hmax)){
        hash->addQuadrant(&q);
      }
    }
  }

  // Free the temporary elements
  delete child0_elems;

  // Free the old elements class
  if (elements){ delete elements; }

  // Cover the hash table to a list and uniquely sort it
  elements = hash->toArray();
  elements->sort();

  delete hash;
}

/*
  Coarsen the quadtree

  This algorithm implements a simple algorithm to obtain a coarser
  quadtree mesh from the current quadtree. For each mesh, if all the
  children of a common parent exist, then the parent quadrant is added
  to the coarse mesh. This is repeated for all quadrants that are
  added to a hash. The final hash is then made into a list and
  returned as a new quadtree object.
*/
TMRQuadtree *TMRQuadtree::coarsen(){
  int size;
  TMRQuadrant *array;
  elements->getArray(&array, &size);

  // Set the offset to be 2**d-1
  const int offset = (1 << 2) - 1;

  TMRQuadrantQueue *queue = new TMRQuadrantQueue();
  int coarsened = 0;

  // Scan through the list, if we have offset quadrants which
  // all share the same parent, then coarsen the parent
  for ( int i = 0; i < size; i++ ){
    int same_parent = 0;
    if (array[i].level > 0 && 
	array[i].childId() == 0 && 
	i+offset < size &&
	array[i+offset].childId() == offset){

      // Test to see if the children have the same parent
      TMRQuadrant p;
      array[i+offset].getSibling(0, &p);
      if (array[i].compare(&p) == 0){
	same_parent = 1;
	array[i].parent(&p);
	queue->push(&p);

	i += offset;
	coarsened++;
      }
    }

    if (!same_parent){
      queue->push(&array[i]);
    }
  }

  // Create the new quadtree from the queue
  TMRQuadtree *tree = new TMRQuadtree(queue->toArray());
  delete queue;

  return tree;
}

/*
  Find the quadrant in the element quadrant list that completely contains
  the provided quadrant.

  This code can be used to find quadrants that completely contain
  another quadrant. This can be useful for finding relationships between
  two quadtrees.
  
  input: 
  quad:     the candidate quadrant

  returns: the occtant that was found or NULL if no such quadrant exists
*/
TMRQuadrant* TMRQuadtree::findEnclosing( TMRQuadrant *quad ){
  // Retrieve the array of elements
  int size = 0;
  TMRQuadrant *elems = NULL;
  elements->getArray(&elems, &size);

  // Set the lower and upper bounds for the quadrant
  const int32_t hquad = 1 << (TMR_MAX_LEVEL - quad->level);
  const int32_t x1 = quad->x;
  const int32_t y1 = quad->y;
  const int32_t x2 = quad->x + hquad;
  const int32_t y2 = quad->y + hquad;
  
  // Set the low and high indices to the first and last
  // element of the element array
  int low = 0;
  int high = size-1;
  int mid = low + (int)((high - low)/2);
	
  // Maintain values of low/high and mid such that the
  // quadrant is between (elems[low], elems[high]).
  // Note that if high-low=1, then mid = high
  while (high != mid){
    // Check if elems[mid] contains the provided quadrant
    const int32_t h = 1 << (TMR_MAX_LEVEL - elems[mid].level);
    if ((elems[mid].x <= x1 && x2 <= elems[mid].x+h) &&
	(elems[mid].y <= y1 && y2 <= elems[mid].y+h)){
      return &elems[mid];
    }
    
    // Compare the ordering of the two quadrants - if the
    // quadrant is less than the other, then adjust the mid point 
    if (quad->compare(&elems[mid]) < 0){
      high = mid-1;
    } 
    else {
      low = mid+1;
    }
    
    // Re compute the mid-point and repeat
    mid = high - (int)((high - low)/2);
  }

  // Check if elems[mid] contains the provided quadrant
  const int32_t h1 = 1 << (TMR_MAX_LEVEL - elems[mid].level);
  if ((elems[mid].x <= x1 && x2 <= elems[mid].x+h1) &&
      (elems[mid].y <= y1 && y2 <= elems[mid].y+h1)){
    return &elems[mid];
  }

  // Check if elems[mid] contains the provided quadrant
  const int32_t h2 = 1 << (TMR_MAX_LEVEL - elems[low].level);
  if ((elems[low].x <= x1 && x2 <= elems[low].x+h2) &&
      (elems[low].y <= y1 && y2 <= elems[low].y+h2)){
    return &elems[low];
  }

  // No quadrant was found, return NULL
  return NULL;
}

/*
  Find the range of element indices that enclose an element

  Unlike the findEnclosing, this will always return a well-defined
  range unless the quadrant tree has not bee balanaced or the quadrants
  fall outside the domain.

  input:
  quad:    the quadrant that we are trying to enclose

  output:
  low:    the lower index within the element array
  high:   the higher index within the element array
*/
void TMRQuadtree::findEnclosingRange( TMRQuadrant *quad,
                                      int *low, int *high ){
  // Set the lower index
  *low = 0;
  
  // Get the number of elements
  elements->getArray(NULL, high);

  // Find the maximum level
  int32_t h = 1 << (TMR_MAX_LEVEL - quad->level);

  // Set the node 
  TMRQuadrant p = *quad;
  p.level = TMR_MAX_LEVEL;

  // Find the lower index
  TMRQuadrant *elow = findEnclosing(&p);
  if (elow){
    *low = elow->tag;
  }

  // Switch the quadrant to the upper-most quadrant and set
  // the range to the difference between the two
  p.x += h-1;
  p.y += h-1;

  // Find the upper index
  TMRQuadrant *ehigh = findEnclosing(&p);
  if (ehigh){
    *high = ehigh->tag+1;
  }
}

/*
  Create the mesh using the internal quadtree data

  This algorithm first adds all the nodes from the element to a hash
  object that keeps a unique list of nodes.  After all the nodes are
  added, they are sorted and uniquified. Next, we go through the
  elements and label any possible dependent node.
*/
void TMRQuadtree::createNodes( int _order ){
  order = _order;
  if (order < 2){ order = 2; }
  if (order > 3){ order = 3; }

  // Free the existing node array
  if (nodes){ delete nodes; }

  // Get the current array of quadrants
  int size;
  TMRQuadrant *array;
  elements->getArray(&array, &size);

  // Allocate an array large enough to store all of the nodes
  int index = 0;
  TMRQuadrant *all_nodes = new TMRQuadrant[ order*order*size ];

  // Loop over all of the current quadrants and add the nodes
  for ( int i = 0; i < size; i++ ){
    if (order == 2){
      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);

      // Add all of the nodes to the hash
      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          // Set the node level to the highest level - this will
          // be updated when the nodes are assigned to the elements
          all_nodes[index].level = 0;
            
          // Set a positive tag, this will be replaced with a 
          // negative tag if the node is dependent
          all_nodes[index].tag = 1;
          
          // Set the node coordinates
          all_nodes[index].x = array[i].x + ii*h;
          all_nodes[index].y = array[i].y + jj*h;
          index++;
        }
      }
    }
    else if (order == 3){
      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level - 1);

      // Add all of the nodes to the hash
      for ( int jj = 0; jj < 3; jj++ ){
        for ( int ii = 0; ii < 3; ii++ ){
          // Set the node level to the highest level - this will
          // be updated when the nodes are assigned to the elements
          all_nodes[index].level = 0;
          
          // Set a positive tag, this will be replaced with a 
          // negative tag if the node is dependent
          all_nodes[index].tag = 1;
          
          // Set the node coordinates
          all_nodes[index].x = array[i].x + ii*h;
          all_nodes[index].y = array[i].y + jj*h;
          index++;
        }
      }
    }
  }

  // Create an array of all the quadrants and uniquely sort it
  nodes = new TMRQuadrantArray(all_nodes, index);
  nodes->sort();
}

/*
  Print out the quadtree to a file for visualization
*/
void TMRQuadtree::printQuadtree( const char *filename ){
  // Iterate through and print out everything  
  FILE * fp = fopen(filename, "w");
  if (fp){
    int size;
    TMRQuadrant *array;
    elements->getArray(&array, &size);

    fprintf(fp, "Variables = X, Y\n");
    fprintf(fp, "ZONE T=TMR N=%d E=%d ", 4*size, size);
    fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEQUADRILATERAL\n");

    // Set the length of the cube
    double dh = 1.0/(1 << TMR_MAX_LEVEL);

    for ( int i = 0; i < size; i++ ){
      int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      int x = array[i].x;
      int y = array[i].y;

      fprintf(fp, "%e %e\n", x*dh, y*dh);
      fprintf(fp, "%e %e\n", (x+h)*dh, y*dh);
      fprintf(fp, "%e %e\n", (x+h)*dh, (y+h)*dh);
      fprintf(fp, "%e %e\n", x*dh, (y+h)*dh);
    }

    for ( int i = 0; i < size; i++ ){
      for ( int k = 0; k < 4; k++ ){
	fprintf(fp, "%d ", 4*i+k+1);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}


