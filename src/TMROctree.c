#include "TMROctree.h"

/*
  Refine the initial octree to the specified depth along all
  coordinate directions 
*/
TMROctree::TMROctree( int refine_level ){
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

  // Keep track of the octants and the length of the octant array
  int nocts = 0;
  TMROctant *array = NULL;
  
  // Compute the number of octants along each edge
  int32_t nx = 1 << refine_level;
    
  // Set the number of octants created
  nocts = nx*nx*nx;

  // Create an array of the octants that will be stored
  array = new TMROctant[ nocts ];

  // Create all the octants required
  int index = 0;
  for ( int32_t z = 0; z < hmax; z += h ){
    for ( int32_t y = 0; y < hmax; y += h ){
      for ( int32_t x = 0; x < hmax; x += h ){
        array[index].x = x;
        array[index].y = y;
        array[index].z = z;
        array[index].level = refine_level;
        index++;
      }
    }
  }

  // Sort the array of octants
  elements = new TMROctantArray(array, nocts);
  elements->sort();

  // Zero the array of octant nodes
  nodes = NULL;
  X = NULL;

  // Set the default order to 2
  order = 2;
}

/*
  Generate a random initial octree that can be used for testing.
*/
TMROctree::TMROctree( int nrand, int min_level, int max_level ){
  // Create an array of the octants that will be stored
  TMROctant *array = new TMROctant[ nrand ];

  // Generate a random number of octants along random directions
  for ( int i = 0; i < nrand; i++ ){
    int32_t level = min_level + (rand() % (max_level - min_level + 1));

    const int32_t h = 1 << (TMR_MAX_LEVEL - level);
    int32_t x = h*(rand() % (1 << level));
    int32_t y = h*(rand() % (1 << level));
    int32_t z = h*(rand() % (1 << level));

    array[i].x = x;
    array[i].y = y;
    array[i].z = z;
    array[i].level = level;
  }

  elements = new TMROctantArray(array, nrand);
  elements->sort();

  // Zero the array of octant nodes
  nodes = NULL;
  X = NULL;

  // Zero everything else
  order = 2;
}

/*
  Create the octree tree from the array
*/
TMROctree::TMROctree( TMROctantArray *_elements ){
  elements = _elements;
  elements->sort();

  // Zero the array of octant nodes
  nodes = NULL;
  X = NULL;

  // Zero everything else
  order = 2;
}

/*
  Free the octree
*/
TMROctree::~TMROctree(){
  // Free the elements and nodes
  if (elements){ delete elements; }
  if (nodes){ delete nodes; }
  if (X){ delete [] X; }
}

/*
  Refine the octree by adding/subtracting elements

  This code takes in an array of integers that are either positive,
  negative or zero. A positive value means that the cell is refined,
  and all children are added. A negative index means that only the
  parent of the octant is refined. A zero index indicates that the
  octant is retained.
*/
void TMROctree::refine( int refinement[], 
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
  if (X){
    delete [] X;
    X = NULL;
  }

  // Create a hash table for the refined  t
  TMROctantHash *hash = new TMROctantHash();

  // Get the current array of octants
  int size;
  TMROctant *array;
  elements->getArray(&array, &size);

  // Add the element
  for ( int i = 0; i < size; i++ ){
    if (refinement[i] == 0){
      // Try adding the 0-sibling
      TMROctant q;
      array[i].getSibling(0, &q);
      hash->addOctant(&q);
    }
    else if (refinement[i] < 0){
      if (array[i].level > min_level){
	TMROctant q;
        array[i].getSibling(0, &q);
        q.level = q.level-1;
	hash->addOctant(&q);
      }
      else {
	hash->addOctant(&array[i]);
      }
    }
    else if (refinement[i] > 0){
      if (array[i].level < max_level){
	TMROctant q = array[i];
	q.level += 1;
        hash->addOctant(&q);
      }
      else {
	hash->addOctant(&array[i]);
      }
    }
  }

  // Now convert the new hash
  TMROctantArray *child0_elems = hash->toArray();
  child0_elems->getArray(&array, &size);

  // Loop over all elements and add their siblings
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  for ( int i = 0; i < size; i++ ){
    for ( int j = 0; j < 8; j++ ){
      TMROctant q;
      array[i].getSibling(j, &q);
      if ((q.x >= 0 && q.x < hmax) &&
          (q.y >= 0 && q.y < hmax) &&
          (q.z >= 0 && q.z < hmax)){
        hash->addOctant(&q);
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
  Coarsen the octree

  This algorithm implements a simple algorithm to obtain a coarser
  octree mesh from the current octree. For each mesh, if all the
  children of a common parent exist, then the parent octant is added
  to the coarse mesh. This is repeated for all octants that are added
  to a hash. The final hash is then made into a list and returned
  as a new octree object.
*/
TMROctree *TMROctree::coarsen(){
  int size;
  TMROctant *array;
  elements->getArray(&array, &size);

  // Set the offset to be 2**d-1
  int offset = (1 << 3) - 1;

  TMROctantQueue *queue = new TMROctantQueue();
  int coarsened = 0;

  // Scan through the list, if we have offset octants which
  // all share the same parent, then coarsen the parent
  for ( int i = 0; i < size; i++ ){
    int same_parent = 0;
    if (array[i].level > 0 && 
	array[i].childId() == 0 && 
	i+offset < size &&
	array[i+offset].childId() == offset){

      // Test to see if the children have the same parent
      TMROctant p;
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

  // Create the new octree from the queue
  TMROctree *tree = new TMROctree(queue->toArray());
  delete queue;

  return tree;
}

/*
  Find the octant in the element octant list that completely contains
  the provided octant.

  This code can be used to find octants that completely contain
  another octant. This can be useful for finding relationships between
  two octrees.
  
  input: 
  oct:     the candidate octant

  returns: the occtant that was found or NULL if no such octant exists
*/
TMROctant* TMROctree::findEnclosing( TMROctant *oct ){
  // Retrieve the array of elements
  int size = 0;
  TMROctant *elems = NULL;
  elements->getArray(&elems, &size);

  // Set the lower and upper bounds for the octant
  const int32_t hoct = 1 << (TMR_MAX_LEVEL - oct->level);
  const int32_t x1 = oct->x;
  const int32_t y1 = oct->y;
  const int32_t z1 = oct->z;
  const int32_t x2 = oct->x + hoct;
  const int32_t y2 = oct->y + hoct;
  const int32_t z2 = oct->z + hoct;

  // Set the low and high indices to the first and last
  // element of the element array
  int low = 0;
  int high = size-1;
  int mid = low + (int)((high - low)/2);
	
  // Maintain values of low/high and mid such that the
  // octant is between (elems[low], elems[high]).
  // Note that if high-low=1, then mid = high
  while (high != mid){
    // Check if elems[mid] contains the provided octant
    const int32_t h = 1 << (TMR_MAX_LEVEL - elems[mid].level);
    if ((elems[mid].x <= x1 && x2 <= elems[mid].x+h) &&
	(elems[mid].y <= y1 && y2 <= elems[mid].y+h) &&
	(elems[mid].z <= z1 && z2 <= elems[mid].z+h)){
      return &elems[mid];
    }
    
    // Compare the ordering of the two octants - if the
    // octant is less than the other, then adjust the mid point 
    if (oct->compare(&elems[mid]) < 0){
      high = mid-1;
    } 
    else {
      low = mid+1;
    }
    
    // Re compute the mid-point and repeat
    mid = high - (int)((high - low)/2);
  }

  // Check if elems[mid] contains the provided octant
  const int32_t h1 = 1 << (TMR_MAX_LEVEL - elems[mid].level);
  if ((elems[mid].x <= x1 && x2 <= elems[mid].x+h1) &&
      (elems[mid].y <= y1 && y2 <= elems[mid].y+h1) &&
      (elems[mid].z <= z1 && z2 <= elems[mid].z+h1)){
    return &elems[mid];
  }

  // Check if elems[mid] contains the provided octant
  const int32_t h2 = 1 << (TMR_MAX_LEVEL - elems[low].level);
  if ((elems[low].x <= x1 && x2 <= elems[low].x+h2) &&
      (elems[low].y <= y1 && y2 <= elems[low].y+h2) &&
      (elems[low].z <= z1 && z2 <= elems[low].z+h2)){
    return &elems[low];
  }

  // No octant was found, return NULL
  return NULL;
}

/*
  Find the range of element indices that enclose an element

  Unlike the findEnclosing, this will always return a well-defined
  range unless the octant tree has not bee balanaced or the octants
  fall outside the domain.

  input:
  oct:    the octant that we are trying to enclose

  output:
  low:    the lower index within the element array
  high:   the higher index within the element array
*/
void TMROctree::findEnclosingRange( TMROctant *oct,
				    int *low, int *high ){
  *low = 0;
  elements->getArray(NULL, high);

  // Find the maximum level
  int32_t h = 1 << (TMR_MAX_LEVEL - oct->level);

  // Set the node 
  TMROctant p = *oct;
  p.level = TMR_MAX_LEVEL;

  // Find the lower index
  TMROctant *elow = findEnclosing(&p);
  if (elow){
    *low = elow->tag;
  }

  // Switch the octant to the upper-most octant and set
  // the range to the difference between the two
  p.x += h-1;
  p.y += h-1;
  p.z += h-1;

  // Find the upper index
  TMROctant *ehigh = findEnclosing(&p);
  if (ehigh){
    *high = ehigh->tag+1;
  }
}

/*
  Create the mesh using the internal octree data

  This algorithm first adds all the nodes from the element to a hash
  object that keeps a unique list of nodes.  After all the nodes are
  added, they are sorted and uniquified. Next, we go through the
  elements and label any possible dependent node.
*/
void TMROctree::createNodes( int _order ){
  order = _order;
  if (order < 2){ order = 2; }
  if (order > 3){ order = 3; }

  // Free the existing node array
  if (nodes){ delete nodes; }
  if (X){ delete [] X; }

  // Get the current array of octants
  int size;
  TMROctant *array;
  elements->getArray(&array, &size);

  // Allocate an array large enough to store all of the nodes
  int index = 0;
  TMROctant *all_nodes = new TMROctant[ order*order*order*size ];

  // Loop over all of the current octants and add the nodes
  for ( int i = 0; i < size; i++ ){
    if (order == 2){
      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);

      // Add all of the nodes to the hash
      for ( int kk = 0; kk < 2; kk++ ){
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
            all_nodes[index].z = array[i].z + kk*h;
            index++;
          }
        }
      }
    }
    else if (order == 3){
      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level - 1);

      // Add all of the nodes to the hash
      for ( int kk = 0; kk < 3; kk++ ){
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
            all_nodes[index].z = array[i].z + kk*h;
            index++;
          }
        }
      }
    }
  }

  // Create an array of all the octants and uniquely sort it
  nodes = new TMROctantArray(all_nodes, index);
  nodes->sort();

  // Allocate the nodes
  nodes->getArray(NULL, &size);
  X = new TMRPoint[ size ];
  memset(X, 0, size*sizeof(TMRPoint));
}

/*
  Print out the octree to a file for visualization
*/
void TMROctree::printOctree( const char * filename ){
  // Iterate through and print out everything  
  FILE * fp = fopen(filename, "w");
  if (fp){
    int size;
    TMROctant *array;
    elements->getArray(&array, &size);

    fprintf(fp, "Variables = X, Y, Z\n");
    fprintf(fp, "ZONE T=TMR N=%d E=%d ", 8*size, size);
    fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEBRICK\n");

    // Set the length of the cube
    double dh = 1.0/(1 << TMR_MAX_LEVEL);

    for ( int i = 0; i < size; i++ ){
      int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      int x = array[i].x;
      int y = array[i].y;
      int z = array[i].z;

      fprintf(fp, "%e %e %e\n", x*dh, y*dh, z*dh);
      fprintf(fp, "%e %e %e\n", (x+h)*dh, y*dh, z*dh);
      fprintf(fp, "%e %e %e\n", (x+h)*dh, (y+h)*dh, z*dh);
      fprintf(fp, "%e %e %e\n", x*dh, (y+h)*dh, z*dh);

      fprintf(fp, "%e %e %e\n", x*dh, y*dh, (z+h)*dh);
      fprintf(fp, "%e %e %e\n", (x+h)*dh, y*dh, (z+h)*dh);
      fprintf(fp, "%e %e %e\n", (x+h)*dh, (y+h)*dh, (z+h)*dh);
      fprintf(fp, "%e %e %e\n", x*dh, (y+h)*dh, (z+h)*dh);
    }

    for ( int i = 0; i < size; i++ ){
      for ( int k = 0; k < 8; k++ ){
	fprintf(fp, "%d ", 8*i+k+1);
      }
      fprintf(fp, "\n");
    }

    fclose(fp);
  }
}


