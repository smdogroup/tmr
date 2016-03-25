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

  // Compute the number of octants along each edge for
  // the leaves
  int32_t nx = 1 << refine_level;

  // Set the number of octants created
  int nocts = nx*nx*nx;

  // Compute the refinement level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - refine_level);

  // Create an array of the octants that will be stored
  TMROctant *array = new TMROctant[ nocts ];

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

  // Zero everything else
  order = 2;
  num_elements = 0;
  num_nodes = 0;
  num_dependent_nodes = 0;

  // Zero the mesh connectivity data
  elem_ptr = NULL;
  elem_conn = NULL;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;
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

  // Zero everything else
  order = 2;
  num_elements = 0;
  num_nodes = 0;
  num_dependent_nodes = 0;

  // Zero the mesh connectivity data
  elem_ptr = NULL;
  elem_conn = NULL;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;
}

/*
  Create the octree tree from the array
*/
TMROctree::TMROctree( TMROctantArray *_elements ){
  elements = _elements;
  elements->sort();

  // Zero the array of octant nodes
  nodes = NULL;

  // Zero everything else
  order = 2;
  num_elements = 0;
  num_nodes = 0;
  num_dependent_nodes = 0;

  // Zero the mesh connectivity data
  elem_ptr = NULL;
  elem_conn = NULL;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;
}

/*
  Free the octree
*/
TMROctree::~TMROctree(){
  // Free the elements and nodes
  if (elements){ delete elements; }
  if (nodes){ delete nodes; }

  // Free the connectivity information if it was allocated
  if (elem_ptr){ delete [] elem_ptr; }
  if (elem_conn){ delete [] elem_conn; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
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

  // Create a hash table for the refined  t
  TMROctantHash *hash = new TMROctantHash();

  // Get the current array of octants
  int size;
  TMROctant *array;
  elements->getArray(&array, &size);

  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Add the element
  for ( int i = 0; i < size; i++ ){
    // Try adding all of the children
    TMROctant q;
    array[i].getSibling(0, &q);

    if (refinement[i] == 0){
      hash->addOctant(&q);
    }
    else if (refinement[i] < 0){
      if (q.level > min_level){
	q.level = q.level-1;
	hash->addOctant(&q);
      }
      else {
	hash->addOctant(&q);
      }
    }
    else if (refinement[i] > 0){
      if (q.level < max_level){
	TMROctant c;
	c.level = q.level + 1;
	const int32_t h = 1 << (TMR_MAX_LEVEL - c.level);
	
	for ( int kk = 0; kk < 2; kk++ ){
	  for ( int jj = 0; jj < 2; jj++ ){
	    for ( int ii = 0; ii < 2; ii++ ){
	      c.x = q.x + 2*h*ii;
	      c.y = q.y + 2*h*jj;
	      c.z = q.z + 2*h*kk;
	      hash->addOctant(&c);
	    }
	  }
	}
      }
      else {
	hash->addOctant(&q);
      }
    }
  }

  // Now convert the new hash
  TMROctantArray *child0_elems = hash->toArray();
  child0_elems->getArray(&array, &size);

  // Loop over all elements and add their siblings
  for ( int i = 0; i < size; i++ ){
    for ( int j = 0; j < 8; j++ ){
      TMROctant q;
      array[i].getSibling(j, &q);
      if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
          q.x < hmax && q.y < hmax && q.z < hmax){
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
  Balance the octree in place

  This algorithm uses a hash and a queue to balance the octree. For
  each element in the octree, we add the neighbors that are required
  to balance to the tree. If the element is not in the hash, we add
  them to a queue, which keeps track of recently added elements. After
  the first pass, the algorithm continues popping elements until the
  queue is empty.

  Note that only 0-th siblings are added/popped on the hash/queue.
  Then at the end, all neighboring siblings are added.

  The type of balancing - face/edge balanced or face/edge/corner
  balanced is determined using the balance_corner flag. Face balancing
  is balancing across faces, corner balances across corners of the
  elements and corner balances across corners. The code always
  balances faces and edges (so that there is at most one depdent node
  per edge) and balances across corners optionally.
*/
void TMROctree::balance( int balance_corner ){
  // Create a hash table for the balanced tree
  TMROctantHash *hash = new TMROctantHash();

  // Go through the existing list of octants and add up everything for
  // balancing
  TMROctant p, q, neighbor;

  // Get the current array of octants
  int size;
  TMROctant *array;
  elements->getArray(&array, &size);

  // Create the queue of octants
  TMROctantQueue *queue = new TMROctantQueue();

  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Add the element
  for ( int i = 0; i < size; i++ ){
    // Try adding all of the children
    array[i].getSibling(0, &q);
    hash->addOctant(&q);
    
    // Get the parent of the octant, and add the their
    // face-matched octants from each face, as long 
    // as they fall within the bounds
    if (q.level > 0){
      q.parent(&p);

      // For each face, get the neighbours along
      // each face
      for ( int face = 0; face < 6; face++ ){
	p.faceNeighbor(face, &neighbor);
	neighbor.getSibling(0, &q);

	// If we're in bounds, add the neighbor
	if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	    q.x < hmax && q.y < hmax && q.z < hmax){
	  if (hash->addOctant(&q)){
	    queue->push(&q);
	  }
	}
      }

      // If we are balancing across edges, also 
      // add the edge-adjacent elements
      for ( int edge = 0; edge < 12; edge++ ){
        p.edgeNeighbor(edge, &neighbor);
        neighbor.getSibling(0, &q);
	  
        // If we're in bounds, add the neighbor
        if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
            q.x < hmax && q.y < hmax && q.z < hmax){
          if (hash->addOctant(&q)){
            queue->push(&q);
          }
        }
      }

      // If we're balancing across edges and 
      if (balance_corner){
	for ( int corner = 0; corner < 8; corner++ ){
	  p.cornerNeighbor(corner, &neighbor);
	  neighbor.getSibling(0, &q);
	  
	  // If we're in bounds, add the neighbor
	  if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	      q.x < hmax && q.y < hmax && q.z < hmax){
	    if (hash->addOctant(&q)){
	      queue->push(&q);
	    }
	  }
	}
      }
    }
  }

  // Now continue until the queue of added octants is
  // empty. At each iteration, pop an octant and add 
  // its neighbours until nothing new is added. This code
  // handles the propagation of octants to adjacent octants.
  while (queue->length() > 0){
    q = queue->pop();
   
    if (q.level > 1){
      q.parent(&p);

      // Add the octants across faces
      for ( int face = 0; face < 6; face++ ){
	p.faceNeighbor(face, &neighbor);
	neighbor.getSibling(0, &q);
	
	// If we're in bounds, add the neighbor
	if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	    q.x < hmax && q.y < hmax && q.z < hmax){
	  if (hash->addOctant(&q)){
	    queue->push(&q);
	  }
	}
      }

      // Add the octants across adjacent edges
      for ( int edge = 0; edge < 12; edge++ ){
        p.edgeNeighbor(edge, &neighbor);
        neighbor.getSibling(0, &q);
	
        // If we're in bounds, add the neighbor
        if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
            q.x < hmax && q.y < hmax && q.z < hmax){
          if (hash->addOctant(&q)){
            queue->push(&q);
          }
        }
      }

      // Add the octants across corners
      if (balance_corner){
	for ( int corner = 0; corner < 8; corner++ ){
	  p.cornerNeighbor(corner, &neighbor);
	  neighbor.getSibling(0, &q);
	  
	  // If we're in bounds, add the neighbor
	  if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	      q.x < hmax && q.y < hmax && q.z < hmax){
	    if (hash->addOctant(&q)){
	      queue->push(&q);
	    }
	  }
	}
      }
    }
  }

  // Now convert the new hash
  TMROctantArray *child0_elems = hash->toArray();
  child0_elems->getArray(&array, &size);

  // Loop over all elements and add their siblings
  for ( int i = 0; i < size; i++ ){
    for ( int j = 0; j < 8; j++ ){
      array[i].getSibling(j, &q);
      if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
          q.x < hmax && q.y < hmax && q.z < hmax){
        hash->addOctant(&q);
      }
    }
  }

  // Free the temporary elements
  delete child0_elems;

  // Free the old elements class
  delete elements;

  // Cover the hash table to a list and uniquely sort it
  elements = hash->toArray();
  elements->sort();

  // Free the hash table and the queue
  delete hash;
  delete queue;
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
  *high = num_elements;

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
  Retrieve the mesh information
*/
void TMROctree::getMesh( int *_num_nodes, 
                         int *_num_elements, 
                         const int **_elem_ptr, 
                         const int **_elem_conn ){
  if (_num_nodes){ *_num_nodes = num_nodes; }
  if (_num_elements){ *_num_elements = num_elements; }
  if (_elem_ptr){ *_elem_ptr = elem_ptr; }
  if (_elem_conn){ *_elem_conn = elem_conn; }
}

/*
  Retrieve the depednent node information
*/
void TMROctree::getDependentMesh( int *_num_dep_nodes, 
                                  const int **_dep_ptr,
                                  const int **_dep_conn,
                                  const double **_dep_weights ){
  if (_num_dep_nodes){ *_num_dep_nodes = num_dependent_nodes; }
  if (_dep_ptr){ *_dep_ptr = dep_ptr; }
  if (_dep_conn){ *_dep_conn = dep_conn; }
  if (_dep_weights){ *_dep_weights = dep_weights; }
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

  // Build the information for the dependent nodes that
  // lie along an element edge
  TMROctantQueue *edge_queue = new TMROctantQueue();
  TMROctantQueue *edge_nodes = new TMROctantQueue();

  // Build the information for the dependent nodes that 
  // lie on an element face
  TMROctantQueue *face_queue = new TMROctantQueue();
  TMROctantQueue *face_nodes = new TMROctantQueue();
  
  if (order == 2){
    // For all operations here, we compare the node numbers
    const int use_nodes = 1;

    // Now check for dependent nodes on each edge and face
    for ( int i = 0; i < size; i++ ){
      // The following indices are used to define the
      // unit vectors along the edges/faces of an element
      const int iindex[] = {1, 0, 0};
      const int jindex[] = {2, 2, 1};

      // Get the side length of the element
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);      
      
      // Get the size of the next finest level
      const int32_t hc = 1 << (TMR_MAX_LEVEL - array[i].level - 1);
      
      // Check if this is a dependent edge
      for ( int k = 0; k < 3; k++ ){
	// The constant and variable unit vectors for
	// determining the edge location
	int ec[3] = {0, 0, 0};
	int ie[3] = {0, 0, 0};
	int je[3] = {0, 0, 0};
	
	// Set the indices of the unit vectors
	ec[k] = 1;
	ie[iindex[k]] = 1;
	je[jindex[k]] = 1;

	for ( int jj = 0; jj < 2; jj++ ){
	  for ( int ii = 0; ii < 2; ii++ ){
	    // Check the location for the dependent node
	    TMROctant c;
	    c.x = array[i].x + hc*ec[0] + h*(ii*ie[0] + jj*je[0]);
	    c.y = array[i].y + hc*ec[1] + h*(ii*ie[1] + jj*je[1]);
	    c.z = array[i].z + hc*ec[2] + h*(ii*ie[2] + jj*je[2]);

	    TMROctant *t;
	    if (t = nodes->contains(&c, use_nodes)){
	      if (t->tag == 1){
		// Push the dependent node from the edge
		edge_queue->push(&c);
            
		// Create and push the independent nodes from the edge
		TMROctant n;
		n.x = array[i].x + h*(ii*ie[0] + jj*je[0]);
		n.y = array[i].y + h*(ii*ie[1] + jj*je[1]);
		n.z = array[i].z + h*(ii*ie[2] + jj*je[2]);
		edge_nodes->push(&n);
             
		// Set the other dependent node
		n.x += h*ec[0];
		n.y += h*ec[1];
		n.z += h*ec[2];
		edge_nodes->push(&n);

		// Set the node so it does not get added again
		t->tag = -1;
	      }
	    }
	  }
	}
      }

      // Check if this is a dependent face
      for ( int k = 0; k < 3; k++ ){
	// The constant and variable unit vectors for
	// determining the edge location
	int ec[3] = {0, 0, 0};
	int ie[3] = {0, 0, 0};
	int je[3] = {0, 0, 0};

	// Set the indices of the unit vectors
	ec[k] = 1;
	ie[iindex[k]] = 1;
	je[jindex[k]] = 1;

	// Scan over each face
	for ( int ii = 0; ii < 2; ii++ ){
	  TMROctant c;
	  c.x = array[i].x + h*ii*ec[0] + hc*(ie[0] + je[0]);
	  c.y = array[i].y + h*ii*ec[1] + hc*(ie[1] + je[1]);
	  c.z = array[i].z + h*ii*ec[2] + hc*(ie[2] + je[2]);

	  TMROctant *t;
	  if (t = nodes->contains(&c, use_nodes)){
	    if (t->tag == 1){
	      // Push the dependent node from the face
	      face_queue->push(&c);

	      // Push the independent nodes from the face 
	      for ( int jface = 0; jface < 2; jface++ ){
		for ( int iface = 0; iface < 2; iface++ ){
		  TMROctant n;
		  n.x = array[i].x + h*ii*ec[0] + h*(iface*ie[0] + jface*je[0]);
		  n.y = array[i].y + h*ii*ec[1] + h*(iface*ie[1] + jface*je[1]);
		  n.z = array[i].z + h*ii*ec[2] + h*(iface*ie[2] + jface*je[2]);
		  face_nodes->push(&n);
		}
	      }
	      t->tag = -2;
	    }
	  }
	}
      }     
    }
  }
  else if (order == 3){}

  // Get the current array of octants
  int node_size;
  TMROctant *node_array;
  nodes->getArray(&node_array, &node_size);

  // Order the elements
  num_elements = 0;
  for ( int i = 0; i < size; i++ ){
    array[i].tag = num_elements;
    num_elements++;
  }

  // Count up the number of independent and dependent
  // node numbers in the array
  num_nodes = 0;
  num_dependent_nodes = 0;

  // Now reset the tags as the node numbers
  for ( int i = 0; i < node_size; i++ ){
    if (node_array[i].tag >= 0){
      node_array[i].tag = num_nodes;
      num_nodes++;
    }
    else {
      node_array[i].tag = -(num_dependent_nodes+1);
      num_dependent_nodes++;
    }
  }

  // Covert the queues to arrays
  TMROctantArray *edge_array = edge_queue->toArray();
  TMROctantArray *face_array = face_queue->toArray();

  // Free the queue objects
  delete edge_queue;
  delete face_queue;

  // Get the arrays of edges and faces
  int nedges, nfaces;
  TMROctant *edges, *faces;
  edge_array->getArray(&edges, &nedges);
  face_array->getArray(&faces, &nfaces);

  // Free the dependent node information if it already exists
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
  
  // Allocate the arrays for the dependent nodes
  dep_ptr = new int[ num_dependent_nodes+1 ];
  dep_conn = new int[ order*nedges + order*order*nfaces ];
  dep_weights = new double[ order*nedges + order*order*nfaces ];

  // Find the dependent variables and check them
  memset(dep_ptr, 0, (num_dependent_nodes+1)*sizeof(int));

  // Add the counts from the edge and face dependent node relationships
  for ( int i = 0; i < nedges; i++ ){
    const int use_nodes = 1;
    TMROctant *t = nodes->contains(&edges[i], use_nodes);
    int node = -t->tag-1;
    dep_ptr[node+1] = 2;
  }
  for ( int i = 0; i < nfaces; i++ ){
    const int use_nodes = 1;
    TMROctant *t = nodes->contains(&faces[i], use_nodes);
    int node = -t->tag-1;
    dep_ptr[node+1] = 4;
  }

  // Count up the totals so that dep_ptr points into the
  // correct location in the dep_conn and dep_weight arrays
  for ( int i = 0; i < num_dependent_nodes; i++ ){
    dep_ptr[i+1] += dep_ptr[i];
  }

  // Loop over all the independent nodes in the queues in the
  // same order as before
  for ( int i = 0; i < nedges; i++ ){
    // Get the dependent node number
    const int use_nodes = 1;
    TMROctant *t = nodes->contains(&edges[i], use_nodes);
    int node = -t->tag-1;

    // Retrieve the independent nodes
    for ( int k = 0; k < order; k++ ){
      TMROctant n = edge_nodes->pop();
      TMROctant *t = nodes->contains(&n, use_nodes);
      dep_conn[dep_ptr[node]+k] = t->tag;
      dep_weights[dep_ptr[node]+k] = 0.5;
    }
  }

  // Loop over all the independent nodes in the faces
  for ( int i = 0; i < nfaces; i++ ){    
    // Get the dependent node number
    const int use_nodes = 1;
    TMROctant *t = nodes->contains(&faces[i], use_nodes);
    int node = -t->tag-1;

    // Retrieve the independent nodes
    for ( int k = 0; k < order*order; k++ ){
      TMROctant n = face_nodes->pop();
      TMROctant *t = nodes->contains(&n, use_nodes);
      dep_conn[dep_ptr[node]+k] = t->tag;
      dep_weights[dep_ptr[node]+k] = 0.25;
    }
  }

  // Free the remaining data
  delete edge_array;
  delete edge_nodes;
  delete face_array;
  delete face_nodes;
}

/*
  Create the mesh connectivity
*/
void TMROctree::createMesh( int _order ){
  // Scan through and create the connectivity for all the
  // elements in the Morton order
  if (!nodes){
    if (_order < 2){ _order = 2; }
    if (_order > 3){ _order = 3; }
    
    // Create the mesh and order the nodes
    createNodes(_order);
  }

  // Get the current array of octants
  int size;
  TMROctant *array;
  elements->getArray(&array, &size);

  // We already know how large the connectivity list
  // will be and the size of the number of dependent nodes
  if (elem_ptr){ delete [] elem_ptr; }
  if (elem_conn){ delete [] elem_conn; }
  elem_ptr = new int[ num_elements+1 ];
  elem_conn = new int[ order*order*order*num_elements ];
  
  // Set the connectivity pointer
  int conn_size = 0;
  elem_ptr[0] = 0;
  for ( int i = 0; i < size; i++ ){
    TMROctant p;

    if (order == 2){
      // For all searches/comparisons, we use node numbers
      const int use_nodes = 1;

      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      p.level = array[i].level;

      for ( int kk = 0; kk < 2; kk++ ){
        for ( int jj = 0; jj < 2; jj++ ){
          for ( int ii = 0; ii < 2; ii++ ){
            p.x = array[i].x + ii*h;
            p.y = array[i].y + jj*h;
            p.z = array[i].z + kk*h;

            // Get the index for the node
            TMROctant *t = nodes->contains(&p, use_nodes);
            
            // Reset the node level if the element level is smaller
            if (array[i].level > t->level){
              t->level = array[i].level;
            }

            // Set the node number in the element connectivity
            elem_conn[conn_size] = t->tag;
            conn_size++;
          }
        }
      }
    }
    else if (order == 3){
      // For all searches, use nodes
      const int use_nodes = 1;

      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level - 1);
      p.level = array[i].level+1;

      // Add all of the nodes to the hash
      for ( int kk = 0; kk < 3; kk++ ){
        for ( int jj = 0; jj < 3; jj++ ){
          for ( int ii = 0; ii < 3; ii++ ){
            p.x = array[i].x + ii*h;
            p.y = array[i].y + jj*h;
            p.z = array[i].z + kk*h;

            // Get the index for the node
            TMROctant *t = nodes->contains(&p, use_nodes);

            // Reset the node level if the element level is smaller
            if (array[i].level > t->level){
              t->level = array[i].level;
            }

            // Set the node number in the element connectivity
            elem_conn[conn_size] = t->tag;
            conn_size++;
          }
        }
      }
    }

    // Set the connectivity pointer
    elem_ptr[i+1] = conn_size;
  }
}

/*
  Create the interpolation operator
*/
void TMROctree::createInterpolation( TMROctree *coarse,
                                     int **_interp_ptr,
                                     int **_interp_conn,
                                     double **_interp_weights ){
  // Set the pointers to NULL
  *_interp_ptr = NULL;
  *_interp_conn = NULL;
  *_interp_weights = NULL;

  // If the nodes have not been allocated, then this function should
  // not be called. Return no result.
  if (!nodes || !coarse->nodes){
    return;
  }

  // Get the size of the fine node array
  int fine_size;
  TMROctant *fine;
  nodes->getArray(&fine, &fine_size);
  
  // Allocate the arrays
  int conn_len = 0;
  int max_conn_len = order*order*order*num_nodes;
  int *interp_ptr = new int[ num_nodes+1 ];
  int *interp_conn = new int[ max_conn_len ];
  double *interp_weights = new double[ max_conn_len ];

  // Allocate the data that will be used
  int nnodes = 0;
  interp_ptr[0] = 0;

  // The maximum possible size of the array of weights. Note
  // that this is found if every node is a dependent node (which is
  // impossible) which points to a dependent face node (also
  // impossible). It is an upper bound.
  int max_size = (order*order*order)*(order*order);
  TMRIndexWeight *weights = new TMRIndexWeight[ max_size ];

  // Set pointers to the coarse dependent data
  const int *cdep_ptr = coarse->dep_ptr;
  const int *cdep_conn = coarse->dep_conn;
  const double *cdep_weights = coarse->dep_weights;

  for ( int i = 0; i < fine_size; i++ ){
    // Use a node-based search
    const int use_nodes = 1;

    // Keep track of the number of new weight elements
    int nweights = 0;

    // Loop over all the adjacent nodes
    if (fine[i].tag >= 0){
      TMROctant *t = coarse->nodes->contains(&fine[i], use_nodes);
      if (t){
        if (t->tag >= 0){
          weights[nweights].index = t->tag;
          weights[nweights].weight = 1.0;
          nweights++;
        }
        else {
          // Un-ravel the dependent node connectivity
          int node = -t->tag-1;
          for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
            weights[nweights].index = cdep_conn[jp];
            weights[nweights].weight = cdep_weights[jp];
            nweights++;
          }
        }
      }
      else {
        // Find the child-id of the node. The node uses the
        // local length scale of the element
        int id = fine[i].childId();

        // Set the element level
        const int32_t h = 1 << (TMR_MAX_LEVEL - fine[i].level);

        if (id == 1 || id == 2 || id == 4){
          // Get the root sibling
          TMROctant n;
          fine[i].getSibling(0, &n);

          // Get the node number of the 0-th sibling
          t = coarse->nodes->contains(&n, use_nodes);
          if (t->tag >= 0){
            weights[nweights].index = t->tag;
            weights[nweights].weight = 0.5;
            nweights++;
          }
          else {
            int node = -t->tag-1;
            for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
              weights[nweights].index = cdep_conn[jp];
              weights[nweights].weight = 0.5*cdep_weights[jp];
              nweights++;
            }
          }

          if (id == 1){
            n.x = n.x + 2*h;
          }
          else if (id == 2){
            n.y = n.y + 2*h;
          }
          else if (id == 4){
            n.z = n.z + 2*h;
          }

          // Get the node number of the other node 
          t = coarse->nodes->contains(&n, use_nodes);
          if (t->tag >= 0){
            weights[nweights].index = t->tag;
            weights[nweights].weight = 0.5;
            nweights++;
          }
          else {
            int node = -t->tag-1;
            for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
              weights[nweights].index = cdep_conn[jp];
              weights[nweights].weight = 0.5*cdep_weights[jp];
              nweights++;
            }
          }
        }
        else if (id == 3 || id == 5 || id == 6){
          // Get the root sibling
          TMROctant n;
          fine[i].getSibling(0, &n);

          int ie[] = {0, 0, 0};
          int je[] = {0, 0, 0};
          if (id == 3){
            ie[0] = 1;  je[1] = 1;
          }
          else if (id == 5){
            ie[0] = 1;  je[2] = 1;
          }
          else if (id == 6){
            ie[1] = 1;  je[2] = 1;
          }

          for ( int jj = 0; jj < 2; jj++ ){
            for ( int ii = 0; ii < 2; ii++ ){
              TMROctant p;
              p.x = n.x + 2*h*(ii*ie[0] + jj*je[0]);
              p.y = n.y + 2*h*(ii*ie[1] + jj*je[1]);
              p.z = n.z + 2*h*(ii*ie[2] + jj*je[2]);

              // Get the node number of the other node 
              t = coarse->nodes->contains(&p, use_nodes);
              if (t->tag >= 0){
                weights[nweights].index = t->tag;
                weights[nweights].weight = 0.25;
                nweights++;
              }
              else {
                int node = -t->tag-1;
                for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
                  weights[nweights].index = cdep_conn[jp];
                  weights[nweights].weight = 0.25*cdep_weights[jp];
                  nweights++;
                }
              }
            }
          }
        }
        else if (id == 7){
          // Get the root sibling
          TMROctant n;
          fine[i].getSibling(0, &n);

          // Add the remaining nodes
          for ( int kk = 0; kk < 2; kk++ ){
            for ( int jj = 0; jj < 2; jj++ ){
              for ( int ii = 0; ii < 2; ii++ ){
                TMROctant p;
                p.x = n.x + 2*h*ii;
                p.y = n.y + 2*h*jj;
                p.z = n.z + 2*h*kk;

                // Get the node number of the other node 
                t = coarse->nodes->contains(&p, use_nodes);
                if (t->tag >= 0){
                  weights[nweights].index = t->tag;
                  weights[nweights].weight = 0.125;
                  nweights++;
                }
                else {
                  int node = -t->tag-1;
                  for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
                    weights[nweights].index = cdep_conn[jp];
                    weights[nweights].weight = 0.125*cdep_weights[jp];
                    nweights++;
                  }
                }
              }
            }
          }
        }
      }

      // Sort the dependent weight values
      nweights = TMRIndexWeight::uniqueSort(weights, nweights);
      
      // Check whether adding these will exceed the size
      // of the array
      if (interp_ptr[nnodes] + nweights >= max_conn_len){
        int estimate = ((2.0*interp_ptr[nnodes]/nnodes)*(num_nodes - nnodes));
        max_conn_len += estimate;
        int *temp = new int[ max_conn_len ];
        memcpy(temp, interp_conn, conn_len*sizeof(int));
        delete [] interp_conn;
        interp_conn = temp;
        
        double *tempw = new double[ max_conn_len ];
        memcpy(tempw, interp_weights, conn_len*sizeof(double));
        delete [] interp_weights;
        interp_weights = tempw;
      }

      // Extract the weights from the sorted list
      for ( int k = 0; k < nweights; k++ ){
        interp_conn[conn_len] = weights[k].index;
        interp_weights[conn_len] = weights[k].weight;
        conn_len++;
      }

      // Increment the pointer
      interp_ptr[nnodes+1] = conn_len;
      nnodes++;
    }
  }

  // Free the weights
  delete [] weights;

  // Set the return arrays
  *_interp_ptr = interp_ptr;
  *_interp_conn = interp_conn;
  *_interp_weights = interp_weights;
}

/*
  Create the restriction operator
*/
void TMROctree::createRestriction( TMROctree *tree,
                                   int **_interp_ptr,
                                   int **_interp_conn,
                                   double **_interp_weights ){
  // Set the pointers to NULL
  *_interp_ptr = NULL;
  *_interp_conn = NULL;
  *_interp_weights = NULL;

  // If the nodes have not been allocated, then this function should
  // not be called. Return no result.
  if (!nodes || !tree->nodes){
    return;
  }

  // Get the array of coarse nodes and the size of nodes
  int coarse_size;
  TMROctant *coarse;
  tree->nodes->getArray(&coarse, &coarse_size);
 
  // Set the weights for the full-approximation
  const double wvals[] = {0.5, 1.0, 0.5};

  // Set the number of independent coarse nodes
  int num_coarse_nodes = tree->num_nodes;

  // Allocate the arrays
  int conn_len = 0;
  int max_conn_len = order*order*order*num_coarse_nodes;
  int *interp_ptr = new int[ num_coarse_nodes+1 ];
  int *interp_conn = new int[ max_conn_len ];
  double *interp_weights = new double[ max_conn_len ];

  // Set the initial pointer values
  int nnodes = 0;
  interp_ptr[0] = 0;

  // The maximum possible size of the array of weights. Note
  // that this is found if every node is a dependent node (which is
  // impossible) which points to a dependent face node (also
  // impossible). It is an upper bound.
  int max_size = 27*(order*order);
  TMRIndexWeight *weights = new TMRIndexWeight[ max_size ];

  // Scan through the nodes of the coarse tree and
  // find the next-to-adjacent nodes in the octree
  for ( int i = 0; i < coarse_size; i++ ){
    if (coarse[i].tag >= 0){
      // Use node-based searching
      const int use_nodes = 1;

      // Keep track of the number of weights used
      int nweights = 0;
      double w = 0.0;

      // Get the fine index number corresponding to the
      // coarse index on this mesh
      TMROctant *fine;
      fine = nodes->contains(&coarse[i], use_nodes);

      // Get the coarse node level
      const int32_t h = 1 << (TMR_MAX_LEVEL - fine->level);
      
      // Scan through the node locations for the coarse mesh
      for ( int kk = 0; kk < 3; kk++ ){
        for ( int jj = 0; jj < 3; jj++ ){
          for ( int ii = 0; ii < 3; ii++ ){
            // Set the location for the adjacent nodes
            TMROctant n;
            n.x = coarse[i].x + h*(ii-1);
            n.y = coarse[i].y + h*(jj-1);
            n.z = coarse[i].z + h*(kk-1);

            // Get the element level
            TMROctant *t;
            t = nodes->contains(&n, use_nodes);
            if (t){
              double wk = wvals[ii]*wvals[jj]*wvals[kk];
              w += wk;

              if (t->tag >= 0){
                weights[nweights].index = t->tag;
                weights[nweights].weight = wk;
                nweights++;
              }
              else {
                int node = -t->tag-1;
                for ( int jp = dep_ptr[node]; jp < dep_ptr[node+1]; jp++ ){
                  weights[nweights].index = dep_conn[jp];
                  weights[nweights].weight = wk*dep_weights[jp];
                  nweights++;
                }
              }
            }
          }
        }
      }

      // Sort the dependent weight values
      nweights = TMRIndexWeight::uniqueSort(weights, nweights);
      
      // Check whether adding these will exceed the size
      // of the array
      if (interp_ptr[nnodes] + nweights >= max_conn_len){
        int estimate = ((2.0*interp_ptr[nnodes]/nnodes)*(num_nodes - nnodes));
        max_conn_len += estimate;
        int *temp = new int[ max_conn_len ];
        memcpy(temp, interp_conn, conn_len*sizeof(int));
        delete [] interp_conn;
        interp_conn = temp;
        
        double *tempw = new double[ max_conn_len ];
        memcpy(tempw, interp_weights, conn_len*sizeof(double));
        delete [] interp_weights;
        interp_weights = tempw;
      }

      // Extract the weights from the sorted list. Normalize
      // the weights by the total weights added at this node
      for ( int k = 0; k < nweights; k++ ){
        interp_conn[conn_len] = weights[k].index;
        interp_weights[conn_len] = weights[k].weight/w;
        conn_len++;
      }

      // Increment the pointer
      interp_ptr[nnodes+1] = conn_len;
      nnodes++;
    }
  }

  // Free the weight used
  delete [] weights;

  // Set the return arrays
  *_interp_ptr = interp_ptr;
  *_interp_conn = interp_conn;
  *_interp_weights = interp_weights;
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


