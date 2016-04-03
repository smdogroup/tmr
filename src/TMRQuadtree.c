#include "TMRQuadtree.h"

/*
  Refine the initial quadtree to the specified depth along all
  coordinate directions 
*/
TMRQuadtree::TMRQuadtree( int refine_level, 
                          TMRQuadrant *_domain, int _ndomain ){
  // Set the domain level
  ndomain = _ndomain;
  domain = NULL;
  if (ndomain > 0){
    domain = new TMRQuadrant[ ndomain ];
    memcpy(domain, _domain, ndomain*sizeof(TMRQuadrant));
  }

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

  // Keep track of the quadrants and the length of the quadrant array
  int nquads = 0;
  TMRQuadrant *array = NULL;
  
  if (domain){
    // Loop over each quadrant in the domain
    for ( int k = 0; k < ndomain; k++ ){
      int d = refine_level - domain[k].level;
      if (d >= 0){
        int nx = 1 << d;
        nquads += nx*nx;
      }
    }

    // Create an array of the quadrants that will be stored
    array = new TMRQuadrant[ nquads ];

    // Create all the quadrants required
    int index = 0;
    for ( int k = 0; k < ndomain; k++ ){
      int d = refine_level - domain[k].level;
      if (d >= 0){
        const int32_t hk = 1 << (TMR_MAX_LEVEL - domain[k].level);
        const int32_t xmax = domain[k].x + hk;
        const int32_t ymax = domain[k].y + hk;

        for ( int32_t y = domain[k].y; y < ymax; y += h ){
          for ( int32_t x = domain[k].x; x < xmax; x += h ){
            array[index].x = x;
            array[index].y = y;
            array[index].level = refine_level;
            index++;
          }
        }
      }
    }
  }
  else {
    // Compute the number of quadrants along each edge
    int32_t nx = 1 << refine_level;
    
    // Set the number of quadrants created
    nquads = nx*nx;

    // Create an array of the quadrants that will be stored
    array = new TMRQuadrant[ nquads ];

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
  }

  // Sort the array of quadrants
  elements = new TMRQuadrantArray(array, nquads);
  elements->sort();

  // Zero the array of quadrant nodes
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
  Generate a random initial quadtree that can be used for testing.
*/
TMRQuadtree::TMRQuadtree( int nrand, int min_level, int max_level ){
  // Create an array of the quadrants that will be stored
  TMRQuadrant *array = new TMRQuadrant[ nrand ];

  // Generate a random number of quadrants along random directions
  for ( int i = 0; i < nrand; i++ ){
    int32_t level = min_level + (rand() % (max_level - min_level + 1));

    const int32_t h = 1 << (TMR_MAX_LEVEL - level);
    int32_t x = h*(rand() % (1 << level));
    int32_t y = h*(rand() % (1 << level));

    array[i].x = x;
    array[i].y = y;
    array[i].level = level;
  }

  elements = new TMRQuadrantArray(array, nrand);
  elements->sort();

  // Zero the array of quadrant nodes
  nodes = NULL;

  // Zero the domain
  ndomain = 0;
  domain = NULL;

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
  Create the quadtree tree from the array
*/
TMRQuadtree::TMRQuadtree( TMRQuadrantArray *_elements,
                          TMRQuadrant *_domain, int _ndomain ){
  elements = _elements;
  elements->sort();

  // Set the domain level
  ndomain = _ndomain;
  domain = NULL;
  if (ndomain > 0){
    domain = new TMRQuadrant[ ndomain ];
    memcpy(domain, _domain, ndomain*sizeof(TMRQuadrant));
  }

  // Zero the array of quadrant nodes
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
  Free the quadtree
*/
TMRQuadtree::~TMRQuadtree(){
  // Free the elements and nodes
  if (elements){ delete elements; }
  if (nodes){ delete nodes; }

  // Free the domain if it exists
  if (domain){ delete [] domain; }

  // Free the connectivity information if it was allocated
  if (elem_ptr){ delete [] elem_ptr; }
  if (elem_conn){ delete [] elem_conn; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
}

/*
  Check that the provided element quadrant is within the domain 
*/
int TMRQuadtree::inDomain( TMRQuadrant *p ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // If no domain is specified, use the unit cube
  // as the full domain
  if (!domain){
    if (p->x >= 0 && p->y >= 0 &&
        p->x < hmax && p->y < hmax){
      return 1;
    }
    else {
      return 0;
    }
  }
  else {
    // If a domain is specified, then search each element within
    // the domain. Note that the domain list should be small.
    for ( int k = 0; k < ndomain; k++ ){
      const int h = 1 << (TMR_MAX_LEVEL - domain[k].level);
      const int x1 = domain[k].x;
      const int y1 = domain[k].y;
      const int x2 = domain[k].x + h;
      const int y2 = domain[k].y + h;
      
      if (p->x >= x1 && p->y >= y1 &&
          p->x < x2 && p->y < y2){
        return 1;
      }
    }
  }

  return 0;
}

/*
  Check if the provided quadrant is on the domain boundary
*/
int TMRQuadtree::onBoundary( TMRQuadrant *p ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // If no domain is specified, use the unit cube
  // as the full domain and check whether the node
  // provided is on the boundary
  if (!domain){
    if ((p->x == 0 || p->x == hmax) ||
        (p->y == 0 || p->y == hmax)){
      return 1;
    }
    else {
      return 0;
    }
  }
  else {
    // If a domain is specified, then search each element within
    // the domain and only return true if the node is on one
    // and only one boundary.
    int nbound = 0;
    int corner = 0; // The node is on a corner
    int edge = 0; // The node is on an edge
    for ( int k = 0; k < ndomain; k++ ){
      const int32_t h = 1 << (TMR_MAX_LEVEL - domain[k].level);
      const int32_t x1 = domain[k].x;
      const int32_t y1 = domain[k].y;
      const int32_t x2 = domain[k].x + h;
      const int32_t y2 = domain[k].y + h;
      
      // Compute the range
      int rx = ((p->x >= x1) && (p->x <= x2));
      int ry = ((p->y >= y1) && (p->y <= y2));
      
      // Flag true if any of the nodes are on the bounds
      int ax = ((p->x == x1) || (p->x == x2)) && ry;
      int ay = ((p->y == y1) || (p->y == y2)) && rx;

      // Check if this is on any of the bounds
      if (ax || ay){
        nbound++;
      }

      // Check if this node is on a corner
      corner = corner || (ax && ay);
    }

    // The node is on the true boundary only if it is only
    // on one of the quadrants that define the boundary
    if (nbound == 0){
      return 0;
    }
    else if (corner){
      return (nbound < 4);
    }
    else {
      return (nbound < 2);
    }
  }

  return 0;
}

/*
  Refine the quadtree by adding/subtracting elements

  This code takes in an array of integers that are either positive,
  negative or zero. A positive value means that the cell is refined,
  and all children are added. A negative index means that only the
  parent of the quadrant is refined. A zero index indicates that the
  quadrant is retained.
*/
void TMRQuadtree::refine( int refinement[], 
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
  TMRQuadrantHash *hash = new TMRQuadrantHash();

  // Get the current array of quadrants
  int size;
  TMRQuadrant *array;
  elements->getArray(&array, &size);

  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Add the element
  for ( int i = 0; i < size; i++ ){
    // Try adding all of the children
    TMRQuadrant q;
    array[i].getSibling(0, &q);

    if (refinement[i] == 0){
      hash->addQuadrant(&q);
    }
    else if (refinement[i] < 0){
      if (q.level > min_level){
	q.level = q.level-1;
	hash->addQuadrant(&q);
      }
      else {
	hash->addQuadrant(&q);
      }
    }
    else if (refinement[i] > 0){
      if (q.level < max_level){
	TMRQuadrant c;
	c.level = q.level + 1;
	const int32_t h = 1 << (TMR_MAX_LEVEL - c.level);
	
        for ( int jj = 0; jj < 2; jj++ ){
          for ( int ii = 0; ii < 2; ii++ ){
            c.x = q.x + 2*h*ii;
            c.y = q.y + 2*h*jj;
            hash->addQuadrant(&c);
          }
        }
      }
      else {
	hash->addQuadrant(&q);
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
      if (inDomain(&q)){
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
  Balance the quadtree in place

  This algorithm uses a hash and a queue to balance the quadtree. For
  each element in the quadtree, we add the neighbors that are required
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
void TMRQuadtree::balance( int balance_corner ){
  // Create a hash table for the balanced tree
  TMRQuadrantHash *hash = new TMRQuadrantHash();

  // Go through the existing list of quadrants and add up everything for
  // balancing
  TMRQuadrant p, q, neighbor;

  // Get the current array of quadrants
  int size;
  TMRQuadrant *array;
  elements->getArray(&array, &size);

  // Create the queue of quadrants
  TMRQuadrantQueue *queue = new TMRQuadrantQueue();

  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Add the element
  for ( int i = 0; i < size; i++ ){
    // Try adding all of the children
    array[i].getSibling(0, &q);
    hash->addQuadrant(&q);
    
    // Get the parent of the quadrant, and add the their
    // face-matched quadrants from each face, as long 
    // as they fall within the bounds
    if (q.level > 0){
      q.parent(&p);

      // If we are balancing across edges, also 
      // add the edge-adjacent elements
      for ( int edge = 0; edge < 4; edge++ ){
        p.edgeNeighbor(edge, &neighbor);
        neighbor.getSibling(0, &q);
	  
        // If we're in bounds, add the neighbor
	if (inDomain(&q)){
          if (hash->addQuadrant(&q)){
            queue->push(&q);
          }
        }
      }

      // If we're balancing across edges and 
      if (balance_corner){
	for ( int corner = 0; corner < 4; corner++ ){
	  p.cornerNeighbor(corner, &neighbor);
	  neighbor.getSibling(0, &q);
	  
	  // If we're in bounds, add the neighbor
          if (inDomain(&q)){
	    if (hash->addQuadrant(&q)){
	      queue->push(&q);
	    }
	  }
	}
      }
    }
  }

  // Now continue until the queue of added quadrants is
  // empty. At each iteration, pop an quadrant and add 
  // its neighbours until nothing new is added. This code
  // handles the propagation of quadrants to adjacent quadrants.
  while (queue->length() > 0){
    q = queue->pop();
   
    if (q.level > 1){
      q.parent(&p);

      // Add the quadrants across adjacent edges
      for ( int edge = 0; edge < 4; edge++ ){
        p.edgeNeighbor(edge, &neighbor);
        neighbor.getSibling(0, &q);
	
        // If we're in bounds, add the neighbor
        if (inDomain(&q)){
          if (hash->addQuadrant(&q)){
            queue->push(&q);
          }
        }
      }

      // Add the quadrants across corners
      if (balance_corner){
	for ( int corner = 0; corner < 4; corner++ ){
	  p.cornerNeighbor(corner, &neighbor);
	  neighbor.getSibling(0, &q);
	  
	  // If we're in bounds, add the neighbor
          if (inDomain(&q)){
	    if (hash->addQuadrant(&q)){
	      queue->push(&q);
	    }
	  }
	}
      }
    }
  }

  // Now convert the new hash
  TMRQuadrantArray *child0_elems = hash->toArray();
  child0_elems->getArray(&array, &size);

  // Loop over all elements and add their siblings
  for ( int i = 0; i < size; i++ ){
    for ( int j = 0; j < 4; j++ ){
      array[i].getSibling(j, &q);
      if (inDomain(&q)){
        hash->addQuadrant(&q);
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
  TMRQuadtree *tree = new TMRQuadtree(queue->toArray(), domain, ndomain);
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
  *low = 0;
  *high = num_elements;

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
  Retrieve the mesh information
*/
void TMRQuadtree::getMesh( int *_num_nodes, 
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
void TMRQuadtree::getDependentMesh( int *_num_dep_nodes, 
                                    const int **_dep_ptr,
                                    const int **_dep_conn,
                                    const double **_dep_weights ){
  if (_num_dep_nodes){ *_num_dep_nodes = num_dependent_nodes; }
  if (_dep_ptr){ *_dep_ptr = dep_ptr; }
  if (_dep_conn){ *_dep_conn = dep_conn; }
  if (_dep_weights){ *_dep_weights = dep_weights; }
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

  // Build the information for the dependent nodes that
  // lie along an element edge
  TMRQuadrantQueue *edge_queue = new TMRQuadrantQueue();
  TMRQuadrantQueue *edge_nodes = new TMRQuadrantQueue();
  
  if (order == 2){
    // For all operations here, we compare the node numbers
    const int use_nodes = 1;

    // Now check for dependent nodes on each edge and face
    for ( int i = 0; i < size; i++ ){
      // Get the side length of the element
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);      
      
      // Get the size of the next finest level
      const int32_t hc = 1 << (TMR_MAX_LEVEL - array[i].level - 1);

      // Check if this is a dependent edge
      for ( int k = 0; k < 2; k++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          // Check the location for the dependent node
          TMRQuadrant c;
          if (k == 0){
            c.x = array[i].x + ii*h;
            c.y = array[i].y + hc;
          }
          else {
            c.x = array[i].x + hc;
            c.y = array[i].y + ii*h;
          }

          TMRQuadrant *t;
          if (t = nodes->contains(&c, use_nodes)){
            if (t->tag == 1){
              // Push the dependent node from the edge
              edge_queue->push(&c);
            
              // Create and push the independent nodes from the edge
              TMRQuadrant n;
              if (k == 0){
                n.x = array[i].x + ii*h;
                n.y = array[i].y;
              }
              else {
                n.x = array[i].x;
                n.y = array[i].y + ii*h;
              }
              edge_nodes->push(&n);
             
              // Set the other dependent node
              if (k == 0){
                n.x = array[i].x + ii*h;
                n.y = array[i].y + h;
              }
              else {
                n.x = array[i].x + h;
                n.y = array[i].y + ii*h;
              }
              edge_nodes->push(&n);
              
              // Set the node so it does not get added again
              t->tag = -1;
            }
          }
        }
      }
    }
  }
  else if (order == 3){}

  // Get the current array of quadrants
  int node_size;
  TMRQuadrant *node_array;
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
  TMRQuadrantArray *edge_array = edge_queue->toArray();

  // Free the queue objects
  delete edge_queue;

  // Get the arrays of edges and faces
  int nedges;
  TMRQuadrant *edges;
  edge_array->getArray(&edges, &nedges);
  
  // Free the dependent node information if it already exists
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
  
  // Allocate the arrays for the dependent nodes
  dep_ptr = new int[ num_dependent_nodes+1 ];
  dep_conn = new int[ order*nedges ];
  dep_weights = new double[ order*nedges ];

  // Find the dependent variables and check them
  memset(dep_ptr, 0, (num_dependent_nodes+1)*sizeof(int));

  // Add the counts from the edge and face dependent node relationships
  for ( int i = 0; i < nedges; i++ ){
    const int use_nodes = 1;
    TMRQuadrant *t = nodes->contains(&edges[i], use_nodes);
    int node = -t->tag-1;
    dep_ptr[node+1] = 2;
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
    TMRQuadrant *t = nodes->contains(&edges[i], use_nodes);
    int node = -t->tag-1;

    // Retrieve the independent nodes
    for ( int k = 0; k < order; k++ ){
      TMRQuadrant n = edge_nodes->pop();
      TMRQuadrant *t = nodes->contains(&n, use_nodes);
      dep_conn[dep_ptr[node]+k] = t->tag;
      dep_weights[dep_ptr[node]+k] = 0.5;
    }
  }

  // Free the remaining data
  delete edge_array;
  delete edge_nodes;
}

/*
  Create the mesh connectivity
*/
void TMRQuadtree::createMesh( int _order ){
  // Scan through and create the connectivity for all the
  // elements in the Morton order
  if (!nodes){
    if (_order < 2){ _order = 2; }
    if (_order > 3){ _order = 3; }
    
    // Create the mesh and order the nodes
    createNodes(_order);
  }

  // Get the current array of quadrants
  int size;
  TMRQuadrant *array;
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
    TMRQuadrant p;

    if (order == 2){
      // For all searches/comparisons, we use node numbers
      const int use_nodes = 1;

      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      p.level = array[i].level;

      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          p.x = array[i].x + ii*h;
          p.y = array[i].y + jj*h;

          // Get the index for the node
          TMRQuadrant *t = nodes->contains(&p, use_nodes);
          
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
    else if (order == 3){
      // For all searches, use nodes
      const int use_nodes = 1;

      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level - 1);
      p.level = array[i].level+1;

      // Add all of the nodes to the hash
      for ( int jj = 0; jj < 3; jj++ ){
        for ( int ii = 0; ii < 3; ii++ ){
          p.x = array[i].x + ii*h;
          p.y = array[i].y + jj*h;
          
          // Get the index for the node
          TMRQuadrant *t = nodes->contains(&p, use_nodes);
          
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

    // Set the connectivity pointer
    elem_ptr[i+1] = conn_size;
  }
}

/*
  Create the interpolation operator
*/
void TMRQuadtree::createInterpolation( TMRQuadtree *coarse,
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
  TMRQuadrant *fine;
  nodes->getArray(&fine, &fine_size);
  
  // Allocate the arrays
  int conn_len = 0;
  int max_conn_len = order*order*num_nodes;
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
  int max_size = (order*order)*order;
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
      TMRQuadrant *t = coarse->nodes->contains(&fine[i], use_nodes);
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

        if (id == 1 || id == 2){
          // Get the root sibling
          TMRQuadrant n;
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
        else if (id == 3){
          // Get the root sibling
          TMRQuadrant n;
          fine[i].getSibling(0, &n);

          for ( int jj = 0; jj < 2; jj++ ){
            for ( int ii = 0; ii < 2; ii++ ){
              TMRQuadrant p;
              p.x = n.x + 2*h*ii;
              p.y = n.y + 2*h*jj;
              
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
void TMRQuadtree::createRestriction( TMRQuadtree *tree,
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
  TMRQuadrant *coarse;
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
  // find the next-to-adjacent nodes in the quadtree
  for ( int i = 0; i < coarse_size; i++ ){
    if (coarse[i].tag >= 0){
      // Use node-based searching
      const int use_nodes = 1;

      // Keep track of the number of weights used
      int nweights = 0;
      double w = 0.0;

      // Get the fine index number corresponding to the
      // coarse index on this mesh
      TMRQuadrant *fine;
      fine = nodes->contains(&coarse[i], use_nodes);

      // Get the coarse node level
      const int32_t h = 1 << (TMR_MAX_LEVEL - fine->level);
      
      // Scan through the node locations for the coarse mesh
      for ( int jj = 0; jj < 3; jj++ ){
        for ( int ii = 0; ii < 3; ii++ ){
          // Set the location for the adjacent nodes
          TMRQuadrant n;
          n.x = coarse[i].x + h*(ii-1);
          n.y = coarse[i].y + h*(jj-1);
          
          // Get the element level
          TMRQuadrant *t;
          t = nodes->contains(&n, use_nodes);
          if (t){
            double wk = wvals[ii]*wvals[jj];
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
    fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEQUAD\n");

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


