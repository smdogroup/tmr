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
  int h = 1 << (TMR_MAX_LEVEL - refine_level);

  // Set the number of octants created
  int nocts = h*h*h;

  // Create an array of the octants that will be stored
  TMROctant *array = new TMROctant[ nocts ];

  // Create all the octants required
  int index = 0;
  for ( int z = 0; z < h; z++ ){
    for ( int y = 0; y < h; y++ ){
      for ( int x = 0; x < h; x++ ){
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

  // Zero everything else
  num_elements = 0;
  num_nodes = 0;
  num_dependent_nodes = 0;
  nodes = NULL;
}

/*
  Generate a random initial octree that can be used for testing.
*/
TMROctree::TMROctree( int nrand, int min_level, int max_level ){
  // Create an array of the octants that will be stored
  TMROctant *array = new TMROctant[ nrand ];

  // Generate a random number of octants along random directions
  for ( int i = 0; i < nrand; i++ ){
    int level = min_level + (rand() % (max_level - min_level + 1));

    int h = 1 << (TMR_MAX_LEVEL - level);
    int x = h*(rand() % (1 << level));
    int y = h*(rand() % (1 << level));
    int z = h*(rand() % (1 << level));

    array[i].x = x;
    array[i].y = y;
    array[i].z = z;
    array[i].level = level;
  }

  elements = new TMROctantArray(array, nrand);
  elements->sort();
  nodes = NULL;

  num_elements = 0;
  num_nodes = 0;
  num_dependent_nodes = 0;
}

/*
  Create the octree tree from the array
*/
TMROctree::TMROctree( TMROctantArray *_elements ){
  elements = _elements;
  elements->sort();
  nodes = NULL;

  num_elements = 0;
  num_nodes = 0;
  num_dependent_nodes = 0;
}

/*
  Free the octree
*/
TMROctree::~TMROctree(){
  if (elements){ delete elements; }
  if (nodes){ delete nodes; }
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

  The type of balancing - face, edge or corner balanced - is
  determined using the balance_type argument. Face balancing
  is balancing across faces, corner balances across corners of 
  the elements and corner balances across corners.

  The code always balances faces and balances edges if (balance_type &
  2) and balances corners if (balance_type & 2) && (balance_type & 4).

  In other words:
  balance_type = 1 balances faces
  balance_type = 2 or 3 balances edge 
  balance_type = 6 or 7 balances corners
*/
void TMROctree::balance( int balance_type ){
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
  const int hmax = 1 << TMR_MAX_LEVEL;

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
      if (balance_type & 2){
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
      }

      // If we're balancing across edges and 
      if ((balance_type & 2) && (balance_type & 4)){
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
      if (balance_type & 2){
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
      }

      // Add the octants across corners
      if ((balance_type & 2) && (balance_type & 4)){
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
      hash->addOctant(&q);
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
  elements->sort();

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
  Create the mesh using the internal octree data
*/
void TMROctree::createNodes( int order, 
                             int **_dep_ptr, 
                             int **_dep_conn,
                             double **_dep_weights ){
  // Create the octant hash that we'll add the element nodes
  TMROctantHash *hash = new TMROctantHash();

  // Get the current array of octants
  int size;
  TMROctant *array;
  elements->getArray(&array, &size);

  // Loop over all of the current octants and add the nodes
  for ( int i = 0; i < size; i++ ){
    // The octant that we'll have to add
    TMROctant p;

    if (order == 2){
      // Add all of the nodes from the adjacent elements
      const int h = 1 << (TMR_MAX_LEVEL - array[i].level);
      p.level = array[i].level;
      p.tag = 1;

      // Add all of the nodes to the hash
      for ( int kk = 0; kk < 2; kk++ ){
        for ( int jj = 0; jj < 2; jj++ ){
          for ( int ii = 0; ii < 2; ii++ ){
            p.x = array[i].x + ii*h;
            p.y = array[i].y + jj*h;
            p.z = array[i].z + kk*h;

            // Add the octant
            hash->addOctant(&p);
          }
        }
      }
    }
    else if (order == 3){
      // Add all of the nodes from the adjacent elements
      const int h = 1 << (TMR_MAX_LEVEL - array[i].level - 1);
      p.level = array[i].level+1;
      p.tag = 1;

      // Add all of the nodes to the hash
      for ( int kk = 0; kk < 3; kk++ ){
        for ( int jj = 0; jj < 3; jj++ ){
          for ( int ii = 0; ii < 3; ii++ ){
            p.x = array[i].x + ii*h;
            p.y = array[i].y + jj*h;
            p.z = array[i].z + kk*h;

            // Add the octant
            hash->addOctant(&p);
          }
        }
      }
    }
  }

  // Cover the hash table to a list and uniquely sort it
  nodes = hash->toArray();
  nodes->sort();

  // Free the hash
  delete hash;

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
      const int h = 1 << (TMR_MAX_LEVEL - array[i].level);      
      
      // Get the size of the next finest level
      const int hc = 1 << (TMR_MAX_LEVEL - array[i].level - 1);
      
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
	for ( int ii = 0; ii < 1; ii++ ){
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
  
  // Allocate the arrays for the dependent nodes
  int *dep_ptr = new int[ num_dependent_nodes+1 ];
  int *dep_conn = new int[ order*nedges + order*order*nfaces ];
  double *dep_weights = new double[ order*nedges + order*order*nfaces ];

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

  // Set the pointers to the return data
  *_dep_ptr = dep_ptr;
  *_dep_conn = dep_conn;
  *_dep_weights = dep_weights;

}

/*
  Create the mesh connectivity
*/
void TMROctree::createMesh( int order,
                            int *_num_nodes,
                            int *_num_dep_nodes,
                            int *_num_elements,
                            int **_elem_ptr, 
                            int **_elem_conn,
			    int **_dep_ptr,
			    int **_dep_conn, 
			    double **_dep_weights ){
  // Scan through and create the connectivity for all the
  // elements in the Morton order
  if (order < 2){ order = 2; }
  if (order > 3){ order = 3; }

  // Create the mesh and order the nodes
  createNodes(order, _dep_ptr, 
	      _dep_conn, _dep_weights);

  // Get the current array of octants
  int size;
  TMROctant *array;
  elements->getArray(&array, &size);

  // We already know how large the connectivity list
  // will be and the size of the number of dependent nodes
  int *elem_ptr = new int[ num_elements+1 ];
  int *elem_conn = new int[ order*order*order*num_elements ];
  
  // Set the connectivity pointer
  int conn_size = 0;
  elem_ptr[0] = 0;
  for ( int i = 0; i < size; i++ ){
    TMROctant p;

    if (order == 2){
      // For all searches/comparisons, we use node numbers
      const int use_nodes = 1;

      // Add all of the nodes from the adjacent elements
      const int h = 1 << (TMR_MAX_LEVEL - array[i].level);
      p.level = array[i].level;

      for ( int kk = 0; kk < 2; kk++ ){
        for ( int jj = 0; jj < 2; jj++ ){
          for ( int ii = 0; ii < 2; ii++ ){
            p.x = array[i].x + ii*h;
            p.y = array[i].y + jj*h;
            p.z = array[i].z + kk*h;

            // Get the index for the node
            TMROctant *t = nodes->contains(&p, use_nodes);
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
      const int h = 1 << (TMR_MAX_LEVEL - array[i].level - 1);
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
            elem_conn[conn_size] = t->tag;
            conn_size++;
          }
        }
      }
    }

    // Set the connectivity pointer
    elem_ptr[i+1] = conn_size;
  }

  // Set the nodes, dependent nodes, elements
  *_num_nodes = num_nodes;
  *_num_dep_nodes = num_dependent_nodes;
  *_num_elements = num_elements;

  // Set the connectivity arrays
  *_elem_ptr = elem_ptr;
  *_elem_conn = elem_conn;
}

/*
  Print out the octree to a file for visualization
*/
void TMROctree::printTree( const char * filename ){
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
      int h = 1 << (TMR_MAX_LEVEL - array[i].level);
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


