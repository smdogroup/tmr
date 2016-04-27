#include "TMRForrest.h"

const int face_to_edge_nodes[][2] = {{0, 2},
                                     {1, 3},
                                     {0, 1},
                                     {2, 3}};

/*
  Create the TMRQuadForrest object
*/
TMRQuadForrest::TMRQuadForrest( MPI_Comm _comm ){
  // Set the MPI communicator
  comm = _comm;

  // Zero out the nodes/edges/faces and all data
  num_nodes = 0;
  num_edges = 0;
  num_faces = 0;

  face_conn = NULL;
  face_edge_conn = NULL;
  node_face_ptr = NULL;
  node_face_conn = NULL;
  edge_face_ptr = NULL;
  edge_face_conn = NULL;
  quadtrees = NULL;
}

/*
  Free the data allocated by the TMRQuadForrest object
*/
TMRQuadForrest::~TMRQuadForrest(){
  if (face_conn){ delete [] face_conn; }
  if (face_edge_conn){ delete [] face_edge_conn; }
  if (node_face_ptr){ delete [] node_face_ptr; }
  if (node_face_conn){ delete [] node_face_conn; }
  if (edge_face_ptr){ delete [] edge_face_ptr; }
  if (edge_face_conn){ delete [] edge_face_conn; }
  if (quadtrees){ 
    for ( int i = 0; i < num_faces; i++ ){
      if (quadtrees[i]){ delete quadtrees[i]; }
    }
    delete [] quadtrees;
  }
}

/*
  Set the connectivity of the faces

  This code sets the face connectivity and generates the following
  additional data that are required:

  1. Face to node connectivity (input)
  2. Node to face connectivity (required for corner balancing)
  3. Unique edge ordering
  4. Face to edge and edge to face connectivity

  This information is required for creating quadtree forrests on the
  unstructured mesh.
*/
void TMRQuadForrest::setConnectivity( int _num_nodes,
                                      const int *_face_conn,
                                      int _num_faces ){
  // Free any data if it has already been allocated. 
  // This will erase everything internally.
  if (face_conn){ delete [] face_conn; }
  if (face_edge_conn){ delete [] face_edge_conn; }
  if (node_face_ptr){ delete [] node_face_ptr; }
  if (node_face_conn){ delete [] node_face_conn; }
  if (edge_face_ptr){ delete [] edge_face_ptr; }
  if (edge_face_conn){ delete [] edge_face_conn; }
  if (quadtrees){ 
    for ( int i = 0; i < num_faces; i++ ){
      if (quadtrees[i]){ delete quadtrees[i]; }
    }
    delete [] quadtrees;
  }

  // Copy over the data locally
  num_nodes = _num_nodes;
  num_faces = _num_faces;
  num_edges = 0;

  // Copy over the face connectivity
  face_conn = new int[ 4*num_faces ];
  memcpy(face_conn, _face_conn, 4*num_faces*sizeof(int));

  // Create the data structure for the node to face connectivity
  node_face_ptr = new int[ num_nodes+1 ];
  memset(node_face_ptr, 0, (num_nodes+1)*sizeof(int));

  // Count the number of times each node is referred to
  for ( int i = 0; i < 4*num_faces; i++ ){
    node_face_ptr[face_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_nodes+1; i++ ){
    node_face_ptr[i] += node_face_ptr[i-1];
  }

  // Allocate the full node to face pointer array
  node_face_conn = new int[ node_face_ptr[num_nodes] ];
  for ( int i = 0; i < num_faces; i++ ){
    for ( int j = 4*i; j < 4*(i+1); j++ ){
      int node = face_conn[j];
      node_face_conn[node_face_ptr[node]] = i;
      node_face_ptr[node]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_nodes; i >= 1; i-- ){
    node_face_ptr[i] = node_face_ptr[i-1];
  }
  node_face_ptr[0] = 0;

  // Now establish a unique ordering of the edges along each face
  face_edge_conn = new int[ 4*num_faces ];
  for ( int i = 0; i < 4*num_faces; i++ ){
    face_edge_conn[i] = -1;
  }

  int edge = 0;
  for ( int i = 0; i < num_faces; i++ ){
    // Loop over each edge in the mesh
    for ( int j = 0; j < 4; j++ ){
      if (face_edge_conn[4*i + j] < 0){
        int n1 = face_conn[4*i + face_to_edge_nodes[j][0]];
        int n2 = face_conn[4*i + face_to_edge_nodes[j][1]];

        // Keep track of the number of edges found
        const int max_nedges = 20;
        int edge_index[max_nedges];
        int nedges = 1;
        edge_index[0] = 4*i + j;

        // Set the edge number - if any is found
        int edge_num = -1;

        // Scan through the faces that share the same
        // node and check if any of the edges are also
        // shared 
        for ( int ip = node_face_ptr[n1];
              ip < node_face_ptr[n1+1]; ip++ ){
          int ii = node_face_conn[ip];
          
          // Loop over each edge in the new face
          for ( int jj = 0; jj < 4; jj++ ){
            int nn1 = face_conn[4*ii + face_to_edge_nodes[jj][0]];
            int nn2 = face_conn[4*ii + face_to_edge_nodes[jj][1]];

            // Check if the face matches
            if ((n1 == nn1 && n2 == nn2) ||
                (n1 == nn2 && n2 == nn1)){
              if (face_edge_conn[4*ii + jj] >= 0){
                // If this edge has been ordered, copy over
                // the edge number
                edge_num = face_edge_conn[4*ii + jj];
              }
              else if (nedges < max_nedges){
                // This edge has not yet been ordered, add it
                // to the unordered list if there is still room
                // if not, we will detect and order it during
                // a future iteration
                edge_index[nedges] = 4*ii + jj;
                nedges++;
              }
            }
          }
        }

        // If this edge does not have a edge number, assign
        // a new one to the list
        edge_num = edge;
        edge++;

        // Set the edge numbers for all the edges that we found
        for ( int ii = 0; ii < nedges; ii++ ){
          face_edge_conn[edge_index[ii]] = edge_num;
        }
      }
    }
  }

  // Set the total number of edges
  num_edges = edge;
 
  // Create the data structure for the edge to face connectivity
  edge_face_ptr = new int[ num_edges+1 ];
  memset(edge_face_ptr, 0, (num_edges+1)*sizeof(int));

  // Count the number of times each edge is referred to
  for ( int i = 0; i < 4*num_faces; i++ ){
    edge_face_ptr[face_edge_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_edges+1; i++ ){
    edge_face_ptr[i] += edge_face_ptr[i-1];
  }

  // Allocate the full node to face pointer array
  edge_face_conn = new int[ edge_face_ptr[num_edges] ];
  for ( int i = 0; i < num_faces; i++ ){
    for ( int j = 4*i; j < 4*(i+1); j++ ){
      int e = face_edge_conn[j];
      edge_face_conn[edge_face_ptr[e]] = i;
      edge_face_ptr[e]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_nodes; i >= 1; i-- ){
    edge_face_ptr[i] = edge_face_ptr[i-1];
  }
  edge_face_ptr[0] = 0;
}

/*
  Create a forrest with the specified refinement level
*/
void TMRQuadForrest::createTrees( int refine_level ){
  // Create the quadtrees
  quadtrees = new TMRQuadtree*[ num_faces ];  
  for ( int i = 0; i < num_faces; i++ ){
    // quadtrees[i] = new TMRQuadtree(refine_level);
    quadtrees[i] = new TMRQuadtree(100, 0, 8);
  }
}

/*
  Given the quadrant 
*/
void TMRQuadForrest::addEdgeNeighbors( int tree,
                                       int edge, 
                                       TMRQuadrant p,
                                       TMRQuadrantHash **hash,
                                       TMRQuadrantQueue **queue ){
  // First determine the global edge number 
  int edge_num = face_edge_conn[4*tree + edge];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Compute the x/y coordinate along the edge
  int32_t ucoord = 0;
  if (edge == 0 || edge == 1){
    ucoord = p.y; 
  }
  else {
    ucoord = p.x;
  }

  // Retrieve the first and second edge numbers
  int n1 = face_conn[4*tree + face_to_edge_nodes[edge][0]];
  int n2 = face_conn[4*tree + face_to_edge_nodes[edge][1]];

  // Now, cycle through all the adjacent trees
  for ( int ip = edge_face_ptr[edge_num];
        ip < edge_face_ptr[edge_num+1]; ip++ ){

    // Get the trees that are adjacent across this edge
    int adjacent = edge_face_conn[ip];
    if (adjacent != tree){
      for ( int j = 0; j < 4; j++ ){
        int nn1 = face_conn[4*adjacent + face_to_edge_nodes[j][0]];
        int nn2 = face_conn[4*adjacent + face_to_edge_nodes[j][1]];

        // Add the quadrant to the list
        if (n1 == nn1 && n2 == nn2){
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          if (j == 0){
            neighbor.x = 0;
            neighbor.y = ucoord;
          }
          else if (j == 1){
            neighbor.x = hmax - h;
            neighbor.y = ucoord;
          }
          else if (j == 2){
            neighbor.x = ucoord;
            neighbor.y = 0; 
          }
          else {
            neighbor.x = ucoord;
            neighbor.y = hmax - h;
          }
          
          // Add the quadrant to the list
          if (hash[adjacent]->addQuadrant(&neighbor)){
            queue[adjacent]->push(&neighbor);
          }
        }
        else if (n1 == nn2 && n2 == nn1){
          // The edges have the opposite orientation
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          if (j == 0){
            neighbor.x = 0;
            neighbor.y = hmax - h - ucoord;
          }
          else if (j == 1){
            neighbor.x = hmax - h;
            neighbor.y = hmax - h - ucoord;
          }
          else if (j == 2){
            neighbor.x = hmax - h - ucoord;
            neighbor.y = 0; 
          }
          else {
            neighbor.x = hmax - h - ucoord;
            neighbor.y = hmax - h;
          }

          // Add the quadrant to the list
          if (hash[adjacent]->addQuadrant(&neighbor)){
            queue[adjacent]->push(&neighbor);
          }
        }
      }
    }
  }
}

/*
  Add the corner neighbors
*/
void TMRQuadForrest::addCornerNeighbors( int tree,
                                         int corner, 
                                         TMRQuadrant p,
                                         TMRQuadrantHash **hash,
                                         TMRQuadrantQueue **queue ){
  // First determine the global edge number 
  int node_num = face_conn[4*tree + corner];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent trees
  for ( int ip = node_face_ptr[node_num];
        ip < node_face_ptr[node_num+1]; ip++ ){
      
    // Get the trees that are adjacent across this edge
    int adjacent = node_face_conn[ip];
    if (adjacent != tree){
      for ( int j = 0; j < 4; j++ ){
        if (face_conn[4*adjacent+j] == node_num){
          // Add the corner node to the list
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          if (j == 0){
            neighbor.x = 0;
            neighbor.y = 0;
          }
          else if (j == 1){
            neighbor.x = hmax - h;
            neighbor.y = 0;
          }
          else if (j == 2){
            neighbor.x = 0;
            neighbor.y = hmax - h; 
          }
          else {
            neighbor.x = hmax - h;
            neighbor.y = hmax - h;
          }

          // Add the quadrant to the list
          if (hash[adjacent]->addQuadrant(&neighbor)){
            queue[adjacent]->push(&neighbor);
          }
        }
      }
    }
  }
}

/*
  Balance the forrest
*/
void TMRQuadForrest::balance( int balance_corner ){
  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Create a hash table for the balanced tree
  TMRQuadrantHash **hash = new TMRQuadrantHash*[ num_faces ];
  TMRQuadrantQueue **queue = new TMRQuadrantQueue*[ num_faces ];
  for ( int face = 0; face < num_faces; face++ ){
    hash[face] = new TMRQuadrantHash();
    queue[face] = new TMRQuadrantQueue();
  }

  for ( int face = 0; face < num_faces; face++ ){
    // Go through the existing list of quadrants and add up everything
    // for balancing
    TMRQuadrant p, q, neighbor;

    // Get the current array of quadrants
    TMRQuadrantArray *elements;
    quadtrees[face]->getElements(&elements);

    // Get the array of elements
    int size;
    TMRQuadrant *array;
    elements->getArray(&array, &size);
    
    // Add the element
    for ( int i = 0; i < size; i++ ){
      // Try adding all of the children
      array[i].getSibling(0, &q);
      hash[face]->addQuadrant(&q);
      
      // Get the parent of the quadrant, and add the their
      // face-matched quadrants from each face, as long 
      // as they fall within the bounds
      if (q.level > 0){
        q.parent(&p);
        
        // If we are balancing across edges
        // add the edge-adjacent elements
        for ( int edge = 0; edge < 4; edge++ ){
          p.edgeNeighbor(edge, &neighbor);
          neighbor.getSibling(0, &q);
	  
          // If we're in bounds, add the neighbor
          if (q.x >= 0 && q.y >= 0 &&
              q.x < hmax && q.y < hmax){
            if (hash[face]->addQuadrant(&q)){
              queue[face]->push(&q);
            }
          }
          else {
            // Add the quadrant to the other trees
            addEdgeNeighbors(face, edge, q, hash, queue);
          }
        }
        
        // If we're balancing across edges and 
        if (balance_corner){
          for ( int corner = 0; corner < 4; corner++ ){
            p.cornerNeighbor(corner, &neighbor);
            neighbor.getSibling(0, &q);
            
            // If we're in bounds, add the neighbor
            if (q.x >= 0 && q.y >= 0 &&
                q.x < hmax && q.y < hmax){
              if (hash[face]->addQuadrant(&q)){
                queue[face]->push(&q);
              }
            }
            else if ((q.x >= 0 && q.x < hmax) ||
                     (q.y >= 0 && q.y < hmax)){
              int edge = 0;
              if (q.x >= hmax){ edge = 1; }
              else if (q.y < 0){ edge = 2; }
              else if (q.y >= hmax){ edge = 3; }
              addEdgeNeighbors(face, edge, q, hash, queue);
            }
            else {
              // Add the quadrant to the other trees
              addCornerNeighbors(face, corner, neighbor, hash, queue);
            }
          }
        }
      }
    }
  }

  int queue_length_flag = 1;
  while (queue_length_flag){
    // Now continue until the queue of added quadrants is
    // empty. At each iteration, pop an quadrant and add 
    // its neighbours until nothing new is added. This code
    // handles the propagation of quadrants to adjacent quadrants.
    for ( int face = 0; face < num_faces; face++ ){
      while (queue[face]->length() > 0){
        TMRQuadrant p, q, neighbor;    
        q = queue[face]->pop();
      
        if (q.level > 1){
          q.parent(&p);

          // Add the quadrants across adjacent edges
          for ( int edge = 0; edge < 4; edge++ ){
            p.edgeNeighbor(edge, &neighbor);
            neighbor.getSibling(0, &q);
            
            // If we're in bounds, add the neighbor
            if (q.x >= 0 && q.y >= 0 &&
                q.x < hmax && q.y < hmax){
              if (hash[face]->addQuadrant(&q)){
                queue[face]->push(&q);
              }
            }
            else {
              // Add the quadrant to the other trees
              addEdgeNeighbors(face, edge, q, hash, queue);
            }
          }
          
          // Add the quadrants across corners
          if (balance_corner){
            for ( int corner = 0; corner < 4; corner++ ){
              p.cornerNeighbor(corner, &neighbor);
              neighbor.getSibling(0, &q);
              
              // If we're in bounds, add the neighbor
              if (q.x >= 0 && q.y >= 0 &&
                  q.x < hmax && q.y < hmax){
                if (hash[face]->addQuadrant(&q)){
                  queue[face]->push(&q);
                }
              }
              else if ((q.x >= 0 && q.x < hmax) ||
                       (q.y >= 0 && q.y < hmax)){
                int edge = 0;
                if (q.x >= hmax){ edge = 1; }
                else if (q.y < 0){ edge = 2; }
                else if (q.y >= hmax){ edge = 3; }
                addEdgeNeighbors(face, edge, q, hash, queue);
              }
              else {
                // Add the quadrant to the other trees
                addCornerNeighbors(face, corner, neighbor, hash, queue);
              }
            }
          }
        }
      }
    }

    queue_length_flag = 0;
    for ( int face = 0; face < num_faces; face++ ){
      if (queue[face]->length() > 0){
        queue_length_flag = 1;
        break;
      }
    }
  }

  // Free the queues - they are no longer required
  for ( int face = 0; face < num_faces; face++ ){
    delete queue[face];
  }
  delete [] queue;

  // Convert the hash tables back to the elements
  for ( int face = 0; face < num_faces; face++ ){
    // Now convert the elements from child-0 elements to
    // elements which cover the full mesh
    TMRQuadrantArray *child0_elems = hash[face]->toArray();
    int size;
    TMRQuadrant *array;
    child0_elems->getArray(&array, &size);

    // Loop over all elements and add their siblings
    for ( int i = 0; i < size; i++ ){
      if (array[i].level > 0){
        for ( int j = 0; j < 4; j++ ){
          TMRQuadrant q;
          array[i].getSibling(j, &q);
          hash[face]->addQuadrant(&q);
        }
      }
    }

    // Free the temporary elements
    delete child0_elems;

    // Set the elements into the quadtree
    TMRQuadrantArray *elements = hash[face]->toArray();
    elements->sort();
    quadtrees[face]->setElements(elements);

    // Free the corresponding hash
    delete hash[face];
  }
  delete hash;
}

/*
  Coarsen the entire forrest
*/
TMRQuadForrest* TMRQuadForrest::coarsen(){
  TMRQuadForrest *coarse = NULL;
  if (quadtrees){
    coarse = new TMRQuadForrest(comm);

    // Copy over the connectivity data 
    coarse->num_nodes = num_nodes;
    coarse->num_edges = num_edges;
    coarse->num_faces = num_faces;
    
    // Allocate/copy the face connectivity
    coarse->face_conn = new int[ 4*num_faces ];
    memcpy(coarse->face_conn, face_conn, 4*num_faces*sizeof(int));

    coarse->face_edge_conn = new int[ 4*num_faces ];
    memcpy(coarse->face_edge_conn, face_edge_conn, 4*num_faces*sizeof(int));
    
    // Allocate the remaining data
    coarse->node_face_ptr = new int[ num_nodes+1 ];
    memcpy(coarse->node_face_ptr, node_face_ptr, (num_nodes+1)*sizeof(int));

    coarse->node_face_conn = new int[ node_face_ptr[num_nodes] ];
    memcpy(coarse->node_face_conn, node_face_conn, 
           node_face_ptr[num_nodes]*sizeof(int));
    
    coarse->edge_face_ptr = new int[ num_edges+1 ];
    memcpy(coarse->edge_face_ptr, edge_face_ptr, (num_edges+1)*sizeof(int));
    
    coarse->edge_face_conn = new int[ edge_face_ptr[num_edges] ];
    memcpy(coarse->edge_face_conn, edge_face_conn,
           edge_face_ptr[num_edges]*sizeof(int));

    // Coarsen all the quadtrees
    coarse->quadtrees = new TMRQuadtree*[ num_faces ];
    for ( int i = 0; i < num_faces; i++ ){
      coarse->quadtrees[i] = quadtrees[i]->coarsen();
    }
  }

  return coarse;
}

/*
  Create the nodes within the mesh
*/
void TMRQuadForrest::createNodes( int _order ){
  // Allocate all possible nodes on all of the trees
  for ( int face = 0; face < num_faces; face++ ){
    quadtrees[face]->createNodes(_order);
  }
  /*
  // Order the nodes
  if (order == 2){
    const int use_nodes = 1;

    // Set the dependent nodes for each face
    for ( int face = 0; face < num_faces; face++ ){
      TMRQuadrantArray *elements;
      quadtrees[face]->getElements(&elements);

      int size;
      TMRQuadrant *array;
      elements->getArray(&array, &size);
    
      // Now check for dependent nodes on each edge and face
      for ( int i = 0; i < size; i++ ){
        // Get the side length of the element
        const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);      
      
        













        // Get the size of the next finest level
        const int32_t hc = 1 << (TMR_MAX_LEVEL - array[i].level - 1);

        // Check if the element across an edge exists
        
        for ( int k = 0; k < 2; k++ ){ // Each coordinate direction
          for ( int ii = 0; ii < 2; ii++ ){
            TMRQuadrant c;
            if (k == 0){
              c.x = array[i].x + ii*h;
              c.y = array[i].y + hc;
            }
            else {
              c.x = array[i].x + hc;
              c.y = array[i].y + ii*h;
            }

            const int use_nodes = 1;
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
             
                // Set the other independent node
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
  }
  else if (order == 3){

  }

  int ndep_nodes = 0;
  int nnodes = 0;

  // Now, order the nodes and dependent nodes within
  // the entire quadtree forrest
  for ( int face = 0; face < num_faces; face++ ){
    // First find the face owners for the edges on
    // this block
    int face_owners[4];
    edge_owners[0] = edge_owners[1] = 
      edge_owners[2] = edge_owners[3] = face;
    
    for ( int j = 0; j < 4; j++ ){
      if (face_edge_conn[4*face+j] < edge_owners[j]){
        edge_owners[j] = face_edge_conn[4*face+j];
      }
    }

    // For each block, cycle through the nodes and
    // assign them, if they are located on another tree,
    // then transform
    for ( int i = 0; i < size; i++ ){
      if (node_array[i].x == 0){
        
      }
      else if (node_array[i].y == 0){
        
      }
      else if (node_array[i].x == hmax){
        
      }
      else if (node_array[i].y == hmax){

      }
      else {
        if (node_array[i].tag >= 0){
          node_array[i].tag = nnodes;
          nnodes++;
        }
        else {
          node_array[i].tag = -(ndep_nodes+1);
          ndep_nodes++;
        }
      }
    }
  }
  */
}
