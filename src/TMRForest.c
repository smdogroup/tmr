#include "TMRForest.h"

const int face_to_edge_nodes[][2] = {{0, 2},
                                     {1, 3},
                                     {0, 1},
                                     {2, 3}};

/*
  Create the TMRQuadForest object
*/
TMRQuadForest::TMRQuadForest( MPI_Comm _comm ){
  // Set the MPI communicator
  comm = _comm;

  // Zero out the nodes/edges/faces and all data
  num_nodes = 0;
  num_edges = 0;
  num_faces = 0;

  // NULL out everything
  face_conn = NULL;
  face_edge_conn = NULL;
  node_face_ptr = NULL;
  node_face_conn = NULL;
  edge_face_ptr = NULL;
  edge_face_conn = NULL;
  quadtrees = NULL;

  // Set the size of the mesh
  mesh_order = 0;
  num_mesh_nodes = 0;
  num_mesh_elements = 0;

  // Set the dependent weight information
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;
}

/*
  Free the data allocated by the TMRQuadForest object
*/
TMRQuadForest::~TMRQuadForest(){
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
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
}

/*
  Set the connectivity of the faces

  This code sets the face connectivity and generates the following
  additional data that are required:

  1. Face to node connectivity (input)
  2. Node to face connectivity (required for corner balancing)
  3. Unique edge ordering
  4. Face to edge and edge to face connectivity

  This information is required for creating quadtree forests on the
  unstructured mesh.
*/
void TMRQuadForest::setConnectivity( int _num_nodes,
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

  // Keep track of the edge numbers
  int edge = 0;
  for ( int i = 0; i < num_faces; i++ ){
    // Loop over each edge on this face
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

        // If this edge does not have an edge number, assign
        // a new one to the list
        if (edge_num < 0){
          edge_num = edge;
          edge++;
        }

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
  for ( int i = num_edges; i >= 1; i-- ){
    edge_face_ptr[i] = edge_face_ptr[i-1];
  }
  edge_face_ptr[0] = 0;
}

/*
  Create a forest with the specified refinement level
*/
void TMRQuadForest::createTrees( int refine_level ){
  if (quadtrees){ 
    for ( int i = 0; i < num_faces; i++ ){
      delete quadtrees[i];
    }
    delete [] quadtrees;
  }

  // Create the quadtrees
  quadtrees = new TMRQuadtree*[ num_faces ];  
  for ( int i = 0; i < num_faces; i++ ){
    quadtrees[i] = new TMRQuadtree(refine_level);
  }
}

/*
  Create a random set of trees

  This function is usually used for testing purposes.
*/
void TMRQuadForest::createRandomTrees( int nrand, 
                                       int min_level, int max_level ){
  if (quadtrees){ 
    for ( int i = 0; i < num_faces; i++ ){
      delete quadtrees[i];
    }
    delete [] quadtrees;
  }

  // Create a random set of quadtrees
  quadtrees = new TMRQuadtree*[ num_faces ];  
  for ( int i = 0; i < num_faces; i++ ){
    quadtrees[i] = new TMRQuadtree(nrand, min_level, max_level);
  }
}

/*
  Add the edge neighbors for a adjacent trees

  This function is called to balance the forest across tree edges.
  Given an quadrant p on the specified corner index, this code ensures
  a corner balanced tree, by adding the corresponding corner qudrant
  to all node-adjacent faces. If the quadrant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  face:    the local face index (p is defined on this face)
  corner:  the tree corner index (p must lie on this corner)
  p:       the local quadrant
  hash:    the array of hash objects for each face
  queue:   the array of newly added qudrants for each face
*/
void TMRQuadForest::addEdgeNeighbors( int face,
                                      int edge, 
                                      TMRQuadrant p,
                                      TMRQuadrantHash **hash,
                                      TMRQuadrantQueue **queue ){
  // First determine the global edge number 
  int edge_num = face_edge_conn[4*face + edge];
  
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
  int n1 = face_conn[4*face + face_to_edge_nodes[edge][0]];
  int n2 = face_conn[4*face + face_to_edge_nodes[edge][1]];

  // Now, cycle through all the adjacent faces
  for ( int ip = edge_face_ptr[edge_num];
        ip < edge_face_ptr[edge_num+1]; ip++ ){

    // Get the faces that are adjacent across this edge
    int adjacent = edge_face_conn[ip];
    if (adjacent != face){
      for ( int j = 0; j < 4; j++ ){
        int nn1 = face_conn[4*adjacent + face_to_edge_nodes[j][0]];
        int nn2 = face_conn[4*adjacent + face_to_edge_nodes[j][1]];

        // Add the quadrant to the list
        if (n1 == nn1 && n2 == nn2){
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          if (j == 0 || j == 1){
            neighbor.x = (hmax - h)*j;
            neighbor.y = ucoord;
          }
          else if (j == 2 || j == 3){
            neighbor.x = ucoord;
            neighbor.y = (hmax - h)*(j % 2); 
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
          if (j == 0 || j == 1){
            neighbor.x = (hmax - h)*j;
            neighbor.y = hmax - h - ucoord;
          }
          else if (j == 2 || j == 3){
            neighbor.x = hmax - h - ucoord;
            neighbor.y = (hmax - h)*(j % 2); 
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
  Add the corner neighbors for a given tree

  This function is called to balance the forest across tree corners.
  Given an quadrant p on the specified corner index, this code ensures
  a corner balanced tree, by adding the corresponding corner qudrant
  to all node-adjacent faces. If the quadrant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  face:    the local face index (p is defined on this face)
  corner:  the tree corner index (p must lie on this corner)
  p:       the local quadrant
  hash:    the array of hash objects for each face
  queue:   the array of newly added qudrants for each face
*/
void TMRQuadForest::addCornerNeighbors( int face,
                                        int corner, 
                                        TMRQuadrant p,
                                        TMRQuadrantHash **hash,
                                        TMRQuadrantQueue **queue ){
  // First determine the global edge number 
  int node_num = face_conn[4*face + corner];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = node_face_ptr[node_num];
        ip < node_face_ptr[node_num+1]; ip++ ){
      
    // Get the faces that are adjacent across this edge
    int adjacent = node_face_conn[ip];
    if (adjacent != face){
      for ( int j = 0; j < 4; j++ ){
        if (face_conn[4*adjacent+j] == node_num){
          // Add the corner node to the list
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          neighbor.x = (hmax - h)*(j%2);
          neighbor.y = (hmax - h)*(j/2);

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
  Add the dependent nodes on all adjacent faces that exist due to
  input quadrant on the local face.

  To ensure inter-face compatibility, nodes that are hanging between
  adjacent faces must be labeled as dependent nodes (with a negative
  tag value.) This function labels all such nodes that are adjacent to
  the quadrant p on the local face across the given edge. These
  dependent nodes are appended to the given array of queues that store
  the dependent node indices.
  
  Note that the dependent nodes are duplicated and are always ordered
  locally on each face independently. This simplifies the
  implementation considerably.

  input:
  face:    the local face index (p is defined on this face)
  edge:    the edge number to search (p must be touching this edge)
  p:       the refined quadrant on the local face
  source:  the source quadrant on the local face
  queue:   the array of queues containing the local dependent nodes
  indep:   the independent node queue - nodes for the dependent constraints
*/
void TMRQuadForest::addEdgeDependentNodes( int face,
                                           int edge,
                                           TMRQuadrant p,
                                           TMRQuadrant source,
                                           TMRQuadrantQueue **queue,
                                           TMRQuadrantQueue **indep ){
  // First determine the global edge number
  int edge_num = face_edge_conn[4*face + edge];
  
  // Compute the edge lengths. Note that the p edge length is
  // guaranteed to be precisely one level below the source edge
  // length.  
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t hs = 1 << (TMR_MAX_LEVEL - source.level);
  const int32_t hp = 1 << (TMR_MAX_LEVEL - p.level);

  // Compute the x/y coordinate along the edge
  int32_t up = p.x;
  int32_t us = source.x;
  if (edge == 0 || edge == 1){
    up = p.y;
    us = source.y;
  }

  // Retrieve the first and second edge numbers
  int n1 = face_conn[4*face + face_to_edge_nodes[edge][0]];
  int n2 = face_conn[4*face + face_to_edge_nodes[edge][1]];

  // Now, cycle through all the adjacent faces
  for ( int ip = edge_face_ptr[edge_num];
        ip < edge_face_ptr[edge_num+1]; ip++ ){

    // Get the faces that are adjacent across this edge
    int adjacent = edge_face_conn[ip];
    if (adjacent != face){
      // Get the quadrant elements
      TMRQuadrantArray *elements;
      quadtrees[adjacent]->getElements(&elements);

      for ( int j = 0; j < 4; j++ ){
        int nn1 = face_conn[4*adjacent + face_to_edge_nodes[j][0]];
        int nn2 = face_conn[4*adjacent + face_to_edge_nodes[j][1]];

        // Check for the orientation of the edges
        int dir = (n1 == nn2 && n2 == nn1);
        if ((n1 == nn1 && n2 == nn2) || dir){
          // First, transform the test p element to the adjacent
          // face coordinate system
          int32_t utp = (dir ? hmax - up - hp : up);
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          if (j == 0 || j == 1){
            neighbor.x = (hmax - hp)*j;
            neighbor.y = utp;
          }
          else if (j == 2 || j == 3){
            neighbor.x = utp;
            neighbor.y = (hmax - hp)*(j % 2);
          }
          
          // If the element exists, we have to add the dependent node
          // corresponding to the source element
          const int use_nodes = 0;
          if (elements->contains(&neighbor, use_nodes)){
            int32_t uts = (dir ? hmax - us - hs : us);

            if (mesh_order == 2){
              // Add the single node for the
              TMRQuadrant node;
              node.tag = -1;
              node.level = p.level;
              if (j == 0 || j == 1){
                // Add the dependent node
                node.x = hmax*j;
                node.y = uts + hp;
                queue[adjacent]->push(&node);

                // Push the independent nodes onto the queue
                node.y = uts;
                indep[adjacent]->push(&node);
                node.y = uts + hs;
                indep[adjacent]->push(&node);
              }
              else if (j == 2 || j == 3){
                // Add the dependent node
                node.x = uts + hp;
                node.y = hmax*(j % 2);
                queue[adjacent]->push(&node);

                // Push the independent nodes onto the queue
                node.x = uts;
                indep[adjacent]->push(&node);
                node.x = uts + hs;
                indep[adjacent]->push(&node);
              }
            }
            else if (mesh_order == 3){              
              const int32_t hr = 1 << (TMR_MAX_LEVEL - p.level - 1);

              // Set the dependent nodes along the edge
              TMRQuadrant node;
              node.tag = -1;
              node.level = p.level;
              if (j == 0 || j == 1){
                // Push the two dependent nodes onto the queue
                node.x = hmax*j;
                node.y = uts + hr;
                queue[adjacent]->push(&node);
                node.y = uts + 3*hr;
                queue[adjacent]->push(&node);

                // Push the nodes for the first constraint
                for ( int k = 0; k < 3; k++ ){
                  node.y = uts + k*hp;
                  indep[adjacent]->push(&node);
                }
                
                // Push the nodes for the second constraint
                for ( int k = 2; k >= 0; k-- ){
                  node.y = uts + k*hp;
                  indep[adjacent]->push(&node);
                }
              }
              else if (j == 2 || j == 3){
                // Push the two dependent nodes onto the queue
                node.x = uts + hr;
                node.y = hmax*(j % 2);
                queue[adjacent]->push(&node);
                node.x = uts + 3*hr;
                queue[adjacent]->push(&node);

                // Push the nodes for the first constraint
                for ( int k = 0; k < 3; k++ ){
                  node.x = uts + k*hp;
                  indep[adjacent]->push(&node);
                }
                
                // Push the nodes for the second constraint
                for ( int k = 2; k >= 0; k-- ){
                  node.x = uts + k*hp;
                  indep[adjacent]->push(&node);
                }
              }
            }
          }
        }
      }
    }
  }
}

/*
  Find the equivalent quadrant from an adjacent face.

  This function returns the quadrant from the adjacent face that
  corresponds to quadrant p on the local face. This algorithm requires
  the the current face index, and the adjacent face index, as well as
  the shared edge number. Note that since the adjacent face is given,
  there is no search over all such adjacent faces.

  input:
  face:      the local face index (p is defined on this face)
  adjacent:  the adjacent face index
  edge:      the edge number shared by face/adjacent
  p:         the quadrant corresponding to the node on the edge

  returns:   the quadrant associated with p on the adjacent face
*/
TMRQuadrant* TMRQuadForest::getEdgeNodeNeighbor( int face,
                                                 int adjacent,
                                                 int edge,
                                                 TMRQuadrant p ){
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Compute the x/y coordinate along the edge
  int32_t ucoord = 0;
  if (edge == 0 || edge == 1){
    ucoord = p.y; 
  }
  else {
    ucoord = p.x;
  }

  // Retrieve the first and second edge numbers
  int n1 = face_conn[4*face + face_to_edge_nodes[edge][0]];
  int n2 = face_conn[4*face + face_to_edge_nodes[edge][1]];

  for ( int j = 0; j < 4; j++ ){
    int nn1 = face_conn[4*adjacent + face_to_edge_nodes[j][0]];
    int nn2 = face_conn[4*adjacent + face_to_edge_nodes[j][1]];

    // Check if the edge is oriented in the same direction or
    // not. Transform the coordinates along the edge accordingly.
    if (n1 == nn1 && n2 == nn2){
      TMRQuadrant neighbor;
      neighbor.level = p.level;
      if (j == 0 || j == 1){
        neighbor.x = hmax*j;
        neighbor.y = ucoord;
      }
      else if (j == 2 || j == 3){
        neighbor.x = ucoord;
        neighbor.y = hmax*(j % 2); 
      }

      // Get the node from the adjacent face
      TMRQuadrantArray *nodes;
      quadtrees[adjacent]->getNodes(&nodes);
     
      const int use_nodes = 1;
      return nodes->contains(&neighbor, use_nodes);
    }
    else if (n1 == nn2 && n2 == nn1){
      // The edges have the opposite orientation
      TMRQuadrant neighbor;
      neighbor.level = p.level;
      if (j == 0 || j == 1){
        neighbor.x = hmax*j;
        neighbor.y = hmax - ucoord;
      }
      else if (j == 2 || j == 3){
        neighbor.x = hmax - ucoord;
        neighbor.y = hmax*(j % 2); 
      }

      // Get the node from the adjacent face
      TMRQuadrantArray *nodes;
      quadtrees[adjacent]->getNodes(&nodes);
     
      const int use_nodes = 1;
      return nodes->contains(&neighbor, use_nodes);
    }
  }

  return NULL;
}

/*
  Balance the forest of quadtrees

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
void TMRQuadForest::balance( int balance_corner ){
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
  delete [] hash;
}

/*
  Duplicate the forest

  This function creates a duplicate representation of the current
  forest. This function copies the global connectivity of the forest
  and copies each individual tree.
*/
TMRQuadForest* TMRQuadForest::duplicate(){
  TMRQuadForest *dup = new TMRQuadForest(comm);
  if (quadtrees){
    // Copy over the connectivity data 
    dup->num_nodes = num_nodes;
    dup->num_edges = num_edges;
    dup->num_faces = num_faces;
    
    // Allocate/copy the face connectivity
    dup->face_conn = new int[ 4*num_faces ];
    memcpy(dup->face_conn, face_conn, 4*num_faces*sizeof(int));

    dup->face_edge_conn = new int[ 4*num_faces ];
    memcpy(dup->face_edge_conn, face_edge_conn, 4*num_faces*sizeof(int));
    
    // Allocate the remaining data
    dup->node_face_ptr = new int[ num_nodes+1 ];
    memcpy(dup->node_face_ptr, node_face_ptr, (num_nodes+1)*sizeof(int));

    dup->node_face_conn = new int[ node_face_ptr[num_nodes] ];
    memcpy(dup->node_face_conn, node_face_conn, 
           node_face_ptr[num_nodes]*sizeof(int));
    
    dup->edge_face_ptr = new int[ num_edges+1 ];
    memcpy(dup->edge_face_ptr, edge_face_ptr, (num_edges+1)*sizeof(int));
    
    dup->edge_face_conn = new int[ edge_face_ptr[num_edges] ];
    memcpy(dup->edge_face_conn, edge_face_conn,
           edge_face_ptr[num_edges]*sizeof(int));

    // Dupn all the quadtrees
    dup->quadtrees = new TMRQuadtree*[ num_faces ];
    for ( int i = 0; i < num_faces; i++ ){
      TMRQuadrantArray *elements;
      quadtrees[i]->getElements(&elements);
      dup->quadtrees[i] = new TMRQuadtree(elements->duplicate());
    }
  }

  return dup;
}

/*
  Coarsen the entire forest

  This function creates a coarsened representation of the current
  forest. This is done by copying the global connectivity of the
  forest and coarsening each individual tree. Note that the resulting
  forest is not necessarily balanced.
*/
TMRQuadForest* TMRQuadForest::coarsen(){
  TMRQuadForest *coarse = new TMRQuadForest(comm);
  if (quadtrees){
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
  Create the nodes from the element mesh

  This function first labels the dependent nodes. These nodes are left
  hanging on an edge in the mesh and depend on independent adjacent
  nodes. While the depedent nodes are determined, we also keep a queue
  of their independent counterparts to generate a data structure that
  gives the dependent nodes in terms of independent nodes.
  
  input:
  order:   the order of the mesh
*/
void TMRQuadForest::createNodes( int order ){
  // Check that the order falls within allowable bounds
  mesh_order = order;
  if (order > 3){ mesh_order = 3; }
  if (order < 2){ mesh_order = 2; }

  // Allocate all possible nodes on all of the trees
  for ( int face = 0; face < num_faces; face++ ){
    quadtrees[face]->createNodes(mesh_order);
  }
  
  // Allocate the space for the dependent node queues for each
  TMRQuadrantQueue **queue = new TMRQuadrantQueue*[ num_faces ];
  TMRQuadrantQueue **indep = new TMRQuadrantQueue*[ num_faces ];
  for ( int face = 0; face < num_faces; face++ ){
    queue[face] = new TMRQuadrantQueue();
    indep[face] = new TMRQuadrantQueue();
  }

  // Set the max element size
  const int hmax = 1 << TMR_MAX_LEVEL;

  // Determine the dependent nodes
  for ( int face = 0; face < num_faces; face++ ){
    // Get the quadrant elements
    TMRQuadrantArray *elements;
    quadtrees[face]->getElements(&elements);

    // Get the elements themselves
    int size;
    TMRQuadrant *array;
    elements->getArray(&array, &size);
    
    for ( int i = 0; i < size; i++ ){
      // Get the side length of the element
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      
      // Check the adjacent elements along each element edge
      for ( int edge = 0; edge < 4; edge++ ){
        // Look for the edge neighbor that is at the next level of
        // refinement from the current element. If this element
        // exists, then we have to use a dependent node.
        TMRQuadrant p = array[i];
        p.level += 1;
        if (edge == 1 || edge == 3){
          p.getSibling(3, &p);
        }
        p.edgeNeighbor(edge, &p);
        
        // Check if this element lies along an edge
        if ((p.x < 0 || p.x >= hmax) ||
            (p.y < 0 || p.y >= hmax)){
          // If the element is on the edge of the tree, then we have
          // to check adjacent trees to see if additional constraints
          // areq required
          addEdgeDependentNodes(face, edge, p, array[i], queue, indep);
        }
        else {
          // If the more-refined element exists, then add the
          // dependent nodes required for compatibility
          const int use_nodes = 0;
          if (elements->contains(&p, use_nodes)){
            if (mesh_order == 2){
              // Find the edge length of the more-refined element
              const int32_t hp = 1 << (TMR_MAX_LEVEL - p.level);

              // This is the one dependent node
              TMRQuadrant node;
              node.tag = -1;
              node.level = p.level;
              if (edge == 0 || edge == 1){
                // Push the dependent node onto the queue
                node.x = array[i].x + edge*h;
                node.y = array[i].y + hp;
                queue[face]->push(&node);
                
                // Push the independent nodes onto the queue
                node.y = array[i].y;
                indep[face]->push(&node);
                node.y = array[i].y + h;
                indep[face]->push(&node);
              }
              else {
                // Push the dependent node onto the queue
                node.x = array[i].x + hp;                
                node.y = array[i].y + (edge % 2)*h;
                queue[face]->push(&node);

                // Push the independent nodes onto the queue
                node.x = array[i].x;
                indep[face]->push(&node);
                node.x = array[i].x + h;
                indep[face]->push(&node);
              }
            }
            else if (mesh_order == 3){
              // Find the edge length of the most-refined level - this
              // is one higher since we are working with the nodal
              // mesh in which every element is refiend once.
              const int32_t hp = 1 << (TMR_MAX_LEVEL - p.level);
              const int32_t hr = 1 << (TMR_MAX_LEVEL - p.level - 1);

              // This there are two dependent nodes for a 3rd order
              // mesh. These correspond to the mid-side nodes of each
              // of the finer two quadrants attached to the larger
              // quadrant.
              TMRQuadrant node;
              node.tag = -1;
              node.level = p.level;
              if (edge == 0 || edge == 1){
                // Push the two dependent nodes onto the queue
                node.x = array[i].x + edge*h;
                node.y = array[i].y + hr;
                queue[face]->push(&node);
                node.y = array[i].y + 3*hr;
                queue[face]->push(&node);

                // Push the nodes for the first constraint
                for ( int k = 0; k < 3; k++ ){
                  node.y = array[i].y + k*hp;
                  indep[face]->push(&node);
                }
                
                // Push the nodes for the second constraint
                for ( int k = 2; k >= 0; k-- ){
                  node.y = array[i].y + k*hp;
                  indep[face]->push(&node);
                }
              }
              else {
                // Push the two dependent nodes onto the queue
                node.x = array[i].x + hr;
                node.y = array[i].y + (edge % 2)*h;
                queue[face]->push(&node);
                node.x = array[i].x + 3*hr;
                queue[face]->push(&node);

                // Push the nodes for the first constraint
                for ( int k = 0; k < 3; k++ ){
                  node.x = array[i].x + k*hp;
                  indep[face]->push(&node);
                }
                
                // Push the nodes for the second constraint
                for ( int k = 2; k >= 0; k-- ){
                  node.x = array[i].x + k*hp;
                  indep[face]->push(&node);
                }
              }
            }
          }
        }
      }
    }
  }

  // Now, go through and label the dependent nodes - first convert
  // all the queues into arrays for easier/consistent access
  TMRQuadrantArray **dep = new TMRQuadrantArray*[ num_faces ];

  for ( int face = 0; face < num_faces; face++ ){
    // Create an array representation of the queue
    dep[face] = queue[face]->toArray();
    delete queue[face];

    // Get the nodes from the quadtree
    TMRQuadrantArray *nodes;
    quadtrees[face]->getNodes(&nodes);

    // Search the node list
    int size;
    TMRQuadrant *array;
    dep[face]->getArray(&array, &size);
    for ( int i = 0; i < size; i++ ){
      const int use_nodes = 1;
      TMRQuadrant *t = nodes->contains(&array[i], use_nodes);
      t->tag = -1;
    }
  }

  // Free the array of queues
  delete [] queue;
  
  // Go through and add the dependent nodes
  num_mesh_nodes = 0;
  num_mesh_dep_nodes = 0;
  num_mesh_elements = 0;
  
  // Now, order the nodes and dependent nodes within the entire
  // quadtree forest
  for ( int face = 0; face < num_faces; face++ ){
    // Count up the number of elements
    int nelems;
    TMRQuadrantArray *elements;
    quadtrees[face]->getElements(&elements);
    elements->getArray(NULL, &nelems);
    num_mesh_elements += nelems;

    // Get the nodes for this quadtree
    int nnodes = 0;
    TMRQuadrantArray *nodes;
    TMRQuadrant *node_array;
    quadtrees[face]->getNodes(&nodes);
    nodes->getArray(&node_array, &nnodes);

    // Find the edge owners
    int edge_owners[4];
    edge_owners[0] = edge_owners[1] = 
      edge_owners[2] = edge_owners[3] = face;

    // Loop over all the faces that connect to each edge
    // and pick the owner - the one that has the lowest 
    // face number
    for ( int j = 0; j < 4; j++ ){
      int edge = face_edge_conn[4*face + j];
      for ( int ip = edge_face_ptr[edge]; 
            ip < edge_face_ptr[edge+1]; ip++ ){
        if (edge_face_conn[ip] < edge_owners[j]){
          edge_owners[j] = edge_face_conn[ip];
        }
      }
    }

    // For each face, cycle through the nodes and assign them, if they
    // are located on another tree, then retrieve the node number from
    // that tree to maintain consistency
    for ( int i = 0; i < nnodes; i++ ){
      if (node_array[i].tag >= 0){
        // First, check if this on an edge
        if (edge_owners[0] < face && node_array[i].x == 0){
          TMRQuadrant *t = getEdgeNodeNeighbor(face, edge_owners[0], 
                                               0, node_array[i]);
          node_array[i].tag = t->tag;
        }
        else if (edge_owners[1] < face && node_array[i].x == hmax){
          TMRQuadrant *t = getEdgeNodeNeighbor(face, edge_owners[1], 
                                               1, node_array[i]);
          node_array[i].tag = t->tag;
        }
        else if (edge_owners[2] < face && node_array[i].y == 0){
          TMRQuadrant *t = getEdgeNodeNeighbor(face, edge_owners[2],
                                               2, node_array[i]);
          node_array[i].tag = t->tag;
        }
        else if (edge_owners[3] < face && node_array[i].y == hmax){
          TMRQuadrant *t = getEdgeNodeNeighbor(face, edge_owners[3],
                                               3, node_array[i]);
          node_array[i].tag = t->tag;
        }
        else {
          node_array[i].tag = num_mesh_nodes;
          num_mesh_nodes++;
        }
      }
      else {
        node_array[i].tag = -(num_mesh_dep_nodes+1);
        num_mesh_dep_nodes++;      
      }
    }
  }

  // Free the dependent node information if it already exists
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
  
  // Allocate the arrays for the dependent nodes
  dep_ptr = new int[ num_mesh_dep_nodes+1 ];
  dep_conn = new int[ mesh_order*num_mesh_dep_nodes ];
  dep_weights = new double[ mesh_order*num_mesh_dep_nodes ];

  // Set the dep_ptr values
  for ( int i = 0; i < num_mesh_dep_nodes+1; i++ ){
    dep_ptr[i] = i*mesh_order;
  }

  // Set the weight values
  double wvals[3];
  wvals[0] = wvals[1] = 0.5;
  wvals[2] = 0.0;

  if (mesh_order == 3){
    wvals[0] = 0.375;
    wvals[1] = 0.75;
    wvals[2] = -0.125;
  }

  // Loop over all the constraint variables and set the dependent
  // weights and variable numbers
  for ( int face = 0; face < num_faces; face++ ){
    // Retrieve the nodes
    TMRQuadrantArray *nodes;
    quadtrees[face]->getNodes(&nodes);

    // Retrieve the dependent nodes
    int size;
    TMRQuadrant *array;
    dep[face]->getArray(&array, &size);

    // Retrieve the edge number
    for ( int i = 0; i < size; i++ ){
      // Get the dependent node number
      const int use_nodes = 1;
      TMRQuadrant *t = nodes->contains(&array[i], use_nodes);
      int node = -t->tag-1;

      // Retrieve the independent nodes
      for ( int k = 0; k < mesh_order; k++ ){
        TMRQuadrant n = indep[face]->pop();
        TMRQuadrant *t = nodes->contains(&n, use_nodes);
        dep_conn[dep_ptr[node]+k] = t->tag;
        dep_weights[dep_ptr[node]+k] = wvals[k];
      }
    }
  }

  // Free things that are no longer required
  for ( int face = 0; face < num_faces; face++ ){
    delete dep[face];
    delete indep[face];
  }
  delete [] indep;
  delete [] dep;
}

/*
  Retrieve the connectivity from the forest
*/
void TMRQuadForest::getMesh( int *nnodes,
                             int *nelems,
                             int **_elem_ptr,
                             int **_elem_conn ){
  // Allocate space for the mesh connectivity
  int *elem_ptr = new int[ num_mesh_elements+1 ];
  int *elem_conn = new int[ mesh_order*mesh_order*num_mesh_elements ];

  // Initialize the first entry of the connectivity
  elem_ptr[0] = 0;

  // Keep a count of the number of elements
  int elem_count = 0;

  // Loop over all the quadtrees, adding to the mesh connectivity
  for ( int face = 0; face < num_faces; face++ ){
    // Add the mesh from this portion of the quadtree
    quadtrees[face]->addMesh(&elem_ptr[elem_count], 
                             &elem_conn[elem_ptr[elem_count]]);

    // Get the number of elements from the quadtree
    elem_count += quadtrees[face]->getNumElements();
  }

  // Set the output
  *nnodes = num_mesh_nodes;
  *nelems = num_mesh_elements;
  *_elem_ptr = elem_ptr;
  *_elem_conn = elem_conn;
}

/*
  Retrieve the dependent nodes that are set internally
*/
void TMRQuadForest::getDependentNodes( int *num_dep_nodes,
                                       const int **_dep_ptr,
                                       const int **_dep_conn,
                                       const double **_dep_weights ){
  if (num_dep_nodes){ *num_dep_nodes = num_mesh_dep_nodes; }
  if (_dep_ptr){ *_dep_ptr = dep_ptr; }
  if (_dep_conn){ *_dep_conn = dep_conn; }
  if (_dep_weights){ *_dep_weights = dep_weights; }
} 
