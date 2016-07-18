#include "TMRQuadForest.h"

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

  // Set the range of nodes
  node_range = NULL;

  // Set the number of elements/dependent nodes
  num_elements = 0;
  num_dep_nodes = 0;

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
  
  // The face partition
  face_owners = NULL;
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
  if (face_owners){ delete [] face_owners; }
  if (node_range){ delete [] node_range; }
}

/*
  Set the connectivity of the faces

  This call is collective on all processors. Every processor must make
  a call with the same connectivity information, otherwise the
  inter-quadtree information will be inconsistent. This sets the
  connectivity and then performs a partition of the mesh using METIS.

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
  if (face_owners){ delete [] face_owners; }

  // Copy over the data locally
  num_nodes = _num_nodes;
  num_faces = _num_faces;
  num_edges = 0;

  // Copy over the face connectivity
  face_conn = new int[ 4*num_faces ];
  memcpy(face_conn, _face_conn, 4*num_faces*sizeof(int));

  // Set up the partition using metis
  face_owners = new int[ num_faces ];

  // Compute the partition on the root processor and
  // broadcast the result to all other processors
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);  
  MPI_Comm_size(comm, &mpi_size);  
  if (mpi_rank == 0){
    // For now.. just perform an assignment
    for ( int i = 0; i < num_faces; i++ ){
      face_owners[i] = i % mpi_size;
    }
  }
  
  // Broadcast the face owners to all processors
  MPI_Bcast(face_owners, num_faces, MPI_INT, 0, comm);

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
      if (quadtrees[i]){ delete quadtrees[i]; }
    }
    delete [] quadtrees;
  }

  // Create the quadtrees
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  quadtrees = new TMRQuadtree*[ num_faces ];
  memset(quadtrees, 0, num_faces*sizeof(TMRQuadtree*));
  for ( int i = 0; i < num_faces; i++ ){
    if (face_owners[i] == mpi_rank){
      quadtrees[i] = new TMRQuadtree(refine_level);
    }
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
      if (quadtrees[i]){ delete quadtrees[i]; }
    }
    delete [] quadtrees;
  }

  // Create a random set of quadtrees
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  quadtrees = new TMRQuadtree*[ num_faces ];
  memset(quadtrees, 0, num_faces*sizeof(TMRQuadtree*));
  for ( int i = 0; i < num_faces; i++ ){
    if (face_owners[i] == mpi_rank){
      quadtrees[i] = new TMRQuadtree(nrand, min_level, max_level);
    }
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
      if (!hash[adjacent]){
        hash[adjacent] = new TMRQuadrantHash();
        queue[adjacent] = new TMRQuadrantQueue();
      }
      
      for ( int j = 0; j < 4; j++ ){
        int nn1 = face_conn[4*adjacent + face_to_edge_nodes[j][0]];
        int nn2 = face_conn[4*adjacent + face_to_edge_nodes[j][1]];

        // Add the quadrant to the list
        if (n1 == nn1 && n2 == nn2){
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          if (j == 0 || j == 1){
            neighbor.x = (hmax - 2*h)*j;
            neighbor.y = ucoord;
          }
          else if (j == 2 || j == 3){
            neighbor.x = ucoord;
            neighbor.y = (hmax - 2*h)*(j % 2); 
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
            neighbor.x = (hmax - 2*h)*j;
            neighbor.y = hmax - 2*h - ucoord;
          }
          else if (j == 2 || j == 3){
            neighbor.x = hmax - 2*h - ucoord;
            neighbor.y = (hmax - 2*h)*(j % 2); 
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
          if (!hash[adjacent]){
            hash[adjacent] = new TMRQuadrantHash();
            queue[adjacent] = new TMRQuadrantQueue();
          }

          // Add the corner node to the list
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          neighbor.x = (hmax - 2*h)*(j%2);
          neighbor.y = (hmax - 2*h)*(j/2);

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
  Balance the forest of quadtrees

  This algorithm uses a hash and a queue to balance the quadtree. For
  each element in the quadtree, we add the neighbors that are required
  to balance to the tree. If the element is not in the hash, we add
  them to a queue, which keeps track of recently added elements. After
  the first pass, the algorithm continues popping elements until the
  all the queues are empty.

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
  // Get the mpi rank
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Create a hash table for the balanced tree
  TMRQuadrantHash **hash = new TMRQuadrantHash*[ num_faces ];
  TMRQuadrantQueue **queue = new TMRQuadrantQueue*[ num_faces ];
  
  // Set the hash tables and queues to NULL
  memset(hash, 0, num_faces*sizeof(TMRQuadrantHash*));
  memset(queue, 0, num_faces*sizeof(TMRQuadrantQueue*));

  for ( int face = 0; face < num_faces; face++ ){
    if (face_owners[face] == mpi_rank){
      // Allocate the hash and queue if they are not already
      // allocated on this processor
      if (!hash[face]){
        hash[face] = new TMRQuadrantHash();
        queue[face] = new TMRQuadrantQueue();
      }

      // Local quadrant data
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
        if (q.level > 1){
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
  }

  int queue_length_flag = 1;
  while (queue_length_flag){
    // Now continue until the queue of added quadrants is
    // empty. At each iteration, pop an quadrant and add 
    // its neighbours until nothing new is added. This code
    // handles the propagation of quadrants to adjacent quadrants.
    for ( int face = 0; face < num_faces; face++ ){
      if (queue[face]){
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
    }
      
    queue_length_flag = 0;
    for ( int face = 0; face < num_faces; face++ ){
      if (queue[face] &&
          queue[face]->length() > 0){
        queue_length_flag = 1;
        break;
      }
    }
  }

  // Free the queues - they are no longer required
  for ( int face = 0; face < num_faces; face++ ){
    if (face_owners[face] != mpi_rank && queue[face]){ 
      delete queue[face];
      queue[face] = NULL;
    }
  }

  // Now everything is locally balanced - all the elements 
  // on the current processor are balanced with all the other
  // element on the current processor, but nothing is 
  // inter-processor balanced yet.

  // The number of quadrants that will be sent from this processor
  // to all other processors in the communicator
  int *quad_counts = new int[ mpi_size ];
  memset(quad_counts, 0, mpi_size*sizeof(int));

  // Create a list of new queues - these will only be allocated
  // as required
  TMRQuadrantQueue **send_queues = new TMRQuadrantQueue*[ mpi_size ];
  memset(send_queues, 0, mpi_size*sizeof(TMRQuadrantQueue*));

  for ( int face = 0; face < num_faces; face++ ){
    if (hash[face] && face_owners[face] != mpi_rank){
      // Create a sorted list of local the 0-child quadrants.
      // This can be further reduced to limit the amount of 
      // memory passed between processors
      TMRQuadrantArray *elems0 = hash[face]->toArray();
      elems0->sort();

      // Get the array of 0-quadrants 
      int size;
      TMRQuadrant *array;
      elems0->getArray(&array, &size);

      // Check that the destination queue has been allocated
      int dest_rank = face_owners[face];
      if (!send_queues[dest_rank]){
        send_queues[dest_rank] = new TMRQuadrantQueue();
      }
      
      // Add the quadrants to their destination queue
      TMRQuadrantQueue *dest = send_queues[dest_rank];
      
      if (size > 0){
        // Get the parent of the quadrant
        TMRQuadrant p, s = array[0];
        s.parent(&p);

        // Loop over all 
        for ( int i = 0; i < size; i++ ){
          if (!p.contains(&array[i])){
            s.tag = face;
            dest->push(&s);
          }
          // Get the next quadrant and find its parent
          s = array[i];
          s.parent(&p);
        }

        // Push the last quadrant onto the queue
        s.tag = face;
        dest->push(&s);
      }
      
      // Update the number of quadrants
      quad_counts[dest_rank] = dest->length();
      
      // Free the elements and the hash table
      delete elems0;
      delete hash[face];
      hash[face] = NULL;
    }
  }

  // Now distribute the quadrants to their destination quadtrees and
  // balance the corresponding quadtrees including the new elements.
  int *quad_recv_counts = new int[ mpi_size ];
  MPI_Alltoall(quad_counts, 1, MPI_INT,
               quad_recv_counts, 1, MPI_INT, comm);

  // Add the recieved quadrants into the local queues
  TMRQuadrantArray **send_arrays = new TMRQuadrantArray*[ mpi_size ];
  TMRQuadrantArray **recv_arrays = new TMRQuadrantArray*[ mpi_size ];
  memset(send_arrays, 0, mpi_size*sizeof(TMRQuadrantArray*));
  memset(recv_arrays, 0, mpi_size*sizeof(TMRQuadrantArray*));

  // Allocate space for the requests
  MPI_Request *send_request = new MPI_Request[ mpi_size ];
  MPI_Request *recv_request = new MPI_Request[ mpi_size ];

  // Set the recv requests to NULL 
  for ( int i = 0; i < mpi_size; i++ ){
    recv_request[i] = MPI_REQUEST_NULL;
  }

  // Loop over all the ranks and send 
  int nsend = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (send_queues[i]){
      // Conver the quadrant array to a list
      send_arrays[i] = send_queues[i]->toArray();
    
      // Delete the send queue 
      delete send_queues[i];

      // Get the quadrant array
      int size;
      TMRQuadrant *array;
      send_arrays[i]->getArray(&array, &size);

      // Post the send to the destination
      MPI_Isend(array, size, TMRQuadrant_MPI_type, 
                i, 0, comm, &send_request[nsend]);
      nsend++;
    }
  }

  delete [] send_queues;

  // Loop over the recieve calls
  int nrecv = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (quad_recv_counts[i] > 0){
      int size = quad_recv_counts[i];
      TMRQuadrant *array = new TMRQuadrant[ size ];
      recv_arrays[i] = new TMRQuadrantArray(array, size);
      
      // Post the recieve
      MPI_Irecv(array, size, TMRQuadrant_MPI_type,
                i, 0, comm, &recv_request[i]);
      nrecv++;
    }
  }

  // Wait for any recieve to complete
  for ( int i = 0; i < nrecv; i++ ){
    int index = 0;
    MPI_Waitany(mpi_size, recv_request, &index, MPI_STATUS_IGNORE);

    // Push the quadrants on to their corresponding queues
    int size;
    TMRQuadrant *array;
    recv_arrays[index]->getArray(&array, &size);
    for ( int j = 0; j < size; j++ ){
      int face = array[j].tag;
      if (hash[face]->addQuadrant(&array[j])){
        queue[face]->push(&array[j]);
      }
    }

    // Free the data
    delete recv_arrays[index];
  }

  delete [] recv_arrays;

  // Wait for any remaining sends to complete
  MPI_Waitall(nsend, send_request, MPI_STATUSES_IGNORE);

  // Free the remaining send data
  for ( int i = 0; i < mpi_size; i++ ){
    if (send_arrays[i]){ delete send_arrays[i]; } 
  }
  delete [] send_arrays;

  // Free other data associated with the parallel communication
  delete [] quad_counts;
  delete [] quad_recv_counts;
  delete [] send_request;
  delete [] recv_request;

  // Now all the received quadrants will balance the tree locally
  // without having to worry about off-processor quadrants.
  for ( int face = 0; face < num_faces; face++ ){
    if (queue[face]){
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
            }
          }
        }
      }
    }
  }

  // Free the remaining queues - these will be the ones that are
  // owned by this processor.
  for ( int face = 0; face < num_faces; face++ ){
    if (queue[face]){ delete queue[face]; }
  }
  delete [] queue;

  // Convert the hash tables back to the elements
  for ( int face = 0; face < num_faces; face++ ){
    if (face_owners[face] == mpi_rank){
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

    dup->face_owners = new int[ num_faces ];
    memcpy(dup->face_owners, face_owners, num_faces*sizeof(int));

    // Duplicate all the quadtrees
    dup->quadtrees = new TMRQuadtree*[ num_faces ];
    memset(dup->quadtrees, 0, num_faces*sizeof(TMRQuadtree*));
    for ( int i = 0; i < num_faces; i++ ){
      if (quadtrees[i]){
        TMRQuadrantArray *elements;
        quadtrees[i]->getElements(&elements);
        dup->quadtrees[i] = new TMRQuadtree(elements->duplicate());
      }
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

    // Copy over the face ownership data
    coarse->face_owners = new int[ num_faces ];
    memcpy(coarse->face_owners, face_owners, num_faces*sizeof(int));

    // Coarsen all the quadtrees
    coarse->quadtrees = new TMRQuadtree*[ num_faces ];
    memset(coarse->quadtrees, 0, num_faces*sizeof(TMRQuadtree*));
    for ( int i = 0; i < num_faces; i++ ){
      if (quadtrees[i]){
        coarse->quadtrees[i] = quadtrees[i]->coarsen();
      }
    }
  }

  return coarse;
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
                                           TMRQuadrant source ){
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

      // Get the quadrant nodes
      TMRQuadrantArray *nodes;
      quadtrees[adjacent]->getNodes(&nodes);

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
                TMRQuadrant *t = nodes->contains(&node, 1);
                t->tag = -1;
              }
              else if (j == 2 || j == 3){
                // Add the dependent node
                node.x = uts + hp;
                node.y = hmax*(j % 2);
                TMRQuadrant *t = nodes->contains(&node, 1);
                t->tag = -2;
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
                TMRQuadrant *t = nodes->contains(&node, 1);
                t->tag = -1;

                node.y = uts + 3*hr;
                t = nodes->contains(&node, 1);
                t->tag = -2;
              }
              else if (j == 2 || j == 3){
                // Push the two dependent nodes onto the queue
                node.x = uts + hr;
                node.y = hmax*(j % 2);
                TMRQuadrant *t = nodes->contains(&node, 1);
                t->tag = -3;
                
                node.x = uts + 3*hr;
                t = nodes->contains(&node, 1);
                t->tag = -4;
              }
            }
          }
        }
      }
    }
  }
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

  // Get the MPI information
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Find the face number corresponding to the owner for 
  // each edge and corner so that we know who should be 
  // ordering what!
  int *edge_face_owners = new int[ num_edges ];
  int *node_face_owners = new int[ num_nodes ];

  // The owner is chosen as the connecting face with the
  // lowest face number
  for ( int edge = 0; edge < num_edges; edge++ ){
    edge_face_owners[edge] = num_faces;

    int ipend = edge_face_ptr[edge+1];
    for ( int ip = edge_face_ptr[edge]; ip < ipend; ip++ ){
      if (edge_face_conn[ip] < edge_face_owners[edge]){
        edge_face_owners[edge] = edge_face_conn[ip];
      }
    }
  }

  // Find the face numbers that own each node
  for ( int node = 0; node < num_nodes; node++ ){
    node_face_owners[node] = num_faces;
    
    int ipend = node_face_ptr[node+1];
    for ( int ip = node_face_ptr[node]; ip < ipend; ip++ ){
      if (node_face_conn[ip] < node_face_owners[node]){
        node_face_owners[node] = node_face_conn[ip];
      }
    }
  }

  // Exchange the quadrants that are located along shared
  // edges/corners so that we can label all the dependent nodes
  // that are owned by this face.
  TMRQuadrantQueue **queues = new TMRQuadrantQueue*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      queues[k] = new TMRQuadrantQueue();
    }
    else {
      queues[k] = NULL;
    }
  }

  for ( int edge = 0; edge < num_edges; edge++ ){
    // Set the destination rank for the 
    int dest_face = edge_face_owners[edge];
    int dest_rank = face_owners[dest_face];

    // Loop through all faces that contribute to this edge
    for ( int ip = edge_face_ptr[edge];
          ip < edge_face_ptr[edge+1]; ip++ ){
      // Get the face corresponding to this edge
      int face = edge_face_conn[ip];
      
      // Check if this face is locally owned, and if
      // we need to send some of its quadrants to another processor
      if ((face_owners[face] == mpi_rank) && 
          (dest_rank != mpi_rank)){
        // Get the element array
        TMRQuadrantArray *elements;
        quadtrees[face]->getElements(&elements);

        // Get the actual quadrant array
        int size;
        TMRQuadrant *array;
        elements->getArray(&array, &size);

        // Find out which edge we need to send
        int edge_num = 0;
        for ( ; edge_num < 4; edge_num++ ){
          if (face_edge_conn[4*face + edge_num] == edge){
            break;
          }
        }

        // Loop over all teh elements and check 
        for ( int i = 0; i < size; i++ ){
          const int32_t hmax = 1 << TMR_MAX_LEVEL; 
          const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
        
          // Check if this lies along an edge
          if ((edge_num == 0 && array[i].x == 0) ||
              (edge_num == 1 && array[i].x + h == hmax) ||
              (edge_num == 2 && array[i].y == 0) ||
              (edge_num == 3 && array[i].y + h == hmax)){
            TMRQuadrant q = array[i];
            q.tag = face;
            queues[dest_rank]->push(&q);
          }
        }
      }
    }
  }

  // Allocate the requests
  MPI_Request *send_requests = new MPI_Request[ mpi_size ];

  // Send the quadrants to the other trees
  int nsends = 0;
  TMRQuadrantArray **arrays = new TMRQuadrantArray*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      // Create the arrays
      arrays[k] = queues[k]->toArray();
      delete queues[k];

      // Get the array of quadrants
      int size;
      TMRQuadrant *array;
      arrays[k]->getArray(&array, &size);

      // Set the array of quadrants to their destination
      MPI_Isend(array, size, TMRQuadrant_MPI_type, 
                k, 0, comm, &send_requests[nsends]);
      nsends++;
    }
    else {
      arrays[k] = NULL;
    }
  }
  
  // Free the array of queues
  delete [] queues;

  // Allocate the queues
  TMRQuadrantQueue **qtrees = new TMRQuadrantQueue*[ num_faces ];
  memset(qtrees, 0, num_faces*sizeof(TMRQuadrantQueue*));

  // Receive the arrays of incoming quadrants
  for ( int k = 0; k < mpi_size-1; k++ ){
    // Probe the recieved messages
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    // Retrieve the size and information for the incoming message
    int source = status.MPI_SOURCE;
    int tag = status.MPI_TAG;
    int size = 0;
    MPI_Get_count(&status, TMRQuadrant_MPI_type, &size);

    // Allocate the incoming array
    TMRQuadrant *array = new TMRQuadrant[ size ];
    MPI_Recv(array, size, TMRQuadrant_MPI_type,
             source, tag, comm, MPI_STATUS_IGNORE);

    // Push the quadrants into their corresponding trees
    for ( int i = 0; i < size; i++ ){
      int face = array[i].tag;
      if (!qtrees[face]){
        qtrees[face] = new TMRQuadrantQueue();
      }
      qtrees[face]->push(&array[i]);
    }

    delete [] array;
  }

  // Now that the queues are completed allocate the faces
  for ( int face = 0; face < num_faces; face++ ){
    if (qtrees[face]){
      quadtrees[face] = new TMRQuadtree(qtrees[face]->toArray());
      delete qtrees[face];
    }
  }
  delete [] qtrees;

  // Allocate all possible nodes on all of the trees
  for ( int face = 0; face < num_faces; face++ ){
    if (quadtrees[face]){
      quadtrees[face]->createNodes(mesh_order);
    }
  }

  // Wait for all the sends to complete
  MPI_Waitall(nsends, send_requests, MPI_STATUSES_IGNORE);
  
  // Now free the arrays
  for ( int k = 0; k < mpi_size; k++ ){
    if (arrays[k]){
      delete arrays[k];
    }
  }
  delete [] arrays;

  // Determine the dependent nodes for each face without labeling
  // the dependent nodes on the edges yet.
  for ( int face = 0; face < num_faces; face++ ){
    if (quadtrees[face]){
      // Get the quadrant elements
      TMRQuadrantArray *elements, *nodes;
      quadtrees[face]->getElements(&elements);
      quadtrees[face]->getNodes(&nodes);
      
      // Get the elements themselves
      int size;
      TMRQuadrant *array;
      elements->getArray(&array, &size);

      // Set flags to indicate whether this face owns
      // each of its edges
      int is_edge_owner[4] = {0, 0, 0, 0};
      for ( int j = 0; j < 4; j++ ){
        int edge = face_edge_conn[4*face + j];
        if (mpi_rank == face_owners[edge_face_owners[edge]]){
          is_edge_owner[j] = 1;
        }
      }

      for ( int i = 0; i < size; i++ ){
        // Get the side length of the element
        const int32_t hmax = 1 << TMR_MAX_LEVEL;
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
          if ((p.x < 0 && is_edge_owner[0]) || 
              (p.x >= hmax && is_edge_owner[1]) ||
              (p.y < 0 && is_edge_owner[2]) || 
              (p.y >= hmax && is_edge_owner[3])){
            // If the element is on the edge of the tree, then we have
            // to check adjacent trees to see if additional constraints
            // areq required
            addEdgeDependentNodes(face, edge, p, array[i]);
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
                if (edge == 0 || edge == 1){
                  node.tag = -1;
                  node.x = array[i].x + edge*h;
                  node.y = array[i].y + hp;
                }
                else {
                  node.tag = -2;
                  node.x = array[i].x + hp;                
                  node.y = array[i].y + (edge % 2)*h;
                }
                
                // Search for dependent node and label it
                const int use_node_search = 1;
                TMRQuadrant *t = nodes->contains(&node, use_node_search);
                t->tag = node.tag;
              }
              else if (mesh_order == 3){
                // Find the edge length of the most-refined level - this
                // is one higher since we are working with the nodal
                // mesh in which every element is refiend once.
                const int32_t hp = 1 << (TMR_MAX_LEVEL - array[i].level - 1);
                const int32_t hr = 1 << (TMR_MAX_LEVEL - array[i].level - 2);

                // This there are two dependent nodes for a 3rd order
                // mesh. These correspond to the mid-side nodes of each
                // of the finer two quadrants attached to the larger
                // quadrant.
                TMRQuadrant node1, node2;
                if (edge == 0 || edge == 1){
                  node1.tag = -1;
                  node1.x = array[i].x + edge*h;
                  node1.y = array[i].y + hr;
                  
                  node2.tag = -2;
                  node2.x = array[i].x + edge*h;
                  node2.y = array[i].y + 3*hr;
                }
                else {
                  node1.tag = -3;
                  node1.x = array[i].x + hr;
                  node1.y = array[i].y + (edge % 2)*h;

                  node2.tag = -4;
                  node2.x = array[i].x + 3*hr;
                  node2.y = array[i].y + (edge % 2)*h;
                }

                // Search for dependent node and label it
                const int use_node_search = 1;
                TMRQuadrant *t;
                t = nodes->contains(&node1, use_node_search);
                t->tag = node1.tag;

                t = nodes->contains(&node2, use_node_search);
                t->tag = node2.tag;
              }
            }
          }
        }
      }
    }
  }

  // Count up all the local variables
  int nlocal = 0;

  // Count up the number of nodes owned on each face
  for ( int face = 0; face < num_faces; face++ ){
    if (mpi_rank == face_owners[face]){
      // Set flags to indicate whether this face owns
      // each of its edges
      int is_edge_owner[4] = {0, 0, 0, 0};
      for ( int j = 0; j < 4; j++ ){
        int edge = face_edge_conn[4*face + j];
        int face = edge_face_owners[edge];
        if (mpi_rank == face_owners[face]){
          is_edge_owner[j] = 1;
        }
      }

      // Count up the number of locally owned corner nodes
      for ( int k = 0; k < 4; k++ ){
        int node = face_conn[4*face + k];
        if (node_face_owners[node] == face){
          nlocal++;
        }
      }

      // Get the quadrant array of nodes
      TMRQuadrantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Extract the array of quadrants
      int size;
      TMRQuadrant *array;
      nodes->getArray(&array, &size);

      // Count up the number of interior nodes
      const int32_t hmax = 1 << TMR_MAX_LEVEL;
      for ( int i = 0; i < size; i++ ){
        if ((array[i].x > 0 && array[i].x < hmax) &&
            (array[i].y > 0 && array[i].y < hmax) && 
            array[i].tag >= 0){
          nlocal++;
        }
        else if (((array[i].y > 0 && array[i].y < hmax && 
                   (is_edge_owner[0] || is_edge_owner[1])) ||
                  (array[i].x > 0 && array[i].x < hmax && 
                   (is_edge_owner[2] || is_edge_owner[3]))) &&
                 array[i].tag >= 0){
          nlocal++;
        }
      }
    }
  }

  // Gather the local variable counts from each processor
  node_range = new int[ mpi_size+1 ];
  memset(node_range, 0, (mpi_size+1)*sizeof(int));
  MPI_Allgather(&nlocal, 1, MPI_INT, &node_range[1], 1, MPI_INT, comm);
 
  // Create the offsets to each node
  for ( int i = 0; i < mpi_size; i++ ){
    node_range[i+1] += node_range[i];
  }
  
  // Set the global ordering of the local faces, edges, and nodes
  // using the offsets computed above
  // Set the node number offset from the range of owned nodes
  int node_num = node_range[mpi_rank];
  for ( int face = 0; face < num_faces; face++ ){
    if (face_owners[face] == mpi_rank){
      // Set a flag to check if this processor owns each
      // of the corner nodes
      int corners[4];
      for ( int k = 0; k < 4; k++ ){
        int node = face_conn[4*face + k];
        corners[k] = (node_face_owners[node] == face);
      }

      // Set a flag to check if this processor owns each
      // of the edges
      int edges[4];
      for ( int k = 0; k < 4; k++ ){
        int edge = face_edge_conn[4*face + k];
        edges[k] = (edge_face_owners[edge] == face);
      }

      // Now order the face. Get the quadrant array of nodes
      TMRQuadrantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Extract the array of quadrants
      int size;
      TMRQuadrant *array;
      nodes->getArray(&array, &size);

      // Count up the number nodes that are locally owned
      // in the interior, edges or corners
      const int32_t hmax = 1 << TMR_MAX_LEVEL;

      for ( int i = 0; i < size; i++ ){
        // Get the x/y coordinates for easier access
        int32_t x = array[i].x;
        int32_t y = array[i].y;

        if (array[i].tag >= 0){
          // Check whether the node is in the interior
          if ((x > 0 && x < hmax) &&
              (y > 0 && y < hmax)){
            array[i].tag = node_num;
            node_num++;
          }
          // Check if the node lies on an edge
          else if ((edges[0] && (x == 0    && (y > 0 && y < hmax))) ||
                   (edges[1] && (x == hmax && (y > 0 && y < hmax))) ||
                   (edges[2] && (y == 0    && (x > 0 && x < hmax))) ||
                   (edges[3] && (y == hmax && (x > 0 && x < hmax)))){
            array[i].tag = node_num;
            node_num++;
          }
          else if ((corners[0] && x == 0 && y == 0) ||
                   (corners[1] && x == hmax && y == 0) ||
                   (corners[2] && x == 0 && y == hmax) ||
                   (corners[3] && x == hmax && y == hmax)){
            array[i].tag = node_num;
            node_num++;
          }
          else if (x == 0 || x == hmax){
            // Label the node as a dependent node
            array[i].tag = -1;
          }
          else if (y == 0 || y == hmax){
            // Label the node as a dependent node
            array[i].tag = -2;
          }
        }
      }
    }
  }

  // Send all the nodes that have a non-negative node number
  queues = new TMRQuadrantQueue*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      queues[k] = new TMRQuadrantQueue();
    }
    else {
      queues[k] = NULL;
    }
  }

  for ( int edge = 0, nsends = 0; edge < num_edges; edge++ ){
    // Get the owner (destination) rank for this edge
    int face = edge_face_owners[edge];
    int edge_owner = face_owners[face];

    // If this edge actually needs to be sent anywhere
    if (mpi_rank == edge_owner){
      // Find the local edge number for the source edge
      int edge_num = 0;
      for ( ; edge_num < 4; edge_num++ ){
        if (face_edge_conn[4*face + edge_num] == edge){
          break;
        }
      }

      // Get the edge nodes to determine the edge orientation
      int n1 = face_conn[4*face + face_to_edge_nodes[edge_num][0]];
      int n2 = face_conn[4*face + face_to_edge_nodes[edge_num][1]];

      // Loop over the edges 
      int ipend = edge_face_ptr[edge+1];
      for ( int ip = edge_face_ptr[edge]; ip < ipend; ip++ ){
        // Get the destination face/numbers
        int dest_face = edge_face_conn[ip];
        int dest_face_owner = face_owners[dest_face];

        // The destination face and the current face are the same
        if (dest_face == face){ 
          continue;
        }

        // Find the destination edge number
        int dest_edge_num = 0;
        for ( ; dest_edge_num < 4; dest_edge_num++ ){
          if (face_edge_conn[4*dest_face + dest_edge_num] == edge){
            break;
          }
        }

        // Get the edge nodes to determine the edge orientation
        int d1 = face_conn[4*dest_face + face_to_edge_nodes[dest_edge_num][0]];
        int d2 = face_conn[4*dest_face + face_to_edge_nodes[dest_edge_num][1]];

        // Loop over the nodes finding those on the edge
        TMRQuadrantArray *nodes;
        quadtrees[face]->getNodes(&nodes);

        // Get the quadrant arrays
        int size;
        TMRQuadrant *array;
        nodes->getArray(&array, &size);
        
        // Loop over all the nodes on the face and find those
        // nodes that lie on the edge
        for ( int i = 0; i < size; i++ ){
          const int32_t hmax = 1 << TMR_MAX_LEVEL; 

          int32_t u = -1;
          if ((edge_num == 0 || edge_num == 1) &&
              (array[i].x == edge_num*hmax) && 
              array[i].tag >= 0){
            u = array[i].y;
          }
          else if ((edge_num == 2 || edge_num == 3) &&
                   (array[i].y == (edge_num % 2)*hmax) && 
                   array[i].tag >= 0){
            u = array[i].x;
          }

          // If the new node actually exists
          if (u >= 0){
            TMRQuadrant q;

            // Use the level to indicate the face number
            q.level = dest_face;

            // Set the node number on the destination face
            q.tag = array[i].tag;
              
            // Transform the nodal coordinates to the new
            // reference system on the opposite quadtree
            if (dest_edge_num == 0 || dest_edge_num == 1){
              q.x = hmax*dest_edge_num;
              if (n1 == d1 && n2 == d2){
                q.y = u;
              }
              else {
                q.y = hmax - u; 
              }
            }
            else if (dest_edge_num == 2 || dest_edge_num == 3){
              q.y = hmax*(dest_edge_num % 2);
              if (n1 == d1 && n2 == d2){
                q.x = u;
              }
              else {
                q.x = hmax - u; 
              }
            }

            if (dest_face_owner == mpi_rank){
              // Get the quadtree face
              TMRQuadrantArray *dest_nodes;
              quadtrees[dest_face]->getNodes(&dest_nodes);
              
              // Get the node and set the tag value
              const int use_nodes = 1;
              TMRQuadrant *t = dest_nodes->contains(&q, use_nodes);
              t->tag = q.tag;
            }
            else {
              // Push the face onto the queue to be sent later
              queues[dest_face_owner]->push(&q);
            }
          }
        }
      }
    }
  }

  // Convert the queues to arrays and send them to their
  // destinations
  nsends = 0;
  arrays = new TMRQuadrantArray*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      // Create the arrays
      arrays[k] = queues[k]->toArray();
      delete queues[k];

      // Get the array of quadrants
      int size;
      TMRQuadrant *array;
      arrays[k]->getArray(&array, &size);

      // Set the array of quadrants to their destination
      MPI_Isend(array, size, TMRQuadrant_MPI_type, 
                k, 0, comm, &send_requests[nsends]);
      nsends++;
    }
    else {
      arrays[k] = NULL;
    }
  }

  // Receive the arrays of incoming quadrants
  for ( int k = 0; k < mpi_size-1; k++ ){
    // Probe the recieved messages
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    // Retrieve the size and information for the incoming message
    int source = status.MPI_SOURCE;
    int tag = status.MPI_TAG;
    int size = 0;
    MPI_Get_count(&status, TMRQuadrant_MPI_type, &size);

    // Allocate the incoming array
    TMRQuadrant *array = new TMRQuadrant[ size ];
    MPI_Recv(array, size, TMRQuadrant_MPI_type,
             source, tag, comm, MPI_STATUS_IGNORE);

    // Push the quadrants into their corresponding trees
    for ( int i = 0; i < size; i++ ){
      int face = array[i].level;

      // Get the face nodes
      TMRQuadrantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Search the nodes
      const int use_nodes = 1;
      TMRQuadrant *t = nodes->contains(&array[i], use_nodes);
      t->tag = array[i].tag;
    }

    delete [] array;
  }

  // Wait for all the sends to complete
  MPI_Waitall(nsends, send_requests, MPI_STATUSES_IGNORE);
  delete [] send_requests;
  
  // Free the quadtrees that are no longer required
  for ( int face = 0; face < num_faces; face++ ){
    if (quadtrees[face] && (face_owners[face] != mpi_rank)){
      delete quadtrees[face];
      quadtrees[face] = NULL;
    }
  }

  // Now free the arrays
  for ( int k = 0; k < mpi_size; k++ ){
    if (arrays[k]){
      delete arrays[k];
    }
  }
  delete [] arrays;

  // Label the dependent nodes
  num_elements = 0;
  num_dep_nodes = 0;
  for ( int face = 0; face < num_faces; face++ ){
    if (face_owners[face] == mpi_rank){
      // Get the nodes quadrants
      TMRQuadrantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Get the node array
      int size;
      TMRQuadrant *array;
      nodes->getArray(&array, &size);
      
      // Count up the number of dependent nodes
      for ( int i = 0; i < size; i++ ){
        if (array[i].tag < 0){
          array[i].tag = -(num_dep_nodes+1);
          num_dep_nodes++;
        }
      }
      
      // Order the elements for this quadrant
      TMRQuadrantArray *elements;
      quadtrees[face]->getElements(&elements);
      elements->getArray(&array, &size);
      for ( int i = 0; i < size; i++ ){
        array[i].tag = num_elements;
        num_elements++;
      }
    }
  }

  MPI_Barrier(comm);

  delete [] edge_face_owners;
  delete [] node_face_owners;
}

/*
  Retrieve the connectivity from the forest
*/
void TMRQuadForest::getMesh( int *nnodes,
                             int *ndep_nodes,
                             int *nelems,
                             int **_elem_conn,
                             int **_dep_conn,
                             double **_dep_weights ){
  // Get the MPI rank
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Allocate space for the mesh connectivity
  int *elem_conn = new int[ mesh_order*mesh_order*num_elements ];
  int *dep_conn = new int[ mesh_order*num_dep_nodes ];
  double *dep_weights = new double[ mesh_order*num_dep_nodes ];

  // Loop over all the quadtrees, adding to the mesh connectivity
  for ( int face = 0; face < num_faces; face++ ){
    if (face_owners[face] == mpi_rank){
      // Add the mesh from this portion of the quadtree
      quadtrees[face]->addMesh(elem_conn, 
                               dep_conn, dep_weights);
    }
  }

  // Set the output
  *nnodes = node_range[mpi_rank+1] - node_range[mpi_rank];
  *ndep_nodes = num_dep_nodes;
  *nelems = num_elements;
  *_elem_conn = elem_conn;
  *_dep_conn = dep_conn;
  *_dep_weights = dep_weights;
}

