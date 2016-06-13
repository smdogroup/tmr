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
  Go through the locally owned edges and tag the dependent nodes
*/
int TMRQuadForest::labelDependentNodes( int edge,
                                        TMRQuadrantArray **edge_nodes,
                                        int count ){
  int num_indep = 0;

  // The nodes are independent only if all the edges have the same
  // node number, otherwise it is a dependent node
  if (count == 1){
    int size;
    TMRQuadrant *array;
    edge_nodes[0]->getArray(&array, &size);

    // Label all the nodes as independent
    for ( ; num_indep < size; num_indep++ ){
      array[num_indep].tag = 1;
    }
  }
  else {
    // Get the sizes, the relative directions
    // and the edges that these correspond to
    TMRQuadrant *array[25];
    int size[25];

    // Retrieve the quadrant arrays
    for ( int j = 0; j < count; j++ ){
      edge_nodes[j]->getArray(&array[j], &size[j]);
    }

    // Set up the information about the orientation of the edges
    int dir[25], loc[25], edge_nums[25];

    // The node numbers for the first edge
    int n1 = 0, n2 = 0;

    // Determine the edge directions
    for ( int j = 0, ip = edge_face_ptr[edge];
          ip < edge_face_ptr[edge+1]; j++, ip++ ){
      int face = edge_face_conn[ip];
      // Find the edge number
      for ( int k = 0; k < 4; k++ ){
        if (face_edge_conn[4*face + k] == edge){
          edge_nums[j] = k;
          break;
        }
      }

      // Determine the direction and the parametric location
      dir[j] = 0, loc[j] = 0;
      int d1 = face_conn[4*face + face_to_edge_nodes[edge_nums[j]][0]]; 
      int d2 = face_conn[4*face + face_to_edge_nodes[edge_nums[j]][1]]; 
        
      if (j == 0){
        n1 = d1;  n2 = d2;
      }
      else {
        if (n1 == d1 && n2 == d2){
          dir[j] = 1;
          loc[j] = 0;
        }
        else {
          dir[j] = -1;
          loc[j] = size[j]-1;
        }
      }
    }

    // Keep track of the position along the edges
    int32_t u[25];

    // Set the maximum side length
    const int32_t hmax = 1 << TMR_MAX_LEVEL;

    // Loop over all the edges
    for ( int i = 0; i < size[0]; i++ ){
      // Pick the corresponding x/y coordinate
      if (edge_nums[0] == 0 || edge_nums[0] == 1){
        u[0] = array[0][i].y;
      }
      else {
        u[0] = array[0][i].x;
      }

      // Now cycle through the remaining edges
      for ( int j = 1; j < count; j++ ){
        // Determine whether to keep on incrementing the
        // parameter along the edge
        int flag = ((dir[j] > 0 && (loc[j] < size[j])) ||
                    (dir[j] < 0 && (loc[j] >= 0)));
        while (flag){
          if (edge_nums[j] == 0 || edge_nums[j] == 1){
            if (dir[j] > 0){
              u[j] = array[j][loc[j]].y;
            }
            else {
              u[j] = hmax - array[j][loc[j]].y;
            }
          }
          else {
            if (dir[j] > 0){
              u[j] = array[j][loc[j]].x;
            }
            else {
              u[j] = hmax - array[j][loc[j]].x;
            }
          }
          
          // If the value of u[j] is less and u[0] keep 
          // incrementing the indicies
          if (u[j] < u[0]){
            array[j][loc[j]].tag = -1;
            loc[j] += dir[j];
          }
          else {
            flag = 0;
          }
        }
      }
          
      // Check if all the nodes are equal
      int equal = 1;
      for ( int j = 1; j < count; j++ ){
        if (u[j] != u[0]){
          equal = 0;
          break;
        }
      }
          
      // If the nodes are not all equal, label the root
      // list as a dependent node
      if (equal){
        // The nodes are equal
        array[0][i].tag = 1+num_indep;
        for ( int j = 1; j < count; j++ ){
          array[j][loc[j]].tag = 1+num_indep;
          loc[j] += dir[j];
        }

        // Increase the number of independent nodes
        num_indep++;
      }
      else {
        array[0][i].tag = -1;
      }
    }
  }

  return num_indep;
}

/*
  Get the node numbers from the owning face and set them for an edge
  along the destination face.  Note that the orientations of these two
  edges may not be the same.

  input:
  face:        the face owner
  edge:        the global edge number 
  dest_face:   the destination face for the nodes shared along the edge
  
  output:
  edge_nodes:  the edge node numbers
*/
void TMRQuadForest::getEdgeNodeNums( int face,
                                     int edge,
                                     int dest_face, 
                                     TMRQuadrantArray *edge_nodes ){
  // Determine the edge number and the node numbers for
  // the owner face - we'll extract the node numbers
  // from this face
  int edge_num = 0;
  for ( int k = 0; k < 4; k++ ){
    if (face_edge_conn[4*face + k] == edge){
      edge_num = k;
      break;
    }
  }

  int n1 = face_conn[4*face + face_to_edge_nodes[edge_num][0]];
  int n2 = face_conn[4*face + face_to_edge_nodes[edge_num][1]];

  // Determine the edge number and node numbers for
  // the destination face - this is the order used
  // in the edge_nodes array
  int dest_edge_num = 0;
  for ( int k = 0; k < 4; k++ ){
    if (face_edge_conn[4*dest_face + k] == edge){
      dest_edge_num = k;
      break;
    }
  }

  int d1 = face_conn[4*dest_face + face_to_edge_nodes[dest_edge_num][0]];
  int d2 = face_conn[4*dest_face + face_to_edge_nodes[dest_edge_num][1]];
  
  // Get the node array
  TMRQuadrantArray *nodes;
  quadtrees[face]->getNodes(&nodes);

  // Get the quadrant arrays
  int size;
  TMRQuadrant *array;
  nodes->getArray(&array, &size);

  // Retrieve the array of edges
  int edge_size;
  TMRQuadrant *edge_array;
  edge_nodes->getArray(&edge_array, &edge_size);

  // Check whether the edges have the same orientation
  int loc = edge_size-1;
  int dir = -1;
  if (n1 == d1 && n2 == d2){ 
    dir = 1; 
    loc = 0;
  }

  // Set the max quadrant side length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Loop over all of the nodes on the face
  for ( int i = 0; i < size; i++ ){
    int32_t u = -1;
    if (edge_num == 0 || edge_num == 1){
      if (array[i].x == edge_num*hmax){
        u = array[i].y;
      }
    }
    else {
      if (array[i].y == (edge_num % 2)*hmax){
        u = array[i].x;
      }
    }

    if (u >= 0){
      if (dir > 0 && loc < edge_size){
        int32_t v = edge_array[loc].x;
        if (dest_edge_num == 0 || dest_edge_num == 1){
          v = edge_array[loc].y;
        }

        for ( ; loc < edge_size && v < u; loc++ ){
          v = edge_array[loc].x;
          if (dest_edge_num == 0 || dest_edge_num == 1){
            v = edge_array[loc].y;
          }
        }

        if (u == v && loc < edge_size){
          edge_array[loc].tag = array[i].tag;
          loc++;
        }
      }
      else if (dir < 0 && loc >= 0){
        int32_t v = hmax - edge_array[loc].x;
        if (dest_edge_num == 0 || dest_edge_num == 1){
          v = hmax - edge_array[loc].y;
        }

        for ( ; loc >= 0 && v < u; loc-- ){
          v = hmax - edge_array[loc].x;
          if (dest_edge_num == 0 || dest_edge_num == 1){
            v = hmax - edge_array[loc].y;
          }
        }

        if (u == v && loc >= 0){
          edge_array[loc].tag = array[i].tag;
          loc--;
        }
      }
    }
  }
}

/*
  Set the node values into the quadtree for the given edge
*/
void TMRQuadForest::setEdgeNodeNums( int face, 
                                     TMRQuadrantArray *edge_nodes ){
  // Retrieve all the nodes from the faces
  TMRQuadrantArray *nodes;
  quadtrees[face]->getNodes(&nodes);

  // Get the size of the array
  int size;
  TMRQuadrant *array;
  nodes->getArray(&array, &size);

  // Get the nodes from the edge
  int edge_size;
  TMRQuadrant *edge_array;
  edge_nodes->getArray(&edge_array, &edge_size);

  for ( int i = 0, ii = 0; (i < size && ii < edge_size); i++ ){
    if (array[i].x == edge_array[ii].x  && 
        array[i].y == edge_array[ii].y){
      // Set the tag value into the face
      array[i].tag = edge_array[ii].tag;
      ii++;
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

  // Allocate the range of nodes
  node_range = new int[ mpi_size+1 ];
  memset(node_range, 0, (mpi_size+1)*sizeof(int));

  // Allocate all possible nodes on all of the trees
  for ( int face = 0; face < num_faces; face++ ){
    if (face_owners[face] == mpi_rank){
      quadtrees[face]->createNodes(mesh_order);
    }
  }

  // Find the face number corresponding to the owner for 
  // each edge and corner so that we know who should be 
  // ordering what!
  int *edge_face_owners = new int[ num_edges ];
  int *node_face_owners = new int[ num_nodes ];

  // The owner is chosen as the connecting face with the
  // lowest face number
  int nedge_sends = 0, nedge_recvs = 0;
  for ( int edge = 0; edge < num_edges; edge++ ){
    edge_face_owners[edge] = num_faces;

    int ipend = edge_face_ptr[edge+1];
    for ( int ip = edge_face_ptr[edge]; ip < ipend; ip++ ){
      if (edge_face_conn[ip] < edge_face_owners[edge]){
        edge_face_owners[edge] = edge_face_conn[ip];
      }
    }

    // Now count up the number of edges that will be sent and received
    int edge_owner = face_owners[edge_face_owners[edge]];
    
    // The edge owner will receive the edges that are not 
    // owned on this processor
    for ( int ip = edge_face_ptr[edge]; ip < ipend; ip++ ){
      int face = edge_face_conn[ip];

      // This processor is the owner and will receive the edge
      if (edge_owner == mpi_rank && 
          face_owners[face] != mpi_rank){
        nedge_recvs++;
      }

      // This processor owns the sending face
      if (face_owners[face] == mpi_rank && 
          edge_owner != mpi_rank){
        nedge_sends++;
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
  
  // Determine the dependent nodes for each face without labeling
  // the dependent nodes on the edges yet.
  for ( int face = 0; face < num_faces; face++ ){
    if (face_owners[face] == mpi_rank){
      // Get the quadrant elements
      TMRQuadrantArray *elements, *nodes;
      quadtrees[face]->getElements(&elements);
      quadtrees[face]->getNodes(&nodes);
      
      // Get the elements themselves
      int size;
      TMRQuadrant *array;
      elements->getArray(&array, &size);
      
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
        
          // Check if this query element lies within
          // the bounds of the current domain
          if ((p.x >= 0 && p.x < hmax) &&
              (p.y >= 0 && p.y < hmax)){
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

  // Allocate the array of quadrant arrays
  int edge_ptr_size = edge_face_ptr[num_edges];

  // Allocate the arrays of edge nodes
  TMRQuadrantArray **edge_nodes = new TMRQuadrantArray*[ edge_ptr_size ];
  memset(edge_nodes, 0, edge_ptr_size*sizeof(TMRQuadrantArray*));

  // Find the max of the send/recvs
  int max_send_recvs = nedge_sends;
  if (nedge_recvs > max_send_recvs){ 
    max_send_recvs = nedge_recvs; 
  }
  MPI_Request *requests = new MPI_Request[ max_send_recvs ];

  // Identify which nodes along each edge need to be sent to which
  // processors to determine whether there are dependent nodes at an
  // interface between quadtrees
  for ( int edge = 0, nsends = 0; edge < num_edges; edge++ ){
    // Get the owner (destination) rank for this edge
    int dest = face_owners[edge_face_owners[edge]];

    // Loop over all the faces connected to this edge
    int ipend = edge_face_ptr[edge+1];
    for ( int ip = edge_face_ptr[edge]; ip < ipend; ip++ ){
      // Get the corresponding face 
      int face = edge_face_conn[ip];

      // Check if this processor is the owner of the sending edge
      if (mpi_rank == face_owners[face]){
        // Find the local edge number
        int edge_num = 0;
        for ( ; edge_num < 4; edge_num++ ){
          if (face_edge_conn[4*face + edge_num] == edge){
            break;
          }
        }

        // Retrieve all the nodes from the faces
        TMRQuadrantArray *nodes;
        quadtrees[face]->getNodes(&nodes);
        
        // Retrieve the node array
        int size;
        TMRQuadrant *array;
        nodes->getArray(&array, &size);
        
        // Allocate the queue that will store all the quadrants
        // along this edge
        TMRQuadrantQueue *edge_queue = new TMRQuadrantQueue();

        // Scan through the quadrants and add the ones that
        // are on the edge to the queue
        const int32_t hmax = 1 << TMR_MAX_LEVEL;
        if (edge_num == 0 || edge_num == 1){
          for ( int i = 0; i < size; i++ ){
            if (array[i].x == edge_num*hmax){
              edge_queue->push(&array[i]);
            }
          }
        }
        else if (edge_num == 2 || edge_num == 3){
          for ( int i = 0; i < size; i++ ){
            if (array[i].y == (edge_num % 2)*hmax){
              edge_queue->push(&array[i]);
            }
          }
        }

        // Create the quadrant arrays from the queue
        edge_nodes[ip] = edge_queue->toArray();
        delete edge_queue;

        // Check whether this is a reciever or a sender
        // for this edge 
        if (dest != mpi_rank){
          // Get the data from the array of nodes and send it to the
          // receiving processor
          edge_nodes[ip]->getArray(&array, &size);

          // Send the array - note use an immediate send and tag the
          // message with the index into the receiver arrays
          MPI_Isend(array, size, TMRQuadrant_MPI_type,
                    dest, ip, comm, &requests[nsends]);
          nsends++;
        }
      }
    }
  }

  // Receive the arrays of incoming quadrants
  for ( int i = 0; i < nedge_recvs; i++ ){
    // Probe the recieved messages
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    // Retrieve the size and information for the incoming message
    int source = status.MPI_SOURCE;
    int ip = status.MPI_TAG;
    int size = 0;
    MPI_Get_count(&status, TMRQuadrant_MPI_type, &size);

    // Allocate the incoming array
    TMRQuadrant *array = new TMRQuadrant[ size ];
    MPI_Recv(array, size, TMRQuadrant_MPI_type,
             source, ip, comm, MPI_STATUS_IGNORE);

    // Allocate the received array
    edge_nodes[ip] = new TMRQuadrantArray(array, size);
  }

  // Wait for all the sends to complete
  MPI_Waitall(nedge_sends, requests, MPI_STATUSES_IGNORE);

  // Allocate an array to store the number of independent nodes
  // per edge in the forest
  int *num_indep_edge_nodes = new int[ num_edges ];
  memset(num_indep_edge_nodes, 0, num_edges*sizeof(int));

  // Loop through and label the dependent nodes
  for ( int edge = 0; edge < num_edges; edge++ ){
    int edge_owner = face_owners[edge_face_owners[edge]];
    if (mpi_rank == edge_owner){
      int ip = edge_face_ptr[edge];
      int ipend = edge_face_ptr[edge+1];
      int count = ipend - ip;

      // Label the dependent nodes
      num_indep_edge_nodes[edge] = 
        labelDependentNodes(edge, &edge_nodes[ip], count);

      // Set the labeled edges back into their owning
      // faces
      for ( ; ip < ipend; ip++ ){
        int face = edge_face_conn[ip];
        if (face_owners[face] == mpi_rank){
          setEdgeNodeNums(face, edge_nodes[ip]);
        }
      }
    }
  }

  // Count up all the local variables
  int nlocal = 0;

  // Count up the number of nodes owned on each face
  for ( int face = 0; face < num_faces; face++ ){
    if (mpi_rank == face_owners[face]){
      // Count up the number of locally owned corner nodes
      for ( int k = 0; k < 4; k++ ){
        int node = face_conn[4*face + k];
        if (node_face_owners[node] == face){
          nlocal++;
        }
      }

      // Count up the number of nodes owned on each face 
      // (not including the corner nodes at the end of each edge)
      for ( int k = 0; k < 4; k++ ){
        int edge = face_edge_conn[4*face + k];
        if (edge_face_owners[edge] == face){
          nlocal += num_indep_edge_nodes[edge]-2;
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
      }
    }
  }

  // Gather the local variable counts from each processor
  node_range[0] = 0;
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
        }
      }
    }
  }

  // Now, reverse the information stored in the quadrant
  // arrays sending them from the original recievers to the 
  // original senders. Copy these guys out and assign them
  // to the values in the quadtrees
  int nsends = 0;
  for ( int edge = 0; edge < num_edges; edge++ ){
    // Check if this processor is the owner of 'edge'
    int edge_owner = face_owners[edge_face_owners[edge]];

    // If this is the edge owner, then send all of the
    // ordered edges back to their original owners
    if (mpi_rank == edge_owner){
      // Loop over the faces that are connected to this edge
      for ( int ip = edge_face_ptr[edge]; 
            ip < edge_face_ptr[edge+1]; ip++ ){
        // Determine the face that this will go to and 
        // the MPI destination rank
        int face = edge_face_conn[ip];
        int dest = face_owners[face];

        // Get the edge node numbers for the specified
        // edge from the locally owned face
        getEdgeNodeNums(edge_face_owners[edge],
                        edge, face, edge_nodes[ip]);

        if (dest != mpi_rank){
          // Get the array
          int size;
          TMRQuadrant *array;
          edge_nodes[ip]->getArray(&array, &size);
          
          // Send the array back to the owner
          MPI_Isend(array, size, TMRQuadrant_MPI_type,
                    dest, ip, comm, &requests[nsends]);
          nsends++;
        }
      }
    }
  }

  // Receive the arrays of incoming quadrants
  for ( int i = 0; i < nedge_sends; i++ ){
    // Probe the recieved messages
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    // Retrieve the size and information for the incoming message
    int source = status.MPI_SOURCE;
    int ip = status.MPI_TAG;
    int size = 0;
    MPI_Get_count(&status, TMRQuadrant_MPI_type, &size);

    // Allocate the received array
    TMRQuadrant *array;
    edge_nodes[ip]->getArray(&array, NULL);

    // Recieve the quadrant arrays
    MPI_Recv(array, size, TMRQuadrant_MPI_type,
             source, ip, comm, MPI_STATUS_IGNORE);
  }

  // Set the received node numbers back into the
  // faces that are locally owned
  for ( int edge = 0; edge < num_edges; edge++ ){
    for ( int ip = edge_face_ptr[edge]; 
          ip < edge_face_ptr[edge+1]; ip++ ){
      int face = edge_face_conn[ip];
      int face_owner = edge_face_owners[ip];

      if (face_owners[face] == mpi_rank &&
          face != face_owner){
        setEdgeNodeNums(face, edge_nodes[ip]);
      }
    }
  }

  // Wait for all the sends to complete
  MPI_Waitall(nedge_recvs, requests, MPI_STATUSES_IGNORE);

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

  // Free the MPI requests
  delete [] requests;
  
  // Delete the allocated edges
  for ( int i = 0; i < edge_face_ptr[num_edges]; i++ ){
    if (edge_nodes[i]){ 
      delete edge_nodes[i]; 
    }
  }
  delete [] edge_nodes;

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

