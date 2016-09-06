#include "TMRQuadForest.h"

// Include METIS
extern "C" {
#include "metis.h"
}

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

/*
  Face to edge node connectivity
*/
const int face_to_edge_nodes[][2] = {{0, 2},
                                     {1, 3},
                                     {0, 1},
                                     {2, 3}};

/*
  Compare integers for sorting
*/
static int compare_integers( const void *a, const void *b ){
  return (*(int*)a - *(int*)b);
}

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

  // NULL the dependent node info
  num_owned_faces = 0;
  owned_faces = NULL;

  // Set the dependent info to NULL
  dep_edges = NULL;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;

  // Set the size of the mesh
  mesh_order = 0;
  
  // The face partition
  mpi_face_owners = NULL;
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
  if (dep_edges){
    for ( int i = 0; i < num_owned_faces; i++ ){
      if (dep_edges[i]){ delete dep_edges[i]; }
    }
    delete [] dep_edges;
  }
  if (mpi_face_owners){ delete [] mpi_face_owners; }
  if (owned_faces){ delete [] owned_faces; }
  if (node_range){ delete [] node_range; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
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
                                     int _num_faces,
                                     int partition ){
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
  if (dep_edges){
    for ( int i = 0; i < num_owned_faces; i++ ){
      if (dep_edges[i]){ delete dep_edges[i]; }
    }
    delete [] dep_edges;
  }
  if (mpi_face_owners){ delete [] mpi_face_owners; }
  if (owned_faces){ delete [] owned_faces; }
  if (node_range){ delete [] node_range; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
  node_range = NULL;
  dep_edges = NULL;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;

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

  // Set up the partition using metis
  mpi_face_owners = new int[ num_faces ];
  memset(mpi_face_owners, 0, num_faces*sizeof(int));

  // Compute the partition on the root processor and
  // broadcast the result to all other processors
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);  
  MPI_Comm_size(comm, &mpi_size);  
  if (mpi_rank == 0){
    if (partition){
      computePartition(mpi_size, NULL, mpi_face_owners);
    }
    else {
      for ( int face = 0; face < num_faces; face++ ){
        mpi_face_owners[face] = face % mpi_size;
      }
    }
  }
  
  // Broadcast the face owners to all processors
  MPI_Bcast(mpi_face_owners, num_faces, MPI_INT, 0, comm);

  // Determine the number of face owners
  num_owned_faces = 0;
  for ( int i = 0; i < num_faces; i++ ){
    if (mpi_rank == mpi_face_owners[i]){
      num_owned_faces++;
    }
  }

  owned_faces = new int[ num_owned_faces ];
  for ( int i = 0, k = 0; i < num_faces; i++ ){
    if (mpi_rank == mpi_face_owners[i]){
      owned_faces[k] = i;
      k++;
    }
  }

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
  for ( int face = 0; face < num_faces; face++ ){
    for ( int j = 0; j < 4; j++ ){
      int e = face_edge_conn[4*face + j];
      edge_face_conn[edge_face_ptr[e]] = 4*face + j;
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
    if (mpi_face_owners[i] == mpi_rank){
      quadtrees[i] = new TMRQuadtree(refine_level);
    }
  }
}

/*
  Create a forest with the specified refinement level
*/
void TMRQuadForest::createTrees( int refine_level[] ){
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
    if (mpi_face_owners[i] == mpi_rank){
      quadtrees[i] = new TMRQuadtree(refine_level[i]);
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
    if (mpi_face_owners[i] == mpi_rank){
      quadtrees[i] = new TMRQuadtree(nrand, min_level, max_level);
    }
  }
}

/*
  Get the array of quadtrees 
*/
int TMRQuadForest::getQuadtrees( TMRQuadtree ***_trees ){
  if (_trees){ *_trees = quadtrees; }
  return num_faces;
}

/*
  Get mesh/ownership information
*/
int TMRQuadForest::getOwnedQuadtrees( const int **_owned_faces ){
  if (_owned_faces){ 
    *_owned_faces = owned_faces; 
  }
  return num_owned_faces;
}

/*
  Retrieve information about the connectivity between faces, edges and
  nodes
*/
void TMRQuadForest::getConnectivity( int *_nfaces, 
                                     int *_nedges, int *_nnodes, 
                                     const int **_face_conn, 
                                     const int **_face_edge_conn ){
  if (_nfaces){ *_nfaces = num_faces; }
  if (_nedges){ *_nedges = num_edges; }
  if (_nnodes){ *_nnodes = num_nodes; }
  if (_face_conn){ *_face_conn = face_conn; }
  if (_face_edge_conn){ *_face_edge_conn = face_edge_conn; }
}

/*
  Retrieve the inverse of the connectivity
*/
void TMRQuadForest::getInverseConnectivity( const int **_node_face_conn,
                                            const int **_node_face_ptr,
                                            const int **_edge_face_conn,
                                            const int **_edge_face_ptr ){
  if (_node_face_conn){ *_node_face_conn = node_face_conn; }
  if (_node_face_ptr){ *_node_face_ptr = node_face_ptr; }
  if (_edge_face_conn){ *_edge_face_conn = edge_face_conn; }
  if (_edge_face_ptr){ *_edge_face_ptr = edge_face_ptr; }
}

/*
  Repartition the mesh based on the number of elements per face
*/
void TMRQuadForest::repartition(){
  // Get the communicator rank/size
  int mpi_size, mpi_rank;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // You can't repartition this on a single processor
  if (mpi_size <= 1){
    return;
  }

  // First, this stores the number of elements on quadtrees owned on
  // each processor
  int *elem_counts = new int[ num_owned_faces ];
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];
    elem_counts[owned] = quadtrees[face]->getNumElements();
  }

  // Allocate space for the new partition
  int *new_part = new int[ num_faces ];

  // Gather the element counts to the root processor
  if (mpi_rank != 0){
    MPI_Gatherv(elem_counts, num_owned_faces, MPI_INT,
                NULL, NULL, NULL, MPI_INT, 0, comm);
  }
  else {
    int *all_elem_counts = new int[ num_faces ];
    int *recv_counts = new int[ mpi_size ];
    int *recv_displ = new int[ mpi_size ];

    // Count up the recvs from each processor
    memset(recv_counts, 0, mpi_size*sizeof(int));
    for ( int face = 0; face < num_faces; face++ ){
      recv_counts[mpi_face_owners[face]]++;
    }
    
    // Compute the displacement offsets
    recv_displ[0] = 0;
    for ( int k = 1; k < mpi_size; k++ ){
      recv_displ[k] = recv_displ[k-1] + recv_counts[k-1];
    }
    
    // Receive all the elements
    MPI_Gatherv(elem_counts, num_owned_faces, MPI_INT,
                all_elem_counts, recv_counts, recv_displ, MPI_INT,
                0, comm);

    // Fill in the number of elements per processor back in the
    // original order from the original faces
    int *nelems = new int[ num_faces ];
    for ( int face = 0; face < num_faces; face++ ){
      int mpi_owner = mpi_face_owners[face];
      
      nelems[face] = all_elem_counts[recv_displ[mpi_owner]];
      recv_displ[mpi_owner]++;
    }

    // Free the offsets etc.
    delete [] all_elem_counts;
    delete [] recv_displ;
    delete [] recv_counts;

    // Compute the new partition on the root processor
    computePartition(mpi_size, nelems, new_part);

    // Free the number of elements per processor
    delete [] nelems;
  }

  // Free the element counts from each processor
  delete [] elem_counts;

  // Broadcast the new partition
  MPI_Bcast(new_part, num_faces, MPI_INT, 0, comm);

  // Now create the quadtrees array
  TMRQuadtree **new_quadtrees = new TMRQuadtree*[ num_faces ];
  memset(new_quadtrees, 0, num_faces*sizeof(TMRQuadtree*));

  // Only redistribute the elements, not the mesh
  int send_count = 0;
  MPI_Request *requests = new MPI_Request[ num_owned_faces ];
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];
    int dest = new_part[face];

    if (dest != mpi_rank){
      // Get the element array
      TMRQuadrantArray *elements;
      quadtrees[face]->getElements(&elements);
      
      // Get the actual quadrant array
      int size;
      TMRQuadrant *array;
      elements->getArray(&array, &size);

      // Send the element array to the new owner
      MPI_Isend(array, size, TMRQuadrant_MPI_type,
                dest, face, comm, &requests[send_count]);
      send_count++;
    }
    else {
      new_quadtrees[face] = quadtrees[face];
      quadtrees[face] = NULL;
    }
  }

  // Determine the new number of face owners based on the new
  // partition and the number of expected incoming messages (quadtrees)
  // from all processors
  int num_owned = 0;
  for ( int face = 0; face < num_faces; face++ ){
    if (mpi_rank == new_part[face]){
      num_owned++;
    }
  }

  int incoming_faces = 0;
  int *new_owned = new int[ num_owned ];
  for ( int face = 0, k = 0; face < num_faces; face++ ){
    if (mpi_rank == new_part[face]){
      new_owned[k] = face;
      k++;

      // Check if this face will need to be sent
      if (mpi_face_owners[face] != mpi_rank){
        incoming_faces++;
      }
    }
  }
  
  // Loop over the old senders
  for ( int k = 0; k < incoming_faces; k++ ){
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    // Get the source, tag == face, and size of the message
    int source = status.MPI_SOURCE;
    int face = status.MPI_TAG;
    int size = 0;
    MPI_Get_count(&status, TMRQuadrant_MPI_type, &size);

    // Allocate the incoming array and receive it
    TMRQuadrant *array = new TMRQuadrant[ size ];
    MPI_Recv(array, size, TMRQuadrant_MPI_type,
             source, face, comm, MPI_STATUS_IGNORE);
    
    // Create the new local quadtree
    TMRQuadrantArray *list = new TMRQuadrantArray(array, size);
    new_quadtrees[face] = new TMRQuadtree(list);
  }

  // Wait for all the sends to complete
  MPI_Waitall(send_count, requests, MPI_STATUSES_IGNORE);

  // Free the old quadtrees that were exchanged to another processor and
  // assign the new quadtrees array
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];
    if (quadtrees[face]){ 
      delete quadtrees[face]; 
    }
  }
  delete [] quadtrees;
  quadtrees = new_quadtrees;

  // Reset local information associated with the nodal ordering
  num_elements = 0;
  num_dep_nodes = 0;
  if (dep_edges){
    for ( int i = 0; i < num_owned_faces; i++ ){
      if (dep_edges[i]){ delete dep_edges[i]; }
    }
    delete [] dep_edges;
  }
  dep_edges = NULL;

  // Assign the new partition and the new owner list
  delete [] owned_faces;
  owned_faces = new_owned;
  num_owned_faces = num_owned;

  delete [] mpi_face_owners;
  mpi_face_owners = new_part;
}

/*
  Partition the mesh based on the super mesh connectivity and
  optionally vertex weights (typically the element count per face)

  input:
  num_part:     the number of partitions
  vwgts:        the vertex weights

  output:
  part:         the face-processor assignment
*/
void TMRQuadForest::computePartition( int part_size, int *vwgts,
                                      int *part ){
  // Set the pointer to the face connectivity
  int *ptr = new int[ num_faces+1 ];
  memset(ptr, 0, (num_faces+1)*sizeof(int));

  // Set the default size for the connectivity
  int max_size = 9*num_faces;
  int *conn = new int[ max_size ];
    
  // Allocate space to store flags to indicate whether the face has
  // been added to the array
  int *cols = new int[ num_faces ];
  memset(cols, 0, num_faces*sizeof(int));

  // The adjacency node count
  int min_adj_count = 4;
  
  for ( int face = 0; face < num_faces; face++ ){
    // Set the new pointer
    ptr[face+1] = ptr[face];
    
    // Loop over all the faces
    for ( int j = 0; j < 4; j++ ){
      int node = face_conn[4*face + j];
      
      // Loop over all of the nodes
      for ( int ip = node_face_ptr[node]; ip < node_face_ptr[node+1]; ip++ ){
        int adj_face = node_face_conn[ip];

        if (cols[adj_face]+1 == min_adj_count*(face+1)){
          // Extend the array
          if (ptr[face+1] >= max_size){
            max_size *= 2;
            int *tmp = new int[ max_size ];
            memcpy(tmp, conn, ptr[face+1]*sizeof(int));
            delete [] conn;
            conn = tmp;
          }

          // Set the new element into the connectivity
          conn[ptr[face+1]] = adj_face;
          ptr[face+1]++;

          // Set the flag so that its not added again
          cols[adj_face]++;
        }
        else if (cols[adj_face] <= min_adj_count*face){
          // This is the first time this face has been encountered on
          // this loop. Add it and increment the pointer
          cols[adj_face] = min_adj_count*face+1;
        }
        else {
          // This face has been encountered before
          cols[adj_face]++;
        }
      }
    }
    
    // Sort the array
    int len = ptr[face+1] - ptr[face];
    qsort(&conn[ptr[face]], len, sizeof(int), compare_integers);
  }
  
  // use the default options
  int options[5];
  options[0] = 0; 
  
  // weights are on the verticies
  int wgtflag = 0; 
  int numflag = 0;
  int edgecut = -1; 
  
  // Weights on vertices and edges
  int *adjwgts = NULL;
  
  if (part_size < 8){
    METIS_PartGraphRecursive(&num_faces, ptr, conn, vwgts, adjwgts,
                             &wgtflag, &numflag, &part_size,
                             options, &edgecut, part);
  }
  else {
    METIS_PartGraphKway(&num_faces, ptr, conn, vwgts, adjwgts, 
                        &wgtflag, &numflag, &part_size, 
                        options, &edgecut, part);
  }
  
  // Free the allocated data
  delete [] cols;
  delete [] conn;
  delete [] ptr;
}

/*
  Duplicate the forest

  This function creates a duplicate representation of the current
  forest. This function copies the global connectivity of the forest
  and copies each individual tree.
*/
TMRQuadForest *TMRQuadForest::duplicate(){
  TMRQuadForest *dup = new TMRQuadForest(comm);
  if (quadtrees){
    // Copy over the connectivity data 
    dup->num_nodes = num_nodes;
    dup->num_edges = num_edges;
    dup->num_faces = num_faces;
    dup->num_owned_faces = num_owned_faces;

    // Allocate/copy the face connectivities
    dup->face_conn = new int[ 4*num_faces ];
    dup->face_edge_conn = new int[ 4*num_faces ];
    memcpy(dup->face_conn, face_conn, 4*num_faces*sizeof(int));
    memcpy(dup->face_edge_conn, face_edge_conn, 4*num_faces*sizeof(int));
    
    // Allocate/copy the inverse relationships
    dup->node_face_ptr = new int[ num_nodes+1 ];
    dup->node_face_conn = new int[ node_face_ptr[num_nodes] ];
    memcpy(dup->node_face_ptr, node_face_ptr, 
           (num_nodes+1)*sizeof(int));
    memcpy(dup->node_face_conn, node_face_conn, 
           node_face_ptr[num_nodes]*sizeof(int));

    dup->edge_face_ptr = new int[ num_edges+1 ];
    dup->edge_face_conn = new int[ edge_face_ptr[num_edges] ];
    memcpy(dup->edge_face_ptr, edge_face_ptr, 
           (num_edges+1)*sizeof(int));
    memcpy(dup->edge_face_conn, edge_face_conn, 
           edge_face_ptr[num_edges]*sizeof(int));

    // Allocate/copy the face ownership
    dup->mpi_face_owners = new int[ num_faces ];
    memcpy(dup->mpi_face_owners, mpi_face_owners, 
           num_faces*sizeof(int));

    dup->owned_faces = new int[ num_owned_faces ];
    memcpy(dup->owned_faces, owned_faces, 
           num_owned_faces*sizeof(int));

    // Duplicate all the quadtrees
    dup->quadtrees = new TMRQuadtree*[ num_faces ];
    memset(dup->quadtrees, 0, num_faces*sizeof(TMRQuadtree*));
    for ( int i = 0; i < num_owned_faces; i++ ){
      int face = owned_faces[i];
      TMRQuadrantArray *elements;
      quadtrees[face]->getElements(&elements);
      dup->quadtrees[face] = new TMRQuadtree(elements->duplicate());
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
TMRQuadForest *TMRQuadForest::coarsen(){
  TMRQuadForest *coarse = new TMRQuadForest(comm);
  if (quadtrees){
    // Copy over the connectivity data 
    coarse->num_nodes = num_nodes;
    coarse->num_edges = num_edges;
    coarse->num_faces = num_faces;
    coarse->num_owned_faces = num_owned_faces;

    // Allocate/copy the face connectivities
    coarse->face_conn = new int[ 4*num_faces ];
    coarse->face_edge_conn = new int[ 4*num_faces ];
    memcpy(coarse->face_conn, face_conn, 4*num_faces*sizeof(int));
    memcpy(coarse->face_edge_conn, face_edge_conn, 4*num_faces*sizeof(int));
    
    // Allocate/copy the inverse relationships
    coarse->node_face_ptr = new int[ num_nodes+1 ];
    coarse->node_face_conn = new int[ node_face_ptr[num_nodes] ];
    memcpy(coarse->node_face_ptr, node_face_ptr, 
           (num_nodes+1)*sizeof(int));
    memcpy(coarse->node_face_conn, node_face_conn, 
           node_face_ptr[num_nodes]*sizeof(int));

    coarse->edge_face_ptr = new int[ num_edges+1 ];
    coarse->edge_face_conn = new int[ edge_face_ptr[num_edges] ];
    memcpy(coarse->edge_face_ptr, edge_face_ptr, 
           (num_edges+1)*sizeof(int));
    memcpy(coarse->edge_face_conn, edge_face_conn, 
           edge_face_ptr[num_edges]*sizeof(int));

    // Allocate/copy the face ownership
    coarse->mpi_face_owners = new int[ num_faces ];
    memcpy(coarse->mpi_face_owners, mpi_face_owners, 
           num_faces*sizeof(int));

    coarse->owned_faces = new int[ num_owned_faces ];
    memcpy(coarse->owned_faces, owned_faces,
           num_owned_faces*sizeof(int));

    // Coarsen all the quadtrees
    coarse->quadtrees = new TMRQuadtree*[ num_faces ];
    memset(coarse->quadtrees, 0, num_faces*sizeof(TMRQuadtree*));
    for ( int i = 0; i < num_owned_faces; i++ ){
      int face = owned_faces[i];
      coarse->quadtrees[face] = quadtrees[face]->coarsen();
    }
  }

  return coarse;
}

/*
  Add the edge neighbors for a adjacent trees

  This function is called to balance the forest across tree edges.
  Given an quadant p on the specified corner index, this code ensures
  a corner balanced tree, by adding the corresponding corner qudrant
  to all node-adjacent faces. If the quadant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  face:        the local face index (p is defined on this face)
  edge_index:  the edge index
  p:           the local quadant
  hash:        the array of hash objects for each face
  queue:       the array of newly added qudrants for each face
*/
void TMRQuadForest::addEdgeNeighbors( int face,
                                      int edge_index, 
                                      TMRQuadrant p,
                                      TMRQuadrantHash **hash,
                                      TMRQuadrantQueue **queue ){
  // First determine the global edge number 
  int edge = face_edge_conn[4*face + edge_index];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Store the u coordinate along the edge
  int32_t ucoord = 0;
  if (edge_index < 2){
    ucoord = p.y;
  }
  else {
    ucoord = p.x;
  }

  // Retrieve the first and second node numbers
  int n1 = face_conn[4*face + face_to_edge_nodes[edge_index][0]];
  int n2 = face_conn[4*face + face_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent faces
  for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
    // Get the face that is adjacent across this edge
    int adjacent = edge_face_conn[ip]/4;
    if (adjacent != face){
      if (!hash[adjacent]){
        hash[adjacent] = new TMRQuadrantHash();
        queue[adjacent] = new TMRQuadrantQueue();
      }

      // Get the adjacent edge index
      int adj_index = edge_face_conn[ip] % 4;

      // Get the nodes on the adjacent face
      int nn1 = face_conn[4*adjacent + face_to_edge_nodes[adj_index][0]];
      int nn2 = face_conn[4*adjacent + face_to_edge_nodes[adj_index][1]];

      // Add the quadant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - 2*h - ucoord;
      }
      
      TMRQuadrant neighbor;
      neighbor.level = p.level;
      if (adj_index < 2){
        neighbor.x = (hmax - 2*h)*(adj_index % 2);
        neighbor.y = u;
      }
      else {
        neighbor.x = u;
        neighbor.y = (hmax - 2*h)*(adj_index % 2);
      }
      
      // Add the quadant to the list
      if (hash[adjacent]->addQuadrant(&neighbor)){
        queue[adjacent]->push(&neighbor);
      }
    }
  }
}

/*
  Add the corner neighbors for a given tree

  This function is called to balance the forest across tree corners.
  Given a quadrant p on the specified corner index, this code ensures
  a corner balanced tree, by adding the corresponding corner qudrant
  to all node-adjacent faces. If the quadrant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  face:    the local face index (p is defined on this face)
  corner:  the corner index (p must lie on this corner)
  p:       the local quadant
  hash:    the array of hash objects for each face
  queue:   the array of newly added qudrants for each face
*/
void TMRQuadForest::addCornerNeighbors( int face,
                                        int corner,
                                        TMRQuadrant p,
                                       TMRQuadrantHash **hash,
                                       TMRQuadrantQueue **queue ){
  // First determine the global edge number 
  int node = face_conn[4*face + corner];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = node_face_ptr[node];
        ip < node_face_ptr[node+1]; ip++ ){
      
    // Get the faces that are adjacent across this edge
    int adjacent = node_face_conn[ip];
    if (adjacent != face){

      // Loop over all of the corners in this adjacent face
      for ( int j = 0; j < 4; j++ ){
        if (face_conn[4*adjacent+j] == node){
          // Allocate the quadant hash/queue if not already allocated
          if (!hash[adjacent]){
            hash[adjacent] = new TMRQuadrantHash();
            queue[adjacent] = new TMRQuadrantQueue();
          }

          // Compute the quadant location
          TMRQuadrant neighbor;
          neighbor.level = p.level;
          neighbor.x = (hmax - 2*h)*(j % 2);
          neighbor.y = (hmax - 2*h)*(j/2);

          // Add the quadant to the list
          if (hash[adjacent]->addQuadrant(&neighbor)){
            queue[adjacent]->push(&neighbor);
          }
        }
      }
    }
  }
}

/*
  Balance the quadant on the entire quadtree

  This code finds the 0-parent of all adjacent quadants either within
  the current tree or within an adjacent tree and adds those quadants
  to balance the input quadant 'quad'.

  input:
  face:            index of the face for the quadant
  quad:            the quadant itself
  hash:            the array of hash tables for each face
  queue:           the quadant queues for each face
  balance_corner:  balance across corners 
  balance_tree:    balance on the entire tree
*/
void TMRQuadForest::balanceQuadrant( int face,
                                     TMRQuadrant *quad,
                                     TMRQuadrantHash **hash,
                                     TMRQuadrantQueue **queue,
                                     const int balance_corner,
                                     const int balance_tree ){
  // Local quadant data
  TMRQuadrant p, neighbor, q;

  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  
  // Get the parent of the quadant, and add the their
  // face-matched quadants from each face, as long 
  // as they fall within the bounds
  if (quad->level > 1){
    quad->parent(&p);

    // Add the edge-adjacent elements
    for ( int edge_index = 0; edge_index < 4; edge_index++ ){
      p.edgeNeighbor(edge_index, &neighbor);
      neighbor.getSibling(0, &q);
	  
      // If we're in bounds, add the neighbor
      if ((q.x >= 0 && q.x < hmax) &&
          (q.y >= 0 && q.y < hmax)){
        if (hash[face]->addQuadrant(&q)){
          queue[face]->push(&q);
        }
      }
      else if (balance_tree){
        // The node may lie across an edge or face
        int ex0 = (q.x < 0);
        int ey0 = (q.y < 0); 
        int ex = (ex0 || q.x >= hmax);
        int ey = (ey0 || q.y >= hmax);
        
        if (ex && ey){
          int corner = (ex0 ? 1 : 0) + (ey0 ? 2 : 0);
          addCornerNeighbors(face, corner, neighbor, hash, queue);
        }
        else {
          // The quadant lies along an edge
          addEdgeNeighbors(face, edge_index, q, hash, queue);
        }
      }
    }

    // If we're balancing across edges and 
    if (balance_corner){
      for ( int corner = 0; corner < 4; corner++ ){
        p.cornerNeighbor(corner, &neighbor);
        neighbor.getSibling(0, &q);

        // If we're in bounds, add the neighbor
        if ((q.x >= 0 && q.x < hmax) &&
            (q.y >= 0 && q.y < hmax)){
          if (hash[face]->addQuadrant(&q)){
            queue[face]->push(&q);
          }
        }
        else if (balance_tree){
          // The node may lie across an edge or face
          int ex0 = (q.x < 0);
          int ey0 = (q.y < 0); 
          int ex = (ex0 || q.x >= hmax);
          int ey = (ey0 || q.y >= hmax);
          
          if (ex && ey){
            int corner = (ex0 ? 1 : 0) + (ey0 ? 2 : 0);
            addCornerNeighbors(face, corner, neighbor, hash, queue);
          }
          else {
            // The quadant lies along an edge
            int edge_index = ex*(ex0 ? 0 : 1) + ey*(ey0 ? 2 : 3);
            addEdgeNeighbors(face, edge_index, q, hash, queue);
          }
        }
      }
    }
  }
}

/*
  Balance the forest of quadtrees

  This algorithm uses a hash and a queue to balance the forest of
  quadtrees. For each element in the quadtree, we add the neighbors that
  are required to balance to the tree. If the element is not in the
  hash, we add them to a queue, which keeps track of recently added
  elements. After the first pass, the algorithm continues popping
  elements until the all the queues are empty.

  Note that only 0-th siblings are added/popped on the hash/queue.
  Then at the end, all neighboring siblings are added.

  The type of balancing - face/edge balanced or face/edge/corner
  balanced is determined using the balance_corner flag. Face balancing
  is balancing across faces, edge balancing is balancing across edges
  of the elements and corner balances across corners. The code always
  balances faces and edges (so that there is at most one depdent node
  per edge) and balances across corners optionally.
*/
void TMRQuadForest::balance( int balance_corner ){
  // Get the mpi rank
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Create a hash table for the balanced tree
  TMRQuadrantHash **hash = new TMRQuadrantHash*[ num_faces ];
  TMRQuadrantQueue **queue = new TMRQuadrantQueue*[ num_faces ];
  
  // Set the hash tables and queues to NULL
  memset(hash, 0, num_faces*sizeof(TMRQuadrantHash*));
  memset(queue, 0, num_faces*sizeof(TMRQuadrantQueue*));

  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Allocate the hash and queue if they are not already
    // allocated on this processor
    if (!hash[face]){
      hash[face] = new TMRQuadrantHash();
      queue[face] = new TMRQuadrantQueue();
    }

    // Get the current array of quadrants
    TMRQuadrantArray *elements;
    quadtrees[face]->getElements(&elements);
    
    // Get the array of elements
    int size;
    TMRQuadrant *array;
    elements->getArray(&array, &size);
    
    // Add all the elements
    for ( int i = 0; i < size; i++ ){
      TMRQuadrant quad;
      array[i].getSibling(0, &quad);
      hash[face]->addQuadrant(&quad);
      
      // Balance the quadrant by push the neighbors
      const int balance_tree = 1;
      balanceQuadrant(face, &quad, hash, queue, 
                      balance_corner, balance_tree);
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
          TMRQuadrant quad = queue[face]->pop();
          // Balance the quadrant by push the neighbors
          const int balance_tree = 1;
          balanceQuadrant(face, &quad, hash, queue, 
                          balance_corner, balance_tree);
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
    if (mpi_face_owners[face] != mpi_rank && queue[face]){ 
      delete queue[face];
      queue[face] = NULL;
    }
  }

  // Now everything is locally balanced - all the elements 
  // on the current processor are balanced with all the other
  // elements on the current processor, but nothing is 
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
    if (hash[face] && 
        (mpi_face_owners[face] != mpi_rank)){
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
      int dest_rank = mpi_face_owners[face];
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
        const int balance_tree = 0;
        TMRQuadrant quad = queue[face]->pop();
        balanceQuadrant(face, &quad, hash, queue, 
                        balance_corner, balance_tree);
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
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

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
  Add the quadrant to the processor queues that correspond to 
  the non-local faces that touch the corner
*/
void TMRQuadForest::addCornerQuadrantToQueues( const int node,
                                               const int mpi_rank,
                                               TMRQuadrant *q,
                                               TMRQuadrantQueue **queues ){
  for ( int ip = node_face_ptr[node]; ip < node_face_ptr[node+1]; ip++ ){
    int face = node_face_conn[ip];
    int rank = mpi_face_owners[face];
    if (rank != mpi_rank){
      queues[rank]->push(q);
    }
  }
}

/*
  Add the quadrant to the processor queues corresponding to the
  non-local faces that touch the given edge
*/
void TMRQuadForest::addEdgeQuadrantToQueues( const int edge,
                                             const int mpi_rank,
                                             TMRQuadrant *q,
                                             TMRQuadrantQueue **queues ){
  for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
    int face = edge_face_conn[ip]/4;
    int rank = mpi_face_owners[face];
    if (rank != mpi_rank){
      queues[rank]->push(q);
    }
  }
}

/*
  The following code exchanges the neighboring quadrants for each
  locally owned quadtree within the forest.

  This code exchanges non-local quadrants across processors so that we
  can locally query quadrants on adjacent quadtrees without having to
  perform parallel communication.

  Note that this code creates partial non-local quadtrees that are
  adjacent to the local quadtrees in the forest. These partial local
  quadtrees should be freed after the nodal ordering has been computed.
*/
void TMRQuadForest::recvQuadNeighbors(){
  // Get the MPI information
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Allocate the queues that store the quadrants destined for each of
  // the processors
  TMRQuadrantQueue **queues = new TMRQuadrantQueue*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      queues[k] = new TMRQuadrantQueue();
    }
    else {
      queues[k] = NULL;
    }
  }

  // For each face, determine the edge, face and corner neighbors 
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];
    
    // Flag to indicate whether any edge/face is non-local
    int has_non_local = 0;
    
    // Check if any of the edges has to be sent to another processor
    for ( int k = 0; k < 4; k++ ){
      int edge = face_edge_conn[4*face + k];
      for ( int ip = edge_face_ptr[edge];
            ip < edge_face_ptr[edge+1]; ip++ ){
        int dest_face = edge_face_conn[ip]/4;
        if (dest_face != mpi_rank){
          has_non_local = 1;
          break;
        }
      }
    }
    
    if (has_non_local){
      // Get the element array
      TMRQuadrantArray *elements;
      quadtrees[face]->getElements(&elements);
      
      // Get the actual quadrant array
      int size;
      TMRQuadrant *array;
      elements->getArray(&array, &size);
        
      // If this is a single quadrant, it fills the entire volume
      if (size == 1){
        TMRQuadrant q = array[0];
        q.tag = face;

        // Send the quadrant to all neighbors
        for ( int node_index = 0; node_index < 4; node_index++ ){
          int node = face_conn[4*face + node_index];
          addCornerQuadrantToQueues(node, mpi_rank, &q, queues);
        }

        // Send the quadrant to all adjacent edges
        for ( int edge_index = 0; edge_index < 4; edge_index++ ){
          int edge = face_edge_conn[4*face + edge_index];
          addEdgeQuadrantToQueues(edge, mpi_rank, &q, queues);
        }
      }
      else {
        // Loop over all the elements and check where we need to send
        // the quadrants that are along each edge/face
        for ( int i = 0; i < size; i++ ){
          const int32_t hmax = 1 << TMR_MAX_LEVEL; 
          const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
          
          // Determine which faces the quadrant touches if any
          int fx0 = (array[i].x == 0);
          int fy0 = (array[i].y == 0);
          int fx = (fx0 || array[i].x + h == hmax);
          int fy = (fy0 || array[i].y + h == hmax);
          
          // Check whether the quadrant lies along an edge
          if (fx || fy){
            // Copy the quadrant from the array and set its tag 
            // to the source face index
            TMRQuadrant q = array[i];
            q.tag = face;
            
            // Pass the face to any adjacent quadtrees
            if (fx && fy){
              int node_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2);
              int node = face_conn[4*face + node_index];
              addCornerQuadrantToQueues(node, mpi_rank, &q, queues);
            }
            else if (fx){
              int edge_index = (fx0 ? 0 : 1);
              int edge = face_edge_conn[4*face + edge_index];
              addEdgeQuadrantToQueues(edge, mpi_rank, &q, queues);
            }
            else if (fy){
              int edge_index = (fy0 ? 2 : 3);
              int edge = face_edge_conn[4*face + edge_index];
              addEdgeQuadrantToQueues(edge, mpi_rank, &q, queues);
            }
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

  // Wait for all the sends to complete
  MPI_Waitall(nsends, send_requests, MPI_STATUSES_IGNORE);
  delete [] send_requests;
  
  // Now free the arrays
  for ( int k = 0; k < mpi_size; k++ ){
    if (arrays[k]){
      delete arrays[k];
    }
  }
  delete [] arrays;
}

/*
  Determine if there is an adjacent quadrant on the connecting edge.

  Return true if an adjacent edge is found across a face-edge and
  false if no quadrant is found.
  
  input:
  edge:          the edge number
  edge_index:    the local edge index
  face_owner:   the index of the face owner
  b:             the quadrant
*/
int TMRQuadForest::checkAdjacentDepEdges( int edge,
                                          int edge_index,
                                          int face_owner,
                                          TMRQuadrant *b ){
  // Get the side length of the quadrant
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);

  // Store the u coordinate along the edge
  int32_t ucoord = 0;
  if (edge_index < 2){
    ucoord = b->y;
  }
  else {
    ucoord = b->x;
  }

  // Retrieve the first and second node numbers
  int n1 = face_conn[4*face_owner + face_to_edge_nodes[edge_index][0]];
  int n2 = face_conn[4*face_owner + face_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent edges
  for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
    int face = edge_face_conn[ip]/4;

    if (face_owner != face){
      // Get the adjacent edge index
      int adj_index = edge_face_conn[ip] % 4;

      // Get the nodes on the adjacent face
      int nn1 = face_conn[4*face + face_to_edge_nodes[adj_index][0]];
      int nn2 = face_conn[4*face + face_to_edge_nodes[adj_index][1]];

      // Add the quadrant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - h - ucoord;
      }

      // Search for the neighboring quadrant
      TMRQuadrant quad;
      quad.level = b->level;
      if (adj_index < 2){
        quad.x = (hmax - h)*(adj_index % 2);
        quad.y = u;
      }
      else {
        quad.x = u;
        quad.y = (hmax - h)*(adj_index % 2);
      }
      
      // Get the elements from the array
      TMRQuadrantArray *elements;
      quadtrees[face]->getElements(&elements);
      
      // If the adjacent element exists then label the
      // corresponding nodes as dependent
      const int use_nodes = 0;
      if (elements->contains(&quad, use_nodes)){
        return 1;
      }
    }
  }

  return 0;
}

/*
  Compute the dependent nodes (hanging edge/face nodes) on each face
  and on the interfaces between adjacent faces.

  The hanging face nodes may occur on any face within the face and
  are associated with the 4 parallel hanging edge nodes.  Within this
  code, we only store the corresponding face node and use the
  associated edges. Edge nodes along face interfaces may be hanging
  even if there is no associated hanging face node.

  side effects:
  dep_edges:   a list of dependent edges (aligned with face edges)
*/
void TMRQuadForest::computeDepEdges(){
  if (dep_edges){ 
    for ( int i = 0; i < num_owned_faces; i++ ){
      if (dep_edges[i]){ delete dep_edges[i]; }
    }
    delete [] dep_edges;
    dep_edges = NULL;
  }
  if (dep_ptr){ delete [] dep_ptr;  dep_ptr = NULL; }
  if (dep_conn){ delete [] dep_conn;  dep_conn = NULL; }
  if (dep_weights){ delete [] dep_weights;  dep_weights = NULL; }

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Allocate the list of dependent nodes/faces
  dep_edges = new TMRQuadrantArray*[ num_owned_faces ];
  memset(dep_edges, 0, num_owned_faces*sizeof(TMRQuadrantArray*));

  // Determine and label the dependent nodes on each processor 
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Create the queue that will store the dependent
    // faces/edges
    TMRQuadrantQueue *dedges = new TMRQuadrantQueue();
    
    // Get the quadrant elements
    TMRQuadrantArray *elements;
    quadtrees[face]->getElements(&elements);
            
    // Get the elements themselves
    int size;
    TMRQuadrant *array;
    elements->getArray(&array, &size);
    
    for ( int i = 0; i < size; i++ ){
      // Get the child id - skip any elements that are not
      // either child-0 or child-7
      int child_id = array[i].childId();

      // Get the side length of the element
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      const int32_t hmax = 1 << TMR_MAX_LEVEL;

      if (array[i].level > 0 && 
          (child_id == 0 || child_id == 3)){        
        // Check the adjacent elements along each face
        int edge_index = 0;
        TMRQuadrant p = array[i];
        if (child_id == 3){
          edge_index = 1;
          p.getSibling(0, &p);
        }
        p.level = p.level - 1;

        for ( ; edge_index < 4; edge_index += 2 ){
          // Look for the face neighbor that is at the next level of
          // refinement from the current edge
          TMRQuadrant q;
          p.edgeNeighbor(edge_index, &q);

          // Check if this element is across a face
          int fx0 = (q.x < 0); 
          int fy0 = (q.y < 0);
          int fx = (fx0 || q.x >= hmax);
          int fy = (fy0 || q.y >= hmax);

          // If an adjacent quadrant of twice the size is on a
          // neighboring tree, or on the same tree, add it as
          // a dependent face
          int add_me = 0;
          if (fx || fy){
            int edge = face_edge_conn[4*face + edge_index];
            if (checkAdjacentDepEdges(edge, edge_index, face, &q)){
              add_me = 1;
            }
          }
          else {
            // If the more-refined element exists then label the
            // corresponding nodes as dependent
            const int use_nodes = 0;
            if (elements->contains(&q, use_nodes)){
              add_me = 1;
            }
          }

          // Add the dependent face
          if (add_me){
            p.tag = edge_index;
            dedges->push(&p);
          }
        }
      }
    }

    // Create the arrays of the dependent nodes/faces
    dep_edges[owned] = dedges->toArray();

    // Delete the queues
    delete dedges;
  }
}

/*
  Label the dependent face and edge nodes

  This code is called after all the dependent faces have been
  computed.  Note that this relies on the mesh being edge-balanced
  (which is required). 
*/
void TMRQuadForest::labelDependentNodes(){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Get the nodes associated with this quadtree
    TMRQuadrantArray *nodes;
    quadtrees[face]->getNodes(&nodes);

    // Get the array of dependent faces
    int dep_size;
    TMRQuadrant *dep_array;
    dep_edges[owned]->getArray(&dep_array, &dep_size);

    for ( int i = 0; i < dep_size; i++ ){
      TMRQuadrant *b = &dep_array[i];

      // Find the edge lengths of the element quadrant and
      // the node spacing on the dependent face
      const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);
      const int32_t hc = 1 << (TMR_MAX_LEVEL - b->level - (mesh_order-1));
          
      // Get the edge index
      int edge_index = b->tag;

      // Loop over all the nodes on this edge
      for ( int ii = 0; ii < 2*mesh_order-1; ii++ ){
        // Label only the dependent nodes
        if (ii % 2 == 1){
          TMRQuadrant node;
          if (edge_index < 2){
            node.x = b->x + h*(edge_index % 2);
            node.y = b->y + ii*hc;
          }
          else {
            node.x = b->x + ii*hc;
            node.y = b->y + h*(edge_index % 2);
          }
            
          // Search for dependent node and label it
          const int use_node_search = 1;
          TMRQuadrant *t = nodes->contains(&node, use_node_search);
          t->tag = -1;
        }
      }
    }
  }
}

/*
  Return an array of flags to indicate whether the given face owns
  each of the edges and nodes 
*/
void TMRQuadForest::getOwnerFlags( int face,
                                   const int *edge_face_owners,
                                   const int *node_face_owners,
                                   int *is_edge_owner, 
                                   int *is_node_owner ){
  // Check whether the face owns the node, edge or node
  if (is_node_owner){
    for ( int k = 0; k < 4; k++ ){
      int node = face_conn[4*face + k];
      is_node_owner[k] = (node_face_owners[node] == face);
    }
  }
  
  // Get the edge owners
  if (is_edge_owner){
    for ( int k = 0; k < 4; k++ ){
      int edge = face_edge_conn[4*face + k];
      is_edge_owner[k] = (edge_face_owners[edge] == face);
    }
  }
}

/*
  Find a global ordering for the nodes now that the dependent nodes
  have been labeled with a negative index.

  This code first counts up the number of nodes that are owned locally
  and then performs an all gather to determine the offsets to the node
  numbering on each local processor.

  input:
  edge_owners:  the face index of the edge owner
  node_owners:  the face index of the corner/node owner
*/
void TMRQuadForest::orderGlobalNodes( const int *edge_face_owners,
                                      const int *node_face_owners ){
  // Get the MPI rank/size 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Count up all the locally owned, global variables and the number
  // of dependent variables
  int nlocal = 0;  
  
  // Set the number of elements/dependent nodes
  num_dep_nodes = 0;
  num_elements = 0;

  // Count up the number of nodes owned on each face
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Check whether the face owns the nodes, edges or faces 
    int is_edge_owner[4], is_node_owner[4];
    getOwnerFlags(face, edge_face_owners, node_face_owners,
                  is_edge_owner, is_node_owner);

    // Get the quadrant array of nodes
    TMRQuadrantArray *nodes;
    quadtrees[face]->getNodes(&nodes);

    // Extract the array of quadrants
    int size;
    TMRQuadrant *array;
    nodes->getArray(&array, &size);

    // Count up the number of elements
    num_elements += quadtrees[face]->getNumElements();

    // Count up the number of interior nodes
    const int32_t hmax = 1 << TMR_MAX_LEVEL;
    for ( int i = 0; i < size; i++ ){
      if (array[i].tag >= 0){
        // Determine which faces the quadrant touches if any
        int fx0 = (array[i].x == 0);
        int fy0 = (array[i].y == 0);
        int fx = (fx0 || array[i].x == hmax);
        int fy = (fy0 || array[i].y == hmax);
          
        // Check whether this node is on the face, edge, corner or is
        // an internal node and order it only if it is locally owned
        if (fx && fy){
          int corner_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2);
          if (is_node_owner[corner_index]){ nlocal++; }
        }
        else if (fx){
          int edge_index = (fx0 ? 0 : 1);
          if (is_edge_owner[edge_index]){ nlocal++; }
        }
        else if (fy){
          int edge_index = (fy0 ? 2 : 3);
          if (is_edge_owner[edge_index]){ nlocal++; }
        }
        else {
          nlocal++;
        }
      }
      else {
        num_dep_nodes++;
        array[i].tag = -num_dep_nodes;
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
  // using the offsets computed above.  Set the node number offset
  // from the range of owned nodes
  int node_num = node_range[mpi_rank];
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Check whether the face owns the nodes, edges or faces 
    int is_edge_owner[4], is_node_owner[4];
    getOwnerFlags(face, edge_face_owners, node_face_owners,
                  is_edge_owner, is_node_owner);

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
      if (array[i].tag >= 0){
        // Determine which faces the quadrant touches if any
        int fx0 = (array[i].x == 0);
        int fy0 = (array[i].y == 0);
        int fx = (fx0 || array[i].x == hmax);
        int fy = (fy0 || array[i].y == hmax);
          
        // Check whether this node is on the face, edge, corner or is
        // an internal node and order it only if it is locally owned
        if (fx && fy){
          int corner_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2);
          if (is_node_owner[corner_index]){ 
            array[i].tag = node_num;  node_num++; 
          }
        }
        else if (fx){
          int edge_index = (fx0 ? 0 : 1);
          if (is_edge_owner[edge_index]){ 
            array[i].tag = node_num;  node_num++; 
          }
        }
        else if (fy){
          int edge_index = (fy0 ? 2 : 3);
          if (is_edge_owner[edge_index]){ 
            array[i].tag = node_num;  node_num++; 
          }
        }
        else {
          array[i].tag = node_num;  node_num++; 
        }
      }
    }
  }
}

/*
  Copy the corner node variable number on the corner with the given corner
  index to the adjacent corners.  These values are then communicated
  back to their owners in a second parallel step.

  input:
  node:        the corner number that the node lies on
  node_index:  the local index of the corner 0 <= node index < 4
  face_owner:  the face index that owns the node
  p:           the quadrant with the correct tag value
*/
void TMRQuadForest::copyCornerNodes( int node,
                                     int node_index,
                                     int face_owner,
                                     TMRQuadrant *p ){
  // Set the maximum quadrant side length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  for ( int ip = node_face_ptr[node]; ip < node_face_ptr[node+1]; ip++ ){
    int face = node_face_conn[ip];
    if (face_owner != face){
      // Find the corresponding node index for this face
      for ( int j = 0; j < 4; j++ ){
        if (node == face_conn[4*face + j]){
          TMRQuadrant quad;
          quad.x = hmax*(j % 2);
          quad.y = hmax*(j/2);

          // Get the node array
          TMRQuadrantArray *nodes;
          quadtrees[face]->getNodes(&nodes);
          
          // Search the nodes and set the tag value
          const int use_node_search = 1;
          TMRQuadrant *t = nodes->contains(&quad, use_node_search);
          t->tag = p->tag;
        }
      }
    }
  }
}

/*
  Copy the edge node variable number on the edge with the given edge
  index to the adjacent edges.  These values are then communicated
  back to their owners in a second parallel step.

  input:
  edge:        the edge number that the node lies on
  edge_index:  the local index of the edge 0 <= edge index < 4
  face_owner:  the face index that owns the edge
  p:           the quadrant with the correct tag value
*/
void TMRQuadForest::copyEdgeNodes( int edge,
                                   int edge_index,
                                   int face_owner,
                                   TMRQuadrant *p ){
  // Set the maximum quadrant side length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Store the u coordinate along the edge
  int32_t u = 0;
  if (edge_index < 2){
    u = p->y;
  }
  else {
    u = p->x;
  }

  // Retrieve the first and second node numbers to determine the
  // relative orientation between this edge and each adjacent edge
  int n1 = face_conn[4*face_owner + face_to_edge_nodes[edge_index][0]];
  int n2 = face_conn[4*face_owner + face_to_edge_nodes[edge_index][1]];

  for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
    int face = edge_face_conn[ip]/4;
    if (face_owner != face){
      // Get the adjacent edge index on the opposite face
      int adj_index = edge_face_conn[ip] % 4;
      
      // Get the orientation 
      int nn1 = face_conn[4*face + face_to_edge_nodes[adj_index][0]];
      int nn2 = face_conn[4*face + face_to_edge_nodes[adj_index][1]];

      // Determine whether the edges are in the same direction or
      // are reversed
      int reverse = (n1 == nn2 && n2 == nn1);

      // Set the u-coordinate along the edge
      int32_t uquad = u;
      if (reverse){
        uquad = hmax - u;
      }
      
      // Transform the quadrant to the adjacent coordinate system
      TMRQuadrant quad;
      if (adj_index < 2){
        quad.x = hmax*(adj_index % 2);
        quad.y = uquad;
      }
      else {
        quad.x = uquad;
        quad.y = hmax*(adj_index % 2);
      }

      // Get the node array
      TMRQuadrantArray *nodes;
      quadtrees[face]->getNodes(&nodes);
      
      // Search the nodes and set the tag value
      const int use_node_search = 1;
      TMRQuadrant *t = nodes->contains(&quad, use_node_search);
      t->tag = p->tag;      
    }
  }
}

/*
  Copy the nodes ordered on adjacent face to the locally owned
  adjacent faces.

  This code loops over all the nodes in all the faces that are owned
  by this processor. When a node on a face, edge or corner is
  encountered that is owned by this processor, its node is copied to
  the non-local interface nodes. These nodes are later communicated to
  their source processors.

  Note that the dependent nodes must be labeled correctly before this
  function will work.

  intput:
  edge_owners:  the face index of the edge owner
  node_owners:  the face index of the corner/node owner
*/
void TMRQuadForest::copyAdjacentNodes( const int *edge_face_owners,
                                       const int *node_face_owners ){
  // Copy over the face node numbers
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Get the node array
    TMRQuadrantArray *nodes;
    quadtrees[face]->getNodes(&nodes);
        
    // Get the actual quadrant array
    int size;
    TMRQuadrant *array;
    nodes->getArray(&array, &size);
    
    // Loop over all the nodes and find the ones that are on the
    // desired face
    for ( int i = 0; i < size; i++ ){
      const int32_t hmax = 1 << TMR_MAX_LEVEL; 
      
      // Only copy the independent nodes that have a non-negative
      // tag node number
      if (array[i].tag >= 0){        
        // Determine which faces the quadrant touches if any
        int fx0 = (array[i].x == 0);
        int fy0 = (array[i].y == 0);
        int fx = (fx0 || array[i].x == hmax);
        int fy = (fy0 || array[i].y == hmax);
          
        if (fx && fy){
          int node_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2);
          int node = face_conn[4*face + node_index];
          if (node_face_owners[node] == face){
            copyCornerNodes(node, node_index, face, &array[i]);
          }
        }
        else if (fx || fy){
          int edge_index = 
            fx*(fx0 ? 0 : 1) + fy*(fy0 ? 2 : 3);
          int edge = face_edge_conn[4*face + edge_index];
          if (edge_face_owners[edge] == face){
            copyEdgeNodes(edge, edge_index, face, &array[i]);
          }
        }
      }
    }
  }
}

/*
  The following code sends the node numbers back to the face owners.

  This code sends the nodes that are owned by this processor but are
  also on face, edge or corner of an adjacent quadtree. These nodes have
  been ordered locally, but their node numbers have not yet been
  passed back to the corresponding face owner (that does not own the
  face/edge/node).

  This code scans over faces that are allocated locally, but are not
  owned locally - these are the partial quadtrees allocated during the
  recvQuadNeighbors call. The nodes on edges/faces/corners of these
  trees are in the local quadtree order.

  input:
  edge_owners:  the face index of the edge owner
  node_owners:  the face index of the corner/node owner
*/
void TMRQuadForest::sendNodeNeighbors( const int *edge_face_owners,
                                       const int *node_face_owners ){
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Send all the nodes that have a non-negative node number
  TMRQuadrantQueue **queues = new TMRQuadrantQueue*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      queues[k] = new TMRQuadrantQueue();
    }
    else {
      queues[k] = NULL;
    }
  }

  // Count up the number of nodes owned on each face
  for ( int face = 0; face < num_faces; face++ ){
    if (quadtrees[face] && 
        (mpi_rank != mpi_face_owners[face])){
      // The destination rank
      int rank = mpi_face_owners[face];

      // Get the element array
      TMRQuadrantArray *nodes;
      quadtrees[face]->getNodes(&nodes);
      
      // Get the actual quadrant array
      int size;
      TMRQuadrant *array;
      nodes->getArray(&array, &size);

      // Loop over all the nodes and check where they need to be sent
      for ( int i = 0; i < size; i++ ){
        if (array[i].tag >= 0){
          const int32_t hmax = 1 << TMR_MAX_LEVEL; 
          
          // Determine which faces the quadrant touches if any
          int fx0 = (array[i].x == 0);
          int fy0 = (array[i].y == 0);
          int fx = (fx0 || array[i].x == hmax);
          int fy = (fy0 || array[i].y == hmax);
        
          // Check whether the node lies on a corner or edge
          if (fx || fy){
            // This is cheating -- set the level -- which is not needed
            // for the nodes to the destination face number
            TMRQuadrant q = array[i];
            q.level = face;
            
            // Pass the face to any adjacent quadrees
            int face_owner = -1;
            if (fx && fy){
              int node_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2);
              int node = face_conn[4*face + node_index];
              face_owner = node_face_owners[node];
            }
            else if (fx){
              int edge_index = (fx0 ? 0 : 1);
              int edge = face_edge_conn[4*face + edge_index];
              face_owner = edge_face_owners[edge];
            }
            else if (fy){
              int edge_index = (fy0 ? 2 : 3);
              int edge = face_edge_conn[4*face + edge_index];
              face_owner = edge_face_owners[edge];
            }

            // If this processor owns the node/edge/face then it will
            // have the correct node number and it needs to be sent back
            // to the face who sent this in the first place
            if (face_owner >= 0 && 
                (mpi_rank == mpi_face_owners[face_owner])){            
              queues[rank]->push(&q);
            }
          }
        }
      }
    }
  }

  // Send everything back to the original senders
  MPI_Request *send_requests = new MPI_Request[ mpi_size ];

  // Send the nodes to the other trees
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
  TMRQuadrantQueue **qnodes = new TMRQuadrantQueue*[ num_faces ];
  memset(qnodes, 0, num_faces*sizeof(TMRQuadrantQueue*));

  // Receive the arrays of incoming nodes
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
      // For this one function, the level encodes the face owner for
      // the node. We now zero this so that there is no
      // confusion. This preserves the tag value which stores the node
      // number itself.
      int face = array[i].level;
      array[i].level = 0;
      if (!qnodes[face]){
        qnodes[face] = new TMRQuadrantQueue();
      }
      qnodes[face]->push(&array[i]);
    }
    delete [] array;
  }

  // Loop over all the faces that have been received 
  for ( int face = 0; face < num_faces; face++ ){
    if (qnodes[face]){
      // Convert the queue to an array and delete the old queue
      TMRQuadrantArray *new_nodes = qnodes[face]->toArray();
      delete qnodes[face];

      // Sort the new array so that we can pass the node numbers
      new_nodes->sort();

      // Get the nodes from the tree
      TMRQuadrantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Get the arrays of the actual quadrants
      int size, new_size;
      TMRQuadrant *array, *new_array;
      nodes->getArray(&array, &size);
      new_nodes->getArray(&new_array, &new_size);
      
      // Scan through the two sorted arrays
      for ( int i = 0, inew = 0; (i < size && inew < new_size); inew++ ){
        while (i < size && (array[i].compare(&new_array[inew]) < 0)){
          i++;
        }

        // If the arrays are equal, then copy over the tag value -- the
        // global node number
        if (array[i].compare(&new_array[inew]) == 0){
          if (array[i].tag >= 0){
            array[i].tag = new_array[inew].tag;
          }
        }
      }
      delete new_nodes;
    }
  }
  delete [] qnodes;

  // Wait for all the sends to complete
  MPI_Waitall(nsends, send_requests, MPI_STATUSES_IGNORE);
  delete [] send_requests;
  
  // Now free the arrays
  for ( int k = 0; k < mpi_size; k++ ){
    if (arrays[k]){
      delete arrays[k];
    }
  }
  delete [] arrays;
}

/*
  Create the nodes from the element mesh

  Note that the element mesh must be balanced before the nodes can be
  ordered.

  This function first computes the face that owns of each of the
  faces, edges and corners(/nodes) within the super-mesh. Next, the
  the non-local quadrants that border each quadree are passed back to the
  processors which own the quadree. This creates a layer of quadrants
  that are temporarily stored in a partial quadree. The code then
  creates nodes within each quadrant element for all of the quadrees
  (including the partial quadrees). Next, the dependent nodes (that are
  hanging on a face or an edge of an element)are labeled according to
  whether they are on an edge or face.
  
  After the dependent nodes are labeled, the nodes are ordered on all
  processors to form a complete global ordering. Next, the nodes are
  communicated locally across quadrant and partial quadrant faces, edges
  and corners. Finally, the new node numbers are returned to the
  processors that border the quadree owners. And lastly, the non-local
  partial quadrees are freed.  

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

  // Find the face numbers corresponding to the owner for each face,
  // edge and corner so that we know who should be ordering what!
  int *edge_face_owners = new int[ num_edges ];
  int *node_face_owners = new int[ num_nodes ];

  // Find the owners for the edges and nodes. The owner is
  // chosen as the connecting face with the lowest face number
  for ( int edge = 0; edge < num_edges; edge++ ){
    edge_face_owners[edge] = num_faces;

    int ipend = edge_face_ptr[edge+1];
    for ( int ip = edge_face_ptr[edge]; ip < ipend; ip++ ){
      int face = edge_face_conn[ip]/4;
      if (face < edge_face_owners[edge]){
        edge_face_owners[edge] = face;
      }
    }
  }

  // Find the node owners
  for ( int node = 0; node < num_nodes; node++ ){
    node_face_owners[node] = num_faces;
    
    int ipend = node_face_ptr[node+1];
    for ( int ip = node_face_ptr[node]; ip < ipend; ip++ ){
      if (node_face_conn[ip] < node_face_owners[node]){
        node_face_owners[node] = node_face_conn[ip];
      }
    }
  }

  // Compute the neighboring quadrants 
  recvQuadNeighbors();

  // Label the dependent nodes on all the faces
  computeDepEdges();

  // Allocate all possible nodes on all of the trees, including the
  // partial trees that have just been exchanged.
  for ( int face = 0; face < num_faces; face++ ){
    if (quadtrees[face]){
      quadtrees[face]->createNodes(mesh_order);
    }
  }

  // Label the dependent nodes
  labelDependentNodes();

  // Compute the global ordering of the nodes
  orderGlobalNodes(edge_face_owners, node_face_owners);

  // Copy the adjacent node numbers between local faces
  copyAdjacentNodes(edge_face_owners, node_face_owners);

  // Send the nodes (complete with global number) back to whence they
  // came to finish the global ordering
  sendNodeNeighbors(edge_face_owners, node_face_owners);

  // Free and NULL the trees that are not local - these are not
  // required anymore.
  for ( int face = 0; face < num_faces; face++ ){
    if (quadtrees[face] && 
        (mpi_rank != mpi_face_owners[face])){
      delete quadtrees[face];
      quadtrees[face] = NULL;
    }
  }
  
  delete [] edge_face_owners;
  delete [] node_face_owners;
}

/*
  Create the local mesh connectivity for the portion of the
  finite-element mesh stored on this processor.

  Before calling the create mesh function, you must have already
  allocated and ordered the nodes

  output:
  conn:    the connectivity (using global node numbers)
  nelems:  the number of elements in the mesh
*/
void TMRQuadForest::createMeshConn( int **_conn, int *_nelems ){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // First, count up the number of elements
  int nelems = 0;
  for ( int face = 0; face < num_faces; face++ ){
    if (mpi_rank == mpi_face_owners[face]){
      nelems += quadtrees[face]->getNumElements();
    }
  }

  // Allocate the connectivity
  int conn_size = 0;
  int *elem_conn = new int[ mesh_order*mesh_order*nelems ];

  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Retrieve the element array and the node array
    TMRQuadrantArray *elements, *nodes;
    quadtrees[face]->getElements(&elements);
    quadtrees[face]->getNodes(&nodes);

    // Get the current array of quadrants
    int size;
    TMRQuadrant *array;
    elements->getArray(&array, &size);

    for ( int i = 0; i < size; i++ ){
      // For all searches/comparisons, we use node numbers
      const int use_nodes = 1;
          
      // Add all of the nodes from this element
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level - (mesh_order-2));
      TMRQuadrant p;
      p.level = array[i].level;
          
      for ( int jj = 0; jj < mesh_order; jj++ ){
        for ( int ii = 0; ii < mesh_order; ii++ ){
          p.x = array[i].x + ii*h;
          p.y = array[i].y + jj*h;
            
          // Get the index for the node
          TMRQuadrant *t = nodes->contains(&p, use_nodes);
          
          // Reset the node level if the element is smaller
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

  // Set the output arrays
  *_conn = elem_conn;
  *_nelems = nelems;
}

/*
  Get the dependent connectivity information (create it if it has not
  been allocated previously).

  output:
  ptr:      pointer for each dependent node number
  conn:     connectivity to each (global) independent node
  weights:  the weight values for each dependent node
*/
int TMRQuadForest::getDepNodeConn( const int **ptr, const int **conn,
                                   const double **weights ){
  if (!dep_ptr){
    if (dep_ptr){ delete [] dep_ptr; }
    if (dep_conn){ delete [] dep_conn; }
    if (dep_weights){ delete [] dep_weights; }
    createDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);
  }

  if (ptr){ *ptr = dep_ptr; }
  if (conn){ *conn = dep_conn; }
  if (weights){ *weights = dep_weights; }
  return num_dep_nodes;
}

/*
  Create the dependent mesh information for all local dependent
  nodes. 

  output:
  ptr:      pointer for each dependent node number
  conn:     connectivity to each (global) independent node
  weights:  the weight values for each dependent node
*/
void TMRQuadForest::createDepNodeConn( int **_ptr, int **_conn, 
                                       double **_weights ){
  // Get the MPI rank
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Count up the total number of dependent faces
  int ndep_edges = 0;
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int ndep = 0;
    dep_edges[owned]->getArray(NULL, &ndep);
    ndep_edges += ndep;
  }

  // Allocate the space for the dependent variable information
  const int nodes_per_edge = (2*mesh_order-1)*(2*mesh_order-1);
  int *edge_nodes = new int[ nodes_per_edge*ndep_edges ];

  // Retrieve the face information from the dependent faces
  int face_offset = 0, edge_offset = 0;
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Get the nodes associated with this quadtree
    TMRQuadrantArray *nodes;
    quadtrees[face]->getNodes(&nodes);

    // Get the dependent edges
    int dep_size;
    TMRQuadrant *dep_array;
    dep_edges[owned]->getArray(&dep_array, &dep_size);

    for ( int i = 0; i < dep_size; i++, edge_offset++ ){
      TMRQuadrant *b = &dep_array[i];

      // Get the edge index
      int edge_index = b->tag;

      // Find the edge length of the element quadrant
      const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);
      const int32_t hc = 1 << (TMR_MAX_LEVEL - b->level - (mesh_order-1));
      
      // Loop over all the nodes on this face
      int *en = &edge_nodes[nodes_per_edge*edge_offset];
      for ( int ii = 0; ii < (2*mesh_order-1); ii++ ){
        // Get the node location corresponding to this face
        TMRQuadrant node;
        if (edge_index < 2){
          node.x = b->x + h*(edge_index % 2);
          node.y = b->y + ii*hc;
        }
        else {
          node.x = b->x + ii*hc;
          node.y = b->y + h*(edge_index % 2);
        }

        // Search for node quadrant and record its label
        const int use_node_search = 1;
        TMRQuadrant *t = nodes->contains(&node, use_node_search);
        en[0] = t->tag;
        en++;
      }
    }
  }

  // Allocate the pointer array
  int *ptr = new int[ num_dep_nodes+1 ];
  memset(ptr, 0, (num_dep_nodes+1)*sizeof(int));

  for ( int i = 0; i < ndep_edges; i++ ){
    // Loop over all the nodes on the dependent edges
    const int *en = &edge_nodes[nodes_per_edge*i];
    for ( int ii = 1; ii < 2*mesh_order-1; ii += 2 ){
      int node = -en[ii]-1;
      ptr[node+1] = mesh_order;
    }
  }

  // Count up the offsets to the pointers
  for ( int i = 0; i < num_dep_nodes; i++ ){
    ptr[i+1] += ptr[i];
  }

  // Set the connectivities and weights
  int *conn = new int[ ptr[num_dep_nodes] ];
  double *weights = new double[ ptr[num_dep_nodes] ];

  const double wt2[] = {0.5, 0.5};
  const double wt31[] = {0.375, 0.75, -0.125};
  const double wt32[] = {-0.125, 0.75, 0.375};

  const double *wt[2] = {NULL, NULL};
  if (mesh_order == 2){
    wt[0] = wt[1] = wt2;
  }
  else { // mesh_order == 3
    wt[0] = wt31;
    wt[1] = wt32;
  }

  for ( int k = 0; k < ndep_edges; k++ ){
    // Loop over all the nodes on the dependent edges
    const int *en = &edge_nodes[nodes_per_edge*k];
    for ( int ii = 1; ii < 2*mesh_order-1; ii += 2 ){
      int node = -en[ii]-1;
      for ( int i = 0; i < mesh_order; i++ ){
        conn[ptr[node] + i] = en[2*i];
        weights[ptr[node] + i] = wt[ii/2][i];
      }
    }
  }
  
  // Free the edge node array
  delete [] edge_nodes;

  // Return the weights
  *_ptr = ptr;
  *_conn = conn;
  *_weights = weights;
}

/*
  Add the node weights to the array. This eliminates dependent nodes
  by unrolling the dependency information.

  input:
  t:             the quadrant node
  w:             the weight for this node
  cdep_ptr:      the dependent node pointer for the coarse mesh
  cdep_conn:     the dependent node connectivity for the coarse mesh
  cdep_weights:  the dependent node weights for the coarse mesh

  input/output: 
  weights:       the list of index weights
  nweights:      the number of weights in the list
*/
void TMRQuadForest::addNodeWeights( TMRQuadrant *t, double w,
                                    const int *cdep_ptr,
                                    const int *cdep_conn,
                                    const double *cdep_weights,
                                    TMRIndexWeight *weights, 
                                    int *nweights ){
  if (t->tag >= 0){
    weights[*nweights].index = t->tag;
    weights[*nweights].weight = w;
    (*nweights)++;
  }
  else {
    // Unravel the dependent node connectivity
    int node = -t->tag-1;
    for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
      weights[*nweights].index = cdep_conn[jp];
      weights[*nweights].weight = w*cdep_weights[jp];
      (*nweights)++;
    }
  }
}

/*
  Create the interpolation operator from the coarse to the fine mesh.

  Each processor builds the interpolation for the locally owned nodes
  in the mesh. This interpolation will refer to other nodes, but the
  output is local.
  
  input:
  coarse:   the coarse quadtree forest that has the same layout as this
  
  output:
  ptr:      the pointer into the local rows
  conn:     the connectivity using global numbers
  weights:  the interpolation weights for each point
*/
void TMRQuadForest::createInterpolation( TMRQuadForest *coarse,
                                         int **_interp_ptr,
                                         int **_interp_conn,
                                         double **_interp_weights ){
  // Set the pointers to NULL
  *_interp_ptr = NULL;
  *_interp_conn = NULL;
  *_interp_weights = NULL;

  // Set the node weights
  const double wt2[] = {0.5, 0.5};
  const double wt31[] = {0.375, 0.75, -0.125};
  const double wt32[] = {-0.125, 0.75, 0.375};
  const double *wt[2] = {NULL, NULL};
  if (mesh_order == 2){
    wt[0] = wt[1] = wt2;
  }
  else { // mesh_order == 3
    wt[0] = wt31;
    wt[1] = wt32;
  }

  // Get the dependent node information
  const int *cdep_ptr, *cdep_conn;
  const double *cdep_weights;
  coarse->getDepNodeConn(&cdep_ptr, &cdep_conn, &cdep_weights);

  // Get the MPI rank
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Determine the number of fine nodes
  int nnodes = node_range[mpi_rank+1] - node_range[mpi_rank];
  
  // Set the current and maximum lengths of the connectivity array
  int conn_len = 0;
  int max_conn_len = 
    mesh_order*mesh_order*nnodes;

  // Allocate the arrays that store the interpolation data/weights
  int *pos = new int[ nnodes ];
  int *flags = new int[ nnodes ];
  memset(flags, 0, nnodes*sizeof(int));

  // Record the connectivity for these nodes
  int *ptr = new int[ nnodes+1 ];
  int *conn = new int[ max_conn_len ];
  double *weights = new double[ max_conn_len ];

  // Set the number of interpolated nodes that have been set
  ptr[0] = 0;

  // The maximum possible size of the array of weights. Note
  // that this is found if every node is a dependent node (which is
  // impossible) which points to a dependent face node (also
  // impossible). It is an upper bound.
  int max_size = (mesh_order*mesh_order)*(mesh_order);
  TMRIndexWeight *wlist = new TMRIndexWeight[ max_size ];

  // Loop over all the faces and compute the interpolation
  int count = 0;
  for ( int owned = 0; owned < num_owned_faces; owned++ ){
    int face = owned_faces[owned];

    // Get the nodes from the fine and coarse mesh
    TMRQuadrantArray *nodes, *coarse_nodes;
    quadtrees[face]->getNodes(&nodes);
    coarse->quadtrees[face]->getNodes(&coarse_nodes);

    // Get the size of the fine node array
    int fine_size;
    TMRQuadrant *fine;
    nodes->getArray(&fine, &fine_size);

    for ( int i = 0; i < fine_size; i++ ){
      // Use a node-based search
      const int use_node_search = 1;

      // Keep track of the number of new weight elements
      int nweights = 0;
      
      // Get the node number
      int node_num = fine[i].tag;
      int loc = node_num - node_range[mpi_rank];

      // Loop over all the adjacent nodes
      if ((node_num >= node_range[mpi_rank] && 
           node_num < node_range[mpi_rank+1]) && 
          flags[loc] == 0){
        // Flag that we've reached this node and store the true
        // location of the interpolation based on the node number
        flags[loc] = 1;
        pos[count] = loc;

        // Retrieve the quadrant pointer
        TMRQuadrant *t = coarse_nodes->contains(&fine[i], use_node_search);
        if (t){
          addNodeWeights(t, 1.0, cdep_ptr, cdep_conn, cdep_weights,
                         wlist, &nweights);
        }
        else {
          // Get the element size for coarse element mesh
          const int32_t h = 2*(1 << (TMR_MAX_LEVEL - fine[i].level));

          // The node spacing for the fine mesh
          const int32_t hf = 
            1 << (TMR_MAX_LEVEL - fine[i].level - (mesh_order-2));
          
          // The node spacing for the coarse mesh
          const int32_t hc = 2*hf;

          // Compute the offsets to the node locations
          int32_t px = fine[i].x % h;
          int32_t py = fine[i].y % h;

          // Determine which interpolation to use
          int32_t sx = fine[i].x % hc;
          int32_t sy = fine[i].y % hc;

          // Add to the interpolation depending on the values of px,
          // py, and pz
          if (sx && sy){
            for ( int jj = 0; jj < mesh_order; jj++ ){
              for ( int ii = 0; ii < mesh_order; ii++ ){
                TMRQuadrant node = fine[i];
                node.x = (fine[i].x - px) + hc*ii;
                node.y = (fine[i].y - py) + hc*jj;
                t = coarse_nodes->contains(&node, use_node_search);
                double w = wt[px/hc][ii]*wt[py/hc][jj];
                addNodeWeights(t, w, cdep_ptr, cdep_conn, cdep_weights,
                               wlist, &nweights);
              }
            }
          }
          else if (sx){
            for ( int ii = 0; ii < mesh_order; ii++ ){
              TMRQuadrant node = fine[i];
              node.x = (fine[i].x - px) + hc*ii;
              t = coarse_nodes->contains(&node, use_node_search);
              addNodeWeights(t, wt[px/hc][ii], cdep_ptr, cdep_conn, cdep_weights,
                             wlist, &nweights);
            }
          }
          else if (sy){
            for ( int jj = 0; jj < mesh_order; jj++ ){
              TMRQuadrant node = fine[i];
              node.y = (fine[i].y - py) + hc*jj;
              t = coarse_nodes->contains(&node, use_node_search);
              addNodeWeights(t, wt[py/hc][jj], cdep_ptr, cdep_conn, cdep_weights,
                             wlist, &nweights);
            }
          }
        }

        // Sort the dependent weight values
        nweights = TMRIndexWeight::uniqueSort(wlist, nweights);
      
        // Check whether adding the new weights will exceed the size
        // of the allocated array
        if (ptr[count] + nweights >= max_conn_len){
          // Estimate the final size of the array
          int estimate = (1 + 2*ptr[count]/count)*(nnodes - count);
          max_conn_len += estimate;

          // Allocate a new array of the new estimated size
          int *temp = new int[ max_conn_len ];
          memcpy(temp, conn, conn_len*sizeof(int));
          delete [] conn;
          conn = temp;
        
          // Copy over the weight values as well
          double *tempw = new double[ max_conn_len ];
          memcpy(tempw, weights, conn_len*sizeof(double));
          delete [] weights;
          weights = tempw;
        }

        // Extract the weights from the sorted list
        for ( int k = 0; k < nweights; k++, conn_len++ ){
          conn[conn_len] = wlist[k].index;
          weights[conn_len] = wlist[k].weight;
        }

        // Increment the pointer
        count++;
        ptr[count] = conn_len;
      }
    }
  }

  // Free the weights and flags
  delete [] wlist;
  delete [] flags;

  // Allocate new interpolation arrays with tight bounds
  int *interp_ptr = new int[ nnodes+1 ];
  int *interp_conn = new int[ conn_len ];
  double *interp_weights = new double[ conn_len ];

  // Set the pointer array so that it offsets into the true global
  // ordering of the nodes
  interp_ptr[0] = 0;
  for ( int k = 0; k < nnodes; k++ ){
    interp_ptr[pos[k]+1] = ptr[k+1] - ptr[k];
  }
  for ( int k = 1; k < nnodes+1; k++ ){
    interp_ptr[k] += interp_ptr[k-1];
  }

  // Set the new interpolation values
  for ( int k = 0; k < nnodes; k++ ){
    int loc = pos[k];
    int len = ptr[k+1] - ptr[k];
    for ( int j = 0; j < len; j++ ){
      interp_conn[interp_ptr[loc] + j] = conn[ptr[k] + j];
      interp_weights[interp_ptr[loc] + j] = weights[ptr[k] + j];
    }
  }

  delete [] pos;
  delete [] ptr;
  delete [] conn;
  delete [] weights;

  // Set the return arrays
  *_interp_ptr = interp_ptr;
  *_interp_conn = interp_conn;
  *_interp_weights = interp_weights;
}

