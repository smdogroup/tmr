#include "TMROctForest.h"

/*
  Map from a block edge number to the local node numbers
*/
const int block_to_edge_nodes[][2] = {{0,1}, {2,3}, {4,5}, {6,7}, 
                                      {0,2}, {1,3}, {4,6}, {5,7}, 
                                      {0,4}, {1,5}, {2,6}, {3,7}};

/*
  Map from the face number to the block number
*/
const int block_to_face_nodes[][4] = {{0,1,2,3}, 
                                      {4,5,6,7},
                                      {0,1,4,5},
                                      {2,3,6,7},
                                      {0,2,4,6},
                                      {1,3,5,7}};

/*
  All possible orientations for two connecting faces
*/
const int face_orientations[][4] = {{0,1,2,3},
                                    {2,0,3,1},
                                    {3,2,1,0},
                                    {1,3,0,2}};

/*
  Given the integer face_index -- the index into the face_orientations
  array -- return the relative x/y locations on this face.

  2 --- 3 
  |     |
  0 --- 1
    
  The face can have any of the following orientations
    
  2 --- 3    3 --- 1    1 --- 0    0 --- 2
  |     |    |     |    |     |    |     |  
  0 --- 1    2 --- 0    3 --- 2    1 --- 3
  
  2 --- 0    0 --- 1    1 --- 3    3 --- 2
  |     |    |     |    |     |    |     |  
  3 --- 1    2 --- 3    0 --- 2    1 --- 0
*/
inline void get_face_node_coords( const int face_id, 
                                  const int32_t u, const int32_t v,
                                  int32_t *x, int32_t *y ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  
  // For these orientations, the normals are in the same direction
  if (face_id == 0){
    // {0,1,2,3} - the block is aligned with the face
    //  2 --- 3
    //  |     |
    //  0 --- 1
    *x = u;
    *y = v;
  }
  else if (face_id == 1){
    // {2,0,3,1} - block rotated to the right
    // 3 --- 1
    // |     |
    // 2 --- 0 
    *x = hmax - v;
    *y = u;
  }
  else if (face_id == 2){
    // {3,2,1,0}
    // 1 --- 0
    // |     |
    // 3 --- 2
    *x = hmax - u;
    *y = hmax - v;
  }
  else if (face_id == 3){
    // {1,3,0,2}
    // 0 --- 2
    // |     |
    // 1 --- 3
    *x = v;
    *y = hmax - u;
  }
}

/*
  Create the TMROctForest object
*/
TMROctForest::TMROctForest( MPI_Comm _comm ){
  // Set the MPI communicator
  comm = _comm;

  // Set the range of nodes
  node_range = NULL;

  // Zero out the nodes/edges/faces and all data
  num_nodes = 0;
  num_edges = 0;
  num_faces = 0;
  num_blocks = 0;

  // Set all the unallocated pointers to NULL
  block_conn = NULL;
  block_face_conn = NULL;
  block_edge_conn = NULL;
  node_block_ptr = NULL;
  node_block_conn = NULL;
  edge_block_ptr = NULL;
  edge_block_conn = NULL;
  face_block_ptr = NULL;
  face_block_conn = NULL;
  octrees = NULL;
  mpi_block_owners = NULL;
  node_range = NULL;

  // Set the size of the mesh
  mesh_order = 0;
}

/*
  Free the data allocated by the TMROctForest object
*/
TMROctForest::~TMROctForest(){
  if (block_conn){ delete [] block_conn; }
  if (block_face_conn){ delete [] block_face_conn; }
  if (block_edge_conn){ delete [] block_edge_conn; }
  if (node_block_ptr){ delete [] node_block_ptr; }
  if (node_block_conn){ delete [] node_block_conn; }
  if (edge_block_ptr){ delete [] edge_block_ptr; }
  if (edge_block_conn){ delete [] edge_block_conn; }
  if (face_block_ptr){ delete [] face_block_ptr; }
  if (face_block_conn){ delete [] face_block_conn; }
  if (octrees){ 
    for ( int i = 0; i < num_blocks; i++ ){
      if (octrees[i]){ delete octrees[i]; }
    }
    delete [] octrees;
  }
  if (mpi_block_owners){ delete [] mpi_block_owners; }
  if (node_range){ delete [] node_range; }
}

/*
  Set the connectivity of the blocks

  This call is collective on all processors. Every processor must make
  a call with the same connectivity information, otherwise the
  inter-octree information will be inconsistent. This sets the
  connectivity and then performs a partition of the mesh using METIS.

  This code sets the block connectivity and generates the following
  additional data that are required:

  1. Block to node connectivity (input)
  2. Node to block connectivity (required for corner balancing)
  3. Unique edge ordering
  4. Unique face ordering
  5. Edge to block connectivity
  6. Face to block connectivity

  This information is required for creating octree forests on the
  unstructured super mesh.
*/
void TMROctForest::setConnectivity( int _num_nodes,
                                     const int *_block_conn,
                                     int _num_blocks ){
  // Free any data if it has already been allocated. 
  // This will erase everything internally.
  if (block_conn){ delete [] block_conn; }
  if (block_face_conn){ delete [] block_face_conn; }
  if (block_edge_conn){ delete [] block_edge_conn; }
  if (node_block_ptr){ delete [] node_block_ptr; }
  if (node_block_conn){ delete [] node_block_conn; }
  if (edge_block_ptr){ delete [] edge_block_ptr; }
  if (edge_block_conn){ delete [] edge_block_conn; }
  if (face_block_ptr){ delete [] face_block_ptr; }
  if (face_block_conn){ delete [] face_block_conn; }
  if (octrees){ 
    for ( int i = 0; i < num_blocks; i++ ){
      if (octrees[i]){ delete octrees[i]; }
    }
    delete [] octrees;
  }
  if (mpi_block_owners){ delete [] mpi_block_owners; }
  if (node_range){ delete [] node_range; }

  // Copy over the data locally
  num_nodes = _num_nodes;
  num_edges = 0;
  num_faces = 0;
  num_blocks = _num_blocks;

  // Copy over the block connectivity
  block_conn = new int[ 8*num_blocks ];
  memcpy(block_conn, _block_conn, 8*num_blocks*sizeof(int));

  // Set up the partition using metis
  mpi_block_owners = new int[ num_blocks ];

  // Compute the partition on the root processor and
  // broadcast the result to all other processors
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);  
  MPI_Comm_size(comm, &mpi_size);  
  if (mpi_rank == 0){
    // For now.. just perform an assignment
    for ( int i = 0; i < num_blocks; i++ ){
      mpi_block_owners[i] = i % mpi_size;
    }
  }
  
  // Broadcast the face owners to all processors
  MPI_Bcast(mpi_block_owners, num_blocks, MPI_INT, 0, comm);

  // Create the data structure for the node to block connectivity
  node_block_ptr = new int[ num_nodes+1 ];
  memset(node_block_ptr, 0, (num_nodes+1)*sizeof(int));

  // Count the number of times each node is referred to
  for ( int i = 0; i < 8*num_blocks; i++ ){
    node_block_ptr[block_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_nodes+1; i++ ){
    node_block_ptr[i] += node_block_ptr[i-1];
  }

  // Allocate the full node to face pointer array
  node_block_conn = new int[ node_block_ptr[num_nodes] ];
  for ( int i = 0; i < num_blocks; i++ ){
    for ( int j = 8*i; j < 8*(i+1); j++ ){
      int node = block_conn[j];
      node_block_conn[node_block_ptr[node]] = i;
      node_block_ptr[node]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_nodes; i >= 1; i-- ){
    node_block_ptr[i] = node_block_ptr[i-1];
  }
  node_block_ptr[0] = 0;

  // Now establish a unique ordering of the edges along each block
  // -------------------------------------------------------------
  block_edge_conn = new int[ 12*num_blocks ];
  for ( int i = 0; i < 12*num_blocks; i++ ){
    block_edge_conn[i] = -1;
  }

  int edge = 0;
  for ( int i = 0; i < num_blocks; i++ ){
    // Loop over each edge on this block
    for ( int j = 0; j < 12; j++ ){
      if (block_edge_conn[12*i + j] < 0){
        int n1 = block_conn[8*i + block_to_edge_nodes[j][0]];
        int n2 = block_conn[7*i + block_to_edge_nodes[j][1]];

        // Keep track of the number of edges found
        const int max_nedges = 128;
        int edge_index[max_nedges];
        int nedges = 1;
        edge_index[0] = 12*i + j;

        // Set the edge number - if any is found
        int edge_num = -1;

        // Scan through the blocks that share the same node and check
        // if any of the edges are also shared 
        for ( int ip = node_block_ptr[n1];
              ip < node_block_ptr[n1+1]; ip++ ){
          int ii = node_block_conn[ip];
          
          // Loop over each edge in the new block
          for ( int jj = 0; jj < 12; jj++ ){
            int nn1 = block_conn[12*ii + block_to_edge_nodes[jj][0]];
            int nn2 = block_conn[12*ii + block_to_edge_nodes[jj][1]];

            // Check if the block matches
            if ((n1 == nn1 && n2 == nn2) ||
                (n1 == nn2 && n2 == nn1)){
              if (block_edge_conn[12*ii + jj] >= 0){
                // If this edge has been ordered, copy over
                // the edge number
                edge_num = block_edge_conn[12*ii + jj];
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
          block_edge_conn[edge_index[ii]] = edge_num;
        }
      }
    }
  }

  // Set the total number of edges
  num_edges = edge;
 
  // Create the data structure for the edge to block connectivity
  edge_block_ptr = new int[ num_edges+1 ];
  memset(edge_block_ptr, 0, (num_edges+1)*sizeof(int));

  // Count the number of times each edge is referred to
  for ( int i = 0; i < 12*num_blocks; i++ ){
    edge_block_ptr[block_edge_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_edges+1; i++ ){
    edge_block_ptr[i] += edge_block_ptr[i-1];
  }

  // Allocate the full edge to block pointer array
  edge_block_conn = new int[ edge_block_ptr[num_edges] ];
  for ( int i = 0; i < num_blocks; i++ ){
    for ( int j = 12*i; j < 12*(i+1); j++ ){
      int e = block_edge_conn[j];
      edge_block_conn[edge_block_ptr[e]] = i;
      edge_block_ptr[e]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_edges; i >= 1; i-- ){
    edge_block_ptr[i] = edge_block_ptr[i-1];
  }
  edge_block_ptr[0] = 0;

  // Now establish a unique ordering of the faces for each block
  // -----------------------------------------------------------
  block_face_conn = new int[ 6*num_blocks ];
  for ( int i = 0; i < 6*num_blocks; i++ ){
    block_face_conn[i] = -1;
  }

  // Count up and order the number of faces in the mesh
  int face = 0;
  for ( int i = 0; i < num_blocks; i++ ){
    // Loop over all the face for this block
    for ( int j = 0; j < 6; j++ ){
      if (block_face_conn[6*i + j] < 0){
        // Get the face nodes for the j-th face
        int face_nodes[4];
        for ( int k = 0; k < 4; k++ ){
          face_nodes[k] = block_conn[8*i + block_to_face_nodes[j][k]];
        }

        // Keep track of the number of connecting faces
        const int max_nfaces = 128;
        int face_index[max_nfaces];
        int nfaces = 1;
        face_index[0] = 6*i + j;

        // Set the face number
        int face_num = -1;
        
        // Scan through and find the blocks that share a common
        // node that is also on this face
        int node = block_conn[8*i + block_to_face_nodes[j][0]];
        for ( int ip = node_block_ptr[node];
              ip < node_block_ptr[node+1]; ip++ ){
          int ii = node_block_conn[ip];

          // Loop over all the faces for the ii-th block
          int face_equiv = 0;
          for ( int jj = 0; jj < 6; jj++ ){
            // Get the nodes corresponding to this face
            int new_face_nodes[4];
            for ( int k = 0; k < 4; k++ ){
              new_face_nodes[k] = 
                block_conn[8*ii + block_to_face_nodes[jj][k]];
            }

            // Loop over the relative orientations between this block
            // and the last block
            for ( int ort = 0; ort < 4; ort++ ){
              face_equiv = 1;
              for ( int k = 0; k < 4; k++ ){
                if (face_nodes[k] != 
                    new_face_nodes[face_orientations[ort][k]]){
                  face_equiv = 0;
                  break;
                }
              }
              if (face_equiv){
                break;
              }
            }

            if (face_equiv){
              if (block_face_conn[6*ii + jj] >= 0){
                face_num = block_face_conn[6*ii + jj];
              }
              else if (nfaces < max_nfaces){
                face_index[nfaces] = 6*ii + jj;
                nfaces++;
              }
            }
          }
        }

        // If this face does not have an associated face number
        // give it one
        if (face_num < 0){
          face_num = face;
          face++;
        }

        // Set the face numbers
        for ( int ii = 0; ii < nfaces; ii++ ){
          block_face_conn[face_index[ii]] = face_num;
        }
      }
    }
  }

  // Set the number of faces
  num_faces = face;

  // Create the data structure for the face to block connectivity
  face_block_ptr = new int[ num_faces+1 ];
  memset(face_block_ptr, 0, (num_faces+1)*sizeof(int));

  // Count the number of times each face is referred to
  for ( int i = 0; i < 12*num_blocks; i++ ){
    face_block_ptr[block_face_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_faces+1; i++ ){
    face_block_ptr[i] += face_block_ptr[i-1];
  }

  // Allocate the full face to block pointer array
  face_block_conn = new int[ face_block_ptr[num_faces] ];
  for ( int i = 0; i < num_blocks; i++ ){
    for ( int j = 6*i; j < 6*(i+1); j++ ){
      int f = block_face_conn[j];
      face_block_conn[face_block_ptr[f]] = i;
      face_block_ptr[f]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_faces; i >= 1; i-- ){
    face_block_ptr[i] = face_block_ptr[i-1];
  }
  face_block_ptr[0] = 0;
}

/*
  Create a forest with the specified refinement level
*/
void TMROctForest::createTrees( int refine_level ){
  if (octrees){ 
    for ( int i = 0; i < num_blocks; i++ ){
      if (octrees[i]){ delete octrees[i]; }
    }
    delete [] octrees;
  }

  // Create the octrees
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  octrees = new TMROctree*[ num_blocks ];
  memset(octrees, 0, num_blocks*sizeof(TMROctree*));
  for ( int i = 0; i < num_blocks; i++ ){
    if (mpi_block_owners[i] == mpi_rank){
      octrees[i] = new TMROctree(refine_level);
    }
  }
}

/*
  Create a random set of trees

  This function is usually used for testing purposes.
*/
void TMROctForest::createRandomTrees( int nrand, 
                                      int min_level, int max_level ){
  if (octrees){ 
    for ( int i = 0; i < num_blocks; i++ ){
      if (octrees[i]){ delete octrees[i]; }
    }
    delete [] octrees;
  }

  // Create a random set of octrees
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  octrees = new TMROctree*[ num_faces ];
  memset(octrees, 0, num_faces*sizeof(TMROctree*));
  for ( int i = 0; i < num_faces; i++ ){
    if (mpi_block_owners[i] == mpi_rank){
      octrees[i] = new TMROctree(nrand, min_level, max_level);
    }
  }
}

/*
  Duplicate the forest

  This function creates a duplicate representation of the current
  forest. This function copies the global connectivity of the forest
  and copies each individual tree.
*/
TMROctForest* TMROctForest::duplicate(){
  TMROctForest *dup = new TMROctForest(comm);
  if (octrees){
    // Copy over the connectivity data 
    dup->num_nodes = num_nodes;
    dup->num_edges = num_edges;
    dup->num_faces = num_faces;
    dup->num_blocks = num_blocks;

    // Allocate/copy the block connectivities
    dup->block_conn = new int[ 8*num_blocks ];
    dup->block_face_conn = new int[ 6*num_faces ];
    dup->block_edge_conn = new int[ 12*num_faces ];
    memcpy(dup->block_conn, block_conn, 8*num_blocks*sizeof(int));
    memcpy(dup->block_face_conn, block_face_conn, 6*num_blocks*sizeof(int));
    memcpy(dup->block_edge_conn, block_edge_conn, 12*num_blocks*sizeof(int));
    
    // Allocate/copy the inverse relationships
    dup->node_block_ptr = new int[ num_nodes+1 ];
    dup->node_block_conn = new int[ node_block_ptr[num_nodes] ];
    memcpy(dup->node_block_ptr, node_block_ptr, (num_nodes+1)*sizeof(int));
    memcpy(dup->node_block_conn, node_block_conn, 
           node_block_ptr[num_nodes]*sizeof(int));

    dup->edge_block_ptr = new int[ num_edges+1 ];
    dup->edge_block_conn = new int[ edge_block_ptr[num_edges] ];
    memcpy(dup->edge_block_ptr, edge_block_ptr, (num_edges+1)*sizeof(int));
    memcpy(dup->edge_block_conn, edge_block_conn, 
           edge_block_ptr[num_edges]*sizeof(int));

    dup->face_block_ptr = new int[ num_faces+1 ];
    dup->face_block_conn = new int[ face_block_ptr[num_faces] ];
    memcpy(dup->face_block_ptr, face_block_ptr, (num_faces+1)*sizeof(int));
    memcpy(dup->face_block_conn, face_block_conn, 
           face_block_ptr[num_faces]*sizeof(int));

    // Allocate/copy the block ownership
    dup->mpi_block_owners = new int[ num_blocks ];
    memcpy(dup->mpi_block_owners, mpi_block_owners, num_blocks*sizeof(int));

    // Duplicate all the octrees
    dup->octrees = new TMROctree*[ num_blocks ];
    memset(dup->octrees, 0, num_blocks*sizeof(TMROctree*));
    for ( int i = 0; i < num_blocks; i++ ){
      if (octrees[i]){
        TMROctantArray *elements;
        octrees[i]->getElements(&elements);
        dup->octrees[i] = new TMROctree(elements->duplicate());
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
TMROctForest* TMROctForest::coarsen(){
  TMROctForest *coarse = new TMROctForest(comm);
  if (octrees){
    // Copy over the connectivity data 
    coarse->num_nodes = num_nodes;
    coarse->num_edges = num_edges;
    coarse->num_faces = num_faces;
    coarse->num_blocks = num_blocks;

    // Allocate/copy the block connectivities
    coarse->block_conn = new int[ 8*num_blocks ];
    coarse->block_face_conn = new int[ 6*num_faces ];
    coarse->block_edge_conn = new int[ 12*num_faces ];
    memcpy(coarse->block_conn, block_conn, 8*num_blocks*sizeof(int));
    memcpy(coarse->block_face_conn, block_face_conn, 6*num_blocks*sizeof(int));
    memcpy(coarse->block_edge_conn, block_edge_conn, 12*num_blocks*sizeof(int));
    
    // Allocate/copy the inverse relationships
    coarse->node_block_ptr = new int[ num_nodes+1 ];
    coarse->node_block_conn = new int[ node_block_ptr[num_nodes] ];
    memcpy(coarse->node_block_ptr, node_block_ptr, (num_nodes+1)*sizeof(int));
    memcpy(coarse->node_block_conn, node_block_conn, 
           node_block_ptr[num_nodes]*sizeof(int));

    coarse->edge_block_ptr = new int[ num_edges+1 ];
    coarse->edge_block_conn = new int[ edge_block_ptr[num_edges] ];
    memcpy(coarse->edge_block_ptr, edge_block_ptr, (num_edges+1)*sizeof(int));
    memcpy(coarse->edge_block_conn, edge_block_conn, 
           edge_block_ptr[num_edges]*sizeof(int));

    coarse->face_block_ptr = new int[ num_faces+1 ];
    coarse->face_block_conn = new int[ face_block_ptr[num_faces] ];
    memcpy(coarse->face_block_ptr, face_block_ptr, (num_faces+1)*sizeof(int));
    memcpy(coarse->face_block_conn, face_block_conn, 
           face_block_ptr[num_faces]*sizeof(int));

    // Allocate/copy the block ownership
    coarse->mpi_block_owners = new int[ num_blocks ];
    memcpy(coarse->mpi_block_owners, mpi_block_owners, num_blocks*sizeof(int));

    // Coarselicate all the octrees
    coarse->octrees = new TMROctree*[ num_blocks ];
    memset(coarse->octrees, 0, num_blocks*sizeof(TMROctree*));
    for ( int i = 0; i < num_blocks; i++ ){
      if (octrees[i]){
        coarse->octrees[i] = octrees[i]->coarsen();
      }
    }
  }

  return coarse;
}

/*
  Add the face neighbors for an adjacent tree

  input:
  tree:   the block


*/
void TMROctForest::addFaceNeighbors( int block,
                                     int face_index, 
                                     TMROctant p,
                                     TMROctantHash **hash,
                                     TMROctantQueue **queue ){
  /*
  // Determine the global face number
  int face = block_face_conn[6*block + face_index];
 
  // Get the u/v coordinates of this face
  int32_t u, v;
  
  if (face_index < 2){ // x-face
    get_face_node_coords(face_id, p.y, p.z, &u, &v);
  }
  else if (face_index < 4){ // y-face
    get_face_node_coords(face_id, p.x, p.z, &u, &v);
  }
  else { // z-face
    get_face_node_coords(face_id, p.x, p.y, &u, &v);
  }

  // Loop over all the adjacent faces and add the block
  for ( int ip = face_block_ptr[face];
        ip < _block_ptr[face+1]; ip++ ){

    // Get the adjacent block
    int adjacent = face_block_conn[ip];
    if (adjacent != block){
      if (!hash[adjacent]){
        hash[adjacent] = new TMROctantHash();
        queue[adjacent] = new TMROctantQueue();
      }

      // Find the corresponding face index for this element
      for ( int j = 0; j < 6; j++ ){
        // Check whether this is the matching face
        if (block_face_conn[6*adjacent + j] == face){
          TMROctant neighbor;

          // Set the element base
          if (j < 2){
            neighbor.x = hmax*(j % 2);
          }
          else if (j < 4){
            neighbor.y = hmax*(j % 2);
          }
          else {
            neighbor.z = hmax*(j % 2);
          }

          // Add the octant to the list
          if (hash[adjacent]->addOctant(&neighbor)){
            queue[adjacent]->push(&neighbor);
          }
        }
      }
    }
  }
  */
}

/*
  Add the edge neighbors for a adjacent trees

  This function is called to balance the forest across tree edges.
  Given an octant p on the specified corner index, this code ensures
  a corner balanced tree, by adding the corresponding corner qudrant
  to all node-adjacent faces. If the octant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  face:    the local face index (p is defined on this face)
  corner:  the tree corner index (p must lie on this corner)
  p:       the local octant
  hash:    the array of hash objects for each face
  queue:   the array of newly added qudrants for each face
*/
void TMROctForest::addEdgeNeighbors( int block,
                                     int edge_index, 
                                     TMROctant p,
                                     TMROctantHash **hash,
                                     TMROctantQueue **queue ){
  // First determine the global edge number 
  int edge = block_edge_conn[12*block + edge_index];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Store the u coordinate along the edge
  int32_t ucoord = 0;
  if (edge_index < 4){
    ucoord = p.x;
  }
  else if (edge_index < 8){
    ucoord = p.y;
  }
  else {
    ucoord = p.z;
  }

  // Retrieve the first and second node numbers
  int n1 = block_conn[8*block + block_to_edge_nodes[edge_index][0]];
  int n2 = block_conn[8*block + block_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent faces
  for ( int ip = edge_block_ptr[edge];
        ip < edge_block_ptr[edge+1]; ip++ ){

    // Get the faces that are adjacent across this edge
    int adjacent = edge_block_conn[ip];
    if (adjacent != block){
      if (!hash[adjacent]){
        hash[adjacent] = new TMROctantHash();
        queue[adjacent] = new TMROctantQueue();
      }

      // Loop over all the edges on the adjacent block
      for ( int j = 0; j < 12; j++ ){
        int nn1 = block_conn[8*adjacent + block_to_edge_nodes[j][0]];
        int nn2 = block_conn[8*adjacent + block_to_edge_nodes[j][1]];

        // Add the octant to the list
        int forward = (n1 == nn1 && n2 == nn2);
        int reverse = (n1 == nn2 && n2 == nn1);

        if (forward || reverse){
          int32_t u = ucoord;
          if (reverse){
            u = hmax - 2*h - ucoord;
          }

          TMROctant neighbor;
          neighbor.level = p.level;
          if (j < 4){
            neighbor.x = u;
            neighbor.y = hmax*(j % 2);
            neighbor.z = hmax*(j/2);
          }
          else if (j < 8){
            neighbor.x = hmax*(j % 2);
            neighbor.y = u;
            neighbor.z = hmax*((j-4)/2);
          }
          else {
            neighbor.x = hmax*(j % 2);
            neighbor.y = hmax*((j-8)/2);
            neighbor.z = u;
          }
         
          // Add the octant to the list
          if (hash[adjacent]->addOctant(&neighbor)){
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
  Given an octant p on the specified corner index, this code ensures
  a corner balanced tree, by adding the corresponding corner qudrant
  to all node-adjacent faces. If the octant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  block:   the local block index (p is defined on this face)
  corner:  the corner index (p must lie on this corner)
  p:       the local octant
  hash:    the array of hash objects for each face
  queue:   the array of newly added qudrants for each face
*/
void TMROctForest::addCornerNeighbors( int block,
                                       int corner,
                                       TMROctant p,
                                       TMROctantHash **hash,
                                       TMROctantQueue **queue ){
  // First determine the global edge number 
  int node = block_conn[8*block + corner];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = node_block_ptr[node];
        ip < node_block_ptr[node+1]; ip++ ){
      
    // Get the faces that are adjacent across this edge
    int adjacent = node_block_conn[ip];
    if (adjacent != block){

      // Loop over all of the corners in this adjacent block
      for ( int j = 0; j < 8; j++ ){
        if (block_conn[8*adjacent+j] == node){
          // Allocate the octant hash/queue if not already allocated
          if (!hash[adjacent]){
            hash[adjacent] = new TMROctantHash();
            queue[adjacent] = new TMROctantQueue();
          }

          // Compute the octant location
          TMROctant neighbor;
          neighbor.level = p.level;
          neighbor.x = (hmax - 2*h)*(j%2);
          neighbor.y = (hmax - 2*h)*(j/2);
          neighbor.z = (hmax - 2*h)*(j/4);

          // Add the octant to the list
          if (hash[adjacent]->addOctant(&neighbor)){
            queue[adjacent]->push(&neighbor);
          }
        }
      }
    }
  }
}

/*
  Balance the octant on the entire octree

  This code finds the 0-parent of all adjacent octants either within
  the current tree or within an adjacent tree and adds those octants
  to balance the input octant 'oct'.

  input:
  block:           index of the block for the octant
  oct:             the octant itself
  hash:            the array of hash tables for each block
  queue:           the octant queues for each block
  balance_corner:  balance across corners 
  balance_tree:    balance on the entire tree
*/
void TMROctForest::balanceOctant( int block,
                                  TMROctant *oct,
                                  TMROctantHash **hash,
                                  TMROctantQueue **queue,
                                  const int balance_corner,
                                  const int balance_tree ){
  // Local octant data
  TMROctant p, neighbor, q;

  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  
  // Get the parent of the octant, and add the their
  // block-matched octants from each block, as long 
  // as they fall within the bounds
  if (oct->level > 1){
    oct->parent(&p);
        
    // Balance across each face
    for ( int face = 0; face < 6; face++ ){
      p.faceNeighbor(face, &neighbor);
      neighbor.getSibling(0, &q);
      
      // If we are within bounds, add the neighbor
      if ((q.x >= 0 && q.x < hmax) &&
          (q.y >= 0 && q.y < hmax) &&
          (q.z >= 0 && q.z < hmax)){
        if (hash[block]->addOctant(&q)){
          queue[block]->push(&q);
        }
      }
      else if (balance_tree){
        addFaceNeighbors(block, face, q, hash, queue);
      }
    }

    // Add the edge-adjacent elements
    for ( int edge = 0; edge < 12; edge++ ){
      p.edgeNeighbor(edge, &neighbor);
      neighbor.getSibling(0, &q);
	  
      // If we're in bounds, add the neighbor
      if ((q.x >= 0 && q.x < hmax) &&
          (q.y >= 0 && q.y < hmax) &&
          (q.z >= 0 && q.z < hmax)){
        if (hash[block]->addOctant(&q)){
          queue[block]->push(&q);
        }
      }
      else if (balance_tree){
        // The node may lie across an edge or face
        int ex = (q.x < 0 || q.x >= hmax);
        int ey = (q.y < 0 || q.y >= hmax);
        int ez = (q.z < 0 || q.z >= hmax);
        
        if ((ex && ey) || (ex && ez) || (ey && ez)){
          // The octant lies along a true edge
          addEdgeNeighbors(block, edge, q, hash, queue);
        }
        else {
          // The octant actually lies along a face 
          int face = 0;
          if (ex){ // This is an x-face
            face = (q.x < 0 ? 0 : 1);
          }
          else if (ey){ // This is a y-face
            face = (q.y < 0 ? 2 : 3);
          }
          else { // This is a z-face
            face = (q.z < 0 ? 4 : 5);
          }
          addFaceNeighbors(block, face, q, hash, queue);
        }
      }
    }

    // If we're balancing across edges and 
    if (balance_corner){
      for ( int corner = 0; corner < 8; corner++ ){
        p.cornerNeighbor(corner, &neighbor);
        neighbor.getSibling(0, &q);
        
        if ((q.x >= 0 && q.x < hmax) &&
            (q.y >= 0 && q.y < hmax) &&
            (q.z >= 0 && q.z < hmax)){
          if (hash[block]->addOctant(&q)){
            queue[block]->push(&q);
          }
        }
        else if (balance_tree){
          // The node may lie across a corner, edge or face
          int ex = (q.x < 0 || q.x >= hmax);
          int ey = (q.y < 0 || q.y >= hmax);
          int ez = (q.z < 0 || q.z >= hmax);

          if (ex && ey && ez){
            // Add the octant to the other trees
            addCornerNeighbors(block, corner, neighbor, hash, queue);
          }
          else if ((ex && ey) || (ex && ez) || (ey && ez)){
            // The octant lies along a true edge
            int edge = 0;
            if (ey && ez){
              edge = (q.y < 0 ? 0 : 1);
              edge += (q.z < 0 ? 0 : 2);
            }
            else if (ex && ez){
              edge = (q.x < 0 ? 4 : 5);
              edge += (q.z < 0 ? 0 : 2);
            }
            else { // (ex && ey)
              edge = (q.x < 0 ? 8 : 9);
              edge += (q.y < 0 ? 0 : 2);
            }
            addEdgeNeighbors(block, edge, q, hash, queue);
          }
          else {
            // The octant actually lies along a face 
            int face = 0;
            if (ex){ // This is an x-face
              face = (q.x < 0 ? 0 : 1);
            }
            else if (ey){ // This is a y-face
              face = (q.y < 0 ? 2 : 3);
            }
            else { // This is a z-face
              face = (q.z < 0 ? 4 : 5);
            }
            addFaceNeighbors(block, face, q, hash, queue);
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
void TMROctForest::balance( int balance_corner ){
  // Get the mpi rank
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Create a hash table for the balanced tree
  TMROctantHash **hash = new TMROctantHash*[ num_blocks ];
  TMROctantQueue **queue = new TMROctantQueue*[ num_blocks ];
  
  // Set the hash tables and queues to NULL
  memset(hash, 0, num_blocks*sizeof(TMROctantHash*));
  memset(queue, 0, num_blocks*sizeof(TMROctantQueue*));

  for ( int block = 0; block < num_blocks; block++ ){
    if (mpi_block_owners[block] == mpi_rank){
      // Allocate the hash and queue if they are not already
      // allocated on this processor
      if (!hash[block]){
        hash[block] = new TMROctantHash();
        queue[block] = new TMROctantQueue();
      }

      // Get the current array of octants
      TMROctantArray *elements;
      octrees[block]->getElements(&elements);

      // Get the array of elements
      int size;
      TMROctant *array;
      elements->getArray(&array, &size);
    
      // Add all the elements
      for ( int i = 0; i < size; i++ ){
        TMROctant oct;
        array[i].getSibling(0, &oct);
        hash[block]->addOctant(&oct);

        // Balance the octant by push the neighbors
        const int balance_tree = 1;
        balanceOctant(block, &oct, hash, queue, 
                      balance_corner, balance_tree);
      }
    }
  }

  int queue_length_flag = 1;
  while (queue_length_flag){
    // Now continue until the queue of added octants is
    // empty. At each iteration, pop an octant and add 
    // its neighbours until nothing new is added. This code
    // handles the propagation of octants to adjacent octants.
    for ( int block = 0; block < num_blocks; block++ ){
      if (queue[block]){
        while (queue[block]->length() > 0){
          TMROctant oct = queue[block]->pop();
          // Balance the octant by push the neighbors
          const int balance_tree = 1;
          balanceOctant(block, &oct, hash, queue, 
                        balance_corner, balance_tree);
        }
      }
    }
      
    queue_length_flag = 0;
    for ( int block = 0; block < num_blocks; block++ ){
      if (queue[block] &&
          queue[block]->length() > 0){
        queue_length_flag = 1;
        break;
      }
    }
  }

  // Free the queues - they are no longer required
  for ( int block = 0; block < num_blocks; block++ ){
    if (mpi_block_owners[block] != mpi_rank && queue[block]){ 
      delete queue[block];
      queue[block] = NULL;
    }
  }

  // Now everything is locally balanced - all the elements 
  // on the current processor are balanced with all the other
  // elements on the current processor, but nothing is 
  // inter-processor balanced yet.

  // The number of octants that will be sent from this processor
  // to all other processors in the communicator
  int *oct_counts = new int[ mpi_size ];
  memset(oct_counts, 0, mpi_size*sizeof(int));

  // Create a list of new queues - these will only be allocated
  // as required
  TMROctantQueue **send_queues = new TMROctantQueue*[ mpi_size ];
  memset(send_queues, 0, mpi_size*sizeof(TMROctantQueue*));

  for ( int block = 0; block < num_blocks; block++ ){
    if (hash[block] && 
        (mpi_block_owners[block] != mpi_rank)){
      // Create a sorted list of local the 0-child octants.
      // This can be further reduced to limit the amount of 
      // memory passed between processors
      TMROctantArray *elems0 = hash[block]->toArray();
      elems0->sort();

      // Get the array of 0-octants 
      int size;
      TMROctant *array;
      elems0->getArray(&array, &size);

      // Check that the destination queue has been allocated
      int dest_rank = mpi_block_owners[block];
      if (!send_queues[dest_rank]){
        send_queues[dest_rank] = new TMROctantQueue();
      }
      
      // Add the octants to their destination queue
      TMROctantQueue *dest = send_queues[dest_rank];
      
      if (size > 0){
        // Get the parent of the octant
        TMROctant p, s = array[0];
        s.parent(&p);

        // Loop over all 
        for ( int i = 0; i < size; i++ ){
          if (!p.contains(&array[i])){
            s.tag = block;
            dest->push(&s);
          }
          // Get the next octant and find its parent
          s = array[i];
          s.parent(&p);
        }

        // Push the last octant onto the queue
        s.tag = block;
        dest->push(&s);
      }
      
      // Update the number of octants
      oct_counts[dest_rank] = dest->length();
      
      // Free the elements and the hash table
      delete elems0;
      delete hash[block];
      hash[block] = NULL;
    }
  }

  // Now distribute the octants to their destination octrees and
  // balance the corresponding octrees including the new elements.
  int *oct_recv_counts = new int[ mpi_size ];
  MPI_Alltoall(oct_counts, 1, MPI_INT,
               oct_recv_counts, 1, MPI_INT, comm);

  // Add the recieved octants into the local queues
  TMROctantArray **send_arrays = new TMROctantArray*[ mpi_size ];
  TMROctantArray **recv_arrays = new TMROctantArray*[ mpi_size ];
  memset(send_arrays, 0, mpi_size*sizeof(TMROctantArray*));
  memset(recv_arrays, 0, mpi_size*sizeof(TMROctantArray*));

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
      // Conver the octant array to a list
      send_arrays[i] = send_queues[i]->toArray();
    
      // Delete the send queue 
      delete send_queues[i];

      // Get the octant array
      int size;
      TMROctant *array;
      send_arrays[i]->getArray(&array, &size);

      // Post the send to the destination
      MPI_Isend(array, size, TMROctant_MPI_type, 
                i, 0, comm, &send_request[nsend]);
      nsend++;
    }
  }

  delete [] send_queues;

  // Loop over the recieve calls
  int nrecv = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (oct_recv_counts[i] > 0){
      int size = oct_recv_counts[i];
      TMROctant *array = new TMROctant[ size ];
      recv_arrays[i] = new TMROctantArray(array, size);
      
      // Post the recieve
      MPI_Irecv(array, size, TMROctant_MPI_type,
                i, 0, comm, &recv_request[i]);
      nrecv++;
    }
  }

  // Wait for any recieve to complete
  for ( int i = 0; i < nrecv; i++ ){
    int index = 0;
    MPI_Waitany(mpi_size, recv_request, &index, MPI_STATUS_IGNORE);

    // Push the octants on to their corresponding queues
    int size;
    TMROctant *array;
    recv_arrays[index]->getArray(&array, &size);
    for ( int j = 0; j < size; j++ ){
      int block = array[j].tag;
      if (hash[block]->addOctant(&array[j])){
        queue[block]->push(&array[j]);
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
  delete [] oct_counts;
  delete [] oct_recv_counts;
  delete [] send_request;
  delete [] recv_request;

  // Now all the received octants will balance the tree locally
  // without having to worry about off-processor octants.
  for ( int block = 0; block < num_blocks; block++ ){
    if (queue[block]){
      while (queue[block]->length() > 0){
        const int balance_tree = 0;
        TMROctant oct = queue[block]->pop();
        balanceOctant(block, &oct, hash, queue, 
                      balance_corner, balance_tree);
      }
    }
  }

  // Free the remaining queues - these will be the ones that are
  // owned by this processor.
  for ( int block = 0; block < num_blocks; block++ ){
    if (queue[block]){ delete queue[block]; }
  }
  delete [] queue;

  // Convert the hash tables back to the elements
  for ( int block = 0; block < num_blocks; block++ ){
    if (mpi_block_owners[block] == mpi_rank){
      // Now convert the elements from child-0 elements to
      // elements which cover the full mesh
      TMROctantArray *child0_elems = hash[block]->toArray();
      int size;
      TMROctant *array;
      child0_elems->getArray(&array, &size);

      // Loop over all elements and add their siblings
      for ( int i = 0; i < size; i++ ){
        if (array[i].level > 0){
          for ( int j = 0; j < 8; j++ ){
            TMROctant q;
            array[i].getSibling(j, &q);
            hash[block]->addOctant(&q);
          }
        }
      }

      // Free the temporary elements
      delete child0_elems;

      // Set the elements into the quadtree
      TMROctantArray *elements = hash[block]->toArray();
      elements->sort();
      octrees[block]->setElements(elements);

      // Free the corresponding hash
      delete hash[block];
    }
  }

  delete [] hash;
}
