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
const int block_to_face_nodes[][4] = {{0,2,4,6}, {1,3,5,7},
                                      {0,1,4,5}, {2,3,6,7},
                                      {0,1,2,3}, {4,5,6,7}};

/*
  All possible orientations for two connecting faces
*/
const int face_orientations[][4] = {{0,1,2,3},
                                    {2,0,3,1},
                                    {3,2,1,0},
                                    {1,3,0,2}};

/*
  Given the x/y locations on the face with orientation face_id,
  find/return the corresponding u/v coordinates on the transformed
  face.
*/
inline void get_face_oct_coords( const int face_id, const int32_t h,
                                 const int32_t x, const int32_t y,
                                 int32_t *u, int32_t *v ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;  
  if (face_id == 0){
    *u = x;
    *v = y;
  }
  else if (face_id == 1){
    *u = hmax - h - y;
    *v = x;
  }
  else if (face_id == 2){
    *u = hmax - h - x;
    *v = hmax - h - y;
  }
  else {
    *u = y;
    *v = hmax - h - x;
  }
}

/*
  Given the u/v locations on the root face, set the coordinates on the
  destination face.
*/
inline void set_face_oct_coords( const int face_id, const int32_t h,
                                 const int32_t u, const int32_t v,
                                 int32_t *x, int32_t *y ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;  
  if (face_id == 0){
    *x = u;
    *y = v;
  }
  else if (face_id == 1){
    *x = v;
    *y = hmax - h - u;
  }
  else if (face_id == 2){
    *x = hmax - h - u;
    *y = hmax - h - v; 
  }
  else {
    *x = hmax - h - v;
    *y = u; 
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
  block_face_ids = NULL;
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
  if (block_face_ids){ delete [] block_face_ids; }
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
  if (block_face_ids){ delete [] block_face_ids; }
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
  node_range = NULL;

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
        int n2 = block_conn[8*i + block_to_edge_nodes[j][1]];

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
            int nn1 = block_conn[8*ii + block_to_edge_nodes[jj][0]];
            int nn2 = block_conn[8*ii + block_to_edge_nodes[jj][1]];

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
                edge_index[nedges] = 12*ii + jj;
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
  for ( int i = 0; i < 6*num_blocks; i++ ){
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

  // Compute the block face ids for each block face
  block_face_ids = new int[ 6*num_blocks ];

  // Determine the face block ids. These ids store both the
  // orientation of the face relative to the face owner - the
  // face owner is on the block with the lowest block index.
  for ( int face = 0; face < num_faces; face++ ){
    int block_owner = num_blocks;
    int owner_index = 0;

    // Scan through the blocks pointing to this face to determine
    // the block owner - the block with the lowest index
    for ( int ip = face_block_ptr[face]; 
          ip < face_block_ptr[face+1]; ip++ ){
      int block = face_block_conn[ip];
      if (block < block_owner){
        block_owner = block;

        // Find the new owner index
        owner_index = 0;
        for ( int j = 0; j < 6; j++, owner_index++ ){
          if (block_face_conn[6*block + j] == face){
            break;
          }
        }
      }
    }

    // Get the nodes for the owner face - these will be used 
    // to establish the relative orientation between the owner
    // and each subsequent block
    int owner_nodes[4];
    for ( int k = 0; k < 4; k++ ){
      owner_nodes[k] = block_conn[8*block_owner + 
                                  block_to_face_nodes[owner_index][k]];
    }

    // Now determine the relative orientations of the faces and store
    // them in the index
    for ( int ip = face_block_ptr[face]; 
          ip < face_block_ptr[face+1]; ip++ ){
      // Find the block index
      int block = face_block_conn[ip];
      for ( int j = 0; j < 6; j++ ){
        if (face == block_face_conn[6*block + j]){
          // Now, find the relative orientations between the two faces.
          // First get the nodes corresponding to this face
          int new_face_nodes[4];
          for ( int k = 0; k < 4; k++ ){
            new_face_nodes[k] = 
              block_conn[8*block + block_to_face_nodes[j][k]];
          }
          
          // Loop over the relative orientations between this block
          // and the last block
          for ( int ort = 0; ort < 4; ort++ ){
            int face_equiv = 1;
            for ( int k = 0; k < 4; k++ ){
              if (owner_nodes[k] != 
                  new_face_nodes[face_orientations[ort][k]]){
                face_equiv = 0;
                break;
              }
            }
            
            // Set the orientation and break out of this loop
            if (face_equiv){
              block_face_ids[6*block + j] = ort;
              break;
            }
          }
        }
      }
    }
  }
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
  octrees = new TMROctree*[ num_blocks ];
  memset(octrees, 0, num_blocks*sizeof(TMROctree*));
  for ( int i = 0; i < num_blocks; i++ ){
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
    dup->block_face_ids = new int[ 6*num_faces ];
    dup->block_edge_conn = new int[ 12*num_faces ];
    memcpy(dup->block_conn, block_conn, 8*num_blocks*sizeof(int));
    memcpy(dup->block_face_conn, block_face_conn, 6*num_blocks*sizeof(int));
    memcpy(dup->block_face_ids, block_face_ids, 6*num_blocks*sizeof(int));
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
    coarse->block_face_ids = new int[ 6*num_faces ];
    coarse->block_edge_conn = new int[ 12*num_faces ];
    memcpy(coarse->block_conn, block_conn, 8*num_blocks*sizeof(int));
    memcpy(coarse->block_face_conn, block_face_conn, 6*num_blocks*sizeof(int));
    memcpy(coarse->block_face_ids, block_face_ids, 6*num_blocks*sizeof(int));
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
  // Determine the global face number
  int face = block_face_conn[6*block + face_index];
  int face_id = block_face_ids[6*block + face_index];
 
  // Set the size of the side-length of the octant
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Get the u/v coordinates of this face
  int32_t u, v;
  if (face_index < 2){ // x-face
    get_face_oct_coords(face_id, 2*h, p.y, p.z, &u, &v);
  }
  else if (face_index < 4){ // y-face
    get_face_oct_coords(face_id, 2*h, p.x, p.z, &u, &v);
  }
  else { // z-face
    get_face_oct_coords(face_id, 2*h, p.x, p.y, &u, &v);
  }

  // Loop over all the adjacent faces and add the block
  for ( int ip = face_block_ptr[face]; 
        ip < face_block_ptr[face+1]; ip++ ){

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
          // Get the face id for the matching face
          face_id = block_face_ids[6*adjacent + j];

          // Get the neighboring octant on the face
          TMROctant neighbor;
          neighbor.level = p.level;
          if (j < 2){
            neighbor.x = (hmax - 2*h)*(j % 2);
            set_face_oct_coords(face_id, h, u, v, 
                                &neighbor.y, &neighbor.z);
          }
          else if (j < 4){
            neighbor.y = (hmax - 2*h)*(j % 2);
            set_face_oct_coords(face_id, h, u, v, 
                                &neighbor.x, &neighbor.z);
          }
          else {
            neighbor.z = (hmax - 2*h)*(j % 2);
            set_face_oct_coords(face_id, h, u, v, 
                                &neighbor.x, &neighbor.y);
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
            neighbor.y = (hmax - 2*h)*(j % 2);
            neighbor.z = (hmax - 2*h)*(j/2);
          }
          else if (j < 8){
            neighbor.x = (hmax - 2*h)*(j % 2);
            neighbor.y = u;
            neighbor.z = (hmax - 2*h)*((j-4)/2);
          }
          else {
            neighbor.x = (hmax - 2*h)*(j % 2);
            neighbor.y = (hmax - 2*h)*((j-8)/2);
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
          neighbor.y = (hmax - 2*h)*((j % 4)/2);
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

  This algorithm uses a hash and a queue to balance the forest of
  octrees. For each element in the octree, we add the neighbors that
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


/*
  Add the quadrant to the processor queues corresponding to the
  non-local blocks that touch the given edge
*/
void TMROctForest::addEdgeOctantToQueues( const int edge, 
                                          const int mpi_rank, 
                                          TMROctant *q,
                                          TMROctantQueue **queues ){
  for ( int ip = edge_block_ptr[edge]; ip < edge_block_ptr[edge+1]; ip++ ){
    int block = edge_block_conn[ip];
    int rank = mpi_block_owners[block];
    if (rank != mpi_rank){
      queues[rank]->push(q);
    }
  }
}

/*
  Add the quadrant to the processor queues corresponding to the
  non-local blocks that touch the given face
*/
void TMROctForest::addFaceOctantToQueues( const int face,
                                          const int mpi_rank,
                                          TMROctant *q,
                                          TMROctantQueue **queues ){
  for ( int ip = face_block_ptr[face]; ip < face_block_ptr[face+1]; ip++ ){
    int block = face_block_conn[ip];
    int rank = mpi_block_owners[block];
    if (rank != mpi_rank){
      queues[rank]->push(q);
    }
  }
}

/*
  The following code exchanges the neighboring octants for each
  locally owned octree within the forest.

  This code exchanges non-local octants across each local octree
  information so that we can locally query octants on adjacent octrees
  without having to perform parallel communication.

  Note that this code creates partial non-local octrees that are
  adjacent to the local octrees in the forest. These partial local
  octrees should be freed after the nodal ordering has been computed.

  input:
  face_owners:   the index of the owning block for each face
  edge_owners:   the index of the owning block for each edge
*/
void TMROctForest::exchangeOctNeighbors( const int *face_block_owners,
                                         const int *edge_block_owners ){ 
  // Get the MPI information
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Allocate the queues that store the octants destined for each of
  // the processors
  TMROctantQueue **queues = new TMROctantQueue*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      queues[k] = new TMROctantQueue();
    }
    else {
      queues[k] = NULL;
    }
  }

  // For each block, determine the edge, face and 
  for ( int block = 0; block < num_blocks; block++ ){
    if (mpi_block_owners[block] == mpi_rank){
      // Flag to indicate whether any edge/face is non-local
      int has_non_local = 0;

      // Check if any of the edges has to be sent to another processor
      for ( int k = 0; k < 12; k++ ){
        int edge = block_edge_conn[12*block + k];
        for ( int ip = edge_block_ptr[edge];
              ip < edge_block_ptr[edge+1]; ip++ ){
          int dest_block = edge_block_conn[ip];
          if (dest_block != mpi_rank){
            has_non_local = 1;
            break;
          }
        }
      }

      if (!has_non_local){
        // If necessary, check if any of the faces should be sent to
        // another processor
        for ( int k = 0; k < 6; k++ ){
          int face = block_face_conn[6*block + k];
          for ( int ip = face_block_ptr[face];
                ip < face_block_ptr[face+1]; ip++ ){
            int dest_block = face_block_conn[ip];
            if (dest_block != mpi_rank){
              has_non_local = 1;
              break;
            }
          }
        }
      }

      if (has_non_local){
        // Get the element array
        TMROctantArray *elements;
        octrees[block]->getElements(&elements);
        
        // Get the actual octant array
        int size;
        TMROctant *array;
        elements->getArray(&array, &size);
        
        // Loop over all the elements and check where we need to send
        // the octants that are along each edge/face
        for ( int i = 0; i < size; i++ ){
          const int32_t hmax = 1 << TMR_MAX_LEVEL; 
          const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
          
          // Determine which faces the octant touches if any
          int fx0 = (array[i].x == 0);
          int fy0 = (array[i].y == 0);
          int fz0 = (array[i].z == 0);
          int fx = (fx0 || array[i].x + h == hmax);
          int fy = (fy0 || array[i].y + h == hmax);
          int fz = (fz0 || array[i].z + h == hmax);
          
          // Check whether the octant lies along an edge
          if (fx || fy || fz){
            // Copy the octant from the array and set its tag 
            // to the source block index
            TMROctant q = array[i];
            q.tag = block;
            
            // Pass the block to any adjacent octrees
            if (fy && fz){
              int edge_index = (fy0 ? 0 : 1) + (fz0 ? 0 : 2);
              int edge = block_edge_conn[12*block + edge_index];
              addEdgeOctantToQueues(edge, mpi_rank, &q, queues);
            }
            if (fx && fz){
              int edge_index = (fx0 ? 4 : 5) + (fz0 ? 0 : 2);
              int edge = block_edge_conn[12*block + edge_index];
              addEdgeOctantToQueues(edge, mpi_rank, &q, queues);
            }
            if (fx && fy){
              int edge_index = (fx0 ? 8 : 9) + (fy0 ? 0 : 2);
              int edge = block_edge_conn[12*block + edge_index];
              addEdgeOctantToQueues(edge, mpi_rank, &q, queues);
            }
            if (fx){
              int face_index = (fx0 ? 0 : 1);
              int face = block_face_conn[6*block + face_index];
              addFaceOctantToQueues(face, mpi_rank, &q, queues);
            }
            if (fy){
              int face_index = (fy0 ? 2 : 3);
              int face = block_face_conn[6*block + face_index];
              addFaceOctantToQueues(face, mpi_rank, &q, queues);
            }
            if (fz){
              int face_index = (fz0 ? 4 : 5);
              int face = block_face_conn[6*block + face_index];
              addFaceOctantToQueues(face, mpi_rank, &q, queues);
            }
          }
        }
      }
    }
  }

  // Allocate the requests
  MPI_Request *send_requests = new MPI_Request[ mpi_size ];

  // Send the octants to the other trees
  int nsends = 0;
  TMROctantArray **arrays = new TMROctantArray*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      // Create the arrays
      arrays[k] = queues[k]->toArray();
      delete queues[k];

      // Get the array of octants
      int size;
      TMROctant *array;
      arrays[k]->getArray(&array, &size);

      // Set the array of octants to their destination
      MPI_Isend(array, size, TMROctant_MPI_type, 
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
  TMROctantQueue **qtrees = new TMROctantQueue*[ num_blocks ];
  memset(qtrees, 0, num_blocks*sizeof(TMROctantQueue*));

  // Receive the arrays of incoming octants
  for ( int k = 0; k < mpi_size-1; k++ ){
    // Probe the recieved messages
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    // Retrieve the size and information for the incoming message
    int source = status.MPI_SOURCE;
    int tag = status.MPI_TAG;
    int size = 0;
    MPI_Get_count(&status, TMROctant_MPI_type, &size);

    // Allocate the incoming array
    TMROctant *array = new TMROctant[ size ];
    MPI_Recv(array, size, TMROctant_MPI_type,
             source, tag, comm, MPI_STATUS_IGNORE);

    // Push the octants into their corresponding trees
    for ( int i = 0; i < size; i++ ){
      int block = array[i].tag;
      if (!qtrees[block]){
        qtrees[block] = new TMROctantQueue();
      }
      qtrees[block]->push(&array[i]);
    }

    delete [] array;
  }

  // Now that the queues are completed allocate the faces
  for ( int block = 0; block < num_blocks; block++ ){
    if (qtrees[block]){
      octrees[block] = new TMROctree(qtrees[block]->toArray());
      delete qtrees[block];
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
  Create the nodes from the element mesh

  Note that the element mesh must be balanced before the nodes can be
  ordered.

  input:
  order:   the order of the mesh
*/
void TMROctForest::createNodes( int order ){
  // Check that the order falls within allowable bounds
  mesh_order = order;
  if (order > 3){ mesh_order = 3; }
  if (order < 2){ mesh_order = 2; }

  // Get the MPI information
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Find the block numbers corresponding to the owner for each face,
  // edge and corner so that we know who should be ordering what!
  int *face_block_owners = new int[ num_faces ];
  int *edge_block_owners = new int[ num_edges ];
  int *node_block_owners = new int[ num_nodes ];

  // Find the owners for the faces, edges and nodes. The owner is
  // chosen as the connecting block with the lowest block number
  for ( int face = 0; face < num_faces; face++ ){
    face_block_owners[face] = num_blocks;

    int ipend = face_block_ptr[face+1];
    for ( int ip = face_block_ptr[face]; ip < ipend; ip++ ){
      if (face_block_conn[ip] < face_block_owners[face]){
        face_block_owners[face] = face_block_conn[ip];
      }
    }
  }

  // Find the edge owners
  for ( int edge = 0; edge < num_edges; edge++ ){
    edge_block_owners[edge] = num_blocks;

    int ipend = edge_block_ptr[edge+1];
    for ( int ip = edge_block_ptr[edge]; ip < ipend; ip++ ){
      if (edge_block_conn[ip] < edge_block_owners[edge]){
        edge_block_owners[edge] = edge_block_conn[ip];
      }
    }
  }

  // Find the node owners
  for ( int node = 0; node < num_nodes; node++ ){
    node_block_owners[node] = num_blocks;
    
    int ipend = node_block_ptr[node+1];
    for ( int ip = node_block_ptr[node]; ip < ipend; ip++ ){
      if (node_block_conn[ip] < node_block_owners[node]){
        node_block_owners[node] = node_block_conn[ip];
      }
    }
  }

  // Compute the neighboring octants 
  exchangeOctNeighbors(face_block_owners, edge_block_owners);

  // Allocate all possible nodes on all of the trees, including the
  // partial trees that have just been exchanged.
  for ( int block = 0; block < num_blocks; block++ ){
    if (octrees[block]){
      octrees[block]->createNodes(mesh_order);
    }
  }


  /*

  // Determine the dependent nodes for each face without labeling
  // the dependent nodes on the edges yet.
  for ( int face = 0; face < num_faces; face++ ){
    if (quadtrees[face]){
      // Get the octant elements
      TMROctantArray *elements, *nodes;
      quadtrees[face]->getElements(&elements);
      quadtrees[face]->getNodes(&nodes);
      
      // Get the elements themselves
      int size;
      TMROctant *array;
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
          TMROctant p = array[i];
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
                TMROctant node;
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
                TMROctant *t = nodes->contains(&node, use_node_search);
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
                // of the finer two octants attached to the larger
                // octant.
                TMROctant node1, node2;
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
                TMROctant *t;
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

      // Get the octant array of nodes
      TMROctantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Extract the array of octants
      int size;
      TMROctant *array;
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

      // Now order the face. Get the octant array of nodes
      TMROctantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Extract the array of octants
      int size;
      TMROctant *array;
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
  queues = new TMROctantQueue*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      queues[k] = new TMROctantQueue();
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
        TMROctantArray *nodes;
        quadtrees[face]->getNodes(&nodes);

        // Get the octant arrays
        int size;
        TMROctant *array;
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
            TMROctant q;

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
              TMROctantArray *dest_nodes;
              quadtrees[dest_face]->getNodes(&dest_nodes);
              
              // Get the node and set the tag value
              const int use_nodes = 1;
              TMROctant *t = dest_nodes->contains(&q, use_nodes);
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
  arrays = new TMROctantArray*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      // Create the arrays
      arrays[k] = queues[k]->toArray();
      delete queues[k];

      // Get the array of octants
      int size;
      TMROctant *array;
      arrays[k]->getArray(&array, &size);

      // Set the array of octants to their destination
      MPI_Isend(array, size, TMROctant_MPI_type, 
                k, 0, comm, &send_requests[nsends]);
      nsends++;
    }
    else {
      arrays[k] = NULL;
    }
  }

  // Receive the arrays of incoming octants
  for ( int k = 0; k < mpi_size-1; k++ ){
    // Probe the recieved messages
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    // Retrieve the size and information for the incoming message
    int source = status.MPI_SOURCE;
    int tag = status.MPI_TAG;
    int size = 0;
    MPI_Get_count(&status, TMROctant_MPI_type, &size);

    // Allocate the incoming array
    TMROctant *array = new TMROctant[ size ];
    MPI_Recv(array, size, TMROctant_MPI_type,
             source, tag, comm, MPI_STATUS_IGNORE);

    // Push the octants into their corresponding trees
    for ( int i = 0; i < size; i++ ){
      int face = array[i].level;

      // Get the face nodes
      TMROctantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Search the nodes
      const int use_nodes = 1;
      TMROctant *t = nodes->contains(&array[i], use_nodes);
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
      // Get the nodes octants
      TMROctantArray *nodes;
      quadtrees[face]->getNodes(&nodes);

      // Get the node array
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);
      
      // Count up the number of dependent nodes
      for ( int i = 0; i < size; i++ ){
        if (array[i].tag < 0){
          array[i].tag = -(num_dep_nodes+1);
          num_dep_nodes++;
        }
      }
      
      // Order the elements for this octant
      TMROctantArray *elements;
      quadtrees[face]->getElements(&elements);
      elements->getArray(&array, &size);
      for ( int i = 0; i < size; i++ ){
        array[i].tag = num_elements;
        num_elements++;
      }
    }
  }
  */
  
  MPI_Barrier(comm);

  delete [] face_block_owners;
  delete [] edge_block_owners;
  delete [] node_block_owners;
}
