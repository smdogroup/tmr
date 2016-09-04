#include "TMROctForest.h"

// Include METIS
extern "C" {
#include "metis.h"
}

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

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
  Map from the face to the edge indices
*/
const int face_to_edge_index[][4] = {{4,6,8,10}, {5,7,9,11},
                                     {0,2,8,9}, {1,3,10,11},
                                     {0,1,4,5}, {2,3,6,7}};

/*
  All possible orientations for two connecting faces
*/
const int face_orientations[][4] = {{0,1,2,3},
                                    {2,0,3,1},
                                    {3,2,1,0},
                                    {1,3,0,2},
                                    {0,2,1,3},
                                    {2,3,0,1},
                                    {3,1,2,0},
                                    {1,0,3,2}};

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
  else if (face_id == 3){
    *u = y;
    *v = hmax - h - x;
  }
  else if (face_id == 4){
    *u = y;
    *v = x;
  }
  else if (face_id == 5){
    *u = x;
    *v = hmax - h - y;
  }
  else if (face_id == 6){
    *u = hmax - h - y;
    *v = hmax - h - x;
  }
  else {
    *u = hmax - h - x;
    *v = y;
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
  else if (face_id == 3){
    *x = hmax - h - v;
    *y = u; 
  }
  else if (face_id == 4){
    *x = v;
    *y = u;
  }
  else if (face_id == 5){
    *x = u;
    *y = hmax - h - v;
  }
  else if (face_id == 6){
    *x = hmax - h - v; 
    *y = hmax - h - u;
  }
  else {
    *x = hmax - h - u; 
    *y = v;
  }
}

/*
  Given the x/y locations on a face with orientation face_id,
  find/return the corresponding u/v coordinates on the transformed
  face.
*/
inline void get_face_node_coords( const int face_id,
                                  const int32_t x, const int32_t y,
                                  int32_t *u, int32_t *v ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;  
  if (face_id == 0){
    *u = x;
    *v = y;
  }
  else if (face_id == 1){
    *u = hmax - y;
    *v = x;
  }
  else if (face_id == 2){
    *u = hmax - x;
    *v = hmax - y;
  }
  else if (face_id == 3){
    *u = y;
    *v = hmax - x;
  }
  else if (face_id == 4){
    *u = y;
    *v = x;
  }
  else if (face_id == 5){
    *u = x;
    *v = hmax - y;
  }
  else if (face_id == 6){
    *u = hmax - y;
    *v = hmax - x;
  }
  else {
    *u = hmax - x;
    *v = y;
  }
}

/*
  Given the u/v locations on the root face, set the coordinates on the
  destination face.
*/
inline void set_face_node_coords( const int face_id,
                                  const int32_t u, const int32_t v,
                                  int32_t *x, int32_t *y ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;  
  if (face_id == 0){
    *x = u;
    *y = v;
  }
  else if (face_id == 1){
    *x = v;
    *y = hmax - u;
  }
  else if (face_id == 2){
    *x = hmax - u;
    *y = hmax - v; 
  }
  else if (face_id == 3){
    *x = hmax - v;
    *y = u; 
  }
  else if (face_id == 4){
    *x = v;
    *y = u;
  }
  else if (face_id == 5){
    *x = u;
    *y = hmax - v;
  }
  else if (face_id == 6){
    *x = hmax - v; 
    *y = hmax - u;
  }
  else {
    *x = hmax - u;
    *y = v;
  }
}

/*
  Compare integers for sorting
*/
static int compare_integers( const void *a, const void *b ){
  return (*(int*)a - *(int*)b);
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
  num_owned_blocks = 0;

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
  dep_faces = NULL;
  dep_edges = NULL;
  mpi_block_owners = NULL;
  owned_blocks = NULL;

  // Set the size of the mesh
  mesh_order = 0;

  // Set data for the number of elements/nodes/dependents
  node_range = NULL;
  num_elements = 0;
  num_dep_nodes = 0;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;
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
  if (dep_faces){
    for ( int i = 0; i < num_owned_blocks; i++ ){
      if (dep_faces[i]){ delete dep_faces[i]; }
    }
    delete [] dep_faces;
  }
  if (dep_edges){
    for ( int i = 0; i < num_owned_blocks; i++ ){
      if (dep_edges[i]){ delete dep_edges[i]; }
    }
    delete [] dep_edges;
  }
  if (mpi_block_owners){ delete [] mpi_block_owners; }
  if (owned_blocks){ delete [] owned_blocks; }
  if (node_range){ delete [] node_range; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
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
                                    int _num_blocks,
                                    int partition ){
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
  if (dep_faces){
    for ( int i = 0; i < num_owned_blocks; i++ ){
      if (dep_faces[i]){ delete dep_faces[i]; }
    }
    delete [] dep_faces;
  }
  if (dep_edges){
    for ( int i = 0; i < num_owned_blocks; i++ ){
      if (dep_edges[i]){ delete dep_edges[i]; }
    }
    delete [] dep_edges;
  }
  if (mpi_block_owners){ delete [] mpi_block_owners; }
  if (owned_blocks){ delete [] owned_blocks; }
  if (node_range){ delete [] node_range; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }
  node_range = NULL;
  dep_faces = NULL;
  dep_edges = NULL;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;

  // Copy over the data locally
  num_nodes = _num_nodes;
  num_edges = 0;
  num_faces = 0;
  num_blocks = _num_blocks;
  num_owned_blocks = 0;

  // Copy over the block connectivity
  block_conn = new int[ 8*num_blocks ];
  memcpy(block_conn, _block_conn, 8*num_blocks*sizeof(int));

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

  // Set up the partition using metis
  mpi_block_owners = new int[ num_blocks ];
  memset(mpi_block_owners, 0, num_blocks*sizeof(int));

  // Compute the partition on the root processor and
  // broadcast the result to all other processors
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);  
  MPI_Comm_size(comm, &mpi_size);
  if (mpi_rank == 0){
    if (partition){
      computePartition(mpi_size, NULL, mpi_block_owners);
    }
    else {
      for ( int block = 0; block < num_blocks; block++ ){
        mpi_block_owners[block] = block % mpi_size;
      }
    }
  }
  
  // Broadcast the face owners to all processors
  MPI_Bcast(mpi_block_owners, num_blocks, MPI_INT, 0, comm);

  // Determine the number of block owners
  num_owned_blocks = 0;
  for ( int i = 0; i < num_blocks; i++ ){
    if (mpi_rank == mpi_block_owners[i]){
      num_owned_blocks++;
    }
  }

  owned_blocks = new int[ num_owned_blocks ];
  for ( int i = 0, k = 0; i < num_blocks; i++ ){
    if (mpi_rank == mpi_block_owners[i]){
      owned_blocks[k] = i;
      k++;
    }
  }

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

  // Loop over all edges and determine their relative orientation
  for ( int edge = 0; edge < num_edges; edge++ ){
    int block_owner = num_blocks;
    int owner_index = 0;

    // Scan through the blocks pointing to this edge to determine
    // the block owner - the block with the lowest index
    for ( int ip = edge_block_ptr[edge]; ip < edge_block_ptr[edge+1]; ip++ ){
      int block = edge_block_conn[ip];
      if (block < block_owner){
        block_owner = block;

        // Find the new owner index
        owner_index = 0;
        for ( int j = 0; j < 12; j++, owner_index++ ){
          if (block_edge_conn[12*block + j] == edge){
            break;
          }
        }
      }
    }

    // Retrieve the first and second node numbers
    int n1 = block_conn[8*block_owner + block_to_edge_nodes[owner_index][0]];
    int n2 = block_conn[8*block_owner + block_to_edge_nodes[owner_index][1]];

    // Now determine the local edge index on each block and adjust
    // connectivity data
    for ( int ip = edge_block_ptr[edge]; ip < edge_block_ptr[edge+1]; ip++ ){
      // Find the block index
      int block = edge_block_conn[ip];
      
      for ( int edge_index = 0; edge_index < 12; edge_index++ ){
        int nn1 = block_conn[8*block + block_to_edge_nodes[edge_index][0]];
        int nn2 = block_conn[8*block + block_to_edge_nodes[edge_index][1]];
        
        // Check if the edges now match up
        if ((n1 == nn1 && n2 == nn2) ||
            (n1 == nn2 && n2 == nn1)){
          edge_block_conn[ip] = 12*block + edge_index;
          break;
        }
      }
    }
  }

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
        // Get the face nodes for the j-th face of the i-th block
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
        int node = face_nodes[0];
        for ( int ip = node_block_ptr[node];
              ip < node_block_ptr[node+1]; ip++ ){
          int ii = node_block_conn[ip];

          // Skip this if the blocks are the same
          if (ii == i){ continue; }

          // Loop over all the faces for the ii-th block
          int face_equiv = 0;
          for ( int jj = 0; jj < 6; jj++ ){
            // Get the nodes corresponding to this face
            int adj_face_nodes[4];
            for ( int k = 0; k < 4; k++ ){
              adj_face_nodes[k] = 
                block_conn[8*ii + block_to_face_nodes[jj][k]];
            }

            // Loop over the relative orientations between this block
            // and the last block
            for ( int ort = 0; ort < 8; ort++ ){
              face_equiv = 
                (face_nodes[0] == adj_face_nodes[face_orientations[ort][0]] &&
                 face_nodes[1] == adj_face_nodes[face_orientations[ort][1]] &&
                 face_nodes[2] == adj_face_nodes[face_orientations[ort][2]] &&
                 face_nodes[3] == adj_face_nodes[face_orientations[ort][3]]);
              // We've found a matching face
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
  memset(block_face_ids, 0, 6*num_blocks*sizeof(int));

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
    for ( int ip = face_block_ptr[face]; ip < face_block_ptr[face+1]; ip++ ){
      // Find the block index
      int block = face_block_conn[ip];
      
      // Find the adjacent face and keep the corresponding face index
      int face_index = 0;
      for ( ; face_index < 6; face_index++ ){
        if (face == block_face_conn[6*block + face_index]){
          // Now, find the relative orientations between the two faces.
          // First get the nodes corresponding to this face
          int adj_face_nodes[4];
          for ( int k = 0; k < 4; k++ ){
            adj_face_nodes[k] = 
              block_conn[8*block + block_to_face_nodes[face_index][k]];
          }
          
          // Loop over the relative orientations between this block
          // and the last block
          int face_equiv = 0;
          for ( int ort = 0; ort < 8; ort++ ){
            face_equiv = 
              (owner_nodes[0] == adj_face_nodes[face_orientations[ort][0]] &&
               owner_nodes[1] == adj_face_nodes[face_orientations[ort][1]] &&
               owner_nodes[2] == adj_face_nodes[face_orientations[ort][2]] &&
               owner_nodes[3] == adj_face_nodes[face_orientations[ort][3]]);
            // Set the orientation and break out of this loop
            if (face_equiv){
              block_face_ids[6*block + face_index] = ort;
              break;
            }
          }

          // If we have found the face index, then break from the
          // outer loop
          if (face_equiv){
            break;
          }
        }
      }

      // Now, readjust the connectivity to also store the face index
      // of the local face on the adjacent block
      face_block_conn[ip] = 6*block + face_index;
    }
  }
}

/*
  Get the array of octrees 
*/
int TMROctForest::getOctrees( TMROctree ***_trees ){
  if (_trees){ *_trees = octrees; }
  return num_blocks;
}

/*
  Get mesh/ownership information
*/
int TMROctForest::getOwnedOctrees( const int **_owned_blocks ){
  if (_owned_blocks){ 
    *_owned_blocks = owned_blocks; 
  }
  return num_owned_blocks;
}

/*
  Retrieve information about the connectivity between 
  blocks, faces, edges and nodes
*/
void TMROctForest::getConnectivity( int *_nblocks, int *_nfaces, 
                                    int *_nedges, int *_nnodes, 
                                    const int **_block_conn, 
                                    const int **_block_face_conn, 
                                    const int **_block_edge_conn,
                                    const int **_block_face_ids ){
  if (_nblocks){ *_nblocks = num_blocks; }
  if (_nfaces){ *_nfaces = num_faces; }
  if (_nedges){ *_nedges = num_edges; }
  if (_nnodes){ *_nnodes = num_nodes; }
  if (_block_conn){ *_block_conn = block_conn; }
  if (_block_face_conn){ *_block_face_conn = block_face_conn; }
  if (_block_edge_conn){ *_block_edge_conn = block_edge_conn; }
  if (_block_face_ids){ *_block_face_ids = block_face_ids; }
}

/*
  Retrieve the inverse of the connectivity
*/
void TMROctForest::getInverseConnectivity( const int **_node_block_conn,
                                           const int **_node_block_ptr,
                                           const int **_edge_block_conn,
                                           const int **_edge_block_ptr,
                                           const int **_face_block_conn,
                                           const int **_face_block_ptr ){
  if (_node_block_conn){ *_node_block_conn = node_block_conn; }
  if (_node_block_ptr){ *_node_block_ptr = node_block_ptr; }
  if (_edge_block_conn){ *_edge_block_conn = edge_block_conn; }
  if (_edge_block_ptr){ *_edge_block_ptr = edge_block_ptr; }
  if (_face_block_conn){ *_face_block_conn = face_block_conn; }
  if (_face_block_ptr){ *_face_block_ptr = face_block_ptr; }
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
  octrees = new TMROctree*[ num_blocks ];
  memset(octrees, 0, num_blocks*sizeof(TMROctree*));
  for ( int i = 0; i < num_owned_blocks; i++ ){
    int block = owned_blocks[i];
    octrees[block] = new TMROctree(refine_level);
  }
}

/*
  Create a forest with the specified refinement levels
*/
void TMROctForest::createTrees( int refine_levels[] ){
  if (octrees){ 
    for ( int i = 0; i < num_blocks; i++ ){
      if (octrees[i]){ delete octrees[i]; }
    }
    delete [] octrees;
  }

  // Create the octrees
  octrees = new TMROctree*[ num_blocks ];
  memset(octrees, 0, num_blocks*sizeof(TMROctree*));
  for ( int i = 0; i < num_owned_blocks; i++ ){
    int block = owned_blocks[i];
    octrees[block] = new TMROctree(refine_levels[block]);
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
  for ( int i = 0; i < num_owned_blocks; i++ ){
    int block = owned_blocks[i];
    octrees[block] = new TMROctree(nrand, min_level, max_level);
  }
}

/*
  Repartition the mesh based on the number of elements per block
*/
void TMROctForest::repartition(){
  // Get the communicator rank/size
  int mpi_size, mpi_rank;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // You can't repartition that which cannot be repartitioned
  if (mpi_size <= 1){
    return;
  }

  // First, this stores the number of elements on octrees owned on
  // each processor
  int *elem_counts = new int[ num_owned_blocks ];
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];
    elem_counts[owned] = octrees[block]->getNumElements();
  }

  // Allocate space for the new partition
  int *new_part = new int[ num_blocks ];

  // Gather the element counts to the root processor
  if (mpi_rank != 0){
    MPI_Gatherv(elem_counts, num_owned_blocks, MPI_INT,
                NULL, NULL, NULL, MPI_INT, 0, comm);
  }
  else {
    int *all_elem_counts = new int[ num_blocks ];
    int *recv_counts = new int[ mpi_size ];
    int *recv_displ = new int[ mpi_size ];

    // Count up the recvs from each processor
    memset(recv_counts, 0, mpi_size*sizeof(int));
    for ( int block = 0; block < num_blocks; block++ ){
      recv_counts[mpi_block_owners[block]]++;
    }
    
    // Compute the displacement offsets
    recv_displ[0] = 0;
    for ( int k = 1; k < mpi_size; k++ ){
      recv_displ[k] = recv_displ[k-1] + recv_counts[k-1];
    }
    
    // Receive all the elements
    MPI_Gatherv(elem_counts, num_owned_blocks, MPI_INT,
                all_elem_counts, recv_counts, recv_displ, MPI_INT,
                0, comm);

    // Fill in the number of elements per processor back in the
    // original order from the original blocks
    int *nelems = new int[ num_blocks ];
    for ( int block = 0; block < num_blocks; block++ ){
      int mpi_owner = mpi_block_owners[block];
      
      nelems[block] = all_elem_counts[recv_displ[mpi_owner]];
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
  MPI_Bcast(new_part, num_blocks, MPI_INT, 0, comm);

  // Now create the octrees array
  TMROctree **new_octrees = new TMROctree*[ num_blocks ];
  memset(new_octrees, 0, num_blocks*sizeof(TMROctree*));

  // Only redistribute the elements, not the mesh
  int send_count = 0;
  MPI_Request *requests = new MPI_Request[ num_owned_blocks ];
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];
    int dest = new_part[block];

    if (dest != mpi_rank){
      // Get the element array
      TMROctantArray *elements;
      octrees[block]->getElements(&elements);
      
      // Get the actual octant array
      int size;
      TMROctant *array;
      elements->getArray(&array, &size);

      // Send the element array to the new owner
      MPI_Isend(array, size, TMROctant_MPI_type,
                dest, block, comm, &requests[send_count]);
      send_count++;
    }
    else {
      new_octrees[block] = octrees[block];
      octrees[block] = NULL;
    }
  }

  // Determine the new number of block owners based on the new
  // partition and the number of expected incoming messages (octrees)
  // from all processors
  int num_owned = 0;
  for ( int block = 0; block < num_blocks; block++ ){
    if (mpi_rank == new_part[block]){
      num_owned++;
    }
  }

  int incoming_blocks = 0;
  int *new_owned = new int[ num_owned ];
  for ( int block = 0, k = 0; block < num_blocks; block++ ){
    if (mpi_rank == new_part[block]){
      new_owned[k] = block;
      k++;

      // Check if this block will need to be sent
      if (mpi_block_owners[block] != mpi_rank){
        incoming_blocks++;
      }
    }
  }
  
  // Loop over the old senders
  for ( int k = 0; k < incoming_blocks; k++ ){
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

    // Get the source, tag == block, and size of the message
    int source = status.MPI_SOURCE;
    int block = status.MPI_TAG;
    int size = 0;
    MPI_Get_count(&status, TMROctant_MPI_type, &size);

    // Allocate the incoming array and receive it
    TMROctant *array = new TMROctant[ size ];
    MPI_Recv(array, size, TMROctant_MPI_type,
             source, block, comm, MPI_STATUS_IGNORE);
    
    // Create the new local octree
    TMROctantArray *list = new TMROctantArray(array, size);
    new_octrees[block] = new TMROctree(list);
  }

  // Wait for all the sends to complete
  MPI_Waitall(send_count, requests, MPI_STATUSES_IGNORE);

  // Free the old octrees that were exchanged to another processor and
  // assign the new octrees array
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];
    if (octrees[block]){ 
      delete octrees[block]; 
    }
  }
  delete [] octrees;
  octrees = new_octrees;

  // Reset local information associated with the nodal ordering
  num_elements = 0;
  num_dep_nodes = 0;
  if (dep_faces){
    for ( int owned = 0; owned < num_owned_blocks; owned++ ){
      if (dep_faces[owned]){ delete dep_faces[owned]; }
    }
    delete [] dep_faces;
  }
  if (dep_edges){
    for ( int i = 0; i < num_owned_blocks; i++ ){
      if (dep_edges[i]){ delete dep_edges[i]; }
    }
    delete [] dep_edges;
  }
  dep_faces = NULL;
  dep_edges = NULL;

  // Assign the new partition and the new owner list
  delete [] owned_blocks;
  owned_blocks = new_owned;
  num_owned_blocks = num_owned;

  delete [] mpi_block_owners;
  mpi_block_owners = new_part;
}

/*
  Partition the mesh based on the super mesh connectivity and
  optionally vertex weights (typically the element count per block)

  input:
  num_part:     the number of partitions
  vwgts:        the vertex weights

  output:
  part:         the block-processor assignment
*/
void TMROctForest::computePartition( int part_size, int *vwgts,
                                     int *part ){
  // Set the pointer to the block connectivity
  int *ptr = new int[ num_blocks+1 ];
  memset(ptr, 0, (num_blocks+1)*sizeof(int));

  // Set the default size for the connectivity
  int max_size = 9*num_blocks;
  int *conn = new int[ max_size ];
    
  // Allocate space to store flags to indicate whether the block has
  // been added to the array
  int *cols = new int[ num_blocks ];
  memset(cols, 0, num_blocks*sizeof(int));

  // The adjacency node count
  int min_adj_count = 4;
  
  for ( int block = 0; block < num_blocks; block++ ){
    // Set the new pointer
    ptr[block+1] = ptr[block];
    
    // Loop over all the blocks
    for ( int j = 0; j < 8; j++ ){
      int node = block_conn[8*block + j];
      
      // Loop over all of the nodes
      for ( int ip = node_block_ptr[node]; 
            ip < node_block_ptr[node+1]; ip++ ){
        int adj_block = node_block_conn[ip];

        if (cols[adj_block]+1 == min_adj_count*(block+1)){
          // Extend the array
          if (ptr[block+1] >= max_size){
            max_size *= 2;
            int *tmp = new int[ max_size ];
            memcpy(tmp, conn, ptr[block+1]*sizeof(int));
            delete [] conn;
            conn = tmp;
          }

          // Set the new element into the connectivity
          conn[ptr[block+1]] = adj_block;
          ptr[block+1]++;

          // Set the flag so that its not added again
          cols[adj_block]++;
        }
        else if (cols[adj_block] <= min_adj_count*block){
          // This is the first time this block has been encountered on
          // this loop. Add it and increment the pointer
          cols[adj_block] = min_adj_count*block+1;
        }
        else {
          // This block has been encountered before
          cols[adj_block]++;
        }
      }
    }
    
    // Sort the array
    int len = ptr[block+1] - ptr[block];
    qsort(&conn[ptr[block]], len, sizeof(int), compare_integers);
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
    METIS_PartGraphRecursive(&num_blocks, ptr, conn, vwgts, adjwgts,
                             &wgtflag, &numflag, &part_size,
                             options, &edgecut, part);
  }
  else {
    METIS_PartGraphKway(&num_blocks, ptr, conn, vwgts, adjwgts, 
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
TMROctForest *TMROctForest::duplicate(){
  TMROctForest *dup = new TMROctForest(comm);
  if (octrees){
    // Copy over the connectivity data 
    dup->num_nodes = num_nodes;
    dup->num_edges = num_edges;
    dup->num_faces = num_faces;
    dup->num_blocks = num_blocks;
    dup->num_owned_blocks = num_owned_blocks;

    // Allocate/copy the block connectivities
    dup->block_conn = new int[ 8*num_blocks ];
    dup->block_face_conn = new int[ 6*num_faces ];
    dup->block_face_ids = new int[ 6*num_faces ];
    dup->block_edge_conn = new int[ 12*num_faces ];
    memcpy(dup->block_conn, block_conn, 8*num_blocks*sizeof(int));
    memcpy(dup->block_face_conn, block_face_conn, 
           6*num_blocks*sizeof(int));
    memcpy(dup->block_face_ids, block_face_ids, 
           6*num_blocks*sizeof(int));
    memcpy(dup->block_edge_conn, block_edge_conn, 
           12*num_blocks*sizeof(int));
    
    // Allocate/copy the inverse relationships
    dup->node_block_ptr = new int[ num_nodes+1 ];
    dup->node_block_conn = new int[ node_block_ptr[num_nodes] ];
    memcpy(dup->node_block_ptr, node_block_ptr, 
           (num_nodes+1)*sizeof(int));
    memcpy(dup->node_block_conn, node_block_conn, 
           node_block_ptr[num_nodes]*sizeof(int));

    dup->edge_block_ptr = new int[ num_edges+1 ];
    dup->edge_block_conn = new int[ edge_block_ptr[num_edges] ];
    memcpy(dup->edge_block_ptr, edge_block_ptr, 
           (num_edges+1)*sizeof(int));
    memcpy(dup->edge_block_conn, edge_block_conn, 
           edge_block_ptr[num_edges]*sizeof(int));

    dup->face_block_ptr = new int[ num_faces+1 ];
    dup->face_block_conn = new int[ face_block_ptr[num_faces] ];
    memcpy(dup->face_block_ptr, face_block_ptr, 
           (num_faces+1)*sizeof(int));
    memcpy(dup->face_block_conn, face_block_conn, 
           face_block_ptr[num_faces]*sizeof(int));

    // Allocate/copy the block ownership
    dup->mpi_block_owners = new int[ num_blocks ];
    memcpy(dup->mpi_block_owners, mpi_block_owners, 
           num_blocks*sizeof(int));

    dup->owned_blocks = new int[ num_owned_blocks ];
    memcpy(dup->owned_blocks, owned_blocks, 
           num_owned_blocks*sizeof(int));

    // Duplicate all the octrees
    dup->octrees = new TMROctree*[ num_blocks ];
    memset(dup->octrees, 0, num_blocks*sizeof(TMROctree*));
    for ( int i = 0; i < num_owned_blocks; i++ ){
      int block = owned_blocks[i];
      TMROctantArray *elements;
      octrees[block]->getElements(&elements);
      dup->octrees[block] = new TMROctree(elements->duplicate());
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
TMROctForest *TMROctForest::coarsen(){
  TMROctForest *coarse = new TMROctForest(comm);
  if (octrees){
    // Copy over the connectivity data 
    coarse->num_nodes = num_nodes;
    coarse->num_edges = num_edges;
    coarse->num_faces = num_faces;
    coarse->num_blocks = num_blocks;
    coarse->num_owned_blocks = num_owned_blocks;

    // Allocate/copy the block connectivities
    coarse->block_conn = new int[ 8*num_blocks ];
    coarse->block_face_conn = new int[ 6*num_faces ];
    coarse->block_face_ids = new int[ 6*num_faces ];
    coarse->block_edge_conn = new int[ 12*num_faces ];
    memcpy(coarse->block_conn, block_conn, 8*num_blocks*sizeof(int));
    memcpy(coarse->block_face_conn, block_face_conn, 
           6*num_blocks*sizeof(int));
    memcpy(coarse->block_face_ids, block_face_ids, 
           6*num_blocks*sizeof(int));
    memcpy(coarse->block_edge_conn, block_edge_conn, 
           12*num_blocks*sizeof(int));
    
    // Allocate/copy the inverse relationships
    coarse->node_block_ptr = new int[ num_nodes+1 ];
    coarse->node_block_conn = new int[ node_block_ptr[num_nodes] ];
    memcpy(coarse->node_block_ptr, node_block_ptr, 
           (num_nodes+1)*sizeof(int));
    memcpy(coarse->node_block_conn, node_block_conn, 
           node_block_ptr[num_nodes]*sizeof(int));

    coarse->edge_block_ptr = new int[ num_edges+1 ];
    coarse->edge_block_conn = new int[ edge_block_ptr[num_edges] ];
    memcpy(coarse->edge_block_ptr, edge_block_ptr, 
           (num_edges+1)*sizeof(int));
    memcpy(coarse->edge_block_conn, edge_block_conn, 
           edge_block_ptr[num_edges]*sizeof(int));

    coarse->face_block_ptr = new int[ num_faces+1 ];
    coarse->face_block_conn = new int[ face_block_ptr[num_faces] ];
    memcpy(coarse->face_block_ptr, face_block_ptr, 
           (num_faces+1)*sizeof(int));
    memcpy(coarse->face_block_conn, face_block_conn, 
           face_block_ptr[num_faces]*sizeof(int));

    // Allocate/copy the block ownership
    coarse->mpi_block_owners = new int[ num_blocks ];
    memcpy(coarse->mpi_block_owners, mpi_block_owners, 
           num_blocks*sizeof(int));

    coarse->owned_blocks = new int[ num_owned_blocks ];
    memcpy(coarse->owned_blocks, owned_blocks,
           num_owned_blocks*sizeof(int));

    // Coarselicate all the octrees
    coarse->octrees = new TMROctree*[ num_blocks ];
    memset(coarse->octrees, 0, num_blocks*sizeof(TMROctree*));
    for ( int i = 0; i < num_owned_blocks; i++ ){
      int block = owned_blocks[i];
      coarse->octrees[block] = octrees[block]->coarsen();
    }
  }

  return coarse;
}

/*
  Add the face neighbors for an adjacent tree
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
  for ( int ip = face_block_ptr[face]; ip < face_block_ptr[face+1]; ip++ ){
    // Get the adjacent block and the face index on the adjacent block
    int adjacent = face_block_conn[ip]/6;
    if (adjacent != block){
      if (!hash[adjacent]){
        hash[adjacent] = new TMROctantHash();
        queue[adjacent] = new TMROctantQueue();
      }

      // Get the face index on the adjacent block
      int adj_index = face_block_conn[ip] % 6;

      // Get the face id for the matching face
      face_id = block_face_ids[6*adjacent + adj_index];

      // Get the neighboring octant on the face
      TMROctant neighbor;
      neighbor.level = p.level;
      if (adj_index < 2){
        neighbor.x = (hmax - 2*h)*(adj_index % 2);
        set_face_oct_coords(face_id, 2*h, u, v, 
                            &neighbor.y, &neighbor.z);
      }
      else if (adj_index < 4){
        neighbor.y = (hmax - 2*h)*(adj_index % 2);
        set_face_oct_coords(face_id, 2*h, u, v, 
                            &neighbor.x, &neighbor.z);
      }
      else {
        neighbor.z = (hmax - 2*h)*(adj_index % 2);
        set_face_oct_coords(face_id, 2*h, u, v, 
                            &neighbor.x, &neighbor.y);
      }
      
      // Add the octant to the list
      if (hash[adjacent]->addOctant(&neighbor)){
        queue[adjacent]->push(&neighbor);
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
  for ( int ip = edge_block_ptr[edge]; ip < edge_block_ptr[edge+1]; ip++ ){
    // Get the block that is adjacent across this edge
    int adjacent = edge_block_conn[ip]/12;
    if (adjacent != block){
      if (!hash[adjacent]){
        hash[adjacent] = new TMROctantHash();
        queue[adjacent] = new TMROctantQueue();
      }

      // Get the adjacent edge index
      int adj_index = edge_block_conn[ip] % 12;

      // Get the nodes on the adjacent block
      int nn1 = block_conn[8*adjacent + block_to_edge_nodes[adj_index][0]];
      int nn2 = block_conn[8*adjacent + block_to_edge_nodes[adj_index][1]];

      // Add the octant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - 2*h - ucoord;
      }
      
      TMROctant neighbor;
      neighbor.level = p.level;
      if (adj_index < 4){
        neighbor.x = u;
        neighbor.y = (hmax - 2*h)*(adj_index % 2);
        neighbor.z = (hmax - 2*h)*(adj_index/2);
      }
      else if (adj_index < 8){
        neighbor.x = (hmax - 2*h)*(adj_index % 2);
        neighbor.y = u;
        neighbor.z = (hmax - 2*h)*((adj_index-4)/2);
      }
      else {
        neighbor.x = (hmax - 2*h)*(adj_index % 2);
        neighbor.y = (hmax - 2*h)*((adj_index-8)/2);
        neighbor.z = u;
      }
      
      // Add the octant to the list
      if (hash[adjacent]->addOctant(&neighbor)){
        queue[adjacent]->push(&neighbor);
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
          neighbor.x = (hmax - 2*h)*(j % 2);
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
              edge = (q.y < 0 ? 0 : 1) + (q.z < 0 ? 0 : 2);
            }
            else if (ex && ez){
              edge = (q.x < 0 ? 4 : 5) + (q.z < 0 ? 0 : 2);
            }
            else { // (ex && ey)
              edge = (q.x < 0 ? 8 : 9) + (q.y < 0 ? 0 : 2);
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
  Balance the forest of octrees

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

  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

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
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

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
    
    // Set the elements into the octree
    TMROctantArray *elements = hash[block]->toArray();
    elements->sort();
    octrees[block]->setElements(elements);
    
    // Free the corresponding hash
    delete hash[block];
  }

  delete [] hash;
}

/*
  Add the octant to the processor queues that correspond to 
  the non-local blocks that touch the corner
*/
void TMROctForest::addCornerOctantToQueues( const int node,
                                            const int mpi_rank,
                                            TMROctant *q,
                                            TMROctantQueue **queues ){
  for ( int ip = node_block_ptr[node]; ip < node_block_ptr[node+1]; ip++ ){
    int block = node_block_conn[ip];
    int rank = mpi_block_owners[block];
    if (rank != mpi_rank){
      queues[rank]->push(q);
    }
  }
}

/*
  Add the octant to the processor queues corresponding to the
  non-local blocks that touch the given edge
*/
void TMROctForest::addEdgeOctantToQueues( const int edge,
                                          const int mpi_rank,
                                          TMROctant *q,
                                          TMROctantQueue **queues ){
  for ( int ip = edge_block_ptr[edge]; ip < edge_block_ptr[edge+1]; ip++ ){
    int block = edge_block_conn[ip]/12;
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
    int block = face_block_conn[ip]/6;
    int rank = mpi_block_owners[block];
    if (rank != mpi_rank){
      queues[rank]->push(q);
    }
  }
}

/*
  The following code exchanges the neighboring octants for each
  locally owned octree within the forest.

  This code exchanges non-local octants across processors so that we
  can locally query octants on adjacent octrees without having to
  perform parallel communication.

  Note that this code creates partial non-local octrees that are
  adjacent to the local octrees in the forest. These partial local
  octrees should be freed after the nodal ordering has been computed.
*/
void TMROctForest::recvOctNeighbors(){
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

  // For each block, determine the edge, face and corner neighbors 
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];
    
    // Flag to indicate whether any edge/face is non-local
    int has_non_local = 0;
    
    // Check if any of the edges has to be sent to another processor
    for ( int k = 0; k < 12; k++ ){
      int edge = block_edge_conn[12*block + k];
      for ( int ip = edge_block_ptr[edge];
            ip < edge_block_ptr[edge+1]; ip++ ){
        int dest_block = edge_block_conn[ip]/12;
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
          int dest_block = face_block_conn[ip]/6;
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
        
      // If this is a single octant, it fills the entire volume
      if (size == 1){
        TMROctant q = array[0];
        q.tag = block;

        // Send the octant to all neighbors
        for ( int node_index = 0; node_index < 8; node_index++ ){
          int node = block_conn[8*block + node_index];
          addCornerOctantToQueues(node, mpi_rank, &q, queues);
        }

        // Send the octant to all adjacent edges
        for ( int edge_index = 0; edge_index < 12; edge_index++ ){
          int edge = block_edge_conn[12*block + edge_index];
          addEdgeOctantToQueues(edge, mpi_rank, &q, queues);
        }

        // Send the octant to all adjacent faces
        for ( int face_index = 0; face_index < 6; face_index++ ){
          int face = block_face_conn[6*block + face_index];
          addFaceOctantToQueues(face, mpi_rank, &q, queues);
        }
      }
      else {
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
            if (fx && fy && fz){
              int node_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2) + (fz0 ? 0 : 4);
              int node = block_conn[8*block + node_index];
              addCornerOctantToQueues(node, mpi_rank, &q, queues);
            }
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
  Determine if there is an adjacent octant on the connecting face.

  Return true if an adjacent face is found across a block-face and
  false if no octant is found.
  
  input:
  face:          the face number
  face_index:    the local face index
  block_owner:   the index of the block owner
  b:             the quadrant
*/
int TMROctForest::checkAdjacentDepFaces( int face,
                                         int face_index,
                                         int block_owner,
                                         TMROctant *b ){
  // Get the side length of the octant
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);

  // Get the face id number
  int face_id = block_face_ids[6*block_owner + face_index];

  // Get the u/v coordinates for this node on the owner face
  int32_t u, v;
  if (face_index < 2){ // x-face
    get_face_oct_coords(face_id, h, b->y, b->z, &u, &v);
  }
  else if (face_index < 4){ // y-face
    get_face_oct_coords(face_id, h, b->x, b->z, &u, &v);
  }
  else { // z-face
    get_face_oct_coords(face_id, h, b->x, b->y, &u, &v);
  }

  // Loop over all the adjacent faces/blocks 
  for ( int ip = face_block_ptr[face]; ip < face_block_ptr[face+1]; ip++ ){
    int block = face_block_conn[ip]/6;

    if (block_owner != block){
      int adj_index = face_block_conn[ip] % 6;

      // Get the face_id corresponding to the orientation of this
      // adjacent face
      face_id = block_face_ids[6*block + adj_index];
      
      // Transform the octant p to the local octant coordinates
      TMROctant oct;
      oct.level = b->level;
      if (adj_index < 2){
        oct.x = (hmax - h)*(adj_index % 2);
        set_face_oct_coords(face_id, h, u, v, &oct.y, &oct.z);
      }
      else if (adj_index < 4){
        oct.y = (hmax - h)*(adj_index % 2);
        set_face_oct_coords(face_id, h, u, v, &oct.x, &oct.z);
      }
      else {
        oct.z = (hmax - h)*(adj_index % 2);
        set_face_oct_coords(face_id, h, u, v, &oct.x, &oct.y);
      }
      
      // Get the nodes from the array
      TMROctantArray *elements;
      octrees[block]->getElements(&elements);
      
      // If the more-refined element exists then label the
      // corresponding nodes as dependent
      const int use_nodes = 0;
      if (elements->contains(&oct, use_nodes)){
        return 1;
      }
    }
  }

  return 0;
}

/*
  Determine if there is an adjacent octant on the connecting edge.

  Return true if an adjacent edge is found across a block-edge and
  false if no octant is found.
  
  input:
  edge:          the edge number
  edge_index:    the local edge index
  block_owner:   the index of the block owner
  b:             the quadrant
*/
int TMROctForest::checkAdjacentDepEdges( int edge,
                                         int edge_index,
                                         int block_owner,
                                         TMROctant *b ){
  // Get the side length of the octant
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);

  // Store the u coordinate along the edge
  int32_t ucoord = 0;
  if (edge_index < 4){
    ucoord = b->x;
  }
  else if (edge_index < 8){
    ucoord = b->y;
  }
  else {
    ucoord = b->z;
  }

  // Retrieve the first and second node numbers
  int n1 = block_conn[8*block_owner + block_to_edge_nodes[edge_index][0]];
  int n2 = block_conn[8*block_owner + block_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent edges
  for ( int ip = edge_block_ptr[edge]; ip < edge_block_ptr[edge+1]; ip++ ){
    int block = edge_block_conn[ip]/12;

    if (block_owner != block){
      // Get the adjacent edge index
      int adj_index = edge_block_conn[ip] % 12;

      // Get the nodes on the adjacent block
      int nn1 = block_conn[8*block + block_to_edge_nodes[adj_index][0]];
      int nn2 = block_conn[8*block + block_to_edge_nodes[adj_index][1]];

      // Add the octant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - h - ucoord;
      }

      // Search for the neighboring octant
      TMROctant oct;
      oct.level = b->level;
      if (adj_index < 4){
        oct.x = u;
        oct.y = (hmax - h)*(adj_index % 2);
        oct.z = (hmax - h)*(adj_index/2);
      }
      else if (adj_index < 8){
        oct.x = (hmax - h)*(adj_index % 2);
        oct.y = u;
        oct.z = (hmax - h)*((adj_index-4)/2);
      }
      else {
        oct.x = (hmax - h)*(adj_index % 2);
        oct.y = (hmax - h)*((adj_index-8)/2);
        oct.z = u;
      }
      
      // Get the elements from the array
      TMROctantArray *elements;
      octrees[block]->getElements(&elements);
      
      // If the adjacent element exists then label the
      // corresponding nodes as dependent
      const int use_nodes = 0;
      if (elements->contains(&oct, use_nodes)){
        return 1;
      }
    }
  }

  return 0;
}

/*
  Compute the dependent nodes (hanging edge/face nodes) on each block
  and on the interfaces between adjacent blocks.

  The hanging face nodes may occur on any face within the block and
  are associated with the 4 parallel hanging edge nodes.  Within this
  code, we only store the corresponding face node and use the
  associated edges. Edge nodes along block interfaces may be hanging
  even if there is no associated hanging face node.

  side effects:
  dep_faces:   a list of the dependent faces
  dep_edges:   a list of dependent edges (aligned with block edges)
*/
void TMROctForest::computeDepFacesAndEdges(){
  if (dep_edges){ 
    for ( int i = 0; i < num_owned_blocks; i++ ){
      if (dep_edges[i]){ delete dep_edges[i]; }
    }
    delete [] dep_edges;
    dep_edges = NULL;
  }
  if (dep_faces){ 
    for ( int i = 0; i < num_owned_blocks; i++ ){
      if (dep_faces[i]){ delete dep_faces[i]; }
    }
    delete [] dep_faces;
    dep_faces = NULL;
  }
  if (dep_ptr){ delete [] dep_ptr;  dep_ptr = NULL; }
  if (dep_conn){ delete [] dep_conn;  dep_conn = NULL; }
  if (dep_weights){ delete [] dep_weights;  dep_weights = NULL; }

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Allocate the list of dependent nodes/faces
  dep_faces = new TMROctantArray*[ num_owned_blocks ];
  dep_edges = new TMROctantArray*[ num_owned_blocks ];
  memset(dep_faces, 0, num_owned_blocks*sizeof(TMROctantArray*));
  memset(dep_edges, 0, num_owned_blocks*sizeof(TMROctantArray*));

  // Determine and label the dependent nodes on each processor 
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

    // Create the queue that will store the dependent
    // faces/edges
    TMROctantQueue *dfaces = new TMROctantQueue();
    TMROctantQueue *dedges = new TMROctantQueue();
    
    // Get the octant elements
    TMROctantArray *elements;
    octrees[block]->getElements(&elements);
            
    // Get the elements themselves
    int size;
    TMROctant *array;
    elements->getArray(&array, &size);
    
    for ( int i = 0; i < size; i++ ){
      // Get the child id - skip any elements that are not
      // either child-0 or child-7
      int child_id = array[i].childId();

      // Get the side length of the element
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      const int32_t hmax = 1 << TMR_MAX_LEVEL;

      if (array[i].level > 0 && 
          (child_id == 0 || child_id == 7)){        
        // Check the adjacent elements along each face
        int face_index = 0;
        TMROctant p = array[i];
        if (child_id == 7){
          face_index = 1;
          p.getSibling(0, &p);
        }
        p.level = p.level - 1;

        for ( ; face_index < 6; face_index += 2 ){
          // Look for the face neighbor that is at the next level of
          // refinement from the current edge
          TMROctant q;
          p.faceNeighbor(face_index, &q);

          // Check if this element is across a face
          int fx = (q.x < 0 || q.x >= hmax);
          int fy = (q.y < 0 || q.y >= hmax);
          int fz = (q.z < 0 || q.z >= hmax);

          // If an adjacent octant of twice the size is on a
          // neighboring tree, or on the same tree, add it as
          // a dependent face
          int add_me = 0;
          if (fx || fy || fz){
            int face = block_face_conn[6*block + face_index];
            if (checkAdjacentDepFaces(face, face_index, block, &q)){
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
            p.tag = face_index;
            dfaces->push(&p);
          }
        }
      }

      // Determine which edges the octant touches if any
      int fx0 = (array[i].x == 0);
      int fy0 = (array[i].y == 0);
      int fz0 = (array[i].z == 0);
      int fx = (fx0 || array[i].x + h == hmax);
      int fy = (fy0 || array[i].y + h == hmax);
      int fz = (fz0 || array[i].z + h == hmax);

      if ((fx && fy) || (fy && fz) || (fx && fz)){
        // Convert the octant to the parent and check all 
        // relevant edges - note we may add more than one
        // dependent edge for each octant
        TMROctant p = array[i];
        p.getSibling(0, &p);
        p.level = p.level - 1;

        if (fy && fz){
          int edge_index = (fy0 ? 0 : 1) + (fz0 ? 0 : 2);
          int edge = block_edge_conn[12*block + edge_index];
          if (checkAdjacentDepEdges(edge, edge_index, block, &p)){
            p.tag = edge_index;
            dedges->push(&p);
          }
        }
        if (fx && fz){
          int edge_index = (fx0 ? 4 : 5) + (fz0 ? 0 : 2);
          int edge = block_edge_conn[12*block + edge_index];
          if (checkAdjacentDepEdges(edge, edge_index, block, &p)){
            p.tag = edge_index;
            dedges->push(&p);
          }
        }
        if (fx && fy){
          int edge_index = (fx0 ? 8 : 9) + (fy0 ? 0 : 2);
          int edge = block_edge_conn[12*block + edge_index];
          if (checkAdjacentDepEdges(edge, edge_index, block, &p)){
            p.tag = edge_index;
            dedges->push(&p);
          }
        }
      }
    }

    // Create the arrays of the dependent nodes/faces
    dep_faces[owned] = dfaces->toArray();
    dep_edges[owned] = dedges->toArray();

    // Delete the queues
    delete dfaces;
    delete dedges;
  }
}

/*
  Label the dependent face and edge nodes

  This code is called after all the dependent faces have been
  computed.  Note that this relies on the mesh being edge-balanced
  (which is required). 
*/
void TMROctForest::labelDependentNodes(){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

    // Get the nodes associated with this octree
    TMROctantArray *nodes;
    octrees[block]->getNodes(&nodes);

    // Get the array of dependent faces
    int dep_size;
    TMROctant *dep_array;
    dep_faces[owned]->getArray(&dep_array, &dep_size);

    for ( int i = 0; i < dep_size; i++ ){
      TMROctant *b = &dep_array[i];

      // Find the edge lengths of the element octant and
      // the node spacing on the dependent face
      const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);
      const int32_t hc = 1 << (TMR_MAX_LEVEL - b->level - (mesh_order-1));
          
      // Get the face index
      int face_index = b->tag;

      // Loop over all the nodes on this face
      for ( int jj = 0; jj < 2*mesh_order-1; jj++ ){
        for ( int ii = 0; ii < 2*mesh_order-1; ii++ ){
          // Label only the dependent nodes
          if ((ii % 2 == 1) || (jj % 2 == 1)){
            TMROctant node;
            if (face_index < 2){
              node.x = b->x + (face_index % 2)*h;
              node.y = b->y + ii*hc;
              node.z = b->z + jj*hc; 
            }
            else if (face_index < 4){
              node.x = b->x + ii*hc;
              node.y = b->y + (face_index % 2)*h;
              node.z = b->z + jj*hc; 
            }
            else {
              node.x = b->x + ii*hc;
              node.y = b->y + jj*hc;
              node.z = b->z + (face_index % 2)*h;
            }
            
            // Search for dependent node and label it
            const int use_node_search = 1;
            TMROctant *t = nodes->contains(&node, use_node_search);
            t->tag = -1;
          }
        }
      }
    }

    // Get the array of dependent edges
    dep_edges[owned]->getArray(&dep_array, &dep_size);

    for ( int i = 0; i < dep_size; i++ ){
      TMROctant *b = &dep_array[i];

      // Find the edge lengths of the element octant and
      // the node spacing on the dependent face
      const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);
      const int32_t hc = 1 << (TMR_MAX_LEVEL - b->level - (mesh_order-1));
          
      // Get the edge index
      int edge_index = b->tag;

      // Loop over all the nodes on this edge
      for ( int ii = 0; ii < 2*mesh_order-1; ii++ ){
        // Label only the dependent nodes
        if (ii % 2 == 1){
          TMROctant node;
          if (edge_index < 4){
            node.x = b->x + ii*hc;
            node.y = b->y + h*(edge_index % 2);
            node.z = b->z + h*(edge_index/2);
          }
          else if (edge_index < 8){
            node.x = b->x + h*(edge_index % 2);
            node.y = b->y + ii*hc;
            node.z = b->z + h*((edge_index-4)/2);
          }
          else {
            node.x = b->x + h*(edge_index % 2);
            node.y = b->y + h*((edge_index-8)/2);
            node.z = b->z + ii*hc;
          }
            
          // Search for dependent node and label it
          const int use_node_search = 1;
          TMROctant *t = nodes->contains(&node, use_node_search);
          t->tag = -1;
        }
      }
    }
  }
}

/*
  Return an array of flags to indicate whether the given block owns
  each of the faces, edges and nodes
*/
void TMROctForest::getOwnerFlags( int block,
                                  const int *face_block_owners,
                                  const int *edge_block_owners,
                                  const int *node_block_owners,
                                  int *is_face_owner, 
                                  int *is_edge_owner, 
                                  int *is_node_owner ){
  // Check whether the block owns the node, edge or node
  if (is_node_owner){
    for ( int k = 0; k < 8; k++ ){
      int node = block_conn[8*block + k];
      is_node_owner[k] = (node_block_owners[node] == block);
    }
  }
  
  // Get the edge owners
  if (is_edge_owner){
    for ( int k = 0; k < 12; k++ ){
      int edge = block_edge_conn[12*block + k];
      is_edge_owner[k] = (edge_block_owners[edge] == block);
    }
  }
  
  // Get the face owners
  if (is_face_owner){
    for ( int k = 0; k < 6; k++ ){
      int face = block_face_conn[6*block + k];
      is_face_owner[k] = (face_block_owners[face] == block);
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
  face_owners:  the block index of the face owner
  edge_owners:  the block index of the edge owner
  node_owners:  the block index of the corner/node owner
*/
void TMROctForest::orderGlobalNodes( const int *face_block_owners,
                                     const int *edge_block_owners,
                                     const int *node_block_owners ){
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
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

    // Check whether the block owns the nodes, edges or faces 
    int is_face_owner[6], is_edge_owner[12], is_node_owner[8];
    getOwnerFlags(block, face_block_owners, 
                  edge_block_owners, node_block_owners,
                  is_face_owner, is_edge_owner, is_node_owner);

    // Get the octant array of nodes
    TMROctantArray *nodes;
    octrees[block]->getNodes(&nodes);

    // Extract the array of octants
    int size;
    TMROctant *array;
    nodes->getArray(&array, &size);

    // Count up the number of elements
    num_elements += octrees[block]->getNumElements();

    // Count up the number of interior nodes
    const int32_t hmax = 1 << TMR_MAX_LEVEL;
    for ( int i = 0; i < size; i++ ){
      if (array[i].tag >= 0){
        // Determine which faces the octant touches if any
        int fx0 = (array[i].x == 0);
        int fy0 = (array[i].y == 0);
        int fz0 = (array[i].z == 0);
        int fx = (fx0 || array[i].x == hmax);
        int fy = (fy0 || array[i].y == hmax);
        int fz = (fz0 || array[i].z == hmax);
          
        // Check whether this node is on the face, edge, corner or is
        // an internal node and order it only if it is locally owned
        if (fx && fy && fz){
          int node_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2) + (fz0 ? 0 : 4);
          if (is_node_owner[node_index]){ nlocal++; }
        }
        else if (fy && fz){
          int edge_index = (fy0 ? 0 : 1) + (fz0 ? 0 : 2);
          if (is_edge_owner[edge_index]){ nlocal++; }
        }
        else if (fx && fz){
          int edge_index = (fx0 ? 4 : 5) + (fz0 ? 0 : 2);
          if (is_edge_owner[edge_index]){ nlocal++; }
        }
        else if (fx && fy){
          int edge_index = (fx0 ? 8 : 9) + (fy0 ? 0 : 2);
          if (is_edge_owner[edge_index]){ nlocal++; }
        }
        else if (fx){
          int face_index = (fx0 ? 0 : 1);
          if (is_face_owner[face_index]){ nlocal++; }
        }
        else if (fy){
          int face_index = (fy0 ? 2 : 3);
          if (is_face_owner[face_index]){ nlocal++; }
        }
        else if (fz){
          int face_index = (fz0 ? 4 : 5);
          if (is_face_owner[face_index]){ nlocal++; }
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
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

    // Check whether the block owns the nodes, edges or faces 
    int is_face_owner[6], is_edge_owner[12], is_node_owner[8];
    getOwnerFlags(block, face_block_owners, 
                  edge_block_owners, node_block_owners,
                  is_face_owner, is_edge_owner, is_node_owner);

    // Now order the face. Get the octant array of nodes
    TMROctantArray *nodes;
    octrees[block]->getNodes(&nodes);

    // Extract the array of octants
    int size;
    TMROctant *array;
    nodes->getArray(&array, &size);

    // Count up the number nodes that are locally owned
    // in the interior, edges or corners
    const int32_t hmax = 1 << TMR_MAX_LEVEL;

    for ( int i = 0; i < size; i++ ){
      if (array[i].tag >= 0){
        // Determine which faces the octant touches if any
        int fx0 = (array[i].x == 0);
        int fy0 = (array[i].y == 0);
        int fz0 = (array[i].z == 0);
        int fx = (fx0 || array[i].x == hmax);
        int fy = (fy0 || array[i].y == hmax);
        int fz = (fz0 || array[i].z == hmax);
          
        // Check whether this node is on the face, edge, corner or is
        // an internal node and order it only if it is locally owned
        if (fx && fy && fz){
          int node_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2) + (fz0 ? 0 : 4);
          if (is_node_owner[node_index]){ 
            array[i].tag = node_num; node_num++; 
          }
        }
        else if (fy && fz){
          int edge_index = (fy0 ? 0 : 1) + (fz0 ? 0 : 2);
          if (is_edge_owner[edge_index]){ 
            array[i].tag = node_num; node_num++; 
          }
        }
        else if (fx && fz){
          int edge_index = (fx0 ? 4 : 5) + (fz0 ? 0 : 2);
          if (is_edge_owner[edge_index]){ 
            array[i].tag = node_num; node_num++; 
          }
        }
        else if (fx && fy){
          int edge_index = (fx0 ? 8 : 9) + (fy0 ? 0 : 2);
          if (is_edge_owner[edge_index]){ 
            array[i].tag = node_num; node_num++; 
          }
        }
        else if (fx){
          int face_index = (fx0 ? 0 : 1);
          if (is_face_owner[face_index]){ 
            array[i].tag = node_num; node_num++; 
          }
        }
        else if (fy){
          int face_index = (fy0 ? 2 : 3);
          if (is_face_owner[face_index]){ 
            array[i].tag = node_num; node_num++; 
          }
        }
        else if (fz){
          int face_index = (fz0 ? 4 : 5);
          if (is_face_owner[face_index]){ 
            array[i].tag = node_num; node_num++; 
          }
        }
        else {
          array[i].tag = node_num; node_num++;
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
  node_index:  the local index of the corner 0 <= node index < 12
  block_owner: the block index that owns the node
  p:           the octant with the correct tag value
*/
void TMROctForest::copyCornerNodes( int node,
                                    int node_index,
                                    int block_owner,
                                    TMROctant *p ){
  // Set the maximum octant side length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  for ( int ip = node_block_ptr[node]; ip < node_block_ptr[node+1]; ip++ ){
    int block = node_block_conn[ip];
    if (block_owner != block){
      // Find the corresponding node index for this block
      for ( int j = 0; j < 8; j++ ){
        if (node == block_conn[8*block + j]){
          TMROctant oct;
          oct.x = hmax*(j % 2);
          oct.y = hmax*((j % 4)/2);
          oct.z = hmax*(j/4);

          // Get the node array
          TMROctantArray *nodes;
          octrees[block]->getNodes(&nodes);
          
          // Search the nodes and set the tag value
          const int use_node_search = 1;
          TMROctant *t = nodes->contains(&oct, use_node_search);
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
  edge_index:  the local index of the edge 0 <= edge index < 12
  block_owner: the block index that owns the edge
  p:           the octant with the correct tag value
*/
void TMROctForest::copyEdgeNodes( int edge,
                                  int edge_index,
                                  int block_owner,
                                  TMROctant *p ){
  // Set the maximum octant side length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Store the u coordinate along the edge
  int32_t u = 0;
  if (edge_index < 4){
    u = p->x;
  }
  else if (edge_index < 8){
    u = p->y;
  }
  else {
    u = p->z;
  }

  // Retrieve the first and second node numbers to determine the
  // relative orientation between this edge and each adjacent edge
  int n1 = block_conn[8*block_owner + block_to_edge_nodes[edge_index][0]];
  int n2 = block_conn[8*block_owner + block_to_edge_nodes[edge_index][1]];

  for ( int ip = edge_block_ptr[edge]; ip < edge_block_ptr[edge+1]; ip++ ){
    int block = edge_block_conn[ip]/12;
    if (block_owner != block){
      // Get the adjacent edge index on the opposite block
      int adj_index = edge_block_conn[ip] % 12;
      
      // Get the orientation 
      int nn1 = block_conn[8*block + block_to_edge_nodes[adj_index][0]];
      int nn2 = block_conn[8*block + block_to_edge_nodes[adj_index][1]];
      
      // Determine whether the edges are in the same direction or
      // are reversed
      int reverse = (n1 == nn2 && n2 == nn1);

      // Set the u-coordinate along the edge
      int32_t uoct = u;
      if (reverse){
        uoct = hmax - u;
      }
      
      // Transform the octant to the adjacent coordinate system
      TMROctant oct;
      if (adj_index < 4){
        oct.x = uoct;
        oct.y = hmax*(adj_index % 2);
        oct.z = hmax*(adj_index/2);
      }
      else if (adj_index < 8){
        oct.x = hmax*(adj_index % 2);
        oct.y = uoct;
        oct.z = hmax*((adj_index-4)/2);
      }
      else {
        oct.x = hmax*(adj_index % 2);
        oct.y = hmax*((adj_index-8)/2);
        oct.z = uoct;
      }

      // Get the node array
      TMROctantArray *nodes;
      octrees[block]->getNodes(&nodes);
      
      // Search the nodes and set the tag value
      const int use_node_search = 1;
      TMROctant *t = nodes->contains(&oct, use_node_search);
      t->tag = p->tag;      
    }
  }
}

/*
  Copy the face node variable number on the given face with the
  provided face index to all adjacent faces with the provided indices.
  These values are then communicated back to their owners in a second
  parallel step.

  input:
  face:        the face number that the node lies on
  face_index:  the local index of the face 0 <= face index < 6
  block_owner: the block index that owns the face
  p:           the octant with the correct tag value
*/
void TMROctForest::copyFaceNodes( int face,
                                  int face_index,
                                  int block_owner, 
                                  TMROctant *p ){
  // Set the maximum octant side length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Get the face id number
  int face_id = block_face_ids[6*block_owner + face_index];

  // Get the u/v coordinates for this node on the owner face
  int32_t u, v;
  if (face_index < 2){ // x-face
    get_face_node_coords(face_id, p->y, p->z, &u, &v);
  }
  else if (face_index < 4){ // y-face
    get_face_node_coords(face_id, p->x, p->z, &u, &v);
  }
  else { // z-face
    get_face_node_coords(face_id, p->x, p->y, &u, &v);
  }

  // Loop over all adjacent faces and copy over the node numbers
  // associated with each co-incident node on the adjacent face
  for ( int ip = face_block_ptr[face]; ip < face_block_ptr[face+1]; ip++ ){
    int block = face_block_conn[ip]/6;

    if (block_owner != block){
      // Get the index of the adjacent face
      int adj_index = face_block_conn[ip] % 6;

      // Get the face_id corresponding to the orientation of this
      // adjacent face
      face_id = block_face_ids[6*block + adj_index];

      // Transform the octant p to the local octant coordinates
      TMROctant oct;
      if (adj_index < 2){
        oct.x = hmax*(adj_index % 2);
        set_face_node_coords(face_id, u, v, &oct.y, &oct.z);
      }
      else if (adj_index < 4){
        oct.y = hmax*(adj_index % 2);
        set_face_node_coords(face_id, u, v, &oct.x, &oct.z);
      }
      else {
        oct.z = hmax*(adj_index % 2);
        set_face_node_coords(face_id, u, v, &oct.x, &oct.y);
      }
      
      // Search for the local octant in the list of nodes
      TMROctantArray *nodes;
      octrees[block]->getNodes(&nodes);
      
      // Search the nodes and set the tag value
      const int use_node_search = 1;
      TMROctant *t = nodes->contains(&oct, use_node_search);
      t->tag = p->tag;
    }
  }
}

/*
  Copy the nodes ordered on adjacent face to the locally owned
  adjacent faces.

  This code loops over all the nodes in all the blocks that are owned
  by this processor. When a node on a face, edge or corner is
  encountered that is owned by this processor, its node is copied to
  the non-local interface nodes. These nodes are later communicated to
  their source processors.

  Note that the dependent nodes must be labeled correctly before this
  function will work.

  intput:
  face_owners:  the block index of the face owner
  edge_owners:  the block index of the edge owner
  node_owners:  the block index of the corner/node owner
*/
void TMROctForest::copyAdjacentNodes( const int *face_block_owners,
                                      const int *edge_block_owners,
                                      const int *node_block_owners ){
  // Copy over the block node numbers
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

    // Get the node array
    TMROctantArray *nodes;
    octrees[block]->getNodes(&nodes);
        
    // Get the actual octant array
    int size;
    TMROctant *array;
    nodes->getArray(&array, &size);
    
    // Loop over all the nodes and find the ones that are on the
    // desired face
    for ( int i = 0; i < size; i++ ){
      const int32_t hmax = 1 << TMR_MAX_LEVEL; 
      
      // Only copy the independent nodes that have a non-negative
      // tag node number
      if (array[i].tag >= 0){        
        // Determine which faces the octant touches if any
        int fx0 = (array[i].x == 0);
        int fy0 = (array[i].y == 0);
        int fz0 = (array[i].z == 0);
        int fx = (fx0 || array[i].x == hmax);
        int fy = (fy0 || array[i].y == hmax);
        int fz = (fz0 || array[i].z == hmax);
          
        if (fx && fy && fz){
          int node_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2) + (fz0 ? 0 : 4);
          int node = block_conn[8*block + node_index];
          if (node_block_owners[node] == block){
            copyCornerNodes(node, node_index, block, &array[i]);
          }
        }
        else if ((fy && fz) || (fx && fz) || (fx && fy)){
          int edge_index = 0;
          if (fy && fz){      
            edge_index = (fy0 ? 0 : 1) + (fz0 ? 0 : 2); 
          }
          else if (fx && fz){ 
            edge_index = (fx0 ? 4 : 5) + (fz0 ? 0 : 2); 
          }
          else {
            edge_index = (fx0 ? 8 : 9) + (fy0 ? 0 : 2); 
          }
          int edge = block_edge_conn[12*block + edge_index];
          if (edge_block_owners[edge] == block){
            copyEdgeNodes(edge, edge_index, block, &array[i]);
          }
        }
        else if (fx || fy || fz){
          int face_index = 
            fx*(fx0 ? 0 : 1) + fy*(fy0 ? 2 : 3) + fz*(fz0 ? 4 : 5);
          int face = block_face_conn[6*block + face_index];
          if (face_block_owners[face] == block){
            copyFaceNodes(face, face_index, block, &array[i]);
          }
        }
      }
    }
  }
}

/*
  The following code sends the node numbers back to the block owners.

  This code sends the nodes that are owned by this processor but are
  also on face, edge or corner of an adjacent octree. These nodes have
  been ordered locally, but their node numbers have not yet been
  passed back to the corresponding block owner (that does not own the
  face/edge/node).

  This code scans over blocks that are allocated locally, but are not
  owned locally - these are the partial octrees allocated during the
  recvOctNeighbors call. The nodes on edges/faces/corners of these
  trees are in the local octree order.

  input:
  face_owners:  the block index of the face owner
  edge_owners:  the block index of the edge owner
  node_owners:  the block index of the corner/node owner
*/
void TMROctForest::sendNodeNeighbors( const int *face_block_owners,
                                      const int *edge_block_owners,
                                      const int *node_block_owners ){
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Send all the nodes that have a non-negative node number
  TMROctantQueue **queues = new TMROctantQueue*[ mpi_size ];
  for ( int k = 0; k < mpi_size; k++ ){
    if (k != mpi_rank){
      queues[k] = new TMROctantQueue();
    }
    else {
      queues[k] = NULL;
    }
  }

  // Count up the number of nodes owned on each face
  for ( int block = 0; block < num_blocks; block++ ){
    if (octrees[block] && 
        (mpi_rank != mpi_block_owners[block])){
      // The destination rank
      int rank = mpi_block_owners[block];

      // Get the element array
      TMROctantArray *nodes;
      octrees[block]->getNodes(&nodes);
      
      // Get the actual octant array
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);

      // Loop over all the nodes and check where they need to be sent
      for ( int i = 0; i < size; i++ ){
        if (array[i].tag >= 0){
          const int32_t hmax = 1 << TMR_MAX_LEVEL; 
          
          // Determine which faces the octant touches if any
          int fx0 = (array[i].x == 0);
          int fy0 = (array[i].y == 0);
          int fz0 = (array[i].z == 0);
          int fx = (fx0 || array[i].x == hmax);
          int fy = (fy0 || array[i].y == hmax);
          int fz = (fz0 || array[i].z == hmax);
        
          // Check whether the node lies on a corner, edge or face
          if (fx || fy || fz){
            // This is cheating -- set the level -- which is not needed
            // for the nodes to the destination block number
            TMROctant q = array[i];
            q.level = block;
            
            // Pass the block to any adjacent octrees
            int block_owner = -1;
            if (fx && fy && fz){
              int node_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2) + (fz0 ? 0 : 4);
              int node = block_conn[8*block + node_index];
              block_owner = node_block_owners[node];
            }
            else if (fy && fz){
              int edge_index = (fy0 ? 0 : 1) + (fz0 ? 0 : 2);
              int edge = block_edge_conn[12*block + edge_index];
              block_owner = edge_block_owners[edge];
            }
            else if (fx && fz){
              int edge_index = (fx0 ? 4 : 5) + (fz0 ? 0 : 2);
              int edge = block_edge_conn[12*block + edge_index];
              block_owner = edge_block_owners[edge];
            }
            else if (fx && fy){
              int edge_index = (fx0 ? 8 : 9) + (fy0 ? 0 : 2);
              int edge = block_edge_conn[12*block + edge_index];
              block_owner = edge_block_owners[edge];
            }
            else if (fx){
              int face_index = (fx0 ? 0 : 1);
              int face = block_face_conn[6*block + face_index];
              block_owner = face_block_owners[face];
            }
            else if (fy){
              int face_index = (fy0 ? 2 : 3);
              int face = block_face_conn[6*block + face_index];
              block_owner = face_block_owners[face];
            }
            else if (fz){
              int face_index = (fz0 ? 4 : 5);
              int face = block_face_conn[6*block + face_index];
              block_owner = face_block_owners[face];
            }

            // If this processor owns the node/edge/face then it will
            // have the correct node number and it needs to be sent back
            // to the block who sent this in the first place
            if (block_owner >= 0 && 
                (mpi_rank == mpi_block_owners[block_owner])){            
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
  TMROctantQueue **qnodes = new TMROctantQueue*[ num_blocks ];
  memset(qnodes, 0, num_blocks*sizeof(TMROctantQueue*));

  // Receive the arrays of incoming nodes
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
      // For this one function, the level encodes the block owner for
      // the node. We now zero this so that there is no
      // confusion. This preserves the tag value which stores the node
      // number itself.
      int block = array[i].level;
      array[i].level = 0;
      if (!qnodes[block]){
        qnodes[block] = new TMROctantQueue();
      }
      qnodes[block]->push(&array[i]);
    }
    delete [] array;
  }

  // Loop over all the blocks that have been received 
  for ( int block = 0; block < num_blocks; block++ ){
    if (qnodes[block]){
      // Convert the queue to an array and delete the old queue
      TMROctantArray *new_nodes = qnodes[block]->toArray();
      delete qnodes[block];

      // Sort the new array so that we can pass the node numbers
      new_nodes->sort();

      // Get the nodes from the tree
      TMROctantArray *nodes;
      octrees[block]->getNodes(&nodes);

      // Get the arrays of the actual octants
      int size, new_size;
      TMROctant *array, *new_array;
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

  This function first computes the block that owns of each of the
  faces, edges and corners(/nodes) within the super-mesh. Next, the
  the non-local octants that border each octree are passed back to the
  processors which own the octree. This creates a layer of octants
  that are temporarily stored in a partial octree. The code then
  creates nodes within each octant element for all of the octrees
  (including the partial octrees). Next, the dependent nodes (that are
  hanging on a face or an edge of an element)are labeled according to
  whether they are on an edge or face.
  
  After the dependent nodes are labeled, the nodes are ordered on all
  processors to form a complete global ordering. Next, the nodes are
  communicated locally across octant and partial octant faces, edges
  and corners. Finally, the new node numbers are returned to the
  processors that border the octree owners. And lastly, the non-local
  partial octrees are freed.  

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
      int block = face_block_conn[ip]/6;
      if (block < face_block_owners[face]){
        face_block_owners[face] = block;
      }
    }
  }

  // Find the edge owners
  for ( int edge = 0; edge < num_edges; edge++ ){
    edge_block_owners[edge] = num_blocks;

    int ipend = edge_block_ptr[edge+1];
    for ( int ip = edge_block_ptr[edge]; ip < ipend; ip++ ){
      int block = edge_block_conn[ip]/12;
      if (block < edge_block_owners[edge]){
        edge_block_owners[edge] = block;
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
  recvOctNeighbors();

  // Label the dependent nodes on all the blocks
  computeDepFacesAndEdges();

  // Allocate all possible nodes on all of the trees, including the
  // partial trees that have just been exchanged.
  for ( int block = 0; block < num_blocks; block++ ){
    if (octrees[block]){
      octrees[block]->createNodes(mesh_order);
    }
  }

  // Label the dependent nodes
  labelDependentNodes();

  // Compute the global ordering of the nodes
  orderGlobalNodes(face_block_owners, edge_block_owners,
                   node_block_owners);

  // Copy the adjacent node numbers between local blocks
  copyAdjacentNodes(face_block_owners, edge_block_owners,
                    node_block_owners);

  // Send the nodes (complete with global number) back to whence they
  // came to finish the global ordering
  sendNodeNeighbors(face_block_owners, edge_block_owners,
                    node_block_owners);

  // Free and NULL the trees that are not local - these are not
  // required anymore.
  for ( int block = 0; block < num_blocks; block++ ){
    if (octrees[block] && 
        (mpi_rank != mpi_block_owners[block])){
      delete octrees[block];
      octrees[block] = NULL;
    }
  }
  
  delete [] face_block_owners;
  delete [] edge_block_owners;
  delete [] node_block_owners;
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
void TMROctForest::createMeshConn( int **_conn, int *_nelems ){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // First, count up the number of elements
  int nelems = 0;
  for ( int block = 0; block < num_blocks; block++ ){
    if (mpi_rank == mpi_block_owners[block]){
      nelems += octrees[block]->getNumElements();
    }
  }

  // Allocate the connectivity
  int conn_size = 0;
  int *elem_conn = new int[ mesh_order*mesh_order*mesh_order*nelems ];

  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

    // Retrieve the element array and the node array
    TMROctantArray *elements, *nodes;
    octrees[block]->getElements(&elements);
    octrees[block]->getNodes(&nodes);

    // Get the current array of octants
    int size;
    TMROctant *array;
    elements->getArray(&array, &size);

    for ( int i = 0; i < size; i++ ){
      // For all searches/comparisons, we use node numbers
      const int use_nodes = 1;
          
      // Add all of the nodes from this element
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level - (mesh_order-2));
      TMROctant p;
      p.level = array[i].level;
          
      for ( int kk = 0; kk < mesh_order; kk++ ){
        for ( int jj = 0; jj < mesh_order; jj++ ){
          for ( int ii = 0; ii < mesh_order; ii++ ){
            p.x = array[i].x + ii*h;
            p.y = array[i].y + jj*h;
            p.z = array[i].z + kk*h;
            
            // Get the index for the node
            TMROctant *t = nodes->contains(&p, use_nodes);
            
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
int TMROctForest::getDepNodeConn( const int **ptr, const int **conn,
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
void TMROctForest::createDepNodeConn( int **_ptr, int **_conn, 
                                      double **_weights ){
  // Get the MPI rank
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Count up the total number of dependent faces
  int ndep_faces = 0, ndep_edges = 0;
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int ndep = 0;
    dep_faces[owned]->getArray(NULL, &ndep);
    ndep_faces += ndep;
    dep_edges[owned]->getArray(NULL, &ndep);
    ndep_edges += ndep;
  }

  // Allocate the space for the dependent variable information
  const int nodes_per_edge = (2*mesh_order-1)*(2*mesh_order-1);
  const int nodes_per_face = (2*mesh_order-1)*(2*mesh_order-1);
  int *face_nodes = new int[ nodes_per_face*ndep_faces ];
  int *edge_nodes = new int[ nodes_per_edge*ndep_edges ];

  // Retrieve the block information from the dependent faces
  int face_offset = 0, edge_offset = 0;
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

    // Get the nodes associated with this octree
    TMROctantArray *nodes;
    octrees[block]->getNodes(&nodes);

    // Get the dependent faces
    int dep_size;
    TMROctant *dep_array;
    dep_faces[owned]->getArray(&dep_array, &dep_size);

    for ( int i = 0; i < dep_size; i++, face_offset++ ){
      TMROctant *b = &dep_array[i];

      // Get the face index
      int face_index = b->tag;

      // Find the edge length of the element octant
      const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);

      // Find the edge length/node separation along the face
      const int32_t hc = 1 << (TMR_MAX_LEVEL - b->level - (mesh_order-1));
      
      // Loop over all the nodes on this face
      int *fn = &face_nodes[nodes_per_face*face_offset];
      for ( int jj = 0; jj < (2*mesh_order-1); jj++ ){
        for ( int ii = 0; ii < (2*mesh_order-1); ii++ ){
          // Get the node location corresponding to this face
          TMROctant node;
          if (face_index < 2){
            node.x = b->x + (face_index % 2)*h;
            node.y = b->y + ii*hc;
            node.z = b->z + jj*hc; 
          }
          else if (face_index < 4){
            node.x = b->x + ii*hc;
            node.y = b->y + (face_index % 2)*h;
            node.z = b->z + jj*hc; 
          }
          else {
            node.x = b->x + ii*hc;
            node.y = b->y + jj*hc;
            node.z = b->z + (face_index % 2)*h;
          }

          // Search for node octant and record its label
          const int use_node_search = 1;
          TMROctant *t = nodes->contains(&node, use_node_search);
          fn[0] = t->tag;
          fn++;
        }
      }
    }

    // Get the dependent edges
    dep_edges[owned]->getArray(&dep_array, &dep_size);

    for ( int i = 0; i < dep_size; i++, edge_offset++ ){
      TMROctant *b = &dep_array[i];

      // Get the edge index
      int edge_index = b->tag;

      // Find the edge length of the element octant
      const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);
      const int32_t hc = 1 << (TMR_MAX_LEVEL - b->level - (mesh_order-1));
      
      // Loop over all the nodes on this face
      int *en = &edge_nodes[nodes_per_edge*edge_offset];
      for ( int ii = 0; ii < (2*mesh_order-1); ii++ ){
        // Get the node location corresponding to this face
        TMROctant node;
        if (edge_index < 4){
          node.x = b->x + ii*hc;
          node.y = b->y + h*(edge_index % 2);
          node.z = b->z + h*(edge_index/2);
        }
        else if (edge_index < 8){
          node.x = b->x + h*(edge_index % 2);
          node.y = b->y + ii*hc;
          node.z = b->z + h*((edge_index-4)/2);
        }
        else {
          node.x = b->x + h*(edge_index % 2);
          node.y = b->y + h*((edge_index-8)/2);
          node.z = b->z + ii*hc;
        }

        // Search for node octant and record its label
        const int use_node_search = 1;
        TMROctant *t = nodes->contains(&node, use_node_search);
        en[0] = t->tag;
        en++;
      }
    }
  }

  // Allocate the pointer to 
  int *ptr = new int[ num_dep_nodes+1 ];
  memset(ptr, 0, (num_dep_nodes+1)*sizeof(int));

  for ( int i = 0; i < ndep_faces; i++ ){
    // Loop over all the nodes on the dependent face
    const int *fn = &face_nodes[nodes_per_face*i];
    for ( int jj = 0; jj < 2*mesh_order-1; jj++ ){
      for ( int ii = 0; ii < 2*mesh_order-1; ii++ ){
        // Nodes where both indices are odd dependent on mesh_order**2
        // nodes, while if only one index is odd, then it depends on only
        // mesh_order nodes
        int node = -fn[ii + (2*mesh_order-1)*jj]-1;
        if ((ii % 2 == 1) && (jj % 2 == 1)){
          ptr[node+1] = mesh_order*mesh_order;
        }
        else if ((ii % 2 == 1) || (jj % 2 == 1)){
          ptr[node+1] = mesh_order;
        }
      }
    }
  }

  for ( int i = 0; i < ndep_edges; i++ ){
    // Loop over all the nodes on the dependent edges
    const int *en = &edge_nodes[nodes_per_edge*i];
    for ( int ii = 0; ii < 2*mesh_order-1; ii++ ){
      int node = -en[ii]-1;
      if (ii % 2 == 1){
        ptr[node+1] = mesh_order;
      }
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

  // Go through and assign the connectivity and weights
  for ( int k = 0; k < ndep_faces; k++ ){
    const int n = 2*mesh_order-1;
    const int m = mesh_order-1;

    // Loop over all the nodes on the dependent face
    for ( int jj = 0; jj < n; jj++ ){
      for ( int ii = 0; ii < n; ii++ ){
        const int *fn = &face_nodes[nodes_per_face*k];
        
        // Get the dependent node number
        int node = -fn[ii + n*jj]-1;
        if ((ii % 2 == 1) && (jj % 2 == 1)){
          for ( int j = 0, kk = 0; j < mesh_order; j++ ){
            for ( int i = 0; i < mesh_order; kk++, i++ ){
              conn[ptr[node] + kk] = fn[2*i + (2*j)*n];
              weights[ptr[node] + kk] = wt[ii/2][i]*wt[jj/2][j];
            }
          }
        }
        else if (ii % 2 == 1){
          for ( int i = 0; i < mesh_order; i++ ){
            conn[ptr[node] + i] = fn[2*i + jj*n];
            weights[ptr[node] + i] = wt[ii/2][i];
          }
        }
        else if (jj % 2 == 1){
          for ( int j = 0; j < mesh_order; j++ ){
            conn[ptr[node] + j] = fn[ii + (2*j)*n];
            weights[ptr[node] + j] = wt[jj/2][j];
          }
        }
      }
    }
  }

  for ( int k = 0; k < ndep_edges; k++ ){
    // Loop over all the nodes on the dependent edges
    const int *en = &edge_nodes[nodes_per_edge*k];
    for ( int ii = 0; ii < 2*mesh_order-1; ii++ ){
      int node = -en[ii]-1;
      if (ii % 2 == 1){
        for ( int i = 0; i < mesh_order; i++ ){
          conn[ptr[node] + i] = en[2*i];
          weights[ptr[node] + i] = wt[ii/2][i];
        }
      }
    }
  }
  
  // Free the temporary memory
  delete [] edge_nodes;
  delete [] face_nodes;

  // Return the weights
  *_ptr = ptr;
  *_conn = conn;
  *_weights = weights;
}

/*
  Add the node weights to the array. This eliminates dependent nodes
  by unrolling the dependency information.

  input:
  t:             the octant node
  w:             the weight for this node
  cdep_ptr:      the dependent node pointer for the coarse mesh
  cdep_conn:     the dependent node connectivity for the coarse mesh
  cdep_weights:  the dependent node weights for the coarse mesh

  input/output: 
  weights:       the list of index weights
  nweights:      the number of weights in the list
*/
void TMROctForest::addNodeWeights( TMROctant *t, double w,
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
  coarse:   the coarse octree forest that has the same layout as this
  
  output:
  ptr:      the pointer into the local rows
  conn:     the connectivity using global numbers
  weights:  the interpolation weights for each point
*/
void TMROctForest::createInterpolation( TMROctForest *coarse,
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
    mesh_order*mesh_order*mesh_order*nnodes;

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
  int max_size = (mesh_order*mesh_order*mesh_order)*(mesh_order*mesh_order);
  TMRIndexWeight *wlist = new TMRIndexWeight[ max_size ];

  // Loop over all the blocks and compute the interpolation
  int count = 0;
  for ( int owned = 0; owned < num_owned_blocks; owned++ ){
    int block = owned_blocks[owned];

    // Get the nodes from the fine and coarse mesh
    TMROctantArray *nodes, *coarse_nodes;
    octrees[block]->getNodes(&nodes);
    coarse->octrees[block]->getNodes(&coarse_nodes);

    // Get the size of the fine node array
    int fine_size;
    TMROctant *fine;
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

        // Retrieve the octant pointer
        TMROctant *t = coarse_nodes->contains(&fine[i], use_node_search);
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
          int32_t pz = fine[i].z % h;

          // Determine which interpolation to use
          int32_t sx = fine[i].x % hc;
          int32_t sy = fine[i].y % hc;
          int32_t sz = fine[i].z % hc;

          // Add to the interpolation depending on the values of px,
          // py, and pz
          if (sx && sy && sz){
            for ( int kk = 0; kk < mesh_order; kk++ ){
              for ( int jj = 0; jj < mesh_order; jj++ ){
                for ( int ii = 0; ii < mesh_order; ii++ ){
                  TMROctant node = fine[i];
                  node.x = (fine[i].x - px) + hc*ii;
                  node.y = (fine[i].y - py) + hc*jj;
                  node.z = (fine[i].z - pz) + hc*kk;
                  t = coarse_nodes->contains(&node, use_node_search);
                  double w = wt[px/hc][ii]*wt[py/hc][jj]*wt[pz/hc][kk];
                  addNodeWeights(t, w, cdep_ptr, cdep_conn, cdep_weights,
                                 wlist, &nweights);
                }
              }
            }
          }
          else if (sy && sz){
            for ( int kk = 0; kk < mesh_order; kk++ ){
              for ( int jj = 0; jj < mesh_order; jj++ ){
                TMROctant node = fine[i];
                node.y = (fine[i].y - py) + hc*jj;
                node.z = (fine[i].z - pz) + hc*kk;
                t = coarse_nodes->contains(&node, use_node_search);
                double w = wt[py/hc][jj]*wt[pz/hc][kk];
                addNodeWeights(t, w, cdep_ptr, cdep_conn, cdep_weights,
                               wlist, &nweights);
              }
            }          
          }
          else if (sx && sz){
            for ( int kk = 0; kk < mesh_order; kk++ ){
              for ( int ii = 0; ii < mesh_order; ii++ ){
                TMROctant node = fine[i];
                node.x = (fine[i].x - px) + hc*ii;
                node.z = (fine[i].z - pz) + hc*kk;
                t = coarse_nodes->contains(&node, use_node_search);
                double w = wt[px/hc][ii]*wt[pz/hc][kk];
                addNodeWeights(t, w, cdep_ptr, cdep_conn, cdep_weights,
                               wlist, &nweights);
              }            
            }
          }
          else if (sx && sy){
            for ( int jj = 0; jj < mesh_order; jj++ ){
              for ( int ii = 0; ii < mesh_order; ii++ ){
                TMROctant node = fine[i];
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
              TMROctant node = fine[i];
              node.x = (fine[i].x - px) + hc*ii;
              t = coarse_nodes->contains(&node, use_node_search);
              addNodeWeights(t, wt[px/hc][ii], cdep_ptr, cdep_conn, cdep_weights,
                             wlist, &nweights);
            }
          }
          else if (sy){
            for ( int jj = 0; jj < mesh_order; jj++ ){
              TMROctant node = fine[i];
              node.y = (fine[i].y - py) + hc*jj;
              t = coarse_nodes->contains(&node, use_node_search);
              addNodeWeights(t, wt[py/hc][jj], cdep_ptr, cdep_conn, cdep_weights,
                             wlist, &nweights);
            }
          }
          else if (sz){
            for ( int kk = 0; kk < mesh_order; kk++ ){
              TMROctant node = fine[i];
              node.z = (fine[i].z - pz) + hc*kk;
              t = coarse_nodes->contains(&node, use_node_search);
              addNodeWeights(t, wt[pz/hc][kk], cdep_ptr, cdep_conn, cdep_weights,
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
