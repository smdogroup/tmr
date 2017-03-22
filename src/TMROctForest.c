#include "TMROctForest.h"

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
  Compare tags for sorting
*/
static int compare_octant_tags( const void *a, const void *b ){
  const TMROctant *A = static_cast<const TMROctant*>(a);
  const TMROctant *B = static_cast<const TMROctant*>(b);
  return A->tag - B->tag;
}

/*
  Create the TMROctForest object
*/
TMROctForest::TMROctForest( MPI_Comm _comm ){
  // Initialize the TMR-specific MPI data types
  if (!TMRIsInitialized){
    TMRInitialize();
  }

  // Set the MPI communicator
  comm = _comm;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Set the topology object to NULL to begin with
  topo = NULL;

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
  face_block_owners = NULL;
  edge_block_owners = NULL;
  node_block_owners = NULL;

  // Null the octant owners/octant list
  owners = NULL;
  octants = NULL;
  adjacent = NULL;
  nodes = NULL;
  dep_faces = NULL;
  dep_edges = NULL;
  X = NULL;

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
  // Free the topology object (if one exists)
  if (topo){ topo->decref(); }

  freeData();
}

/*
  Free any data that has been allocated
*/
void TMROctForest::freeData(){
  // Free the connectivity data
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

  // Free the ownership data
  if (face_block_owners){ delete [] face_block_owners; }
  if (edge_block_owners){ delete [] edge_block_owners; }
  if (node_block_owners){ delete [] node_block_owners; }

  // Free the octants/adjacency/dependency data
  if (owners){ delete [] owners; }
  if (octants){ delete octants; }
  if (adjacent){ delete adjacent; }
  if (nodes){ delete nodes; }
  if (dep_faces){ delete dep_faces; }
  if (dep_edges){ delete dep_edges; }
  if (X){ delete [] X; }

  if (node_range){ delete [] node_range; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }

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
  face_block_owners = NULL;
  edge_block_owners = NULL;
  node_block_owners = NULL;

  // Null the octant owners/octant list
  owners = NULL;
  octants = NULL;
  adjacent = NULL;
  nodes = NULL;
  dep_faces = NULL;
  dep_edges = NULL;
  X = NULL;

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
  Copy the connectivity data, but not the octants/nodes
*/
void TMROctForest::copyData( TMROctForest *copy ){
  // Copy over the connectivity data
  copy->num_nodes = num_nodes;
  copy->num_edges = num_edges;
  copy->num_faces = num_faces;
  copy->num_blocks = num_blocks;

  // Allocate/copy the block connectivities
  copy->block_conn = new int[ 8*num_blocks ];
  copy->block_face_conn = new int[ 6*num_faces ];
  copy->block_face_ids = new int[ 6*num_faces ];
  copy->block_edge_conn = new int[ 12*num_faces ];
  memcpy(copy->block_conn, block_conn, 8*num_blocks*sizeof(int));
  memcpy(copy->block_face_conn, block_face_conn, 
         6*num_blocks*sizeof(int));
  memcpy(copy->block_face_ids, block_face_ids, 
         6*num_blocks*sizeof(int));
  memcpy(copy->block_edge_conn, block_edge_conn, 
         12*num_blocks*sizeof(int));
    
  // Allocate/copy the inverse relationships
  copy->node_block_ptr = new int[ num_nodes+1 ];
  copy->node_block_conn = new int[ node_block_ptr[num_nodes] ];
  memcpy(copy->node_block_ptr, node_block_ptr, 
         (num_nodes+1)*sizeof(int));
  memcpy(copy->node_block_conn, node_block_conn, 
         node_block_ptr[num_nodes]*sizeof(int));

  copy->edge_block_ptr = new int[ num_edges+1 ];
  copy->edge_block_conn = new int[ edge_block_ptr[num_edges] ];
  memcpy(copy->edge_block_ptr, edge_block_ptr, 
         (num_edges+1)*sizeof(int));
  memcpy(copy->edge_block_conn, edge_block_conn, 
         edge_block_ptr[num_edges]*sizeof(int));

  copy->face_block_ptr = new int[ num_faces+1 ];
  copy->face_block_conn = new int[ face_block_ptr[num_faces] ];
  memcpy(copy->face_block_ptr, face_block_ptr, 
         (num_faces+1)*sizeof(int));
  memcpy(copy->face_block_conn, face_block_conn, 
         face_block_ptr[num_faces]*sizeof(int));

  // Copy the ownership information
  copy->face_block_owners = new int[ num_faces ];
  copy->edge_block_owners = new int[ num_edges ];
  copy->node_block_owners = new int[ num_nodes ];
  memcpy(copy->face_block_owners, face_block_owners,
         num_faces*sizeof(int));
  memcpy(copy->edge_block_owners, edge_block_owners,
         num_edges*sizeof(int));
  memcpy(copy->node_block_owners, node_block_owners,
         num_nodes*sizeof(int));

  // Copy over the topology object
  copy->topo = topo;
  if (copy->topo){
    copy->topo->incref();
  }
}

/*
  Set the mesh topology object

  This sets to internal topology object (used to define the geometry)
  and resets the internal connectivity (if any has been defined).
*/
void TMROctForest::setTopology( TMRTopology *_topo ){
  if (_topo){
    _topo->incref();
    if (topo){ topo->decref(); }
    topo = _topo;

    // Compute the connectivity
    int _num_nodes, _num_edges, _num_faces, _num_blocks;
    const int *_block_conn, *_block_edge_conn;
    const int *_block_face_conn;
    topo->getConnectivity(&_num_nodes, &_num_edges, 
                          &_num_faces, &_num_blocks,
                          &_block_conn, &_block_edge_conn,
                          &_block_face_conn);

    // Set the connectivity internally
    setFullConnectivity(_num_nodes, _num_edges,
                        _num_faces, _num_blocks,
                        _block_conn, _block_edge_conn,
                        _block_face_conn);
  }
}

/*
  Set the connectivity of the blocks

  This call is collective on all processors. Every processor must make
  a call with the same connectivity information, otherwise the
  inter-octree information will be inconsistent.

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
  freeData();

  // Copy over the data locally
  num_nodes = _num_nodes;
  num_edges = 0;
  num_faces = 0;
  num_blocks = _num_blocks;
  
  // Copy over the block connectivity
  block_conn = new int[ 8*num_blocks ];
  memcpy(block_conn, _block_conn, 8*num_blocks*sizeof(int));

  // Compute the node to block information
  computeNodesToBlocks();

  // Compute the edge connectivity from the block data
  computeEdgesFromNodes();
  computeEdgesToBlocks();

  // Compute the face connectivity from the block data
  computeFacesFromNodes();
  computeFacesToBlocks();

  // Compute the block owners based on the node, edge and face data
  computeBlockOwners();  
}

/*
  Set the full connectivity, specifying the node, edge, face and block
  numbers independently.
*/
void TMROctForest::setFullConnectivity( int _num_nodes, 
                                        int _num_edges,
                                        int _num_faces, 
                                        int _num_blocks,
                                        const int *_block_conn,
                                        const int *_block_edge_conn,
                                        const int *_block_face_conn ){
  // Free any data allocated for other connectivities
  freeData();

  // Copy over the number of geometric entities
  num_nodes = _num_nodes;
  num_edges = _num_edges;
  num_faces = _num_faces;
  num_blocks = _num_blocks;

  // Copy over the block connectivity
  block_conn = new int[ 8*num_blocks ];
  memcpy(block_conn, _block_conn, 8*num_blocks*sizeof(int));

  // Compute the node to block information
  computeNodesToBlocks();

  // Copy over the edge information
  block_edge_conn = new int[ 12*num_blocks ];
  memcpy(block_edge_conn, _block_edge_conn, 12*num_blocks*sizeof(int));
  
  // Compute the edge to block information
  computeEdgesToBlocks();

  // Compute the face connectivity from the block data
  block_face_conn = new int[ 6*num_blocks ];
  memcpy(block_face_conn, _block_face_conn, 6*num_blocks*sizeof(int));

  // Compute the face to block information
  computeFacesToBlocks();

  // Compute the block owners based on the node, edge and face data
  computeBlockOwners();  
}

/*
  Compute the connectivity storing the node to block information 
*/
void TMROctForest::computeNodesToBlocks(){
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

  // Loop over all the blocks and reset node->block connectivity
  // to store both the adjacent block and the corresponding
  // node index into that array
  for ( int node = 0; node < num_nodes; node++ ){
    for ( int ip = node_block_ptr[node];
          ip < node_block_ptr[node+1]; ip++ ){
      int adj = node_block_conn[ip];
      int adj_index = 0;
      for ( ; adj_index < 8; adj_index++ ){
        if (block_conn[8*adj + adj_index] == node){
          break;
        }
      }

      node_block_conn[ip] = 8*adj + adj_index;
    }
  }
}

/*
  Compute the reverse relationship from the edges to the blocks
*/
void TMROctForest::computeEdgesToBlocks(){
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
}

/*
  Establish a unique ordering of the edges along each block
*/
void TMROctForest::computeEdgesFromNodes(){
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
          int ii = node_block_conn[ip]/8;
          
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
}

/*
  Establish a unique ordering of the faces for each block
*/
void TMROctForest::computeFacesFromNodes(){
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
          int ii = node_block_conn[ip]/8;

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
}

/*
  Compute the inverse connectivity from the faces to the blocks 
*/
void TMROctForest::computeFacesToBlocks(){
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
    for ( int ip = face_block_ptr[face]; 
          ip < face_block_ptr[face+1]; ip++ ){
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
  Once the connectivity for the blocks/faces/edges have been set,
  compute the owners of each object.  
*/
void TMROctForest::computeBlockOwners(){
  // Find the block numbers corresponding to the owner for each face,
  // edge and corner so that we know who should be ordering what!
  face_block_owners = new int[ num_faces ];
  edge_block_owners = new int[ num_edges ];
  node_block_owners = new int[ num_nodes ];

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
      int block = node_block_conn[ip]/8;
      if (block < node_block_owners[node]){
        node_block_owners[node] = block;
      }
    }
  }
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
  Allocate the trees for each element within the mesh
*/
void TMROctForest::createTrees( int refine_level ){
  int32_t level = refine_level;
  if (level < 0){
    level = 0;
  }
  else if (level >= TMR_MAX_LEVEL){
    level = TMR_MAX_LEVEL-1;
  }

  // Set who owns what blocks
  int nblocks = num_blocks/mpi_size;
  int remain = num_blocks % mpi_size;
  int start = mpi_rank*nblocks;
  int end = (mpi_rank+1)*nblocks;
  if (mpi_rank < remain){
    nblocks += 1;
    start += mpi_rank;
    end += mpi_rank+1;
  }
  else {
    start += remain;
    end += remain;
  }
  
  // Create an array of the octants that will be stored
  int nelems = 1 << (TMR_MAX_LEVEL - refine_level);
  int size = nelems*nelems*nelems*nblocks;
  TMROctant *array = new TMROctant[ size ];

  // Generate all of the octants on the associated blocks
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - level);
  for ( int count = 0, block = start; block < end; block++ ){
    for ( int32_t x = 0; x < hmax; x += h ){
      for ( int32_t y = 0; y < hmax; y += h ){
        for ( int32_t z = 0; z < hmax; z += h ){
          array[count].tag = 0;
          array[count].block = block;
          array[count].level = level;
          array[count].x = x;
          array[count].y = y;
          array[count].z = z;
          count++;
        }
      }
    }
  }

  // Create the array of octants
  octants = new TMROctantArray(array, size);
  octants->sort();

  // Set the last octant
  TMROctant p;
  p.tag = -1;
  p.block = num_blocks-1;
  p.x = p.y = p.z = hmax;
  if (size > 0){
    p = array[0];
  }

  owners = new TMROctant[ mpi_size ];
  MPI_Allgather(&p, 1, TMROctant_MPI_type, 
                owners, 1, TMROctant_MPI_type, comm);

  // Set the offsets if some of the processors have zero
  // octants
  for ( int k = 1; k < mpi_size; k++ ){
    if (owners[k].tag == -1){
      owners[k] = owners[k-1];
    }
  }
}

/*
  Create a set of random trees
*/
void TMROctForest::createRandomTrees( int nrand,
                                      int min_level, int max_level ){
  // Set who owns what blocks
  int nblocks = num_blocks/mpi_size;
  int remain = num_blocks % mpi_size;
  int start = mpi_rank*nblocks;
  int end = (mpi_rank+1)*nblocks;
  if (mpi_rank < remain){
    nblocks += 1;
    start += mpi_rank;
    end += mpi_rank+1;
  }
  else {
    start += remain;
    end += remain;
  }
  
  // Create an array of the octants that will be stored
  int size = nrand*nblocks;
  TMROctant *array = new TMROctant[ size ];

  // Generate a random number of octants along random directions
  for ( int count = 0, block = start; block < end; block++ ){
    for ( int i = 0; i < nrand; i++, count++ ){
      int32_t level = min_level + (rand() % (max_level - min_level + 1));

      const int32_t h = 1 << (TMR_MAX_LEVEL - level);
      int32_t x = h*(rand() % (1 << level));
      int32_t y = h*(rand() % (1 << level));
      int32_t z = h*(rand() % (1 << level));
      
      array[count].tag = 0;
      array[count].block = block;
      array[count].level = level;
      array[count].x = x;
      array[count].y = y;
      array[count].z = z;
    }
  }

  // Create the array of octants
  octants = new TMROctantArray(array, size);
  octants->sort();

  // Set the last octant
  TMROctant p;
  p.tag = -1;
  p.block = num_blocks-1;
  p.x = p.y = p.z = 1 << TMR_MAX_LEVEL;
  if (size > 0){
    p = array[0];
  }

  owners = new TMROctant[ mpi_size ];
  MPI_Allgather(&p, 1, TMROctant_MPI_type, 
                owners, 1, TMROctant_MPI_type, comm);

  // Set the offsets if some of the processors have zero
  // octants
  for ( int k = 1; k < mpi_size; k++ ){
    if (owners[k].tag == -1){
      owners[k] = owners[k-1];
    }
  }
}

/*
  Repartition the octants across all processors
*/
void TMROctForest::repartition(){
  // First, this stores the number of elements on octrees owned on
  // each processor
  int *ptr = new int[ mpi_size+1 ];
  int size;
  TMROctant *array;
  octants->getArray(&array, &size);

  // Gather the sizes from all the arrays
  MPI_Allgather(&size, 1, MPI_INT, &ptr[1], 1, MPI_INT, comm);

  // Set the pointers
  ptr[0] = 0;
  for ( int k = 0; k < mpi_size; k++ ){
    ptr[k+1] += ptr[k];
  }

  // Compute the average size of the new counts
  int average_count = ptr[mpi_size]/mpi_size;
  int remain = ptr[mpi_size] - average_count*mpi_size;
  
  // Figure out what goes where on the new distribution of octants
  int *new_ptr = new int[ mpi_size+1 ];
  new_ptr[0] = 0;
  for ( int k = 0; k < mpi_size; k++ ){
    new_ptr[k+1] = new_ptr[k] + average_count;
    if (k < remain){
      new_ptr[k+1] += 1;
    }
  }

  // Allocate the new array of octants
  int new_size = new_ptr[mpi_rank+1] - new_ptr[mpi_rank];
  TMROctant *new_array = new TMROctant[ new_size ];

  // Ptr:      |----|---|--------------------|-| 
  // New ptr:  |-------|-------|-------|-------|

  // Count up the number of sends/recvs
  int send_count = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank && 
        new_ptr[i+1] >= ptr[mpi_rank] &&
        new_ptr[i] < ptr[mpi_rank+1]){
      send_count++;
    }
  }

  // Allocate arrays send/recv arrays
  MPI_Request *send_requests = new MPI_Request[ send_count ];

  // Zero the send/recv counts
  send_count = 0;

  // Match up the intervals for sends
  for ( int i = 0; i < mpi_size; i++ ){
    // If we're on the correct interval, 
    if (new_ptr[i+1] >= ptr[mpi_rank] &&
        new_ptr[i] < ptr[mpi_rank+1]){
      // The start and end locations within the local array
      int start = new_ptr[i] - ptr[mpi_rank];
      if (start < 0){ start = 0; }

      int max_count = ptr[mpi_rank+1] - ptr[mpi_rank];
      int end = new_ptr[i+1] - ptr[mpi_rank];
      if (end > max_count){ end = max_count; }

      // Set the count
      int count = end - start;
      
      if (i == mpi_rank){
        int new_start = ptr[i] - new_ptr[i];
        if (new_start < 0){ new_start = 0; }
        memcpy(&new_array[new_start], &array[start], 
               count*sizeof(TMROctant)); 
      }
      else if (count > 0){
        // Send the element array to the new owner
        MPI_Isend(&array[start], count, TMROctant_MPI_type,
                  i, 0, comm, &send_requests[send_count]);
        send_count++;
      }
    }
  }

  // Match up the intervals for recvs
  for ( int i = 0; i < mpi_size; i++ ){
    // If we're on the correct interval, 
    if (i != mpi_rank && 
        ptr[i+1] >= new_ptr[mpi_rank] && 
        ptr[i] < new_ptr[mpi_rank+1]){
      // The start location within the local array
      int start = ptr[i] - new_ptr[mpi_rank];
      if (start < 0){ start = 0; }

      // Set the end location
      int max_count = new_ptr[mpi_rank+1] - new_ptr[mpi_rank];
      int end = ptr[i+1] - new_ptr[mpi_rank];
      if (end > max_count){ end = max_count; }

      // Set the count
      int count = end - start;

      // Send the element array to the new owner
      if (count > 0){
        MPI_Recv(&new_array[start], count, TMROctant_MPI_type,
                 i, 0, comm, MPI_STATUS_IGNORE);
      }
    }
  }

  // Free the memory
  delete [] ptr;
  delete [] new_ptr;

  // Wait for any remaining sends to complete
  MPI_Waitall(send_count, send_requests, MPI_STATUSES_IGNORE);
  delete [] send_requests;

  // Free the octant arrays
  delete octants;
  octants = new TMROctantArray(new_array, new_size);

  owners = new TMROctant[ mpi_size ];
  MPI_Allgather(&new_array[0], 1, TMROctant_MPI_type, 
                owners, 1, TMROctant_MPI_type, comm);
}

/*
  Duplicate the forest

  This function creates a duplicate representation of the current
  forest. This function copies the global connectivity of the forest
  and copies each individual tree.
*/
TMROctForest *TMROctForest::duplicate(){
  TMROctForest *dup = new TMROctForest(comm);
  if (block_conn){
    copyData(dup);

    // Copy the octants
    dup->octants = octants->duplicate();
    dup->owners = new TMROctant[ mpi_size ];
    memcpy(dup->owners, owners, sizeof(TMROctant)*mpi_size);
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
  if (block_conn){
    copyData(coarse);

    // Coarsen the octants in the array
    int size;
    TMROctant *array;
    octants->getArray(&array, &size);

    // Set the offset to be 2**d-1
    int offset = (1 << 3) - 1;

    // Create a new queue of octants
    TMROctantQueue *queue = new TMROctantQueue();
 
    // Scan through the list and add the coarse 0-child of each octant
    // to the new queue of octants
    for ( int i = 0; i < size; i++ ){
      if (array[i].level > 0){
        if (array[i].childId() == 0){
          TMROctant p;
          array[i].parent(&p);
          queue->push(&p);
        }
      }
      else {
        queue->push(&array[i]);
      }
    }

    // Create the coarse octants
    coarse->octants = queue->toArray();

    // Set the owner array
    coarse->octants->getArray(&array, &size);
    coarse->owners = new TMROctant[ mpi_size ];
    MPI_Allgather(&array[0], 1, TMROctant_MPI_type, 
                  coarse->owners, 1, TMROctant_MPI_type, comm);
  }

  return coarse;
}

/*
  Refine the octree mesh based on the input refinement levels
*/
void TMROctForest::refine( const int refinement[],
                           int min_level, int max_level ){
  // Adjust the min and max levels to ensure consistency
  if (min_level < 0){ min_level = 0; }
  if (max_level > TMR_MAX_LEVEL){ max_level = TMR_MAX_LEVEL; }

  // This is just a sanity check
  if (min_level > max_level){ min_level = max_level; }

  // Free memory if it has been allocated
  if (adjacent){ delete adjacent; }
  if (nodes){ delete nodes; } 
  if (dep_faces){ delete dep_faces; }
  if (dep_edges){ delete dep_edges; }
  if (X){ delete X; }
  adjacent = NULL;
  nodes = NULL;
  dep_faces = NULL;
  dep_edges = NULL;
  X = NULL;
  
  // Create a hash table for the refined octants and the octants
  // that are external (on other processors)
  TMROctantHash *hash = new TMROctantHash();
  TMROctantHash *ext_hash = new TMROctantHash();

  // Get the current array of octants
  int size;
  TMROctant *array;
  octants->getArray(&array, &size);

  if (refinement){
    for ( int i = 0; i < size; i++ ){
      if (refinement[i] == 0){
        // We know that this octant is locally owned
        hash->addOctant(&array[i]);
      }
      else if (refinement[i] < 0){
        // Coarsen this octant
        if (array[i].level > min_level){
          TMROctant q;
          array[i].getSibling(0, &q);
          q.level = q.level-1;
          if (mpi_rank == getOctantMPIOwner(&q)){
            hash->addOctant(&q);
          }
          else {
            ext_hash->addOctant(&q);
          }
        }
        else {
          // If it is already at the min level, just add it
          hash->addOctant(&array[i]);
        }
      }
      else if (refinement[i] > 0){
        // Refine this octant
        if (array[i].level < max_level){
          TMROctant q = array[i];
          q.level += 1;
          q.getSibling(0, &q);
          if (mpi_rank == getOctantMPIOwner(&q)){
            hash->addOctant(&q);
          }
          else {
            ext_hash->addOctant(&q);
          }
        }
        else {
          // If the octant is at the max level add it without
          // refinement
          hash->addOctant(&array[i]);
        }
      }
    }
  }
  else {
    // No refinement array is provided. Just go ahead and refine
    // everything...
    for ( int i = 0; i < size; i++ ){
      if (array[i].level < max_level){
        TMROctant q = array[i];
        q.level += 1;
        q.getSibling(0, &q);
        if (mpi_rank == getOctantMPIOwner(&q)){
          hash->addOctant(&q);
        }
        else {
          ext_hash->addOctant(&q);
        }
      }
    }
  }
  
  // Free the old octants class
  delete octants;

  // Sort the list of external octants
  TMROctantArray *list = ext_hash->toArray();
  list->sort();
  delete ext_hash;

  // Get the local list octants added from other processors
  // and add them to the local hash table
  TMROctantArray *local = distributeOctants(list);
  delete list;

  // Get the local octants and add them to the hash table
  local->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    hash->addOctant(&array[i]);
  }
  delete local;

  // Cover the hash table to a list and uniquely sort it
  octants = hash->toArray();
  octants->sort();

  delete hash;
}

/*
  Transform the node from a local coordinate system into the global
  node numbers 

  This transforms the given octant to the coordinate system of the
  lowest octant touching this node if it is on an octree boundary.

  input/output:
  oct:  the octant representing a node in the local coordinate system
*/
void TMROctForest::transformNode( TMROctant *oct ){
  // Get the maximum octant length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Check if this node lies on an octree boundary
  int fx0 = (oct->x == 0);
  int fy0 = (oct->y == 0);
  int fz0 = (oct->z == 0);
  int fx = (fx0 || oct->x == hmax);
  int fy = (fy0 || oct->y == hmax);
  int fz = (fz0 || oct->z == hmax);

  if (fx || fy || fz){
    // Get the original block index
    int block = oct->block;

    if (fx && fy && fz){
      // This node lies on a corner
      int corner = (fx0 ? 0 : 1) + (fy0 ? 0 : 2) + (fz0 ? 0 : 4);
      int node = block_conn[8*block + corner];

      // Transform the octant to each other octree frame
      // and check which processor owns it
      int owner = node_block_owners[node];

      if (block != owner){
        for ( int ip = node_block_ptr[node];
              ip < node_block_ptr[node+1]; ip++ ){
          int adjacent = node_block_conn[ip]/8;
          if (adjacent == owner){
            int adj_index = node_block_conn[ip] % 8;

            oct->block = adjacent;
            oct->x = hmax*(adj_index % 2);
            oct->y = hmax*((adj_index % 4)/2);
            oct->z = hmax*(adj_index/4);
            break;
          }
        }
      }
    }
    else if ((fy && fz) || (fx && fz) || (fx && fy)){
      // This node lies on an edge
      int edge_index = 0;
      int32_t u = 0;
      if (fy && fz){
        edge_index = (fy0 ? 0 : 1) + (fz0 ? 0 : 2);
        u = oct->x;
      }
      else if (fx && fz){
        edge_index = (fx0 ? 4 : 5) + (fz0 ? 0 : 2);
        u = oct->y;
      }
      else {
        edge_index = (fx0 ? 8 : 9) + (fy0 ? 0 : 2);
        u = oct->z;
      }

      // Get the edge number
      int edge = block_edge_conn[12*block + edge_index];

      // Get the edge owner
      int owner = edge_block_owners[edge];
      
      if (block != owner){
        // Retrieve the first and second node numbers to determine the
        // relative orientation between this edge and each adjacent edge
        int n1 = block_conn[8*block + block_to_edge_nodes[edge_index][0]];
        int n2 = block_conn[8*block + block_to_edge_nodes[edge_index][1]];

        for ( int ip = edge_block_ptr[edge]; 
              ip < edge_block_ptr[edge+1]; ip++ ){
          // Get the adjacent edge index on the opposite block
          int adjacent = edge_block_conn[ip]/12;

          if (owner == adjacent){
            int adj_index = edge_block_conn[ip] % 12;
      
            // Get the orientation
            int nn1 = block_conn[8*adjacent + block_to_edge_nodes[adj_index][0]];
            int nn2 = block_conn[8*adjacent + block_to_edge_nodes[adj_index][1]];
            
            // Determine whether the edges are in the same direction
            // or are reversed
            int reverse = (n1 == nn2 && n2 == nn1);
            
            // Set the u-coordinate along the edge
            int32_t uoct = u;
            if (reverse){
              uoct = hmax - u;
            }
          
            // Transform the octant to the adjacent coordinate system
            oct->block = adjacent;
            if (adj_index < 4){
              oct->x = uoct;
              oct->y = hmax*(adj_index % 2);
              oct->z = hmax*(adj_index/2);
            }
            else if (adj_index < 8){
              oct->x = hmax*(adj_index % 2);
              oct->y = uoct;
              oct->z = hmax*((adj_index-4)/2);
            }
            else {
              oct->x = hmax*(adj_index % 2);
              oct->y = hmax*((adj_index-8)/2);
              oct->z = uoct;
            }
            
            break;
          }
        }
      }
    }
    else {
      // Which face index are we dealing with?
      int face_index =
        fx*(fx0 ? 0 : 1) + fy*(fy0 ? 2 : 3) + fz*(fz0 ? 4 : 5);
      
      // Get the face owner
      int face = block_face_conn[6*block + face_index];

      // Get the face owner
      int owner = face_block_owners[face];

      if (block != owner){
        // Get the face id number
        int face_id = block_face_ids[6*block + face_index];
        
        // Get the u/v coordinates for this node on the owner face
        int32_t u, v;
        if (face_index < 2){ // x-face
          get_face_node_coords(face_id, oct->y, oct->z, &u, &v);
        }
        else if (face_index < 4){ // y-face
          get_face_node_coords(face_id, oct->x, oct->z, &u, &v);
        }
        else { // z-face
          get_face_node_coords(face_id, oct->x, oct->y, &u, &v);
        }

        // Loop over the blocks adjacent to this face to find
        // the block owner and its orientation
        for ( int ip = face_block_ptr[face]; 
              ip < face_block_ptr[face+1]; ip++ ){
        
          // Get the adjacent block index
          int adjacent = face_block_conn[ip]/6;
          if (adjacent == owner){
            int adj_index = face_block_conn[ip] % 6;

            // Get the face_id corresponding to the orientation of
            // this adjacent face
            face_id = block_face_ids[6*adjacent + adj_index];
          
            // Transform the octant p to the local octant coordinates
            oct->block = adjacent;
            if (adj_index < 2){
              oct->x = hmax*(adj_index % 2);
              set_face_node_coords(face_id, u, v, &oct->y, &oct->z);
            }
            else if (adj_index < 4){
              oct->y = hmax*(adj_index % 2);
              set_face_node_coords(face_id, u, v, &oct->x, &oct->z);
            }
            else {
              oct->z = hmax*(adj_index % 2);
              set_face_node_coords(face_id, u, v, &oct->x, &oct->y);
            }
          
            break;
          }
        }
      }
    }

    // Truncate the node back into the domain if it is on any of the
    // outer boundaries
    if (oct->x == hmax){ oct->x = hmax-1; }
    if (oct->y == hmax){ oct->y = hmax-1; }
    if (oct->z == hmax){ oct->z = hmax-1; }
  }
}

/*
  Get the owner of the octant
*/
int TMROctForest::getOctantMPIOwner( TMROctant *oct ){
  int rank = 0; 

  // while (owners[rank+1] <= oct) rank++
  for ( ; (rank < mpi_size-1 && 
           owners[rank+1].compareEncoding(oct) <= 0); rank++ );

  return rank;
}

/*
  Match the octant intervals
*/
void TMROctForest::matchOctantIntervals( TMROctant *array,
                                         int size, 
                                         int *ptr ){
  ptr[0] = 0;

  int index = 0;
  for ( int rank = 0; rank < mpi_size-1; rank++ ){
    while (index < size &&
           owners[rank+1].compareEncoding(&array[index]) > 0){
      index++;
    }
    ptr[rank+1] = index;
  }
  ptr[mpi_size] = size;
}

/*
  Match the MPI intervals
*/
void TMROctForest::matchMPIIntervals( TMROctant *array,
                                      int size, int *ptr ){
  ptr[0] = 0;
  for ( int i = 0, rank = 0; rank < mpi_size; rank++ ){
    while (i < size && array[i].tag <= rank){
      i++;
    }
    ptr[rank+1] = i;
  }
  ptr[mpi_size] = size;
}

/*
  Send a distributed list of octants to their owner processors
*/
TMROctantArray *TMROctForest::distributeOctants( TMROctantArray *list,
                                                 int use_tags,
                                                 int **_oct_ptr, 
                                                 int **_oct_recv_ptr,
                                                 int include_local ){
  // Get the array itself
  int size;
  TMROctant *array;
  list->getArray(&array, &size);

  // The number of octants that will be sent from this processor
  // to all other processors in the communicator
  int *oct_ptr = new int[ mpi_size+1 ];
  int *oct_recv_ptr = new int[ mpi_size+1 ];

  // Match the octant intervals to determine how mnay octants
  // need to be sent to each processor
  if (use_tags){
    matchMPIIntervals(array, size, oct_ptr);
  }
  else {
    matchOctantIntervals(array, size, oct_ptr);
  }

  // Count up the number of octants
  int *oct_counts = new int[ mpi_size ];
  for ( int i = 0; i < mpi_size; i++ ){
    if (!include_local && i == mpi_rank){
      oct_counts[i] = 0;
    }
    else {
      oct_counts[i] = oct_ptr[i+1] - oct_ptr[i];
    }
  }

  // Now distribute the octants to their destination octrees and
  // balance the corresponding octrees including the new elements.
  int *oct_recv_counts = new int[ mpi_size ];
  MPI_Alltoall(oct_counts, 1, MPI_INT,
               oct_recv_counts, 1, MPI_INT, comm);

  // Now use oct_ptr to point into the recv array
  oct_recv_ptr[0] = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    oct_recv_ptr[i+1] = oct_recv_ptr[i] + oct_recv_counts[i];
  }

  delete [] oct_counts;
  delete [] oct_recv_counts;

  // Create the distributed array
  TMROctantArray *dist = sendOctants(list, oct_ptr, oct_recv_ptr);

  // Free other data associated with the parallel communication
  if (_oct_ptr){ 
    *_oct_ptr = oct_ptr;
  }
  else {
    delete [] oct_ptr;
  }
  if (_oct_recv_ptr){
    *_oct_recv_ptr = oct_recv_ptr;
  }
  else {
    delete [] oct_recv_ptr;
  }

  return dist;
}

/*
  Send the octants to the processors designated by the pointer arrays
*/
TMROctantArray *TMROctForest::sendOctants( TMROctantArray *list,
                                           const int *oct_ptr,
                                           const int *oct_recv_ptr ){
  // Get the array itself
  int size;
  TMROctant *array;
  list->getArray(&array, &size);

  // Count up the number of recvs
  int nsends = 0, nrecvs = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank){
      if (oct_ptr[i+1] - oct_ptr[i] > 0){
        nsends++;
      }
      if (oct_recv_ptr[i+1] - oct_recv_ptr[i] > 0){
        nrecvs++;
      }
    }
  }

  // Allocate the space for the recving array
  int recv_size = oct_recv_ptr[mpi_size];
  TMROctant *recv_array = new TMROctant[ recv_size ];

  // Allocate space for the requests
  MPI_Request *send_request = new MPI_Request[ nsends ];

  // Loop over all the ranks and send 
  for ( int i = 0, j = 0; i < mpi_size; i++ ){
    if (i != mpi_rank && oct_ptr[i+1] - oct_ptr[i] > 0){
      // Post the send to the destination
      int count = oct_ptr[i+1] - oct_ptr[i];
      MPI_Isend(&array[oct_ptr[i]], count, TMROctant_MPI_type, 
                i, 0, comm, &send_request[j]);
      j++;
    }
    else if (i == mpi_rank){
      int count = oct_recv_ptr[i+1] - oct_recv_ptr[i];
      if (count > 0 &&
          (count == oct_ptr[i+1] - oct_ptr[i])){
        memcpy(&recv_array[oct_recv_ptr[i]], 
               &array[oct_ptr[i]], count*sizeof(TMROctant));
      }
    }
  }

  // Loop over the recieve calls
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank && oct_recv_ptr[i+1] > oct_recv_ptr[i]){
      int recv_count = oct_recv_ptr[i+1] - oct_recv_ptr[i];
     
      MPI_Recv(&recv_array[oct_recv_ptr[i]], recv_count, 
               TMROctant_MPI_type,
               i, 0, comm, MPI_STATUS_IGNORE);
    }
  }

  // Wait for any remaining sends to complete
  MPI_Waitall(nsends, send_request, MPI_STATUSES_IGNORE);
  delete [] send_request;

  return new TMROctantArray(recv_array, recv_size);
}

/*
  Add the face neighbors for an adjacent tree

  This function balances the given octant p across the face with the
  given face_index. The adjacent octants are added to the hash object.
  If the octants are not in the hash, they are also pushed to the queue
  for balancing.

  input:
  face_index:  index for the given face
  p:           octant to balance
  hash:        the array of hash objects for each processor
  queue:       the array of newly added qudrants for each processor
*/
void TMROctForest::addFaceNeighbors( int face_index, 
                                     TMROctant p,
                                     TMROctantHash *hash,
                                     TMROctantHash *ext_hash,
                                     TMROctantQueue *queue ){
  // Determine the global face number
  int block = p.block;
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
      // Get the face index on the adjacent block
      int adj_index = face_block_conn[ip] % 6;

      // Get the face id for the matching face
      face_id = block_face_ids[6*adjacent + adj_index];

      // Get the neighboring octant on the face
      TMROctant neighbor;
      neighbor.block = adjacent;
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
      
      // Find the octant owner and add the octant to the hash
      // tables and possibly queue
      int owner = getOctantMPIOwner(&neighbor);
      if (owner == mpi_rank){
        if (hash->addOctant(&neighbor)){
          queue->push(&neighbor);
        }
      }
      else if (ext_hash && ext_hash->addOctant(&neighbor)){
        queue->push(&neighbor);
      }
    }
  }
}

/*
  Add the edge neighbors for adjacent trees

  This function is called to balance the forest across tree edges.
  Given an octant p on the specified edge index, this code ensures a
  edge balanced tree, by adding the corresponding edge octants to all
  edge-adjacent octrees. If the octant is added to the hash object, it
  is then appended to the local face queue to ensure that it is also
  balanced locally.

  input:
  edge_index:  index for the given edge
  p:           octant to balance
  hash:        the array of hash objects for each processor
  queue:       the array of newly added qudrants for each processor
*/
void TMROctForest::addEdgeNeighbors( int edge_index, 
                                     TMROctant p,
                                     TMROctantHash *hash,
                                     TMROctantHash *ext_hash,
                                     TMROctantQueue *queue ){
  // First determine the global edge number
  int block = p.block;
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

  // Now, cycle through all the adjacent blocks
  for ( int ip = edge_block_ptr[edge]; ip < edge_block_ptr[edge+1]; ip++ ){
    // Get the block that is adjacent across this edge
    int adjacent = edge_block_conn[ip]/12;
    if (adjacent != block){
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
      neighbor.block = adjacent;
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

      // Find the octant owner and add the octant to the hash
      // tables and possibly queue
      int owner = getOctantMPIOwner(&neighbor);
      if (owner == mpi_rank){
        if (hash->addOctant(&neighbor)){
          queue->push(&neighbor);
        }
      }
      else if (ext_hash && ext_hash->addOctant(&neighbor)){
        queue->push(&neighbor);
      }
    }
  }
}

/*
  Add the corner neighbors for a given tree

  This function is called to balance the forest across tree corners.
  Given an octant p on the specified corner index, this code ensures a
  corner balanced tree, by adding the corresponding corner octants to
  all node-adjacent octrees. If the octant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  corner:  the corner index (p must lie on this corner)
  p:       the local octant
  hash:    the array of hash objects for each processor
  queue:   the array of newly added qudrants for each processor
*/
void TMROctForest::addCornerNeighbors( int corner,
                                       TMROctant p,
                                       TMROctantHash *hash,
                                       TMROctantHash *ext_hash,
                                       TMROctantQueue *queue ){
  // First determine the global edge number
  int block = p.block;
  int node = block_conn[8*block + corner];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = node_block_ptr[node];
        ip < node_block_ptr[node+1]; ip++ ){
      
    // Get the faces that are adjacent across this edge
    int adjacent = node_block_conn[ip]/8;
    if (adjacent != block){
      int adj_index = node_block_conn[ip] % 8;

      // Compute the octant location
      TMROctant neighbor;
      neighbor.block = p.block;
      neighbor.level = p.level;
      neighbor.x = (hmax - 2*h)*(adj_index % 2);
      neighbor.y = (hmax - 2*h)*((adj_index % 4)/2);
      neighbor.z = (hmax - 2*h)*(adj_index/4);
      
      // Find the octant owner and add the octant to the hash
      // tables and possibly queue
      int owner = getOctantMPIOwner(&neighbor);
      if (owner == mpi_rank){
        if (hash->addOctant(&neighbor)){
          queue->push(&neighbor);
        }
      }  
      else if (ext_hash && ext_hash->addOctant(&neighbor)){
        queue->push(&neighbor);
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
  ext_hash:        the array of hash tables for each block
  queue:           the octant queues for each block
  balance_corner:  balance across corners 
  balance_tree:    balance on the entire tree
*/
void TMROctForest::balanceOctant( TMROctant *oct,
                                  TMROctantHash *hash,
                                  TMROctantHash *ext_hash,
                                  TMROctantQueue *queue,
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
        // Get the MPI rank of the octant owner
        int owner = getOctantMPIOwner(&q);
        if (owner == mpi_rank){
          if (hash->addOctant(&q)){
            queue->push(&q);
          }
        }
        else if (ext_hash && ext_hash->addOctant(&q)){
          queue->push(&q);
        }
      }
      else if (balance_tree){
        addFaceNeighbors(face, q, hash, ext_hash, queue);
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
        // Get the MPI rank of the octant owner
        int owner = getOctantMPIOwner(&q);
        if (owner == mpi_rank){
          if (hash->addOctant(&q)){
            queue->push(&q);
          }
        }
        else if (ext_hash && ext_hash->addOctant(&q)){
          queue->push(&q);
        }
      }
      else if (balance_tree){
        // The node may lie across an edge or face
        int ex = (q.x < 0 || q.x >= hmax);
        int ey = (q.y < 0 || q.y >= hmax);
        int ez = (q.z < 0 || q.z >= hmax);
        
        if ((ex && ey) || (ex && ez) || (ey && ez)){
          // The octant lies along a true edge
          addEdgeNeighbors(edge, q, hash, ext_hash, queue);
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
          addFaceNeighbors(face, q, hash, ext_hash, queue);
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
          // Get the MPI rank of the octant owner
          int owner = getOctantMPIOwner(&q);
          if (owner == mpi_rank){
            if (hash->addOctant(&q)){
              queue->push(&q);
            }
          }
          else if (ext_hash && ext_hash->addOctant(&q)){
            queue->push(&q);
          }
        }
        else if (balance_tree){
          // The node may lie across a corner, edge or face
          int ex = (q.x < 0 || q.x >= hmax);
          int ey = (q.y < 0 || q.y >= hmax);
          int ez = (q.z < 0 || q.z >= hmax);

          if (ex && ey && ez){
            // Add the octant to the other trees
            addCornerNeighbors(corner, neighbor, hash, ext_hash, queue);
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
            addEdgeNeighbors(edge, q, hash, ext_hash, queue);
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
            addFaceNeighbors(face, q, hash, ext_hash, queue);
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
  // Create a hash table for the balanced tree
  TMROctantHash *hash = new TMROctantHash();
  TMROctantHash *ext_hash = new TMROctantHash();
  TMROctantQueue *queue = new TMROctantQueue();
  
  // Get the array of octants
  int oct_size;
  TMROctant *oct_array;
  octants->getArray(&oct_array, &oct_size);
    
  // Add all the elements
  for ( int i = 0; i < oct_size; i++ ){
    TMROctant oct;
    oct_array[i].getSibling(0, &oct);

    // Get the octant owner
    int owner = getOctantMPIOwner(&oct);

    // Add the owner
    if (owner == mpi_rank){
      hash->addOctant(&oct);
    }
    else {
      ext_hash->addOctant(&oct);
    }
      
    // Balance the octants locally
    const int balance_tree = 1;
    balanceOctant(&oct, hash, ext_hash, queue, 
                  balance_corner, balance_tree);
  }

  // Free the original octant array and set it to NULL
  delete octants;

  while (queue->length() > 0){
    // Now continue until the queue of added octants is
    // empty. At each iteration, pop an octant and add 
    // its neighbours until nothing new is added. This code
    // handles the propagation of octants to adjacent octants.
    TMROctant oct = queue->pop();
    const int balance_tree = 1;
    balanceOctant(&oct, hash, ext_hash, queue, 
                  balance_corner, balance_tree);
  }

  // Now everything is locally balanced - all the elements on the
  // current processor are balanced with all the other elements on the
  // current processor, but nothing is inter-processor balanced yet.
  // Create a sorted list of local the 0-child octants. This can be
  // further reduced to limit the amount of memory passed between
  // processors
  TMROctantArray *elems0 = ext_hash->toArray();
  elems0->sort();

  // Get the array of 0-octants 
  int size;
  TMROctant *array;
  elems0->getArray(&array, &size);

  if (size > 0){
    // Get the parent of the octant
    TMROctant p, s = array[0];
    s.parent(&p);

    // Loop over all of the local octants 
    for ( int i = 0; i < size; i++ ){
      if (!p.contains(&array[i])){
        queue->push(&s);
      }
      // Get the next octant and find its parent
      s = array[i];
      s.parent(&p);
    }
    
    // Push the last octant onto the queue
    queue->push(&s);
  }
  
  // Free the elements and the hash table
  delete elems0;

  // Send the octants to their destination...
  TMROctantArray *list = queue->toArray();
  delete queue;

  // Get the local octants from the list
  TMROctantArray *local = distributeOctants(list);
  delete list;

  // Allocate a new queue
  queue = new TMROctantQueue();

  // Get the local array of octants and add them to the
  // hash table
  local->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    if (hash->addOctant(&array[i])){
      queue->push(&array[i]);
    }
  }

  // Now all the received octants will balance the tree locally
  // without having to worry about off-processor octants.
  while (queue->length() > 0){
    const int balance_tree = 1;
    TMROctant oct = queue->pop();
    balanceOctant(&oct, hash, NULL, queue, 
                  balance_corner, balance_tree);
  }

  // Now convert the elements from child-0 elements to
  // elements which cover the full mesh
  TMROctantArray *child0_elems = hash->toArray();
  child0_elems->getArray(&oct_array, &oct_size);

  // Loop over all elements and add their siblings
  for ( int i = 0; i < oct_size; i++ ){
    if (oct_array[i].level > 0){
      for ( int j = 0; j < 8; j++ ){
        TMROctant q;
        oct_array[i].getSibling(j, &q);
        int owner = getOctantMPIOwner(&q);
        if (mpi_rank == owner){
          hash->addOctant(&q);
        }
        else {
          queue->push(&q);
        }
      }
    }
  }

  // Free the temporary elements
  delete child0_elems;

  // Turn the queue into an array
  list = queue->toArray();
  delete queue;

  // Sort the list before distributing it
  list->sort();

  // Get the local list octants added from other processors
  // and add them to the local hash table
  local = distributeOctants(list);
  delete list;

  // Get the local octants and add them to the hash table
  local->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    hash->addOctant(&array[i]);
  }
  delete local;
  
  // Set the elements into the octree
  octants = hash->toArray();
  octants->sort();

  // Free the hash
  delete hash;
}

/*
  Add the octant to the processor queues corresponding to the
  non-local blocks that touch the given face
*/
void TMROctForest::addAdjacentFaceToQueue( int face_index, 
                                           TMROctant p,
                                           TMROctantQueue *queue, 
                                           TMROctant orig ){
  // Determine the global face number
  int block = p.block;
  int face = block_face_conn[6*block + face_index];
  int face_id = block_face_ids[6*block + face_index];
 
  // Set the size of the side-length of the octant
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Get the u/v coordinates of this face
  int32_t u, v;
  if (face_index < 2){ // x-face
    get_face_oct_coords(face_id, h, p.y, p.z, &u, &v);
  }
  else if (face_index < 4){ // y-face
    get_face_oct_coords(face_id, h, p.x, p.z, &u, &v);
  }
  else { // z-face
    get_face_oct_coords(face_id, h, p.x, p.y, &u, &v);
  }

  // Loop over all the adjacent faces and add the block
  for ( int ip = face_block_ptr[face]; ip < face_block_ptr[face+1]; ip++ ){
    // Get the adjacent block and the face index on the adjacent block
    int adjacent = face_block_conn[ip]/6;
    if (adjacent != block){
      // Get the face index on the adjacent block
      int adj_index = face_block_conn[ip] % 6;

      // Get the face id for the matching face
      face_id = block_face_ids[6*adjacent + adj_index];

      // Get the neighboring octant on the face
      TMROctant neighbor;
      neighbor.block = adjacent;
      neighbor.level = p.level;
      if (adj_index < 2){
        neighbor.x = (hmax - h)*(adj_index % 2);
        set_face_oct_coords(face_id, h, u, v, 
                            &neighbor.y, &neighbor.z);
      }
      else if (adj_index < 4){
        neighbor.y = (hmax - h)*(adj_index % 2);
        set_face_oct_coords(face_id, h, u, v, 
                            &neighbor.x, &neighbor.z);
      }
      else {
        neighbor.z = (hmax - h)*(adj_index % 2);
        set_face_oct_coords(face_id, h, u, v, 
                            &neighbor.x, &neighbor.y);
      }
      
      // Find the octant owner
      int owner = getOctantMPIOwner(&neighbor);
      if (owner != mpi_rank){
        orig.tag = owner;
        queue->push(&orig);
      }
    }
  }
}

/*
  Add the octant to the processor queue corresponding to the non-local
  blocks that touch the given edge 
*/
void TMROctForest::addAdjacentEdgeToQueue( int edge_index, 
                                           TMROctant p,
                                           TMROctantQueue *queue, 
                                           TMROctant orig ){
  // First determine the global edge number
  int block = p.block;
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
      // Get the adjacent edge index
      int adj_index = edge_block_conn[ip] % 12;

      // Get the nodes on the adjacent block
      int nn1 = block_conn[8*adjacent + block_to_edge_nodes[adj_index][0]];
      int nn2 = block_conn[8*adjacent + block_to_edge_nodes[adj_index][1]];

      // Add the octant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - h - ucoord;
      }
      
      TMROctant neighbor;
      neighbor.block = adjacent;
      neighbor.level = p.level;
      if (adj_index < 4){
        neighbor.x = u;
        neighbor.y = (hmax - h)*(adj_index % 2);
        neighbor.z = (hmax - h)*(adj_index/2);
      }
      else if (adj_index < 8){
        neighbor.x = (hmax - h)*(adj_index % 2);
        neighbor.y = u;
        neighbor.z = (hmax - h)*((adj_index-4)/2);
      }
      else {
        neighbor.x = (hmax - h)*(adj_index % 2);
        neighbor.y = (hmax - h)*((adj_index-8)/2);
        neighbor.z = u;
      }
      
      // Find the octant owner
      int owner = getOctantMPIOwner(&neighbor);
      if (owner != mpi_rank){
        orig.tag = owner;
        queue->push(&orig);
      }
    }
  }
}

/*
  Add the octant to the queue that correspond to the non-local blocks
  that touch the corner 
*/
void TMROctForest::addAdjacentCornerToQueue( int corner,
                                             TMROctant p,
                                             TMROctantQueue *queue, 
                                             TMROctant orig ){
  // First determine the global edge number
  int block = p.block;
  int node = block_conn[8*block + corner];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = node_block_ptr[node];
        ip < node_block_ptr[node+1]; ip++ ){
      
    // Get the faces that are adjacent across this edge
    int adjacent = node_block_conn[ip]/8;
    if (adjacent != block){
      int adj_index = node_block_conn[ip] % 8;

      // Compute the octant location
      TMROctant neighbor;
      neighbor.block = adjacent;
      neighbor.level = p.level;
      neighbor.x = (hmax - h)*(adj_index % 2);
      neighbor.y = (hmax - h)*((adj_index % 4)/2);
      neighbor.z = (hmax - h)*(adj_index/4);
      
      // Find the octant owner
      int owner = getOctantMPIOwner(&neighbor);
      if (owner != mpi_rank){
        orig.tag = owner;
        queue->push(&orig);
      }
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
void TMROctForest::computeAdjacentOctants(){
  if (adjacent){ delete adjacent; }

  // Allocate the queue that stores the octants destined for each of
  // the processors
  TMROctantQueue *queue = new TMROctantQueue();
  
  // Get the actual octant array
  int size;
  TMROctant *array;
  octants->getArray(&array, &size);

  // Sibling ids adjacent to each local face index
  const int face_ids[][4] = {{0,2,4,6}, {1,3,5,7},
                             {0,1,4,5}, {2,3,6,7},
                             {0,1,2,3}, {4,5,6,7}};

  // Sibling ids adjacet to each local edge index
  const int edge_ids[][2] =
    {{0,1}, {2,3}, {4,5}, {6,7},
     {0,2}, {1,3}, {4,6}, {5,7},
     {0,4}, {1,5}, {2,6}, {3,7}};
  
  // Loop over all the elements and check where we need to send
  // the octants that are along each edge/face
  for ( int i = 0; i < size; i++ ){
    const int32_t hmax = 1 << TMR_MAX_LEVEL;
    const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
    
    // Look across each face
    for ( int face_index = 0; face_index < 6; face_index++ ){
      for ( int k = 0; k < 4; k++ ){
        // Get the neighbor across each face and check if the
        // neighboring octant is locally owned. If not, add array[i]
        // to queue destined for the adjacent processor.
        TMROctant q = array[i];
        q.level += 1;
        q.getSibling(face_ids[face_index][k], &q);
        q.faceNeighbor(face_index, &q);

        if ((q.x >= 0 && q.x < hmax) &&
            (q.y >= 0 && q.y < hmax) &&
            (q.z >= 0 && q.z < hmax)){
          // Get the MPI rank of the octant owner
          int owner = getOctantMPIOwner(&q);
          if (owner != mpi_rank){
            TMROctant p = array[i];
            p.tag = owner;
            queue->push(&p);
          }
        }
        else {
          // Transform the octant across the boundary
          // and perform the same check on each
          addAdjacentFaceToQueue(face_index, q, queue, array[i]);
        }
      }
    }

    // Add the edge-adjacent octant across the boundary
    for ( int edge_index = 0; edge_index < 12; edge_index++ ){
      for ( int k = 0; k < 2; k++ ){
        // Get the next-lowest level octant
        TMROctant q = array[i];
        q.level += 1;
        q.getSibling(edge_ids[edge_index][k], &q);
        q.edgeNeighbor(edge_index, &q);
	  
        // If we're in bounds, add the neighbor
        if ((q.x >= 0 && q.x < hmax) &&
            (q.y >= 0 && q.y < hmax) &&
            (q.z >= 0 && q.z < hmax)){
          // Get the MPI rank of the octant owner
          int owner = getOctantMPIOwner(&q);
          if (owner != mpi_rank){
            TMROctant p = array[i];
            p.tag = owner;
            queue->push(&p);
          }
        }
        else {
          // The node may lie across an edge or face
          int ex = (q.x < 0 || q.x >= hmax);
          int ey = (q.y < 0 || q.y >= hmax);
          int ez = (q.z < 0 || q.z >= hmax);
          
          if ((ex && ey) || (ex && ez) || (ey && ez)){
            // The octant lies along a true edge
            addAdjacentEdgeToQueue(edge_index, q, queue, array[i]);
        }
          else {
            // The octant actually lies along a face 
            int face_index = 0;
            if (ex){ // This is an x-face
              face_index = (q.x < 0 ? 0 : 1);
            }
            else if (ey){ // This is a y-face
              face_index = (q.y < 0 ? 2 : 3);
            }
            else { // This is a z-face
              face_index = (q.z < 0 ? 4 : 5);
            }
            addAdjacentFaceToQueue(face_index, q, queue, array[i]);
          }
        }
      }
    }

    // Add corner-adjacent octants
    for ( int corner = 0; corner < 8; corner++ ){
      TMROctant q = array[i];
      q.level += 1;
      q.getSibling(corner, &q);
      q.cornerNeighbor(corner, &q);
        
      if ((q.x >= 0 && q.x < hmax) &&
          (q.y >= 0 && q.y < hmax) &&
          (q.z >= 0 && q.z < hmax)){
        // Get the MPI rank of the octant owner
        int owner = getOctantMPIOwner(&q);
        if (owner != mpi_rank){
          TMROctant p = array[i];
          p.tag = owner;
          queue->push(&p);
        }
      }
      else {
        // The node may lie across a corner, edge or face
        int ex = (q.x < 0 || q.x >= hmax);
        int ey = (q.y < 0 || q.y >= hmax);
        int ez = (q.z < 0 || q.z >= hmax);
        
        if (ex && ey && ez){
          // Add the octant to the other trees
          addAdjacentCornerToQueue(corner, q, queue, array[i]);
        }
        else if ((ex && ey) || (ex && ez) || (ey && ez)){
          // The octant lies along a true edge
          int edge_index = 0;
          if (ey && ez){
            edge_index = (q.y < 0 ? 0 : 1) + (q.z < 0 ? 0 : 2);
          }
          else if (ex && ez){
            edge_index = (q.x < 0 ? 4 : 5) + (q.z < 0 ? 0 : 2);
          }
          else { // (ex && ey)
            edge_index = (q.x < 0 ? 8 : 9) + (q.y < 0 ? 0 : 2);
          }
          addAdjacentEdgeToQueue(edge_index, q, queue, array[i]);
        }
        else {
          // The octant actually lies along a face 
          int face_index = 0;
          if (ex){ // This is an x-face
            face_index = (q.x < 0 ? 0 : 1);
          }
          else if (ey){ // This is a y-face
            face_index = (q.y < 0 ? 2 : 3);
          }
          else { // This is a z-face
            face_index = (q.z < 0 ? 4 : 5);
          }
          addAdjacentFaceToQueue(face_index, q, queue, array[i]);
        }
      }
    }
  }

  // Convert the local adjacency non-local list of octants
  TMROctantArray *list = queue->toArray();
  list->getArray(&array, &size);
  qsort(array, size, sizeof(TMROctant), compare_octant_tags);

  // Distribute the octants
  int use_tags = 1;
  adjacent = distributeOctants(list, use_tags);
  adjacent->sort();

  delete list;
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
int TMROctForest::checkAdjacentDepFaces( int face_index,
                                         TMROctant *b,
                                         TMROctantArray *adjocts ){
  // Get the side length of the octant
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);

  // Get the face id number
  int block_owner = b->block;
  int face = block_face_conn[6*block_owner + face_index];
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
      oct.block = block;
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
      
      // If the more-refined element exists then label the
      // corresponding nodes as dependent
      const int use_node_search = 0;
      if (octants->contains(&oct, use_node_search) || 
          (adjocts && adjocts->contains(&oct, use_node_search))){
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
  edge_index:    the local edge index
  b:             the quadrant
*/
int TMROctForest::checkAdjacentDepEdges( int edge_index,
                                         TMROctant *b,
                                         TMROctantArray *adjocts ){
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
  int block_owner = b->block;
  int edge = block_edge_conn[12*block_owner + edge_index];
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
      oct.block = block;
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
      
      // If the more-refined element exists then label the
      // corresponding nodes as dependent
      const int use_node_search = 0;
      if (octants->contains(&oct, use_node_search) || 
          (adjocts && adjocts->contains(&oct, use_node_search))){
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
    delete dep_edges;
    dep_edges = NULL;
  }
  if (dep_faces){
    delete dep_faces;
    dep_faces = NULL;
  }
  if (dep_ptr){ 
    delete [] dep_ptr;  
    delete [] dep_conn;
    delete [] dep_weights; 
    dep_ptr = NULL;  
    dep_conn = NULL; 
    dep_weights = NULL; 
  }

  // Create the queue that will store the dependent faces/edges
  TMROctantQueue *dfaces = new TMROctantQueue();
  TMROctantQueue *dedges = new TMROctantQueue();
            
  for ( int iter = 0; iter < 2; iter++ ){
    // Get the elements either in the regular octant array or
    // in the adjacent element array
    int size = 0;
    TMROctant *array = NULL;
    TMROctantArray *adjocts = NULL;
    if (iter == 0){
      octants->getArray(&array, &size);
      adjocts = adjacent;
    }
    else if (adjacent){
      adjacent->getArray(&array, &size);
      adjocts = NULL;
    }

    // Face sibling ids for each face index
    const int face_ids[][4] = {{0,2,4,6}, {1,3,5,7},
                               {0,1,4,5}, {2,3,6,7},
                               {0,1,2,3}, {4,5,6,7}};

    for ( int i = 0; i < size; i++ ){
      // Get the side length of the element
      const int32_t hmax = 1 << TMR_MAX_LEVEL;
    
      // Check the adjacent elements along each face
      for ( int face_index = 0; face_index < 6; face_index++ ){
        // Scan through the adjacent octants
        int add_me = 0;

        for ( int k = 0; k < 4; k++ ){
          // Get the octant and increase the level
          TMROctant p = array[i];
          p.level += 1;

          // Get the sibling id for each octant along the
          // face that we're on right now
          TMROctant q;
          p.getSibling(face_ids[face_index][k], &q);

          // Get the face neighbor
          q.faceNeighbor(face_index, &q);
      
          // Check if this element is across a face
          int fx = (q.x < 0 || q.x >= hmax);
          int fy = (q.y < 0 || q.y >= hmax);
          int fz = (q.z < 0 || q.z >= hmax);
          
          // If this face lies across an octree face, then
          // check the adjacent octree
          if (fx || fy || fz){
            if (checkAdjacentDepFaces(face_index, &q, adjocts)){
              add_me = 1;
              break;
            }
          }
          else {
            // If the more-refined element exists then label the
            // corresponding nodes as dependent
            const int use_node_search = 0;
            if (octants->contains(&q, use_node_search) || 
                (adjocts && adjocts->contains(&q, use_node_search))){
              add_me = 1;
              break;
            }
          }
        }

        if (add_me){
          // Add the dependent face
          TMROctant t = array[i];
          t.tag = face_index;
          dfaces->push(&t);
        }
      }

      // Enumerate the sibling-ids for each edge. Along each edge
      // there are two siblings that touch the edge. 
      //
      //       /----7----/
      //      3|        2|
      //    /--|-6----/  |
      //    | 11      |  9
      //   10  |      8  |
      //    |  /----5-|--/
      //    | 1       | 0
      //    /----4----/   / x
      //                 /
      //        y ------/ 

      // Edge->sibling id data...
      const int edge_ids[][2] =
        {{0,1}, {2,3}, {4,5}, {6,7},
         {0,2}, {1,3}, {4,6}, {5,7},
         {0,4}, {1,5}, {2,6}, {3,7}};

      // Check whether the next-level refined element exists over an
      // adjacent edge
      for ( int edge_index = 0; edge_index < 12; edge_index++ ){
        // Keep a flag to check if we need to add the edge as a
        // dependent or not
        int add_me = 0;

        for ( int k = 0; k < 2; k++ ){
          // Get the octant and increase the level
          TMROctant p = array[i];
          p.level += 1;

          // Get the sibling id for each octant along the
          // face that we're on right now
          TMROctant q;
          p.getSibling(edge_ids[edge_index][k], &q);

          // Get the edge neighbor
          q.edgeNeighbor(edge_index, &q);

          // Check if the adjacent octant q lies over an octree edge
          // or face and if so check for the corresponding octant
          int fx0 = (q.x < 0);
          int fy0 = (q.y < 0);
          int fz0 = (q.z < 0);
          int fx = (fx0 || q.x >= hmax);
          int fy = (fy0 || q.y >= hmax);
          int fz = (fz0 || q.z >= hmax);
          
          if ((fx && fy) || (fy && fz) || (fx && fz)){
            if (checkAdjacentDepEdges(edge_index, &q, adjocts)){
              add_me = 1; break;
            }
          }
          else if (fx || fy || fz){
            int face_index =
              fx*(fx0 ? 0 : 1) + fy*(fy0 ? 2 : 3) + fz*(fz0 ? 4 : 5);
            if (checkAdjacentDepFaces(face_index, &q, adjocts)){
              add_me = 1; break;
            }
          }
          else {
            const int use_node_search = 0;
            if (octants->contains(&q, use_node_search) || 
                (adjocts && adjocts->contains(&q, use_node_search))){
              add_me = 1; break;
            }
          }
        }
        
        if (add_me){
          TMROctant t = array[i];
          t.tag = edge_index;
          dedges->push(&t);
        }
      }
    }
  }

  // Create the arrays of the dependent nodes/faces
  dep_faces = dfaces->toArray();
  dep_edges = dedges->toArray();

  // Free the data
  delete dfaces;
  delete dedges;
}

/*
  Label the dependent face and edge nodes

  This code is called after all the dependent faces have been
  computed.  Note that this relies on the mesh being edge-balanced
  (which is required).
*/
void TMROctForest::labelDependentNodes(){
  // Get the array of dependent faces
  int dep_size;
  TMROctant *dep_array;
  dep_faces->getArray(&dep_array, &dep_size);

  // Loop over the labeled dependent faces
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
          node.block = b->block;
          node.level = 0;
          node.tag = -1;
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

          // Transform the node to the global ordering
          transformNode(&node);

          // Search for dependent node and label it
          const int use_node_search = 1;
          TMROctant *t = nodes->contains(&node, use_node_search);
          if (t){
            t->tag = -1;
          }
        }
      }
    }
  }

  // Get the array of dependent edges
  dep_edges->getArray(&dep_array, &dep_size);

  for ( int i = 0; i < dep_size; i++ ){
    TMROctant *b = &dep_array[i];
    
    // Find the edge lengths of the element octant and
    // the node spacing on the dependent face
    const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);
    const int32_t hc = 1 << (TMR_MAX_LEVEL - b->level - (mesh_order-1));
    
    // Get the edge index
    int edge_index = b->tag;

    // Loop over all the nodes on this edge
    for ( int ii = 1; ii < 2*mesh_order-1; ii += 2 ){
      // Label only the dependent nodes
      TMROctant node;
      node.block = b->block;
      node.level = 0;
      node.tag = -1;
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

      // Transform the node to the global ordering
      transformNode(&node);
      
      // Search for dependent node and label it
      const int use_node_search = 1;
      TMROctant *t = nodes->contains(&node, use_node_search);
      if (t){ 
        t->tag = -1;
      }
    }
  }
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

  // Send/recv the adjacent octants
  computeAdjacentOctants();

  // Compute the dependent face nodes
  computeDepFacesAndEdges();
  
  // Free the node data if it does not already exist
  if (nodes){ delete nodes; }

  // Get the current array of octants
  int size;
  TMROctant *array;
  octants->getArray(&array, &size);

  // Allocate an array large enough to store all of the nodes
  int index = 0;
  int max_nodes = mesh_order*mesh_order*mesh_order*size;
  TMROctant *all_nodes = new TMROctant[ max_nodes ];

  // Loop over all of the current octants and add the nodes
  for ( int i = 0; i < size; i++ ){
    if (mesh_order == 2){
      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);

      // Add all of the nodes to the hash
      for ( int kk = 0; kk < 2; kk++ ){
        for ( int jj = 0; jj < 2; jj++ ){
          for ( int ii = 0; ii < 2; ii++ ){
            // Set the node level to the highest level - this will
            // be updated when the nodes are assigned to the elements
            all_nodes[index].block = array[i].block;
            all_nodes[index].level = 0;
            
            // Set a positive tag, this will be replaced with a 
            // negative tag if the node is dependent
            all_nodes[index].tag = 1;

            // Set the node coordinates
            all_nodes[index].x = array[i].x + ii*h;
            all_nodes[index].y = array[i].y + jj*h;
            all_nodes[index].z = array[i].z + kk*h;

            // Transform the node to the global frame
            transformNode(&all_nodes[index]);
            index++;
          }
        }
      }
    }
    else if (mesh_order == 3){
      // Add all of the nodes from the adjacent elements
      const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level - 1);

      // Add all of the nodes to the hash
      for ( int kk = 0; kk < 3; kk++ ){
        for ( int jj = 0; jj < 3; jj++ ){
          for ( int ii = 0; ii < 3; ii++ ){
            // Set the node level to the highest level - this will
            // be updated when the nodes are assigned to the elements
            all_nodes[index].block = array[i].block;
            all_nodes[index].level = 0;
            
            // Set a positive tag, this will be replaced with a 
            // negative tag if the node is dependent
            all_nodes[index].tag = 1;

            // Set the node coordinates
            all_nodes[index].x = array[i].x + ii*h;
            all_nodes[index].y = array[i].y + jj*h;
            all_nodes[index].z = array[i].z + kk*h;

            // Transform the node to the global frame
            transformNode(&all_nodes[index]);
            index++;
          }
        }
      }
    }
  }

  // Create an array of all the octants and uniquely sort it
  nodes = new TMROctantArray(all_nodes, index);
  nodes->sort();

  // Count up all the locally owned, global variables and the number
  // of dependent variables
  int nlocal = 0;

  // Label the dependent nodes
  labelDependentNodes();

  // Extract the array of nodes
  nodes->getArray(&array, &size);

  // Count up the number of nodes that are owned by this processor
  // and append the nodes that are not owned by this processor to
  // a series of queues destined for other processors.
  for ( int i = 0; i < size; i++ ){
    // Get the owner of this node
    int owner = getOctantMPIOwner(&array[i]);

    // If this is the owner, then it is either dependent
    if (owner == mpi_rank && array[i].tag >= 0){
      nlocal++;
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

  // Count up the number of interior nodes
  for ( int i = 0; i < size; i++ ){
    // Get the owner of this node
    TMROctant p = array[i];
    int owner = getOctantMPIOwner(&p);
    
    // If this is the owner, set the node number
    if (owner == mpi_rank && array[i].tag >= 0){
      array[i].tag = node_num;
      node_num++;
    }
  }

  // Distribute the nodes
  int use_tags = 0;
  int *send_ptr, *recv_ptr;
  TMROctantArray *dist = distributeOctants(nodes, use_tags,
                                           &send_ptr, &recv_ptr);

  // Loop over the off-processor nodes and search for them in the
  // sorted node list and assign them the correct tag
  dist->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    if (array[i].tag >= 0){
      const int use_node_search = 1;
      TMROctant *t = nodes->contains(&array[i], use_node_search);
      array[i].tag = t->tag;
    }
  }

  // Adjust the pointer size
  int local_offset = send_ptr[mpi_rank+1] - send_ptr[mpi_rank]; 
  for ( int i = mpi_rank+1; i <= mpi_size; i++ ){
    send_ptr[i] -= local_offset;
  }

  // Send the nodes back to the original processors
  TMROctantArray *ext_nodes = sendOctants(dist, recv_ptr, send_ptr);
  delete dist;
  delete [] recv_ptr;

  // These are now the external octants with the correct
  // node numbers. Copy the node numbers into the local array
  nodes->getArray(&array, &size);
  int ext_size;
  TMROctant *ext_array;
  ext_nodes->getArray(&ext_array, &ext_size);
  int j = 0;
  for ( int i = 0; i < send_ptr[mpi_rank]; i++, j++ ){
    if (array[i].tag >= 0){
      array[i].tag = ext_array[j].tag;
    }
  }

  for ( int i = send_ptr[mpi_rank+1], ii = send_ptr[mpi_rank+1] + local_offset;
        i < send_ptr[mpi_size]; i++, j++, ii++ ){
    if (array[ii].tag >= 0){
      array[ii].tag = ext_array[j].tag;
    }
  }
  delete [] send_ptr;
  delete ext_nodes;
  
  // Count up the dependent nodes
  num_dep_nodes = 0;
  for ( int i = 0; i < size; i++ ){
    if (array[i].tag < 0){
      num_dep_nodes++;
      array[i].tag = -num_dep_nodes;
    }
  }

  // Allocate the positions for all of the nodes within the mesh
  nodes->getArray(&array, &size);
  X = new TMRPoint[ size ];
  memset(X, 0, size*sizeof(TMRPoint));

  const int32_t hmax = 1<< TMR_MAX_LEVEL;
  if (topo){    
    for ( int i = 0; i < size; i++ ){
      double u = 0.0, v = 0.0, w = 0.0;
      if (array[i].x == hmax-1){ 
        u = 1.0;
      }
      else {
        u = 1.0*array[i].x/hmax;
      }
      if (array[i].y == hmax-1){
        v = 1.0;
      }
      else {
        v = 1.0*array[i].y/hmax;
      }
      if (array[i].z == hmax-1){
        w = 1.0;
      }
      else {
        w = 1.0*array[i].z/hmax;
      }

      TMRVolume *volume;
      topo->getVolume(array[i].block, &volume);
      if (volume){
        volume->evalPoint(u, v, w, &X[i]);
      }
    }
  }
}

/*
  Get the elements that either lie in a volume, on a face or on a
  curve with a given attribute.

  This code loops over all octants that are locally owned on this
  processor and checks the attribute of the underlying topological
  entry. If the volume attribute matches, the octant is added
  directly, otherwise the local face or volume index is set as the
  tag. 

  input:
  attr:   string attribute associated with the geometric feature

  returns:
  list:   an array of octants satisfying the attribute
*/
TMROctantArray* TMROctForest::getOctsWithAttribute( const char *attr ){
  if (!topo){
    return NULL;
  }

  // Create a queue to store the elements that we find
  TMROctantQueue *queue = new TMROctantQueue();

  // Get the quadrants
  int size;
  TMROctant *array;
  octants->getArray(&array, &size);

  // Loop over the quadrants and find out whether it touches
  // a face or edge with the prescribed attribute
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  for ( int i = 0; i < size; i++ ){
    const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);

    // Get the surface quadrant
    TMRVolume *vol;
    topo->getVolume(array[i].block, &vol);  
    const char *vol_attr = vol->getAttribute();
    if (vol_attr && strcmp(vol_attr, attr) == 0){
      queue->push(&array[i]);
    }
    else {
      // Check if this octant lies on a face
      int fx0 = (array[i].x == 0);
      int fy0 = (array[i].y == 0);
      int fz0 = (array[i].z == 0);
      int fx = (fx0 || array[i].x + h == hmax);
      int fy = (fy0 || array[i].y + h == hmax);
      int fz = (fz0 || array[i].z + h == hmax);

      if (fx || fy || fz){
        TMRFace *face = NULL;
        TMROctant oct = array[i];
        if (fx){
          int face_index = (fx0 ? 0 : 1);
          int face_num = block_face_conn[6*array[i].block + face_index];
          topo->getFace(face_num, &face);
          const char *face_attr = face->getAttribute();
          if (face_attr && strcmp(face_attr, attr) == 0){
            oct.tag = face_index;
            queue->push(&oct);
          }
        }
        if (fy){
          int face_index = (fy0 ? 2 : 3);
          int face_num = block_face_conn[6*array[i].block + face_index];
          topo->getFace(face_num, &face);
          const char *face_attr = face->getAttribute();
          if (face_attr && strcmp(face_attr, attr) == 0){
            oct.tag = face_index;
            queue->push(&oct);
          }
        }
        if (fz){
          int face_index = (fz0 ? 4 : 5);
          int face_num = block_face_conn[6*array[i].block + face_index];
          topo->getFace(face_num, &face);
          const char *face_attr = face->getAttribute();
          if (face_attr && strcmp(face_attr, attr) == 0){
            oct.tag = face_index;
            queue->push(&oct);
          }
        }
      }
    }
  }

  TMROctantArray *list = queue->toArray();
  delete queue;
  return list;
}

/*
  Create an array of the nodes that are lie on a surface, edge or
  corner with a given attribute

  This code loops over all nodes and check whether they lie on a
  geometric entity that has the given attribute.

  input:
  attr:   the string of the attribute to search

  returns:
  list:   the nodes matching the specified attribute
*/
TMROctantArray* TMROctForest::getNodesWithAttribute( const char *attr ){
  if (!topo){
    return NULL;
  }

  // Create a queue to store the nodes that we find
  TMROctantQueue *queue = new TMROctantQueue();

  // Get the local nodal quadrants
  int size;
  TMROctant *array;
  nodes->getArray(&array, &size);

  // Loop over the quadrants and find out whether it touches
  // a face or edge with the prescribed attribute
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  for ( int i = 0; i < size; i++ ){
    // Check if this node lies on an octree boundary
    int fx0 = (array[i].x == 0);
    int fy0 = (array[i].y == 0);
    int fz0 = (array[i].z == 0);
    int fx = (fx0 || array[i].x == hmax);
    int fy = (fy0 || array[i].y == hmax);
    int fz = (fz0 || array[i].z == hmax);

    if (fx && fy && fz){
      // This node lies on a corner
      TMRVertex *vert;
      int vert_index = (fx0 ? 0 : 1) + (fy0 ? 0 : 2) + (fz0 ? 0 : 4);
      int vert_num = block_conn[8*array[i].block + vert_index];
      topo->getVertex(vert_num, &vert);
      const char *vert_attr = vert->getAttribute();
      if (vert_attr && strcmp(vert_attr, attr) == 0){
        queue->push(&array[i]);
      }
    }
    else if ((fy && fz) || (fx && fz) || (fx && fy)){
      // This node lies on an edge
      TMREdge *edge;
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
      int edge_num = block_edge_conn[12*array[i].block + edge_index];
      topo->getEdge(edge_num, &edge);
      const char *edge_attr = edge->getAttribute();
      if (edge_attr && strcmp(edge_attr, attr) == 0){
        queue->push(&array[i]);
      }
    }
    else if (fx || fy || fz){
      // Which face index are we dealing with?
      TMRFace *face;
      int face_index =
        fx*(fx0 ? 0 : 1) + fy*(fy0 ? 2 : 3) + fz*(fz0 ? 4 : 5);
      int face_num = block_face_conn[6*array[i].block + face_index];
      topo->getFace(face_num, &face);  
      const char *face_attr = face->getAttribute();
      if (face_attr && strcmp(face_attr, attr) == 0){
        queue->push(&array[i]);
      }
    }
    else {
      TMRVolume *volume;
      topo->getVolume(array[i].block, &volume);
      const char *volume_attr = volume->getAttribute();
      if (volume_attr && strcmp(volume_attr, attr) == 0){
        queue->push(&array[i]);
      }
    }
  }

  TMROctantArray *list = queue->toArray();
  delete queue;
  return list;
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
  // Get the element octants from this processor
  int nelems = 0;
  TMROctant *array;
  octants->getArray(&array, &nelems);
  
  // Allocate the array to store the connectivity
  int conn_size = 0;
  int *elem_conn = new int[ mesh_order*mesh_order*mesh_order*nelems ];

  for ( int i = 0; i < nelems; i++ ){
    // For all searches/comparisons, we use node numbers
    const int use_node_search = 1;
          
    // Compute the node separation distance
    const int32_t h = 
      1 << (TMR_MAX_LEVEL - array[i].level - (mesh_order-2));
          
    // Loop over all element nodes and determine the global
    // node number
    for ( int kk = 0; kk < mesh_order; kk++ ){
      for ( int jj = 0; jj < mesh_order; jj++ ){
        for ( int ii = 0; ii < mesh_order; ii++ ){
          TMROctant node;
          node.block = array[i].block;
          node.level = 0;
          node.tag = 1;
          node.x = array[i].x + ii*h;
          node.y = array[i].y + jj*h;
          node.z = array[i].z + kk*h;

          // Transform the node to the global ordering
          transformNode(&node);
            
          // Get the node from the sorted list
          TMROctant *t = nodes->contains(&node, use_node_search);
            
          // Set the node number in the element connectivity
          elem_conn[conn_size] = t->tag;
          conn_size++;
        }
      }
    }
  }

  // Set the output arrays
  if (_conn){ *_conn = elem_conn; }
  else { delete [] elem_conn; }
  if (_nelems){ *_nelems = nelems; }
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
  // Go through and label the dependent edges/nodes
  int ndep_faces, ndep_edges;
  TMROctant *face_array, *edge_array;
  dep_faces->getArray(&face_array, &ndep_faces);
  dep_edges->getArray(&edge_array, &ndep_edges);

  // Keep a unique hash of the nodes that are on another
  // processor
  TMROctantHash *ext_hash = new TMROctantHash();

  // Allocate the pointer into the dependent edges
  int *ptr = new int[ num_dep_nodes+1 ];
  memset(ptr, 0, (num_dep_nodes+1)*sizeof(int));

  TMROctantQueue *queue = new TMROctantQueue();

  // Go through the dependent faces and determine the dependent
  // node information
  for ( int i = 0; i < ndep_faces; i++ ){
    // Find the edge length of the element octant
    const int32_t h = 1 << (TMR_MAX_LEVEL - face_array[i].level);
    
    // Find the edge length/node separation along the face
    const int32_t hc = 
      1 << (TMR_MAX_LEVEL - face_array[i].level - (mesh_order-1));
      
    // Get the constraint face index
    int face_index = face_array[i].tag;

    // Loop over the dependent nodes on this face and set the
    // number of associated independent nodes
    for ( int jj = 0; jj < (2*mesh_order-1); jj++ ){
      for ( int ii = 0; ii < (2*mesh_order-1); ii++ ){
        // Get the node location corresponding to this face
        TMROctant node;
        node.block = face_array[i].block;
        node.level = 0;
        node.tag = 0;

        if (face_index < 2){
          node.x = face_array[i].x + (face_index % 2)*h;
          node.y = face_array[i].y + ii*hc;
          node.z = face_array[i].z + jj*hc;
        }
        else if (face_index < 4){
          node.x = face_array[i].x + ii*hc;
          node.y = face_array[i].y + (face_index % 2)*h;
          node.z = face_array[i].z + jj*hc;
        }
        else {
          node.x = face_array[i].x + ii*hc;
          node.y = face_array[i].y + jj*hc;
          node.z = face_array[i].z + (face_index % 2)*h;
        }

        // Convert the node to the global encoding
        transformNode(&node);

        // Determine if this is a dependent or independent ndoe
        int conn_size = 0;
        if ((ii % 2 == 1) && (jj % 2 == 1)){
          conn_size = mesh_order*mesh_order;
        }
        else if ((ii % 2 == 1) || (jj % 2 == 1)){
          conn_size = mesh_order;
        }
        
        // Search for node octant and record its label
        const int use_node_search = 1;
        TMROctant *t = nodes->contains(&node, use_node_search);
        if (conn_size > 0){
          if (t){
            // Get the dependent node number
            int dep_node = -t->tag-1;
            ptr[dep_node+1] = conn_size;
          }
        }
        else if (!t){
          // conn_size == 0 -> this should be an independent node
          // but it does not exist in the local node list (it is
          // not referred to in a local element
          ext_hash->addOctant(&node);
        }
      }
    }
  }

  // Go through the dependent edges and determine the dependent
  // node information
  for ( int i = 0; i < ndep_edges; i++ ){
    // Find the edge length of the element octant
    const int32_t h = 1 << (TMR_MAX_LEVEL - edge_array[i].level);
    
    // Find the edge length/node separation
    const int32_t hc = 
      1 << (TMR_MAX_LEVEL - edge_array[i].level - (mesh_order-1));
      
    // Get the constraint edge index
    int edge_index = edge_array[i].tag;
      
    // Loop over all the nodes on this edge
    for ( int ii = 0; ii < (2*mesh_order-1); ii++ ){
      // Get the node location
      TMROctant node;
      node.block = edge_array[i].block;
      node.level = 0;
      node.tag = 0;

      if (edge_index < 4){
        node.x = edge_array[i].x + ii*hc;
        node.y = edge_array[i].y + h*(edge_index % 2);
        node.z = edge_array[i].z + h*(edge_index/2);
      }
      else if (edge_index < 8){
        node.x = edge_array[i].x + h*(edge_index % 2);
        node.y = edge_array[i].y + ii*hc;
        node.z = edge_array[i].z + h*((edge_index-4)/2);
      }
      else {
        node.x = edge_array[i].x + h*(edge_index % 2);
        node.y = edge_array[i].y + h*((edge_index-8)/2);
        node.z = edge_array[i].z + ii*hc;
      }

      // Convert the node to the global encoding
      transformNode(&node);

      // Search for node octant and record its label
      const int use_node_search = 1;
      TMROctant *t = nodes->contains(&node, use_node_search);
      if (ii % 2 == 1){
        if (t){
          // Get the dependent node number
          int dep_node = -t->tag-1;
          ptr[dep_node+1] = mesh_order;
        }
      }
      else if (!t){
        ext_hash->addOctant(&node);
      }
    }
  }

  // Count up the offsets to the pointers
  for ( int i = 0; i < num_dep_nodes; i++ ){
    ptr[i+1] += ptr[i];
  }

  // Create an array of all the independent nodes that are
  // owned by other processors and not referenced on the
  // elements on this processor.
  TMROctantArray *ext_array = ext_hash->toArray();
  delete ext_hash;
  ext_array->sort();

  // Distribute the non-local nodes back to their owning processors
  // to determine their node numbers
  int use_tags = 0;
  int *send_ptr, *recv_ptr;
  TMROctantArray *dist = distributeOctants(ext_array, use_tags,
                                           &send_ptr, &recv_ptr);
  delete ext_array;

  // Loop over the off-processor nodes and search for them in the
  // sorted node list and assign them the correct tag
  int size;
  TMROctant *array;
  dist->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    const int use_node_search = 1;
    TMROctant *t = nodes->contains(&array[i], use_node_search);
    array[i].tag = t->tag;
  }

  // Send the nodes back to the original processors
  TMROctantArray *ext_nodes = sendOctants(dist, recv_ptr, send_ptr);
  delete dist;
  delete [] recv_ptr;
  delete [] send_ptr;

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
  for ( int i = 0; i < ndep_faces; i++ ){
    // Get the face node numbers
    int fn[25];

    // Find the edge length of the element octant
    const int32_t h = 1 << (TMR_MAX_LEVEL - face_array[i].level);
    
    // Find the edge length/node separation along the face
    const int32_t hc = 
      1 << (TMR_MAX_LEVEL - face_array[i].level - (mesh_order-1));
      
    // Get the constraint face index
    int face_index = face_array[i].tag;

    // Loop over the dependent nodes on this face and set the
    // number of associated independent nodes
    for ( int jj = 0; jj < (2*mesh_order-1); jj++ ){
      for ( int ii = 0; ii < (2*mesh_order-1); ii++ ){
        // Get the node location corresponding to this face
        TMROctant node;
        node.block = face_array[i].block;
        node.level = 0;
        node.tag = 0;

        if (face_index < 2){
          node.x = face_array[i].x + (face_index % 2)*h;
          node.y = face_array[i].y + ii*hc;
          node.z = face_array[i].z + jj*hc;
        }
        else if (face_index < 4){
          node.x = face_array[i].x + ii*hc;
          node.y = face_array[i].y + (face_index % 2)*h;
          node.z = face_array[i].z + jj*hc;
        }
        else {
          node.x = face_array[i].x + ii*hc;
          node.y = face_array[i].y + jj*hc;
          node.z = face_array[i].z + (face_index % 2)*h;
        }

        // Convert the node to the global encoding
        transformNode(&node);

        // Determine if this is a dependent or independent ndoe
        int conn_size = 0;
        if ((ii % 2 == 1) && (jj % 2 == 1)){
          conn_size = mesh_order*mesh_order;
        }
        else if ((ii % 2 == 1) || (jj % 2 == 1)){
          conn_size = mesh_order;
        }

        // Search the node list
        const int use_node_search = 1;
        TMROctant *t = nodes->contains(&node, use_node_search);

        // Record the face nodes
        fn[ii + (2*mesh_order-1)*jj] = -num_dep_nodes-1;
        if (t){
          fn[ii + (2*mesh_order-1)*jj] = t->tag;
        }
        else if ((conn_size == 0) && !t){
          t = ext_nodes->contains(&node, use_node_search);
          fn[ii + (2*mesh_order-1)*jj] = t->tag;
        }
      }
    }
 
    // Loop over all the nodes on the dependent face
    const int n = 2*mesh_order-1;
    const int m = mesh_order-1;

    for ( int jj = 0; jj < n; jj++ ){
      for ( int ii = 0; ii < n; ii++ ){        
        // Get the dependent node number
        int node = -fn[ii + n*jj]-1;

        if (node != num_dep_nodes){
          if ((ii % 2 == 1) && (jj % 2 == 1)){
            for ( int jk = 0, kk = 0; jk < mesh_order; jk++ ){
              for ( int ik = 0; ik < mesh_order; ik++, kk++ ){
                conn[ptr[node] + kk] = fn[2*ik + (2*jk)*n];
                weights[ptr[node] + kk] = wt[ii/2][ik]*wt[jj/2][jk];
              }
            }
          }
          else if (ii % 2 == 1){
            for ( int ik = 0; ik < mesh_order; ik++ ){
              conn[ptr[node] + ik] = fn[2*ik + jj*n];
              weights[ptr[node] + ik] = wt[jj/2][ik];
            }
          }
          else if (jj % 2 == 1){
            for ( int jk = 0; jk < mesh_order; jk++ ){
              conn[ptr[node] + jk] = fn[ii + (2*jk)*n];
              weights[ptr[node] + jk] = wt[ii/2][jk];
            }
          }
        }
      }
    }
  }

  // Loop over the dependent edges
  for ( int i = 0; i < ndep_edges; i++ ){
    // The edge node numbers
    int en[5];

    // Find the edge length of the element octant
    const int32_t h = 1 << (TMR_MAX_LEVEL - edge_array[i].level);
    
    // Find the edge length/node separation
    const int32_t hc = 
      1 << (TMR_MAX_LEVEL - edge_array[i].level - (mesh_order-1));
      
    // Get the constraint edge index
    int edge_index = edge_array[i].tag;
      
    // Loop over all the nodes on this edge
    for ( int ii = 0; ii < (2*mesh_order-1); ii++ ){
      // Get the node location
      TMROctant node;
      node.block = edge_array[i].block;
      node.level = 0;
      node.tag = 0;

      if (edge_index < 4){
        node.x = edge_array[i].x + ii*hc;
        node.y = edge_array[i].y + h*(edge_index % 2);
        node.z = edge_array[i].z + h*(edge_index/2);
      }
      else if (edge_index < 8){
        node.x = edge_array[i].x + h*(edge_index % 2);
        node.y = edge_array[i].y + ii*hc;
        node.z = edge_array[i].z + h*((edge_index-4)/2);
      }
      else {
        node.x = edge_array[i].x + h*(edge_index % 2);
        node.y = edge_array[i].y + h*((edge_index-8)/2);
        node.z = edge_array[i].z + ii*hc;
      }

      // Convert the node to the global encoding
      transformNode(&node);

      // Compute the independent node size
      int conn_size = 0;
      if (ii % 2 == 1){
        conn_size = mesh_order;
      }

      const int use_node_search = 1;
      TMROctant *t = nodes->contains(&node, use_node_search);
      en[ii] = -num_dep_nodes-1;
      if (t){
        en[ii] = t->tag;
      }
      else if ((conn_size == 0) && !t){
        t = ext_nodes->contains(&node, use_node_search);
        en[ii] = t->tag;
      }
    }

    // Loop over all the nodes on the dependent edges
    for ( int jj = 1; jj < 2*mesh_order-1; jj += 2 ){
      int node = -en[jj]-1;
      if (node != num_dep_nodes){
        if (jj % 2 == 1){
          for ( int j = 0; j < mesh_order; j++ ){
            conn[ptr[node] + j] = en[2*j];
            weights[ptr[node] + j] = wt[jj/2][j];
          }
        }
      }
    }
  }
  
  // Free the external nodes
  delete ext_nodes;

  // Return the weights
  *_ptr = ptr;
  *_conn = conn;
  *_weights = weights;
}

/*
  Given a node, find the enclosing octant

  This code is used to find the octant in the octant array that
  encloses the given node.
*/
TMROctant* TMROctForest::findEnclosing( TMROctant *node ){
  // Retrieve the array of elements
  int size = 0;
  TMROctant *array = NULL;
  octants->getArray(&array, &size);

  // Set the lower and upper bounds for the octant
  const int32_t block = node->block;
  const int32_t x = node->x;
  const int32_t y = node->y;
  const int32_t z = node->z;

  // Set the low and high indices to the first and last
  // element of the element array
  int low = 0;
  int high = size-1;
  int mid = low + (int)((high - low)/2);

  // Maintain values of low/high and mid such that the
  // octant is between (elems[low], elems[high]).
  // Note that if high-low=1, then mid = high
  while (high != mid){
    // Check if array[mid] contains the provided octant
    const int32_t h = 1 << (TMR_MAX_LEVEL - array[mid].level);
    if ((array[mid].block == block) &&
        (array[mid].x <= x && x <= array[mid].x+h) &&
	(array[mid].y <= y && y <= array[mid].y+h) &&
	(array[mid].z <= z && z <= array[mid].z+h)){
      return &array[mid];
    }
    
    // Compare the ordering of the two octants - if the
    // octant is less than the other, then adjust the mid point 
    if (node->compareEncoding(&array[mid]) < 0){
      high = mid-1;
    } 
    else {
      low = mid+1;
    }
    
    // Re compute the mid-point and repeat
    mid = high - (int)((high - low)/2);
  }

  // Check if array[mid] contains the provided octant
  const int32_t h1 = 1 << (TMR_MAX_LEVEL - array[mid].level);
  if ((array[mid].block == block) &&
      (array[mid].x <= x && x <= array[mid].x+h1) &&
      (array[mid].y <= y && y <= array[mid].y+h1) &&
      (array[mid].z <= z && z <= array[mid].z+h1)){
    return &array[mid];
  }

  // Check if elems[mid] contains the provided octant
  const int32_t h2 = 1 << (TMR_MAX_LEVEL - array[low].level);
  if ((array[low].block == block) &&
      (array[low].x <= x && x <= array[low].x+h2) &&
      (array[low].y <= y && y <= array[low].y+h2) &&
      (array[low].z <= z && z <= array[low].z+h2)){
   return &array[low];
 }

 // No octant was found, return NULL
 return NULL;
}

/*
  Compute the interpolation
*/
int TMROctForest::computeInterpWeights( const int order,
                                        const int32_t u, const int32_t h,
                                        double Nu[] ){
  if (u == 0 || u == h){
    Nu[0] = 1.0;
    return 1;
  }
  else if (order == 2){
    double ud = 1.0*u/h;
    Nu[0] = 1.0 - ud;
    Nu[1] = ud;
    return 2;
  }
  else {
    double ud = 1.0*u/h;
    Nu[0] = 2.0*(0.5 - ud)*(1.0 - ud);
    Nu[1] = 4.0*ud*(1.0 - ud);
    Nu[2] = 2.0*ud*(ud - 0.5);
    return 3;
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
                                        TACSBVecInterp *interp ){
  // Get the dependent node information
  const int *cdep_ptr, *cdep_conn;
  const double *cdep_weights;
  coarse->getDepNodeConn(&cdep_ptr, &cdep_conn, &cdep_weights);

  // Get the node array
  int node_size;
  TMROctant *node_array;
  nodes->getArray(&node_array, &node_size);
  
  // Find the size of the local node array
  int local_size = node_range[mpi_rank+1] - node_range[mpi_rank];
  TMROctant *local_array = new TMROctant[ local_size ];

  // Read out the nodes that are locally owned
  for ( int i = 0, j = 0; i < node_size; i++ ){
    if (node_array[i].tag >= node_range[mpi_rank] &&
        node_array[i].tag < node_range[mpi_rank+1]){
      local_array[j] = node_array[i];
      j++;
    }
  }

  // Copy the locally owned nodes to the allocated array
  TMROctantArray *local = new TMROctantArray(local_array, local_size);

  // Distribute the octants to the owners - include the local octants
  // in the new array since everything has to be interpolated
  int use_tags = 0; // Use the octant ownership to distribute (not the tags)
  int include_local = 1; // Include the locally owned octants
  TMROctantArray *fine_nodes = 
    coarse->distributeOctants(local, use_tags, NULL, NULL, include_local);
  delete local;

  // Get the number of locally owned nodes on this processor
  int fine_size;
  TMROctant *fine;
  fine_nodes->getArray(&fine, &fine_size);

  // Loop over the nodes in the fine mesh that are owned by octants
  // on the coarse mesh stored on this processor
  for ( int i = 0; i < fine_size; i++ ){
    // Find an octant that encloses the node - this is not unique, but
    // does produce a unique interpolation (since edges/face/corners
    // will be treated the same when adjacent octants that both share
    // a common node location touch)
    TMROctant *oct = coarse->findEnclosing(&fine[i]);
    
    if (oct){
      // The maximum possible size of the array of weights. Note
      // that this is found if every node is a dependent node (which is
      // impossible) which points to a dependent face node (also
      // impossible). It is an upper bound.
      const int max_size = (3*3*3)*(3*3);
      TMRIndexWeight weights[max_size];

      // Get the element size for coarse element
      const int32_t h = 1 << (TMR_MAX_LEVEL - oct->level);
      const int32_t hc = 
        1 << (TMR_MAX_LEVEL - oct->level - (coarse->mesh_order - 2));

      // Compute the parametric location
      int32_t u = fine[i].x - oct->x;
      int32_t v = fine[i].y - oct->y;
      int32_t w = fine[i].z - oct->z;

      // Set the base node location
      int32_t x = oct->x + (u == h ? h : 0);
      int32_t y = oct->y + (v == h ? h : 0);
      int32_t z = oct->z + (w == h ? h : 0);

      // Compute the interpolation weights
      double Nu[3], Nv[3], Nw[3];
      int nu = computeInterpWeights(coarse->mesh_order, u, h, Nu);
      int nv = computeInterpWeights(coarse->mesh_order, v, h, Nv);
      int nw = computeInterpWeights(coarse->mesh_order, w, h, Nw);
    
      // Loop over the nodes that are within this octant
      int nweights = 0;
      for ( int kk = 0; kk < nw; kk++ ){
        for ( int jj = 0; jj < nv; jj++ ){
          for ( int ii = 0; ii < nu; ii++ ){
            // Compute the interpolation weight
            double weight = Nu[ii]*Nv[jj]*Nw[kk];

            // Set the node locations
            TMROctant node;
            node.block = oct->block;
            node.x = x + hc*ii;
            node.y = y + hc*jj;
            node.z = z + hc*kk;
            
            // Transform the node using the coarse transform
            coarse->transformNode(&node);

            // Find the coarse mesh
            int use_node_search = 1;
            TMROctant *t = coarse->nodes->contains(&node, use_node_search);

            if (t->tag >= 0){
              weights[nweights].index = t->tag;
              weights[nweights].weight = weight;
              nweights++;
            }
            else {
              // Unravel the dependent node connectivity
              int node = -t->tag-1;
              for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
                weights[nweights].index = cdep_conn[jp];
                weights[nweights].weight = weight*cdep_weights[jp];
                nweights++;
              }
            }
          }
        }
      }

      // Sort the dependent weight values
      nweights = TMRIndexWeight::uniqueSort(weights, nweights);

      // The interpolation variables/weights on the coarse mesh
      int vars[27];
      double wvals[27];
      for ( int k = 0; k < nweights; k++ ){
        vars[k] = weights[k].index;
        wvals[k] = weights[k].weight;
      }

      // Add the weights/indices to the interpolation object
      interp->addInterp(fine[i].tag, wvals, vars, nweights);
    }
  }

  delete fine_nodes;
}

/*
  Create a sorted, unique array of the external node numbers that are
  referenced on this processor, but are not local.

  This function can be used to determine a local order for the nodes
  on this processor.

  in/out:
  ext_nodes:   the external nodes
  
  returns:     the number of external nodes
*/
int TMROctForest::getExtNodeNums( int **_ext_nodes ){
  // Determine the number of fine nodes
  int node_size;
  TMROctant *node_array;
  nodes->getArray(&node_array, &node_size);

  // Number of local nodes
  int local_size = node_range[mpi_rank+1] - node_range[mpi_rank];
    
  // The maximum number of external nodes
  int max_ext_nodes = node_size - local_size;
  int num_ext = 0;
  int *ext_nodes = new int[ max_ext_nodes ];

  // Scan through and add any external, independent nodes
  for ( int i = 0; i < node_size; i++ ){
    if (node_array[i].tag >= 0 &&
        (node_array[i].tag < node_range[mpi_rank] ||
         node_array[i].tag >= node_range[mpi_rank+1])){
      ext_nodes[num_ext] = node_array[i].tag;
      num_ext++;
    }
  }

  // Sort the array using quicksort
  qsort(ext_nodes, num_ext, sizeof(int), compare_integers);

  // Create a smaller array to store the result
  *_ext_nodes = new int[ num_ext ];
  memcpy(*_ext_nodes, ext_nodes, num_ext*sizeof(int));

  // Free the larger array
  delete [] ext_nodes;

  return num_ext;
}
