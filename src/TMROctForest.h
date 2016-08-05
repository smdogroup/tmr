#ifndef TMR_OCTANT_FOREST_H
#define TMR_OCTANT_FOREST_H

#include "TMROctree.h"

class TMROctForest {
 public:
  TMROctForest( MPI_Comm _comm );
  ~TMROctForest();

  // Set the connectivity
  // --------------------
  void setConnectivity( int _num_nodes,
                        const int *_block_conn,
                        int _num_blocks,
                        int partition=0 );

  // Re-partition the mesh based on element count
  // --------------------------------------------
  void repartition();

  // Create the forest of octrees
  // ----------------------------
  void createTrees( int refine_level );
  void createTrees( int refine_levels[] );
  void createRandomTrees( int nrand=10, 
                          int min_level=0, int max_level=8 );

  // Duplicate or coarsen the forest
  // -------------------------------
  TMROctForest* duplicate();
  TMROctForest *coarsen();

  // Balance the octree meshes
  // -------------------------
  void balance( int balance_corner=0 );

  // Create and order the nodes
  // --------------------------
  void createNodes( int order=2 );

  // Create the mesh connectivity
  // ----------------------------
  void createMeshConn( int **_conn, int *_nelems );

  // Retrieve the dependent mesh nodes
  // ---------------------------------
  int createDependentNodes( int **_ptr, int **_conn,
                            double **_weights );

  // Get the array of octrees
  // ------------------------
  int getOctrees( TMROctree ***_trees ){
    if (_trees){ *_trees = octrees; }
    return num_blocks;
  }

 private:
  // Compute the partition using METIS
  // ---------------------------------
  void computePartition( int part_size, int *vwgts, int *part );

  // Balance-related routines
  // ------------------------
  // Balance the octant across the local tree and the forest
  void balanceOctant( int block, TMROctant *oct,
                      TMROctantHash **hash, TMROctantQueue **queue,
                      const int balance_corner,
                      const int balance_tree );

  // Add adjacent quadrants to the hashes/queues for balancing
  void addFaceNeighbors( int block, int face_index, 
                         TMROctant p,
                         TMROctantHash **hash,
                         TMROctantQueue **queue );
  void addEdgeNeighbors( int block, int edge_index, 
                         TMROctant p,
                         TMROctantHash **hash,
                         TMROctantQueue **queue );
  void addCornerNeighbors( int block, int corner, 
                           TMROctant p,
                           TMROctantHash **hash,
                           TMROctantQueue **queue );

  // Nodal ordering routines
  // -----------------------
  // Add octants to adjacent non-owner processor queues
  void addCornerOctantToQueues( const int node, 
                                const int mpi_rank, 
                                TMROctant *q,
                                TMROctantQueue **queues );
  void addEdgeOctantToQueues( const int edge, 
                              const int mpi_rank, 
                              TMROctant *q,
                              TMROctantQueue **queues );
  void addFaceOctantToQueues( const int face, 
                              const int mpi_rank, 
                              TMROctant *q,
                              TMROctantQueue **queues );

  // Exchange non-local octant neighbors
  void recvOctNeighbors();
    
  // Label the dependent nodes on the locally owned blocks
  void labelDependentFaceNodes( const int order );
  int checkAdjacentDepFaces( int face, int face_index,
                             int block_owner, TMROctant *b );
  void computeDependentFaces();

  // Get the owner flags
  void getOwnerFlags( int block, int mpi_rank,
                      const int *face_block_owners,
                      const int *edge_block_owners,
                      const int *node_block_owners,
                      int *is_face_owner, int *is_edge_owner, 
                      int *is_node_owner );

  // Order the global nodes
  void orderGlobalNodes( const int *face_block_owners,
                         const int *edge_block_owners,
                         const int *node_block_owners );

  // Send the ordered nodes back to their owners
  void sendNodeNeighbors( const int *face_block_owners,
                          const int *edge_block_owners,
                          const int *node_block_owners );

  // Copy the nodes numbers from one block to an adjacent block
  void copyCornerNodes( int node, int node_index,
                        int block_owner, TMROctant *p );
  void copyEdgeNodes( int edge, int edge_index,
                      int block_owner, TMROctant *p );
  void copyFaceNodes( int face, int face_index,
                      int block_owner, TMROctant *p );
  void copyAdjacentNodes( const int *face_block_owners,
                          const int *edge_block_owners,
                          const int *node_block_owners );

  // The communicator
  MPI_Comm comm;

  // The following data is the same across all processors
  // ----------------------------------------------------
  // Set the nodes/edges/faces/blocks
  int num_nodes, num_edges, num_faces, num_blocks;

  // Information for the face/edge/node connectivity
  int *block_conn, *block_face_conn, *block_edge_conn;
  int *node_block_ptr, *node_block_conn;
  int *edge_block_ptr, *edge_block_conn;
  int *face_block_ptr, *face_block_conn;

  // Information to enable transformations between faces
  int *block_face_ids;

  // Information about the mesh
  int mesh_order;

  // Set the range of nodes owned by each processor
  int *node_range;

  // The mpi rank of the block owners
  int *mpi_block_owners;

  // The following data is processor-local
  // -------------------------------------
  // Set the elements, nodes and dependent nodes
  int num_elements;
  int num_dep_nodes;

  // Keep a pointer to the forest of quadtrees
  TMROctree **octrees;

  // Pointers to the dependent faces/edges
  TMROctantArray **dep_faces;

  // A short cut to the owned blocks
  int num_owned_blocks;
  int *owned_blocks;
};

#endif // TMR_OCTANT_FOREST_H
