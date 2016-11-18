#ifndef TMR_OCTANT_FOREST_H
#define TMR_OCTANT_FOREST_H

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

#include "TMROctant.h"

/*
  TMR Forest class

  This class defines a forest of octrees. The octrees within the
  forest can be distributed across processors. The connectivity
  between octrees is defined on all processors by setting a octree to
  node connectivity.
  
  The octrees can be redistributed across processors by using the
  repartition function. This destroys the nodes that may have been
  created (but can easily be recomputed).
  
  The duplicate() and coarsen() functions create a forest that is
  aligned with the parallel distribution of octrees. This facilitates
  the construction of the interpolation operators that can be used for
  multigrid solution algorithms.
*/

class TMROctForest {
 public:
  TMROctForest( MPI_Comm _comm );
  ~TMROctForest();

  // Get the MPI communicator
  // ------------------------
  MPI_Comm getMPIComm(){ return comm; }

  // Set the connectivity
  // --------------------
  void setConnectivity( int _num_nodes,
                        const int *_block_conn,
                        int _num_blocks,
                        int partition=0 );

  // Re-partition the octrees based on element count
  // -----------------------------------------------
  void repartition();
  
  /*
  // Create the forest of octrees
  // ----------------------------
  void createTrees( int refine_level );
  void createTrees( int refine_levels[] );
  */
  void createRandomTrees( int nrand=10, 
                          int min_level=0, int max_level=8 );

  // Duplicate or coarsen the forest
  // -------------------------------
  TMROctForest *duplicate();
  TMROctForest *coarsen();

  // Balance the octree meshes
  // -------------------------
  void balance( int balance_corner=0 );

  // Create and order the nodes
  // --------------------------
  void createNodes( int order=2 );

  // Create the mesh connectivity
  // ----------------------------
  // void createMeshConn( int **_conn, int *_nelems );

  // Retrieve the dependent mesh nodes
  // ---------------------------------
  // int getDepNodeConn( const int **_ptr, const int **_conn,
  //                    const double **_weights );
 
  // Create interpolation/restriction operators
  // ------------------------------------------
  // void createInterpolation( TMROctForest *coarse, 
  //                          int **_interp_ptr, int **_interp_conn,
  //                          double **_interp_weights );

  // Get the external node numbers
  // -----------------------------
  // int getExtNodeNums( int **_extNodes );
  
  // Get the mesh order
  // ------------------
  // int getMeshOrder(){ return mesh_order; }

  // Get the node-processor ownership range
  // --------------------------------------
  /*
    void getOwnedNodeRange( const int **_node_range ){
    if (_node_range){
      *_node_range = node_range;
    }
  }
  */

  // Get the octants
  // ---------------
  void getOctants( TMROctantArray **_octants ){
    if (_octants){
      *_octants = octants;
    }
  }

  // Retrieve the connectivity information
  // -------------------------------------
  void getConnectivity( int *_nblocks, int *_nfaces, 
                        int *_nedges, int *_nnodes, 
                        const int **_block_conn, 
                        const int **_block_face_conn, 
                        const int **_block_edge_conn,
                        const int **_block_face_ids );
  void getInverseConnectivity( const int **_node_block_conn,
                               const int **_node_block_ptr,
                               const int **_edge_block_conn,
                               const int **_edge_block_ptr,
                               const int **_face_block_conn,
                               const int **_face_block_ptr );
  
  // private:
  // Compute the partition using METIS
  // ---------------------------------
  // void computePartition( int part_size, int *vwgts, int *part );

  // Free/copy the allocated data
  void freeData();
  void copyData( TMROctForest *copy );

  // Transform the octant to the global order
  void transformNode( TMROctant *oct );

  // Get the octant owner
  int getOctantMPIOwner( TMROctant *oct );

  // match the ownership intervals
  void matchOctantIntervals( TMROctant *array,
                             int size, int *ptr );
  void matchMPIIntervals( TMROctant *array,
                          int size, int *ptr );

  // Distribute the octant array
  TMROctantArray *distributeOctants( TMROctantArray *list,
                                     int use_tags=0,
                                     int **oct_ptr=NULL, 
                                     int **oct_recv_ptr=NULL );
  TMROctantArray *sendOctants( TMROctantArray *list,
                               const int *oct_ptr,
                               const int *oct_recv_ptr );

  // Balance-related routines
  // ------------------------
  // Balance the octant across the local tree and the forest
  void balanceOctant( TMROctant *oct,
                      TMROctantHash *hash, TMROctantHash *ext_hash,
                      TMROctantQueue *queue,
                      const int balance_corner,
                      const int balance_tree );

  // Add adjacent quadrants to the hashes/queues for balancing
  void addFaceNeighbors( int face_index, 
                         TMROctant p,
                         TMROctantHash *hash,
                         TMROctantHash *ext_hash,
                         TMROctantQueue *queue );
  void addEdgeNeighbors( int edge_index, 
                         TMROctant p,
                         TMROctantHash *hash,
                         TMROctantHash *ext_hash,
                         TMROctantQueue *queue );
  void addCornerNeighbors( int corner, 
                           TMROctant p,
                           TMROctantHash *hash,
                           TMROctantHash *ext_hash,
                           TMROctantQueue *queue );

  // Nodal ordering routines
  // -----------------------
  // Add octants to adjacent non-owner processor queues
  void addAdjacentFaceToQueue( int face_index,
                               TMROctant p,
                               TMROctantQueue *queue, 
                               TMROctant orig );
  void addAdjacentEdgeToQueue( int edge_index,
                               TMROctant p,
                               TMROctantQueue *queue, 
                               TMROctant orig );
  void addAdjacentCornerToQueue( int corner,
                                 TMROctant p,
                                 TMROctantQueue *queue, 
                                 TMROctant orig );

  // Exchange non-local octant neighbors
  void computeAdjacentOctants();

  // Find the dependent faces and edges in the mesh
  void computeDepFacesAndEdges();
  int checkAdjacentDepFaces( int face_index, TMROctant *b,
                             TMROctantArray *adjocts );
  int checkAdjacentDepEdges( int edge_index, TMROctant *b,
                             TMROctantArray *adjocts );
  
  // Label the dependent nodes on the locally owned blocks
  void labelDependentNodes();

  /*
  // Create the dependent node connectivity
  void createDepNodeConn( int **_ptr, int **_conn,
                          double **_weights );

  // Add the nodal weighting values to an interpolant
  void addNodeWeights( TMROctant *t, double w,
                       const int *cdep_ptr, const int *cdep_conn,
                       const double *cdep_weights,
                       TMRIndexWeight *weights, int *nweights );
  */

  // The communicator
  MPI_Comm comm;
  int mpi_rank, mpi_size;

  // The following data is the same across all processors
  // ----------------------------------------------------
  // Set the nodes/edges/faces/blocks
  int num_nodes, num_edges, num_faces, num_blocks;

  // Information for the face/edge/node connectivity
  int *block_conn, *block_face_conn, *block_edge_conn;
  int *node_block_ptr, *node_block_conn;
  int *edge_block_ptr, *edge_block_conn;
  int *face_block_ptr, *face_block_conn;

  // Store the face/edge/node owners
  int *face_block_owners;
  int *edge_block_owners;
  int *node_block_owners;

  // Information to enable transformations between faces
  int *block_face_ids;

  // Information about the mesh
  int mesh_order;

  // Set the range of nodes owned by each processor
  int *node_range;

  // The owner octants 
  TMROctant *owners;

  // The array of all octants
  TMROctantArray *octants;
  
  // The octants that are adjacent to this processor
  TMROctantArray *adjacent;

  // The array of all the nodes
  TMROctantArray *nodes;

  // The following data is processor-local
  // -------------------------------------
  // The number of elements and dependent nodes
  int num_elements;
  int num_dep_nodes;
  int *dep_ptr, *dep_conn;
  double *dep_weights;

  // Pointers to the dependent faces/edges
  TMROctantArray *dep_faces;
  TMROctantArray *dep_edges;
};

#endif // TMR_OCTANT_FOREST_H
