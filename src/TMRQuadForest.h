#ifndef TMR_QUADTREE_FOREST_H
#define TMR_QUADTREE_FOREST_H

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

#include "TMRQuadtree.h"

/*
  A parallel forest of quadtrees
*/

class TMRQuadForest {
 public:
  TMRQuadForest( MPI_Comm _comm );
  ~TMRQuadForest();

  // Get the MPI communicator
  // ------------------------
  MPI_Comm getMPIComm(){ return comm; }

  // Set the connectivity
  // --------------------
  void setConnectivity( int _num_nodes,
                        const int *_face_conn,
                        int _num_faces,
                        int partition=0 );

  // Re-partition the mesh based on element count
  // --------------------------------------------
  void repartition();

  // Create the forest of quadtrees
  // ----------------------------
  void createTrees( int refine_level );
  void createTrees( int refine_levels[] );
  void createRandomTrees( int nrand=10, 
                          int min_level=0, int max_level=8 );

  // Duplicate or coarsen the forest
  // -------------------------------
  TMRQuadForest *duplicate();
  TMRQuadForest *coarsen();

  // Balance the quadtree meshes
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
  int getDepNodeConn( const int **_ptr, const int **_conn,
                      const double **_weights );

  // Create interpolation/restriction operators
  // ------------------------------------------
  void createInterpolation( TMRQuadForest *coarse, 
                            int **_interp_ptr, int **_interp_conn,
                            double **_interp_weights );

  // Get the array of quadtrees - careful, many are NULL
  // -------------------------------------------------
  int getQuadtrees( TMRQuadtree ***_trees );
  
  // Get mesh/ownership information - short cut to the non-NULL quadtrees
  // ------------------------------------------------------------------
  int getOwnedQuadtrees( const int **_owned_faces );

  // Get the node-processor ownership range
  // --------------------------------------
  void getOwnedNodeRange( const int **_node_range ){
    if (_node_range){
      *_node_range = node_range;
    }
  }

  // Retrieve the connectivity information
  // -------------------------------------
  void getConnectivity( int *_nfaces, int *_nedges, int *_nnodes, 
                        const int **_face_conn, 
                        const int **_face_edge_conn );
  void getInverseConnectivity( const int **_node_face_conn,
                               const int **_node_face_ptr,
                               const int **_edge_face_conn,
                               const int **_edge_face_ptr );

 private:
  // Compute the partition using METIS
  // ---------------------------------
  void computePartition( int part_size, int *vwgts, int *part );

  // Balance-related routines
  // ------------------------
  // Balance the quadrant across the local tree and the forest
  void balanceQuadrant( int face, TMRQuadrant *quad,
                        TMRQuadrantHash **hash, 
                        TMRQuadrantQueue **queue,
                        const int balance_corner,
                        const int balance_tree );

  // Add adjacent quadrants to the hashes/queues for balancing
  void addEdgeNeighbors( int face, int edge_index, 
                         TMRQuadrant p,
                         TMRQuadrantHash **hash,
                         TMRQuadrantQueue **queue );
  void addCornerNeighbors( int face, int corner, 
                           TMRQuadrant p,
                           TMRQuadrantHash **hash,
                           TMRQuadrantQueue **queue );

  // Nodal ordering routines
  // -----------------------
  // Add quadrants to adjacent non-owner processor queues
  void addCornerQuadrantToQueues( const int node, 
                                  const int mpi_rank, 
                                  TMRQuadrant *q,
                                  TMRQuadrantQueue **queues );
  void addEdgeQuadrantToQueues( const int edge, 
                                const int mpi_rank, 
                                TMRQuadrant *q,
                                TMRQuadrantQueue **queues );

  // Exchange non-local quadrant neighbors
  void recvQuadNeighbors();
    
  // Label the dependent nodes on the locally owned faces
  void labelDependentNodes();
  int checkAdjacentDepEdges( int edge, int edge_index,
                             int face_owner, TMRQuadrant *b );
  void computeDepEdges();

  // Get the owner flags
  void getOwnerFlags( int face,
                      const int *edge_face_owners,
                      const int *node_face_owners,
                      int *is_edge_owner, int *is_node_owner );

  // Order the global nodes
  void orderGlobalNodes( const int *edge_face_owners,
                         const int *node_face_owners );

  // Send the ordered nodes back to their owners
  void sendNodeNeighbors( const int *edge_face_owners,
                          const int *node_face_owners );

  // Copy the nodes numbers from one face to an adjacent face
  void copyCornerNodes( int node, int node_index,
                        int face_owner, TMRQuadrant *p );
  void copyEdgeNodes( int edge, int edge_index,
                      int face_owner, TMRQuadrant *p );
  void copyAdjacentNodes( const int *edge_face_owners,
                          const int *node_face_owners );

  // Create the dependent node connectivity
  void createDepNodeConn( int **_ptr, int **_conn,
                          double **_weights );

  // Add the nodal weighting values to an interpolant
  void addNodeWeights( TMRQuadrant *t, double w,
                       const int *cdep_ptr, const int *cdep_conn,
                       const double *cdep_weights,
                       TMRIndexWeight *weights, int *nweights );
  // The communicator
  MPI_Comm comm;

  // The following data is the same across all processors
  // ----------------------------------------------------
  // Set the nodes/edges/faces
  int num_nodes, num_edges, num_faces;

  // Information for the face/edge/node connectivity
  int *face_conn, *face_edge_conn;
  int *node_face_ptr, *node_face_conn;
  int *edge_face_ptr, *edge_face_conn;

  // Information about the mesh
  int mesh_order;

  // Set the range of nodes owned by each processor
  int *node_range;

  // The mpi rank of the face owners
  int *mpi_face_owners;

  // The following data is processor-local
  // -------------------------------------
  // The number of elements and dependent nodes
  int num_elements;
  int num_dep_nodes;
  int *dep_ptr, *dep_conn;
  double *dep_weights;

  // Keep a pointer to the forest of quadtrees
  TMRQuadtree **quadtrees;

  // Pointers to the dependent faces/edges
  TMRQuadrantArray **dep_edges;

  // A short cut to the owned faces
  int num_owned_faces;
  int *owned_faces;
};

#endif // TMR_QUADTREE_FOREST_H
