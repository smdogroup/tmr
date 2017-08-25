#ifndef TMR_OCTANT_FOREST_H
#define TMR_OCTANT_FOREST_H

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

#include "TMRTopology.h"
#include "TMROctant.h"
#include "BVecInterp.h"

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
class TMROctForest : public TMREntity {
 public:
  TMROctForest( MPI_Comm _comm );
  ~TMROctForest();

  // Get the MPI communicator
  // ------------------------
  MPI_Comm getMPIComm(){ return comm; }

  // Set the topology (and determine the connectivity)
  // -------------------------------------------------
  void setTopology( TMRTopology *_topo );
  TMRTopology *getTopology();

  // Set the connectivity
  // --------------------
  void setConnectivity( int _num_nodes,
                        const int *_block_conn,
                        int _num_blocks );
  void setFullConnectivity( int _num_nodes, int _num_edges,
                            int _num_faces, int _num_blocks,
                            const int *_block_conn,
                            const int *_block_edge_conn,
                            const int *_block_face_conn );

  // Re-partition the octrees based on element count
  // -----------------------------------------------
  void repartition();
  
  // Create the forest of octrees
  // ----------------------------
  void createTrees( int refine_level );
  void createRandomTrees( int nrand=10, 
                          int min_level=0, int max_level=8 );

  // Duplicate or coarsen the forest
  // -------------------------------
  TMROctForest *duplicate();
  TMROctForest *coarsen();

  // Refine the mesh
  // ---------------
  void refine( const int refinement[]=NULL,
               int min_level=0, int max_level=TMR_MAX_LEVEL );

  // Balance the octree meshes
  // -------------------------
  void balance( int balance_corner=0 );

  // Create and order the nodes
  // --------------------------
  void createNodes( int order=2 );

  // Expert routines for first allocating then ordering nodes (needed
  // when the nodes have non-uniform dof/node)
  // ----------------------------------------------------------------
  void allocateNodes( int order=2 );
  void orderNodes();

  // Get the nodes or elements with certain attributes
  // -------------------------------------------------
  TMROctantArray* getOctsWithAttribute( const char *attr );
  TMROctantArray* getNodesWithAttribute( const char *attr, 
                                         int intersect=1 );

  // Create the mesh connectivity
  // ----------------------------
  void createMeshConn( int **_conn, int *_nelems );

  // Retrieve the dependent mesh nodes
  // ---------------------------------
  void createDepNodeConn();
  int getDepNodeConn( const int **_ptr, const int **_conn,
                      const double **_weights );
 
  // Create interpolation/restriction operators
  // ------------------------------------------
  void createInterpolation( TMROctForest *coarse,
                            TACSBVecInterp *interp );

  // Get the external node numbers
  // -----------------------------
  int getExtNodeNums( int **_extNodes );
  
  // Get the mesh order
  // ------------------
  int getMeshOrder(){ return mesh_order; }

  // Get the node-processor ownership range
  // --------------------------------------
  int getOwnedNodeRange( const int **_node_range ){
    if (_node_range){
      *_node_range = node_range;
    }
    return mpi_size;
  }

  // Get the octants and the nodes
  // -----------------------------
  void getOctants( TMROctantArray **_octants ){
    if (_octants){ *_octants = octants; }
  }
  void getNodes( TMROctantArray **_nodes ){
    if (_nodes){ *_nodes = nodes; }
  }
  void getPoints( TMRPoint **_X ){
    if (_X){ *_X = X; }
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

  // Find the octant enclosing the given node
  // ----------------------------------------
  TMROctant* findEnclosing( TMROctant *node );

  // Transform the octant to the global order
  // ----------------------------------------
  void transformNode( TMROctant *oct );

  // Distribute the octant array in parallel to other processors
  // -----------------------------------------------------------
  TMROctantArray *distributeOctants( TMROctantArray *list,
                                     int use_tags=0,
                                     int **oct_ptr=NULL, 
                                     int **oct_recv_ptr=NULL,
                                     int include_local=0 );

  // Send the octants back to their original processors (dual of distribute)
  // -----------------------------------------------------------------------
  TMROctantArray *sendOctants( TMROctantArray *list,
                               const int *oct_ptr,
                               const int *oct_recv_ptr );

  // Write out files showing the connectivity
  // ----------------------------------------
  void writeToVTK( const char *filename );
  void writeToTecplot( const char *filename );
  void writeForestToVTK( const char *filename );

 private:
  // Compute the node connectivity information
  void computeNodesToBlocks();

  // Compute the connectivity information
  void computeEdgesFromNodes();
  void computeFacesFromNodes();

  // Compute the inverse connectivities
  void computeEdgesToBlocks();
  void computeFacesToBlocks();
  
  // Set the owners - this determines how the mesh will be ordered
  void computeBlockOwners();

  // Free/copy the allocated data
  void freeData();
  void copyData( TMROctForest *copy );

  // Get the octant owner
  int getOctantMPIOwner( TMROctant *oct );

  // match the ownership intervals
  void matchOctantIntervals( TMROctant *array,
                             int size, int *ptr );
  void matchMPIIntervals( TMROctant *array,
                          int size, int *ptr );

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
  
  // Set the dependent node locations
  void setDepNodeLocations();

  // Label the dependent nodes on the locally owned blocks
  void labelDependentNodes();

  // Create the dependent node connectivity
  void createDepNodeConn( int **_ptr, int **_conn,
                          double **_weights );

  // Compute the interpolation weights
  int computeInterpWeights( const int order,
                            const int32_t u, const int32_t h,
                            double Nu[] );

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
  TMRPoint *X;

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

  // The topology of the underlying model (if any)
  TMRTopology *topo;
};

#endif // TMR_OCTANT_FOREST_H
