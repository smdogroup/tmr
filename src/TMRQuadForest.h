#ifndef TMR_QUADTREE_FOREST_H
#define TMR_QUADTREE_FOREST_H

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

#include "TMRTopology.h"
#include "TMRQuadrant.h"
#include "BVecInterp.h"

/*
  A parallel forest of quadtrees

  This class defines a parallel forest of quadrtrees. The connectivity
  between quadtrees is defined at a global level. The quadrants can
  easily be redistributed across processors using the repartition()
  call.

  The duplicate() and coarsen() calls can be used to create a nested
  sequence of meshes that can be used in conjunction with multigrid
  methods.
*/
class TMRQuadForest : public TMREntity {
 public:
  TMRQuadForest( MPI_Comm _comm );
  ~TMRQuadForest();

  // Get the MPI communicator
  // ------------------------
  MPI_Comm getMPIComm(){ return comm; }

  // Set the topology (and determine the connectivity)
  // -------------------------------------------------
  void setTopology( TMRTopology *_topo );

  // Set the connectivity directly
  // -----------------------------
  void setConnectivity( int _num_nodes,
                        const int *_face_conn,
                        int _num_faces );
  void setFullConnectivity( int _num_nodes, 
                            int _num_edges,
                            int _num_faces, 
                            const int *_face_conn,
                            const int *_face_edge_conn );

  // Re-partition the quadtrees based on element count
  // -------------------------------------------------
  void repartition();

  // Create the forest of quadtrees
  // ----------------------------
  void createTrees( int refine_level );
  void createRandomTrees( int nrand=10, 
                          int min_level=0, int max_level=8 );

  // Refine the mesh
  // ---------------
  void refine( const int refinement[],
               int min_level=0, int max_level=TMR_MAX_LEVEL );

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
                            TACSBVecInterp *interp );

  // Get the external node numbers
  // -----------------------------
  int getExtNodeNums( int **_extNodes );
  
  // Get the mesh order
  // ------------------
  int getMeshOrder(){ return mesh_order; }

  // Get the node-processor ownership range
  // --------------------------------------
  void getOwnedNodeRange( const int **_node_range ){
    if (_node_range){
      *_node_range = node_range;
    }
  }

  // Get the quadrants and the nodes
  // -----------------------------
  void getQuadrants( TMRQuadrantArray **_quadrants ){
    if (_quadrants){ *_quadrants = quadrants; }
  }
  void getNodes( TMRQuadrantArray **_nodes ){
    if (_nodes){ *_nodes = nodes; }
  }
  void getPoints( TMRPoint **_X ){
    if (_X){ *_X = X; }
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

  // Transform the quadrant to the global order
  // ------------------------------------------
  void transformNode( TMRQuadrant *quad );

 private:
  // Free the internally stored data and zero things
  void freeData();
  
  // Duplicate data
  void copyData( TMRQuadForest *copy );

  // Set up the connectivity from nodes -> faces
  void computeNodesToFaces();

  // Set up the edge connectivity
  void computeEdgesFromNodes();
  void computeEdgesToFaces();

  // Compute the faces that own the edges and nodes
  void computeFaceOwners();

  // Get the quadrant owner
  int getQuadrantMPIOwner( TMRQuadrant *quad );

  // match the ownership intervals
  void matchQuadrantIntervals( TMRQuadrant *array,
                               int size, int *ptr );
  void matchMPIIntervals( TMRQuadrant *array,
                          int size, int *ptr );

  // Distribute the quadrant array
  TMRQuadrantArray *distributeQuadrants( TMRQuadrantArray *list,
                                         int use_tags=0,
                                         int **quad_ptr=NULL, 
                                         int **quad_recv_ptr=NULL,
                                         int include_local=0 );
  TMRQuadrantArray *sendQuadrants( TMRQuadrantArray *list,
                                   const int *quad_ptr,
                                   const int *quad_recv_ptr );

  // Balance-related routines
  // ------------------------
  // Balance the quadrant across the local tree and the forest
  void balanceQuadrant( TMRQuadrant *quad,
                        TMRQuadrantHash *hash, TMRQuadrantHash *ext_hash,
                        TMRQuadrantQueue *queue,
                        const int balance_corner,
                        const int balance_tree );

  // Add adjacent quadrants to the hashes/queues for balancing
  void addEdgeNeighbors( int edge_index, 
                         TMRQuadrant p,
                         TMRQuadrantHash *hash,
                         TMRQuadrantHash *ext_hash,
                         TMRQuadrantQueue *queue );
  void addCornerNeighbors( int corner, 
                           TMRQuadrant p,
                           TMRQuadrantHash *hash,
                           TMRQuadrantHash *ext_hash,
                           TMRQuadrantQueue *queue );

  // Nodal ordering routines
  // -----------------------
  // Add quadrants to adjacent non-owner processor queues
  void addAdjacentEdgeToQueue( int edge_index,
                               TMRQuadrant p,
                               TMRQuadrantQueue *queue, 
                               TMRQuadrant orig );
  void addAdjacentCornerToQueue( int corner,
                                 TMRQuadrant p,
                                 TMRQuadrantQueue *queue, 
                                 TMRQuadrant orig );

  // Exchange non-local quadrant neighbors
  void computeAdjacentQuadrants();

  // Find the dependent faces and edges in the mesh
  void computeDepEdges();
  int checkAdjacentDepEdges( int edge_index, TMRQuadrant *b,
                             TMRQuadrantArray *adjquads );
  
  // Label the dependent nodes on the locally owned blocks
  void labelDependentNodes();

  // Create the dependent node connectivity
  void createDepNodeConn( int **_ptr, int **_conn,
                          double **_weights );

  // Find the quadrant 
  TMRQuadrant* findEnclosing( TMRQuadrant *node );

  // Compute the interpolation weights
  int computeInterpWeights( const int order,
                            const int32_t u, const int32_t h,
                            double Nu[] );
  // The communicator
  MPI_Comm comm;
  int mpi_rank, mpi_size;

  // The owner ranges
  TMRQuadrant *owners;

  // The following data is the same across all processors
  // ----------------------------------------------------
  // Set the nodes/edges/faces
  int num_nodes, num_edges, num_faces;

  // Information for the face/edge/node connectivity
  int *face_conn, *face_edge_conn;
  int *node_face_ptr, *node_face_conn;
  int *edge_face_ptr, *edge_face_conn;

  // Set the node/edge owners
  int *node_face_owners, *edge_face_owners;

  // Information about the mesh
  int mesh_order;

  // Set the range of nodes owned by each processor
  int *node_range;

  // The array of all quadrants
  TMRQuadrantArray *quadrants;
  
  // The quadrants that are adjacent to this processor
  TMRQuadrantArray *adjacent;

  // The array of all the nodes
  TMRQuadrantArray *nodes;
  TMRPoint *X;

  // The following data is processor-local
  // -------------------------------------
  // The number of elements and dependent nodes
  int num_elements;
  int num_dep_nodes;
  int *dep_ptr, *dep_conn;
  double *dep_weights;

  // Pointers to the dependent edges
  TMRQuadrantArray *dep_edges;

  // The topology of the underlying model (if any)
  TMRTopology *topo;
};

#endif // TMR_QUADTREE_FOREST_H
