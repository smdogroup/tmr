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
  TMRQuadForest( MPI_Comm _comm, int mesh_order=2,
                 TMRInterpolationType interp_type=TMR_GAUSS_LOBATTO_POINTS );
  ~TMRQuadForest();

  // Get the MPI communicator
  // ------------------------
  MPI_Comm getMPIComm(){ return comm; }

  // Set the topology (and determine the connectivity)
  // -------------------------------------------------
  void setTopology( TMRTopology *_topo );
  TMRTopology *getTopology();

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
  void refine( const int refinement[]=NULL,
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
  void createNodes();

  // Retrieve the dependent mesh nodes
  // ---------------------------------
  void getNodeConn( const int **_conn=NULL, 
                    int *_num_elements=NULL,
                    int *_num_owned_nodes=NULL,
                    int *_num_local_nodes=NULL );
  int getDepNodeConn( const int **_ptr, const int **_conn,
                      const double **_weights );

  // Create interpolation/restriction operators
  // ------------------------------------------
  void createInterpolation( TMRQuadForest *coarse,
                            TACSBVecInterp *interp );

  // Get the nodes or elements with certain attributes
  // -------------------------------------------------
  TMRQuadrantArray* getQuadsWithAttribute( const char *attr );
  int getNodesWithAttribute( const char *attr, int **nodes );
  
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

  // Get the quadrants and the nodes
  // -----------------------------
  void getQuadrants( TMRQuadrantArray **_quadrants ){
    if (_quadrants){ *_quadrants = quadrants; }
  }
  int getNodeNumbers( const int **_node_numbers ){
    if (_node_numbers){ *_node_numbers = node_numbers; }
    return num_local_nodes;
  }
  int getPoints( TMRPoint **_X ){
    if (_X){ *_X = X; }
    return num_local_nodes;
  }

  // Get the index of the given element or node
  // ------------------------------------------
  int getElementIndex( TMRQuadrant *element );

  // Retrieve the connectivity information
  // -------------------------------------
  void getConnectivity( int *_nfaces, int *_nedges, int *_nnodes, 
                        const int **_face_conn, 
                        const int **_face_edge_conn );
  void getInverseConnectivity( const int **_node_face_conn,
                               const int **_node_face_ptr,
                               const int **_edge_face_conn,
                               const int **_edge_face_ptr );

  // Write out files showing the connectivity
  // ----------------------------------------
  void writeToVTK( const char *filename );
  void writeToTecplot( const char *filename );
  void writeForestToVTK( const char *filename );
  void writeAdjacentToVTK( const char *filename );

  // Distribute the quadrant array
  // -----------------------------
  TMRQuadrantArray *distributeQuadrants( TMRQuadrantArray *list,
                                         int use_tags=0,
                                         int **quad_ptr=NULL, 
                                         int **quad_recv_ptr=NULL,
                                         int include_local=0,
                                         int use_node_index=0 );
  TMRQuadrantArray *sendQuadrants( TMRQuadrantArray *list,
                                   const int *quad_ptr,
                                   const int *quad_recv_ptr,
                                   int use_node_index=0 );

  // Find the quadrant 
  // -----------------
  TMRQuadrant* findEnclosing( TMRQuadrant *node );

 private:
  // Labels for the nodes
  static const int TMR_QUAD_NODE_LABEL = 0;
  static const int TMR_QUAD_EDGE_LABEL = 1;
  static const int TMR_QUAD_FACE_LABEL = 2;

  // Free the internally stored data and zero things
  void freeData();
  void freeMeshData( int free_quads=1, int free_owners=1 );
  
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
  void matchTagIntervals( TMRQuadrant *array,
                          int size, int *ptr );
 
  // Balance-related routines
  // ------------------------
  // Balance the quadrant across the local tree and the forest
  void balanceQuadrant( TMRQuadrant *quad,
                        TMRQuadrantHash *hash, 
                        TMRQuadrantHash *ext_hash,
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
  void computeAdjacentDepEdges( int edge_index, TMRQuadrant *b,
                                TMRQuadrantArray *adjquads );

  // Transform the node to the global numbering. Note that the node
  // may be labeled a node, edge or face
  void transformNode( TMRQuadrant *quad, int *edge_reversed=NULL );

  // Label the dependent nodes in the dependent node list
  void labelDependentNodes( int *nodes );

  // Create the global node ownership data
  TMRQuadrantArray* createLocalNodes();

  // Create the local connectivity based on the input node array
  void createLocalConn( TMRQuadrantArray *nodes, const int *node_offset );

  // Create the dependent node connectivity
  void createDependentConn( const int *node_nums,
                            TMRQuadrantArray *nodes, 
                            const int *node_offset );

  // Compute the node locations
  void evaluateNodeLocations();

  // The communicator 
  MPI_Comm comm;
  int mpi_rank, mpi_size;

  // Information about the type of interpolation
  TMRInterpolationType interp_type;
  double *interp_knots;

  // The owner quadrant ranges for each processor. Note that this is
  // in the quadrant space not the node space
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

  // The mesh order/connectivity information
  int mesh_order;
  int *conn;

  // Set the range of node numbers owned by each processor
  int *node_range;

  // The nodes are organized as follows
  // |--- dependent nodes -- | ext_pre | -- owned local -- | - ext_post -|

  // The following data is processor-local
  int *node_numbers; // All the local node numbers ref'd on this proc
  int num_local_nodes; // Total number of locally ref'd nodes
  int num_dep_nodes; // Number of dependent nodes
  int num_owned_nodes; // Number of nodes that are owned by me
  int ext_pre_offset; // Number of nodes before pre

  // The dependent node information
  int *dep_ptr, *dep_conn;
  double *dep_weights;

  // The array of all quadrants
  TMRQuadrantArray *quadrants;
  
  // The quadrants that are adjacent to this processor
  TMRQuadrantArray *adjacent;

  // The array of all the nodes
  TMRPoint *X;

  // The topology of the underlying model (if any)
  TMRTopology *topo;
};

#endif // TMR_QUADTREE_FOREST_H
