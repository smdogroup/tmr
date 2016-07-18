#ifndef TMR_QUADTREE_FOREST_H
#define TMR_QUADTREE_FOREST_H

#include "TMRQuadtree.h"

class TMRQuadForest {
 public:
  TMRQuadForest( MPI_Comm _comm );
  ~TMRQuadForest();

  // Functions for setting the connecitivity
  // ---------------------------------------
  void setConnectivity( int _num_nodes,
                        const int *_face_conn, 
                        int _num_faces );

  // Create the forest of quadrants
  // ------------------------------
  void createTrees( int refine_level );
  void createRandomTrees( int nrand=10, 
                          int min_level=0, int max_level=8 );

  // Balance the quadrant meshes
  // ---------------------------
  void balance( int balance_corner=0 );

  // Duplicate or coarsen the forest
  // -------------------------------
  TMRQuadForest* duplicate();
  TMRQuadForest *coarsen();

  // Create the mesh connectivity
  // ----------------------------
  void createNodes( int order );

  // Retrieve the mesh connectivity
  // ------------------------------
  void getMesh( int *nnodes, int *ndep_nodes, int *nelems,
                int **_elem_conn, int **_dep_conn,
                double **_dep_weights );

  // Retrieve the individual quadtrees 
  // ---------------------------------
  int getQuadtrees( TMRQuadtree ***_quadtrees ){
    if (_quadtrees){ *_quadtrees = quadtrees; }
    return num_faces;
  }

 private:
  // Add the adjacent quadrant to the hash/queues for balancing
  void addCornerNeighbors( int tree, int corner, 
                           TMRQuadrant p,
                           TMRQuadrantHash **hash,
                           TMRQuadrantQueue **queue );
  void addEdgeNeighbors( int tree, int edge, 
                         TMRQuadrant p,
                         TMRQuadrantHash **hash,
                         TMRQuadrantQueue **queue );
  
  // Label the dependent nodes on adjacent edges
  int labelDependentNodes( int edge,
                           TMRQuadrantArray **edge_nodes,
                           int count );

  // Get/set the node numbers from a specified edge
  void getEdgeNodeNums( int face, int edge, int dest_face, 
                        TMRQuadrantArray *edge_nodes );
  void setEdgeNodeNums( int face, TMRQuadrantArray *edge_nodes );

  // Label the dependent edge nodes
  void addEdgeDependentNodes( int face, int edge,
                              TMRQuadrant p, TMRQuadrant source );

  // The communicator
  MPI_Comm comm;

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
  int num_elements;
  int num_dep_nodes;

  // Keep a pointer to the forest of quadtrees
  TMRQuadtree **quadtrees; 

  // The mpi rank of the face owner
  int *face_owners;
};

#endif // TMR_QUADTREE_FOREST_H
