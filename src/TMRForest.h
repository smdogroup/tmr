#ifndef TMR_FORREST_H
#define TMR_FORREST_H

#include "mpi.h"
#include "TMRQuadtree.h"
#include "TMROctree.h"

class TMRQuadForest {
 public:
  TMRQuadForest( MPI_Comm _comm );
  ~TMRQuadForest();

  // Functions for setting the connecitivity
  // ---------------------------------------
  void setConnectivity( int _num_nodes,
                        const int *_face_conn, 
                        int _num_faces );

  // Create the forrest of quadrants
  // -------------------------------
  void createTrees( int refine_level );
  void createRandomTrees( int nrand=10, 
                          int min_level=0, int max_level=8 );


  // Balance the quadrant meshes
  // ---------------------------
  void balance( int balance_corner=0 );

  // Duplicate or coarsen the forrest
  // --------------------------------
  TMRQuadForest* duplicate();
  TMRQuadForest *coarsen();

  // Create the mesh connectivity
  // ----------------------------
  void createNodes( int order );

  // Retrieve the mesh connectivity
  // ------------------------------
  void getMesh( int *nnodes, int *nelems,
                int **_elem_ptr, int **_elem_conn );
  void getDependentNodes( int *num_dep_nodes,
                          const int **_dep_ptr,
                          const int **_dep_conn,
                          const double **_dep_weights );  

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

  // Get the adjacent quadrant for node ordering
  TMRQuadrant* getEdgeNodeNeighbor( int face, 
                                    int adjacent,
                                    int edge, 
                                    TMRQuadrant p );

  // Label dependent nodes
  void addEdgeDependentNodes( int face, int edge,
                              TMRQuadrant p, 
                              TMRQuadrant source,
                              TMRQuadrantQueue **queue,
                              TMRQuadrantQueue **indep );

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
  int num_mesh_nodes;
  int num_mesh_dep_nodes;
  int num_mesh_elements;

  // Keep the dependent weights stored locally
  int *dep_ptr;
  int *dep_conn;
  double *dep_weights;

  // Keep a pointer to the forrest of quadtrees
  TMRQuadtree **quadtrees; 
};

#endif // TMR_FORREST_H
