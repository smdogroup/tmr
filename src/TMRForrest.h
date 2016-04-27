#ifndef TMR_FORREST_H
#define TMR_FORREST_H

#include "mpi.h"
#include "TMRQuadtree.h"
#include "TMROctree.h"

class TMRQuadForrest {
 public:
  TMRQuadForrest( MPI_Comm _comm );
  ~TMRQuadForrest();

  // Functions for setting the connecitivity
  // ---------------------------------------
  void setConnectivity( int _num_nodes,
                        const int *_face_conn, 
                        int _num_faces );

  // Create the forrest of quadrants
  // -------------------------------
  void createTrees( int refine_level );
  
  // Add the adjacent quadrant to the hash/queues
  // --------------------------------------------
  void addCornerNeighbors( int tree, int corner, 
                           TMRQuadrant p,
                           TMRQuadrantHash **hash,
                           TMRQuadrantQueue **queue );
  void addEdgeNeighbors( int tree, int edge, 
                         TMRQuadrant p,
                         TMRQuadrantHash **hash,
                         TMRQuadrantQueue **queue );

  // Balance the quadrant meshes
  // ---------------------------
  void balance( int balance_corner=0 );

  // Coarsen the forrest
  // -------------------
  TMRQuadForrest *coarsen();

  // Create the mesh connectivity
  // ----------------------------
  void createNodes( int _order );
  void createMesh( int _order );

  // Retrieve the individual quadtrees 
  // ---------------------------------
  int getQuadtrees( TMRQuadtree ***_quadtrees ){
    if (_quadtrees){ *_quadtrees = quadtrees; }
    return num_faces;
  }

 private:
  // The communicator
  MPI_Comm comm;

  // Set the nodes/edges/faces
  int num_nodes, num_edges, num_faces;

  // Information for the face/edge/node connectivity
  int *face_conn, *face_edge_conn;
  int *node_face_ptr, *node_face_conn;
  int *edge_face_ptr, *edge_face_conn;

  // Keep a pointer to the forrest of quadtrees
  TMRQuadtree **quadtrees; 
};

#endif // TMR_FORREST_H
