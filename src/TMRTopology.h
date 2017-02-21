#ifndef TMR_TOPOLOGY_H
#define TMR_TOPOLOGY_H

#include "TMRBase.h"
#include "TMRGeometry.h"

/*
  The following header file includes the definitions of the
  topological entries within a multiblock mesh.

  These entries refer to the underlying geometry definitions and
  provide connectivity for edges, faces and blocks within a
  multi-block mesh.
*/

/*
  The main topology class that contains the objects used to build the
  underlying mesh.
*/
class TMRTopology : public TMREntity {
 public:
  TMRTopology( TMRGeometry *geo );
  ~TMRTopology();
  
  // Retrieve the face/edge/node information
  void getSurface( int face_num, TMRSurface **face );
  void getFaceCurve( int face_num, int edge_index, TMRCurve **curve );

  // Retrive the connectivity from the topology object
  void getConnectivity( int *nnodes, int *nedges, int *nfaces,
                        const int **face_nodes, const int **face_edges );

 private:
  // The connectivity information
  int *edge_to_vertices;
  int *face_to_edges;
  int *face_to_vertices;

  // The geometry class
  TMRGeometry *geo;
};

#endif // TMR_TOPOLOGY_H
