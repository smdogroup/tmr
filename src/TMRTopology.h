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
  Topological curve entity between two vertices

  The vertices themselves are extracted from the curve. Note that the
  vertices must be set in the curve.
*/
class TMREdge : public TMREntity {
 public:
  TMREdge( MPI_Comm _comm, TMRCurve *_curve,
           TMRVertex *_v1, TMRVertex *_v2 );
  ~TMREdge();

  // Get the underlying curve associated with this edge
  void getCurve( TMRCurve **_curve );

  // Retrieve the vertices at the beginning/end of the curve
  void getVertices( TMRVertex **_v1, TMRVertex **_v2 );

 private:
  // Record the communicator
  MPI_Comm comm;

  // The edge identifier - the combination of the edge
  // identifier and the communicator are unique
  int edge_id; 
  
  // The underlying curve associated with this edge
  TMRCurve *curve;
};

/*
  Topological surface entity bounded by four non-degenerate edges

  The TMRFace object represents the topology of a non-degenerate
  quadrilateral surface with four edges. Since a general surface
  object can be bounded by three or more edges, and we require a
  specified edge ordering, the edge objects are required arguments.
  
  The edges are ordered as follows, where the orientation must match

  \--------3--------\
  |                 |
  |      v ^        |
  |        |        |
  0        ----> u  1
  |                 |
  |                 |
  |                 |
  \--------2--------\

  Where the normal of the surface is out of plane in this case. Note
  that u/v locations are notional, and the edges need not run along
  constant u/v values on the surface.
*/
class TMRFace : public TMREntity {
 public:
  TMRFace( MPI_Comm comm, TMRSurface *_surface,
           TMREdge *_edges[], int _edge_dir[] );
  ~TMRFace();

  // Get the four edges that bound this face
  void getEdges( TMREdge ***edges, const int **_edge_dir );

 private:
  // The communicator
  MPI_Comm comm;

  // The face identifier -- unique identifier for the
  // face
  int face_id;

  // Pointers to the edges bounding this face
  TMRSurface *surface;
  TMREdge *edges[4];
  int edge_dir[4];
};

/*
  The main topology class that contains the objects used to build the
  underlying mesh.
*/
class TMRTopology : public TMREntity {
 public:
  TMRTopology( MPI_Comm _comm,
               TMRFace **_faces, int num_faces );
  ~TMRTopology();
  
  // Retrieve the face/edge/node information
  void getFace( int face_num, TMRFace **face );

  void getConnectivity( int *nnodes, int *nedges, int *nfaces,
                        const int **face_nodes, const int **edge_nodes );

 private:
  // The communicator associated 
  MPI_Comm comm;

  // Build the underlying connectivity for the mesh topology
  void buildConnectivity();

  // The face information
  int num_faces;
  TMRFace **faces;

  // The edge information
  int num_edges;
  TMREdge **edges;

  // The node information
  int num_nodes;
  TMRVertex **nodes;
};

#endif // TMR_TOPOLOGY_H
