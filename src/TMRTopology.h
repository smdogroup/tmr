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

// Declarations of the topology information
class TMREdge;
class TMRFace;
class TMRBlock;

/*
  Topological curve entity between two vertices

  The vertices themselves are extracted from the curve. Note that the
  vertices must be set in the curve.
*/
class TMREdge : public TMREntity {
 public:
  TMREdge( TMRCurve *_curve );
  TMREdge( TMRCurve *_curve, TMRVertex *_v1, TMRVertex *_v2 );
  ~TMREdge();

 private:
  // The underlying curve associated with this edge
  TMRCurve *curve;

  // The list of faces that refer to this edge
  class FaceList {
  public:
    TMRFace *surf;
    FaceList *next;
  } *faces;
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
  TMRFace( TMRSurface *_surface, TMREdge *_edges[] );
  ~TMRFace();

 private:
  // Pointers to the edges bounding this face
  TMRSurface *surface;
  TMREdge *edges[4];

  // The list of faces that refer to this edge
  class BlockList {
  public:
    TMRBlock *block;
    BlockList *next;
  } *blocks;
};

/*
  Topological volume entity bounded by six faces.
*/
class TMRBlock : public TMREntity {
 public:
  TMRBlock( TMRFace *_faces[] );
  ~TMRBlock();

 private:
  // Pointers to the faces bounding the volume
  TMRFace *faces[6];
};

/*
  The main topology class that contains the objects used to build the
  underlying mesh.
*/
class TMRTopology : public TMREntity {
 public:
  TMRTopology( MPI_Comm _comm );
  ~TMRTopology();

  void addBlock( TMRBlock *block );
  void buildTopology();

 private:
  MPI_Comm comm;

  // Set the block connectivity
  int num_blocks;
  TMRBlock *blocks;
};

#endif // TMR_TOPOLOGY_H
