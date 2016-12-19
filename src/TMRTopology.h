#ifndef TMR_TOPOLOGY_H
#define TMR_TOPOLOGY_H

#include "TMRBase.h"
#include "TMRTopology.h"

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
*/
class TMREdge : public TMREntity {
 public:
  TMREdge( TMRCurve *_curve );
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
  Topological surface entity bounded by four edges
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
  Topological volume entity bounded by six faces
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
  TMRTopology( MPI_Comm _comm ){}
  ~TMRTopology(){}

  void addBlock( TMRBlock *block );
  void buildTopology();

 private:
  int num_blocks;
  TMRBlock *blocks;
};

#endif // TMR_TOPOLOGY_H
