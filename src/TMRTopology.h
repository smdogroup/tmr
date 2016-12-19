#ifndef TMR_TOPOLOGY_H
#define TMR_TOPOLOGY_H

#include "TMRBase.h"
#include "TMRTopology.h"

/*
  The following header file includes the definitions of the
  topological entries within a block mesh.

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
  TMREdge( TMRVertex *_v1, TMRVertex *_v2 ){
    v1 = _v1;  v1->incref();
    v2 = _v2;  v2->incref();
  }
  ~TMREdge(){
    v1->decref();
    v2->decref();
  }

  // Associate the edge to the curve
  void toCurve( TMRCurve *_curve ){
    curve = _curve;
    curve->incref();
  }

 private:
  // The vertices
  TMRVertex *v1, *v2;
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
  TMRFace( TMREdge *_edges[] ){
    for ( int k = 0; k < 4; k++ ){
      edges[k] = _edge[k];
      edges[k]->incref();
    }
  }
  ~TMRFace(){
    for ( int k = 0; k < 4; k++ ){
      edges[k]->decref();
    }
  }

  // Associate the face with a surface
  void toSurface( TMREdge *_face ){
    face = _face;
    face->incref();
  }

 private:
  // Pointers to the edges bounding this face
  TMREdge *edges[4];
  TMRSurface *surf;

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
  TMRBlock( TMRFace *_faces[] ){
    for ( int k = 0; k < 6; k++ ){
      faces[k] = _faces[k];
      faces[k]->incref();
    }
  }
  ~TMRBlock(){
    for ( int k = 0; k < 6; k++ ){
      faces[k]->decref();
    }
  }

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

  // Add all of the blocks to the mesh
  void addBlocks( TMRBlock *blocks, int num_blocks );

  // Build the underlying topology
  void buildForest();

  
  
};

#endif // TMR_TOPOLOGY_H
