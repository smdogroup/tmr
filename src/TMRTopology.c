#include "TMRTopology.h"

/*
  Topological curve entity between two vertices
*/
TMREdge::TMREdge( TMRCurve *_curve ){
  curve = _curve;
  curve->incref();
}
 
TMREdge::~TMREdge(){
  curve->incref();
}

/*
  Topological surface entity bounded by four edges
*/
TMRFace::TMRFace( TMRSurface *_surface, 
                  TMREdge *_edges[] ){
  // Set the surface object
  surface = _surface;
  surface->incref();

  for ( int k = 0; k < 4; k++ ){
    edges[k] = _edge[k];
    edges[k]->incref();
  }
}

TMRFace::~TMRFace(){
  surface->decref();
  for ( int k = 0; k < 4; k++ ){
    edges[k]->decref();
  }
}

/*
  Topological volume entity bounded by six faces
*/
TMRBlock::TMRBlock( TMRFace *_faces[] ){
  for ( int k = 0; k < 6; k++ ){
    faces[k] = _faces[k];
    faces[k]->incref();
  }
}
 
TMRBlock::~TMRBlock(){
  for ( int k = 0; k < 6; k++ ){
    faces[k]->decref();
  }
}

/*
  The main topology class that contains the objects used to build the
  underlying mesh.  
*/
TMRTopology::TMRTopology( MPI_Comm _comm ){
  comm = _comm;

}

TMRTopology::~TMRTopology(){

}

void TMRTopology::addBlock( TMRBlock *block ){


}

/*
  Build the topology which associates the nodes/edges/faces and blocks
  with node numbers, creates the underlying TMROctForest object.

  This can be used to generate the underlying mesh.
*/
void TMRTopology::buildTopology(){
  

}
