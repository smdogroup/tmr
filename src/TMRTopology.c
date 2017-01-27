#include "TMRTopology.h"

/*
  Topological curve entity between two vertices
*/
TMREdge::TMREdge( TMRCurve *_curve ){
  curve = _curve;
  curve->incref();
}
 
TMREdge::TMREdge( TMRCurve *_curve, TMRVertex *v1, TMRVertex *v2 ){
  curve = _curve;
  curve->incref();
  curve->setVertices(v1, v2);
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
    edges[k] = _edges[k];
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

  num_faces = 0;
  faces = NULL;
  num_edges = 0;
  edges = NULL;
  num_nodes = 0;
  nodes = NULL;
}

/*
  Free the topology data
*/
TMRTopology::~TMRTopology(){
  // Free the faces, edges and nodes
  if (faces){
    for ( int i = 0; i < num_faces; i++ ){
      faces[i]->decref();
    }
    delete [] faces;
  }
  if (edges){
    for ( int i = 0; i < num_edges; i++ ){
      edges[i]->decref();
    }
    delete [] edges;
  }
  if (nodes){
    for ( int i = 0; i < num_nodes; i++ ){
      nodes[i]->decref();
    }
    delete [] nodes;
  }
}

/*
  Set all the edges in the topology
*/
void TMRTopology::setEdges( TMREdge **_edges, int _num_edges ){
  num_edges = _num_edges;
  edges = new TMREdge*[ num_edges ];
  for ( int i = 0; i < num_edges; i++ ){
    edges[i] = _edges[i];
    edges[i]->incref();
  }
}

/*
  Set all the faces into the mesh
*/
void TMRTopology::setFaces( TMRFace **_faces, int _num_faces ){
  num_faces = _num_faces;
  faces = new TMRFace*[ num_faces ];
  for ( int i = 0; i < num_faces; i++ ){
    faces[i] = _faces[i];
    faces[i]->incref();
  }
}

/*
  Build the topology which associates the nodes/edges/faces and blocks
  with node numbers, creates the underlying TMROctForest object.

  This can be used to generate the underlying mesh.
*/
void TMRTopology::buildConnectivity(){
  // First, uniquely order the vertices within the objects


  // Next uniquely order the edges



  
}

/*
  Retrieve the node associate with the given node number
 */
void TMRTopology::getNode( int node_num, TMRVertex **node ){
  *node = NULL;
  if (nodes && (node_num >= 0 && node_num < num_nodes)){
    *node = nodes[node_num];
  }
}

/*
  Retrieve the edge
*/
void TMRTopology::getEdge( int edge_num, 
                           TMREdge **edge, int *_n1, int *_n2 ){

}

/*
  Retrieve the face
*/
void TMRTopology::getFace( int face_num, TMRFace **face ){
  *face = NULL;
  if (faces && (face_num >= 0 && face_num < num_faces)){
    *face = faces[face_num];
  }
}


/*
  
 */
void TMRTopology::getConnectivity( int *nnodes, 
                                   int *nedges, int *nfaces,
                                   const int *face_nodes, 
                                   const int *edge_nodes ){


}
