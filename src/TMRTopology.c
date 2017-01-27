#include "TMRTopology.h"

/*
  Topological curve entity between two vertices
*/
TMREdge::TMREdge( MPI_Comm _comm, TMRCurve *_curve,
                  TMRVertex *v1, TMRVertex *v2 ){
  comm = _comm;
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
TMRFace::TMRFace( MPI_Comm _comm, TMRSurface *_surface, 
                  TMREdge *_edges[], int _edge_dir[] ){
  // Set the surface object
  comm = _comm;
  surface = _surface;
  surface->incref();

  for ( int k = 0; k < 4; k++ ){
    edges[k] = _edges[k];
    edges[k]->incref();
    edge_dir[k] = _edge_dir[k];
  }
}

TMRFace::~TMRFace(){
  surface->decref();
  for ( int k = 0; k < 4; k++ ){
    edges[k]->decref();
  }
}

/*
  The main topology class that contains the objects used to build the
  underlying mesh.  
*/
TMRTopology::TMRTopology( MPI_Comm _comm,
                          TMRFace **_faces, int _num_faces ){
  comm = _comm;

  // Record the faces
  num_faces = _num_faces;
  faces = new TMRFace*[ num_faces ];
  for ( int i = 0; i < num_faces; i++ ){
    faces[i] = _faces[i];
    faces[i]->incref();
  }

  // Keep track of the underlying edges/nodes still to be counted
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
  Build the topology which associates the nodes/edges/faces and blocks
  with node numbers, creates the underlying TMROctForest object.

  This can be used to generate the underlying mesh.
*/
void TMRTopology::buildConnectivity(){
  // Establish an ordering for the edges
  TMREdge **tmp_edges = new TMREdge*[ 4*num_faces ];

  // Loop over all of the edges and compute 


  
}

/*
  Retrieve the face object
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
                                   const int **face_nodes, 
                                   const int **edge_nodes ){
  
}
