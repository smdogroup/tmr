#include "TMRTopology.h"
#include <stdio.h>

/*
  Topological curve entity between two vertices
*/
TMREdge::TMREdge( TMRCurve *_curve ){
  curve = _curve;
  curve->incref();
  curve->getVertices(&v1, &v2);
}

TMREdge::TMREdge( TMRCurve *_curve,
                  TMRVertex *_v1, TMRVertex *_v2 ){
  curve = _curve;
  curve->incref();
  v1 = _v1;
  v2 = _v2;
}

TMREdge::~TMREdge(){
  curve->incref();
}

/*m
  Retrive the curve
*/
void TMREdge::getCurve( TMRCurve **_curve ){
  *_curve = curve;
}

/*
  Retrieve the vertices
*/
void TMREdge::getVertices( TMRVertex **_v1, TMRVertex **_v2 ){
  if (_v1){ *_v1 = v1; }
  if (_v2){ *_v2 = v2; }
}

/*
  Topological surface entity bounded by four edges
*/
TMRFace::TMRFace( TMRSurface *_surface, 
                  TMREdge *_edges[], const int _edge_dir[] ){
  // Set the surface object
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
  Retrieve the underlying surface
*/
void TMRFace::getSurface( TMRSurface **surf ){
  *surf = surface;
}

/*
  Retrieve the connected edges and directions
*/
void TMRFace::getEdges( TMREdge ***_edges, const int **_edge_dir ){
  if (_edges){
    *_edges = edges;
  }
  if (_edge_dir){
    *_edge_dir = edge_dir;
  }
}

/*
  The main topology class that contains the objects used to build the
  underlying mesh.  
*/
TMRTopology::TMRTopology( TMRVertex **_vertices, int _num_vertices,
                          TMREdge **_edges, int _num_edges,
                          TMRFace **_faces, int _num_faces ){
  // Record the nodes
  num_vertices = _num_vertices;
  vertices = new TMRVertex*[ num_vertices ];
  for ( int i = 0; i < num_vertices; i++ ){
    vertices[i] = _vertices[i];
    vertices[i]->incref();
  }

  // Record the edges
  num_edges = _num_edges;
  edges = new TMREdge*[ num_edges ];
  for ( int i = 0; i < num_edges; i++ ){
    edges[i] = _edges[i];
    edges[i]->incref();
  }

  // Record the faces
  num_faces = _num_faces;
  faces = new TMRFace*[ num_faces ];
  for ( int i = 0; i < num_faces; i++ ){
    faces[i] = _faces[i];
    faces[i]->incref();
  }

  // Build the topology information based on the inputs
  face_to_vertices = new int[ 4*num_faces ];
  face_to_edges = new int[ 4*num_faces ];
  edge_to_vertices = new int[ 2*num_edges ];

  // Sort the addresses to make the array
  for ( int i = 0; i < num_faces; i++ ){
    TMREdge **e;
    faces[i]->getEdges(&e, NULL);

    // Search for the faces
    for ( int j = 0; j < 4; j++ ){
      for ( int k = 0; k < num_edges; k++ ){
        if (e[j] == edges[k]){
          face_to_edges[4*i+j] = k;
          break;
        }
      }
    }
  }

  // Sort the addresses to make the array 
  for ( int i = 0; i < num_edges; i++ ){
    TMRVertex *v1, *v2;
    edges[i]->getVertices(&v1, &v2);

    // Search for the vertices
    for ( int k = 0; k < num_vertices; k++ ){
      if (v1 == vertices[k]){
        edge_to_vertices[2*i] = k;
        break;
      }
    }

    // Search for the vertices
    for ( int k = 0; k < num_vertices; k++ ){
      if (v2 == vertices[k]){
        edge_to_vertices[2*i+1] = k;
        break;
      }
    }
  }

  // Get 
  for ( int i = 0; i < num_faces; i++ ){
    const int *edge_dir;
    faces[i]->getEdges(NULL, &edge_dir);

    int edge = face_to_edges[4*i];
    if (edge_dir[0] > 0){
      face_to_vertices[4*i] = edge_to_vertices[2*edge];
      face_to_vertices[4*i+2] = edge_to_vertices[2*edge+1];
    }
    else {
      face_to_vertices[4*i] = edge_to_vertices[2*edge+1];
      face_to_vertices[4*i+2] = edge_to_vertices[2*edge];
    }

    edge = face_to_edges[4*i+1];
    if (edge_dir[1] > 0){
      face_to_vertices[4*i+1] = edge_to_vertices[2*edge];
      face_to_vertices[4*i+3] = edge_to_vertices[2*edge+1];
    }
    else {
      face_to_vertices[4*i+1] = edge_to_vertices[2*edge+1];
      face_to_vertices[4*i+3] = edge_to_vertices[2*edge];
    }
  }
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
  if (vertices){
    for ( int i = 0; i < num_vertices; i++ ){
      vertices[i]->decref();
    }
    delete [] vertices;
  }

  delete [] face_to_vertices;
  delete [] face_to_edges;
  delete [] edge_to_vertices;
}

/*
  Retrieve the face object
*/
void TMRTopology::getSurface( int face_num, TMRSurface **face ){
  *face = NULL;
  if (faces && (face_num >= 0 && face_num < num_faces)){
    faces[face_num]->getSurface(face);
  }
}

/*
  Retrieve the curve object associated with the given face/edge index
*/
void TMRTopology::getFaceCurve( int face_num, int edge_index, TMRCurve **curve ){
  *curve = NULL;
  if (faces && (face_num >= 0 && face_num < num_faces)){
    TMREdge **edges;
    faces[face_num]->getEdges(&edges, NULL);
    if (edge_index >= 0 && edge_index < 4){
      edges[edge_index]->getCurve(curve);
    }
  }
}

/*
  Retrieve the connectivity information
*/
void TMRTopology::getConnectivity( int *nnodes, 
                                   int *nedges, int *nfaces,
                                   const int **face_nodes,
                                   const int **face_edges ){
  *nnodes = num_vertices;
  *nedges = num_edges;
  *nfaces = num_faces;
  *face_nodes = face_to_vertices;
  *face_edges = face_to_edges;
}
