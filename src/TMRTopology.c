#include "TMRTopology.h"
#include <stdio.h>

/*
  The main topology class that contains the objects used to build the
  underlying mesh.  
*/
TMRTopology::TMRTopology( TMRGeometry *_geo ){
  // Increase the ref. count to the geometry object
  geo = _geo;
  geo->incref();

  // Get the geometry objects
  int num_vertices, num_edges, num_faces;
  TMRVertex **vertices;
  TMRCurve **edges;
  TMRSurface **faces;
  geo->getVertices(&num_vertices, &vertices);
  geo->getCurves(&num_edges, &edges);
  geo->getSurfaces(&num_faces, &faces);

  // Build the topology information based on the inputs
  face_to_vertices = new int[ 4*num_faces ];
  face_to_edges = new int[ 4*num_faces ];
  edge_to_vertices = new int[ 2*num_edges ];

  // Sort the addresses to make the array
  const int coordinate_to_edge[] = {3, 1, 0, 2};

  // Create the face -> edge information
  for ( int i = 0; i < num_faces; i++ ){
    TMRCurve **e;
    faces[i]->getCurveSegment(0, NULL, &e, NULL);

    // Search for the face indices
    for ( int jp = 0; jp < 4; jp++ ){
      int j = coordinate_to_edge[jp];
      face_to_edges[4*i+j] = geo->getCurveIndex(e[j]);
    }
  }

  // Create the edge -> vertex information
  for ( int i = 0; i < num_edges; i++ ){
    TMRVertex *v1, *v2;
    edges[i]->getVertices(&v1, &v2);
    
    edge_to_vertices[2*i] = geo->getVertexIndex(v1);
    edge_to_vertices[2*i+1] = geo->getVertexIndex(v2);
  }

  // Create the face -> vertex information
  for ( int i = 0; i < num_faces; i++ ){
    const int *edge_dir;
    faces[i]->getCurveSegment(0, NULL, NULL, &edge_dir);

    int edge = face_to_edges[4*i];
    if (edge_dir[3] > 0){
      face_to_vertices[4*i] = edge_to_vertices[2*edge+1];
      face_to_vertices[4*i+2] = edge_to_vertices[2*edge];
    }
    else {
      face_to_vertices[4*i] = edge_to_vertices[2*edge];
      face_to_vertices[4*i+2] = edge_to_vertices[2*edge+1];
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
  // Free the geometry object
  geo->decref();
  delete [] face_to_vertices;
  delete [] face_to_edges;
  delete [] edge_to_vertices;
}

/*
  Retrieve the face object
*/
void TMRTopology::getSurface( int face_num, TMRSurface **face ){
  int num_faces;
  TMRSurface **faces;
  geo->getSurfaces(&num_faces, &faces);
  *face = NULL;
  if (faces && (face_num >= 0 && face_num < num_faces)){
    *face = faces[face_num];
  }
}

/*
  Retrieve the curve object associated with the given face/edge index
*/
void TMRTopology::getFaceCurve( int face_num, int edge_index, TMRCurve **curve ){
  int num_faces;
  TMRSurface **faces;
  geo->getSurfaces(&num_faces, &faces);

  // Coordinate edge number to actual edge number
  const int coordinate_to_edge[] = {3, 1, 0, 2};

  *curve = NULL;
  if (faces && (face_num >= 0 && face_num < num_faces)){
    TMRCurve **edges;
    faces[face_num]->getCurveSegment(0, NULL, &edges, NULL);
    if (edge_index >= 0 && edge_index < 4){
      edge_index = coordinate_to_edge[edge_index];
      *curve = edges[edge_index];
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
  int num_vertices, num_edges, num_faces;
  geo->getVertices(&num_vertices, NULL);
  geo->getCurves(&num_edges, NULL);
  geo->getSurfaces(&num_faces, NULL);
  *nnodes = num_vertices;
  *nedges = num_edges;
  *nfaces = num_faces;
  *face_nodes = face_to_vertices;
  *face_edges = face_to_edges;
}
