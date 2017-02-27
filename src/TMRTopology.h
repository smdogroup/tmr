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
  The abstract edge class
*/
class TMREdge : public TMREntity {
 public:
  
  // Parametrize the curve on the given surface
  virtual int getParamsOnSurface( TMRSurface *surface, double t, 
                                  int dir, double *u, double *v );

  // Set/retrive the vertices at the beginning and end of the curve
  void setVertices( TMRVertex *_v1, TMRVertex *_v2 );
  void getVertices( TMRVertex **_v1, TMRVertex **_v2 );

  // Integrate along the edge and return an array containing
  // the parametric locations to provide an even spacing
  double integrate( double t1, double t2, double tol,
                    double **tvals, double **dist, int *nvals );

  // Set/retrieve the mesh
  void setMesh( TMRCurveMesh *_mesh );
  void getMesh( TMRCurveMesh **_mesh );

 private:
  // The start/end vertices of the curve
  TMRVertex *v1, *v2;

  // The mesh for the curve - if it exists
  TMRCurveMesh *mesh;

};

/*
  The edge loop class
*/
class TMREdgeLoop : public TMREntity {
 public:

 private:

};



/*
  The abstract face class
*/
class TMRFace : public TMREntity {
 public:
  
  // Set/retrieve the mesh
  void setMesh( TMRSurfaceMesh *_mesh );
  void getMesh( TMRSurfaceMesh **_mesh );

 private:

  // The mesh for the curve - if it exists
  TMRSurfaceMesh *mesh;

};


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
