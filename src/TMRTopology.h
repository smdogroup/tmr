#ifndef TMR_TOPOLOGY_H
#define TMR_TOPOLOGY_H

/*
  The following header file includes the definitions of the
  topological entries within a multiblock mesh.

  These entries refer to the underlying geometry definitions and
  provide connectivity for edges, faces and blocks within a
  multi-block mesh.
*/

#include "TMRBase.h"
#include "TMRGeometry.h"

class TMREdge;
class TMRFace;
class TMREdgeMesh;
class TMRFaceMesh;

/*
  The vertex class: Note that this is used to store both the
  point and to represent the underlying geometry
*/
class TMRVertex : public TMREntity {
 public:
  TMRVertex(){ var = -1; }
  virtual ~TMRVertex(){}

  // Evalue the point
  virtual int evalPoint( TMRPoint *p ) = 0;

  // Get the parameters on the associated curve/surface
  virtual int getParamOnEdge( TMREdge *edge, double *t );
  virtual int getParamsOnFace( TMRFace *face,
                               double *u, double *v );

  // Set/retrieve the node numbers
  int setNodeNum( int *num );
  int getNodeNum( int *num );

 private:
  int var;
};

/*
  The abstract edge class
*/
class TMREdge : public TMREntity {
 public:
  TMREdge();
  virtual ~TMREdge();

  // Get the parameter range for this surface
  virtual void getRange( double *tmin, double *tmax ) = 0;
 
  // Given the parametric point, compute the x,y,z location
  virtual int evalPoint( double t, TMRPoint *X ) = 0;
  
  // Perform the inverse evaluation
  virtual int invEvalPoint( TMRPoint p, double *t );

  // Given the parametric point, evaluate the first derivative 
  virtual int evalDeriv( double t, TMRPoint *Xt );
  
  // Parametrize the curve on the given surface
  virtual int getParamsOnFace( TMRFace *face, double t, 
                               int dir, double *u, double *v );

  // Set/retrive the vertices at the beginning and end of the curve
  void setVertices( TMRVertex *_v1, TMRVertex *_v2 );
  void getVertices( TMRVertex **_v1, TMRVertex **_v2 );

  // Integrate along the edge and return an array containing
  // the parametric locations to provide an even spacing
  double integrate( double t1, double t2, double tol,
                    double **tvals, double **dist, int *nvals );

  // Set/retrieve the mesh
  void setMesh( TMREdgeMesh *_mesh );
  void getMesh( TMREdgeMesh **_mesh );

  // Write the object to the VTK file
  void writeToVTK( const char *filename );
 private:
  // The start/end vertices of the curve
  TMRVertex *v1, *v2;

  // The mesh for the curve - if it exists
  TMREdgeMesh *mesh;

  // Derivative step size
  static double deriv_step_size;
};

/*
  The edge loop class

  The curve segments must form a closed loop which is checked 
  by the code. The boundary must lie counterclockwise around a 
  surface while holes/cutouts must run clockwise so that the 
  domain always lies to the left of the curve.
*/
class TMREdgeLoop : public TMREntity {
 public:
  TMREdgeLoop( int _ncurves, TMREdge *_edges[], 
               const int _dir[] );
  ~TMREdgeLoop();

  void getEdgeLoop( int *_ncurves, TMREdge **_edges[], 
                    const int *_dir[] );

 private:
  int nedges;
  TMREdge **edges;
  int *dir;
};

/*
  The abstract face class
*/
class TMRFace : public TMREntity {
 public:
  TMRFace();
  virtual ~TMRFace();

  // Get the parameter range for this surface
  virtual void getRange( double *umin, double *vmin,
                         double *umax, double *vmax ) = 0;
 
  // Given the parametric point, compute the x,y,z location
  virtual int evalPoint( double u, double v, TMRPoint *X ) = 0;
  
  // Perform the inverse evaluation
  virtual int invEvalPoint( TMRPoint p, double *u, double *v );

  // Given the parametric point, evaluate the first derivative 
  virtual int evalDeriv( double u, double v, 
                         TMRPoint *Xu, TMRPoint *Xv );

  // Add an edge loop to the face
  int getNumEdgeLoops();
  void addEdgeLoop( TMREdgeLoop *loop );
  void getEdgeLoop( int k, TMREdgeLoop **loop );

  // Set/retrieve the mesh
  void setMesh( TMRFaceMesh *_mesh );
  void getMesh( TMRFaceMesh **_mesh );

  // Write the object to the VTK file
  void writeToVTK( const char *filename );

 private:
  // The mesh for the curve - if it exists
  TMRFaceMesh *mesh;

  // Store the loop information
  int num_loops, max_num_loops;
  TMREdgeLoop **loops;

  // Derivative step size
  static double deriv_step_size;
};

/*
  The TMR Geometry class. 

  This contains the geometry objects -- vertices, curves and surfaces
  -- that are used to define the geometry of a model.
*/
class TMRModel : public TMREntity {
 public:
  TMRModel( int _num_vertices, TMRVertex **_vertices, 
            int _num_edges, TMREdge **_edges,
            int _num_faces, TMRFace **_faces );
  ~TMRModel();

  // Retrieve the underlying vertices/curves/surfaces
  void getVertices( int *_num_vertices, TMRVertex ***_vertices );
  void getEdges( int *_num_edges, TMREdge ***_edges );
  void getFaces( int *_num_faces, TMRFace ***_faces );

  // Query geometric objects based on pointer values
  int getVertexIndex( TMRVertex *vertex );
  int getEdgeIndex( TMREdge *edge );
  int getFaceIndex( TMRFace *face );

 private:
  // Verify that everything is more or less well defined. Print
  // out error messages if something doesn't make sense.
  int verify();

  // The verticies, curves and surfaces that define a BRep
  int num_vertices, num_edges, num_faces;
  TMRVertex **vertices;
  TMREdge **edges;
  TMRFace **faces;

  // This keeps track of the ordering of the geometry objects and
  // enables a fast object -> object index lookup
  template <class ctype>
  class OrderedPair {
  public:
    int num;
    ctype *obj;
  };

  template <class ctype>
  static int compare_ordered_pairs( const void *avoid, const void *bvoid );

  OrderedPair<TMRVertex> *ordered_verts;
  OrderedPair<TMREdge> *ordered_edges;
  OrderedPair<TMRFace> *ordered_faces;
};

/*
  The main topology class that contains the objects used to build the
  underlying mesh.
*/
class TMRTopology : public TMREntity {
 public:
  TMRTopology( TMRModel *geo );
  ~TMRTopology();
  
  // Retrieve the face/edge/node information
  void getSurface( int face_num, TMRFace **face );
  void getFaceCurve( int face_num, int edge_index, TMREdge **curve );

  // Retrive the connectivity from the topology object
  void getConnectivity( int *nnodes, int *nedges, int *nfaces,
                        const int **face_nodes, const int **face_edges );

 private:
  // The connectivity information
  int *edge_to_vertices;
  int *face_to_edges;
  int *face_to_vertices;

  // The geometry class
  TMRModel *geo;
};

#endif // TMR_TOPOLOGY_H
