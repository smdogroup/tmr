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
class TMRVolume;
class TMREdgeMesh;
class TMRFaceMesh;
class TMRVolumeMesh;

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

  // Is this edge degenerate
  virtual int isDegenerate(){ return 0; }

  // Set/retrieve the vertices at the beginning and end of the curve
  void setVertices( TMRVertex *_v1, TMRVertex *_v2 );
  void getVertices( TMRVertex **_v1, TMRVertex **_v2 );

  // Set/retrieve the source edge
  void setSource( TMREdge *_edge );
  void getSource( TMREdge **_edge );

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
  TMREdge *source; // Source edge (may be NULL)

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

  // Retrieve the edges in the loop and their orientations
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
  TMRFace( int _normal_orient=1 );
  virtual ~TMRFace();

  // Get the underlying normal direction
  int getOrientation();

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

  // Set/retrieve the source face
  void setSource( TMRVolume *_volume, TMRFace *_face );
  void getSource( int *dir, TMRVolume **_volume, TMRFace **_face );

  // Set/retrieve the mesh
  void setMesh( TMRFaceMesh *_mesh );
  void getMesh( TMRFaceMesh **_mesh );

  // Write the object to the VTK file
  void writeToVTK( const char *filename );

 private:
  // Relative orientation of the normal direction and
  // the parametric normal direction
  const int normal_orient;

  // The mesh for the curve - if it exists
  TMRFaceMesh *mesh;

  int source_dir; // The relative source face direction
  TMRVolume *source_volume; // Source volume
  TMRFace *source; // Source face (may be NULL)

  // Store the loop information
  int num_loops, max_num_loops;
  TMREdgeLoop **loops;

  // Derivative step size
  static double deriv_step_size;
};

/*
  The TMR volume object.

  This object is used primarily as a container class and does not
  neccessarily parametrize the interior of a volume, except in some
  instances - for instance when using transfinite interpolation from
  surfaces.
*/
class TMRVolume : public TMREntity {
 public:
  TMRVolume( int _nfaces, TMRFace **_faces, const int *_dir=NULL );
  virtual ~TMRVolume();

  // Get the parameter range for this volume
  virtual void getRange( double *umin, double *vmin, double *wmin,
                         double *umax, double *vmax, double *wmax );

  // Given the parametric point u,v,w compute the physical location x,y,z
  virtual int evalPoint( double u, double v, double w, TMRPoint *X );

  // Get the faces that enclose this volume
  void getFaces( int *_num_faces, TMRFace ***_faces, const int **_dir );

  // Set/retrieve the mesh
  void setMesh( TMRVolumeMesh *_mesh );
  void getMesh( TMRVolumeMesh **_mesh );

  // Write the object to the VTK file
  void writeToVTK( const char *filename );

 private:
  // Store the face information
  int num_faces;
  TMRFace **faces;
  int *dir;

  // Set the volume mesh
  TMRVolumeMesh *mesh;
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
            int _num_faces, TMRFace **_faces,
            int _num_volumes=0, TMRVolume **_volumes=NULL );
  ~TMRModel();

  // Retrieve the underlying vertices/curves/surfaces
  void getVertices( int *_num_vertices, TMRVertex ***_vertices );
  void getEdges( int *_num_edges, TMREdge ***_edges );
  void getFaces( int *_num_faces, TMRFace ***_faces );
  void getVolumes( int *_num_volumes, TMRVolume ***_volumes );

  // Query geometric objects based on pointer values
  int getVertexIndex( TMRVertex *vertex );
  int getEdgeIndex( TMREdge *edge );
  int getFaceIndex( TMRFace *face );
  int getVolumeIndex( TMRVolume *volume );

 private:
  // Verify that everything is more or less well defined. Print out
  // error messages if something doesn't make sense.
  int verify();

  // The verticies, curves and surfaces that define a BRep
  int num_vertices, num_edges, num_faces, num_volumes;
  TMRVertex **vertices;
  TMREdge **edges;
  TMRFace **faces;
  TMRVolume **volumes;

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
  OrderedPair<TMRVolume> *ordered_volumes;
};

/*
  The main topology class that contains the objects used to build the
  underlying mesh.
  
  This class takes in a general TMRModel, but there are additional
  requirements that are placed on the model to create a proper
  topology object. These requirements are as follows:

  1. No edge can degenerate to a vertex.
  2. No face can degenerate to an edge or vertex.
  3. All faces must be surrounded by a single edge loop with 4
  non-degenerate edges.
  4. All volumes must contain 6 non-degenerate faces that are
  ordered in coordinate ordering as shown below. Furthermore, all
  volumes must be of type TMRTFIVolume.

  Coordinate ordering for the faces of a volume:

       6 ---------------------- 7
      /|                       /|
     / |                      / |
    /  |                 5 --/  |
   /   |                /   /   |
  4 ---------------------- 5    |
  |    |                   |    |
  |    |                   |3 --|
  |    |                   ||   |
  |    2 ------------------|--- 3
  | 0 /                    | 1 /
  |/|/                     |/|/
  | /--4               2 --| /
  |/  /                |   |/
  0 ---------------------- 1
  
  w
  |  / v
  | /
  |/
  . --- u
*/
class TMRTopology : public TMREntity {
 public:
  TMRTopology( MPI_Comm _comm, TMRModel *geo );
  ~TMRTopology();
  
  // Retrieve the face/edge/node information
  int getNumVolumes();
  int getNumFaces();
  int getNumEdges();
  int getNumVertices();
  void getVolume( int vol_num, TMRVolume **volume );
  void getFace( int face_num, TMRFace **face );
  void getEdge( int edge_num, TMREdge **edge );
  void getVertex( int vertex_num, TMRVertex **vertex );

  // Retrive the quadtree connectivity from the topology object
  void getConnectivity( int *nnodes, int *nedges, int *nfaces,
                        const int **face_nodes, const int **face_edges );

  // Retrieve the octree connectivity from the topology object
  void getConnectivity( int *nnodes, int *nedges, int *nfaces, int *nvolumes,
                        const int **volume_nodes, const int **volume_edges,
                        const int **volume_faces );

 private:
  // Compute the face connectivity
  void computeConnectivty( int num_entities, int num_edges, 
                           int num_faces, const int ftoedges[],
                           int **_face_to_face_ptr,
                           int **_face_to_face );
  void reorderEntities( int num_entities, int num_edges, int num_faces,
                        const int *ftoedges, int *entity_to_new_num, 
                        int *new_num_to_entity, int use_rcm=1 );

  // Connectivity for face -> edge, face -> vertex and edge -> vertex
  void computeFaceConn();

  // Connectivity for volume -> face/edge/vertex
  void computeVolumeConn();

  // Get the MPI communicator
  MPI_Comm comm;

  // The connectivity information
  int *edge_to_vertices;
  int *face_to_edges;
  int *face_to_vertices;
  int *volume_to_vertices;
  int *volume_to_edges;
  int *volume_to_faces;

  // The reordering for the faces
  int *face_to_new_num;
  int *new_num_to_face;

  // Reordering of the volumes
  int *volume_to_new_num;
  int *new_num_to_volume;

  // The geometry class
  TMRModel *geo;
};

#endif // TMR_TOPOLOGY_H
