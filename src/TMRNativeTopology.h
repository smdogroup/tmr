#ifndef TMR_NATIVE_TOPOLOGY_H
#define TMR_NATIVE_TOPOLOGY_H

#include "TMRTopology.h"

/*
  TMREdge from curve class
*/
class TMREdgeFromCurve : public TMREdge {
 public:
  TMREdgeFromCurve( TMRCurve *_curve ){
    curve = _curve;
    curve->incref();
  }
  ~TMREdgeFromCurve(){
    curve->decref();
  }
  void getRange( double *tmin, double *tmax ){
    curve->getRange(tmin, tmax);
  } 
  int evalPoint( double t, TMRPoint *X ){
    return curve->evalPoint(t, X);
  } 
  int invEvalPoint( TMRPoint p, double *t ){
    return curve->invEvalPoint(p, t);
  }
  int evalDeriv( double t, TMRPoint *Xt ){
    return curve->evalDeriv(t, Xt);
  }
 private:
  TMRCurve *curve;
};

/*
  TMRFace object from a surface class
*/
class TMRFaceFromSurface : public TMRFace {
 public:
  TMRFaceFromSurface( TMRSurface *_surf ){
    surf = _surf;
    surf->incref();
  }
  ~TMRFaceFromSurface(){
    surf->decref();
  }
  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax ){
    surf->getRange(umin, vmin, umax, vmax);
  } 
  int evalPoint( double u, double v, TMRPoint *X ){
    return surf->evalPoint(u, v, X);
  } 
  int invEvalPoint( TMRPoint p, double *u, double *v ){
    return surf->invEvalPoint(p, u, v);
  }
  int evalDeriv( double u, double v, 
                 TMRPoint *Xu, TMRPoint *Xv ){
    return surf->evalDeriv(u, v, Xu, Xv);
  }
 private:
  TMRSurface *surf;
};

/*
  Set the TMRVertex from a point
*/
class TMRVertexFromPoint : public TMRVertex {
 public:
  TMRVertexFromPoint( TMRPoint p );
  ~TMRVertexFromPoint(){}
  int evalPoint( TMRPoint *p );
 private:
  TMRPoint pt;
};

/*
  Set the TMRVertex location based on a parametric location along a
  curve.  

  This takes either a parametric point or does an inverse evaluation
  first to determine the parametric location.
*/
class TMRVertexFromEdge : public TMRVertex {
 public:
  TMRVertexFromEdge( TMREdge *_edge, double _t );
  TMRVertexFromEdge( TMREdge *_edge, TMRPoint p );
  ~TMRVertexFromEdge();
  int evalPoint( TMRPoint *p );
  int getParamOnEdge( TMREdge *edge, double *t );
  int getParamsOnFace( TMRFace *face, double *u, double *v );
  TMREdge* getEdge();

 private:
  double t;
  TMREdge *edge;
};

/*
  Evaluate a vertex location based on its parametric location on
  a surface.
*/
class TMRVertexFromFace : public TMRVertex {
 public:
  TMRVertexFromFace( TMRFace *_face, double _u, double _v );
  TMRVertexFromFace( TMRFace *_face, TMRPoint p );
  ~TMRVertexFromFace();
  int evalPoint( TMRPoint *p );
  int getParamsOnFace( TMRFace *surf,
                       double *u, double *v );
  TMRFace* getFace();

 private:
  double u, v;
  TMRFace *face;
};

/*
  The curve parametrized on the surface C(t) = S(u(t), v(t))
*/
class TMREdgeFromFace : public TMREdge {
 public:
  TMREdgeFromFace( TMRFace *_face, TMRPcurve *_pcurve, int _is_degen=0 );
  ~TMREdgeFromFace();
  void getRange( double *tmin, double *tmax );
  int evalPoint( double t, TMRPoint *X );
  int getParamsOnFace( TMRFace *face, double t, 
                       int dir, double *u, double *v );
  int invEvalPoint( TMRPoint X, double *t );
  int evalDeriv( double t, TMRPoint *Xt );
  int isDegenerate(){ return is_degen; }
  void addEdgeFromFace( TMRFace *_face, TMRPcurve *_pcurve );
 private:
  int is_degen;
  int nfaces;
  TMRFace **faces;
  TMRPcurve **pcurves;
};

/*
  Split the curve
*/
class TMRSplitEdge : public TMREdge {
 public:
  TMRSplitEdge( TMREdge *_edge, double _t1, double _t2 );
  TMRSplitEdge( TMREdge *_edge, TMRPoint *p1, TMRPoint *p2 );
  TMRSplitEdge( TMREdge *_edge, TMRVertex *p1, TMRVertex *p2 );
  TMRSplitEdge();

  void getRange( double *tmin, double *tmax );
  int evalPoint( double t, TMRPoint *X );
  int getParamsOnFace( TMRFace *face, double t, 
                       int dir, double *u, double *v );

 private:
  // Evaluate the bspline curve
  double t1, t2;
  TMREdge *edge;
};

/*
  Transfinite-interpolation edge

  This is just a line between the vertices v1 and v2
*/
class TMRTFIEdge : public TMREdge {
 public:
  TMRTFIEdge( TMRVertex *_v1, TMRVertex *_v2 );
  ~TMRTFIEdge();
  void getRange( double *tmin, double *tmax );
  int evalPoint( double t, TMRPoint *X );
};

/*
  Transfinite-interpolation surface

  v3-------e2-------v2
  |                 |
  |      v ^        |
  |        |        |
  e3       ----> u  e1
  |                 |
  |                 |
  |                 |
  v0-------e0-------v1
*/
class TMRTFIFace : public TMRFace {
 public:
  TMRTFIFace( TMREdge *_edges[], const int dir[],
              TMRVertex *verts[] );
  ~TMRTFIFace();

  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax ); 
  int evalPoint( double u, double v, TMRPoint *X ); 
  int invEvalPoint( TMRPoint p, double *u, double *v );
  int evalDeriv( double u, double v, 
                 TMRPoint *Xu, TMRPoint *Xv );

 private:
  // The edges within the transfinite-interpolation
  TMREdge *edges[4];
  int dir[4];

  // The vertex locations
  TMRPoint c[4];
};

/*
  Parametric TFI class

  This face class defines a segment of a surface defined in parameter
  space. This class is used to create a surface-conforming patch. This
  class is used by TMRMesh to create a surface topology object needed
  for TMRQuadForest.

  The parametric curves and the input vertices are used to define a
  segment on the surface in parameter space that looks like this:

  v3-------e2-------v2
  |                 |
  |      v ^        |
  |        |        |
  e3       ----> u  e1
  |                 |
  |                 |
  |                 |
  v0-------e0-------v1

  Note that all parameter curves must have the range [0,
  1]. Furthermore, note that the constructor for this object
  automatically sets the vertices into the curves and adds the curve
  segments (using the TMRSurface::addEdgeLoop) for convenience.
*/
class TMRParametricTFIFace : public TMRFace {
 public:
  TMRParametricTFIFace( TMRFace *_face, 
                        TMREdge *_edges[], const int dir[],
                        TMRVertex *verts[] );
  ~TMRParametricTFIFace();
  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax ); 
  int evalPoint( double u, double v, TMRPoint *X ); 
  int invEvalPoint( TMRPoint p, double *u, double *v );
  int evalDeriv( double u, double v, 
                 TMRPoint *Xu, TMRPoint *Xv );

 private:
  TMRFace *face;
  TMREdge *edges[4];
  int dir[4];

  // Parametric coordinates of the vertex points on the surface
  double vupt[4], vvpt[4];
};

/*
  The following class performs a transfinite interpolation over a
  brick volume element. The arguments to the volume are arranged in a
  coordinate ordering. This ordering is designed to be easily
  compatible with the TMROctForest ordering. 

  The input faces, edges, and vertices are ordered in the following
  manner:
  
       v6----------e3----------v7
      /|                       /|
     e6|                      e7|
    /  |                 f5--/  |
   /   |                /   /   |
  v4----------e2-----------v5   |
  |   e11                  |   e10
  |    |                   |f3--|
  |    |                   ||   |
  e8   v2---------e1-------e9--v3
  |f0 /                    |f1 /
  |/|e4                    |/|e5
  | /--f4              f2--| /
  |/  /                |   |/
  v0----------e0----------v1
  
  The orient argument captures the root node information for each face
  such that faces that have different orientations will enclose a
  properly oriented volume. 
*/
class TMRTFIVolume : public TMRVolume {
 public:
  TMRTFIVolume( TMRFace *_faces[], const int orient[],
                TMREdge *_edges[], const int dir[],
                TMRVertex *verts[] );
  ~TMRTFIVolume();
  void getRange( double *umin, double *vmin, double *wmin,
                 double *umax, double *vmax, double *wmax );
  int evalPoint( double u, double v, double w, TMRPoint *X ); 
  void getEntities( TMRFace ***_faces, TMREdge ***_edges, 
                    TMRVertex ***_verts );

 private:
  // Faces surrounding the volume: coordinate ordered
  TMRFace *faces[6];
  int orient[6];

  // Edges defining the volume intersections
  TMREdge *edges[12];
  int edge_dir[12];
  
  // Locations of the vertices/corners
  TMRVertex *verts[8];
  TMRPoint c[8];
};

#endif // TMR_NATIVE_TOPOLOGY_H
