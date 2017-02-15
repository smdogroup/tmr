#ifndef TMR_GEOMETRY_H
#define TMR_GEOMETRY_H

#include "TMRBase.h"

/*
  The following header file contains the interface for the geometry/
  topology for the TMR objects. These vertex/edge/surface and volume
  objects are used to map the 

  The vertex, edge, face and volume classes are used in conjunction
  with the TMROct(Quad)Forest class to evaluate nodal locations with
  the mesh. These interfaces are designed to be overriden with an
  external geometry engine.
*/

/*
  The vertex class: Note that this is used to store both the
  point and to represent the underlying geometry
*/
class TMRVertex : public TMREntity {
 public:
  TMRVertex(){}
  virtual ~TMRVertex(){}

  // Evalue the point
  virtual int evalPoint( TMRPoint *p ) = 0;
};

/*
  The parametrization for a curve
*/
class TMRCurve : public TMREntity {
 public:
  TMRCurve();
  TMRCurve( TMRVertex *v1, TMRVertex *v2 );
  virtual ~TMRCurve();

  // Get the parameter range for this edge
  virtual void getRange( double *tmin, double *tmax ) = 0;
  
  // Given the parametric point, evaluate the x,y,z location
  virtual int evalPoint( double t, TMRPoint *X ) = 0;

  // Given the point, find the parametric location
  virtual int invEvalPoint( TMRPoint X, double *t );

  // Given the parametric point, evaluate the derivative 
  virtual int evalDeriv( double t, TMRPoint *Xt );
  
  // Set/retrive the vertices at the beginning and end of the curve
  void setVertices( TMRVertex *_v1, TMRVertex *_v2 );
  void getVertices( TMRVertex **_v1, TMRVertex **_v2 );

  // Integrate along the edge and return an array containing
  // the parametric locations to provide an even spacing
  double integrate( double t1, double t2, double tol,
                    double **tvals, double **dist, int *nvals );

  // Write the object to the VTK file
  void writeToVTK( const char *filename );
  
 private:
  // Set the step size
  static double deriv_step_size;

  // The start/end vertices of the curve
  TMRVertex *v1, *v2;
};

/*
  The parametrization of a surface
*/
class TMRSurface : public TMREntity {
 public:
  TMRSurface();
  virtual ~TMRSurface();

  // Get the parameter range for this surface
  virtual void getRange( double *umin, double *vmin,
                         double *umax, double *vmax ) = 0;
 
  // Given the parametric point, compute the x,y,z location
  virtual int evalPoint( double u, double v, TMRPoint *X ) = 0;
  
  // Perform the inverse evaluation
  virtual int invEvalPoint( TMRPoint p, double *u, double *v ) = 0;

  // Given the parametric point, evaluate the first derivative 
  virtual int evalDeriv( double u, double v, 
                         TMRPoint *Xu, TMRPoint *Xv ) = 0;

  // Add a curve segment. The curve segments must form a closed
  // loop which is checked by the code. The boundary must lie
  // counterclockwise around the surface while holes/cutouts
  // must run clockwise so that the domain always lies to the
  // left of the curve.
  int addCurveSegment( TMRCurve **_curves, 
                       const int _dir[], int ncurves );
  void getCurves( TMRCurve ***_curves, 
                  const int **_dir, int *_num_curves );

  // Write the object to the VTK file
  void writeToVTK( const char *filename );

 private:
  // Pointers to the curves that enclose the object. 
  // Note  to the list of curves
  int num_curves;
  TMRCurve **curves;
  int *dir;

  // Set the step size
  static double deriv_step_size;
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
class TMRVertexFromCurve : public TMRVertex {
 public:
  TMRVertexFromCurve( TMRCurve *_curve, double _t );
  TMRVertexFromCurve( TMRCurve *_curve, TMRPoint p );
  ~TMRVertexFromCurve();
  int evalPoint( TMRPoint *p );
  TMRCurve* getCurve();
  double getParamPoint();

 private:
  double t;
  TMRCurve *curve;
};

/*
  Evaluate a vertex location based on its parametric location on
  a surface.
*/
class TMRVertexFromSurface : public TMRVertex {
 public:
  TMRVertexFromSurface( TMRSurface *_surface, double _u, double _v );
  TMRVertexFromSurface( TMRSurface *_surface, TMRPoint p );
  ~TMRVertexFromSurface();
  int evalPoint( TMRPoint *p );

 private:
  double u, v;
  TMRSurface *surface;
};

/*
  Project a curve onto a surface and evaluate the surface location
*/
class TMRCurveFromSurfaceProjection : public TMRCurve {
 public:
  TMRCurveFromSurfaceProjection( TMRSurface *_surface, 
                                 TMRCurve *_curve );
  ~TMRCurveFromSurfaceProjection();  
  void getRange( double *tmin, double *tmax );
  int evalPoint( double t, TMRPoint *p );

 private:
  TMRCurve *curve;
  TMRSurface *surface;
};

/*
  Split the curve
*/
class TMRSplitCurve : public TMRCurve {
 public:
  TMRSplitCurve( TMRCurve *_curve, double _t1, double _t2 );
  TMRSplitCurve( TMRCurve *_curve, TMRPoint *p1, TMRPoint *p2 );
  TMRSplitCurve( TMRCurve *_curve, TMRVertex *p1, TMRVertex *p2 );
  TMRSplitCurve();

  // Get the parameter range
  void getRange( double *tmin, double *tmax );

  // Evaluate the point
  int evalPoint( double t, TMRPoint *X );

 private:
  // Evaluate the bspline curve
  double t1, t2;
  TMRCurve *curve;
};

#endif // TMR_GEOMETRY_H
