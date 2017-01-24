#ifndef TMR_BSPLINE_H
#define TMR_BSPLINE_H

#include "TMRGeometry.h"

/*
  Defines a general transformation: projection/translation/rotation
  for a B-spline object.
*/
class TMRBsplineTransform : public TMREntity {
 public:
  TMRBsplineTransform();
  ~TMRBsplineTransform();

  // Apply the transformation
  void applyTransform( TMRPoint *in, double *win,
                       TMRPoint *out, double *wout,
                       int npts );

 private:
  // Transformation coefficients
  double C[16]; 
};

/*
  This file contains the TMRBsplineCurve and TMRBsplineSurface
  classes that define B-spline and NURBS curves/surfaces for 
  simple geometries. 

  These classes implement the TMRCurve and TMRSurface interfaces,
  respectively.  
*/
class TMRBsplineCurve : public TMRCurve {
 public:
  TMRBsplineCurve( int n, int k, TMRPoint *pts );
  TMRBsplineCurve( int n, int k, const double *Tu, TMRPoint *pts );
  TMRBsplineCurve( int n, int k, const double *Tu, const double *wts, 
                   TMRPoint *pts );
  ~TMRBsplineCurve();

  // Get the parameter range for this edge
  void getRange( double *tmin, double *tmax );
  
  // Given the parametric point, evaluate the x,y,z location
  int evalPoint( double t, TMRPoint *X );

  // Given the x,y,z location, find the parametric coordinates
  int invEvalPoint( TMRPoint X, double *t );

  // Given the parametric point, evaluate the derivative 
  int evalDeriv( double t, TMRPoint *Xt );

  // Refine the knot vector using knot insertion
  TMRBsplineCurve* refineKnots( const double *Tnew, int nnew );

  // Retrieve the underlying data
  void getData( int *_n, int *_k, const double **_Tu,
                const double **_wts, const TMRPoint **_pts ){
    if (_n){ *_n = nctl; }
    if (_k){ *_k = ku; }
    if (_Tu){ *_Tu = Tu; }
    if (_wts){ *_wts = wts; }
    if (_pts){ *_pts = pts; }
  }
  
 private:
  // The number of control points and b-spline order
  int nctl, ku;

  // The knot locations
  double *Tu;

  // The control point locations
  TMRPoint *pts;

  // The weighs (when it is a NURBS curve)
  double *wts;

  // Maximum number of newton iterations for the B-spline inverse
  // point code
  static int max_newton_iters;
};

/*
  The tensor-product B-spline/NURBS surface class
*/
class TMRBsplineSurface : public TMRSurface {
 public:
  TMRBsplineSurface( int nu, int nv, int ku, int kv, 
                     TMRPoint *pts );
  TMRBsplineSurface( int nu, int nv, int ku, int kv, 
                     const double *Tu, const double *Tv, 
                     TMRPoint *pts );
  TMRBsplineSurface( int nu, int nv, int ku, int kv, 
                     const double *Tu, const double *Tv,
                     const double *wts, TMRPoint *pts );
  ~TMRBsplineSurface();

  // Get the parameter range for this surface
  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax );
 
  // Given the parametric point, compute the x,y,z location
  int evalPoint( double u, double v, TMRPoint *X );
  
  // Perform the inverse evaluation
  int invEvalPoint( TMRPoint p, double *u, double *v );

  // Given the parametric point, evaluate the first derivative 
  int evalDeriv( double u, double v, 
                 TMRPoint *Xu, TMRPoint *Xv );
 private:
  // The number of control points and b-spline order
  int ku, kv;
  int nu, nv;

  // The knot locations
  double *Tu, *Tv;

  // The control point locations
  TMRPoint *pts;

  // The weighs (when it is a NURBS curve)
  double *wts;

  // Maximum number of newton iterations for the B-spline inverse
  // point code
  static int max_newton_iters;
};

/*
  The TMRCurveInterpolation enables the creation of a b-spline curve
  using a number of different interpolation strategies.

  The object is created, then internal parameters/options can be set,
  and finally, the interpolation b-spline object can be created.  
*/
class TMRCurveInterpolation : public TMREntity {
 public:
  TMRCurveInterpolation( const TMRPoint *interp, int ninterp );
  ~TMRCurveInterpolation();
  
  // Set the number of control points (must be < ninterp)
  void setNumControlPoints( int _nctl ){
    nctl = _nctl;
  }

  // Create the interpolation
  TMRBsplineCurve *createCurve( int ku );

 private:
  // Get the interpolation points
  void getInterpLoc( double *ubar );

  int nctl;
  int ninterp;
  TMRPoint *interp;
};

/*
  Create a b-spline surface by lofting a spline over b-spline curves
 
  The loft requires 3 steps:
  1. Converting the knot space to span the same to span the same 
  2. Inserting knots to create a common knot space
  3. Creating an interpolating spline to pass through all curves
*/
class TMRCurveLofter : public TMREntity {
 public:
  TMRCurveLofter( TMRBsplineCurve **_curves, int _num_curves );
  ~TMRCurveLofter();

  // Create the surface interpolation
  TMRBsplineSurface* createSurface( int kv );

 private:
  // The curves used for surface lofting
  int num_curves;
  TMRBsplineCurve **curves;
  TMRBsplineCurve **consist;
};

#endif // TMR_BSPLINE_H
