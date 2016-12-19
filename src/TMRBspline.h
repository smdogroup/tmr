#ifndef TMR_BSPLINE_H
#define TMR_BSPLINE_H

#include "TMRGeometry.h"

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

#endif // TMR_BSPLINE_H
