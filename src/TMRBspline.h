/*
  This file is part of the package TMR for adaptive mesh refinement.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

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
  TMRBsplineCurve( int n, int k, const TMRPoint *pts );
  TMRBsplineCurve( int n, int k, const double *Tu, const TMRPoint *pts );
  TMRBsplineCurve( int n, int k, const double *Tu, const double *wts, 
                   const TMRPoint *pts );
  ~TMRBsplineCurve();

  // Get the parameter range for this edge
  void getRange( double *tmin, double *tmax );
  
  // Given the parametric point, evaluate the x,y,z location
  int evalPoint( double t, TMRPoint *X );

  // Given the x,y,z location, find the parametric coordinates
  int invEvalPoint( TMRPoint X, double *t );

  // Given the parametric point, evaluate the derivative 
  int evalDeriv( double t, TMRPoint *Xt );

  // Given the parametric point, evaluate the second derivative 
  int eval2ndDeriv( double t, TMRPoint *Xt );

  // Split the curve at the given parametric point
  void split( double tsplit, TMRBsplineCurve **c1, TMRBsplineCurve **c2 );

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
                     const TMRPoint *pts );
  TMRBsplineSurface( int nu, int nv, int ku, int kv, 
                     const double *Tu, const double *Tv, 
                     const TMRPoint *pts );
  TMRBsplineSurface( int nu, int nv, int ku, int kv, 
                     const double *Tu, const double *Tv,
                     const double *wts, const TMRPoint *pts );
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

  // Given the parametric point, evaluate the second derivatives
  int eval2ndDeriv( double u, double v,
                    TMRPoint *Xuu, TMRPoint *Xuv, TMRPoint *Xvv );

  void getData( int *_nu, int *_nv, int *_ku, int *_kv,
                const double **_Tu, const double **_Tv,
                const double **_wts, const TMRPoint **_pts ){
    if (_nu){ *_nu = nu; }
    if (_nv){ *_nv = nv; }
    if (_ku){ *_ku = ku; }
    if (_kv){ *_kv = kv; }
    if (_Tu){ *_Tu = Tu; }
    if (_Tv){ *_Tv = Tv; }
    if (_wts){ *_wts = wts; }
    if (_pts){ *_pts = pts; }
  }
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
  TMRBsplinePcurve: Parametric curve

  This class defines a parametric curve such that u(t) and v(t) are 
  functions of the parameter t.
*/
class TMRBsplinePcurve : public TMRPcurve {
 public:
  TMRBsplinePcurve( int n, int k, const double *pts );
  TMRBsplinePcurve( int n, int k, const double *Tu, const double *pts );
  TMRBsplinePcurve( int n, int k, const double *Tu, const double *wts, 
                    const double *pts );
  ~TMRBsplinePcurve();

  // Get the parameter range for this edge
  void getRange( double *tmin, double *tmax );
  
  // Given the parametric point, evaluate the x,y,z location
  int evalPoint( double t, double *u, double *v );

  // Given the parametric point, evaluate the derivative 
  int evalDeriv( double t, double *ut, double *vt );

  // Evaluate the second derivative
  int eval2ndDeriv( double t, double *utt, double *vtt );

  // Refine the knot vector using knot insertion
  TMRBsplinePcurve* refineKnots( const double *Tnew, int nnew );
  
 private:
  // The number of control points and b-spline order
  int nctl, ku;

  // The knot locations
  double *Tu;

  // The control point locations
  double *pts;

  // The weighs (when it is a NURBS curve)
  double *wts;
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
