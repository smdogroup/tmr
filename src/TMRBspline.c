#include "TMRBspline.h"
#include "tmrlapack.h"
#include <math.h>
#include <stdio.h>

// The maximum order for the spline
static const int MAX_BSPLINE_ORDER = 6;

/*
  Helper function for comparing integers for qsort/bsearch
*/
static int compare_double( const void *a, const void *b ){
  const double *aa = static_cast<const double*>(a);
  const double *bb = static_cast<const double*>(b);
  if (*aa > *bb){
    return 1;
  }
  else if (*aa < *bb){
    return -1;
  }
  else {
    return 0;
  }
}

/*
  Set up a knot interval vector with uniform spacing

  input:
  n:     number of control points
  k:     order of the b-spline
  umin:  minimum knot value
  umax:  maximum knot value
  
  output:
  T:     knot vector with uniform spacing
*/
void bspline_knots( int n, int k, double umin, double umax, 
                    double *T ){
  // Set the knot vectors in the u-direction
  for ( int i = 0; i < k; i++ ){
    T[i] = umin;
    T[n+k-1-i] = umax;
  }

  // Set the uniform spacing
  for ( int i = k; i < n; i++ ){
    T[i] = umin + ((umax - umin)*(i-k+1))/(n-k+1);
  }
}

/*
  Find the interval for the computing the basis function

  input:
  u: the parametric location
  T: the knot locations
  n: the number of control points
  k: the order of the b-spline

  returns:
  the knot interval
*/
int bspline_interval( double u, const double *T, int n, int k ){
  if (u >= T[n]){
    return n-1;
  }
  else if (u < T[k-1]){
    return k-1;
  }
  else {
    // Use a binary search to locate the interval required
    int low = k-1;
    int high = n;
    int mid = low + (high - low)/2;

    while (u < T[mid] || u >= T[mid+1]){
      if (u < T[mid]){
        high = mid;
      }
      else {
        low = mid;
      }

      mid = low + (high - low)/2;
    }

    return mid;
  }
}

/*
  Evaluate the B-spline basis functions

  output:
  N:  the basis functions - an array of size k

  input:
  idx:  the index of the interval for u
  u:    the parametric position
  Tu:   the knot vector
  ku:   the order of the basis functions
  work: a temporary array of size 2*k

  u is on the idx-th knot span such that u is in the interval
  u \in [Tu[idx], Tu[idx+1]) 
*/
void bspline_basis( double *N, const int idx, const double u, 
                    const double *Tu, 
                    const int ku, double *work ){
  N[0] = 1.0;
  
  // Set the pointers for the temporary work arrays
  // Note that left[j] = u - Tu[i+1 - j]
  // and right[j] = Tu[i+j] - u
  double *left = &work[0];
  double *right = &work[ku];
  
  for ( int j = 1; j < ku; j++ ){
    left[j] = u - Tu[idx+1-j];
    right[j] = Tu[idx+j] - u;
    
    N[j] = 0.0;
    for ( int i = 0; i < j; i++ ){
      double temp = N[i]/(right[i+1] + left[j-i]);
      N[i] = N[j] + right[i+1]*temp;
      N[j] = left[j-i]*temp;
    }
  }
}

/*
  Compute the derivative of the b-spline basis

  output:
  N:  the basis functions - an array of size k

  input:
  idx:  the index of the interval for u
  u:    the parametric position
  Tu:   the knot vector
  ku:   the order of the basis functions
  work: a temporary array of size 2*k + k*k

  u is on the idx-th knot span such that u is in the interval
  u \in [Tu[idx], Tu[idx+1]) 
*/
void bspline_basis_derivative( double *N, const int idx, const double u,
                               int ideriv, const double *Tu,
                               const int ku, double *work ){

  if (ideriv >= ku){
    // Set all higher-order derivatives to zero
    memset(&N[ku*ku], 0, (ideriv+1-ku)*ku*sizeof(double));

    // set the maximum derivative to ideriv
    ideriv = ku-1;
  }

  // Set the pointers for the temporary work arrays
  // Note that left[j] = u - Tu[i+1 - j]
  // and right[j] = Tu[i+j] - u
  double *left = &work[0];
  double *right = &work[ku];

  // ndu is a matrix of total size ku*ku
  // The basis functions are stored in the upper triangular part of
  // the matrix such that N_{idx-j, i} = ndu[i + j*ku] if i >= j is
  // the basis function. The knot differences are stored in the lower
  // portion of the matrix such that
  // ndu[i + j*ku] = u_{idx+i+1} - u_{idx+j-i}
  double *ndu = &work[2*ku];
  ndu[0] = 1.0;

  // Compute all orders of the basis functions
  for ( int j = 1; j < ku; j++ ){
    left[j] = u - Tu[idx+1-j];
    right[j] = Tu[idx+j] - u;

    double njj = 0.0;
    for ( int i = 0; i < j; i++ ){
      // Compute the difference Tu[idx+i+1] - Tu[idx+j-i]
      ndu[i + j*ku] = right[i+1] + left[j-i];
      double temp = ndu[j-1 + i*ku]/ndu[i + j*ku];

      // Compute the basis function
      ndu[j + i*ku] = njj + right[i+1]*temp;
      njj = left[j-i]*temp;
    }

    // Store the basis function
    ndu[j*(ku+1)] = njj; 
  }

  // Set the basis functions
  for ( int i = 0; i < ku; i++ ){
    N[i] = ndu[(ku-1) + i*ku]; 
  }

  // Set the temporary arrays for the a-coefficients
  double *a0 = &work[0];
  double *a1 = &work[ku];

  for ( int i = 0; i < ku; i++ ){
    a0[0] = 1.0;

    for ( int k = 1; k <= ideriv; k++ ){ 
      double d = 0.0;

      // Compute the first of the a-terms
      // a_{k,0} = a_{k-1,0}/(u_{i+ku-k} - u_{i})
      if (i >= k){
        a1[0] = a0[0]/ndu[(i-k) + (ku-k)*ku];
        d = a1[0]*ndu[(ku-k-1) + (i-k)*ku];
      }

      int jstart = k-i;
      if (i >= k-1){ jstart = 1; }

      int jend = ku-i-1;
      if (i <= ku-k){ jend = k-1; }

      for ( int j = jstart; j <= jend; j++ ){
        // Compute a_{k,j} = (a_{k-1,j} - a_{k-1,j-1})
        // /(u_{i+ku+j-k} -u_{i+j})
        a1[j] = (a0[j] - a0[j-1])/ndu[(i-k+j) + (ku-k)*ku];
        d += a1[j]*ndu[(ku-k-1) + (i-k+j)*ku];
      }

      // Compute the term 
      // a_{k,k} = -a_{k-1}/(u_{i+ku} - u_{i+k})
      if (i <= ku-k-1){
        a1[k] = -a0[k-1]/ndu[i + (ku-k)*ku];
        d += a1[k]*ndu[(ku-k-1) + i*ku];
      }

      // Set the basis function
      N[i + k*ku] = d;

      // Swap the rows of a
      double *t = a0;
      a0 = a1;  a1 = t;
    }
  }

  // Multiply the basis by the factorial term
  int r = ku-1;
  for ( int k = 1; k <= ideriv; k++ ){
    for ( int j = 0; j < ku; j++ ){
      N[j + k*ku] *= r;
    }
    r *= (ku-1 - k);
  }
}

/*
  Create a B-spline curve with a uniform knot vector
*/
TMRBsplineCurve::TMRBsplineCurve( int _n, int _k,
                                  TMRPoint *_pts ){
  nctl = _n;

  // Check the bounds on the bspline order
  ku = _k;
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }

  // Set a uniform spacing of the knot vector
  Tu = new double[ nctl+ku ];
  bspline_knots(nctl, ku, 0.0, 1.0, Tu);

  // Allocate and copy over the points
  pts = new TMRPoint[ nctl ];
  memcpy(pts, _pts, nctl*sizeof(TMRPoint));
  
  // Uniform weights
  wts = NULL;
}

/*
  Create a B-spline curve with a specified knot vector

  n:    the number of control points
  k:    the order of the B-spline
  Tu:   the knot vector (n+k)
  pts:  the TMRPoint interpolation points
*/
TMRBsplineCurve::TMRBsplineCurve( int _n, int _k, 
                                  const double *_Tu,
                                  TMRPoint *_pts ){
  nctl = _n;

  // Check the bounds on the bspline order
  ku = _k;
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }

  // Allocate the knots and control points
  Tu = new double[ nctl+ku ];
  if (_Tu){
    memcpy(Tu, _Tu, (nctl+ku)*sizeof(double));
  }
  else {
    bspline_knots(nctl, ku, 0.0, 1.0, Tu);
  }
  
  pts = new TMRPoint[ nctl ];
  memcpy(pts, _pts, nctl*sizeof(TMRPoint));
  
  // Uniform weights
  wts = NULL;
}

/*
  Create a B-spline curve with a specified knot vector

  n:    the number of control points
  k:    the order of the B-spline
  Tu:   the knot vector (n+k)
  pts:  the TMRPoint interpolation points
*/
TMRBsplineCurve::TMRBsplineCurve( int _n, int _k, 
                                  const double *_Tu,
                                  const double *_wts, 
                                  TMRPoint *_pts ){
  nctl = _n;

  // Check the bounds on the bspline order
  ku = _k;
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }

  // Set the knots
  Tu = new double[ nctl+ku ];
  if (_Tu){
    memcpy(Tu, _Tu, (nctl+ku)*sizeof(double));
  }
  else {
    bspline_knots(nctl, ku, 0.0, 1.0, Tu);
  }
  
  // Allocate the weights
  wts = new double[ nctl ];
  memcpy(wts, _wts, nctl*sizeof(double));

  // Copy over the points
  pts = new TMRPoint[ nctl ];
  memcpy(pts, _pts, nctl*sizeof(TMRPoint));
}

/*
  Free the curve
*/
TMRBsplineCurve::~TMRBsplineCurve(){
  delete [] Tu;
  delete [] pts;
  if (wts){ delete [] wts; }
}

// Set the maximum number of newton iterations
int TMRBsplineCurve::max_newton_iters = 25;

/*
  Get the range of the knot vectors
*/
void TMRBsplineCurve::getRange( double *tmin, double *tmax ){
  *tmin = Tu[0];
  *tmax = Tu[nctl+ku-1];
}
  
/*
  Given the parametric point, evaluate the x,y,z location
*/  
int TMRBsplineCurve::evalPoint( double t, TMRPoint *X ){
  double Nu[MAX_BSPLINE_ORDER];
  double work[2*MAX_BSPLINE_ORDER];

  // Compute the knot span
  int intu = bspline_interval(t, Tu, nctl, ku);

  // Compute the basis functions
  bspline_basis(Nu, intu, t, Tu, ku, work);

  // Set the interval to the initial control point
  intu = intu - ku + 1;

  // Zero the point location
  X->zero();

  // If this is a NURBS curve add the effect of the weights
  if (wts){
    // Evaluate the b-spline
    double w = 0.0;
    for ( int i = 0; i < ku; i++ ){
      X->x += wts[intu + i]*Nu[i]*pts[intu + i].x;
      X->y += wts[intu + i]*Nu[i]*pts[intu + i].y;
      X->z += wts[intu + i]*Nu[i]*pts[intu + i].z;
      w += Nu[i]*wts[intu + i];
    }

    // Divide through by the weights
    if (w != 0.0){
      w = 1.0/w;
      X->x *= w;
      X->y *= w;
      X->z *= w;
    }
  }
  else {
    // Evaluate the b-spline
    for ( int i = 0; i < ku; i++ ){
      X->x += Nu[i]*pts[intu + i].x;
      X->y += Nu[i]*pts[intu + i].y;
      X->z += Nu[i]*pts[intu + i].z;
    }
  }

  // Success - we didn't fail
  return 0;
}

/*
  Perform the inverse point evaluation
*/
int TMRBsplineCurve::invEvalPoint( TMRPoint point, double *tf ){
  // Did this evaluation fail or not? 
  int fail = 1;

  // Set a default value for the output parameter value
  *tf = 0.0;

  // Temporary internal data
  double Nu[3*MAX_BSPLINE_ORDER];
  double work[2*MAX_BSPLINE_ORDER + MAX_BSPLINE_ORDER*MAX_BSPLINE_ORDER];
  double t = 0.5*(Tu[0] + Tu[nctl+ku-1]);

  // Perform a newton iteration until convergence
  for ( int j = 0; j < max_newton_iters; j++ ){
    // Compute the knot span
    int intu = bspline_interval(t, Tu, nctl, ku);
    
    // Compute the basis functions
    int idu = 2;
    bspline_basis_derivative(Nu, intu, t, idu, Tu, ku, work);

    // Set the interval to the initial control point
    intu = intu - ku + 1;

    // The positions and their derivatives on the b-spline
    // surface
    TMRPoint X, Xt, Xtt;
    X.zero();  Xt.zero();  Xtt.zero();

    if (wts){
      // Evaluate the weighted b-spline
      TMRPoint p, pt, ptt;
      p.zero();  pt.zero();  ptt.zero();

      double w = 0.0, wt = 0.0, wtt = 0.0;
      for ( int i = 0; i < ku; i++ ){
        // Evaluate the numerator and its derivative
        p.x += wts[intu + i]*Nu[i]*pts[intu + i].x;
        p.y += wts[intu + i]*Nu[i]*pts[intu + i].y;
        p.z += wts[intu + i]*Nu[i]*pts[intu + i].z;
        
        pt.x += wts[intu + i]*Nu[ku + i]*pts[intu + i].x;
        pt.y += wts[intu + i]*Nu[ku + i]*pts[intu + i].y;
        pt.z += wts[intu + i]*Nu[ku + i]*pts[intu + i].z;

        ptt.x += wts[intu + i]*Nu[2*ku + i]*pts[intu + i].x;
        ptt.y += wts[intu + i]*Nu[2*ku + i]*pts[intu + i].y;
        ptt.z += wts[intu + i]*Nu[2*ku + i]*pts[intu + i].z;
      
        // Compute the weights and their derivatives
        w += Nu[i]*wts[intu + i];
        wt += Nu[ku + i]*wts[intu + i];
        wtt += Nu[2*ku + i]*wts[intu + i];
      }

      // Divide through by the weights
      if (w != 0.0){
        // Compute the first and second derivatives of the inverse of
        // the weights
        w = 1.0/w;
        wtt = 2.0*w*w*w*wt*wt - w*w*wtt;
        wt = -w*w*wt;

        // Evaluate the derivatives 
        X.x = w*p.x;
        X.y = w*p.y;
        X.z = w*p.z;
        
        Xt.x = w*pt.x + wt*p.x;
        Xt.y = w*pt.y + wt*p.y;
        Xt.z = w*pt.z + wt*p.z;
        
        Xtt.x = 2.0*wt*pt.x + w*ptt.x + wtt*p.x;
        Xtt.y = 2.0*wt*pt.y + w*ptt.y + wtt*p.y;
        Xtt.z = 2.0*wt*pt.z + w*ptt.z + wtt*p.z;
      }
      else {
        fail = 1;
        return fail;
      }
    }
    else {
      // Evaluate the b-spline
      for ( int i = 0; i < ku; i++ ){
        X.x += Nu[i]*pts[intu + i].x;
        X.y += Nu[i]*pts[intu + i].y;
        X.z += Nu[i]*pts[intu + i].z;
        
        Xt.x += Nu[ku + i]*pts[intu + i].x;
        Xt.y += Nu[ku + i]*pts[intu + i].y;
        Xt.z += Nu[ku + i]*pts[intu + i].z;
        
        Xtt.x += Nu[2*ku + i]*pts[intu + i].x;
        Xtt.y += Nu[2*ku + i]*pts[intu + i].y;
        Xtt.z += Nu[2*ku + i]*pts[intu + i].z;
      }
    }

    // Compute the dot product of the tangent w.r.t. the difference
    // between the point on the curve and the projection point
    TMRPoint r;
    r.x = (X.x - point.x);
    r.y = (X.y - point.y);
    r.z = (X.z - point.z);

    // Compute the residual
    double res = Xt.dot(r);

    // Compute the derivative of the dot product
    double deriv = Xtt.dot(r) + Xt.dot(Xt);

    // Compute the update for t
    double tnew = t;
    if (deriv != 0.0){
      tnew = t - res/deriv;
    }
    
    // Truncate the new value to the boundary
    // if (open curve){
    if (tnew < Tu[0]){
      tnew = Tu[0];
    }
    else if (tnew > Tu[nctl+ku-1]){
      tnew = Tu[nctl+ku-1];
    }
    // }
    // else { }
  
    // Check if the convergence test satisfied
    if (fabs(r.x) < eps_dist && 
        fabs(r.y) < eps_dist && 
        fabs(r.z) < eps_dist){
      *tf = t;
      return 0;
    }
    if (fabs((tnew - t)*Xt.x) < eps_dist &&
        fabs((tnew - t)*Xt.y) < eps_dist &&
        fabs((tnew - t)*Xt.z) < eps_dist){
      *tf = t;
      return 0;
    }

    // Perform the cosine check
    double dot1 = Xt.dot(Xt);
    double dot2 = r.dot(r);
    if (res*res < eps_cosine*eps_cosine*dot1*dot2){
      *tf = t;
      return 0;
    }

    // Update the new parameter value
    t = tnew;
  }

  // The newton method has failed!!
  fail = 1;
  return fail;
}

/*
  Given the parametric point, evaluate the derivative 
*/  
int TMRBsplineCurve::evalDeriv( double t, TMRPoint *Xt ){
  double Nu[2*MAX_BSPLINE_ORDER];
  double work[2*MAX_BSPLINE_ORDER + MAX_BSPLINE_ORDER*MAX_BSPLINE_ORDER];

  // Compute the knot span
  int intu = bspline_interval(t, Tu, nctl, ku);

  // Compute the basis functions
  int idu = 1;
  bspline_basis_derivative(Nu, intu, t, idu, Tu, ku, work);

  // Set the interval to the initial control point
  intu = intu - ku + 1;

  // Evaluate the b-spline
  Xt->x = 0.0;
  Xt->y = 0.0;
  Xt->z = 0.0;

  if (wts){
    // Evaluate the weighted b-spline
    TMRPoint p, pt;
    p.x = p.y = p.z = 0.0;
    pt.x = pt.y = pt.z = 0.0;

    double w = 0.0, wt = 0.0;
    for ( int i = 0; i < ku; i++ ){
      // Evaluate the numerator and its derivative
      p.x += wts[intu + i]*Nu[i]*pts[intu + i].x;
      p.y += wts[intu + i]*Nu[i]*pts[intu + i].y;
      p.z += wts[intu + i]*Nu[i]*pts[intu + i].z;

      pt.x += wts[intu + i]*Nu[ku + i]*pts[intu + i].x;
      pt.y += wts[intu + i]*Nu[ku + i]*pts[intu + i].y;
      pt.z += wts[intu + i]*Nu[ku + i]*pts[intu + i].z;

      // Compute the weights
      w += Nu[i]*wts[intu + i];
      wt += Nu[i + ku]*wts[intu + i];
    }

    // Divide through by the weights
    if (w != 0.0){
      w = 1.0/w;
      wt = -w*w*wt;
      Xt->x = w*pt.x + wt*p.x;
      Xt->y = w*pt.y + wt*p.y;
      Xt->z = w*pt.z + wt*p.z;
    }
  }
  else {
    // Evaluate the b-spline
    for ( int i = 0; i < ku; i++ ){
      Xt->x += Nu[ku + i]*pts[intu + i].x;
      Xt->y += Nu[ku + i]*pts[intu + i].y;
      Xt->z += Nu[ku + i]*pts[intu + i].z;
    }
  }

  return 0;
}

/*
  Refine the knot vector by inserting the knots specified in the input
  Tnew. This function returns a new, equivalent curve with the new
  knots at the specified locations.  
*/
TMRBsplineCurve* TMRBsplineCurve::refineKnots( const double *_Tnew, 
                                               int nnew ){
  // Allocate the new vectors 
  int kbar = ku;
  int nbar = nctl + nnew;
  double *Tbar = new double[ nbar+ku ];
  TMRPoint *pbar = new TMRPoint[ nbar ];

  // Allocate the new Tnew and
  double *Tnew = new double[ nnew ];
  memcpy(Tnew, _Tnew, nnew*sizeof(double));
  qsort(Tnew, nnew, sizeof(double), compare_double);

  // Allocate the control point weights
  double *wbar = NULL;
  if (wts){ wbar = new double[ nbar ]; }

  // Locate the first/end interval
  int a = bspline_interval(Tnew[0], Tu, nctl, ku);
  int b = bspline_interval(Tnew[nnew-1], Tu, nctl, ku) + 1;

  // Copy over the points that will not change
  for ( int j = 0; j <= a-ku+1; j++ ){
    pbar[j] = pts[j];
  }
  for ( int j = b-1; j < nctl; j++ ){
    pbar[j+nnew] = pts[j];
  }

  if (wts){
    // Copy over the NURBS weights which will not change
    for ( int j = 0; j <= a-ku+1; j++ ){
      wbar[j] = wts[j];
    }
    for ( int j = b-1; j < nctl; j++ ){
      wbar[j+nnew] = wts[j];
    }
  }
  
  // Copy over the knots that do not change
  for ( int j = 0; j <= a; j++ ){
    Tbar[j] = Tu[j];
  }
  for ( int j = b+ku-1; j < nctl+ku; j++ ){
    Tbar[j+nnew] = Tu[j];
  }

  // Add the new values
  for ( int i = b+ku-2, j = nnew-1, k = b+ku+nnew-2; j >= 0; j-- ){
    while (Tnew[j] <= Tu[i] && i > a){
      // Read over points and weights
      pbar[k-ku] = pts[i-ku];
      if (wts){
        wbar[k-ku] = wts[i-ku];
      }
      // Read in the new knots
      Tbar[k] = Tu[i];
      k--; 
      i--;
    }
    // Read in the last pt/weight pair
    pbar[k-ku] = pbar[k-ku+1];
    if (wts){
      wbar[k-ku] = wbar[k-ku+1];
    }

    // Form the linear combinations of the point/weight pairs
    for ( int l = 1; l < ku; l++ ){
      int ind = k-(ku-1)+l;     
      double alpha = Tbar[k+l] - Tnew[j];
     
      if (alpha == 0.0){
        pbar[ind-1] = pbar[ind];
      }
      else {
        alpha = alpha/(Tbar[k+l] - Tu[i-ku+1+l]);       
        if (wts){
          pbar[ind-1].x = (alpha*wbar[ind-1]*pbar[ind-1].x + 
                           (1.0 - alpha)*wbar[ind]*pbar[ind].x);
          pbar[ind-1].y = (alpha*wbar[ind-1]*pbar[ind-1].y + 
                           (1.0 - alpha)*wbar[ind]*pbar[ind].y);
          pbar[ind-1].z = (alpha*wbar[ind-1]*pbar[ind-1].z + 
                           (1.0 - alpha)*wbar[ind]*pbar[ind].z);
          wbar[ind-1] = (alpha*wbar[ind-1] + (1.0 - alpha)*wbar[ind]);
        }
        else {
          pbar[ind-1].x = alpha*pbar[ind-1].x + (1.0 - alpha)*pbar[ind].x;
          pbar[ind-1].y = alpha*pbar[ind-1].y + (1.0 - alpha)*pbar[ind].y;
          pbar[ind-1].z = alpha*pbar[ind-1].z + (1.0 - alpha)*pbar[ind].z;
        }
      }
    }

    // Copy over the knot
    Tbar[k] = Tnew[j];
    k--;
  }

  // Now we have the new inputs
  TMRBsplineCurve *refined = NULL;

  if (wbar){
    refined = new TMRBsplineCurve(nbar, ku, Tbar, wbar, pbar);
  }
  else {
    refined = new TMRBsplineCurve(nbar, ku, Tbar, pbar);
  }

  // Free the allocated data
  delete [] Tnew;
  delete [] Tbar;
  delete [] pbar;
  if (wbar){
    delete [] wbar;
  }

  // return the new object
  return refined;
}

/*
  Create a B-spline curve with a uniform knot vector
*/
TMRBsplinePcurve::TMRBsplinePcurve( int _n, int _k, 
                                    const double *_pts ){
  nctl = _n;

  // Check the bounds on the bspline order
  ku = _k;
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }

  // Set a uniform spacing of the knot vector
  Tu = new double[ nctl+ku ];
  bspline_knots(nctl, ku, 0.0, 1.0, Tu);

  // Allocate and copy over the points
  pts = new double[ 2*nctl ];
  memcpy(pts, _pts, 2*nctl*sizeof(double));
  
  // Uniform weights
  wts = NULL;
}

/*
  Create a B-spline curve with a specified knot vector

  n:    the number of control points
  k:    the order of the B-spline
  Tu:   the knot vector (n+k)
  pts:  the TMRPoint interpolation points
*/
TMRBsplinePcurve::TMRBsplinePcurve( int _n, int _k, 
                                    const double *_Tu,
                                    const double *_pts ){
  nctl = _n;

  // Check the bounds on the bspline order
  ku = _k;
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }

  // Allocate the knots and control points
  Tu = new double[ nctl+ku ];
  if (_Tu){
    memcpy(Tu, _Tu, (nctl+ku)*sizeof(double));
  }
  else {
    bspline_knots(nctl, ku, 0.0, 1.0, Tu);
  }
  
  pts = new double[ 2*nctl ];
  memcpy(pts, _pts, 2*nctl*sizeof(double));
  
  // Uniform weights
  wts = NULL;
}

/*
  Create a B-spline curve with a specified knot vector

  n:    the number of control points
  k:    the order of the B-spline
  Tu:   the knot vector (n+k)
  pts:  the TMRPoint interpolation points
*/
TMRBsplinePcurve::TMRBsplinePcurve( int _n, int _k, 
                                    const double *_Tu,
                                    const double *_wts, 
                                    const double *_pts ){
  nctl = _n;

  // Check the bounds on the bspline order
  ku = _k;
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }

  // Set the knots
  Tu = new double[ nctl+ku ];
  if (_Tu){
    memcpy(Tu, _Tu, (nctl+ku)*sizeof(double));
  }
  else {
    bspline_knots(nctl, ku, 0.0, 1.0, Tu);
  }
  
  // Allocate the weights
  wts = new double[ nctl ];
  memcpy(wts, _wts, nctl*sizeof(double));

  // Copy over the points
  pts = new double[ 2*nctl ];
  memcpy(pts, _pts, 2*nctl*sizeof(double));
}

/*
  Free the curve
*/
TMRBsplinePcurve::~TMRBsplinePcurve(){
  delete [] Tu;
  delete [] pts;
  if (wts){ delete [] wts; }
}

/*
  Get the range of the knot vectors
*/
void TMRBsplinePcurve::getRange( double *tmin, double *tmax ){
  *tmin = Tu[0];
  *tmax = Tu[nctl+ku-1];
}
  
/*
  Given the parametric point, evaluate the x,y,z location
*/  
int TMRBsplinePcurve::evalPoint( double t, double *u, double *v ){
  double Nu[MAX_BSPLINE_ORDER];
  double work[2*MAX_BSPLINE_ORDER];

  // Compute the knot span
  int intu = bspline_interval(t, Tu, nctl, ku);

  // Compute the basis functions
  bspline_basis(Nu, intu, t, Tu, ku, work);

  // Set the interval to the initial control point
  intu = intu - ku + 1;

  // Zero the point location
  *u = *v = 0.0;

  // If this is a NURBS curve add the effect of the weights
  if (wts){
    // Evaluate the b-spline
    double w = 0.0;
    for ( int i = 0; i < ku; i++ ){
      *u += wts[intu + i]*Nu[i]*pts[2*(intu + i)];
      *v += wts[intu + i]*Nu[i]*pts[2*(intu + i)+1];
      w += Nu[i]*wts[intu + i];
    }

    // Divide through by the weights
    if (w != 0.0){
      w = 1.0/w;
      *u *= w;
      *v *= w;
    }
  }
  else {
    // Evaluate the b-spline
    for ( int i = 0; i < ku; i++ ){
      *u += Nu[i]*pts[2*(intu + i)];
      *v += Nu[i]*pts[2*(intu + i)+1];
    }
  }

  // Success - we didn't fail
  return 0;
}

/*
  Given the parametric point, evaluate the derivative 
*/  
int TMRBsplinePcurve::evalDeriv( double t, double *ut, double *vt ){
  double Nu[2*MAX_BSPLINE_ORDER];
  double work[2*MAX_BSPLINE_ORDER + MAX_BSPLINE_ORDER*MAX_BSPLINE_ORDER];

  // Compute the knot span
  int intu = bspline_interval(t, Tu, nctl, ku);

  // Compute the basis functions
  int idu = 1;
  bspline_basis_derivative(Nu, intu, t, idu, Tu, ku, work);

  // Set the interval to the initial control point
  intu = intu - ku + 1;

  // Evaluate the b-spline
  *ut = *vt = 0.0;

  if (wts){
    // Evaluate the weighted b-spline
    double u = 0.0, v = 0.0;
    double w = 0.0, wt = 0.0;
    for ( int i = 0; i < ku; i++ ){
      // Evaluate the numerator and its derivative
      u += wts[intu + i]*Nu[i]*pts[2*(intu + i)];
      v += wts[intu + i]*Nu[i]*pts[2*(intu + i)+1];
      
      *ut += wts[intu + i]*Nu[ku + i]*pts[2*(intu + i)];
      *vt += wts[intu + i]*Nu[ku + i]*pts[2*(intu + i)+1];
      
      // Compute the weights
      w += Nu[i]*wts[intu + i];
      wt += Nu[i + ku]*wts[intu + i];
    }

    // Divide through by the weights
    if (w != 0.0){
      w = 1.0/w;
      wt = -w*w*wt;
      *ut = w*(*ut) + wt*u;
      *vt = w*(*vt) + wt*v;
    }
  }
  else {
    // Evaluate the b-spline
    for ( int i = 0; i < ku; i++ ){
      *ut += Nu[ku + i]*pts[2*(intu + i)];
      *vt += Nu[ku + i]*pts[2*(intu + i)+1];
    }
  }

  return 0;
}

/*
  Refine the knot vector by inserting the knots specified in the input
  Tnew. This function returns a new, equivalent curve with the new
  knots at the specified locations.  
*/
TMRBsplinePcurve* TMRBsplinePcurve::refineKnots( const double *_Tnew, 
                                                 int nnew ){
  // Allocate the new vectors 
  int kbar = ku;
  int nbar = nctl + nnew;
  double *Tbar = new double[ nbar+ku ];
  double *pbar = new double[ 2*nbar ];

  // Allocate the new Tnew and
  double *Tnew = new double[ nnew ];
  memcpy(Tnew, _Tnew, nnew*sizeof(double));
  qsort(Tnew, nnew, sizeof(double), compare_double);

  // Allocate the control point weights
  double *wbar = NULL;
  if (wts){ wbar = new double[ nbar ]; }

  // Locate the first/end interval
  int a = bspline_interval(Tnew[0], Tu, nctl, ku);
  int b = bspline_interval(Tnew[nnew-1], Tu, nctl, ku) + 1;

  // Copy over the points that will not change
  for ( int j = 0; j <= a-ku+1; j++ ){
    pbar[j] = pts[j];
  }
  for ( int j = b-1; j < nctl; j++ ){
    pbar[j+nnew] = pts[j];
  }

  if (wts){
    // Copy over the NURBS weights which will not change
    for ( int j = 0; j <= a-ku+1; j++ ){
      wbar[j] = wts[j];
    }
    for ( int j = b-1; j < nctl; j++ ){
      wbar[j+nnew] = wts[j];
    }
  }
  
  // Copy over the knots that do not change
  for ( int j = 0; j <= a; j++ ){
    Tbar[j] = Tu[j];
  }
  for ( int j = b+ku-1; j < nctl+ku; j++ ){
    Tbar[j+nnew] = Tu[j];
  }

  // Add the new values
  for ( int i = b+ku-2, j = nnew-1, k = b+ku+nnew-2; j >= 0; j-- ){
    while (Tnew[j] <= Tu[i] && i > a){
      // Read over points and weights
      pbar[k-ku] = pts[i-ku];
      if (wts){
        wbar[k-ku] = wts[i-ku];
      }
      // Read in the new knots
      Tbar[k] = Tu[i];
      k--; 
      i--;
    }
    // Read in the last pt/weight pair
    pbar[k-ku] = pbar[k-ku+1];
    if (wts){
      wbar[k-ku] = wbar[k-ku+1];
    }

    // Form the linear combinations of the point/weight pairs
    for ( int l = 1; l < ku; l++ ){
      int ind = k-(ku-1)+l;     
      double alpha = Tbar[k+l] - Tnew[j];
     
      if (alpha == 0.0){
        pbar[ind-1] = pbar[ind];
      }
      else {
        alpha = alpha/(Tbar[k+l] - Tu[i-ku+1+l]);       
        if (wts){
          pbar[2*(ind-1)] = (alpha*wbar[ind-1]*pbar[2*(ind-1)] + 
                            (1.0 - alpha)*wbar[ind]*pbar[2*ind]);
          pbar[2*(ind-1)+1] = (alpha*wbar[ind-1]*pbar[2*(ind-1)+1] + 
                              (1.0 - alpha)*wbar[ind]*pbar[2*ind]);
          wbar[ind-1] = (alpha*wbar[ind-1] + (1.0 - alpha)*wbar[ind]);
        }
        else {
          pbar[2*(ind-1)] = alpha*pbar[2*(ind-1)] + (1.0 - alpha)*pbar[2*ind];
          pbar[2*(ind-1)+1] = alpha*pbar[2*(ind-1)+1] + (1.0 - alpha)*pbar[2*ind+1];
        }
      }
    }

    // Copy over the knot
    Tbar[k] = Tnew[j];
    k--;
  }

  // Now we have the new inputs
  TMRBsplinePcurve *refined = NULL;

  if (wbar){
    refined = new TMRBsplinePcurve(nbar, ku, Tbar, wbar, pbar);
  }
  else {
    refined = new TMRBsplinePcurve(nbar, ku, Tbar, pbar);
  }

  // Free the allocated data
  delete [] Tnew;
  delete [] Tbar;
  delete [] pbar;
  if (wbar){
    delete [] wbar;
  }

  // return the new object
  return refined;
}

/*
  Evaluate the TMRSurface
*/
TMRBsplineSurface::TMRBsplineSurface( int _nu, int _nv, 
                                      int _ku, int _kv, 
                                      TMRPoint *_pts ){
  // Record the number of control points and the basis order
  nu = _nu;
  nv = _nv;
  ku = _ku;
  kv = _kv;
  
  // Check the bounds on the bspline order
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }
  if (kv < 2){ kv = 2; }
  if (kv > MAX_BSPLINE_ORDER){ kv = MAX_BSPLINE_ORDER; }

  // Allocate the knot vectors
  Tu = new double[ nu+ku ];
  Tv = new double[ nv+kv ];
  bspline_knots(nu, ku, 0.0, 1.0, Tu);
  bspline_knots(nv, kv, 0.0, 1.0, Tv);

  // Set the wts = 1.0
  wts = NULL;

  // Allocate the points array
  pts = new TMRPoint[ nu*nv ];
  memcpy(pts, _pts, nu*nv*sizeof(TMRPoint));
}

/*
  Allocate the B-spline surface and set the knot locations
*/
TMRBsplineSurface::TMRBsplineSurface( int _nu, int _nv, 
                                      int _ku, int _kv, 
                                      const double *_Tu, 
                                      const double *_Tv, 
                                      TMRPoint *_pts ){
  // Record the number of control points and the basis order
  nu = _nu;
  nv = _nv;
  ku = _ku;
  kv = _kv;

  // Check the bounds on the bspline order
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }
  if (kv < 2){ kv = 2; }
  if (kv > MAX_BSPLINE_ORDER){ kv = MAX_BSPLINE_ORDER; }
  
  // Allocate the knot vectors
  Tu = new double[ nu+ku ];
  Tv = new double[ nv+kv ];
  if (_Tu){
    memcpy(Tu, _Tu, (nu+ku)*sizeof(double));
  }
  else {
    bspline_knots(nu, ku, 0.0, 1.0, Tu);
  }
  if (_Tv){
    memcpy(Tv, _Tv, (nv+kv)*sizeof(double));
  }
  else {
    bspline_knots(nv, kv, 0.0, 1.0, Tv);
  }
  
  // Set the wts = 1.0
  wts = NULL;

  // Allocate the points array
  pts = new TMRPoint[ nu*nv ];
  memcpy(pts, _pts, nu*nv*sizeof(TMRPoint));
}

TMRBsplineSurface::TMRBsplineSurface( int _nu, int _nv, 
                                      int _ku, int _kv, 
                                      const double *_Tu, const double *_Tv,
                                      const double *_wts, 
                                      TMRPoint *_pts ){
  // Record the number of control points and the basis order
  nu = _nu;
  nv = _nv;
  ku = _ku;
  kv = _kv;

  // Check the bounds on the bspline order
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }
  if (kv < 2){ kv = 2; }
  if (kv > MAX_BSPLINE_ORDER){ kv = MAX_BSPLINE_ORDER; }
  
  // Allocate the knot vectors
  Tu = new double[ nu+ku ];
  Tv = new double[ nv+kv ];
  if (_Tu){
    memcpy(Tu, _Tu, (nu+ku)*sizeof(double));
  }
  else {
    bspline_knots(nu, ku, 0.0, 1.0, Tu);
  }
  if (_Tv){
    memcpy(Tv, _Tv, (nv+kv)*sizeof(double));
  }
  else {
    bspline_knots(nv, kv, 0.0, 1.0, Tv);
  }

  // Set the weights for the rational part of the NURBS
  wts = new double[ nu*nv ];
  memcpy(wts, _wts, nu*nv*sizeof(double));

  // Allocate the points array
  pts = new TMRPoint[ nu*nv ];
  memcpy(pts, _pts, nu*nv*sizeof(TMRPoint));
}

/*
  Free the B-spline/NURBS surface
*/
TMRBsplineSurface::~TMRBsplineSurface(){
  delete [] Tu;
  delete [] Tv;
  if (wts){ delete [] wts; }
  delete [] pts;
}

/*
  Set the maximum number of newton iterations
*/
int TMRBsplineSurface::max_newton_iters = 25;

/* 
   Get the parameter range for this surface from the knot vectors
*/
void TMRBsplineSurface::getRange( double *umin, double *vmin,
                                  double *umax, double *vmax ){
  *umin = Tu[0];
  *umax = Tu[nu+ku-1];
  *vmin = Tv[0];
  *vmax = Tv[nv+kv-1];
}
 
/*
  Given the parametric point, compute the x,y,z location
*/
int TMRBsplineSurface::evalPoint( double u, double v, TMRPoint *X ){
  // The basis functions/work arrays
  double Nu[MAX_BSPLINE_ORDER], Nv[MAX_BSPLINE_ORDER];
  double work[2*MAX_BSPLINE_ORDER];

  // Compute the knot intervals
  int intu = bspline_interval(u, Tu, nu, ku);
  int intv = bspline_interval(v, Tv, nv, kv);

  // Evaluate the basis functions
  bspline_basis(Nu, intu, u, Tu, ku, work);
  bspline_basis(Nv, intv, v, Tv, kv, work);
  
  // Set the interval to the initial control point
  intu = intu - ku + 1;
  intv = intv - kv + 1;

  // Zero the point location
  X->zero();

  // If this is a NURBS surface add the effect of the weights
  if (wts){
    // Evaluate the b-spline
    double w = 0.0;
    for ( int j = 0; j < kv; j++ ){
      for ( int i = 0; i < ku; i++ ){
        int index = intu+i + (intv+j)*nu;
        X->x += wts[index]*Nu[i]*Nv[j]*pts[index].x;
        X->y += wts[index]*Nu[i]*Nv[j]*pts[index].y;
        X->z += wts[index]*Nu[i]*Nv[j]*pts[index].z;
        w += wts[index]*Nu[i]*Nv[j];
      }
    }

    // Divide through by the weights
    if (w != 0.0){
      w = 1.0/w;
      X->x *= w;
      X->y *= w;
      X->z *= w;
    }
  }
  else {
    // Evaluate the b-spline
    for ( int j = 0; j < kv; j++ ){
      for ( int i = 0; i < ku; i++ ){
        int index = intu+i + (intv+j)*nu;
        X->x += Nu[i]*Nv[j]*pts[index].x;
        X->y += Nu[i]*Nv[j]*pts[index].y;
        X->z += Nu[i]*Nv[j]*pts[index].z;
      }
    }
  }
  
  // Success
  return 0;
}
  
/*
  Perform the inverse evaluation
*/
int TMRBsplineSurface::invEvalPoint( TMRPoint point, 
                                     double *uf, double *vf ){
  // Did this evaluation fail or not?
  int fail = 1;

  // Set the default parameter values
  *uf = *vf = 0.0;
  
  // The basis functions/work arrays
  double Nu[3*MAX_BSPLINE_ORDER], Nv[3*MAX_BSPLINE_ORDER];
  double work[2*MAX_BSPLINE_ORDER + MAX_BSPLINE_ORDER*MAX_BSPLINE_ORDER];
  double u = 0.5*(Tu[0] + Tu[nu+ku-1]);
  double v = 0.5*(Tv[0] + Tv[nv+kv-1]);

  // Get the bounds
  double umin, vmin, umax, vmax;
  getRange(&umin, &vmin, &umax, &vmax);

  // Perform a newton iteration until convergence
  for ( int k = 0; k < max_newton_iters; k++ ){
    // Compute the knot intervals
    int intu = bspline_interval(u, Tu, nu, ku);
    int intv = bspline_interval(v, Tv, nv, kv);

    // Evaluate the basis functions
    int idu = 2, idv = 2;
    bspline_basis_derivative(Nu, intu, u, idu, Tu, ku, work);
    bspline_basis_derivative(Nv, intv, v, idv, Tv, kv, work);
  
    // Set the interval to the initial control point
    intu = intu - ku + 1;
    intv = intv - kv + 1;

    // Set the point and their derivatives
    TMRPoint X, Xu, Xv, Xuu, Xuv, Xvv;
    X.zero();  Xu.zero();  Xv.zero();
    Xuu.zero();  Xuv.zero();  Xvv.zero();
    
    // If this is a NURBS surface add the effect of the weights
    if (wts){
      // Evaluate the weighted b-spline
      TMRPoint p, pu, pv, puu, puv, pvv;
      p.zero();  pu.zero();  pv.zero();
      puu.zero();  puv.zero();  pvv.zero();
      
      double w = 0.0, wu = 0.0, wv = 0.0;
      double wuu = 0.0, wuv = 0.0, wvv = 0.0;
      for ( int j = 0; j < kv; j++ ){
        for ( int i = 0; i < ku; i++ ){
          // Evaluate the numerator and its derivative
          int index = intu+i + (intv+j)*nu;
          p.x += wts[index]*Nu[i]*Nv[j]*pts[index].x;
          p.y += wts[index]*Nu[i]*Nv[j]*pts[index].y;
          p.z += wts[index]*Nu[i]*Nv[j]*pts[index].z;
          
          pu.x += wts[index]*Nu[ku+i]*Nv[j]*pts[index].x;
          pu.y += wts[index]*Nu[ku+i]*Nv[j]*pts[index].y;
          pu.z += wts[index]*Nu[ku+i]*Nv[j]*pts[index].z;
          
          pv.x += wts[index]*Nu[i]*Nv[kv+j]*pts[index].x;
          pv.y += wts[index]*Nu[i]*Nv[kv+j]*pts[index].y;
          pv.z += wts[index]*Nu[i]*Nv[kv+j]*pts[index].z;

          puu.x += wts[index]*Nu[2*ku+i]*Nv[j]*pts[index].x;
          puu.y += wts[index]*Nu[2*ku+i]*Nv[j]*pts[index].y;
          puu.z += wts[index]*Nu[2*ku+i]*Nv[j]*pts[index].z;

          puv.x += wts[index]*Nu[ku+i]*Nv[kv+j]*pts[index].x;
          puv.y += wts[index]*Nu[ku+i]*Nv[kv+j]*pts[index].y;
          puv.z += wts[index]*Nu[ku+i]*Nv[kv+j]*pts[index].z;
          
          pvv.x += wts[index]*Nu[i]*Nv[2*kv+j]*pts[index].x;
          pvv.y += wts[index]*Nu[i]*Nv[2*kv+j]*pts[index].y;
          pvv.z += wts[index]*Nu[i]*Nv[2*kv+j]*pts[index].z;
          
          // Compute the weights and the derivatives of weights
          w += wts[index]*Nu[i]*Nv[j];
          wu += wts[index]*Nu[ku+i]*Nv[j];
          wv += wts[index]*Nu[i]*Nv[kv+j];
          wuu += wts[index]*Nu[2*ku+i]*Nv[j];
          wuv += wts[index]*Nu[ku+i]*Nv[kv+j];
          wvv += wts[index]*Nu[i]*Nv[2*kv+j];
        }
      }
      
      // Divide through by the weights
      if (w != 0.0){
        w = 1.0/w;
        wuu = 2.0*w*w*w*wu*wu - w*w*wuu;
        wuv = 2.0*w*w*w*wu*wv - w*w*wuv;
        wvv = 2.0*w*w*w*wv*wv - w*w*wvv;
        wu = -w*w*wu;
        wv = -w*w*wv;

        X.x = w*p.x;
        X.y = w*p.y;
        X.z = w*p.z;
        
        Xu.x = w*pu.x + wu*p.x;
        Xu.y = w*pu.y + wu*p.y;
        Xu.z = w*pu.z + wu*p.z;
        
        Xv.x = w*pv.x + wv*p.x;
        Xv.y = w*pv.y + wv*p.y;
        Xv.z = w*pv.z + wv*p.z;

        Xuu.x = 2.0*wu*pu.x + w*puu.x + wuu*p.x;
        Xuu.y = 2.0*wu*pu.y + w*puu.y + wuu*p.y;
        Xuu.z = 2.0*wu*pu.z + w*puu.z + wuu*p.z;
        
        Xuv.x = wu*pv.x + w*puv.x + wuv*p.x + wv*pu.x;
        Xuv.y = wu*pv.y + w*puv.y + wuv*p.y + wv*pu.y;
        Xuv.z = wu*pv.z + w*puv.z + wuv*p.z + wv*pu.z;

        Xvv.x = 2.0*wv*pv.x + w*pvv.x + wvv*p.x;
        Xvv.y = 2.0*wv*pv.y + w*pvv.y + wvv*p.y;
        Xvv.z = 2.0*wv*pv.z + w*pvv.z + wvv*p.z;
      }
    }
    else {
      for ( int j = 0; j < kv; j++ ){
        for ( int i = 0; i < ku; i++ ){
          // Evaluate the numerator and its derivative
          int index = intu+i + (intv+j)*nu;
          X.x += Nu[i]*Nv[j]*pts[index].x;
          X.y += Nu[i]*Nv[j]*pts[index].y;
          X.z += Nu[i]*Nv[j]*pts[index].z;

          Xu.x += Nu[ku+i]*Nv[j]*pts[index].x;
          Xu.y += Nu[ku+i]*Nv[j]*pts[index].y;
          Xu.z += Nu[ku+i]*Nv[j]*pts[index].z;
        
          Xv.x += Nu[i]*Nv[kv+j]*pts[index].x;
          Xv.y += Nu[i]*Nv[kv+j]*pts[index].y;
          Xv.z += Nu[i]*Nv[kv+j]*pts[index].z;

          Xuu.x += Nu[2*ku+i]*Nv[j]*pts[index].x;
          Xuu.y += Nu[2*ku+i]*Nv[j]*pts[index].y;
          Xuu.z += Nu[2*ku+i]*Nv[j]*pts[index].z;
        
          Xuv.x += Nu[ku+i]*Nv[kv+j]*pts[index].x;
          Xuv.y += Nu[ku+i]*Nv[kv+j]*pts[index].y;
          Xuv.z += Nu[ku+i]*Nv[kv+j]*pts[index].z;

          Xvv.x += Nu[i]*Nv[2*kv+j]*pts[index].x;
          Xvv.y += Nu[i]*Nv[2*kv+j]*pts[index].y;
          Xvv.z += Nu[i]*Nv[2*kv+j]*pts[index].z;
        }
      }
    }
  
    // Compute the difference between the position on the surface
    // and the point
    TMRPoint r;
    r.x = (X.x - point.x);
    r.y = (X.y - point.y);
    r.z = (X.z - point.z);

    // Compute the residual
    double ru = Xu.dot(r);
    double rv = Xv.dot(r);

    // Compute the elements of the Jacobian matrix
    double Juu = Xuu.dot(r) + Xu.dot(Xu);
    double Juv = Xuv.dot(r) + Xu.dot(Xv);
    double Jvv = Xvv.dot(r) + Xv.dot(Xv);

    double du = 0.0, dv = 0.0;

    // Check for the bounds on u
    if (u <= umin && ru >= 0.0){
      ru = 0.0;
      Juu = 1.0; 
      Juv = 0.0;
    }
    else if (u >= umax && ru <= 0.0){
      ru = 0.0;
      Juu = 1.0;
      Juv = 0.0;
    }
    
    // Check for the bounds on v
    if (v <= vmin && rv >= 0.0){
      rv = 0.0;
      Jvv = 1.0;
      Juv = 0.0;
    }
    else if (v >= vmax && rv <= 0.0){
      rv = 0.0;
      Jvv = 1.0;
      Juv = 0.0;
    }
                  
    // Solve the 2x2 system
    double det = Juu*Jvv - Juv*Juv;
    if (det != 0.0){
      du = (Jvv*ru - Juv*rv)/det;
      dv = (Juu*rv - Juv*ru)/det;
    }
    
    // Compute the updates
    double unew = u - du;
    double vnew = v - dv;

    // Truncate the u/v values to the bounds
    if (unew < umin){ unew = umin; }
    else if (unew > umax){ unew = umax; }
    if (vnew < vmin){ vnew = vmin; }
    else if (vnew > vmax){ vnew = vmax; }

    // Check if the convergence test is satisfied
    if (fabs(r.x) < eps_dist && 
        fabs(r.y) < eps_dist && 
        fabs(r.z) < eps_dist){
      *uf = u;
      *vf = v;
      return 0;
    }

    // Perform the cosine check
    double dotr = r.dot(r);
    double dotu = Xu.dot(Xu);
    double dotv = Xv.dot(Xv);
    if (ru*ru < eps_cosine*eps_cosine*dotu*dotr && 
        rv*rv < eps_cosine*eps_cosine*dotv*dotr){
      *uf = u;
      *vf = v;
      return 0;
    }

    // Update the new parameter values
    u = unew;
    v = vnew;
  }

  // The newton method failed
  return 0;
}

/*
  Given the parametric location, evaluate the first derivative 
*/
int TMRBsplineSurface::evalDeriv( double u, double v, 
                                  TMRPoint *Xu, TMRPoint *Xv ){
  // The basis functions/work arrays
  double Nu[2*MAX_BSPLINE_ORDER], Nv[2*MAX_BSPLINE_ORDER];
  double work[2*MAX_BSPLINE_ORDER + MAX_BSPLINE_ORDER*MAX_BSPLINE_ORDER];

  // Compute the knot intervals
  int intu = bspline_interval(u, Tu, nu, ku);
  int intv = bspline_interval(v, Tv, nv, kv);

  // Evaluate the basis functions
  int idu = 1, idv = 1;
  bspline_basis_derivative(Nu, intu, u, idu, Tu, ku, work);
  bspline_basis_derivative(Nv, intv, v, idv, Tv, kv, work);
  
  // Set the interval to the initial control point
  intu = intu - ku + 1;
  intv = intv - kv + 1;

  // Zero the point location
  Xu->x = Xu->y = Xu->z = 0.0;
  Xv->x = Xv->y = Xv->z = 0.0;

  // If this is a NURBS surface add the effect of the weights
  if (wts){
    // Evaluate the weighted b-spline
    TMRPoint p, pu, pv;
    p.x = p.y = p.z = 0.0;
    pu.x = pu.y = pu.z = 0.0;
    pv.x = pv.y = pv.z = 0.0;

    double w = 0.0, wu = 0.0, wv = 0.0;
    for ( int j = 0; j < kv; j++ ){
      for ( int i = 0; i < ku; i++ ){
        // Evaluate the numerator and its derivative
        int index = intu+i + (intv+j)*nu;
        p.x += wts[index]*Nu[i]*Nv[j]*pts[index].x;
        p.y += wts[index]*Nu[i]*Nv[j]*pts[index].y;
        p.z += wts[index]*Nu[i]*Nv[j]*pts[index].z;

        pu.x += wts[index]*Nu[ku+i]*Nv[j]*pts[index].x;
        pu.y += wts[index]*Nu[ku+i]*Nv[j]*pts[index].y;
        pu.z += wts[index]*Nu[ku+i]*Nv[j]*pts[index].z;
        
        pv.x += wts[index]*Nu[i]*Nv[kv+j]*pts[index].x;
        pv.y += wts[index]*Nu[i]*Nv[kv+j]*pts[index].y;
        pv.z += wts[index]*Nu[i]*Nv[kv+j]*pts[index].z;

        // Compute the weights and the derivatives of weights
        w += wts[index]*Nu[i]*Nv[j];
        wu += wts[index]*Nu[ku+i]*Nv[j];
        wv += wts[index]*Nu[i]*Nv[kv+j];
      }
    }

    // Divide through by the weights
    if (w != 0.0){
      w = 1.0/w;
      wu = -w*w*wu;
      wv = -w*w*wv;

      Xu->x = w*pu.x + wu*p.x;
      Xu->y = w*pu.y + wu*p.y;
      Xu->z = w*pu.z + wu*p.z;

      Xv->x = w*pv.x + wv*p.x;
      Xv->y = w*pv.y + wv*p.y;
      Xv->z = w*pv.z + wv*p.z;
    }
  }
  else {
    for ( int j = 0; j < kv; j++ ){
      for ( int i = 0; i < ku; i++ ){
        // Evaluate the numerator and its derivative
        int index = intu+i + (intv+j)*nu;
        Xu->x += Nu[ku+i]*Nv[j]*pts[index].x;
        Xu->y += Nu[ku+i]*Nv[j]*pts[index].y;
        Xu->z += Nu[ku+i]*Nv[j]*pts[index].z;
        
        Xv->x += Nu[i]*Nv[kv+j]*pts[index].x;
        Xv->y += Nu[i]*Nv[kv+j]*pts[index].y;
        Xv->z += Nu[i]*Nv[kv+j]*pts[index].z;
      }
    }
  }

  // Success
  return 0;
}

/* 
  Create the interpolation curve generating object.

  This is designed so that different parameters can be set internally
  before the interpolated b-spline curve is created.

  input:
  interp:  the interpolation points
  n:       the number of points
*/
TMRCurveInterpolation::TMRCurveInterpolation( const TMRPoint *_interp,
                                              int _ninterp ){
  ninterp = _ninterp;
  nctl = _ninterp;
  interp = new TMRPoint[ ninterp ];
  memcpy(interp, _interp, ninterp*sizeof(TMRPoint));
}

/*
  Deallocate the interpolation object
*/
TMRCurveInterpolation::~TMRCurveInterpolation(){
  delete [] interp;
}

/*
  Determine the knot intervals at which we'll fit the data
  points. This is based on the chord length between interpolation
  points.
*/
void TMRCurveInterpolation::getInterpLoc( double *ubar ){
  // First, find the chord length of the entire set of points
  double d = 0.0;
  for ( int i = 0; i < ninterp-1; i++ ){
    TMRPoint p;
    p.x = interp[i+1].x - interp[i].x;
    p.y = interp[i+1].y - interp[i].y;
    p.z = interp[i+1].z - interp[i].z;
    d += sqrt(p.dot(p));
  }

  // Now compute ubar inteprolation points based on the 
  // chord length distance between the control points
  ubar[0] = 0.0;
  for ( int i = 0; i < ninterp-1; i++ ){
    TMRPoint p;
    p.x = interp[i+1].x - interp[i].x;
    p.y = interp[i+1].y - interp[i].y;
    p.z = interp[i+1].z - interp[i].z;
    ubar[i+1] = ubar[i] + sqrt(p.dot(p))/d;
  }

  // Ensure that the final interpolation point occurs at the end
  // of the interval
  ubar[ninterp-1] = 1.0;
}

/*
  Create the interpolation and return the curve
*/
TMRBsplineCurve* TMRCurveInterpolation::createCurve( int ku ){
  if (ku < 2){ ku = 2; }
  if (ku > MAX_BSPLINE_ORDER){ ku = MAX_BSPLINE_ORDER; }

  // Pick the ubar values - the parametric locations on the b-spline
  // where the specified values will be interpolated
  double *ubar = new double[ ninterp ];
  getInterpLoc(ubar);

  // Set knot/point data that will be allocated
  double *Tu = NULL;
  TMRPoint *pts = NULL;

  if (nctl < ninterp){
    // Allocate the knot vector
    Tu = new double[ nctl+ku ];
    for ( int i = 0; i < ku; i++ ){
      Tu[i] = 0.0;
      Tu[nctl+ku-1-i] = 1.0;
    }

    // Specify the internal knot vectors
    double d = 1.0*(ninterp-1)/(nctl-ku+1);
    for ( int i = ku; i < nctl; i++ ){
      int j = int(floor((i-ku+1)*d));
      double alpha = (i-ku+1)*d - j;
      Tu[i] = alpha*ubar[j+1] + (1.0 - alpha)*ubar[j];
    }

    // The bandwidth of the linear equation
    int upp_bnd = ku-1, low_bnd = ku-1;
    int bnd = upp_bnd + low_bnd + 1;
    int lda = 2*low_bnd + upp_bnd + 1;
    
    // Allocate the banded matrix
    double *A = new double[ lda*nctl ];
    memset(A, 0, lda*nctl*sizeof(double));

    double *rhs = new double[ 3*nctl ];
    memset(rhs, 0, 3*nctl*sizeof(double));

    // Set the values into the A matrix
    for ( int p = 0; p < ninterp; p++ ){
      // Evaluate the B-spline basis functions
      double Nu[MAX_BSPLINE_ORDER];
      double work[2*MAX_BSPLINE_ORDER];
      
      // Compute the b-spline interval
      int intu = bspline_interval(ubar[p], Tu, nctl, ku);
      
      // Compute the b-spline basis
      bspline_basis(Nu, intu, ubar[p], Tu, ku, work);
      
      // Convert to the initial control point
      intu = intu - ku + 1;

      // Add the contributions to the right-hand side
      for ( int jp = 0; jp < ku; jp++ ){
        int i = intu+jp;
        rhs[i] += Nu[jp]*interp[p].x;
        rhs[nctl+i] += Nu[jp]*interp[p].y;
        rhs[2*nctl+i] += Nu[jp]*interp[p].z;
      }
      
      // Add the contribution to the matrix
      for ( int kp = 0; kp < ku; kp++ ){
        int i = intu+kp;
        for ( int jp = 0; jp < ku; jp++ ){
          int j = intu+jp;
          A[bnd-1+i-j + lda*j] += Nu[jp]*Nu[kp];
        }
      }
    }
    
    // Set the first and last point to be fixed to the
    // initial/final location
    for ( int jp = 0; jp < ku; jp++ ){
      double valu = 0.0;
      if (jp == 0){ valu = 1.0; }

      int i = 0;
      int j = jp;
      A[bnd-1+i-j + lda*j] = valu;

      i = nctl-1;
      j = nctl-1-jp;
      A[bnd-1+i-j + lda*j] = valu;
    }

    // Set the right-hand side to fix the initial point
    rhs[0] = interp[0].x;
    rhs[nctl] = interp[0].y;
    rhs[2*nctl] = interp[0].z;

    // Fix the final location
    rhs[nctl-1] = interp[ninterp-1].x;
    rhs[2*nctl-1] = interp[ninterp-1].y;
    rhs[3*nctl-1] = interp[ninterp-1].z;
    
    // Solve the linear system
    int nrhs = 3;
    int *ipiv = new int[ nctl ];
    int info = 0;
    LAPACKdgbsv(&nctl, &upp_bnd, &low_bnd, &nrhs, 
                A, &lda, ipiv, rhs, &nctl, &info);
    
    if (info != 0){
      fprintf(stderr, 
              "TMRCurveInterpolation error failed with LAPACK error %d\n",
              info);
      return NULL;
    }
    
    // Record/store the control point locations
    pts = new TMRPoint[ nctl ];
    for ( int i = 0; i < nctl; i++ ){
      pts[i].x = rhs[i];
      pts[i].y = rhs[nctl+i];
      pts[i].z = rhs[2*nctl+i];
    }
    
    // Free data associated with the linear system
    delete [] A;
    delete [] rhs;
    delete [] ipiv;
  }
  else { // curve_type == INTERPOLATION){
    nctl = ninterp;

    // Allocate the knot vector
    Tu = new double[ nctl+ku ];
    for ( int i = 0; i < ku; i++ ){
      Tu[i] = 0.0;
      Tu[nctl+ku-1-i] = 1.0;
    }
    
    // Specify the interior knot locations
    for ( int i = ku; i < nctl; i++ ){
      Tu[i] = 0.0;
      
      // Average the ubars to obtain the knot locations
      for ( int j = i; j < i+ku-1; j++ ){
        Tu[i] += ubar[j-ku+1];
      }
      Tu[i] = Tu[i]/(ku-1.0);
    }

    // The bandwidth of the linear equation
    int upp_bnd = ku-1, low_bnd = ku-1;
    int bnd = upp_bnd + low_bnd+1;
    int lda = 2*low_bnd + upp_bnd + 1;
    
    // Allocate the banded matrix
    double *A = new double[ lda*nctl ];
    memset(A, 0, lda*nctl*sizeof(double));
    
    // Set the values into the A matrix
    for ( int i = 0; i < nctl; i++ ){
      // Evaluate the B-spline basis functions
      double Nu[MAX_BSPLINE_ORDER];
      double work[2*MAX_BSPLINE_ORDER];
      
      // Compute the b-spline interval
      int intu = bspline_interval(ubar[i], Tu, nctl, ku);
      
      // Compute the b-spline basis
      bspline_basis(Nu, intu, ubar[i], Tu, ku, work);
      
      // Convert to the initial control point
      intu = intu - ku + 1;
      
      for ( int jp = 0; jp < ku; jp++ ){
        // A(KL+KU+1+i-j,j) = A(i,j)
        int j = intu+jp;
        A[bnd-1+i-j + lda*j] = Nu[jp];
      }
    }
    
    // Compute the right-hand-side values - the x/y/z locations
    // of the interpolation points
    double *rhs = new double[ 3*nctl ];
    for ( int i = 0; i < nctl; i++ ){
      rhs[i] = interp[i].x;
      rhs[nctl+i] = interp[i].y;
      rhs[2*nctl+i] = interp[i].z;
    }
    
    // Solve the linear system
    int nrhs = 3;
    int *ipiv = new int[ nctl ];
    int info = 0;
    LAPACKdgbsv(&nctl, &upp_bnd, &low_bnd, &nrhs, 
                A, &lda, ipiv, rhs, &nctl, &info);
    
    if (info != 0){
      fprintf(stderr, 
              "TMRCurveInterpolation error failed with LAPACK error %d\n",
              info);
      return NULL;
    }
    
    // Record/store the control point locations
    pts = new TMRPoint[ nctl ];
    for ( int i = 0; i < nctl; i++ ){
      pts[i].x = rhs[i];
      pts[i].y = rhs[nctl+i];
      pts[i].z = rhs[2*nctl+i];
    }

    // Free data associated with the linear system
    delete [] A;
    delete [] rhs;
    delete [] ipiv;
  }

  // Create the B-spline and return it
  TMRBsplineCurve *curve = new TMRBsplineCurve(nctl, ku, Tu, pts);

  // Free the data that was allocated
  delete [] ubar;
  delete [] Tu;
  delete [] pts;

  // Return the b-spline object
  return curve;
}

/*
  The following function defines a lofting operation
*/
TMRCurveLofter::TMRCurveLofter( TMRBsplineCurve **_curves, 
                                int _num_curves ){
  num_curves = _num_curves;
  curves = new TMRBsplineCurve*[ num_curves ];
  consist = NULL;

  for ( int i = 0; i < num_curves; i++ ){
    curves[i] = _curves[i];
    curves[i]->incref();
  }
}

/*
  Free the curves
*/
TMRCurveLofter::~TMRCurveLofter(){
  for ( int i = 0; i < num_curves; i++ ){
    curves[i]->decref();
  }
  delete [] curves;

  if (consist){
    for ( int i = 0; i < num_curves; i++ ){
      consist[i]->decref();
    }
    delete [] consist;
  }
}

/*
  Create the lofted surface
*/
TMRBsplineSurface* TMRCurveLofter::createSurface( int kv ){
  // Tuncate the kv value to the permissible range
  if (kv < 2){ kv = 2; }
  if (kv > MAX_BSPLINE_ORDER){ kv = MAX_BSPLINE_ORDER; }

  // Allocate data to store the indices and  
  int *index = new int[ num_curves ];
  int *ptr = new int[ num_curves ];
  int *curve_mult = new int[ num_curves ];

  // Retrieve the underlying data from the curves
  int max_size = 0;
  for ( int i = 0; i < num_curves; i++ ){
    int nu, ku;
    curves[i]->getData(&nu, &ku, NULL, NULL, NULL);
    ptr[i] = 0;
    max_size += nu+ku;
  }

  // Set the values for the candidate smallest knot vector
  int count = 0; // Number of knots on the "stack"
  double tmin = 0.0;
  int max_mult = 0;
  
  // Compute the initial curve multiplicities
  for ( int i = 0; i < num_curves; i++ ){
    int k = ptr[i];
    int nu, ku;
    const double *curve_Tu;
    curves[i]->getData(&nu, &ku, &curve_Tu, NULL, NULL);
    while (k < (nu+ku) && 
           curve_Tu[ptr[i]] == curve_Tu[k]){
      k++;
    }

    // Record the knot multiplicity
    curve_mult[i] = k;

    if (i == 0){
      tmin = curve_Tu[0];
      max_mult = curve_mult[i];
      index[0] = 0;
      count = 1;
    }
  }

  // Allocate the maximum possible size for the knot vector
  int len = 0;
  double *Tu = new double[ max_size ];

  // Find a unified knot vector across all input curves
  int done = 0;
  while (!done){
    // Loop over all curves, finding the min knot value and its
    // maximum multiplicity across all curves
    for ( int i = 0; i < num_curves; i++ ){
      // Skip the first entry - we've already add this guy..
      if (i == index[0]){
        continue;
      }

      // Extract the info for this curve
      int nu, ku;
      const double *curve_Tu;
      curves[i]->getData(&nu, &ku, &curve_Tu, NULL, NULL);

      if (ptr[i] < nu+ku){
        if (curve_Tu[ptr[i]] == tmin){
          // Check if we need to update the max multiplicity
          if (max_mult < curve_mult[i]){
            max_mult = curve_mult[i];
          }
          index[count] = i;
          count++;
        }
        else if (curve_Tu[ptr[i]] < tmin){
          tmin = curve_Tu[ptr[i]];
          max_mult = curve_mult[i];
          index[0] = i;
          count = 1;
        }
      }
    }

    // Add the pointer to the unified knot vector
    for ( int i = 0; i < max_mult; i++, len++ ){
      Tu[len] = tmin;
    }

    // Update the pointers/curve multiplicities
    for ( int j = 0; j < count; j++ ){
      int i = index[j];
      ptr[i] += curve_mult[i];

      // Extract the knot information
      int nu, ku;
      const double *curve_Tu;
      curves[i]->getData(&nu, &ku, &curve_Tu, NULL, NULL);

      // Find the multiplicity of the new knot
      int k = ptr[i];
      while (k < nu+ku &&
             curve_Tu[ptr[i]] == curve_Tu[k]){
        k++;
      }

      // Record the curve multiplicity
      curve_mult[i] = k - ptr[i];
    }

    // Check if we're all done or not
    done = 1;
    for ( int i = 0; i < num_curves; i++ ){
      // Extract the knot information
      int nu, ku;
      const double *curve_Tu;
      curves[i]->getData(&nu, &ku, &curve_Tu, NULL, NULL);

      if (ptr[i] < nu+ku){
        // Reset the data for the next candidate smallest knot value
        tmin = curve_Tu[ptr[i]];
        max_mult = curve_mult[i];
        index[0] = i;
        count = 1;

        // We're not done yet...
        done = 0;
        break;
      }
    }
  }

  // Free the data that is no longer required
  delete [] ptr;
  delete [] index;
  delete [] curve_mult;

  // Now, traverse each curve, and match the difference between the
  // unified knot vector and its derivative
  double *Tnew = new double[ len ];
  consist = new TMRBsplineCurve*[ num_curves ];

  // Loop over all of the curves, determining the knots which need
  // to be added to make it consistent with the unifed Tu knot vector
  for ( int i = 0; i < num_curves; i++ ){
    int nu, ku;
    const double *curve_Tu;
    curves[i]->getData(&nu, &ku, &curve_Tu, NULL, NULL);

    // Keep track of the number of new knots required for each curve
    int nnew = 0;
    for ( int j = 0, k = 0; j < len; j++ ){
      if (Tu[j] == curve_Tu[k]){
        k++;
      }
      else {
        Tnew[nnew] = Tu[j];
        nnew++;
      }
    }

    // Refine the knot vector for each curve to make it consistent
    consist[i] = curves[i]->refineKnots(Tnew, nnew);
    consist[i]->incref();
  }

  // Keep track of the maximum value of kv - this must be consistent
  // with the number of specified curves
  int nv = num_curves;
  if (kv > num_curves){ kv = num_curves; }

  // Create a knot vector for the v-direction
  double *Tv = new double[ nv+kv ];
  bspline_knots(nv, kv, 0.0, 1.0, Tv);

  // Check whether any of the curves have the rational part
  int is_nurbs = 0;
  for ( int i = 0; i < num_curves; i++ ){
    const double *curve_wts;
    consist[i]->getData(NULL, NULL, NULL, &curve_wts, NULL);
    if (curve_wts){
      is_nurbs = 1;
      break;
    }
  }

  // Create a TMRPoint array containing all of the poitns
  int ku = 4;
  int nu = len-ku;
  TMRPoint *pts = new TMRPoint[ nu*nv ];
  double *wts = NULL;
  if (is_nurbs){
    wts = new double[ nu*nv ];
  }

  // Get the order of the curve
  consist[0]->getData(NULL, &ku, NULL, NULL, NULL);

  // Go through and extract the points
  for ( int j = 0; j < num_curves; j++ ){
    const double *curve_wts;
    const TMRPoint *curve_pts;
    consist[j]->getData(NULL, NULL, NULL, &curve_wts, &curve_pts);

    // Copy over the weights - if any
    if (wts){
      if (curve_wts){
        for ( int i = 0; i < nu; i++ ){
          wts[j*nu + i] = curve_wts[i];
        }
      }
      else {
        for ( int i = 0; i < nu; i++ ){
          wts[j*nu + i] = 1.0;
        }
      }
    }

    // Copy over the points
    for ( int i = 0; i < nu; i++ ){
      pts[j*nu + i] = curve_pts[i];
    }
  }

  // Create the surface
  TMRBsplineSurface *surf = NULL;
  if (wts){
    surf = new TMRBsplineSurface(nu, nv, ku, kv, Tu, Tv, wts, pts);
  }
  else {
    surf = new TMRBsplineSurface(nu, nv, ku, kv, Tu, Tv, pts);
  }

  // Free the data that is no longer required
  delete [] Tu;
  delete [] Tv;
  delete [] Tnew;
  delete [] pts;
  if (wts){ delete [] wts; }

  return surf;
}
