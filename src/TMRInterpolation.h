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

#ifndef TMR_INTERPOLATION_FUNCTIONS_H
#define TMR_INTERPOLATION_FUNCTIONS_H

/*
  The following file defines the inline interpolation functions used
  by TMR.
*/

/*
  Evaluate the shape functions at the given parametric point

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  knots:  the interpolation knots in parameter space

  output:
  N:      the values of the shape functions at u
*/
inline void lagrange_shape_functions( const int order,
                                      const double u,
                                      const double *knots,
                                      double *N ){
  // Loop over the shape functions
  for ( int i = 0; i < order; i++ ){
    N[i] = 1.0;
    for ( int j = 0; j < order; j++ ){
      if (i != j){
        double d = 1.0/(knots[i] - knots[j]);
        N[i] *= (u - knots[j])*d;
      }
    }
  }
}

/*
  Evaluate the shape functions and the derivative of the shape functions
  with respect to the parameter coordinate

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  knots:  the interpolation knots in parameter space

  output:
  N:      the values of the shape functions at u
  Nd:     the derivative of the shape functions at u
*/
inline void lagrange_shape_func_derivative( const int order,
                                            const double u,
                                            const double *knots,
                                            double *N,
                                            double *Nd ){
  // Loop over the shape function knot locations
  for ( int i = 0; i < order; i++ ){
    N[i] = 1.0;
    Nd[i] = 0.0;

    // Loop over each point again, except for the current control
    // point, adding the contribution to the shape function
    for ( int j = 0; j < order; j++ ){
      if (i != j){
        double d = 1.0/(knots[i] - knots[j]);
        N[i] *= (u - knots[j])*d;

        // Now add up the contribution to the derivative
        for ( int k = 0; k < order; k++ ){
          if (k != i && k != j){
            d *= (u - knots[k])/(knots[i] - knots[k]);
          }
        }

        // Add the derivative contribution
        Nd[i] += d;
      }
    }
  }
}

/*
  Evaluate the shape functions and their first and second derivatives
  with respect to the parameter coordinate

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  knots:  the interpolation knots in parameter space

  output:
  N:      the values of the shape functions at u
  Nd:     the derivative of the shape functions at u
  Ndd:    the second derivative of the shape functions at u
*/
inline void lagrange_shape_func_second_derivative( const int order,
                                                   const double u,
                                                   const double *knots,
                                                   double *N,
                                                   double *Nd,
                                                   double *Ndd ){
 // Loop over the shape function control points
  for ( int i = 0; i < order; i++ ){
    N[i] = 1.0;
    Nd[i] = 0.0;
    Ndd[i] = 0.0;

    // Loop over each point again, except for the current control
    // point, adding the contribution to the shape function
    for ( int j = 0; j < order; j++ ){
      if (i != j){
        double tj = 1.0/(knots[i] - knots[j]);
        double dj = tj;
        N[i] = N[i]*(u - knots[j])*dj;

        // Loop over all the knots again to determine the
        // contribution to the derivative of the shape function
        for ( int k = 0; k < order; k++ ){
          if (k != i && k != j){
            double dk = 1.0/(knots[i] - knots[k]);
            dj *= (u - knots[k])*dk;
            dk *= tj;

            for ( int m = 0; m < order; m++ ){
              if (m != i && m != j && m != k){
                dk *= (u - knots[m])/(knots[i] - knots[m]);
              }
            }
            Ndd[i] += dk;
          }
        }
        Nd[i] += dj;
      }
    }
  }
}

/*
  Evaluate the Bernstein shape functions at the given parametric point

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  idx:    the index of the interval for u
  knots:  the interpolation knots in parameter space
  work:   a temporary array of size 2*k

  output:
  N:      the values of the shape functions at u
*/
inline void bernstein_shape_functions( const int order,
                                       const double u,
                                       const double *knots,
                                       double *N ){
  N[0] = 1.0;
  int idx = order-1;
  // Set the pointers for the temporary work arrays
  // Note that left[j] = u - knots[i+1 - j]
  // and right[j] = knots[i+j] - u
  double *left = new double[ order ];
  double *right = new double[order];
  
  for ( int j = 1; j < order; j++ ){
    left[j] = u - knots[idx+1-j];
    right[j] = knots[idx+j] - u;
    
    N[j] = 0.0;
    for ( int i = 0; i < j; i++ ){
      double temp = N[i]/(right[i+1] + left[j-i]);
      N[i] = N[j] + right[i+1]*temp;
      N[j] = left[j-i]*temp;
    }
  }    
  delete [] left;
  delete [] right;
}

/*
  Evaluate the shape functions and the derivative of the shape functions
  with respect to the parameter coordinate

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  knots:  the interpolation knots in parameter space

  output:
  N:      the values of the shape functions at u
  Nd:     the derivative of the shape functions at u
*/
inline void bernstein_shape_func_derivative( const int order,
                                             const double u,
                                             const double *knots,
                                             double *N,
                                             double *Nd ){
  int ideriv = 1;
  int ku = order;
  const int idx = order-1;
  double *work = new double[4*ku];
  N[0] = 1.0;

  // Set the pointers for the temporary work arrays
  // Note that left[j] = u - knots[i+1 - j]
  // and right[j] = knots[i+j] - u
  double *left = &work[0];
  double *right = &work[ku];

  // ndu is a matrix of total size ku*ku
  // The basis functions are stored in the upper triangular part of the matrix
  // such that N_{idx-j, i} = ndu[i + j*ku] if i >= j is the basis function.
  // The knot differences are stored in the lower portion of the matrix such that
  // ndu[i + j*ku] = u_{idx+i+1} - u_{idx+j-i}
  double *ndu = &work[2*ku];
  ndu[0] = 1.0;

  // Compute all orders of the basis functions
  for ( int j = 1; j < ku; j++ ){
    left[j] = u - knots[idx+1-j];
    right[j] = knots[idx+j] - u;
    // Compute the shape functions
    N[j] = 0.0;
    for ( int i = 0; i < j; i++ ){
      double temp = N[i]/(right[i+1] + left[j-i]);
      N[i] = N[j] + right[i+1]*temp;
      N[j] = left[j-i]*temp;
    }
    // Compute the terms needed for the 1st and 2nd derivative
    double njj = 0.0;
    for ( int i = 0; i < j; i++ ){
      // Compute the difference knots[idx+i+1] - knots[idx+j-i]
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
    Nd[i] = ndu[(ku-1) + i*ku]; 
  }

  // Set the temporary arrays for the a-coefficients
  double *a0 = &work[0];
  double *a1 = &work[ku];

  for ( int i = 0; i < ku; i++ ){
    a0[0] = 1.0;
    double d = 0.0;

    // Compute the first of the a-terms
    // a_{k,0} = a_{k-1,0}/(u_{i+ku-k} - u_{i})
    if (i >= ideriv){
      a1[0] = a0[0]/ndu[(i-ideriv) + (ku-ideriv)*ku];
      d = a1[0]*ndu[(ku-ideriv-1) + (i-ideriv)*ku];
    }
    
    int jstart = ideriv-i;
    if (i >= ideriv-1){ jstart = 1; }
    
    int jend = ku-i-1;
    if (i <= ku-ideriv){ jend = ideriv-1; }
    
    for ( int j = jstart; j <= jend; j++ ){
      // Compute a_{k,j} = (a_{k-1,j} - a_{k-1,j-1})/(u_{i+ku+j-k} - u_{i+j})
      a1[j] = (a0[j] - a0[j-1])/ndu[(i-ideriv+j) + (ku-ideriv)*ku];
      d += a1[j]*ndu[(ku-ideriv-1) + (i-ideriv+j)*ku];
    }
    
    // Compute the term 
    // a_{k,k} = -a_{k-1}/(u_{i+ku} - u_{i+k})
    if (i <= ku-ideriv-1){
      a1[ideriv] = -a0[ideriv-1]/ndu[i + (ku-ideriv)*ku];
      d += a1[ideriv]*ndu[(ku-ideriv-1) + i*ku];
    }
    
    // Set the basis function
    Nd[i + ideriv*ku] = d;
    
    // Swap the rows of a
    double *t = a0;
    a0 = a1;  a1 = t;
  }

  // Multiply the basis by the factorial term
  int r = ku-1;
  for ( int j = 0; j < ku; j++ ){
    Nd[j + ideriv*ku] *= r;
  }
}

/*
  Evaluate the shape functions and their first and second derivatives
  with respect to the parameter coordinate

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  knots:  the interpolation knots in parameter space

  output:
  N:      the values of the shape functions at u
  Nd:     the derivative of the shape functions at u
  Ndd:    the second derivative of the shape functions at u
*/
inline void bernstein_shape_func_second_derivative( const int order,
                                                    const double u,
                                                    const double *knots,
                                                    double *N,
                                                    double *Nd,
                                                    double *Ndd ){
  int ku = order;
  const int idx = order-1;
  double *work = new double[4*ku];
  N[0] = 1.0;

  // Set the pointers for the temporary work arrays
  // Note that left[j] = u - knots[i+1 - j]
  // and right[j] = knots[i+j] - u
  double *left = &work[0];
  double *right = &work[ku];

  // ndu is a matrix of total size ku*ku
  // The basis functions are stored in the upper triangular part of the matrix
  // such that N_{idx-j, i} = ndu[i + j*ku] if i >= j is the basis function.
  // The knot differences are stored in the lower portion of the matrix such that
  // ndu[i + j*ku] = u_{idx+i+1} - u_{idx+j-i}
  double *ndu = &work[2*ku];
  ndu[0] = 1.0;

  // Compute all orders of the basis functions
  for ( int j = 1; j < ku; j++ ){
    left[j] = u - knots[idx+1-j];
    right[j] = knots[idx+j] - u;
    // Compute the shape functions
    N[j] = 0.0;
    for ( int i = 0; i < j; i++ ){
      double temp = N[i]/(right[i+1] + left[j-i]);
      N[i] = N[j] + right[i+1]*temp;
      N[j] = left[j-i]*temp;
    }
    // Compute the second derivative
    double njj = 0.0;
    for ( int i = 0; i < j; i++ ){
      // Compute the difference knots[idx+i+1] - knots[idx+j-i]
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
    Nd[i] = ndu[(ku-1) + i*ku];
    Ndd[i] = ndu[(ku-1) + i*ku];
  }

  // Set the temporary arrays for the a-coefficients
  double *a0 = &work[0];
  double *a1 = &work[ku];
  int n_deriv[] = {1,2};
  // Compute the nth derivative
  for (int n = 0; n < 2; n++){
    for (int i = 0; i < ku; i++){
      a0[0] = 1.0;
      
      for ( int k = 1; k <= n_deriv[n]; k++ ){
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
          // Compute a_{k,j} = (a_{k-1,j} - a_{k-1,j-1})/(u_{i+ku+j-k} - u_{i+j})
          a1[j] = (a0[j] - a0[j-1])/ndu[(i-k+j) + (ku-k)*ku];
          d += a1[j]*ndu[(ku-k-1) + (i-k+j)*ku];
        }

        // Compute the term
        // a_{k,k} = -a_{k-1}/(u_{i+ku} - u_{i+k})
        if (i <= ku-k-1){
          a1[k] = -a0[k-1]/ndu[i + (ku-k)*ku];
          d += a1[k]*ndu[(ku-k-1) + i*ku];
        }
        if (n_deriv[n] == 1){
          // Set the basis function
          Nd[i + k*ku] = d;
        }
        else{
          // Set the basis function
          Ndd[i + k*ku] = d;
        }
        // Swap the rows of a
        double *t = a0;
        a0 = a1;  a1 = t;
      }
      // Multiply the basis by the factorial term
      int r = ku-1;
      for ( int k = 1; k <= n_deriv[n]; k++ ){
        for ( int j = 0; j < ku; j++ ){
          if (n_deriv[n] == 1){
          // Set the basis function
            Nd[j + k*ku] *= r;
          }
          else{
            Ndd[j + k*ku] *= r;
          }
        }
        r *= (ku-1 - k);
      }
    }
  } 
}

#endif // TMR_INTERPOLATION_FUNCTIONS_H
