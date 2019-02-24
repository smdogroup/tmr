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
                                       double *N ){
  double u1 = 0.5*(1.0 - u);
  double u2 = 0.5*(u + 1.0);

  N[0] = 1.0;
  for ( int j = 1; j < order; j++ ){
    double s = 0.0;
    for ( int k = 0; k < j; k++ ){
      double t = N[k];
      N[k] = s + u1*t;
      s = u2*t;
    }
    N[j] = s;
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
inline void bernstein_shape_func_derivative( const int order,
                                             const double u,
                                             double *N,
                                             double *Nd ){
  double u1 = 0.5*(1.0 - u);
  double u2 = 0.5*(u + 1.0);

  // Compute the basis for the order-1 bernstein polynomial
  N[0] = 1.0;
  for ( int j = 1; j < order-1; j++ ){
    double s = 0.0;
    for ( int k = 0; k < j; k++ ){
      double t = N[k];
      N[k] = s + u1*t;
      s = u2*t;
    }
    N[j] = s;
  }

  // Add the contributions to the derivative
  for ( int j = 0; j < order; j++ ){
    Nd[j] = 0.0;
    if (j > 0){
      Nd[j] += 0.5*(order-1)*N[j-1];
    }
    if (j < order-1){
      Nd[j] -= 0.5*(order-1)*N[j];
    }
  }

  // Now compute the full order basis
  bernstein_shape_functions(order, u, N);
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
                                                    double *N,
                                                    double *Nd,
                                                    double *Ndd ){
  double u1 = 0.5*(1.0 - u);
  double u2 = 0.5*(u + 1.0);

  // Compute the basis for the order-1 bernstein polynomial
  N[0] = 1.0;
  for ( int j = 1; j < order-2; j++ ){
    double s = 0.0;
    for ( int k = 0; k < j; k++ ){
      double t = N[k];
      N[k] = s + u1*t;
      s = u2*t;
    }
    N[j] = s;
  }

  // Add the contributions to the derivative
  for ( int j = 0; j < order-1; j++ ){
    Nd[j] = 0.0;
    if (j > 0){
      Nd[j] += 0.5*(order-2)*N[j-1];
    }
    if (j < order-2){
      Nd[j] -= 0.5*(order-2)*N[j];
    }
  }

  for ( int j = 0; j < order; j++ ){
    Ndd[j] = 0.0;
    if (j > 0){
      Ndd[j] += 0.5*(order-1)*Nd[j-1];
    }
    if (j < order-1){
      Ndd[j] -= 0.5*(order-1)*Nd[j];
    }
  }

  bernstein_shape_func_derivative(order, u, N, Nd);
}

/*
  Evaluate the constraint weights for a Bernstein polynomial along an
  edge for the given mesh order.

  This code relies on the
  o---------------o
  |               |
  |               |
  |               |
  | <<<< Ec >>>>> |
  o-------x-------o
  |       |       |
  |       |       |
  o-------o-------o

  The basis functions on Ec are indexed from -(mesh_order-1) to
  (mesh_order-1). The input argument u gives the index of the basis 
  function on the two adjacent refined elements.

  input:
  mesh_order:  the order of the mesh
  u:           the index of the basis function on the refined edge

  output:
  N:           the weights on the basis functions on the coarse edge 
*/
inline void eval_bernstein_weights( const int mesh_order,
                                    const int u,
                                    double *N ){
  // Weights evaluated from the subdivison matrix
  if (mesh_order == 2){
    if (u == -1){
      N[0] = 1.0;
      N[1] = 0.0;
    }
    else if (u == 0){
      N[0] = 0.5;
      N[1] = 0.5;
    }
    else { // u == 1
      N[0] = 0.0;
      N[1] = 1.0;
    }
  }
  else if (mesh_order == 3){
    if (u == -2){
      N[0] = 1.0;
      N[1] = 0.0;
      N[2] = 0.0;
    }
    else if (u == -1){
      N[0] = 0.5;
      N[1] = 0.5;
      N[2] = 0.0;
    }
    else if (u == 0){
      N[0] = 0.25;
      N[1] = 0.5;
      N[2] = 0.25;
    }
    else if (u == 1){
      N[0] = 0.0;
      N[1] = 0.5;
      N[2] = 0.5;
    }
    else { // u == 2
      N[0] = 0.0;
      N[1] = 0.0;
      N[2] = 1.0;
    }
  }
  else if (mesh_order == 4){
    if (u == -3){
      N[0] = 1.0;
      N[1] = 0.0;
      N[2] = 0.0;
      N[3] = 0.0;
    }
    else if (u == -2){
      N[0] = 0.5;
      N[1] = 0.5;
      N[2] = 0.0;
      N[3] = 0.0;
    }
    else if (u == -1){
      N[0] = 0.25;
      N[1] = 0.5;
      N[2] = 0.25;
      N[3] = 0.0;
    }
    else if (u == 0){
      N[0] = 0.125;
      N[1] = 0.375;
      N[2] = 0.375;
      N[3] = 0.125;
    }
    else if (u == 1){
      N[0] = 0.0;
      N[1] = 0.25;
      N[2] = 0.5;
      N[3] = 0.25;
    }
    else if (u == 2){
      N[0] = 0.0;
      N[1] = 0.0;
      N[2] = 0.5;
      N[3] = 0.5;
    }
    else { // u == 3
      N[0] = 0.0;
      N[1] = 0.0;
      N[2] = 0.0;
      N[3] = 1.0;
    }
  }
  else if (mesh_order == 5){
    if (u == -4){
      N[0] = 1.0;
      N[1] = 0.0;
      N[2] = 0.0;
      N[3] = 0.0;
      N[4] = 0.0;
    }
    else if (u == -3){
      N[0] = 0.5;
      N[1] = 0.5;
      N[2] = 0.0;
      N[3] = 0.0;
      N[4] = 0.0;
    }
    else if (u == -2){
      N[0] = 0.25;
      N[1] = 0.5;
      N[2] = 0.25;
      N[3] = 0.0;
      N[4] = 0.0;
    }
    else if (u == -1){
      N[0] = 0.125;
      N[1] = 0.375;
      N[2] = 0.375;
      N[3] = 0.125;
      N[4] = 0.0;
    }
    else if (u == 0){
      N[0] = 0.0625;
      N[1] = 0.25;
      N[2] = 0.375;
      N[3] = 0.25;
      N[4] = 0.0625;
    }
    else if (u == 1){
      N[0] = 0.0;
      N[1] = 0.125;
      N[2] = 0.375;
      N[3] = 0.375;
      N[4] = 0.125;
    }
    else if (u == 2){
      N[0] = 0.0;
      N[1] = 0.0;
      N[2] = 0.25;
      N[3] = 0.5;
      N[4] = 0.25;
    }
    else if (u == 3){
      N[0] = 0.0;
      N[1] = 0.0;
      N[2] = 0.0;
      N[3] = 0.5;
      N[4] = 0.5;
    }
    else { // u == 4
      N[0] = 0.0;
      N[1] = 0.0;
      N[2] = 0.0;
      N[3] = 0.0;
      N[4] = 1.0;
    }
  }
}

/*
  Evaluate the weights on the coarse mesh that interpolate the 
  Bernstein polynomial from the fine mesh.

  Note: This only works for cases in which there is at most one-order
  difference between the meshes.

  input:
  mesh_order:         the order on the fine space
  coarse_mesh_order:  the order on the coarse space
  u:                  the basis function index (on the fine space)

  output:
  N:                  the weights in the coarse space
*/
inline int eval_bernstein_interp_weights( const int mesh_order,
                                          const int coarse_mesh_order,
                                          const int u,
                                          double *N ){
  // Weights evaluated from the subdivison matrix
  if (mesh_order - coarse_mesh_order == 1){
    if (coarse_mesh_order == 2){
      if (u == 0){
        N[0] = 1.0;
        N[1] = 0.0;
      }
      else if (u == 1){
        N[0] = 0.5;
        N[1] = 0.5;
      }
      else if (u == 2){
        N[0] = 0.0;
        N[1] = 1.0;
      }
    }
    else if (coarse_mesh_order == 3){
      if (u == 0){
        N[0] = 1.0;
        N[1] = 0.0;
        N[2] = 0.0;
      }
      else if (u == 1){
        N[0] = 1./3;
        N[1] = 2./3;
        N[2] = 0.0;
      }
      else if (u == 2){
        N[0] = 0.0;
        N[1] = 2./3;
        N[2] = 1./3;
      }
      else if (u == 3){
        N[0] = 0.0;
        N[1] = 0.0;
        N[2] = 1.0;
      }
    }
    else if (coarse_mesh_order == 4){
      if (u == 0){
        N[0] = 1.0;
        N[1] = 0.0;
        N[2] = 0.0;
        N[3] = 0.0;        
      }
      else if (u == 1){
        N[0] = 0.25;
        N[1] = 0.75;
        N[2] = 0.0;
        N[3] = 0.0;
      }
      else if (u == 2){
        N[0] = 0.0;
        N[1] = 0.0;
        N[2] = 0.75;
        N[3] = 0.25;
      }
      else if (u == 3){
        N[0] = 0.0;
        N[1] = 0.5;
        N[2] = 0.5;
        N[3] = 0.0;
      }
      else if (u == 4){
        N[0] = 0.0;
        N[1] = 0.0;
        N[2] = 0.0;
        N[3] = 1.0;        
      }
    }
    else if (coarse_mesh_order == 5){
      if (u == 0){
        N[0] = 1.0;
        N[1] = 0.0;
        N[2] = 0.0;
        N[3] = 0.0;
        N[4] = 0.0;
      }
      else if (u == 1){
        N[0] = 0.2;
        N[1] = 0.8;
        N[2] = 0.0;
        N[3] = 0.0;
        N[4] = 0.0;
      }
      else if (u == 2){
        N[0] = 0.0;
        N[1] = 0.4;
        N[2] = 0.6;
        N[3] = 0.0;
        N[4] = 0.0;
      }
      else if (u == 3){
        N[0] = 0.0;
        N[1] = 0.0;
        N[2] = 0.6;
        N[3] = 0.4;
        N[4] = 0.0;
      }
      else if (u == 4){
        N[0] = 0.0;
        N[1] = 0.0;
        N[2] = 0.0;
        N[3] = 0.8;
        N[4] = 0.2;      
      }
      else if (u == 5){
        N[0] = 0.0;
        N[1] = 0.0;
        N[2] = 0.0;
        N[3] = 0.0;
        N[4] = 1.0;
      }
    }

    return 0;
  }
  else if (mesh_order == coarse_mesh_order){
    memset(N, 0, coarse_mesh_order*sizeof(double));
    N[u] = 1.0;
    
    return 0;
  }

  // Failed
  return 1;
}

#endif // TMR_INTERPOLATION_FUNCTIONS_H
