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

#endif // TMR_INTERPOLATION_FUNCTIONS_H
