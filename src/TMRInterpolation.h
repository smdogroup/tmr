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
