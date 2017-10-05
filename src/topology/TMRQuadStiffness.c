#include "TMRQuadStiffness.h"

/*
  Create the octree stiffness object based on an interpolation from
  the filter variables
*/
TMRQuadStiffness::TMRQuadStiffness( TMRIndexWeight *_weights,
                                    int _nweights,
                                    TacsScalar _density,
                                    TacsScalar _E, TacsScalar _nu,
                                    double _q, double _eps ){
  // Record the density, Poisson ratio
  density = _density;
  nu = _nu;
  q = _q;
  eps = _eps;
  E = _E;
  // Set the weights/local indices
  nweights = _nweights;
  memcpy(weights, _weights, nweights*sizeof(TMRIndexWeight));
  
  // Set the initial value for the densities
  rho = 0.95;
}

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
void TMRQuadStiffness::setDesignVars( const TacsScalar x[], int numDVs ){
  rho = 0.0;
  for ( int i = 0; i < nweights; i++ ){
    rho += weights[i].weight*x[weights[i].index];
  }
}

/*
  Get the design variable values

  This is not possible to back out once the density is computed, so
  instead we use a constant value of 0.95 as the output.
*/
void TMRQuadStiffness::getDesignVars( TacsScalar x[], int numDVs ){
  for ( int i = 0; i < nweights; i++ ){
    if (weights[i].index >= 0 && weights[i].index < numDVs){
      x[weights[i].index] = 0.95;
    }
  }
}

/*
  Retrieve the range of the design variable values
*/
void TMRQuadStiffness::getDesignVarRange( TacsScalar lb[], 
                                         TacsScalar ub[], int numDVs ){
  for ( int i = 0; i < nweights; i++ ){
    if (weights[i].index >= 0 && weights[i].index < numDVs){
      lb[weights[i].index] = 0.0;
      ub[weights[i].index] = 1.0;
    }
  }
}

/*
  Given the values of everything, compute the stress
*/
void TMRQuadStiffness::calculateStress( const double pt[],
                                        const TacsScalar e[], 
                                        TacsScalar s[] ){
  // Compute the penalized stiffness
  TacsScalar w = rho/(1.0 + q*(1.0 - rho));
  TacsScalar D = E/(1.0-nu*nu)*w;
  
  s[0] = D*e[0]+nu*D*e[1];
  s[1] = nu*D*e[0]+D*e[1];
  s[2] = 0.5*(1.0-nu)*D*e[2];
}

/*
  Add the derivative of the product of the stress with the vector psi
  times alpha to the design variable array dvSens
*/
void TMRQuadStiffness::addStressDVSens( const double pt[], 
                                        const TacsScalar e[], 
                                        TacsScalar alpha, 
                                        const TacsScalar psi[], 
                                        TacsScalar fdvSens[], int dvLen ){
  TacsScalar dxw = 1.0+q*(1.0-rho);
  TacsScalar w = ((1.0+q)/(dxw*dxw));
  TacsScalar s[3];
  TacsScalar D = E/(1.0-nu*nu)*w;

  // Compute the product dD/dx*(B*u)
  s[0] = D*e[0]+nu*D*e[1];
  s[1] = nu*D*e[0]+D*e[1];
  s[2] = 0.5*(1.0-nu)*D*e[2];
  // Compute the term psi^{T}*B^{T}*dD/dx*B*u
  for ( int i = 0; i < nweights; i++ ){
    fdvSens[weights[i].index] += 
      alpha*weights[i].weight*(s[0]*psi[0] + s[1]*psi[1] 
                               + s[2]*psi[2]); 
  }
}

/*
  Compute the density of the material as a linear function of the
  element-wise density
*/
void TMRQuadStiffness::getPointwiseMass( const double pt[], 
                                         TacsScalar mass[] ){
  mass[0] = rho*density;
}

/*
  Compute the derivative of the pointwise mass as function of the design
  variable values
*/
void TMRQuadStiffness::addPointwiseMassDVSens( const double pt[], 
                                               const TacsScalar alpha[],
                                               TacsScalar dvSens[], 
                                               int dvLen ){
  TacsScalar scale = density*alpha[0];
  for ( int i = 0; i < nweights; i++ ){
    dvSens[weights[i].index] += scale*weights[i].weight;
  }
}
