#include "TMROctStiffness.h"

/*
  Create the octree stiffness object based on an interpolation from
  the filter variables
*/
TMROctStiffness::TMROctStiffness( TMRIndexWeight *_weights,
                                  int _nweights,
                                  TacsScalar _density,
                                  TacsScalar E, TacsScalar _nu,
                                  double _q ){
  // Record the density, Poisson ratio, D and the shear modulus
  density = _density;
  nu = _nu;
  D = E/((1.0 + nu)*(1.0 - 2.0*nu));
  G = 0.5*E/(1.0 + nu);
  q = _q;

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
void TMROctStiffness::setDesignVars( const TacsScalar x[], int numDVs ){
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
void TMROctStiffness::getDesignVars( TacsScalar x[], int numDVs ){
  for ( int i = 0; i < nweights; i++ ){
    if (weights[i].index >= 0 && weights[i].index < numDVs){
      x[weights[i].index] = 0.95;
    }
  }
}

/*
  Retrieve the range of the design variable values
*/
void TMROctStiffness::getDesignVarRange( TacsScalar lb[], 
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
void TMROctStiffness::calculateStress( const double pt[],
                                       const TacsScalar e[], 
                                       TacsScalar s[] ){
  // Set the minimum relative stiffness value
  const double eps = 1e-3;

  // Compute the penalized stiffness
  TacsScalar penalty = rho/(1.0 + q*(1.0 - rho));
  TacsScalar Dp = (penalty + eps)*D;
  TacsScalar Gp = (penalty + eps)*G;
  s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
  s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
  s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
  s[3] = Gp*e[3];
  s[4] = Gp*e[4];
  s[5] = Gp*e[5];
}

/*
  Add the derivative of the product of the stress with the vector psi
  times alpha to the design variable array dvSens
*/
void TMROctStiffness::addStressDVSens( const double pt[], 
                                       const TacsScalar e[], 
                                       TacsScalar alpha, 
                                       const TacsScalar psi[], 
                                       TacsScalar fdvSens[], int dvLen ){
  TacsScalar penalty = (q + 1.0)/((1.0 + q*(1.0 - rho))*(1.0 + q*(1.0 - rho)));
  TacsScalar Dp = alpha*penalty*D;
  TacsScalar Gp = alpha*penalty*G;
  TacsScalar s[6];
  s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
  s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
  s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
  s[3] = Gp*e[3];
  s[4] = Gp*e[4];
  s[5] = Gp*e[5];

  for ( int i = 0; i < nweights; i++ ){
    fdvSens[weights[i].index] += 
      weights[i].weight*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] + 
			 s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
  }
}

/*
  Compute the density of the material as a linear function of the
  element-wise density
*/
void TMROctStiffness::getPointwiseMass( const double pt[], 
                                        TacsScalar mass[] ){
  mass[0] = rho*density;
}

/*
  Compute the derivative of the pointwise mass as function of the design
  variable values
*/
void TMROctStiffness::addPointwiseMassDVSens( const double pt[], 
                                              const TacsScalar alpha[],
                                              TacsScalar dvSens[], 
                                              int dvLen ){
  TacsScalar scale = density*alpha[0];
  for ( int i = 0; i < nweights; i++ ){
    dvSens[weights[i].index] += scale*weights[i].weight;
  }
}
