#include "TMRQuadStiffness.h"
#include "FElibrary.h"
#include "YSlibrary.h"
/*
  Create the quadtree stiffness object based on an interpolation from
  the filter variables
*/
TMRQuadStiffness::TMRQuadStiffness( TMRIndexWeight *_weights,
                                    int _nweights,
                                    TacsScalar _density,
                                    TacsScalar _E,
                                    TacsScalar _nu,
                                    TacsScalar _ys,
                                    double _q, double _eps ){
  // Record the density, Poisson ratio, D and the shear modulus
  density = _density;
  nu = _nu;
  q = _q;
  eps = _eps;
  E = _E;
  ys = _ys;
  // Set the weights/local indices
  nweights = _nweights;
  memcpy(weights, _weights, nweights*sizeof(TMRIndexWeight));
 
  // Set the initial value for the densities
  xw = 0.50;
}

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
void TMRQuadStiffness::setDesignVars( const TacsScalar x[], int numDVs ){
  xw = 0.0;
  for ( int i = 0; i < nweights; i++ ){
    xw += weights[i].weight*x[weights[i].index];
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
                                          TacsScalar ub[], 
                                          int numDVs ){
  for ( int i = 0; i < nweights; i++ ){
    if (weights[i].index >= 0 && weights[i].index < numDVs){
      lb[weights[i].index] = 1e-3;
      ub[weights[i].index] = 1.0;
    }
  }
}
/*
  Compute the density of the material as a linear function of the
  element-wise density
*/
void TMRQuadStiffness::getPointwiseMass( const double pt[], 
                                         TacsScalar mass[] ){
  mass[0] = xw*density;
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

/*
  Compute the matrix vector product of D*e
*/
void TMRQuadStiffness::calculateStress(const double pt[],
                                       const TacsScalar e[], 
                                       TacsScalar s[]){
  // Compute the penalized stiffness
  TacsScalar w = xw/(1.0 + q*(1.0 - xw));
  TacsScalar D = E/(1.0-nu*nu)*w;

  // Compute the resulting matrix-vector product D*e
  s[0] = D*e[0]+nu*D*e[1];
  s[1] = nu*D*e[0]+D*e[1];
  s[2] = 0.5*(1.0-nu)*D*e[2];
}
void TMRQuadStiffness::calculateThermalStress(const double pt[],
                                              const TacsScalar e[], 
                                              TacsScalar s[]){
  // Compute the penalized stiffness
  TacsScalar w = xw*1.0;
  TacsScalar D = E/(1.0-nu*nu)*w;

  // Compute the resulting matrix-vector product D*e
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
                                        TacsScalar fdvSens[],
                                        int dvLen ){
  
  TacsScalar dxw = 1.0+q*(1.0-xw);
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

// Evaluate the failure criteria
void TMRQuadStiffness::failure( const double pt[], 
                                const TacsScalar e[],
                                TacsScalar * fail ){
  // Use the von Mises failure criterion
  // Compute the relaxation factor
  TacsScalar r_factor = 1.0;
  if (eps > 0.0){
    r_factor = xw/(eps*(1.0-xw)+xw);
  }
  else {
    TacsScalar p = -1.0*eps;
    r_factor = pow(xw,p);
  }
  TacsScalar s[3];
  TacsScalar D = E/(1.0-nu*nu);
  // Compute the resulting matrix-vector product D*e
  s[0] = D*e[0]+nu*D*e[1];
  s[1] = nu*D*e[0]+D*e[1];
  s[2] = 0.5*(1.0-nu)*D*e[2];
 
  *fail = r_factor*VonMisesFailurePlaneStress(s,ys);
}
// Evaluate the failure criteria w.r.t. design variables
void TMRQuadStiffness::addFailureDVSens( const double pt[], 
                                         const TacsScalar e[],
                                         TacsScalar alpha,
                                         TacsScalar dvSens[], 
                                         int dvLen ){
  // Compute the relaxation factor
  TacsScalar r_factor_sens = 0.0;
  if (eps > 0.0){
    TacsScalar d = 1.0/(eps*(1.0-xw)+xw);
    r_factor_sens = eps*d*d;
  }
  else {
    TacsScalar p = -1.0*eps;
    r_factor_sens = pow(xw,p-1)*p;
  }
  TacsScalar s[3];
  TacsScalar D = E/(1.0-nu*nu);
  
  // Compute the resulting stress
  s[0] = D*e[0]+nu*D*e[1];
  s[1] = nu*D*e[0]+D*e[1];
  s[2] = 0.5*(1.0-nu)*D*e[2];
  
  TacsScalar fail = VonMisesFailurePlaneStress(s,ys);
  TacsScalar inner = alpha*r_factor_sens*fail;
  for ( int i = 0; i < nweights; i++ ){
    dvSens[weights[i].index] += weights[i].weight*inner;
  } 
}
void TMRQuadStiffness::failureStrainSens( const double pt[], 
                                          const TacsScalar e[],
                                          TacsScalar sens[]){
 
  TacsScalar s[3], ps_sens[3];
  // Use the von Mises failure criterion
  // Compute the relaxation factor
  TacsScalar r_factor = 1.0;
  if (eps > 0.0){
    r_factor = xw/(eps*(1.0-xw)+xw);
  }
  else {
    TacsScalar p = -1.0*eps;
    r_factor = pow(xw,p);
    //r_factor = sqrt(xw);
    
  }
  TacsScalar D = E/(1.0-nu*nu);
  // Compute the resulting stress
  s[0] = D*e[0]+nu*D*e[1];
  s[1] = nu*D*e[0]+D*e[1];
  s[2] = 0.5*(1.0-nu)*D*e[2];
  
  TacsScalar fail = VonMisesFailurePlaneStressSens(ps_sens, s, ys);
  sens[0] = r_factor*D*(ps_sens[0]+nu*ps_sens[1]);
  sens[1] = r_factor*D*(nu*ps_sens[0]+ps_sens[1]);
  sens[2] = r_factor*0.5*(1.0-nu)*D*ps_sens[2];
  
}

