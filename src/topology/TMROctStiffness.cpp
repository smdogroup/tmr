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

#include "TMROctStiffness.h"

/*
  Create the octree stiffness object based on an interpolation from
  the filter variables
*/
TMROctStiffness::TMROctStiffness( TMRIndexWeight *_weights,
                                  int _nweights,
                                  TMRStiffnessProperties *_props,
                                  double _q, double _eps ){
  // Record the density, Poisson ratio, D and the shear modulus
  props = _props;
  props->incref();
  q = _q;
  eps = _eps;

  // Set the weights/local indices
  nweights = _nweights;
  weights = new TMRIndexWeight[ nweights ];
  memcpy(weights, _weights, nweights*sizeof(TMRIndexWeight));
  
  // Set the initial value for the densities
  nvars = 1;
  rho[0] = 0.95;
  if (props->nmats > 1){
    nvars = props->nmats+1;
    rho[0] = 1.0;
    for ( int j = 1; j < nvars; j++ ){
      rho[j] = 1.0/(nvars-1);
    }
  }  
}

/*
  Decref the props
*/
TMROctStiffness::~TMROctStiffness(){
  props->decref();
  delete [] weights;
}

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
void TMROctStiffness::setDesignVars( const TacsScalar x[], int numDVs ){
  for ( int j = 0; j < nvars; j++ ){
    rho[j] = 0.0;
    for ( int i = 0; i < nweights; i++ ){
      rho[j] += weights[i].weight*x[nvars*weights[i].index + j];
    }
  }
}

/*
  Get the design variable values

  This is not possible to back out once the density is computed, so
  instead we use a constant value of 0.95 as the output.
*/
void TMROctStiffness::getDesignVars( TacsScalar x[], int numDVs ){
  if (nvars == 1){
    for ( int i = 0; i < nweights; i++ ){
      if (weights[i].index >= 0 && weights[i].index < numDVs){
        x[weights[i].index] = 0.95;
      }
    }
  }
  else {
    for ( int j = 0; j < nvars; j++ ){
      TacsScalar value = 0.5;
      if (j > 1){
        value = 0.5/(nvars-1);
      }
      for ( int i = 0; i < nweights; i++ ){
        if (weights[i].index >= 0 && weights[i].index < numDVs){
          x[nvars*weights[i].index + j] = value;
        }
      }
    }
  }
}

/*
  Retrieve the range of the design variable values
*/
void TMROctStiffness::getDesignVarRange( TacsScalar lb[], 
                                         TacsScalar ub[], int numDVs ){
  if (nvars == 1){
    for ( int i = 0; i < nweights; i++ ){
      if (weights[i].index >= 0 && weights[i].index < numDVs){
        lb[weights[i].index] = 0.0;
        ub[weights[i].index] = 1.0;
      }
    }
  }
  else {
    for ( int j = 0; j < nvars; j++ ){
      for ( int i = 0; i < nweights; i++ ){
        if (weights[i].index >= 0 && weights[i].index < numDVs){
          lb[nvars*weights[i].index + j] = 0.0;
          ub[nvars*weights[i].index + j] = 1e30;
        }
      }
    }
  }
}

/*
  Given the values of everything, compute the stress
*/
void TMROctStiffness::calculateStress( const double pt[],
                                       const TacsScalar e[], 
                                       TacsScalar s[] ){
  if (nvars == 1){
    // Compute the penalized stiffness
    TacsScalar penalty = rho[0]/(1.0 + q*(1.0 - rho[0]));

    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];

    TacsScalar Dp = (penalty + eps)*D;
    TacsScalar Gp = (penalty + eps)*G;
    s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
    s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
    s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
    s[3] = Gp*e[3];
    s[4] = Gp*e[4];
    s[5] = Gp*e[5];
  }
  else {
    // Compute the penalized stiffness
    s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      // Compute the penalty
      TacsScalar penalty = rho[j]/(1.0 + q*(1.0 - rho[j]));

      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];

      // Add the penalized value
      TacsScalar Dp = (penalty + eps)*D;
      TacsScalar Gp = (penalty + eps)*G;
      s[0] += Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
      s[1] += Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
      s[2] += Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
      s[3] += Gp*e[3];
      s[4] += Gp*e[4];
      s[5] += Gp*e[5];
    }
  }
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
  if (nvars == 1){
    // Compute the derivative of the penalty
    TacsScalar penalty = 
      (q + 1.0)/((1.0 + q*(1.0 - rho[0]))*(1.0 + q*(1.0 - rho[0])));

    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];

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
  else {
    for ( int j = 1; j < nvars; j++ ){
      // Compute the penalty
      TacsScalar penalty = 
        (q + 1.0)/((1.0 + q*(1.0 - rho[j]))*(1.0 + q*(1.0 - rho[j])));

      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];

      // Add the result to the derivative
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
        fdvSens[nvars*weights[i].index + j] += 
          weights[i].weight*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] + 
                             s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
      }
    }
  }
}

/*
  Compute the density of the material as a linear function of the
  element-wise density
*/
void TMROctStiffness::getPointwiseMass( const double pt[], 
                                        TacsScalar mass[] ){
  if (nvars == 1){
    mass[0] = rho[0]*props->density[0];
  }
  else {
    mass[0] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      mass[0] += rho[j]*props->density[j-1];
    }
  }
}

/*
  Compute the derivative of the pointwise mass as function of the design
  variable values
*/
void TMROctStiffness::addPointwiseMassDVSens( const double pt[], 
                                              const TacsScalar alpha[],
                                              TacsScalar dvSens[], 
                                              int dvLen ){
  if (nvars == 1){
    TacsScalar scale = props->density[0]*alpha[0];
    for ( int i = 0; i < nweights; i++ ){
      dvSens[weights[i].index] += scale*weights[i].weight;
    }
  }
  else {
    for ( int j = 1; j < nvars; j++ ){
      TacsScalar scale = props->density[j-1]*alpha[0];
      for ( int i = 0; i < nweights; i++ ){
        dvSens[nvars*weights[i].index + j] += scale*weights[i].weight;
      }
    }
  }
}

/*
  Create the octree stiffness object based on an interpolation from
  the filter variables
*/
TMRLinearOctStiffness::TMRLinearOctStiffness( TMRIndexWeight *_weights,
                                              int _nweights,
                                              TacsScalar _x_init,
                                              TacsScalar _density,
                                              TacsScalar E, 
                                              TacsScalar _nu,
                                              double _q,
                                              PenaltyType _type,
                                              double _eps ){
  // Record the density, Poisson ratio, D and the shear modulus
  density = _density;
  nu = _nu;
  D = E/((1.0 + nu)*(1.0 - 2.0*nu));
  G = 0.5*E/(1.0 + nu);
  q = _q;
  penalty_type = _type;
  eps = _eps;

  // Set the initial values of the linearization
  rho_const = eps;
  rho_linear = 1.0;

  // Set the weights/local indices
  nweights = _nweights;
  memcpy(weights, _weights, nweights*sizeof(TMRIndexWeight));

  // Set the lower bound for the design variables
  for ( int i = 0; i < nweights; i++ ){
    x_vals[i] = _x_init;
    x_lb[i] = 0.0;
  }

  // Set the initial value for the density
  rho = 0.95;
}

// Set the default penalty type
TMRLinearOctStiffness::PenaltyType 
TMRLinearOctStiffness::penalty_type = TMRLinearOctStiffness::SIMP;

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
void TMRLinearOctStiffness::setDesignVars( const TacsScalar x[], 
                                           int numDVs ){
  rho = 0.0;
  for ( int i = 0; i < nweights; i++ ){
    x_vals[i] = x[weights[i].index];
    rho += weights[i].weight*x[weights[i].index];
  }
}

/*
  Get the design variable values

  This is not possible to back out once the density is computed, so
  instead we use a constant value of 0.95 as the output.
*/
void TMRLinearOctStiffness::getDesignVars( TacsScalar x[], int numDVs ){
  for ( int i = 0; i < nweights; i++ ){
    if (weights[i].index >= 0 && weights[i].index < numDVs){
      x[weights[i].index] = x_vals[i];
    }
  }
}

/*
  Retrieve the range of the design variable values
*/
void TMRLinearOctStiffness::getDesignVarRange( TacsScalar lb[], 
                                               TacsScalar ub[], 
                                               int numDVs ){
  for ( int i = 0; i < nweights; i++ ){
    if (weights[i].index >= 0 && weights[i].index < numDVs){
      lb[weights[i].index] = x_lb[i];
      ub[weights[i].index] = 1.0;
    }
  }
}

/*
  Set the linearization
*/
void TMRLinearOctStiffness::setLinearization( TacsScalar _q,
                                              const TacsScalar x[], 
                                              int numDVs ){ 
  // Set the value of the RAMP penalization
  q = _q;

  // Compute the value of rho
  rho = 0.0;
  for ( int i = 0; i < nweights; i++ ){
    x_vals[i] = x[weights[i].index];
    rho += weights[i].weight*x[weights[i].index];
  }
 
  if (penalty_type == SIMP){
    // Set the linearization for the material weights
    if (q > 1.0){
      rho_const = pow(rho, q) + eps;
      rho_linear = q*pow(rho, q-1.0);
    
      // Set the lower bound
      for ( int i = 0; i < nweights; i++ ){
        x_lb[i] = ((q - 1.0)/q)*x_vals[i];
      }
    }
    else {
      rho_const = eps;
      rho_linear = 1.0;

      // Set the lower bound
      for ( int i = 0; i < nweights; i++ ){
        x_lb[i] = 0.0;
      }
    }
  }
  else {
    if (q > 0.0){
      // Set the linearization for the material weights
      rho_const = rho/(1.0 + q*(1.0 - rho)) + eps;
      rho_linear = (q + 1.0)/((1.0 + q*(1.0 - rho))*(1.0 + q*(1.0 - rho)));
    
      // Set the lower bound
      for ( int i = 0; i < nweights; i++ ){
        x_lb[i] = (q/(1.0 + q))*x_vals[i];
      }
    }
    else {
      rho_const = eps;
      rho_linear = 1.0;

      // Set the lower bound
      for ( int i = 0; i < nweights; i++ ){
        x_lb[i] = 0.0;
      }
    }
  }
}

/*
  Given the values of everything, compute the stress
*/
void TMRLinearOctStiffness::calculateStress( const double pt[],
                                             const TacsScalar e[], 
                                             TacsScalar s[] ){
  // Compute the penalized stiffness
  TacsScalar Dp = (rho_const + rho_linear*rho)*D;
  TacsScalar Gp = (rho_const + rho_linear*rho)*G;
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
void TMRLinearOctStiffness::addStressDVSens( const double pt[], 
                                             const TacsScalar e[], 
                                             TacsScalar alpha, 
                                             const TacsScalar psi[], 
                                             TacsScalar fdvSens[], 
                                             int dvLen ){
  TacsScalar Dp = alpha*rho_linear*D;
  TacsScalar Gp = alpha*rho_linear*G;
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
void TMRLinearOctStiffness::getPointwiseMass( const double pt[], 
                                              TacsScalar mass[] ){
  mass[0] = rho*density;
}

/*
  Compute the derivative of the pointwise mass as function of the design
  variable values
*/
void TMRLinearOctStiffness::addPointwiseMassDVSens( const double pt[], 
                                                    const TacsScalar alpha[],
                                                    TacsScalar dvSens[], 
                                                    int dvLen ){
  TacsScalar scale = density*alpha[0];
  for ( int i = 0; i < nweights; i++ ){
    dvSens[weights[i].index] += scale*weights[i].weight;
  }
}

/*
  Compute the derivative of the stress projected onto the provided
  design variable vector.
*/
void TMRLinearOctStiffness::calcStressDVProject( const double pt[],
                                                 const TacsScalar e[],
                                                 const TacsScalar px[],
                                                 int dvLen, 
                                                 TacsScalar s[] ){
  TacsScalar alpha = 0.0;
  for ( int i = 0; i < nweights; i++ ){
    alpha += px[i]*weights[i].weight;
  }

  // Compute the product of the stress derivative times px
  TacsScalar Dp = alpha*rho_linear*D;
  TacsScalar Gp = alpha*rho_linear*G;
  s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
  s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
  s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
  s[3] = Gp*e[3];
  s[4] = Gp*e[4];
  s[5] = Gp*e[5];
}
