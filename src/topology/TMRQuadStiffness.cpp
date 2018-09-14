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

#include "TMRQuadStiffness.h"

/*
  Create the quadtree stiffness object based on an interpolation from
  the filter variables
*/
TMRQuadStiffness::TMRQuadStiffness( TMRIndexWeight *_weights,
                                    int _nweights,
                                    TMRQuadStiffnessProperties *_props ){
  // Record the density, Poisson ratio, D and the shear modulus
  props = _props;
  props->incref();

  // Set the weights/local indices
  nweights = _nweights;
  weights = new TMRIndexWeight[ nweights ];
  memcpy(weights, _weights, nweights*sizeof(TMRIndexWeight));
  // Set the initial value for the densities
  nvars = 1;
  x[0] = 0.95;
  rho[0] = 0.95;
  if (props->nmats > 1){
    nvars = props->nmats+1;
    x[0] = 1.0;
    rho[0] = 1.0;
    for ( int j = 1; j < nvars; j++ ){
      x[j] = 1.0/(nvars-1);
      rho[j] = 1.0/(nvars-1);
    }
  }
}

/*
  Decref the props
*/
TMRQuadStiffness::~TMRQuadStiffness(){
  props->decref();
  delete [] weights;
}

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
void TMRQuadStiffness::setDesignVars( const TacsScalar xdv[], int numDVs ){
  const double beta = props->beta;
  const double xoffset = props->xoffset;
  const int use_project = props->use_project;

  for ( int j = 0; j < nvars; j++ ){
    x[j] = 0.0;
    for ( int i = 0; i < nweights; i++ ){
      x[j] += weights[i].weight*xdv[nvars*weights[i].index + j];
    }

    // Apply the projection to obtain the projected value
    // of the density
    if (use_project){
      rho[j] = 1.0/(1.0 + exp(-beta*(x[j] - xoffset)));
    }
    else{
      rho[j] = x[j];
    }
  }
}

/*
  Get the design variable values

  This is not possible to back out once the density is computed, so
  instead we use a constant value of 0.95 as the output.
*/
void TMRQuadStiffness::getDesignVars( TacsScalar xdv[], int numDVs ){
  if (nvars == 1){
    for ( int i = 0; i < nweights; i++ ){
      if (weights[i].index >= 0 && weights[i].index < numDVs){
        xdv[weights[i].index] = 0.5;
      }
    }
  }
  else {
    for ( int j = 0; j < nvars; j++ ){
      TacsScalar value = 0.5;
      if (j >= 1){
        value = 0.5/(nvars-1);
      }
      for ( int i = 0; i < nweights; i++ ){
        if (weights[i].index >= 0 && weights[i].index < numDVs){
          xdv[nvars*weights[i].index + j] = value;
        }
      }
    }
  }
}

/*
  Retrieve the range of the design variable values
*/
void TMRQuadStiffness::getDesignVarRange( TacsScalar lb[], 
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
void TMRQuadStiffness::calculateStress( const double pt[],
                                        const TacsScalar e[], 
                                        TacsScalar s[] ){
  const double k0 = props->k0;
  const double q = props->q;
  if (nvars == 1){
    // Compute the penalized stiffness
    TacsScalar penalty = rho[0]/(1.0 + q*(1.0 - rho[0]));

    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];

    TacsScalar Dp = (penalty + k0)*D;
    TacsScalar Gp = (penalty + k0)*G;
    
    s[0] = Dp*e[0]+nu*Dp*e[1];
    s[1] = nu*Dp*e[0]+Dp*e[1];
    s[2] = Gp*e[2];
  }
  else {
    // Compute the penalized stiffness
    s[0] = s[1] = s[2] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      // Compute the penalty
      TacsScalar penalty = rho[j]/(1.0 + q*(1.0 - rho[j]));

      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];

      // Add the penalized value
      TacsScalar Dp = (penalty + k0)*D;
      TacsScalar Gp = (penalty + k0)*G;

      s[0] += Dp*e[0]+nu*Dp*e[1];
      s[1] += nu*Dp*e[0]+Dp*e[1];
      s[2] += Gp*e[2];
    }
  }  
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
  const double q = props->q;
  const double beta = props->beta;
  const double xoffset = props->xoffset;
  const int use_project = props->use_project;

  if (nvars == 1){
    // Compute the derivative of the penalization with respect to
    // the projected density
    TacsScalar penalty =
      (q + 1.0)/((1.0 + q*(1.0 - rho[0]))*(1.0 + q*(1.0 - rho[0])));

    // Add the derivative of the projection
    if (use_project){
      penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
    }

    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];

    TacsScalar Dp = alpha*penalty*D;
    TacsScalar Gp = alpha*penalty*G;
    TacsScalar s[3];
    s[0] = Dp*e[0]+nu*Dp*e[1];
    s[1] = nu*Dp*e[0]+Dp*e[1];
    s[2] = Gp*e[2];
    
    TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
    for ( int i = 0; i < nweights; i++ ){
      fdvSens[weights[i].index] += weights[i].weight*product;
    }
  }
  else {
    for ( int j = 1; j < nvars; j++ ){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (q + 1.0)/((1.0 + q*(1.0 - rho[j]))*(1.0 + q*(1.0 - rho[j])));

      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[j] - xoffset))*rho[j]*rho[j];
      }

      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];

      // Add the result to the derivative
      TacsScalar Dp = alpha*penalty*D;
      TacsScalar Gp = alpha*penalty*G;
      TacsScalar s[3];
      s[0] = Dp*e[0]+nu*Dp*e[1];
      s[1] = nu*Dp*e[0]+Dp*e[1];
      s[2] = Gp*e[2];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[nvars*weights[i].index + j] += weights[i].weight*product;
      }
    }
  }
}

/*
  Compute the density of the material as a linear function of the
  element-wise density
*/
void TMRQuadStiffness::getPointwiseMass( const double pt[], 
                                         TacsScalar mass[] ){
  if (nvars == 1){
    mass[0] = x[0]*props->density[0];
  }
  else {
    mass[0] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      mass[0] += x[j]*props->density[j-1];
    }
  }  
}

/*
  Compute the derivative of the pointwise mass as function of the design
  variable values
*/
void TMRQuadStiffness::addPointwiseMassDVSens( const double pt[], 
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
