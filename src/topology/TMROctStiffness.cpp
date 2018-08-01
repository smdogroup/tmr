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
                                  TMRStiffnessProperties *_props ){
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
TMROctStiffness::~TMROctStiffness(){
  props->decref();
  delete [] weights;
}

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
void TMROctStiffness::setDesignVars( const TacsScalar xdv[], int numDVs ){
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
void TMROctStiffness::getDesignVars( TacsScalar xdv[], int numDVs ){
  if (nvars == 1){
    for ( int i = 0; i < nweights; i++ ){
      if (weights[i].index >= 0 && weights[i].index < numDVs){
        xdv[weights[i].index] = 0.95;
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
      TacsScalar Dp = (penalty + k0)*D;
      TacsScalar Gp = (penalty + k0)*G;
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
    TacsScalar s[6];
    s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
    s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
    s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
    s[3] = Gp*e[3];
    s[4] = Gp*e[4];
    s[5] = Gp*e[5];

    TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                          s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
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
      TacsScalar s[6];
      s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
      s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
      s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
      s[3] = Gp*e[3];
      s[4] = Gp*e[4];
      s[5] = Gp*e[5];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                            s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
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
void TMROctStiffness::getPointwiseMass( const double pt[],
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

// Evaluate the failure criteria
void TMROctStiffness::failure( const double pt[],
                               const TacsScalar e[],
                               TacsScalar *fail ){
  const double eps = props->eps;

  if (nvars == 1){
    // Use the von Mises failure criterion
    // Compute the relaxation factor
    TacsScalar r_factor = 1.0;
    TacsScalar xw = rho[0];

    if (eps > 0.0){
      r_factor = xw/(eps*(1.0 - xw) + xw);
    }
    else {
      TacsScalar p = -1.0*eps;
      r_factor = pow(xw, p);
    }

    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];
    TacsScalar ys = props->ys[0];
    TacsScalar s[6];

    // Compute the resulting stress
    s[0] = D*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
    s[1] = D*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
    s[2] = D*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
    s[3] = G*e[3];
    s[4] = G*e[4];
    s[5] = G*e[5];

    *fail = r_factor*VonMisesFailure3D(s, ys);
  }
}

/*
  Evaluate the derivative of the failure criteria w.r.t. design
  variables
*/
void TMROctStiffness::addFailureDVSens( const double pt[],
                                        const TacsScalar e[],
                                        TacsScalar alpha,
                                        TacsScalar dvSens[],
                                        int dvLen ){
  const double eps = props->eps;
  const double beta = props->beta;
  const double xoffset = props->xoffset;
  const int use_project = props->use_project;

  if (nvars == 1){
    // Compute the relaxation factor
    TacsScalar r_factor_sens = 0.0;
    TacsScalar xw = rho[0];

    if (eps > 0.0){
      TacsScalar d = 1.0/(eps*(1.0-xw) + xw);
      r_factor_sens = eps*d*d;
    }
    else {
      TacsScalar p = -eps;
      r_factor_sens = pow(xw, p-1.0)*p;
    }

    // Add the contribution from the derivative of the projection
    if (use_project){
      r_factor_sens *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
    }

    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];
    TacsScalar ys = props->ys[0];
    TacsScalar s[6];

    // Compute the resulting stress
    s[0] = D*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
    s[1] = D*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
    s[2] = D*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
    s[3] = G*e[3];
    s[4] = G*e[4];
    s[5] = G*e[5];

    TacsScalar fail = VonMisesFailure3D(s, ys);
    TacsScalar inner = alpha*r_factor_sens*fail;
    for ( int i = 0; i < nweights; i++ ){
      dvSens[weights[i].index] += weights[i].weight*inner;
    }
  }
  else {
    printf("Not implemented yet \n");
  }
}

void TMROctStiffness::failureStrainSens( const double pt[],
                                         const TacsScalar e[],
                                         TacsScalar sens[]){
  const double eps = props->eps;
  if (nvars == 1){
    TacsScalar s[6], ps_sens[6];
    // Use the von Mises failure criterion
    // Compute the relaxation factor
    TacsScalar r_factor = 1.0;
    TacsScalar xw = 1.0*rho[0];
    if (eps > 0.0){
      r_factor = xw/(eps*(1.0 - xw) + xw);
    }
    else {
      TacsScalar p = -1.0*eps;
      r_factor = pow(xw, p);
    }

    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];
    TacsScalar ys = props->ys[0];

    // Compute the resulting stress
    s[0] = D*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
    s[1] = D*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
    s[2] = D*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
    s[3] = G*e[3];
    s[4] = G*e[4];
    s[5] = G*e[5];
    VonMisesFailure3DStressSens(ps_sens, s, ys);
    sens[0] = r_factor*D*((1.0-nu)*ps_sens[0] + nu*(ps_sens[1] +
                                                    ps_sens[2]));
    sens[1] = r_factor*D*((1.0-nu)*ps_sens[1] + nu*(ps_sens[0] +
                                                    ps_sens[2]));
    sens[2] = r_factor*D*((1.0-nu)*ps_sens[2] + nu*(ps_sens[0] +
                                                    ps_sens[1]));

    sens[3] = r_factor*G*ps_sens[3];
    sens[4] = r_factor*G*ps_sens[4];
    sens[5] = r_factor*G*ps_sens[5];
  }
  else {
    printf("Not implemented yet \n");
  }
}

/*
  Create the octree stiffness object based on an interpolation from
  the filter variables
*/
TMRAnisotropicStiffness::TMRAnisotropicStiffness( TMRIndexWeight *_weights,
                                                  int _nweights,
                                                  TMRAnisotropicProperties
                                                    *_props ){
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
TMRAnisotropicStiffness::~TMRAnisotropicStiffness(){
  props->decref();
  delete [] weights;
}

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
void TMRAnisotropicStiffness::setDesignVars( const TacsScalar xdv[],
                                             int numDVs ){
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
void TMRAnisotropicStiffness::getDesignVars( TacsScalar xdv[],
                                             int numDVs ){
  if (nvars == 1){
    for ( int i = 0; i < nweights; i++ ){
      if (weights[i].index >= 0 && weights[i].index < numDVs){
        xdv[weights[i].index] = 0.95;
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
void TMRAnisotropicStiffness::getDesignVarRange( TacsScalar lb[],
                                                 TacsScalar ub[],
                                                 int numDVs ){
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
void TMRAnisotropicStiffness::calculateStress( const double pt[],
                                               const TacsScalar e[],
                                               TacsScalar s[] ){
  const double k0 = props->k0;
  const double q = props->q;

  if (nvars == 1){
    // Compute the penalized stiffness
    TacsScalar penalty = rho[0]/(1.0 + q*(1.0 - rho[0])) + k0;

    // Compute the stress
    s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = 0.0;
    addStress(penalty, props->C, e, s);
  }
  else {
    // Compute the penalized stiffness
    s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      // Compute the penalty
      TacsScalar penalty = rho[j]/(1.0 + q*(1.0 - rho[j])) + k0;

      addStress(penalty, &props->C[21*(j-1)], e, s);
    }
  }
}

/*
  Add the derivative of the product of the stress with the vector psi
  times alpha to the design variable array dvSens
*/
void TMRAnisotropicStiffness::addStressDVSens( const double pt[],
                                               const TacsScalar e[],
                                               TacsScalar alpha,
                                               const TacsScalar psi[],
                                               TacsScalar fdvSens[],
                                               int dvLen ){
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

    // Compute the product
    TacsScalar product = penalty*evalStressProduct(props->C, e, psi);
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

      // Compute the product
      TacsScalar product =
        penalty*evalStressProduct(&props->C[21*(j-1)], e, psi);
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
void TMRAnisotropicStiffness::getPointwiseMass( const double pt[],
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
void TMRAnisotropicStiffness::addPointwiseMassDVSens( const double pt[],
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

// Evaluate the failure criteria
void TMRAnisotropicStiffness::failure( const double pt[],
                                       const TacsScalar e[],
                                       TacsScalar *fail ){}

/*
  Evaluate the derivative of the failure criteria w.r.t. design
  variables
*/
void TMRAnisotropicStiffness::addFailureDVSens( const double pt[],
                                                const TacsScalar e[],
                                                TacsScalar alpha,
                                                TacsScalar dvSens[],
                                                int dvLen ){}

void TMRAnisotropicStiffness::failureStrainSens( const double pt[],
                                                 const TacsScalar e[],
                                                 TacsScalar sens[]){}
