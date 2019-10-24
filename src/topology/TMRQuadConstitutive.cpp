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

#include "TMRQuadConstitutive.h"

/*
  Create the octree stiffness object based on an interpolation from
  the filter variables
*/
TMRQuadConstitutive::TMRQuadConstitutive( TMRStiffnessProperties *_props,
                                          TMRQuadForest *_forest ):
  TACSPlaneStressConstitutive(NULL){
  // Record the density, Poisson ratio, D and the shear modulus
  props = _props;
  props->incref();

  forest = _forest;
  forest->incref();

  nmats = props->nmats;
  nvars = 1;
  if (nmats >= 2){
    nvars += nmats;
  }

  // Get the connectivity information
  int num_elements;
  const int order = forest->getMeshOrder();
  forest->getNodeConn(NULL, &num_elements);
  int nconn = order*order*num_elements;

  // Allocate space for the shape functions
  N = new double[ order*order ];
  temp_array = new TacsScalar[ 2*nmats ];

  // Initialize the design vector
  x = new TacsScalar[ nvars*nconn ];
  if (nvars == 1){
    for ( int i = 0; i < nconn; i++ ){
      x[i] = 0.95;
    }
  }
  else {
    for ( int i = 0; i < nvars*nconn; i++ ){
      if (i % nvars == 0){
        x[i] = 0.95;
      }
      else {
        x[i] = 0.95/nmats;
      }
    }
  }
}

/*
  Deallocate the constitutive object
*/
TMRQuadConstitutive::~TMRQuadConstitutive(){
  props->decref();
  forest->decref();
  delete [] x;
  delete [] N;
  delete [] temp_array;
}

/*
  Retrieve the design variable values
*/
int TMRQuadConstitutive::getDesignVarNums( int elemIndex,
                                           int dvLen,
                                           int dvNums[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;
  if (dvNums){
    const int *conn;
    forest->getNodeConn(&conn);
    conn = &conn[len*elemIndex];
    for ( int i = 0; (i < len && i < dvLen); i++ ){
      dvNums[i] = conn[i];
    }
  }
  return len;
}

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
int TMRQuadConstitutive::setDesignVars( int elemIndex,
                                        int dvLen,
                                        const TacsScalar dvs[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  TacsScalar *xptr = &x[nvars*len*elemIndex];
  for ( int i = 0; i < nvars*len; i++ ){
    xptr[i] = dvs[i];
  }

  return len;
}

/*
  Get the design variable values

  This is not possible to back out once the density is computed, so
  instead we use a constant value of 0.95 as the output.
*/
int TMRQuadConstitutive::getDesignVars( int elemIndex,
                                        int dvLen,
                                        TacsScalar dvs[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  TacsScalar *xptr = &x[nvars*len*elemIndex];
  for ( int i = 0; i < nvars*len; i++ ){
    dvs[i] = xptr[i];
  }

  return len;
}

/*
  Retrieve the range of the design variable values
*/
int TMRQuadConstitutive::getDesignVarRange( int elemIndex,
                                            int dvLen,
                                            TacsScalar lb[],
                                            TacsScalar ub[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  if (nmats == 1){
    for ( int i = 0; i < len; i++ ){
      lb[i] = 0.0;
      ub[i] = 1.0;
    }
  }
  else {
    for ( int i = 0; i < nvars*len; i++ ){
      lb[i] = 0.0;
      ub[i] = 1e30;
    }
  }

  return len;
}

/*
  Evaluate the material density: Linear combination of the density for
  each material times the design density
*/
TacsScalar TMRQuadConstitutive::evalDensity( int elemIndex,
                                             const double pt[],
                                             const TacsScalar X[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  // Evaluate the density
  TacsScalar density = 0.0;
  if (nvars == 1){
    TacsScalar rho = 0.0;
    for ( int i = 0; i < len; i++ ){
      rho += N[i]*xptr[i];
    }

    density = rho*props->props[0]->getDensity();
  }
  else {
    for ( int j = 0; j < nmats; j++ ){
      TacsScalar rho = 0.0;
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }

      density += rho*props->props[j]->getDensity();
    }
  }

  return density;
}

// Add the derivative of the density
void TMRQuadConstitutive::addDensityDVSens( int elemIndex,
                                            const double pt[],
                                            const TacsScalar X[],
                                            const TacsScalar scale,
                                            int dvLen,
                                            TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Add the derivative of the density
  if (nvars == 1){
    TacsScalar density = props->props[0]->getDensity();
    for ( int i = 0; i < len; i++ ){
      dfdx[i] += scale*N[i]*density;
    }
  }
  else {
    for ( int j = 0; j < nmats; j++ ){
      TacsScalar density = props->props[j]->getDensity();
      for ( int i = 0; i < len; i++ ){
        dfdx[nvars*i + j+1] += scale*N[i]*density;
      }
    }
  }
}

/*
  Evaluate the specific heat: Rule of mixtures for the specific heat
*/
TacsScalar TMRQuadConstitutive::evalSpecificHeat( int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  // Evaluate the density
  TacsScalar specific_heat = 0.0;
  if (nvars == 1){
    TacsScalar rho = 0.0;
    for ( int i = 0; i < len; i++ ){
      rho += N[i]*xptr[i];
    }

    specific_heat = rho*props->props[0]->getSpecificHeat();
  }
  else {
    for ( int j = 0; j < nmats; j++ ){
      TacsScalar rho = 0.0;
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }

      specific_heat += rho*props->props[j]->getSpecificHeat();
    }
  }

  return specific_heat;
}

/*
  Add the derivative of the specific heat
*/
void TMRQuadConstitutive::addSpecificHeatDVSens( int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 const TacsScalar scale,
                                                 int dvLen,
                                                 TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Add the derivative of the density
  if (nvars == 1){
    TacsScalar specific_heat = props->props[0]->getSpecificHeat();
    for ( int i = 0; i < len; i++ ){
      dfdx[i] += scale*N[i]*specific_heat;
    }
  }
  else {
    for ( int j = 0; j < nmats; j++ ){
      TacsScalar specific_heat = props->props[j]->getSpecificHeat();
      for ( int i = 0; i < len; i++ ){
        dfdx[nvars*i + j+1] += scale*N[i]*specific_heat;
      }
    }
  }
}

/*
  Evaluate the stresss
*/
void TMRQuadConstitutive::evalStress( int elemIndex,
                                      const double pt[],
                                      const TacsScalar X[],
                                      const TacsScalar e[],
                                      TacsScalar s[] ){
  TacsScalar C[6];
  evalTangentStiffness(elemIndex, pt, X, C);
  s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
  s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
  s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);
}

/*
  Evaluate the tangent stiffness matrix
*/
void TMRQuadConstitutive::evalTangentStiffness( int elemIndex,
                                                const double pt[],
                                                const TacsScalar X[],
                                                TacsScalar C[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;
  const double q = props->q;
  const double beta = props->beta;
  const double xoffset = props->xoffset;

  memset(C, 0, 6*sizeof(TacsScalar));

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  // Add the derivative of the density
  for ( int j = 0; j < nmats; j++ ){
    TacsScalar rho = 0.0;
    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[i];
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }
    }

    // Use projection
    if (props->use_project){
      rho = 1.0/(1.0 + exp(-beta*(rho - xoffset)));
    }

    // Evaluate the tangent stiffness
    TacsScalar Cj[6];
    props->props[j]->evalTangentStiffness2D(Cj);

    // Compute the penalty
    TacsScalar penalty = rho/(1.0 + q*(1.0 - rho));

    // Add up the contribution to the tangent stiffness matrix
    for ( int i = 0; i < 6; i++ ){
      C[i] += penalty*Cj[i];
    }
  }
}

/*
  Add the derivative of the product of the stress with the vector psi
  times alpha to the design variable array dvSens
*/
void TMRQuadConstitutive::addStressDVSens( int elemIndex,
                                           const double pt[],
                                           const TacsScalar X[],
                                           const TacsScalar e[],
                                           TacsScalar scale,
                                           const TacsScalar psi[],
                                           int dvLen,
                                           TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;
  const double q = props->q;
  const double beta = props->beta;
  const double xoffset = props->xoffset;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  // Add the derivative of the density
  for ( int j = 0; j < nmats; j++ ){
    TacsScalar rho = 0.0;
    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[i];
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }
    }

    // Compute the derivative of the penalization with respect to
    // the projected density
    TacsScalar dpenalty =
      (q + 1.0)/((1.0 + q*(1.0 - rho))*(1.0 + q*(1.0 - rho)));

    // Add the derivative of the projection
    if (props->use_project){
      TacsScalar rho_proj = 1.0/(1.0 + exp(-beta*(rho - xoffset)));
      dpenalty *= beta*exp(-beta*(rho - xoffset))*rho_proj*rho_proj;
    }

    TacsScalar C[6], s[3];
    props->props[j]->evalTangentStiffness2D(C);
    s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    TacsScalar product = scale*dpenalty*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);

    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        dfdx[i] += N[i]*product;
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        dfdx[nvars*i + j+1] += N[i]*product;
      }
    }
  }
}

// Evaluate the thermal strain
void TMRQuadConstitutive::evalThermalStrain( int elemIndex,
                                             const double pt[],
                                             const TacsScalar X[],
                                             TacsScalar theta,
                                             TacsScalar e[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  // Add the derivative of the density
  for ( int j = 0; j < nmats; j++ ){
    TacsScalar rho = 0.0;
    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[i];
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }
    }

    TacsScalar et[3];
    props->props[j]->evalThermalStrain2D(et);
    e[0] += rho*theta*et[0];
    e[1] += rho*theta*et[1];
    e[2] += rho*theta*et[2];
  }
}

void TMRQuadConstitutive::addThermalStrainDVSens( int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  TacsScalar theta,
                                                  const TacsScalar psi[],
                                                  int dvLen,
                                                  TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  // Add the derivative of the density
  for ( int j = 0; j < nmats; j++ ){
    TacsScalar rho = 0.0;
    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[i];
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }
    }

    TacsScalar et[3];
    props->props[j]->evalThermalStrain2D(et);

    TacsScalar product = theta*(et[0]*psi[0] + et[1]*psi[1] + et[2]*psi[2]);

    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        dfdx[i] += product*N[i];
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        dfdx[nvars*i + j+1] += product*N[i];
      }
    }
  }
}

// Evaluate the heat flux, given the thermal gradient
void TMRQuadConstitutive::evalHeatFlux( int elemIndex,
                                        const double pt[],
                                        const TacsScalar X[],
                                        const TacsScalar grad[],
                                        TacsScalar flux[] ){
  TacsScalar C[3];
  evalTangentHeatFlux(elemIndex, pt, X, C);
  flux[0] = C[0]*grad[0] + C[1]*grad[1];
  flux[1] = C[1]*grad[0] + C[2]*grad[1];
}

// Evaluate the tangent of the heat flux
void TMRQuadConstitutive::evalTangentHeatFlux( int elemIndex,
                                               const double pt[],
                                               const TacsScalar X[],
                                               TacsScalar C[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;
  const double qcond = props->qcond;
  const double beta = props->beta;
  const double xoffset = props->xoffset;

  memset(C, 0, 3*sizeof(TacsScalar));

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  // Add the derivative of the density
  for ( int j = 0; j < nmats; j++ ){
    TacsScalar rho = 0.0;
    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[i];
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }
    }

    // Use projection
    if (props->use_project){
      rho = 1.0/(1.0 + exp(-beta*(rho - xoffset)));
    }

    // Evaluate the tangent stiffness
    TacsScalar Cj[3];
    props->props[j]->evalTangentHeatFlux2D(Cj);

    // Compute the penalty
    TacsScalar penalty = rho/(1.0 + qcond*(1.0 - rho));

    // Add up the contribution to the tangent stiffness matrix
    for ( int i = 0; i < 3; i++ ){
      C[i] += penalty*Cj[i];
    }
  }
}

// Add the derivative of the heat flux
void TMRQuadConstitutive::addHeatFluxDVSens( int elemIndex,
                                             const double pt[],
                                             const TacsScalar X[],
                                             const TacsScalar grad[],
                                             TacsScalar scale,
                                             const TacsScalar psi[],
                                             int dvLen,
                                             TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;
  const double qcond = props->qcond;
  const double beta = props->beta;
  const double xoffset = props->xoffset;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  // Add the derivative of the density
  for ( int j = 0; j < nmats; j++ ){
    TacsScalar rho = 0.0;
    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[i];
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }
    }

    // Compute the derivative of the penalization with respect to
    // the projected density
    TacsScalar dpenalty =
      (qcond + 1.0)/((1.0 + qcond*(1.0 - rho))*(1.0 + qcond*(1.0 - rho)));

    // Add the derivative of the projection
    if (props->use_project){
      TacsScalar rho_proj = 1.0/(1.0 + exp(-beta*(rho - xoffset)));
      dpenalty *= beta*exp(-beta*(rho - xoffset))*rho_proj*rho_proj;
    }

    TacsScalar C[3], flux[2];
    props->props[j]->evalTangentHeatFlux2D(C);
    flux[0] = C[0]*grad[0] + C[1]*grad[1];
    flux[1] = C[1]*grad[0] + C[2]*grad[1];

    TacsScalar product = scale*dpenalty*(flux[0]*psi[0] +
                                         flux[1]*psi[1]);

    if (nvars == 1){
      for ( int i = 0; i < len; i++ ){
        dfdx[i] += N[i]*product;
      }
    }
    else {
      for ( int i = 0; i < len; i++ ){
        dfdx[nvars*i + j+1] += N[i]*product;
      }
    }
  }
}

// Evaluate the failure criteria
TacsScalar TMRQuadConstitutive::evalFailure( int elemIndex,
                                             const double pt[],
                                             const TacsScalar X[],
                                             const TacsScalar e[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;
  const double eps = props->eps;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  if (nvars == 1){
    TacsScalar C[6];
    props->props[0]->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    TacsScalar rho = 0.0;
    for ( int i = 0; i < len; i++ ){
      rho += N[i]*xptr[i];
    }

    TacsScalar r_factor = rho/(eps*(1.0 - rho) + rho);

    TacsScalar fail = props->props[0]->vonMisesFailure2D(s);

    return r_factor*fail;
  }
  else {
    const double ksWeight = props->ksWeight;
    TacsScalar *fail = temp_array;
    TacsScalar max_fail = -1e20;

    for ( int j = 0; j < nmats; j++ ){
      TacsScalar C[6];
      props->props[j]->evalTangentStiffness2D(C);

      TacsScalar s[3];
      s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
      s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
      s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

      TacsScalar rho = 0.0;
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }

      TacsScalar r_factor = rho/(eps*(1.0 - rho) + rho);

      fail[j] = r_factor*props->props[j]->vonMisesFailure2D(s);
      if (TacsRealPart(fail[j]) > max_fail){
        max_fail = fail[j];
      }
    }

    // Compute the KS aggregate of the failure
    TacsScalar ksSum = 0.0;
    for ( int j = 0; j < nmats; j++ ){
      ksSum += exp(ksWeight*(fail[j] - max_fail));
    }

    TacsScalar ksFail = max_fail + log(ksSum)/ksWeight;

    return ksFail;
  }

  return 0.0;
}

/*
  Evaluate the derivative of the failure criteria w.r.t. design
  variables
*/
void TMRQuadConstitutive::addFailureDVSens( int elemIndex,
                                            const double pt[],
                                            const TacsScalar X[],
                                            const TacsScalar e[],
                                            TacsScalar scale,
                                            int dvLen,
                                            TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;
  const double eps = props->eps;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  if (nvars == 1){
    TacsScalar C[6];
    props->props[0]->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    TacsScalar rho = 0.0;
    for ( int i = 0; i < len; i++ ){
      rho += N[i]*xptr[i];
    }

    // Compute the derivative of the stress relaxation factor
    TacsScalar d = 1.0/(eps*(1.0 - rho) + rho);
    TacsScalar r_factor_sens = eps*d*d;

    // Compute the contribution to the von Mises failure
    TacsScalar fail = props->props[0]->vonMisesFailure2D(s);

    // Add the derivative
    TacsScalar contrib = scale*r_factor_sens*fail;
    for ( int i = 0; i < len; i++ ){
      dfdx[i] += contrib*N[i];
    }
  }
  else {
    const double ksWeight = props->ksWeight;
    TacsScalar *fail = &temp_array[0];
    TacsScalar *rho = &temp_array[nmats];
    TacsScalar max_fail = -1e20;

    for ( int j = 0; j < nmats; j++ ){
      TacsScalar C[6];
      props->props[j]->evalTangentStiffness2D(C);

      TacsScalar s[3];
      s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
      s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
      s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

      // Compute the density for the j-th material
      rho[j] = 0.0;
      for ( int i = 0; i < len; i++ ){
        rho[j] += N[i]*xptr[nvars*i + j+1];
      }

      // Compute the stress relaxation factor
      TacsScalar r_factor = rho[j]/(eps*(1.0 - rho[j]) + rho[j]);

      // Compute the failure value
      fail[j] = r_factor*props->props[j]->vonMisesFailure2D(s);
      if (TacsRealPart(fail[j]) > max_fail){
        max_fail = fail[j];
      }
    }

    // Compute the KS aggregate of the failure
    TacsScalar ksSum = 0.0;
    for ( int j = 0; j < nmats; j++ ){
      fail[j] = exp(ksWeight*(fail[j] - max_fail));
      ksSum += fail[j];
    }

    for ( int j = 0; j < nmats; j++ ){
      TacsScalar C[6];
      props->props[j]->evalTangentStiffness2D(C);

      TacsScalar s[3];
      s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
      s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
      s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

      // Compute the derivative of the stress relaxation factor
      TacsScalar d = 1.0/(eps*(1.0 - rho[j]) + rho[j]);
      TacsScalar r_factor_sens = eps*d*d;

      // Compute the contribution to the von Mises failure
      TacsScalar fail_val = props->props[j]->vonMisesFailure2D(s);

      // Add the derivative
      TacsScalar contrib = scale*r_factor_sens*fail_val*(fail[j]/ksSum);
      for ( int i = 0; i < len; i++ ){
        dfdx[i*nvars + j+1] += contrib*N[i];
      }
    }
  }
}

TacsScalar TMRQuadConstitutive::evalFailureStrainSens( int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[],
                                                       const TacsScalar e[],
                                                       TacsScalar dfde[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order;
  const double eps = props->eps;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  if (nvars == 1){
    TacsScalar C[6];
    props->props[0]->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    TacsScalar t[3];
    TacsScalar fail = props->props[0]->vonMisesFailure2DStressSens(s, t);

    dfde[0] = (C[0]*t[0] + C[1]*t[1] + C[2]*t[2]);
    dfde[1] = (C[1]*t[0] + C[3]*t[1] + C[4]*t[2]);
    dfde[2] = (C[2]*t[0] + C[4]*t[1] + C[5]*t[2]);

    TacsScalar rho = 0.0;
    for ( int i = 0; i < len; i++ ){
      rho += N[i]*xptr[i];
    }

    TacsScalar r_factor = rho/(eps*(1.0 - rho) + rho);

    for ( int i = 0; i < 3; i++ ){
      dfde[i] *= r_factor;
    }

    return r_factor*fail;
  }
  else {
    const double ksWeight = props->ksWeight;
    TacsScalar *fail = &temp_array[0];
    TacsScalar *rho = &temp_array[nmats];
    TacsScalar max_fail = -1e20;

    for ( int j = 0; j < nmats; j++ ){
      TacsScalar C[6];
      props->props[j]->evalTangentStiffness2D(C);

      TacsScalar s[3];
      s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
      s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
      s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

      // Compute the density for the j-th material
      rho[j] = 0.0;
      for ( int i = 0; i < len; i++ ){
        rho[j] += N[i]*xptr[nvars*i + j+1];
      }

      // Compute the stress relaxation factor
      TacsScalar r_factor = rho[j]/(eps*(1.0 - rho[j]) + rho[j]);

      // Compute the failure value
      fail[j] = r_factor*props->props[j]->vonMisesFailure2D(s);
      if (TacsRealPart(fail[j]) > max_fail){
        max_fail = fail[j];
      }
    }

    // Compute the KS aggregate of the failure
    TacsScalar ksSum = 0.0;
    for ( int j = 0; j < nmats; j++ ){
      fail[j] = exp(ksWeight*(fail[j] - max_fail));
      ksSum += fail[j];
    }

    TacsScalar ksFail = max_fail + log(ksSum)/ksWeight;

    memset(dfde, 0, 3*sizeof(TacsScalar));
    for ( int j = 0; j < nmats; j++ ){
      TacsScalar C[6];
      props->props[j]->evalTangentStiffness2D(C);

      TacsScalar s[3];
      s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
      s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
      s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

      TacsScalar t[3];
      props->props[j]->vonMisesFailure2DStressSens(s, t);

      // Compute the stress relaxation factor
      TacsScalar r_factor = rho[j]/(eps*(1.0 - rho[j]) + rho[j]);

      TacsScalar scale = r_factor*fail[j]/ksSum;
      dfde[0] += scale*(C[0]*t[0] + C[1]*t[1] + C[2]*t[2]);
      dfde[1] += scale*(C[1]*t[0] + C[3]*t[1] + C[4]*t[2]);
      dfde[2] += scale*(C[2]*t[0] + C[4]*t[1] + C[5]*t[2]);
    }

    return ksFail;
  }

  return 0.0;
}

const char* TMRQuadConstitutive::getObjectName(){
  return "TMRQuadConstitutive";
}
