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

#include "TMROctConstitutive.h"

/*
  Store the stiffness properties for the constitutive classes
*/
TMRStiffnessProperties::TMRStiffnessProperties( int _nmats,
                                                TACSMaterialProperties **_props,
                                                double _q,
                                                double _eps,
                                                double _k0,
                                                double _ksWeight,
                                                double _qtemp,
                                                double _qcond,
                                                double _beta,
                                                double _xoffset,
                                                int _use_project ){
  nmats = _nmats;
  props = new TACSMaterialProperties*[ nmats ];
  for ( int i = 0; i < nmats; i++ ){
    props[i] = _props[i];
    props[i]->incref();
  }

  q = _q;
  k0 = _k0;
  eps = _eps;

  // Set the KS weight for the multimaterial problem
  ksWeight = _ksWeight;

  // Penalization for thermoelastic problem
  qtemp = _qtemp;
  qcond = _qcond;

  // Projection
  beta = _beta;
  xoffset = _xoffset;
  use_project = _use_project;
}

/*
  Create the octree stiffness object based on an interpolation from
  the filter variables
*/
TMROctConstitutive::TMROctConstitutive( TMRStiffnessProperties *_props,
                                        TMROctForest *_forest ):
  TACSSolidConstitutive(NULL){
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
  int nconn = order*order*order*num_elements;

  // Allocate space for the shape functions
  N = new double[ order*order*order ];
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
TMROctConstitutive::~TMROctConstitutive(){
  props->decref();
  forest->decref();
  delete [] x;
  delete [] N;
  delete [] temp_array;
}

/*
  Retrieve the design variable values
*/
int TMROctConstitutive::getDesignVarNums( int elemIndex,
                                          int dvLen,
                                          int dvNums[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;
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
int TMROctConstitutive::setDesignVars( int elemIndex,
                                       int dvLen,
                                       const TacsScalar dvs[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

  TacsScalar *xptr = &x[nmats*len*elemIndex];
  for ( int i = 0; i < nmats*len; i++ ){
    xptr[i] = dvs[i];
  }

  return len;
}

/*
  Get the design variable values

  This is not possible to back out once the density is computed, so
  instead we use a constant value of 0.95 as the output.
*/
int TMROctConstitutive::getDesignVars( int elemIndex,
                                       int dvLen,
                                       TacsScalar dvs[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

  TacsScalar *xptr = &x[nvars*len*elemIndex];
  for ( int i = 0; i < nmats*len; i++ ){
    dvs[i] = xptr[i];
  }

  return len;
}

/*
  Retrieve the range of the design variable values
*/
int TMROctConstitutive::getDesignVarRange( int elemIndex,
                                           int dvLen,
                                           TacsScalar lb[],
                                           TacsScalar ub[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

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
TacsScalar TMROctConstitutive::evalDensity( int elemIndex,
                                            const double pt[],
                                            const TacsScalar X[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

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
void TMROctConstitutive::addDensityDVSens( int elemIndex,
                                           const double pt[],
                                           const TacsScalar X[],
                                           const TacsScalar scale,
                                           int dvLen,
                                           TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

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
TacsScalar TMROctConstitutive::evalSpecificHeat( int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

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
void TMROctConstitutive::addSpecificHeatDVSens( int elemIndex,
                                                const double pt[],
                                                const TacsScalar X[],
                                                const TacsScalar scale,
                                                int dvLen,
                                                TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

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
void TMROctConstitutive::evalStress( int elemIndex,
                                     const double pt[],
                                     const TacsScalar X[],
                                     const TacsScalar e[],
                                     TacsScalar s[] ){
  TacsScalar C[21];
  evalTangentStiffness(elemIndex, pt, X, C);
  s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
  s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
  s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
  s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
  s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
  s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];
}

/*
  Evaluate the tangent stiffness matrix
*/
void TMROctConstitutive::evalTangentStiffness( int elemIndex,
                                               const double pt[],
                                               const TacsScalar X[],
                                               TacsScalar C[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;
  const double q = props->q;
  const double beta = props->beta;
  const double xoffset = props->xoffset;

  memset(C, 0, 21*sizeof(TacsScalar));

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
    TacsScalar Cj[21];
    props->props[j]->evalTangentStiffness3D(Cj);

    // Compute the penalty
    TacsScalar penalty = rho/(1.0 + q*(1.0 - rho));

    // Add up the contribution to the tangent stiffness matrix
    for ( int i = 0; i < 21; i++ ){
      C[i] += penalty*Cj[i];
    }
  }
}

/*
  Add the derivative of the product of the stress with the vector psi
  times alpha to the design variable array dvSens
*/
void TMROctConstitutive::addStressDVSens( int elemIndex,
                                          const double pt[],
                                          const TacsScalar X[],
                                          const TacsScalar e[],
                                          TacsScalar scale,
                                          const TacsScalar psi[],
                                          int dvLen,
                                          TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;
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

    TacsScalar C[21], s[6];
    props->props[j]->evalTangentStiffness3D(C);
    s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
    s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
    s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
    s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
    s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
    s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

    TacsScalar product = scale*dpenalty*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                                         s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);

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
void TMROctConstitutive::evalThermalStrain( int elemIndex,
                                            const double pt[],
                                            const TacsScalar X[],
                                            TacsScalar theta,
                                            TacsScalar e[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

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

    TacsScalar et[6];
    props->props[j]->evalThermalStrain3D(et);
    e[0] += rho*theta*et[0];
    e[1] += rho*theta*et[1];
    e[2] += rho*theta*et[2];
    e[3] += rho*theta*et[3];
    e[4] += rho*theta*et[4];
    e[5] += rho*theta*et[5];
  }
}

void TMROctConstitutive::addThermalStrainDVSens( int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 TacsScalar theta,
                                                 const TacsScalar psi[],
                                                 int dvLen,
                                                 TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;

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

    TacsScalar et[6];
    props->props[j]->evalThermalStrain3D(et);

    TacsScalar product = theta*(et[0]*psi[0] + et[1]*psi[1] + et[2]*psi[2] +
                                et[3]*psi[3] + et[4]*psi[4] + et[5]*psi[5]);

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
void TMROctConstitutive::evalHeatFlux( int elemIndex,
                                       const double pt[],
                                       const TacsScalar X[],
                                       const TacsScalar grad[],
                                       TacsScalar flux[] ){
  TacsScalar C[6];
  evalTangentHeatFlux(elemIndex, pt, X, C);
  flux[0] = C[0]*grad[0] + C[1]*grad[1] + C[2]*grad[2];
  flux[1] = C[1]*grad[0] + C[3]*grad[1] + C[4]*grad[2];
  flux[2] = C[2]*grad[0] + C[4]*grad[1] + C[5]*grad[2];
}

// Evaluate the tangent of the heat flux
void TMROctConstitutive::evalTangentHeatFlux( int elemIndex,
                                              const double pt[],
                                              const TacsScalar X[],
                                              TacsScalar C[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;
  const double qcond = props->qcond;
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
    props->props[j]->evalTangentHeatFlux3D(Cj);

    // Compute the penalty
    TacsScalar penalty = rho/(1.0 + qcond*(1.0 - rho));

    // Add up the contribution to the tangent stiffness matrix
    for ( int i = 0; i < 6; i++ ){
      C[i] += penalty*Cj[i];
    }
  }
}

// Add the derivative of the heat flux
void TMROctConstitutive::addHeatFluxDVSens( int elemIndex,
                                            const double pt[],
                                            const TacsScalar X[],
                                            const TacsScalar grad[],
                                            TacsScalar scale,
                                            const TacsScalar psi[],
                                            int dvLen,
                                            TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;
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

    TacsScalar C[6], flux[3];
    props->props[j]->evalTangentHeatFlux3D(C);
    flux[0] = C[0]*grad[0] + C[1]*grad[1] + C[2]*grad[2];
    flux[1] = C[1]*grad[0] + C[3]*grad[1] + C[4]*grad[2];
    flux[2] = C[2]*grad[0] + C[4]*grad[1] + C[5]*grad[2];

    TacsScalar product = scale*dpenalty*(flux[0]*psi[0] +
                                         flux[1]*psi[1] +
                                         flux[2]*psi[2]);

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
TacsScalar TMROctConstitutive::evalFailure( int elemIndex,
                                            const double pt[],
                                            const TacsScalar X[],
                                            const TacsScalar e[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;
  const double eps = props->eps;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  if (nvars == 1){
    TacsScalar C[21];
    props->props[0]->evalTangentStiffness3D(C);

    TacsScalar s[6];
    s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
    s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
    s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
    s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
    s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
    s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

    TacsScalar rho = 0.0;
    for ( int i = 0; i < len; i++ ){
      rho += N[i]*xptr[i];
    }

    TacsScalar r_factor = rho/(eps*(1.0 - rho) + rho);

    TacsScalar fail = props->props[0]->vonMisesFailure3D(s);

    return r_factor*fail;
  }
  else {
    const double ksWeight = props->ksWeight;
    TacsScalar *fail = temp_array;
    TacsScalar max_fail = -1e20;

    for ( int j = 0; j < nmats; j++ ){
      TacsScalar C[21];
      props->props[j]->evalTangentStiffness3D(C);

      TacsScalar s[6];
      s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
      s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
      s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
      s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
      s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
      s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

      TacsScalar rho = 0.0;
      for ( int i = 0; i < len; i++ ){
        rho += N[i]*xptr[nvars*i + j+1];
      }

      TacsScalar r_factor = rho/(eps*(1.0 - rho) + rho);

      fail[j] = r_factor*props->props[j]->vonMisesFailure3D(s);
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
void TMROctConstitutive::addFailureDVSens( int elemIndex,
                                           const double pt[],
                                           const TacsScalar X[],
                                           const TacsScalar e[],
                                           TacsScalar scale,
                                           int dvLen,
                                           TacsScalar dfdx[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;
  const double eps = props->eps;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  if (nvars == 1){
    TacsScalar C[21];
    props->props[0]->evalTangentStiffness3D(C);

    TacsScalar s[6];
    s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
    s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
    s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
    s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
    s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
    s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

    TacsScalar rho = 0.0;
    for ( int i = 0; i < len; i++ ){
      rho += N[i]*xptr[i];
    }

    // Compute the derivative of the stress relaxation factor
    TacsScalar d = 1.0/(eps*(1.0 - rho) + rho);
    TacsScalar r_factor_sens = eps*d*d;

    // Compute the contribution to the von Mises failure
    TacsScalar fail = props->props[0]->vonMisesFailure3D(s);

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
      TacsScalar C[21];
      props->props[j]->evalTangentStiffness3D(C);

      TacsScalar s[6];
      s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
      s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
      s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
      s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
      s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
      s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

      // Compute the density for the j-th material
      rho[j] = 0.0;
      for ( int i = 0; i < len; i++ ){
        rho[j] += N[i]*xptr[nvars*i + j+1];
      }

      // Compute the stress relaxation factor
      TacsScalar r_factor = rho[j]/(eps*(1.0 - rho[j]) + rho[j]);

      // Compute the failure value
      fail[j] = r_factor*props->props[j]->vonMisesFailure3D(s);
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
      TacsScalar C[21];
      props->props[j]->evalTangentStiffness3D(C);

      TacsScalar s[6];
      s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
      s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
      s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
      s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
      s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
      s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

      // Compute the derivative of the stress relaxation factor
      TacsScalar d = 1.0/(eps*(1.0 - rho[j]) + rho[j]);
      TacsScalar r_factor_sens = eps*d*d;

      // Compute the contribution to the von Mises failure
      TacsScalar fail_val = props->props[j]->vonMisesFailure3D(s);

      // Add the derivative
      TacsScalar contrib = scale*r_factor_sens*fail_val*(fail[j]/ksSum);
      for ( int i = 0; i < len; i++ ){
        dfdx[i*nvars + j+1] += contrib*N[i];
      }
    }
  }
}

TacsScalar TMROctConstitutive::evalFailureStrainSens( int elemIndex,
                                                      const double pt[],
                                                      const TacsScalar X[],
                                                      const TacsScalar e[],
                                                      TacsScalar dfde[] ){
  const int order = forest->getMeshOrder();
  const int len = order*order*order;
  const double eps = props->eps;

  // Evaluate the shape functions
  forest->evalInterp(pt, N);

  // Get the design variable values
  TacsScalar *xptr = &x[nvars*len*elemIndex];

  if (nvars == 1){
    TacsScalar C[21];
    props->props[0]->evalTangentStiffness3D(C);

    TacsScalar s[6];
    s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
    s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
    s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
    s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
    s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
    s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

    TacsScalar t[6];
    TacsScalar fail = props->props[0]->vonMisesFailure3DStressSens(s, t);

    dfde[0] = C[0]*t[0] + C[1]*t[1]  + C[2]*t[2]  + C[3]*t[3]  + C[4]*t[4]  + C[5]*t[5];
    dfde[1] = C[1]*t[0] + C[6]*t[1]  + C[7]*t[2]  + C[8]*t[3]  + C[9]*t[4]  + C[10]*t[5];
    dfde[2] = C[2]*t[0] + C[7]*t[1]  + C[11]*t[2] + C[12]*t[3] + C[13]*t[4] + C[14]*t[5];
    dfde[3] = C[3]*t[0] + C[8]*t[1]  + C[12]*t[2] + C[15]*t[3] + C[16]*t[4] + C[17]*t[5];
    dfde[4] = C[4]*t[0] + C[9]*t[1]  + C[13]*t[2] + C[16]*t[3] + C[18]*t[4] + C[19]*t[5];
    dfde[5] = C[5]*t[0] + C[10]*t[1] + C[14]*t[2] + C[17]*t[3] + C[19]*t[4] + C[20]*t[5];

    TacsScalar rho = 0.0;
    for ( int i = 0; i < len; i++ ){
      rho += N[i]*xptr[i];
    }

    TacsScalar r_factor = rho/(eps*(1.0 - rho) + rho);

    for ( int i = 0; i < 6; i++ ){
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
      TacsScalar C[21];
      props->props[j]->evalTangentStiffness3D(C);

      TacsScalar s[6];
      s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
      s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
      s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
      s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
      s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
      s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

      // Compute the density for the j-th material
      rho[j] = 0.0;
      for ( int i = 0; i < len; i++ ){
        rho[j] += N[i]*xptr[nvars*i + j+1];
      }

      // Compute the stress relaxation factor
      TacsScalar r_factor = rho[j]/(eps*(1.0 - rho[j]) + rho[j]);

      // Compute the failure value
      fail[j] = r_factor*props->props[j]->vonMisesFailure3D(s);
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

    memset(dfde, 0, 6*sizeof(TacsScalar));
    for ( int j = 0; j < nmats; j++ ){
      TacsScalar C[21];
      props->props[j]->evalTangentStiffness3D(C);

      TacsScalar s[6];
      s[0] = C[0]*e[0] + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
      s[1] = C[1]*e[0] + C[6]*e[1]  + C[7]*e[2]  + C[8]*e[3]  + C[9]*e[4]  + C[10]*e[5];
      s[2] = C[2]*e[0] + C[7]*e[1]  + C[11]*e[2] + C[12]*e[3] + C[13]*e[4] + C[14]*e[5];
      s[3] = C[3]*e[0] + C[8]*e[1]  + C[12]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
      s[4] = C[4]*e[0] + C[9]*e[1]  + C[13]*e[2] + C[16]*e[3] + C[18]*e[4] + C[19]*e[5];
      s[5] = C[5]*e[0] + C[10]*e[1] + C[14]*e[2] + C[17]*e[3] + C[19]*e[4] + C[20]*e[5];

      TacsScalar t[6];
      props->props[j]->vonMisesFailure3DStressSens(s, t);

      // Compute the stress relaxation factor
      TacsScalar r_factor = rho[j]/(eps*(1.0 - rho[j]) + rho[j]);

      TacsScalar scale = r_factor*fail[j]/ksSum;
      dfde[0] += scale*(C[0]*t[0] + C[1]*t[1]  + C[2]*t[2]  + C[3]*t[3]  + C[4]*t[4]  + C[5]*t[5]);
      dfde[1] += scale*(C[1]*t[0] + C[6]*t[1]  + C[7]*t[2]  + C[8]*t[3]  + C[9]*t[4]  + C[10]*t[5]);
      dfde[2] += scale*(C[2]*t[0] + C[7]*t[1]  + C[11]*t[2] + C[12]*t[3] + C[13]*t[4] + C[14]*t[5]);
      dfde[3] += scale*(C[3]*t[0] + C[8]*t[1]  + C[12]*t[2] + C[15]*t[3] + C[16]*t[4] + C[17]*t[5]);
      dfde[4] += scale*(C[4]*t[0] + C[9]*t[1]  + C[13]*t[2] + C[16]*t[3] + C[18]*t[4] + C[19]*t[5]);
      dfde[5] += scale*(C[5]*t[0] + C[10]*t[1] + C[14]*t[2] + C[17]*t[3] + C[19]*t[4] + C[20]*t[5]);
    }

    return ksFail;
  }

  return 0.0;
}

const char* TMROctConstitutive::getObjectName(){
  return "TMROctConstitutive";
}
