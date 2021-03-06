/*
  This file is part of the package TMR for adaptive mesh refinement.

  Copyright (C) 2018 Georgia Tech Research Corporation.
  Additional copyright (C) 2018 Graeme Kennedy.
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

#include "TMRCoupledThermoOctStiffness.h"
#include "FElibrary.h"
#include "YSlibrary.h"
/*
  Create the octree stiffness object based on an interpolation from
  the filter variables
*/
TMRCoupledThermoOctStiffness::TMRCoupledThermoOctStiffness( TMRIndexWeight *_weights,
                                                            int _nweights,
                                                            TMRStiffnessProperties *_props,
                                                            TMROctForest *_filter, 
                                                            int *_index ){
  // Record the density, Poisson ratio, D and the shear modulus
  props = _props;
  props->incref();

  // Set the weights/local indices
  nweights = _nweights;
  weights = NULL;
  if (_weights){
    weights = new TMRIndexWeight[ nweights ];
    memcpy(weights, _weights, nweights*sizeof(TMRIndexWeight));
  }
  index = NULL;
  if (_index){
    index = new int[ nweights ];
    memcpy(index, _index, nweights*sizeof(int));
  }  
  if (weights && index){
    printf("Define only using TMRIndexWeight or Index\n");
  }
  // Copy over the filter if provided
  filter = NULL;
  if (_filter){
    filter = _filter;
    filter->incref();
  }
  
  if (weights && filter){
    index = new int[ nweights ];
    for (int i = 0; i < nweights; i++){
      index[i] = weights[i].index;
    }    
  }
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
  // Set the dv vector to NULL
  dv = new TacsScalar[ nvars*nweights ];
  for (int i = 0; i < nvars*nweights; i++){
    dv[i] = 0.95;
  }
}
/*
  Decref the props
*/
TMRCoupledThermoOctStiffness::~TMRCoupledThermoOctStiffness(){
  props->decref();
  if (weights){
    delete [] weights;
  }
  if (index){
    delete [] index;
  }
  if (filter){
    filter->decref();
  }
  if (dv){
    delete [] dv;
  }
}

/*
  Loop over the design variable inputs and compute the local value of
  the density
*/
void TMRCoupledThermoOctStiffness::setDesignVars( const TacsScalar xdv[], 
                                                  int numDVs ){
  const double beta = props->beta;
  const double xoffset = props->xoffset;
  const int use_project = props->use_project;
  
  if (weights){
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += weights[i].weight*xdv[nvars*weights[i].index + j];
        dv[nvars*i+j] = xdv[nvars*weights[i].index + j];
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
  else { // Check for dependent nodes
    for ( int j = 0; j < nvars; j++ ){
      for ( int i = 0; i < nweights; i++ ){
        dv[nvars*i+j] = xdv[nvars*index[i] + j];
      }
    }
  }
}

/*
  Get the design variable values

  This is not possible to back out once the density is computed, so
  instead we use a constant value of 0.95 as the output.
*/
void TMRCoupledThermoOctStiffness::getDesignVars( TacsScalar xdv[], 
                                                  int numDVs ){
  if (weights){
    if (nvars == 1){
      for ( int i = 0; i < nweights; i++ ){
        if (weights[i].index >= 0 && weights[i].index < numDVs){
          xdv[weights[i].index] = 0.95;
        }
      }
    }
    else {
      for ( int j = 0; j < nvars; j++ ){
        TacsScalar value = 0.05;
        if (j >= 1){
          value = 0.95/(nvars-1);
        }
        for ( int i = 0; i < nweights; i++ ){
          if (weights[i].index >= 0 && weights[i].index < numDVs){
            xdv[nvars*weights[i].index + j] = value;
          }
        }
      }
    }
  }
  else { // Using Bernstein polynomial; check for dependent nodes
    if (nvars == 1){
      for ( int i = 0; i < nweights; i++ ){
        xdv[index[i]] = 0.95;
      }
    }
    else {
      for ( int j = 0; j < nvars; j++ ){
        TacsScalar value = 0.05;
        if (j >= 1){
          value = 0.95/(nvars-1);
        }
        for ( int i = 0; i < nweights; i++ ){
          xdv[nvars*index[i] + j] = value;          
        }
      }
    }
  }
}

/*
  Retrieve the range of the design variable values
*/
void TMRCoupledThermoOctStiffness::getDesignVarRange( TacsScalar lb[], 
                                                      TacsScalar ub[], 
                                                      int numDVs ){
  if (weights){
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
  else { // Using Bernstein polynomail; check for dependent nodes
    if (nvars == 1){
      for ( int i = 0; i < nweights; i++ ){
        lb[index[i]] = 0.0;
        ub[index[i]] = 1.0;
      }
    }
    else {
      for ( int j = 0; j < nvars; j++ ){
        for ( int i = 0; i < nweights; i++ ){
          lb[nvars*index[i] + j] = 0.0;
          ub[nvars*index[i] + j] = 1e30;
        }
      }
    }
  }
}
/*
  Compute the density of the material as a linear function of the
  element-wise density
*/
void TMRCoupledThermoOctStiffness::getPointwiseMass( const double pt[], 
                                                     TacsScalar mass[] ){
  // Get the density value at isoparametric point pt
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];     
      }
    }
  }
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
void TMRCoupledThermoOctStiffness::addPointwiseMassDVSens( const double pt[], 
                                                           const TacsScalar alpha[],
                                                           TacsScalar dvSens[], 
                                                           int dvLen ){
  // Get the density value at isoparametric point pt
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    if (nvars == 1){
      TacsScalar scale = props->density[0]*alpha[0];
      for ( int i = 0; i < nweights; i++ ){
        dvSens[index[i]] += (scale*N[i]);
      }
    }
    else {
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar scale = props->density[j-1]*alpha[0];
        for ( int i = 0; i < nweights; i++ ){
          dvSens[nvars*index[i] + j] += (scale*N[i]);
        }
      }
    }
  }
  else{
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
}

/*
  Compute the matrix vector product of D*e
*/
void TMRCoupledThermoOctStiffness::calculateStress( const double pt[],
                                                    const TacsScalar e[], 
                                                    TacsScalar s[] ){
  const double k0 = props->k0;
  const double qs = props->q;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];
    }
  }
  
  if (nvars == 1){
    TacsScalar penalty = rho[0]/(1.0+qs*(1.0-rho[0]));
    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];
    // Add the penalized value
    TacsScalar Dp = (penalty+k0)*D;
    TacsScalar Gp = (penalty+k0)*G;
    // Compute the resulting matrix-vector product D*e
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
      TacsScalar penalty = rho[j]/(1.0+qs*(1.0-rho[j]));
      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];

      // Add the penalized value
      TacsScalar Dp = (penalty+k0)*D;
      TacsScalar Gp = (penalty+k0)*G;
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
void TMRCoupledThermoOctStiffness::addStressDVSens( const double pt[], 
                                                    const TacsScalar e[], 
                                                    TacsScalar alpha, 
                                                    const TacsScalar psi[], 
                                                    TacsScalar fdvSens[],
                                                    int dvLen ){
  const double q = props->q;
  const double beta = props->beta;
  const double xoffset = props->xoffset;
  const int use_project = props->use_project;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];        
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];      
    }
    if (nvars == 1){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (q + 1.0)/((1.0 + q*(1.0 - rho[0]))*(1.0 + q*(1.0 - rho[0])));
      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[6];
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      // Compute the penalized value
      TacsScalar Dp = alpha*penalty*D;
      TacsScalar Gp = alpha*penalty*G;
      // Compute the product dD/dx*(B*u)
      s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
      s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
      s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
      s[3] = Gp*e[3];
      s[4] = Gp*e[4];
      s[5] = Gp*e[5];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                            s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[index[i]] += N[i]*product;
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

        TacsScalar s[6];
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];
        // Add the result to the derivative
        TacsScalar Dp = alpha*penalty*D;
        TacsScalar Gp = alpha*penalty*G;
        // Compute the product dD/dx*(B*u)
        s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
        s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
        s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
        s[3] = Gp*e[3];
        s[4] = Gp*e[4];
        s[5] = Gp*e[5];
        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                              s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*index[i] + j] += N[i]*product;
        }
      }
    }
  }
  else {
    if (nvars == 1){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (q + 1.0)/((1.0 + q*(1.0 - rho[0]))*(1.0 + q*(1.0 - rho[0])));
      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[6];
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      // Compute the penalized value
      TacsScalar Dp = alpha*penalty*D;
      TacsScalar Gp = alpha*penalty*G;
      // Compute the product dD/dx*(B*u)
      s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
      s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
      s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
      s[3] = Gp*e[3];
      s[4] = Gp*e[4];
      s[5] = Gp*e[5];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                            s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
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

        TacsScalar s[6];
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];
        // Add the result to the derivative
        TacsScalar Dp = alpha*penalty*D;
        TacsScalar Gp = alpha*penalty*G;
        // Compute the product dD/dx*(B*u)
        s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
        s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
        s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
        s[3] = Gp*e[3];
        s[4] = Gp*e[4];
        s[5] = Gp*e[5];
        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                              s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*weights[i].index + j] += weights[i].weight*product;
        }
      }
    }
  }
}
/*
  Compute the matrix vector product of L*e
*/
void TMRCoupledThermoOctStiffness::calculateThermal(const double pt[],
                                                    const TacsScalar e[], 
                                                    TacsScalar s[]){
  const double k0 = props->k0;
  const double qtemp = props->qtemp;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      rho[j] = x[j];
    }  
  }
  if (nvars == 1){
    TacsScalar penalty = rho[0]/(1.0+qtemp*(1.0-rho[0]));
    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];
    TacsScalar aT = props->aT[0];
    // Add the penalized value
    TacsScalar Dp = (penalty+k0)*D*aT;
    TacsScalar Gp = (penalty+k0)*G;
    // Compute the resulting matrix-vector product D*e
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
      TacsScalar penalty = rho[j]/(1.0+qtemp*(1.0-rho[j]));
      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];
      TacsScalar aT = props->aT[j-1];
      // Add the penalized value
      TacsScalar Dp = (penalty+k0)*D*aT;
      TacsScalar Gp = (penalty+k0)*G;

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
void TMRCoupledThermoOctStiffness::addThermalDVSens( const double pt[], 
                                                     const TacsScalar e[], 
                                                     TacsScalar alpha, 
                                                     const TacsScalar psi[], 
                                                     TacsScalar fdvSens[],
                                                     int dvLen ){
  const double qtemp = props->qtemp;
  const double beta = props->beta;
  const double xoffset = props->xoffset;
  const int use_project = props->use_project;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];
    }
    if (nvars == 1){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (qtemp + 1.0)/((1.0 + qtemp*(1.0 - rho[0]))*(1.0 + qtemp*(1.0 - rho[0])));

      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[6];
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      TacsScalar aT = props->aT[0];
      // Compute the penalized value
      TacsScalar Dp = alpha*penalty*D*aT;
      TacsScalar Gp = alpha*penalty*G;
      // Compute the product dD/dx*(B*u)
      s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
      s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
      s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
      s[3] = Gp*e[3];
      s[4] = Gp*e[4];
      s[5] = Gp*e[5];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                            s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
    
      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[index[i]] += N[i]*product;
      }
    }
    else {
      // Need to include aT[j] here
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar penalty =
          (qtemp + 1.0)/((1.0 + qtemp*(1.0 - rho[j]))*(1.0 + qtemp*(1.0 - rho[j])));

        // Add the derivative of the projection
        if (use_project){
          penalty *= beta*exp(-beta*(x[j] - xoffset))*rho[j]*rho[j];
        }
        TacsScalar s[6];
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];
        TacsScalar aT = props->aT[j-1];
        // Add the result to the derivative
        TacsScalar Dp = alpha*penalty*D*aT;
        TacsScalar Gp = alpha*penalty*G;
        // Compute the product dD/dx*(B*u)
        s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
        s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
        s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
        s[3] = Gp*e[3];
        s[4] = Gp*e[4];
        s[5] = Gp*e[5];

        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                              s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
      
        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*index[i] + j] += N[i]*product;
        }
      }
    }
  }
  else {
    if (nvars == 1){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (qtemp + 1.0)/((1.0 + qtemp*(1.0 - rho[0]))*(1.0 + qtemp*(1.0 - rho[0])));

      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[6];
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      TacsScalar aT = props->aT[0];
      // Compute the penalized value
      TacsScalar Dp = alpha*penalty*D*aT;
      TacsScalar Gp = alpha*penalty*G;
      // Compute the product dD/dx*(B*u)
      s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
      s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
      s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
      s[3] = Gp*e[3];
      s[4] = Gp*e[4];
      s[5] = Gp*e[5];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                            s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
    
      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[weights[i].index] += weights[i].weight*product;
      }
    }
    else {
      // Need to include aT[j] here
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar penalty =
          (qtemp + 1.0)/((1.0 + qtemp*(1.0 - rho[j]))*(1.0 + qtemp*(1.0 - rho[j])));

        // Add the derivative of the projection
        if (use_project){
          penalty *= beta*exp(-beta*(x[j] - xoffset))*rho[j]*rho[j];
        }
        TacsScalar s[6];
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];
        TacsScalar aT = props->aT[j-1];
        // Add the result to the derivative
        TacsScalar Dp = alpha*penalty*D*aT;
        TacsScalar Gp = alpha*penalty*G;
        // Compute the product dD/dx*(B*u)
        s[0] = Dp*((1.0 - nu)*e[0] + nu*(e[1] + e[2]));
        s[1] = Dp*((1.0 - nu)*e[1] + nu*(e[0] + e[2]));
        s[2] = Dp*((1.0 - nu)*e[2] + nu*(e[0] + e[1]));
        s[3] = Gp*e[3];
        s[4] = Gp*e[4];
        s[5] = Gp*e[5];

        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
                              s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5]);
      
        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*weights[i].index + j] += weights[i].weight*product;
        }
      }
    }
  }  
}

/*
  Compute the conduction terms D_{th}*u
*/
void TMRCoupledThermoOctStiffness::calculateConduction( const double pt[],
                                                        const TacsScalar e[], 
                                                        TacsScalar s[] ){
  const double k0 = props->k0;
  const double qcond = props->qcond;

  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      rho[j] = x[j];
    }  
  }

  // Assuming isotropic conduction i.e. only diagonal terms
  if (nvars == 1){
    TacsScalar penalty = rho[0]/(1.0+qcond*(1.0-rho[0]));
    TacsScalar kcond = props->kcond[0];
    TacsScalar Cp = (kcond+k0)*penalty;
    s[0] = Cp*e[0];
    s[1] = Cp*e[1];
    s[2] = Cp*e[2];
  }
  else {
    s[0] = s[1] = s[2] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      TacsScalar penalty = rho[j]/(1.0+qcond*(1.0-rho[j]));
      // Extract the properties
      TacsScalar kcond = props->kcond[j-1];
      // Add the penalized value
      TacsScalar Cp = (penalty+k0)*kcond;
      s[0] += Cp*e[0];
      s[1] += Cp*e[1];
      s[2] += Cp*e[2];
    }
  }  
}
/*
  Add derivative w.r.t. the design variables due to the conduction term
*/
void TMRCoupledThermoOctStiffness::addConductionDVSens( const double pt[], 
                                                        const TacsScalar e[], 
                                                        TacsScalar alpha, 
                                                        const TacsScalar psi[], 
                                                        TacsScalar fdvSens[],
                                                        int dvLen ){
  const double qcond = props->qcond;
  const double beta = props->beta;
  const double xoffset = props->xoffset;
  const int use_project = props->use_project;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];      
    }
    if (nvars == 1){
      TacsScalar penalty =
        (qcond + 1.0)/((1.0 + qcond*(1.0 - rho[0]))*(1.0 + qcond*(1.0 - rho[0])));

      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[3];
      // Compute the product dD/dx*(B*u)
      TacsScalar kcond = props->kcond[0];
      TacsScalar Cp = alpha*kcond*penalty;
      s[0] = Cp*e[0];
      s[1] = Cp*e[1];
      s[2] = Cp*e[2];
    
      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[index[i]] +=  N[i]*product;
      }
    }
    else {
      for ( int j = 1; j < nvars; j++ ){
        // Compute the derivative of the penalization with respect to
        // the projected density
        TacsScalar penalty =
          (qcond + 1.0)/((1.0 + qcond*(1.0 - rho[j]))*(1.0 + qcond*(1.0 - rho[j])));

        // Add the derivative of the projection
        if (use_project){
          penalty *= beta*exp(-beta*(x[j] - xoffset))*rho[j]*rho[j];
        }
      
        // Extract the properties
        TacsScalar kcond = props->kcond[j-1];
        // Add the penalized value
        TacsScalar Cp = alpha*penalty*kcond;
        TacsScalar s[3];
        s[0] = Cp*e[0];
        s[1] = Cp*e[1];
        s[2] = Cp*e[2];
      
        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
      
        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*index[i] + j] += N[i]*product;
        }
      }
    }
  }
  else {
    if (nvars == 1){
      TacsScalar penalty =
        (qcond + 1.0)/((1.0 + qcond*(1.0 - rho[0]))*(1.0 + qcond*(1.0 - rho[0])));

      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[3];
      // Compute the product dD/dx*(B*u)
      TacsScalar kcond = props->kcond[0];
      TacsScalar Cp = alpha*kcond*penalty;
      s[0] = Cp*e[0];
      s[1] = Cp*e[1];
      s[2] = Cp*e[2];
    
      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[weights[i].index] +=  weights[i].weight*product;
      }
    }
    else {
      for ( int j = 1; j < nvars; j++ ){
        // Compute the derivative of the penalization with respect to
        // the projected density
        TacsScalar penalty =
          (qcond + 1.0)/((1.0 + qcond*(1.0 - rho[j]))*(1.0 + qcond*(1.0 - rho[j])));

        // Add the derivative of the projection
        if (use_project){
          penalty *= beta*exp(-beta*(x[j] - xoffset))*rho[j]*rho[j];
        }
      
        // Extract the properties
        TacsScalar kcond = props->kcond[j-1];
        // Add the penalized value
        TacsScalar Cp = alpha*penalty*kcond;
        TacsScalar s[3];
        s[0] = Cp*e[0];
        s[1] = Cp*e[1];
        s[2] = Cp*e[2];
      
        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
      
        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*weights[i].index + j] += weights[i].weight*product;
        }
      }
    }
  }
}

/*
  Returns the unpenalized term alpha for candidate material j
*/
TacsScalar TMRCoupledThermoOctStiffness::getEffThermalAlpha( int vars_j ){
  // Let calculateThermal take care of aT[j]
  if (vars_j > 0){
    // Extract the properties
    TacsScalar aT = props->aT[vars_j-1];    
    return aT;
  }
  else {
    return props->aT[0];
  }
}

// Evaluate the failure criteria
// Change in temperature for the element is appended to the strain vector
void TMRCoupledThermoOctStiffness::failure( const double pt[], 
                                            const TacsScalar T[],
                                            const TacsScalar e[],
                                            TacsScalar * fail ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];
    }
  }
  if (nvars == 1){
    // Use the von Mises failure criterion
    // Compute the relaxation factor
    TacsScalar r_factor = 1.0;
    if (eps > 0.0){
      r_factor = rho[0]/(eps*(1.0-rho[0])+rho[0]);
    }
    else {
      TacsScalar p = -1.0*eps;
      r_factor = pow(rho[0],p);
    }
    TacsScalar s[6], eff_e[6];
    // Extract the properties
    TacsScalar nu = props->nu[0];
    TacsScalar D = props->D[0];
    TacsScalar G = props->G[0];
    TacsScalar ys = props->ys[0];
    TacsScalar aT = props->aT[0];

    eff_e[0] = e[0]-aT*x[0]*T[0];
    eff_e[1] = e[1]-aT*x[0]*T[0];
    eff_e[2] = e[2]-aT*x[0]*T[0];
    eff_e[3] = e[3];
    eff_e[4] = e[4];
    eff_e[5] = e[5];
     
    // Compute the resulting stress
    s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
    s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
    s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
    s[3] = G*eff_e[3];
    s[4] = G*eff_e[4];
    s[5] = G*eff_e[5];
 
    *fail = r_factor*VonMisesFailure3D(s,ys);
  }
  else {
    *fail = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = rho[j]/(eps*(1.0-rho[j])+rho[j]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(rho[j],p);
      }      
      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];
      TacsScalar ys = props->ys[j-1];
      TacsScalar aT = props->aT[j-1];

      // Compute effective strain = B*u-alpha*dT
      TacsScalar s[6], eff_e[6];
      eff_e[0] = e[0]-aT*T[0];
      eff_e[1] = e[1]-aT*T[0];
      eff_e[2] = e[2]-aT*T[0];
      eff_e[3] = e[3]*1.0;
      eff_e[4] = e[4]*1.0;
      eff_e[5] = e[5]*1.0;
      // Compute the contribution of failure from each material
      s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
      s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
      s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
      s[3] = G*eff_e[3];
      s[4] = G*eff_e[4];
      s[5] = G*eff_e[5];
      
      *fail += r_factor*VonMisesFailure3D(s,ys);
    }
  }
}
// Evaluate the failure criteria w.r.t. design variables
void TMRCoupledThermoOctStiffness::addFailureDVSens( const double pt[], 
                                                     const TacsScalar T[],
                                                     const TacsScalar e[],
                                                     TacsScalar alpha,
                                                     TacsScalar dvSens[], 
                                                     int dvLen ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];        
      }
    }
    if (nvars == 1){
      // Compute the relaxation factor
      TacsScalar r_factor_sens = 0.0;
      if (eps > 0.0){
        TacsScalar d = 1.0/(eps*(1.0-x[0])+x[0]);
        r_factor_sens = eps*d*d;
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor_sens = pow(x[0],p-1)*p;
      }
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      TacsScalar ys = props->ys[0];
      TacsScalar aT = props->aT[0];

      TacsScalar s[6], eff_e[6];
      // Compute effective strain = B*u-alpha*dT
      eff_e[0] = e[0]-aT*x[0]*T[0];
      eff_e[1] = e[1]-aT*x[0]*T[0];
      eff_e[2] = e[2]-aT*x[0]*T[0];
      eff_e[3] = e[3]*1.0;
      eff_e[4] = e[4]*1.0;
      eff_e[5] = e[5]*1.0;
      // Compute the contribution of failure from each material
      s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
      s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
      s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
      s[3] = G*eff_e[3];
      s[4] = G*eff_e[4];
      s[5] = G*eff_e[5];

      TacsScalar fail = VonMisesFailure3D(s,ys);
      TacsScalar inner = alpha*r_factor_sens*fail;
      for ( int i = 0; i < nweights; i++ ){
        dvSens[index[i]] += N[i]*inner;
      }
    }
    else {
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar r_factor_sens = 0.0;
        if (eps > 0.0){
          TacsScalar d = 1.0/(eps*(1.0-x[j])+x[j]);
          r_factor_sens = eps*d*d;
        }
        else {
          TacsScalar p = -1.0*eps;
          r_factor_sens = pow(x[j],p-1)*p;
        }
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];
        TacsScalar ys = props->ys[j-1];
        TacsScalar aT = props->aT[j-1];
      
        TacsScalar s[6], eff_e[6];
        // Compute effective strain = B*u-alpha*dT
        eff_e[0] = e[0]-aT*T[0];
        eff_e[1] = e[1]-aT*T[0];
        eff_e[2] = e[2]-aT*T[0];
        eff_e[3] = e[3]*1.0;
        eff_e[4] = e[4]*1.0;
        eff_e[5] = e[5]*1.0;
        // Compute the contribution of failure from each material
        s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
        s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
        s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
        s[3] = G*eff_e[3];
        s[4] = G*eff_e[4];
        s[5] = G*eff_e[5];
        
        // For the derivative dr/dx*sigma_vm/ys
        TacsScalar fail = VonMisesFailure3D(s,ys);
          
        // Compute the inner product
        TacsScalar inner = alpha*(r_factor_sens*fail);
        
        for ( int i = 0; i < nweights; i++ ){
          dvSens[nvars*index[i] + j] += N[i]*inner;
        }
      }
    }    
  }
  else{
    if (nvars == 1){
      // Compute the relaxation factor
      TacsScalar r_factor_sens = 0.0;
      if (eps > 0.0){
        TacsScalar d = 1.0/(eps*(1.0-rho[0])+rho[0]);
        r_factor_sens = eps*d*d;
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor_sens = pow(rho[0],p-1)*p;
      }
    
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      TacsScalar ys = props->ys[0];
      TacsScalar aT = props->aT[0];
      TacsScalar s[6], eff_e[6];
      // Compute effective strain = B*u-alpha*dT
      eff_e[0] = e[0]-aT*T[0];
      eff_e[1] = e[1]-aT*T[0];
      eff_e[2] = e[2]-aT*T[0];
      eff_e[3] = e[3]*1.0;
      eff_e[4] = e[4]*1.0;
      eff_e[5] = e[5]*1.0;
      // Compute the contribution of failure from each material
      s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
      s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
      s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
      s[3] = G*eff_e[3];
      s[4] = G*eff_e[4];
      s[5] = G*eff_e[5];
    
      TacsScalar fail = VonMisesFailure3D(s,ys);
      TacsScalar inner = alpha*r_factor_sens*fail;
      for ( int i = 0; i < nweights; i++ ){
        dvSens[weights[i].index] += weights[i].weight*inner;
      }
    }
    else{
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar r_factor_sens = 0.0;
        if (eps > 0.0){
          TacsScalar d = 1.0/(eps*(1.0-rho[j])+rho[j]);
          r_factor_sens = eps*d*d;
        }
        else {
          TacsScalar p = -1.0*eps;
          r_factor_sens = pow(rho[j],p-1)*p;
        }
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];
        TacsScalar ys = props->ys[j-1];
        TacsScalar aT = props->aT[j-1];
        TacsScalar s[6], eff_e[6];
        // Compute effective strain = B*u-alpha*dT
        eff_e[0] = e[0]-aT*T[0];
        eff_e[1] = e[1]-aT*T[0];
        eff_e[2] = e[2]-aT*T[0];
        eff_e[3] = e[3]*1.0;
        eff_e[4] = e[4]*1.0;
        eff_e[5] = e[5]*1.0;
        // Compute the contribution of failure from each material
        s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
        s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
        s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
        s[3] = G*eff_e[3];
        s[4] = G*eff_e[4];
        s[5] = G*eff_e[5];

        // For the derivative dr/dx*sigma_vm/ys
        TacsScalar fail = VonMisesFailure3D(s,ys);
        TacsScalar inner = alpha*(r_factor_sens*fail);
        for ( int i = 0; i < nweights; i++ ){
          dvSens[nvars*weights[i].index + j] += weights[i].weight*inner;
        }
      }
    }
  }
}
void TMRCoupledThermoOctStiffness::failureStrainSens( const double pt[], 
                                                      const TacsScalar T[],
                                                      const TacsScalar e[],
                                                      TacsScalar sens[], 
                                                      int vars_j ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];        
      }  
    }
    if (nvars == 1){
      TacsScalar s[6], ps_sens[6], eff_e[6];
      // Use the von Mises failure criterion
      // Compute the relaxation factor
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = rho[0]/(eps*(1.0-x[0])+x[0]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(x[0],p);
      }
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      TacsScalar ys = props->ys[0];
      TacsScalar aT = props->aT[0];
      // Compute effective strain = B*u-alpha*dT
      eff_e[0] = e[0]-aT*T[0];
      eff_e[1] = e[1]-aT*T[0];
      eff_e[2] = e[2]-aT*T[0];
      eff_e[3] = e[3]*1.0;
      eff_e[4] = e[4]*1.0;
      eff_e[5] = e[5]*1.0;
      // Compute the contribution of failure from each material
      s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
      s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
      s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
      s[3] = G*eff_e[3];
      s[4] = G*eff_e[4];
      s[5] = G*eff_e[5];

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
    else{
      int j = vars_j*1;
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = x[j]/(eps*(1.0-x[j])+x[j]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(x[j],p);
      }

      TacsScalar s[6], ps_sens[6], eff_e[6];
      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];
      TacsScalar ys = props->ys[j-1];
      TacsScalar aT = props->aT[j-1];
      // Compute effective strain = B*u-alpha*dT
      eff_e[0] = e[0]-aT*T[0];
      eff_e[1] = e[1]-aT*T[0];
      eff_e[2] = e[2]-aT*T[0];
      eff_e[3] = e[3]*1.0;
      eff_e[4] = e[4]*1.0;
      eff_e[5] = e[5]*1.0;
    
      // Compute the contribution of failure from each material
      s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
      s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
      s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
      s[3] = G*eff_e[3];
      s[4] = G*eff_e[4];
      s[5] = G*eff_e[5];

      VonMisesFailure3DStressSens(ps_sens, s, ys);
      sens[0] = r_factor*D*((1.0-nu)*ps_sens[0] + 
                            nu*(ps_sens[1] + ps_sens[2]));
      sens[1] = r_factor*D*((1.0-nu)*ps_sens[1] + 
                            nu*(ps_sens[0] + ps_sens[2]));
      sens[2] = r_factor*D*((1.0-nu)*ps_sens[2] + 
                            nu*(ps_sens[0] + ps_sens[1]));
      
      sens[3] = r_factor*G*ps_sens[3];
      sens[4] = r_factor*G*ps_sens[4];
      sens[5] = r_factor*G*ps_sens[5];
    }
  }
  else {
    if (nvars == 1){
      TacsScalar s[6], ps_sens[6], eff_e[6];
      // Use the von Mises failure criterion
      // Compute the relaxation factor
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = rho[0]/(eps*(1.0-rho[0])+rho[0]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(rho[0],p);
      }
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      TacsScalar ys = props->ys[0];
      TacsScalar aT = props->aT[0];
      // Compute effective strain = B*u-alpha*dT
      eff_e[0] = e[0]-aT*T[0];
      eff_e[1] = e[1]-aT*T[0];
      eff_e[2] = e[2]-aT*T[0];
      eff_e[3] = e[3]*1.0;
      eff_e[4] = e[4]*1.0;
      eff_e[5] = e[5]*1.0;
      // Compute the contribution of failure from each material
      s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
      s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
      s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
      s[3] = G*eff_e[3];
      s[4] = G*eff_e[4];
      s[5] = G*eff_e[5];

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
    else{
      int j = vars_j*1;
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = rho[j]/(eps*(1.0-rho[j])+rho[j]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(rho[j],p);
      }
      TacsScalar s[6], ps_sens[6], eff_e[6];
      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];
      TacsScalar ys = props->ys[j-1];
      TacsScalar aT = props->aT[j-1];
      // Compute effective strain = B*u-alpha*dT
      eff_e[0] = e[0]-aT*T[0];
      eff_e[1] = e[1]-aT*T[0];
      eff_e[2] = e[2]-aT*T[0];
      eff_e[3] = e[3]*1.0;
      eff_e[4] = e[4]*1.0;
      eff_e[5] = e[5]*1.0;
    
      // Compute the contribution of failure from each material
      s[0] = D*((1.0 - nu)*eff_e[0] + nu*(eff_e[1] + eff_e[2]));
      s[1] = D*((1.0 - nu)*eff_e[1] + nu*(eff_e[0] + eff_e[2]));
      s[2] = D*((1.0 - nu)*eff_e[2] + nu*(eff_e[0] + eff_e[1]));
      s[3] = G*eff_e[3];
      s[4] = G*eff_e[4];
      s[5] = G*eff_e[5];

      VonMisesFailure3DStressSens(ps_sens, s, ys);
      sens[0] = r_factor*D*((1.0-nu)*ps_sens[0] + 
                            nu*(ps_sens[1] + ps_sens[2]));
      sens[1] = r_factor*D*((1.0-nu)*ps_sens[1] + 
                            nu*(ps_sens[0] + ps_sens[2]));
      sens[2] = r_factor*D*((1.0-nu)*ps_sens[2] + 
                            nu*(ps_sens[0] + ps_sens[1]));
      
      sens[3] = r_factor*G*ps_sens[3];
      sens[4] = r_factor*G*ps_sens[4];
      sens[5] = r_factor*G*ps_sens[5];
    }
  }
}
TacsScalar TMRCoupledThermoOctStiffness::getDVOutputValue( int dvIndex, 
                                                           const double pt[] ){
  // Compute the filter shape functions
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];
    }
  }
  if (dvIndex == 0){
    if (nvars > 1){
      return rho[1];
    }
    return rho[0];
  }
  else if (dvIndex == 1){
    if (nvars > 1){
      return rho[2];
    }
    return rho[0];
  }
  else if (dvIndex == 2){
    return rho[0];
  }  
  return 0;
}

// Evaluate the heat flux normal to the surface
// Change in temperature for the element is appended to the strain vector
void TMRCoupledThermoOctStiffness::heatflux( const double pt[], 
                                             const TacsScalar normal[],
                                             const TacsScalar e[],
                                             TacsScalar *qn ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];
    }
  }
  if (nvars == 1){
    // Use the von Mises failure criterion
    // Compute the relaxation factor
    TacsScalar r_factor = 1.0;
    if (eps > 0.0){
      r_factor = rho[0]/(eps*(1.0-rho[0])+rho[0]);
    }
    else {
      TacsScalar p = -1.0*eps;
      r_factor = pow(rho[0],p);
    }
    TacsScalar q[3];
    // Extract the properties
    TacsScalar kcond = props->kcond[0];
    q[0] = kcond*e[0];
    q[1] = kcond*e[1];
    q[2] = kcond*e[2];
 
    *qn = r_factor*(normal[0]*q[0] + normal[1]*q[1] + normal[2]*q[2]);
  }
  else {
    *qn = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = rho[j]/(eps*(1.0-rho[j])+rho[j]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(rho[j],p);
      }      
      // Extract the properties
      TacsScalar kcond = props->kcond[j-1];
      TacsScalar q[3];
      q[0] = kcond*e[0];
      q[1] = kcond*e[1];
      q[2] = kcond*e[2];
      
      *qn += r_factor*(normal[0]*q[0] + normal[1]*q[1] + normal[2]*q[2]);
    }
  }
}

// Evaluate the failure criteria w.r.t. design variables
void TMRCoupledThermoOctStiffness::addHeatFluxDVSens( const double pt[], 
                                                      const TacsScalar normal[],
                                                      const TacsScalar e[],
                                                      TacsScalar alpha,
                                                      TacsScalar dvSens[], 
                                                      int dvLen ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];        
      }
    }
    if (nvars == 1){
      // Compute the relaxation factor
      TacsScalar r_factor_sens = 0.0;
      if (eps > 0.0){
        TacsScalar d = 1.0/(eps*(1.0-x[0])+x[0]);
        r_factor_sens = eps*d*d;
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor_sens = pow(x[0],p-1)*p;
      }
      // Extract the properties
      TacsScalar kcond = props->kcond[0];
      TacsScalar q[3];
      q[0] = kcond*e[0];
      q[1] = kcond*e[1];
      q[2] = kcond*e[2];
      TacsScalar qn = (normal[0]*q[0] + normal[1]*q[1] +
                       normal[2]*q[2]); 

      TacsScalar inner = alpha*r_factor_sens*qn;
      for ( int i = 0; i < nweights; i++ ){
        dvSens[index[i]] += N[i]*inner;
      }
    }
    else {
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar r_factor_sens = 0.0;
        if (eps > 0.0){
          TacsScalar d = 1.0/(eps*(1.0-x[j])+x[j]);
          r_factor_sens = eps*d*d;
        }
        else {
          TacsScalar p = -1.0*eps;
          r_factor_sens = pow(x[j],p-1)*p;
        }
        TacsScalar kcond = props->kcond[j-1];
        TacsScalar q[3];
        q[0] = kcond*e[0];
        q[1] = kcond*e[1];
        q[2] = kcond*e[2];
        TacsScalar qn = (normal[0]*q[0] + normal[1]*q[1] +
                         normal[2]*q[2]); 
        
        TacsScalar inner = alpha*r_factor_sens*qn;
        for ( int i = 0; i < nweights; i++ ){
          dvSens[nvars*index[i] + j] += N[i]*inner;
        }
      }      
    }
  }
}
void TMRCoupledThermoOctStiffness::heatfluxStrainSens( const double pt[], 
                                                       const TacsScalar normal[],
                                                       TacsScalar sens[] ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];        
      }  
    }
    if (nvars == 1){
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = rho[0]/(eps*(1.0-x[0])+x[0]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(x[0],p);
      }
      // Extract the properties
      TacsScalar kcond = props->kcond[0];
      sens[0] = r_factor*kcond*normal[0];
      sens[1] = r_factor*kcond*normal[1];
      sens[2] = r_factor*kcond*normal[2];
    }
    else {
      for (int j = 1; j < nvars; j++){
        TacsScalar r_factor = 1.0;
        if (eps > 0.0){
          r_factor = x[j]/(eps*(1.0-x[j])+x[j]);
        }
        else {
          TacsScalar p = -1.0*eps;
          r_factor = pow(x[j],p);
        }
        // Extract the properties
        TacsScalar kcond = props->kcond[j-1];
        sens[0] += r_factor*kcond*normal[0];
        sens[1] += r_factor*kcond*normal[1];
        sens[2] += r_factor*kcond*normal[2];
      }      
    }    
  }
}
// Evaluate the failure criteria
// Change in temperature for the element is appended to the strain vector
void TMRCoupledThermoOctStiffness::maxtemp( const double pt[], 
                                            const TacsScalar maxtemp,
					    TacsScalar * fail,
					    int mat_j ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];
    }    
  }
  if (nvars == 1){
    // Compute the relaxation factor
    TacsScalar r_factor = 1.0;
    if (eps > 0.0){
      r_factor = rho[0]/(eps*(1.0-rho[0])+rho[0]);
    }
    else {
      TacsScalar p = -1.0*eps;
      r_factor = pow(rho[0],p);
    }    
    *fail = r_factor*maxtemp/props->Tmax[0];
  }
  else {
    *fail = 0.0;
    int j = mat_j+1;
    TacsScalar r_factor = 1.0;
    if (eps > 0.0){
      r_factor = rho[j]/(eps*(1.0-rho[j])+rho[j]);
    }
    else {
      TacsScalar p = -1.0*eps;
      r_factor = pow(rho[j],p);
    }
    
    *fail = r_factor*maxtemp/props->Tmax[mat_j];
    // for ( int j = 1; j < nvars; j++ ){
    //   TacsScalar r_factor = 1.0;
    //   if (eps > 0.0){
    //     r_factor = rho[j]/(eps*(1.0-rho[j])+rho[j]);
    //   }
    //   else {
    //     TacsScalar p = -1.0*eps;
    //     r_factor = pow(rho[j],p);
    //   }
    //   *fail += r_factor*maxtemp;
    // }
  }
}

void TMRCoupledThermoOctStiffness::addMaxTempDVSens( const double pt[], 
						     const TacsScalar maxtemp,
						     TacsScalar alpha,
						     TacsScalar dvSens[], int dvLen,
						     int mat_j ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];        
      }
      rho[j] = x[j];
    }
    if (nvars == 1){
      // Compute the relaxation factor
      TacsScalar r_factor_sens = 0.0;
      if (eps > 0.0){
        TacsScalar d = 1.0/(eps*(1.0-rho[0])+rho[0]);
        r_factor_sens = eps*d*d;
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor_sens = pow(rho[0],p-1)*p;
      }
      TacsScalar inner = alpha*r_factor_sens*maxtemp/props->Tmax[0];
      for ( int i = 0; i < nweights; i++ ){
        dvSens[index[i]] += N[i]*inner;
      }
    }
    else {
      TacsScalar r_factor_sens = 0.0;
      int j = mat_j+1;
      if (eps > 0.0){
	TacsScalar d = 1.0/(eps*(1.0-rho[j])+rho[j]);
	r_factor_sens = eps*d*d;
      }
      else {
	TacsScalar p = -1.0*eps;
	r_factor_sens = pow(rho[j],p-1)*p;
      }
      // Compute the inner product
      TacsScalar inner = alpha*(r_factor_sens*maxtemp/props->Tmax[mat_j]);
        
      for ( int i = 0; i < nweights; i++ ){
	dvSens[nvars*index[i] + j] += N[i]*inner;
      }
      // for ( int j = 1; j < nvars; j++ ){
      //   TacsScalar r_factor_sens = 0.0;
      //   if (eps > 0.0){
      //     TacsScalar d = 1.0/(eps*(1.0-x[j])+x[j]);
      //     r_factor_sens = eps*d*d;
      //   }
      //   else {
      //     TacsScalar p = -1.0*eps;
      //     r_factor_sens = pow(x[j],p-1)*p;
      //   }
      // 	// Compute the inner product
      //   TacsScalar inner = alpha*(r_factor_sens*maxtemp);
        
      //   for ( int i = 0; i < nweights; i++ ){
      //     dvSens[nvars*index[i] + j] += N[i]*inner;
      //   }
      // }
    }
  }
  else{
    if (nvars == 1){
      // Compute the relaxation factor
      TacsScalar r_factor_sens = 0.0;
      if (eps > 0.0){
        TacsScalar d = 1.0/(eps*(1.0-rho[0])+rho[0]);
        r_factor_sens = eps*d*d;
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor_sens = pow(rho[0],p-1)*p;
      }
      TacsScalar inner = alpha*r_factor_sens*maxtemp/props->Tmax[0];
      for ( int i = 0; i < nweights; i++ ){
        dvSens[weights[i].index] += weights[i].weight*inner;
      }
    }
    else {
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar r_factor_sens = 0.0;
        if (eps > 0.0){
          TacsScalar d = 1.0/(eps*(1.0-rho[j])+rho[j]);
          r_factor_sens = eps*d*d;
        }
        else {
          TacsScalar p = -1.0*eps;
          r_factor_sens = pow(rho[j],p-1)*p;
        }
	TacsScalar inner = alpha*(r_factor_sens*maxtemp/props->Tmax[j-1]);
        for ( int i = 0; i < nweights; i++ ){
          dvSens[nvars*weights[i].index + j] += weights[i].weight*inner;
        }
      }
    }
  }

}

void TMRCoupledThermoOctStiffness::maxtempStrainSens( const double pt[], 
						      const TacsScalar maxtemp,
						      TacsScalar sens[],
						      int mat_j ){
  const double eps = props->eps;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*i + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];
    }
    if (nvars == 1){
      // Compute the relaxation factor
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = rho[0]/(eps*(1.0-x[0])+x[0]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(x[0],p);
      }
      sens[0] = r_factor/props->Tmax[0];
    }
    else {
      sens[0] = 0.0;
      int j = mat_j+1;
      TacsScalar r_factor_sens = 0.0;
      if (eps > 0.0){
	TacsScalar d = 1.0/(eps*(1.0-rho[j])+rho[j]);
	r_factor_sens = eps*d*d;
      }
      else {
	TacsScalar p = -1.0*eps;
	r_factor_sens = pow(rho[j],p-1)*p;
      }
      sens[0] += r_factor_sens/props->Tmax[mat_j];
      
      // for ( int j = 1; j < nvars; j++ ){
      //   TacsScalar r_factor_sens = 0.0;
      //   if (eps > 0.0){
      //     TacsScalar d = 1.0/(eps*(1.0-rho[j])+rho[j]);
      //     r_factor_sens = eps*d*d;
      //   }
      //   else {
      //     TacsScalar p = -1.0*eps;
      //     r_factor_sens = pow(rho[j],p-1)*p;
      //   }
      // 	sens[0] += r_factor_sens;
      // }
      
    }
  }
  else {
    if (nvars == 1){
      // Compute the relaxation factor
      TacsScalar r_factor = 1.0;
      if (eps > 0.0){
        r_factor = rho[0]/(eps*(1.0-x[0])+x[0]);
      }
      else {
        TacsScalar p = -1.0*eps;
        r_factor = pow(x[0],p);
      }
      sens[0] = r_factor/props->Tmax[0];
    }
    else {
      sens[0] = 0.0;
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar r_factor_sens = 0.0;
        if (eps > 0.0){
          TacsScalar d = 1.0/(eps*(1.0-rho[j])+rho[j]);
          r_factor_sens = eps*d*d;
        }
        else {
          TacsScalar p = -1.0*eps;
          r_factor_sens = pow(rho[j],p-1)*p;
        }
	sens[0] += r_factor_sens/props->Tmax[j-1];
      }
    }
  }
}
