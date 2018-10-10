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

#include "TMRCoupledThermoQuadStiffness.h"
#include "FElibrary.h"
#include "YSlibrary.h"

/*
  Create the quadtree stiffness object based on an interpolation from
  the filter variables
*/
TMRCoupledThermoQuadStiffness::TMRCoupledThermoQuadStiffness( TMRIndexWeight *_weights,
                                                              int _nweights,
                                                              TMRQuadStiffnessProperties *_props,
 TMRQuadForest *_filter ){
  

  // Record the density, Poisson ratio, D and the shear modulus
  props = _props;
  props->incref();
  
  // Set the weights/local indices
  nweights = _nweights;
  weights = new TMRIndexWeight[ nweights ];
  memcpy(weights, _weights, nweights*sizeof(TMRIndexWeight));
    
  // Copy over the filter if provided
  filter = NULL;
  if (_filter){
    filter = _filter;
    filter->incref();
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
  dv = NULL;
}
/*
  Decref the props
*/
TMRCoupledThermoQuadStiffness::~TMRCoupledThermoQuadStiffness(){
  props->decref();
  delete [] weights;
  if (filter){
    filter->decref();
  }
  if (dv){
    delete [] dv;
  }
}
   
/*
  Loop over the design variable inputs and compute the local value 
  of the density
*/
void TMRCoupledThermoQuadStiffness::setDesignVars( const TacsScalar xdv[], 
                                                   int numDVs ){
  const double beta = props->beta;
  const double xoffset = props->xoffset;
  const int use_project = props->use_project;

  if (!dv){
    dv = new TacsScalar[ numDVs ];
  }
  memcpy(dv, xdv, numDVs*sizeof(TacsScalar));
  
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
void TMRCoupledThermoQuadStiffness::getDesignVars( TacsScalar xdv[], 
                                                   int numDVs ){
  dv = new TacsScalar[ numDVs ];
  memset(dv, 0.0, numDVs*sizeof(TacsScalar));

  if (nvars == 1){
    for ( int i = 0; i < nweights; i++ ){
      if (weights[i].index >= 0 && weights[i].index < numDVs){
        xdv[weights[i].index] = 0.95;
      }
    }
  }
  else {
    for ( int j = 0; j < nvars; j++ ){
      TacsScalar value = 0.25;
      if (j >= 1){
        value = 0.75/(nvars-1);
      }
      for ( int i = 0; i < nweights; i++ ){
        if (weights[i].index >= 0 && weights[i].index < numDVs){
          xdv[nvars*weights[i].index + j] = value;
        }
      }
    }
  }
  memcpy(dv, xdv, numDVs*sizeof(TacsScalar));
}

/*
  Retrieve the range of the design variable values
*/
void TMRCoupledThermoQuadStiffness::getDesignVarRange( TacsScalar lb[], 
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
  Compute the density of the material as a linear function of the
  element-wise density
*/
void TMRCoupledThermoQuadStiffness::getPointwiseMass( const double pt[], 
                                                      TacsScalar mass[] ){

  // Get the density value at isoparametric point pt
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        x[j] += N[i]*dv[nvars*weights[i].index + j];
        //x[j] += N[i]*weights[i].weight*dv[nvars*weights[i].index + j];
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
  Compute the derivative of the pointwise mass as function of the
  design variable values
*/
void TMRCoupledThermoQuadStiffness::addPointwiseMassDVSens( const double pt[], 
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
        dvSens[weights[i].index] += (scale*N[i]);
        //dvSens[weights[i].index] += (scale*weights[i].weight*N[i]);
      }
    }
    else {
      for ( int j = 1; j < nvars; j++ ){
        TacsScalar scale = props->density[j-1]*alpha[0];
        for ( int i = 0; i < nweights; i++ ){
          dvSens[nvars*weights[i].index + j] += (scale*N[i]);
          //dvSens[nvars*weights[i].index + j] += (scale*weights[i].weight*N[i]);
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
void TMRCoupledThermoQuadStiffness::calculateStress( const double pt[],
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
        x[j] += N[i]*dv[nvars*weights[i].index + j];
        //x[j] += N[i]*weights[i].weight*dv[nvars*weights[i].index + j];
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
    s[0] = Dp*(e[0]+nu*e[1]);
    s[1] = Dp*(nu*e[0]+e[1]);
    s[2] = Gp*e[2];
  }
  else {
    // Compute the penalized stiffness
    s[0] = s[1] = s[2] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      TacsScalar penalty = rho[j]/(1.0+qs*(1.0-rho[j]));
      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];

      // Add the penalized value
      TacsScalar Dp = (penalty+k0)*D;
      TacsScalar Gp = (penalty+k0)*G;

      s[0] += Dp*(e[0]+nu*e[1]);
      s[1] += Dp*(nu*e[0]+e[1]);
      s[2] += Gp*e[2];
    }
  }
}
/*
  Compute the matrix vector product of L*e
*/
void TMRCoupledThermoQuadStiffness::calculateThermal( const double pt[],
                                                      const TacsScalar e[], 
                                                      TacsScalar s[] ){
  const double k0 = props->k0;
  const double qtemp = props->qtemp;
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        //x[j] += N[i]*weights[i].weight*dv[nvars*weights[i].index + j];
        x[j] += N[i]*dv[nvars*weights[i].index + j];
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
    s[0] = Dp*(e[0]+nu*e[1]);
    s[1] = Dp*(nu*e[0]+e[1]);
    s[2] = Gp*e[2];
  }
  else {
    // aT is added here
    s[0] = s[1] = s[2] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      TacsScalar penalty = rho[j]/(1.0+qtemp*(1.0-rho[j]));
      // Extract the properties
      TacsScalar nu = props->nu[j-1];
      TacsScalar D = props->D[j-1];
      TacsScalar G = props->G[j-1];
      TacsScalar aT = props->aT[j-1];
      // Add the penalized value
      TacsScalar Dp = (penalty+k0)*D*aT;
      TacsScalar Gp = (penalty+k0)*G;

      s[0] += Dp*(e[0]+nu*e[1]);
      s[1] += Dp*(nu*e[0]+e[1]);
      s[2] += Gp*e[2];
    }
  }
}
/*
  Compute the conduction terms D_{th}*u
*/
void TMRCoupledThermoQuadStiffness::calculateConduction( const double pt[],
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
        x[j] += N[i]*dv[nvars*weights[i].index + j];
        //x[j] += N[i]*weights[i].weight*dv[nvars*weights[i].index + j];
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
  }
  else {
    s[0] = s[1] = 0.0;
    for ( int j = 1; j < nvars; j++ ){
      TacsScalar penalty = rho[j]/(1.0+qcond*(1.0-rho[j]));
      // Extract the properties
      TacsScalar kcond = props->kcond[j-1];
      // Add the penalized value
      TacsScalar Cp = (penalty+k0)*kcond;
      s[0] += Cp*e[0];
      s[1] += Cp*e[1];
    }
  }  
}

/*
  Add the derivative of the product of the stress with the vector psi
  times alpha to the design variable array dvSens
*/
void TMRCoupledThermoQuadStiffness::addStressDVSens( const double pt[], 
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
        x[j] += N[i]*dv[nvars*weights[i].index + j];
        //x[j] += N[i]*weights[i].weight*dv[nvars*weights[i].index + j];
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
      TacsScalar s[3];
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];

      // Compute the penalized value
      TacsScalar Dp = alpha*penalty*D;
      TacsScalar Gp = alpha*penalty*G;
      // Compute the product dD/dx*(B*u)
      s[0] = Dp*(e[0]+nu*e[1]);
      s[1] = Dp*(nu*e[0]+e[1]);
      s[2] = Gp*e[2];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);

      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[weights[i].index] += N[i]*product;
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
        TacsScalar s[3];
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];

        // Add the result to the derivative
        TacsScalar Dp = alpha*penalty*D;
        TacsScalar Gp = alpha*penalty*G;
        // Compute the product dD/dx*(B*u)
        s[0] = Dp*(e[0]+nu*e[1]);
        s[1] = Dp*(nu*e[0]+e[1]);
        s[2] = Gp*e[2];

        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);

        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*weights[i].index + j] += N[i]*product;
        }
      }
    }
  }
  else{
    if (nvars == 1){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (q + 1.0)/((1.0 + q*(1.0 - rho[0]))*(1.0 + q*(1.0 - rho[0])));
      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[3];
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];

      // Compute the penalized value
      TacsScalar Dp = alpha*penalty*D;
      TacsScalar Gp = alpha*penalty*G;
      // Compute the product dD/dx*(B*u)
      s[0] = Dp*(e[0]+nu*e[1]);
      s[1] = Dp*(nu*e[0]+e[1]);
      s[2] = Gp*e[2];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);

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
        TacsScalar s[3];
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];

        // Add the result to the derivative
        TacsScalar Dp = alpha*penalty*D;
        TacsScalar Gp = alpha*penalty*G;
        // Compute the product dD/dx*(B*u)
        s[0] = Dp*(e[0]+nu*e[1]);
        s[1] = Dp*(nu*e[0]+e[1]);
        s[2] = Gp*e[2];

        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);

        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*weights[i].index + j] += weights[i].weight*product;
        }
      }
    }
  }
}

void TMRCoupledThermoQuadStiffness::addThermalDVSens( const double pt[], 
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
        //x[j] += N[i]*weights[i].weight*dv[nvars*weights[i].index + j];
        x[j] += N[i]*dv[nvars*weights[i].index + j];
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
      TacsScalar s[3];
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      TacsScalar aT = props->aT[0];
      // Compute the penalized value
      TacsScalar Dp = alpha*penalty*D*aT;
      TacsScalar Gp = alpha*penalty*G;
      // Compute the product dD/dx*(B*u)
      s[0] = Dp*(e[0]+nu*e[1]);
      s[1] = Dp*(nu*e[0]+e[1]);
      s[2] = Gp*e[2];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);

      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[weights[i].index] += N[i]*product;
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
        TacsScalar s[3];
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];
        TacsScalar aT = props->aT[j-1];
        // Add the result to the derivative
        TacsScalar Dp = alpha*penalty*D*aT;
        TacsScalar Gp = alpha*penalty*G;
        // Compute the product dD/dx*(B*u)
        s[0] = Dp*(e[0]+nu*e[1]);
        s[1] = Dp*(nu*e[0]+e[1]);
        s[2] = Gp*e[2];

        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
        // Compute the term psi^{T}*B^{T}*dD/dx*B*u
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*weights[i].index + j] += N[i]*product;
        }
      }
    }
  }
  else{
    if (nvars == 1){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (qtemp + 1.0)/((1.0 + qtemp*(1.0 - rho[0]))*(1.0 + qtemp*(1.0 - rho[0])));

      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[3];
      // Extract the properties
      TacsScalar nu = props->nu[0];
      TacsScalar D = props->D[0];
      TacsScalar G = props->G[0];
      TacsScalar aT = props->aT[0];
      // Compute the penalized value
      TacsScalar Dp = alpha*penalty*D*aT;
      TacsScalar Gp = alpha*penalty*G;
      // Compute the product dD/dx*(B*u)
      s[0] = Dp*(e[0]+nu*e[1]);
      s[1] = Dp*(nu*e[0]+e[1]);
      s[2] = Gp*e[2];

      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);

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
        TacsScalar s[3];
        // Extract the properties
        TacsScalar nu = props->nu[j-1];
        TacsScalar D = props->D[j-1];
        TacsScalar G = props->G[j-1];
        TacsScalar aT = props->aT[j-1];
        // Add the result to the derivative
        TacsScalar Dp = alpha*penalty*D*aT;
        TacsScalar Gp = alpha*penalty*G;
        // Compute the product dD/dx*(B*u)
        s[0] = Dp*(e[0]+nu*e[1]);
        s[1] = Dp*(nu*e[0]+e[1]);
        s[2] = Gp*e[2];

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
  Add derivative w.r.t. the design variables due to the conduction term
*/
void TMRCoupledThermoQuadStiffness::addConductionDVSens( const double pt[], 
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
        x[j] += N[i]*dv[nvars*weights[i].index + j];
      }
      // Apply the projection to obtain the projected value
      // of the density
      rho[j] = x[j];      
    }
    if (nvars == 1){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (qcond + 1.0)/((1.0 + qcond*(1.0 - rho[0]))*(1.0 + qcond*(1.0 - rho[0])));

      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[2];
      // Compute the product dD/dx*(B*u)
      TacsScalar kcond = props->kcond[0];
      TacsScalar Cp = alpha*kcond*penalty;
      s[0] = Cp*e[0];
      s[1] = Cp*e[1];
    
      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1]);
      // Compute the term psi^{T}*B^{T}*dD/dx*B*u
      for ( int i = 0; i < nweights; i++ ){
        fdvSens[weights[i].index] += N[i]*product;
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
        TacsScalar s[2];
        s[0] = Cp*e[0];
        s[1] = Cp*e[1];
        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1]);
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*weights[i].index + j] += N[i]*product;
        }
      }
    }
  }
  else {
    if (nvars == 1){
      // Compute the derivative of the penalization with respect to
      // the projected density
      TacsScalar penalty =
        (qcond + 1.0)/((1.0 + qcond*(1.0 - rho[0]))*(1.0 + qcond*(1.0 - rho[0])));

      // Add the derivative of the projection
      if (use_project){
        penalty *= beta*exp(-beta*(x[0] - xoffset))*rho[0]*rho[0];
      }
      TacsScalar s[2];
      // Compute the product dD/dx*(B*u)
      TacsScalar kcond = props->kcond[0];
      TacsScalar Cp = alpha*kcond*penalty;
      s[0] = Cp*e[0];
      s[1] = Cp*e[1];
    
      TacsScalar product = (s[0]*psi[0] + s[1]*psi[1]);
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
          (qcond + 1.0)/((1.0 + qcond*(1.0 - rho[j]))*(1.0 + qcond*(1.0 - rho[j])));

        // Add the derivative of the projection
        if (use_project){
          penalty *= beta*exp(-beta*(x[j] - xoffset))*rho[j]*rho[j];
        }
      
        // Extract the properties
        TacsScalar kcond = props->kcond[j-1];
        // Add the penalized value
        TacsScalar Cp = alpha*penalty*kcond;
        TacsScalar s[2];
        s[0] = Cp*e[0];
        s[1] = Cp*e[1];
        TacsScalar product = (s[0]*psi[0] + s[1]*psi[1]);
        for ( int i = 0; i < nweights; i++ ){
          fdvSens[nvars*weights[i].index + j] += weights[i].weight*product;
        }
      }
    }
  }
}

/*
  Returns the unpenalized term alpha
*/
TacsScalar TMRCoupledThermoQuadStiffness::getEffThermalAlpha( int vars_j ){
  if (vars_j > 0){
    // Extract the properties
    TacsScalar aT = props->aT[vars_j-1];    
    return aT;
  }
  else {
    return props->aT[0];
  }
}

/*
  Returns the unpenalized term Tref
*/
TacsScalar TMRCoupledThermoQuadStiffness::getReferenceTemperature(){
  return Tref;
}


// Evaluate the failure criteria
// The strain vector is B*u, need to account for dT here
void TMRCoupledThermoQuadStiffness::failure( const double pt[], 
                                             const TacsScalar T[],
                                             const TacsScalar e[],
                                             TacsScalar *fail ){
  // const double eps = props->eps;

  // if (nvars == 1){
  //   // Use the von Mises failure criterion
  //   // Compute the relaxation factor
  //   TacsScalar r_factor = 1.0;

  //   if (eps > 0.0){
  //     r_factor = rho[0]/(eps*(1.0-rho[0])+rho[0]);
  //   }
  //   else {
  //     TacsScalar p = -1.0*eps;
  //     r_factor = pow(rho[0],p);
  //   }
  //   // Extract the properties
  //   TacsScalar nu = props->nu[0];
  //   TacsScalar D = props->D[0];
  //   TacsScalar G = props->G[0];
  //   TacsScalar aT = props->aT[0];
  //   TacsScalar ys = props->ys[0];
  //   TacsScalar s[3], eff_e[3];

  //   eff_e[0] = e[0]-aT*T[0];
  //   eff_e[1] = e[1]-aT*T[0];
  //   eff_e[2] = e[2]*1.0;

  //   // Compute the resulting stress
  //   s[0] = D*(eff_e[0]+nu*eff_e[1]);
  //   s[1] = D*(nu*eff_e[0]+eff_e[1]);
  //   s[2] = G*eff_e[2];
    
  //   *fail = r_factor*VonMisesFailurePlaneStress(s,ys);
  // }
  // else {    
  //   *fail = 0.0;    
  //   for ( int j = 1; j < nvars; j++ ){
  //     TacsScalar r_factor = 1.0;
  //     if (eps > 0.0){
  //       r_factor = rho[j]/(eps*(1.0-rho[j])+rho[j]);
  //     }
  //     else {
  //       TacsScalar p = -1.0*eps;
  //       r_factor = pow(rho[j],p);
  //     }
  //     // Extract the properties
  //     TacsScalar nu = props->nu[j-1];
  //     TacsScalar D = props->D[j-1];
  //     TacsScalar G = props->G[j-1];
  //     TacsScalar ys = props->ys[j-1];
  //     TacsScalar aT = props->aT[j-1];
  //     // Compute effective strain = B*u-alpha*dT
  //     TacsScalar s[3], eff_e[3];
  //     eff_e[0] = e[0]-aT*T[0];
  //     eff_e[1] = e[1]-aT*T[0];
  //     eff_e[2] = e[2]*1.0;
      
  //     // Compute the contribution of failure from each material
  //     s[0] = D*(eff_e[0]+nu*eff_e[1]);
  //     s[1] = D*(nu*eff_e[0]+eff_e[1]);
  //     s[2] = G*eff_e[2];
      
  //     *fail += r_factor*VonMisesFailurePlaneStress(s,ys);      
  //   }
  // }  
}
// Evaluate the failure criteria w.r.t. design variables
void TMRCoupledThermoQuadStiffness::addFailureDVSens( const double pt[], 
                                                      const TacsScalar T[],
                                                      const TacsScalar e[],
                                                      TacsScalar alpha,
                                                      TacsScalar dvSens[], 
                                                      int dvLen ){
  // const double eps = props->eps;
  // if (nvars == 1){
  //   // Compute the relaxation factor
  //   TacsScalar r_factor_sens = 0.0;
  //   if (eps > 0.0){
  //     TacsScalar d = 1.0/(eps*(1.0-rho[0])+rho[0]);
  //     r_factor_sens = eps*d*d;
  //   }
  //   else {
  //     TacsScalar p = -1.0*eps;
  //     r_factor_sens = pow(rho[0],p-1)*p;
  //   }
  //   // Extract the properties
  //   TacsScalar nu = props->nu[0];
  //   TacsScalar D = props->D[0];
  //   TacsScalar G = props->G[0];
  //   TacsScalar ys = props->ys[0];
  //   TacsScalar aT = props->aT[0];
  //   TacsScalar s[3], eff_e[3];

  //   eff_e[0] = e[0]-aT*T[0];
  //   eff_e[1] = e[1]-aT*T[0];
  //   eff_e[2] = e[2]*1.0;
  //   // Compute the resulting stress
  //   s[0] = D*(eff_e[0]+nu*eff_e[1]);
  //   s[1] = D*(nu*eff_e[0]+eff_e[1]);
  //   s[2] = G*eff_e[2];
  //   TacsScalar fail = VonMisesFailurePlaneStress(s,ys);
  //   TacsScalar inner = alpha*r_factor_sens*fail;
  //   for ( int i = 0; i < nweights; i++ ){
  //     dvSens[weights[i].index] += weights[i].weight*inner;
  //   }
  // }
  // else {
  //   for ( int j = 1; j < nvars; j++ ){
  //     TacsScalar r_factor_sens = 0.0;
  //     if (eps > 0.0){
  //       TacsScalar d = 1.0/(eps*(1.0-rho[j])+rho[j]);
  //       r_factor_sens = eps*d*d;
  //     }
  //     else {
  //       TacsScalar p = -1.0*eps;
  //       r_factor_sens = pow(rho[j],p-1)*p;
  //     }
  //     // Extract the properties
  //     TacsScalar nu = props->nu[j-1];
  //     TacsScalar D = props->D[j-1];
  //     TacsScalar G = props->G[j-1];
  //     TacsScalar ys = props->ys[j-1];
  //     TacsScalar aT = props->aT[j-1];
      
  //     TacsScalar s[3], eff_e[3];
  //     eff_e[0] = e[0]-aT*T[0];
  //     eff_e[1] = e[1]-aT*T[0];
  //     eff_e[2] = e[2]*1.0;
      
  //     // Compute the contribution of failure from each material
  //     s[0] = D*(eff_e[0]+nu*eff_e[1]);
  //     s[1] = D*(nu*eff_e[0]+eff_e[1]);
  //     s[2] = G*eff_e[2];
     
  //     // For the derivative dr/dx*sigma_vm/ys
  //     TacsScalar fail = VonMisesFailurePlaneStress(s,ys);
          
  //     // Compute the inner product
  //     TacsScalar inner = alpha*(r_factor_sens*fail);
      
  //     for ( int i = 0; i < nweights; i++ ){
  //       dvSens[nvars*weights[i].index + j] += weights[i].weight*inner;
  //     }
  //   }
  // }
}
void TMRCoupledThermoQuadStiffness::failureStrainSens( const double pt[],
                                                            const TacsScalar T[],
                                                            const TacsScalar e[],
                                                            TacsScalar sens[],
                                                            int vars_j ){
  // const double eps = props->eps;
  // if (nvars == 1){
  //   TacsScalar s[3], ps_sens[3], eff_e[3];
  //   // Use the von Mises failure criterion
  //   // Compute the relaxation factor
  //   TacsScalar r_factor = 1.0;
  //   if (eps > 0.0){
  //     r_factor = rho[0]/(eps*(1.0-rho[0])+rho[0]);
  //   }
  //   else {
  //     TacsScalar p = -1.0*eps;
  //     r_factor = pow(rho[0],p);
  //   }
  //   // Extract the properties
  //   TacsScalar nu = props->nu[0];
  //   TacsScalar D = props->D[0];
  //   TacsScalar G = props->G[0];
  //   TacsScalar ys = props->ys[0];
  //   TacsScalar aT = props->aT[0];
  //   // Compute the effective strain
  //   eff_e[0] = e[0]-aT*T[0];
  //   eff_e[1] = e[1]-aT*T[0];
  //   eff_e[2] = e[2]*1.0;

  //   // Compute the resulting stress
  //   s[0] = D*(eff_e[0]+nu*eff_e[1]);
  //   s[1] = D*(nu*eff_e[0]+eff_e[1]);
  //   s[2] = G*eff_e[2];
    
  //   VonMisesFailurePlaneStressSens(ps_sens, s, ys);
  //   sens[0] = r_factor*D*(ps_sens[0]+nu*ps_sens[1]);
  //   sens[1] = r_factor*D*(nu*ps_sens[0]+ps_sens[1]);
  //   sens[2] = r_factor*G*ps_sens[2];
  // }
  // else {
  //   int j = vars_j*1;
  //   TacsScalar r_factor = 1.0;
  //   if (eps > 0.0){
  //     r_factor = rho[j]/(eps*(1.0-rho[j])+rho[j]);
  //   }
  //   else {
  //     TacsScalar p = -1.0*eps;
  //     r_factor = pow(rho[j],p);
  //   }

  //   TacsScalar s[3], ps_sens[3], eff_e[3];
  //   // Extract the properties
  //   TacsScalar nu = props->nu[j-1];
  //   TacsScalar D = props->D[j-1];
  //   TacsScalar G = props->G[j-1];
  //   TacsScalar ys = props->ys[j-1];
  //   TacsScalar aT = props->aT[j-1];
  //   // Compute the effective strain B*u-a*dT
  //   eff_e[0] = e[0]-aT*T[0];
  //   eff_e[1] = e[1]-aT*T[0];
  //   eff_e[2] = e[2]*1.0;
  //   // Compute the contribution of failure from each material
  //   s[0] = D*(eff_e[0]+nu*eff_e[1]);
  //   s[1] = D*(nu*eff_e[0]+eff_e[1]);
  //   s[2] = G*eff_e[2];
      
  //   VonMisesFailurePlaneStressSens(ps_sens, s, ys);
  //   sens[0] = r_factor*D*(ps_sens[0]+nu*ps_sens[1]);
  //   sens[1] = r_factor*D*(nu*ps_sens[0]+ps_sens[1]);
  //   sens[2] = r_factor*G*ps_sens[2];
  // }
}
TacsScalar TMRCoupledThermoQuadStiffness::getDVOutputValue( int dvIndex, 
                                                            const double pt[] ){
  // Compute the filter shape functions
  if (filter && dv){
    double N[nweights];
    filter->evalInterp(pt, N);
    for ( int j = 0; j < nvars; j++ ){
      x[j] = 0.0;
      for ( int i = 0; i < nweights; i++ ){
        printf("N[%d]: %f weights: %f %d %d\n",i, N[i],
               weights[i].weight, weights[i].index, weights[i].mask);
        x[j] += N[i]*dv[nvars*weights[i].index + j];
        //x[j] += N[i]*weights[i].weight*dv[nvars*weights[i].index + j];
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
  // if (filter && dv){
  //   TMRIndexWeight *local_weights = new TMRIndexWeight[ nweights ];
  //   computeLocalWeights(pt, local_weights);
  //   double N[nweights];
  //   filter->evalInterp(pt, N);
  //   for ( int j = 0; j < nvars; j++ ){
  //     x[j] = 0.0;
  //     for ( int i = 0; i < nweights; i++ ){
  //       printf("N[%d]: %f weights: %f %d, new: %f %d\n",i, N[i],
  //              weights[i].weight, weights[i].index, local_weights[i].weight,
  //              local_weights[i].index);
  //       x[j] += N[i]*dv[nvars*weights[i].index + j];
  //       //x[j] += N[i]*weights[i].weight*dv[nvars*weights[i].index + j];
  //     }
  //     // Apply the projection to obtain the projected value
  //     // of the density
  //     rho[j] = x[j];
  //   }
  // }
  // if (dvIndex == 0){
  //   if (nvars > 1){
  //     return rho[1];
  //   }
  //   return rho[0];
  // }
  else if (dvIndex == 1){
    if (nvars > 1){
      return rho[2];
    }
    return rho[0];
  }
  return 0;
}

// void TMRCoupledThermoQuadStiffness::computeLocalWeights( const double pt[],
//                                                          TMRIndexWeight *local_weights ){
//   double N[nweights];
//   filter->evalInterp(pt, N);

//   // Get the dependent node information for this mesh
//   const int *dep_ptr, *dep_conn;
//   const double *dep_weights;
//   filter->getDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);

//   // Get the mesh order
//   const int order = filter->getMeshOrder();
  
//   // Get the connectivity
//   const int *conn;
//   filter->getNodeConn(&conn);
//   const int *c = &conn[quad->tag*order*order];

//   // Loop over the adjacent nodes within the filter
//   int local_nweights = 0;
//   for ( int jj = 0; jj < order; jj++ ){
//     for ( int ii = 0; ii < order; ii++ ){
//       // Set the weights
//       int offset = ii + jj*order;
//       double weight = N[offset];

//       // Get the tag number
//       if (c[offset] >= 0){
//         local_weights[local_nweights].index = c[offset];
//         local_weights[local_nweights].weight = weight;
//         local_nweights++;
//       }
//       else {
//         int node = -c[offset]-1;
//         for ( int jp = dep_ptr[node]; jp < dep_ptr[node+1]; jp++ ){
//           local_weights[local_nweights].index = dep_conn[jp];
//           local_weights[local_nweights].weight = weight*dep_weights[jp];
//           local_nweights++;
//         }
//       }
//     }
//   }
//   local_nweights = TMRIndexWeight::uniqueSort(local_weights, local_nweights);

// }

