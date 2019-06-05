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

#ifndef TMR_QUADRANT_STIFFNESS_H
#define TMR_QUADRANT_STIFFNESS_H

/*
  The following class defines a PlaneStressQuadStiffness base class for topology
  optimization using TACS.  
*/

#include "TMRBase.h"
#include "PlaneStressStiffness.h"
#include "YSlibrary.h"
/*
  The TMRStiffnessProperties class
*/
class TMRQuadStiffnessProperties : public TMREntity {
 public:
  static const int MAX_NUM_MATERIALS = 5;

  // Set the stiffness properties
  TMRQuadStiffnessProperties( int _nmats,
                              double _q, double _eps, double _k0,
                              double _beta, double _xoffset,
                              TacsScalar *rho, TacsScalar *_E, TacsScalar *_nu,
                              TacsScalar *_ys=NULL, TacsScalar *_aT=NULL,
                              TacsScalar *_kcond=NULL,TacsScalar *_Tmax=NULL,
			      double _qtemp=0.0,
                              double _qcond=0.0, int _use_project=0 ){
    if (_nmats > MAX_NUM_MATERIALS){
      _nmats = MAX_NUM_MATERIALS;
    }
    nmats = _nmats;
    q = _q;
    k0 = _k0;
    eps = _eps;
    beta = _beta;
    xoffset = _xoffset;
    use_project = _use_project;
    // Penalization for thermoelastic problem
    qtemp = _qtemp;
    qcond = _qcond;

    for ( int i = 0; i < nmats; i++ ){
      density[i] = rho[i];
      E[i] = _E[i];
      nu[i] = _nu[i];
      ys[i] = 1.0e8;
      G[i] = 0.5*E[i]/(1.0 + nu[i]);
      D[i] = E[i]/(1.0 - nu[i]*nu[i]);
      aT[i] = 0.0;
      kcond[i] = 0.0;
      if (_ys){
        ys[i] = _ys[i];
      }
      if (_aT){
        aT[i] = _aT[i];
      }
      if (_kcond){
        kcond[i] = _kcond[i];
      }
      if (_Tmax){
	Tmax[i] = _Tmax[i];
      }
    }
  }

  int nmats;   // Number of materials to use
  double q;    // RAMP penalization factor
  double eps;  // Stress relaxation parameter
  double k0;   // Small stiffness factor >= 0 ~ 1e-6
  double beta; // Parameter for the logistics function
  double xoffset;  // Offset parameter in the logistics function
  int use_project; // Flag to indicate if projection should be used (0, 1)
  TacsScalar density[MAX_NUM_MATERIALS]; // Material density
  TacsScalar E[MAX_NUM_MATERIALS]; // Young's modulus
  TacsScalar D[MAX_NUM_MATERIALS]; // Stiffness
  TacsScalar nu[MAX_NUM_MATERIALS]; // Poisson's ratio
  TacsScalar G[MAX_NUM_MATERIALS]; // Shear stiffness
  TacsScalar ys[MAX_NUM_MATERIALS]; // Yield stress
  TacsScalar aT[MAX_NUM_MATERIALS]; // Heat coefficient of expansion
  TacsScalar kcond[MAX_NUM_MATERIALS]; // Heat conductivity
  TacsScalar Tmax[MAX_NUM_MATERIALS]; // Maximum temperature
  TacsScalar qtemp; // RAMP penalty parameter for temperature
  TacsScalar qcond; // RAMP penalty parameter for conduction
};

/*
  The TMRQuadStiffness class

  This defines the TMRQuadStiffness class which takes the weights from
  up to 4 adjacent vertices. This uses the RAMP method for
  penalization.
*/
class TMRQuadStiffness : public PlaneStressStiffness {
 public: 
  static const int MAX_NUM_MATERIALS = 5;
  TMRQuadStiffness( TMRIndexWeight *_weights, int _nweights,
                    TMRQuadStiffnessProperties *_props );
  ~TMRQuadStiffness();

  // Set the design variable values in the object
  // --------------------------------------------
  void setDesignVars( const TacsScalar x[], int numDVs );
  void getDesignVars( TacsScalar x[], int numDVs );
  void getDesignVarRange( TacsScalar lb[], TacsScalar ub[], int numDVs );

  // Compute the stress
  // ------------------
  void calculateStress( const double pt[],
                        const TacsScalar e[], TacsScalar s[] );
  void addStressDVSens( const double pt[], const TacsScalar strain[], 
                        TacsScalar alpha, const TacsScalar psi[], 
                        TacsScalar dvSens[], int dvLen );

  // Evaluate the pointwise mass
  // ---------------------------
  void getPointwiseMass( const double pt[], 
                         TacsScalar mass[] );
  void addPointwiseMassDVSens( const double pt[], 
                               const TacsScalar alpha[],
                               TacsScalar dvSens[], int dvLen );

  // Return the density as the design variable
  // -----------------------------------------
  TacsScalar getDVOutputValue( int dvIndex, const double pt[] ){ 
    return rho[0]; 
  }
  // Return the failure
  // -------------------
  void failure( const double pt[],
                const TacsScalar strain[],
                TacsScalar * fail );

  void addFailureDVSens( const double pt[],
                         const TacsScalar strain[],
                         TacsScalar alpha,
                         TacsScalar dvSens[], int dvLen );

  void failureStrainSens(const double pt[],
                         const TacsScalar strain[],
                         TacsScalar sens[] );

 private:
  // The stiffness properties
  TMRQuadStiffnessProperties *props;

  // The value of the design-dependent density
  int nvars;
  TacsScalar x[MAX_NUM_MATERIALS+1];
  TacsScalar rho[MAX_NUM_MATERIALS+1];


  // The local density of the
  int nweights;
  TMRIndexWeight *weights;
};

#endif // TMR_QUADRANT_STIFFNESS
