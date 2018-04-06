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

#ifndef TMR_OCTANT_STIFFNESS_H
#define TMR_OCTANT_STIFFNESS_H

/*
  The following class defines a SolidStiffness base class for topology
  optimization using TACS.  
*/

#include "TMRBase.h"
#include "SolidStiffness.h"
#include "YSlibrary.h"

/*
  The TMRStiffnessProperties class
*/
class TMRStiffnessProperties : public TMREntity {
 public:
  static const int MAX_NUM_MATERIALS = 5;
  
  // Set the stiffness properties
  TMRStiffnessProperties( int _nmats, 
                          double _q, double _eps, double _k0,
                          double _beta, double _xoffset,
                          TacsScalar *rho, TacsScalar *_E, TacsScalar *_nu,
                          TacsScalar *_ys=NULL ){
    if (_nmats > MAX_NUM_MATERIALS){
      _nmats = MAX_NUM_MATERIALS;
    }
    nmats = _nmats;
    q = _q;
    k0 = _k0;
    eps = _eps;
    for ( int i = 0; i < nmats; i++ ){
      density[i] = rho[i];
      E[i] = _E[i];
      nu[i] = _nu[i];
      ys[i] = 1.0e8;
      G[i] = 0.5*E[i]/(1.0 + nu[i]);
      D[i] = E[i]/((1.0 + nu[i])*(1.0 - 2.0*nu[i]));
      if (_ys){
        ys[i] = _ys[i];
      }
    }
  }
                   
  int nmats;   // Number of materials to use
  double q;    // RAMP penalization factor
  double eps;  // Stress relaxation parameter
  double k0;   // Small stiffness factor >= 0 ~ 1e-6
  double beta; // Parameter for the logistics function
  double xoffset;  // Offset parameter in the logistics function
  TacsScalar density[MAX_NUM_MATERIALS]; // Material density
  TacsScalar E[MAX_NUM_MATERIALS]; // Young's modulus
  TacsScalar D[MAX_NUM_MATERIALS]; // Stiffness
  TacsScalar nu[MAX_NUM_MATERIALS]; // Poisson's ratio
  TacsScalar G[MAX_NUM_MATERIALS]; // Shear stiffness
  TacsScalar ys[MAX_NUM_MATERIALS]; // Yield stress
};

/*
  The TMROctStiffness class

  This defines the TMROctStiffness class which takes the weights from
  up to 8 adjacent vertices. This class uses the RAMP method for
  penalization.
*/
class TMROctStiffness : public SolidStiffness {
 public:
  static const int MAX_NUM_MATERIALS = 5;
  static const int MAX_NUM_WEIGHTS = 8;

  TMROctStiffness( TMRIndexWeight *_weights, int _nweights,
                   TMRStiffnessProperties *_props );
  ~TMROctStiffness();

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
  
  // Return the density as the design variable
  // -----------------------------------------
  TacsScalar getDVOutputValue( int dvIndex, const double pt[] ){ 
    return rho[0]; 
  }

 private:
  // The stiffness properties
  TMRStiffnessProperties *props;

  // The value of the design-dependent density
  int nvars;
  TacsScalar x[MAX_NUM_MATERIALS+1];
  TacsScalar rho[MAX_NUM_MATERIALS+1];

  // The local density of the
  int nweights;
  TMRIndexWeight *weights;
};

#endif // TMR_OCTANT_STIFFNESS_H
