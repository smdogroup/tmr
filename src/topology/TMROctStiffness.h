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

/*
  The TMRStiffnessProperties class
*/
class TMRStiffnessProperties : public TMREntity {
 public:
  static const int MAX_NUM_MATERIALS = 5;
  
  // Set the stiffness properties
  TMRStiffnessProperties( int _nmats, TacsScalar *rho, 
                          TacsScalar *_E, TacsScalar *_nu ){
    if (_nmats > MAX_NUM_MATERIALS){
      _nmats = MAX_NUM_MATERIALS;
    }
    nmats = _nmats;
    for ( int i = 0; i < nmats; i++ ){
      density[i] = rho[i];
      E[i] = _E[i];
      nu[i] = _nu[i];
      G[i] = 0.5*E[i]/(1.0 + nu[i]);
      D[i] = E[i]/((1.0 + nu[i])*(1.0 - 2.0*nu[i]));
    }
  }
                   
  int nmats;
  TacsScalar density[MAX_NUM_MATERIALS];
  TacsScalar E[MAX_NUM_MATERIALS];
  TacsScalar D[MAX_NUM_MATERIALS];
  TacsScalar nu[MAX_NUM_MATERIALS];
  TacsScalar G[MAX_NUM_MATERIALS];
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
                   TMRStiffnessProperties *_props,
                   double _q, double _eps=1e-3 );
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

  // Return the density as the design variable
  // -----------------------------------------
  TacsScalar getDVOutputValue( int dvIndex, const double pt[] ){ 
    return rho[0]; 
  }

 private:
  // The stiffness properties
  TMRStiffnessProperties *props;

  // The RAMP penalization factor
  TacsScalar q;
  double eps;

  // The value of the design-dependent density
  int nvars;
  TacsScalar rho[MAX_NUM_MATERIALS+1];

  // The local density of the
  int nweights;
  TMRIndexWeight *weights;
};

/*

*/
/*
class TMROctStiff : public SolidStiffness {
 public:
  TMROctStiff( TMROctForest *_forest, int _elem, 
               TMRStiffnessProperties *_props,
               double _q, double _eps=1e-3 ){

  }
  
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

 private:
  // Keep track of the element/forest
  int elem;
  TMROctForest *forest;

  // The stiffness properties
  TMRStiffnessProperties *props;

  // The RAMP penalization factor
  TacsScalar q;
  double eps;

  // Compute the forest values
  TacsScalar *rho;
};
*/
/*
  TMRLinearOctStiffness class

  The linearized stiffness class. This can be used to perform a
  sequential linearized convex optimization.
*/
class TMRLinearOctStiffness : public SolidStiffness {
 public: 
  enum PenaltyType { SIMP, RAMP };
  static const int MAX_NUM_WEIGHTS = 8;

  TMRLinearOctStiffness( TMRIndexWeight *_weights, int _nweights,
                         TacsScalar _x_init,
                         TacsScalar _density, TacsScalar E, 
                         TacsScalar _nu, double _q, 
                         PenaltyType _type=RAMP, double _eps=1e-3);
  
  // Set the linearization coefficients
  // ----------------------------------
  void setLinearization( TacsScalar _q,
                         const TacsScalar dvs[], int numDVs );

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

  // Compute the derivative of the stress
  // ------------------------------------
  void calcStressDVProject( const double pt[],
                            const TacsScalar e[],
                            const TacsScalar px[],
                            int dvLen, TacsScalar s[] );

  // Evaluate the pointwise mass
  // ---------------------------
  void getPointwiseMass( const double pt[], TacsScalar mass[] );
  void addPointwiseMassDVSens( const double pt[], 
                               const TacsScalar alpha[],
                               TacsScalar dvSens[], int dvLen );

  // Return the density as the design variable
  // -----------------------------------------
  TacsScalar getDVOutputValue( int dvIndex, const double pt[] ){ 
    return rho; 
  }

 private:
  // The density and stiffness properties
  TacsScalar density;
  TacsScalar D, G, nu;

  // The RAMP penalization factor
  static PenaltyType penalty_type;
  TacsScalar q;
  double eps;

  // The value of the design-dependent density
  TacsScalar rho;
  TacsScalar rho_const, rho_linear;

  // The lower bounds and design variable values
  TacsScalar x_vals[MAX_NUM_WEIGHTS];
  TacsScalar x_lb[MAX_NUM_WEIGHTS];

  // The local density of the
  int nweights;
  TMRIndexWeight weights[MAX_NUM_WEIGHTS];
};

#endif // TMR_OCTANT_STIFFNESS_H
