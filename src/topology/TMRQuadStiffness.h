#ifndef TMR_QUADRANT_STIFFNESS_H
#define TMR_QUADRANT_STIFFNESS_H

/*
  The following class defines a PlaneStressQuadStiffness base class for topology
  optimization using TACS.  
*/

#include "TMRBase.h"
#include "PlaneStressStiffness.h"

/*
  The TMRQuadStiffness class

  This defines the TMRQuadStiffness class which takes the weights from
  up to 4 adjacent vertices. This uses the RAMP method for
  penalization.
*/
class TMRQuadStiffness : public PlaneStressStiffness {
 public: 
  static const int MAX_NUM_WEIGHTS = 4;

  TMRQuadStiffness( TMRIndexWeight *_weights, int _nweights,
                    TacsScalar _density, TacsScalar E, 
                    TacsScalar _nu, double _q, double _eps=1e-3 );

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
    return rho; 
  }

 private:
  // The density and stiffness properties
  TacsScalar density;
  TacsScalar E, nu;

  // The RAMP penalization factor
  TacsScalar q;
  double eps;

  // The value of the design-dependent density
  TacsScalar rho;

  // The local density of the
  int nweights;
  TMRIndexWeight weights[MAX_NUM_WEIGHTS];
};

#endif // TMR_QUADRANT_STIFFNESS
