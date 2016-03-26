#ifndef TMR_OCTANT_STIFFNESS_H
#define TMR_OCTANT_STIFFNESS_H

#include "TMROctree.h"
#include "TMROctant.h"
#include "SolidStiffness.h"

class OctStiffness : public SolidStiffness {
 public: 
  OctStiffness( TMROctree *filter, TMROctant *oct,
		TacsScalar _density,
		TacsScalar _E, TacsScalar _nu,
		double _q );

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
  static const int FILTER_ORDER = 2;
  static const int MAX_NUM_WEIGHTS = 
    (FILTER_ORDER*FILTER_ORDER*FILTER_ORDER)*(FILTER_ORDER*FILTER_ORDER);

  // The density and stiffness properties
  TacsScalar density;
  TacsScalar D, G, nu;

  // The RAMP penalization factor
  TacsScalar q;

  // The value of the design-dependent density
  TacsScalar rho;

  // The local density of the
  int nweights;
  TMRIndexWeight weights[MAX_NUM_WEIGHTS];
};

#endif // TMR_OCTANT_STIFFNESS
