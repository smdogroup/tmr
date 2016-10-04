#ifndef TMR_OCTANT_STIFFNESS_H
#define TMR_OCTANT_STIFFNESS_H

/*
  The following class defines a SolidStiffness base class for topology
  optimization using TACS.  
*/

#include "TMRBase.h"
#include "SolidStiffness.h"

/*
  The TMROctStiffness class

  This defines the TMROctStiffness class which takes the weights from
  up to 8 adjacent vertices. This uses the RAMP method for
  penalization.
*/
class TMROctStiffness : public SolidStiffness {
 public: 
  static const int MAX_NUM_WEIGHTS = 8;

  TMROctStiffness( TMRIndexWeight *_weights, int _nweights,
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
  TacsScalar D, G, nu;

  // The RAMP penalization factor
  TacsScalar q;
  double eps;

  // The value of the design-dependent density
  TacsScalar rho;

  // The local density of the
  int nweights;
  TMRIndexWeight weights[MAX_NUM_WEIGHTS];
};

/*
  TMRLinearOctStiffness class

  The linearized stiffness class. This can be used to perform a
  sequential linearized convex optimization.
*/
class TMRLinearOctStiffness : public SolidStiffness {
 public: 
  static const int MAX_NUM_WEIGHTS = 8;

  TMRLinearOctStiffness( TMRIndexWeight *_weights, int _nweights,
                         TacsScalar _density, TacsScalar E, 
                         TacsScalar _nu, double _q, double _eps=1e-3);
  
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
  TacsScalar D, G, nu;

  // The RAMP penalization factor
  TacsScalar q;
  double eps;

  // The value of the design-dependent density
  TacsScalar rho;
  TacsScalar rho_const, rho_linear;

  // Set the lower bounds
  TacsScalar x_lb[MAX_NUM_WEIGHTS];

  // The local density of the
  int nweights;
  TMRIndexWeight weights[MAX_NUM_WEIGHTS];
};

#endif // TMR_OCTANT_STIFFNESS
