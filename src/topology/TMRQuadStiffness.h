#ifndef TMR_QUADRANT_STIFFNESS_H
#define TMR_QUADRANT_STIFFNESS_H

/*
  The following class defines a SolidStiffness base class for topology
  optimization using TACS.  
*/

#include "TMRBase.h"
#include "PlaneStressStiffness.h"
#include "TACSConstitutive.h"
/*
  Stiffness properties object
*/
class TMRQuadStiffnessProperties {
 public:
  // Initialize the properties so that they have some value
  TMRQuadStiffnessProperties(){
    rho = 1.0;
    E = 70.0e9;
    nu = 0.3;
    q = 1.0;
    ys = 280.e3;
  }

  TacsScalar rho; // Material density
  TacsScalar E; // Young's modulus
  TacsScalar nu; // Poisson ratio
  TacsScalar q; // RAMP penalty parameter
  TacsScalar ys; // Yield stress
};

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
                    TacsScalar _density, TacsScalar _E, 
                    TacsScalar _nu, TacsScalar _ys,
                    double _q, double _eps=1e-3);
  
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
  void calculateThermalStress( const double pt[],
                               const TacsScalar e[], TacsScalar s[] );
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
                         TacsScalar sens[]);

  // Return the density as the design variable
  // -----------------------------------------
  TacsScalar getDVOutputValue( int dvIndex, const double pt[] ){ 
    return xw*density; 
  }
  
 private:
  // The density and stiffness properties
  double E, nu, ys, epsilon;
  TacsScalar rho, alpha, xw,density;

  // The RAMP penalization factor
  TacsScalar q;
  double eps;

  // The local density of the
  int nweights;
  TMRIndexWeight weights[MAX_NUM_WEIGHTS];
};

#endif // TMR_THERMO_QUADRANT_STIFFNESS_H
