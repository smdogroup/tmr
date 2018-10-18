#ifndef TMR_COUPLED_THERMO_QUADRANT_STIFFNESS_H
#define TMR_COUPLED_THERMO_QUADRANT_STIFFNESS_H

/*
  The following class defines a coupled thermoelastic base class for topology
  optimization using TACS with multi-material.  
*/

#include "TMRBase.h"
#include "CoupledThermoPlaneStressStiffness.h"
#include "TACSConstitutive.h"
#include "TMRQuadStiffness.h"
#include "TMRQuadForest.h"
/*
  The TMRCoupledThermoQuadStiffness class

  This defines the TMRCoupledThermoQuadStiffness class which takes the weights from
  up to 4 adjacent vertices. This uses the RAMP method for
  penalization for TSC and x.
*/
class TMRCoupledThermoQuadStiffness : public CoupledThermoPlaneStressStiffness {
 public: 
  static const int MAX_NUM_MATERIALS = 5;
  TMRCoupledThermoQuadStiffness( TMRIndexWeight *_weights, 
                                 int _nweights,
                                 TMRQuadStiffnessProperties *_props,
                                 TMRQuadForest *_filter=NULL,
                                 int *_index=NULL );
  ~TMRCoupledThermoQuadStiffness();
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
                const TacsScalar T[],
                const TacsScalar strain[],
                TacsScalar * fail );
  void addFailureDVSens( const double pt[], 
                         const TacsScalar T[],
                         const TacsScalar strain[],
                         TacsScalar alpha,
                         TacsScalar dvSens[], int dvLen );
  void failureStrainSens( const double pt[],
                          const TacsScalar T[],
                          const TacsScalar strain[],
                          TacsScalar sens[],
                          int vars_j=0 );

  // Return the density as the design variable
  // -----------------------------------------
  TacsScalar getDVOutputValue( int dvIndex, const double pt[] );

  // Compute the conduction term
  // -------------------------------------------------------------
  TacsScalar getEffThermalAlpha( int vars_j=0 );
  void calculateThermal( const double pt[],
                        const TacsScalar e[], TacsScalar s[] );
  void addThermalDVSens( const double pt[], 
                         const TacsScalar e[], 
                         TacsScalar alpha, 
                         const TacsScalar psi[], 
                         TacsScalar fdvSens[],
                         int dvLen );
  void calculateConduction( const double pt[],
                            const TacsScalar e[], TacsScalar s[] );
  void addConductionDVSens( const double pt[], 
                            const TacsScalar e[], 
                            TacsScalar alpha, 
                            const TacsScalar psi[], 
                            TacsScalar fdvSens[],
                            int dvLen );
  TacsScalar getReferenceTemperature();
  int getVarsPerNode(){
    return nvars;
  }
 private:
  void computeLocalWeights( const double pt[], 
                            TMRIndexWeight *local_weights );

  // The density and stiffness properties
  TMRQuadStiffnessProperties *props;
  int nmats, nvars;
  TacsScalar x[MAX_NUM_MATERIALS+1];
  TacsScalar rho[MAX_NUM_MATERIALS+1];
  
  // The local weights
  int nweights;
  TMRIndexWeight *weights;
  int *index;
  // Local design variable vector
  TacsScalar *dv;
  TMRQuadForest *filter;
};

#endif // TMR_THERMO_QUADRANT_STIFFNESS_H
