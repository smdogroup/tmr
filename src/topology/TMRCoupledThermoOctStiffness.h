#ifndef TMR_COUPLED_THERMO_OCTANT_STIFFNESS_H
#define TMR_COUPLED_THERMO_OCTANT_STIFFNESS_H

/*
  The following class defines a SolidStiffness base class for topology
  optimization using TACS.  
*/

#include "TMRBase.h"
#include "CoupledThermoSolidStiffness.h"
#include "TACSConstitutive.h"
#include "TMROctStiffness.h"
#include "TMROctForest.h"
/*
  The TMRCoupledThermoOctStiffness class

  This defines the TMRCoupledThermoOctStiffness class which takes the weights from
  up to 8 adjacent vertices. This uses the RAMP method for
  penalization.
*/
class TMRCoupledThermoOctStiffness : public CoupledThermoSolidStiffness {
 public: 
  static const int MAX_NUM_MATERIALS = 5;
  TMRCoupledThermoOctStiffness( TMRIndexWeight *_weights, 
                                int _nweights,
                                TMRStiffnessProperties *_props,
                                TMROctForest *_filter=NULL,
                                int *_index=NULL );
  ~TMRCoupledThermoOctStiffness();
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
  int getVarsPerNode(){
    return nvars;
  }
  void heatflux( const double pt[],
                 const TacsScalar normal[],
                 const TacsScalar strain[],
                 TacsScalar * qn );
  void addHeatFluxDVSens( const double pt[],
                          const TacsScalar normal[],
                          const TacsScalar strain[],
                          TacsScalar alpha,
                          TacsScalar dvSens[], int dvLen );
  void heatfluxStrainSens( const double pt[],
                           const TacsScalar normal[],
                           const TacsScalar strain[],
                           TacsScalar sens[], 
                           int vars_j=0 );

  
 private:
  // The density and stiffness properties
  TMRStiffnessProperties *props;
  int nmats, nvars;
  TacsScalar x[MAX_NUM_MATERIALS+1];
  TacsScalar rho[MAX_NUM_MATERIALS+1];

  // The local weights
  int nweights;
  TMRIndexWeight *weights;
  int *index;
  // Local design variable vector
  TacsScalar *dv;
  TMROctForest *filter;
};

#endif // TMR_THERMO_OCTANT_STIFFNESS_H
