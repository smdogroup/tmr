#ifndef TMR_REFINEMENT_TOOLS_H
#define TMR_REFINEMENT_TOOLS_H

#include "TMRQuadForest.h"
#include "TMROctForest.h"
#include "TACSAssembler.h"
#include "TACSMg.h"

/*
  Create a TACS multigrid object
*/
void TMR_CreateTACSMg( int nlevels, TACSAssembler *tacs[],
                       TMROctForest *forest[],
                       TACSMg **_mg, int use_pairs=0 );
void TMR_CreateTACSMg( int nlevels, TACSAssembler *tacs[],
                       TMRQuadForest *forest[],
                       TACSMg **_mg );

/*
  Compute the reconstructed solution on a uniformly refined mesh
*/
void TMR_ComputeReconSolution( TACSAssembler *tacs,
                               TMRQuadForest *forest,
                               TACSAssembler *tacs_refined,
                               TACSBVec *_uvec=NULL,
                               TACSBVec *_uvec_refined=NULL );
void TMR_ComputeReconSolution( TACSAssembler *tacs,
                               TMROctForest *forest,
                               TACSAssembler *tacs_refined,
                               TACSBVec *_uvec=NULL,
                               TACSBVec *_uvec_refined=NULL );

/*
  Perform a mesh refinement based on the strain engery refinement
  criteria.
*/
TacsScalar TMR_StrainEnergyRefine( TACSAssembler *tacs,
                                   TMRQuadForest *forest,
                                   double target_err,
                                   int min_level=0, 
                                   int max_level=TMR_MAX_LEVEL );
TacsScalar TMR_StrainEnergyRefine( TACSAssembler *tacs,
                                   TMROctForest *forest,
                                   double target_err,
                                   int min_level=0, 
                                   int max_level=TMR_MAX_LEVEL );

/*
  Perform adjoint-based mesh refinement on the forest of quadtrees
*/
TacsScalar TMR_AdjointRefine( TACSAssembler *tacs,
                              TACSAssembler *refine,
                              TACSBVec *adjvec,
                              TMRQuadForest *forest,
                              double target_err,
                              int min_level=0, 
                              int max_level=TMR_MAX_LEVEL,
                              TacsScalar *adj_corr=NULL );

/*
  Evaluate a stress constraint based on a higher-order interpolation
  of the stresses in the problem.

  This makes strong assumptions about the element type and constitutive
  matrix. Be careful when using this method
*/
class TMRStressConstraint : public TMREntity {
 public:
  TMRStressConstraint( int order, TACSAssembler *tacs, 
                       TacsScalar _ks_weight=30.0 );
  ~TMRStressConstraint();

  // Evaluate the aggregated stress constraint across all processors
  TacsScalar evalConstraint( TACSBVec *_uvec );

  // Evaluate the terms required for the total derivative
  void evalConDeriv( TacsScalar *dfdxx, int size, TACSBVec *dfdu );

 private:
  // Evaluate the element strain at the given point
  TacsScalar evalStrain( const double pt[],
                         const TacsScalar *Xpts,
                         const TacsScalar *vars, 
                         const TacsScalar *ubar,
                         TacsScalar J[], TacsScalar e[] );

  // Evaluate the derivatives of the element strain
  TacsScalar addStrainDeriv( const double pt[], const TacsScalar J[],
                             const TacsScalar alpha, const TacsScalar dfde[],
                             TacsScalar dfdu[], TacsScalar dfdubar[] );

  // The mesh order
  int order;

  // The values used to compute the KS function
  TacsScalar ks_weight;
  TacsScalar ks_max_fail, ks_fail_sum;

  // The TACSAssembler 
  TACSAssembler *tacs;

  // The solution vector and the derivatives at the nodes
  TACSBVec *uvec, *uderiv;

  // The temporary derivative vector
  TACSBVec *dfduderiv;

  // The weights on the local 
  TACSBVec *weights;

  // Arrays to store the element-wise variables
  TacsScalar *Xpts;
  TacsScalar *vars, *dvars, *ddvars;

  // The values of the element-wise reconstructed derivatives
  // at each node
  TacsScalar *varderiv;

  // Values of the variables for each enrichment function
  TacsScalar *ubar;

  // Temporary vector
  TacsScalar *tmp;
};

#endif // TMR_REFINEMENT_TOOLS_H
