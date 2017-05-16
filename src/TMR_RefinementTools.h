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
void TMR_ComputeReconSolution( const int order,
                               TACSAssembler *tacs,
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

#endif // TMR_REFINEMENT_TOOLS_H
