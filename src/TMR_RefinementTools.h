#ifndef TMR_REFINEMENT_TOOLS_H
#define TMR_REFINEMENT_TOOLS_H

#include "TMRQuadForest.h"
#include "TACSAssembler.h"

/*
  Perform a mesh refinement based on the strain engery refinement
  criteria.
*/
TacsScalar TMR_StrainEnergyRefine( TACSAssembler *tacs,
                                   TMRQuadForest *forest,
                                   double target_err,
                                   int min_level=0, int max_level=TMR_MAX_LEVEL );

#endif // TMR_REFINEMENT_TOOLS
