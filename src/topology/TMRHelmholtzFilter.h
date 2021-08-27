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

#ifndef TMR_HELMHOLTZ_FILTER_H
#define TMR_HELMHOLTZ_FILTER_H

#include "TMRConformFilter.h"
#include "TACSAssembler.h"
#include "TACSMg.h"

/*
  Create a Helmholtz filter object
*/
class TMRHelmholtzFilter : public TMRConformFilter {
 public:
  TMRHelmholtzFilter( double helmholtz_radius,
                      int _nlevels,
                      TACSAssembler *_assembler[],
                      TMROctForest *_filter[] );
  TMRHelmholtzFilter( double helmholtz_radius,
                      int _nlevels,
                      TACSAssembler *_assembler[],
                      TMRQuadForest *_filter[] );
  ~TMRHelmholtzFilter();

  // Set the design variable values (including all local values)
  void setDesignVars( TACSBVec *x );

  // Set values/add values to the vector
  void addValues( TACSBVec *vec );

 private:
  void initialize_helmholtz( double helmholtz_radius );

  // Apply the Helmholtz filter to a vector
  void applyFilter( TACSBVec *vec );
  void applyTranspose( TACSBVec *input, TACSBVec *output );

  // Data (that may be NULL) for the Helmholtz-based PDE filter
  TACSAssembler **helmholtz_assembler;
  TACSKsm *helmholtz_ksm;
  TACSMg *helmholtz_mg;
  TACSBVec *helmholtz_rhs, *helmholtz_psi;
  TACSBVec *helmholtz_vec;

  // Temporary vector
  TACSBVec *temp;
};

#endif // TMR_HELMHOLTZ_FILTER_H
