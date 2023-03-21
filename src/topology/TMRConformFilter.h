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

#ifndef TMR_CONFORM_FILTER_H
#define TMR_CONFORM_FILTER_H

#include "TACSAssembler.h"
#include "TMROctForest.h"
#include "TMRQuadForest.h"
#include "TMRTopoFilter.h"
#include "TMR_STLTools.h"

/*
  Build a conforming interpolation filter
*/
class TMRConformFilter : public TMRTopoFilter {
 public:
  TMRConformFilter(int _nlevels, TACSAssembler *_tacs[],
                   TMROctForest *_filter[]);
  TMRConformFilter(int _nlevels, TACSAssembler *_tacs[],
                   TMRQuadForest *_filter[]);
  ~TMRConformFilter();

  // Get the MPI communicator
  MPI_Comm getMPIComm();

  // Get the TACSAssembler instance (on the finest mesh level)
  TACSAssembler *getAssembler();

  // Get the Quad or OctForest on the finest mesh level
  TMRQuadForest *getFilterQuadForest();
  TMROctForest *getFilterOctForest();

  // Set the design variable values (including all local values)
  void setDesignVars(TACSBVec *x);

  // Set values/add values to the vector
  void addValues(TACSBVec *vec);
  void setValues(TACSBVec *vec);

  void writeSTLFile(int k, double cutoff, const char *filename) {
    if (oct_filter) {
      TMR_GenerateBinFile(filename, oct_filter[0], x[0], k, cutoff);
    }
  }

 protected:
  // The number of multigrid levels
  int nlevels;
  TACSAssembler **assembler;

  // The maximum number of local design variables
  int max_local_vars;

  // Set the information about the filter at each level
  TMROctForest **oct_filter;
  TMRQuadForest **quad_filter;
  TACSBVecInterp **filter_interp;

  // Create the design variable values at each level
  TACSBVec **x;

 private:
  // Initialize the problem
  void initialize(int _nlevels, TACSAssembler *_tacs[],
                  TMROctForest *_oct_filter[], TMRQuadForest *_quad_filter[]);
};

#endif  // TMR_CONFORM_FILTER_H
