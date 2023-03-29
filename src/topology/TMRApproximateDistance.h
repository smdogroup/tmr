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

#ifndef TMR_APPROXIMATE_DISTANCE_H
#define TMR_APPROXIMATE_DISTANCE_H

#include "TACSAssembler.h"
#include "TMROctForest.h"
#include "TMRQuadForest.h"

void TMRApproximateDistance(TMRQuadForest *filter, int index, double cutoff,
                            double t, TACSBVec *rho, const char *filename,
                            double *min_dist);
void TMRApproximateDistance(TMROctForest *filter, int index, double cutoff,
                            double t, TACSBVec *rho, const char *filename,
                            double *min_dist);

#endif  // TMR_APPROXIMATE_DISTANCE_H
