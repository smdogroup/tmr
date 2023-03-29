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

#ifndef TMR_TACS_TOPO_CREATOR_H
#define TMR_TACS_TOPO_CREATOR_H

#include "TACSAssembler.h"
#include "TMR_TACSCreator.h"

/*
  This is an abstract base class used to create octforests specialized
  for topology optimization. This class can be overriden with the
  createElement function to create different types of topolgy
  optimization problems for 3D structures.
*/
class TMROctTACSTopoCreator : public TMROctTACSCreator {
 public:
  TMROctTACSTopoCreator(TMRBoundaryConditions *_bcs, int _design_vars_per_node,
                        TMROctForest *_filter);
  ~TMROctTACSTopoCreator();

  // Create the elements
  void createElements(int order, TMROctForest *forest, int num_elements,
                      TACSElement **elements);

  // Create the element
  virtual TACSElement *createElement(int order, TMROctant *oct, int nweights,
                                     TMRIndexWeight *weights) = 0;

 private:
  // Compute the weights for a given point
  void computeWeights(const int mesh_order, const double *knots,
                      TMROctant *node, TMROctant *oct, TMRIndexWeight *weights,
                      double *tmp);

  // The filter indices. This defines the relationship between the
  // local design variable numbers and the global design variable
  // numbers.
  TACSBVecIndices *filter_indices;
};

/*
  This class is an abstract base class for setting up
  topology optimzation problems using TMRQuadForest objects
*/
class TMRQuadTACSTopoCreator : public TMRQuadTACSCreator {
 public:
  TMRQuadTACSTopoCreator(TMRBoundaryConditions *_bcs, int _design_vars_per_node,
                         TMRQuadForest *_filter);
  ~TMRQuadTACSTopoCreator();

  // Create the elements
  void createElements(int order, TMRQuadForest *forest, int num_elements,
                      TACSElement **elements);

  // Create the element
  virtual TACSElement *createElement(int order, TMRQuadrant *oct, int nweights,
                                     TMRIndexWeight *weights) = 0;

 private:
  // Compute the weights for a given point
  void computeWeights(const int mesh_order, const double *knots,
                      TMRQuadrant *node, TMRQuadrant *quad,
                      TMRIndexWeight *weights, double *tmp, int sort = 1);

  // The filter indices. This defines the relationship between the
  // local design variable numbers and the global design variable
  // numbers.
  TACSBVecIndices *filter_indices;
};

/*
  This class is an abstract base class for setting up topology optimzation
  problems using TMROctForest objects. This is simplified for when the
  underlying forest for the analysis and design mesh is identical
*/
class TMROctConformTACSTopoCreator : public TMROctTACSCreator {
 public:
  TMROctConformTACSTopoCreator(
      TMRBoundaryConditions *_bcs, int _design_vars_per_node,
      TMROctForest *_forest, int order = -1,
      TMRInterpolationType interp_type = TMR_UNIFORM_POINTS);
  ~TMROctConformTACSTopoCreator();

  // Create the elements
  void createElements(int order, TMROctForest *forest, int num_elements,
                      TACSElement **elements);

  // Create the element
  virtual TACSElement *createElement(int order, TMROctant *oct, int nconn,
                                     const int *index,
                                     TMROctForest *filter) = 0;
};

/*
  This class is an abstract base class for setting up topology optimzation
  problems using TMRQuadForest objects. This is simplified for when the
  underlying forest for the analysis and design mesh is identical
*/
class TMRQuadConformTACSTopoCreator : public TMRQuadTACSCreator {
 public:
  TMRQuadConformTACSTopoCreator(
      TMRBoundaryConditions *_bcs, int _design_vars_per_node,
      TMRQuadForest *_forest, int order = -1,
      TMRInterpolationType interp_type = TMR_UNIFORM_POINTS);
  ~TMRQuadConformTACSTopoCreator();

  // Create the elements
  void createElements(int order, TMRQuadForest *forest, int num_elements,
                      TACSElement **elements);

  // Create the element
  virtual TACSElement *createElement(int order, TMRQuadrant *oct, int nconn,
                                     const int *index,
                                     TMRQuadForest *filter) = 0;
};

#endif  // TMR_TACS_TOPO_CREATOR_H
