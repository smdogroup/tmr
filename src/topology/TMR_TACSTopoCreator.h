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

#include "TMR_TACSCreator.h"
#include "TACSAssembler.h"
#include "TMROctStiffness.h"

/*
  This is an abstract base class used to create octforests specialized
  for topology optimization. This class can be overriden with the
  createElement function to create different types of topolgy
  optimization problems for 3D structures.
*/
class TMROctTACSTopoCreator : public TMROctTACSCreator {
 public:
  TMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                         TMROctForest *_filter );
  ~TMROctTACSTopoCreator();

  // Create the elements
  void createElements( int order, 
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements );

  // Create the element
  virtual TACSElement *createElement( int order, 
                                      TMROctant *oct,
                                      TMRIndexWeight *weights, 
                                      int nweights ) = 0;

  // Get the underlying objects that define the filter
  void getFilter( TMROctForest **filter );
  void getMap( TACSVarMap **_map );
  void getIndices( TACSBVecIndices **_indices );  

 private:
  // Compute the weights for a given point
  void computeWeights( const int mesh_order, const double *knots,
                       TMROctant *node, TMROctant *oct,
                       TMRIndexWeight *weights, double *tmp,
                       int sort=1);

  // The forest that defines the filter
  TMROctForest *filter;

  // The filter map for this object. This defines how the design
  // variables are distributed across all of the processors.
  TACSVarMap *filter_map;

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
  TMRQuadTACSTopoCreator( TMRBoundaryConditions *_bcs,
                          TMRQuadForest *_filter );
  ~TMRQuadTACSTopoCreator();

  // Create the elements
  void createElements( int order,
                       TMRQuadForest *forest,
                       int num_elements,
                       TACSElement **elements );

  void createElements( TMRQuadForest *forest,
                       int num_elements,
                       TACSElement **elements );

  // Create the element
  virtual TACSElement *createElement( int order, 
                                      TMRQuadrant *oct,
                                      TMRIndexWeight *weights, 
                                      int nweights ) = 0;

  // Get the underlying objects that define the filter
  void getFilter( TMRQuadForest **filter );
  void getMap( TACSVarMap **_map );
  void getIndices( TACSBVecIndices **_indices );  

 private:
  // Compute the weights for a given point
  void computeWeights( const int mesh_order, const double *knots,
                       TMRQuadrant *node, TMRQuadrant *quad,
                       TMRIndexWeight *weights, double *tmp );

  // The forest that defines the filter
  TMRQuadForest *filter;

  // The filter map for this object. This defines how the design
  // variables are distributed across all of the processors.
  TACSVarMap *filter_map;

  // The filter indices. This defines the relationship between the
  // local design variable numbers and the global design variable
  // numbers.
  TACSBVecIndices *filter_indices;
};

#endif // TMR_TACS_TOPO_CREATOR_H

