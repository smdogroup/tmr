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

#ifndef TMR_LAGRANGE_FILTER_H
#define TMR_LAGRANGE_FILTER_H

#include "TMRTopoFilter.h"
#include "TMROctForest.h"
#include "TMRQuadForest.h"
#include "TACSAssembler.h"
#include "TMR_STLTools.h"

/*
  Form a filter for the Lagrange interpolation filter
*/
class TMRLagrangeFilter : public TMRTopoFilter {
 public:
  TMRLagrangeFilter( int _nlevels,
                     TACSAssembler *_tacs[],
                     TMROctForest *_filter[],
                     TACSVarMap *_filter_maps[],
                     TACSBVecIndices *filter_indices[],
                     int _vars_per_node=1 );
  TMRLagrangeFilter( int _nlevels,
                     TACSAssembler *_tacs[],
                     TMRQuadForest *_filter[],
                     TACSVarMap *_filter_maps[],
                     TACSBVecIndices *filter_indices[],
                     int _vars_per_node=1 );
  ~TMRLagrangeFilter();

  // Get the MPI communicator
  MPI_Comm getMPIComm();

  // Get the design variable map
  TACSVarMap *getDesignVarMap();

  // Get the TACSAssembler instance (on the finest mesh level)
  TACSAssembler *getAssembler();

  // Get problem definitions maximum local size of the design variable values
  int getVarsPerNode();
  int getNumLocalVars();
  int getMaxNumLocalVars();

  // Create a design vector on the finest mesh level
  TACSBVec *createVec();

  // Set the design variable values (including all local values)
  void setDesignVars( TACSBVec *x );

  // Set values/add values to the vector
  void addValues( TacsScalar *in, TACSBVec *out );
  void setValues( TacsScalar *in, TACSBVec *out );

  // Write the STL file
  void writeSTLFile( int k, double cutoff, const char *filename ){
    if (oct_filter){
      TMR_GenerateBinFile(filename, oct_filter[0], x[0], k, cutoff);
    }
  }
 private:
  // Initialize the problem
  void initialize( int _nlevels,
                   TACSAssembler *_tacs[],
                   TMROctForest *_oct_filter[],
                   TMRQuadForest *_quad_filter[],
                   TACSVarMap *_filter_maps[],
                   TACSBVecIndices *filter_indices[],
                   int _vars_per_node=1 );

  // Get/set values from the TACSBVec object
  int getLocalValuesFromBVec( int level, TACSBVec *vec, TacsScalar *xloc );
  void setBVecFromLocalValues( int level, const TacsScalar *xloc, TACSBVec *vec,
                               TACSBVecOperation op );

  // The number of multigrid levels
  int nlevels;
  TACSAssembler **tacs;

  // The number of variables per node
  int vars_per_node;

  // The maximum number of local design variables
  int max_local_vars;

  // Set the information about the filter at each level
  TMROctForest **oct_filter;
  TMRQuadForest **quad_filter;
  TACSVarMap **filter_maps;
  TACSBVecDistribute **filter_dist;
  TACSBVecInterp **filter_interp;

  // Create the design variable values at each level
  TACSBVec **x;
};

#endif // TMR_LAGRANGE_FILTER_H
