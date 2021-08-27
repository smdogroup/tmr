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

#include "TMRLagrangeFilter.h"

TMRLagrangeFilter::TMRLagrangeFilter( int _nlevels,
                                      TACSAssembler *_assembler[],
                                      TMROctForest *_filter[] ){
  initialize(_nlevels, _assembler, _filter, NULL);
}

TMRLagrangeFilter::TMRLagrangeFilter( int _nlevels,
                                      TACSAssembler *_assembler[],
                                      TMRQuadForest *_filter[] ){
  initialize(_nlevels, _assembler, NULL, _filter);
}

void TMRLagrangeFilter::initialize( int _nlevels,
                                    TACSAssembler *_assembler[],
                                    TMROctForest *_oct_filter[],
                                    TMRQuadForest *_quad_filter[] ){
  // Set the number of multigrid levels and the number of variables
  // per node
  nlevels = _nlevels;
  
  // Allocate arrays to store the assembler objects/forests
  assembler = new TACSAssembler*[ nlevels ];

  oct_filter = NULL;
  quad_filter = NULL;
  if (_oct_filter){
    oct_filter = new TMROctForest*[ nlevels ];
  }
  else {
    quad_filter = new TMRQuadForest*[ nlevels ];
  }

  // The design variable vector for each level
  x = new TACSBVec*[ nlevels ];

  for ( int k = 0; k < nlevels; k++ ){
    // Set the TACSAssembler objects for each level
    assembler[k] = _assembler[k];
    assembler[k]->incref();

    x[k] = assembler[k]->createDesignVec();
    x[k]->incref();

    // Set the filter object
    if (_oct_filter){
      oct_filter[k] = _oct_filter[k];
      oct_filter[k]->incref();
    }
    else {
      quad_filter[k] = _quad_filter[k];
      quad_filter[k]->incref();
    }
  }

  // Now create the interpolation between filter levels
  filter_interp = new TACSBVecInterp*[ nlevels-1 ];

  for ( int k = 1; k < nlevels; k++ ){
    // Create the interpolation object
    filter_interp[k-1] = new TACSBVecInterp(assembler[k]->getDesignNodeMap(),
                                            assembler[k-1]->getDesignNodeMap(),
                                            assembler[0]->getDesignVarsPerNode());
    filter_interp[k-1]->incref();

    // Create the interpolation on the TMR side
    if (oct_filter){
      oct_filter[k-1]->createInterpolation(oct_filter[k], filter_interp[k-1]);
    }
    else {
      quad_filter[k-1]->createInterpolation(quad_filter[k], filter_interp[k-1]);
    }
    filter_interp[k-1]->initialize();
  }
}

/*
  Free all the data for the filter object
*/
TMRLagrangeFilter::~TMRLagrangeFilter(){
  // Decrease the reference counts
  for ( int k = 0; k < nlevels; k++ ){
    assembler[k]->decref();
    if (quad_filter && quad_filter[k]){
      quad_filter[k]->decref();
    }
    if (oct_filter && oct_filter[k]){
      oct_filter[k]->decref();
    }
    x[k]->decref();
  }
  delete [] assembler;
  if (quad_filter){
    delete [] quad_filter;
  }
  if (oct_filter){
    delete [] oct_filter;
  }
  delete [] x;

  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->decref();
  }
  delete [] filter_interp;
}

/*
  Get the root assembler object
*/
TACSAssembler* TMRLagrangeFilter::getAssembler(){
  return assembler[0];
}

/*
  Get the root TMRQuadForest object (if any)
*/
TMRQuadForest* TMRLagrangeFilter::getFilterQuadForest(){
  if (quad_filter){
    return quad_filter[0];
  }
  return NULL;
}

/*
  Get the root TMROctForest object (if any)
*/
TMROctForest* TMRLagrangeFilter::getFilterOctForest(){
  if (oct_filter){
    return oct_filter[0];
  }
  return NULL;
}

/*
  Set the design variables for each level
*/
void TMRLagrangeFilter::setDesignVars( TACSBVec *xvec ){
  // Copy the values to the local design variable vector
  x[0]->copyValues(xvec);
  assembler[0]->setDesignVars(x[0]);

  // Set the design variable values on all processors
  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->multWeightTranspose(x[k], x[k+1]);
    assembler[k+1]->setDesignVars(x[k+1]);
  }
}

/*
  Add values to the output TACSBVec
*/
void TMRLagrangeFilter::addValues( TACSBVec *vec ){
  vec->beginSetValues(TACS_ADD_VALUES);
  vec->endSetValues(TACS_ADD_VALUES);
}

/*
  Set values to the output vector
*/
void TMRLagrangeFilter::setValues( TACSBVec *vec ){
  vec->beginSetValues(TACS_INSERT_VALUES);
  vec->endSetValues(TACS_INSERT_VALUES);
}
