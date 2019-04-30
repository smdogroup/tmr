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
                                      TACSAssembler *_tacs[],
                                      TMROctForest *_filter[],
                                      TACSVarMap *_filter_maps[],
                                      TACSBVecIndices *filter_indices[],
                                      int _vars_per_node ){
  initialize(_nlevels, _tacs, _filter, NULL,
             _filter_maps, filter_indices, _vars_per_node);
}

TMRLagrangeFilter::TMRLagrangeFilter( int _nlevels,
                                      TACSAssembler *_tacs[],
                                      TMRQuadForest *_filter[],
                                      TACSVarMap *_filter_maps[],
                                      TACSBVecIndices *filter_indices[],
                                      int _vars_per_node ){
  initialize(_nlevels, _tacs, NULL, _filter,
             _filter_maps, filter_indices, _vars_per_node);
}

void TMRLagrangeFilter::initialize( int _nlevels,
                                    TACSAssembler *_tacs[],
                                    TMROctForest *_oct_filter[],
                                    TMRQuadForest *_quad_filter[],
                                    TACSVarMap *_filter_maps[],
                                    TACSBVecIndices *filter_indices[],
                                    int _vars_per_node ){
  // Set the number of multigrid levels and the number of variables
  // per node
  nlevels = _nlevels;
  vars_per_node = _vars_per_node;
  max_local_vars = 0;
  
  // Allocate arrays to store the assembler objects/forests
  tacs = new TACSAssembler*[ nlevels ];

  oct_filter = NULL;
  quad_filter = NULL;
  if (_oct_filter){
    oct_filter = new TMROctForest*[ nlevels ];
  }
  else {
    quad_filter = new TMRQuadForest*[ nlevels ];
  }

  for ( int k = 0; k < nlevels; k++ ){
    // Set the TACSAssembler objects for each level
    tacs[k] = _tacs[k];
    tacs[k]->incref();

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

  // Allocate arrays to store the filter data
  filter_maps = new TACSVarMap*[ nlevels ];
  filter_dist = new TACSBVecDistribute*[ nlevels ];

  // The design variable vector for each level
  x = new TACSBVec*[ nlevels ];

  // Copy over the assembler objects and filters
  for ( int k = 0; k < nlevels; k++ ){
    filter_maps[k] = _filter_maps[k];
    filter_maps[k]->incref();

    // Create the distribution object for the design variables
    filter_dist[k] = new TACSBVecDistribute(filter_maps[k],
                                            filter_indices[k]);
    filter_dist[k]->incref();

    x[k] = new TACSBVec(filter_maps[k], vars_per_node, filter_dist[k]);
    x[k]->incref();
  }

  // Now create the interpolation between filter levels
  filter_interp = new TACSBVecInterp*[ nlevels-1 ];

  for ( int k = 1; k < nlevels; k++ ){
    // Create the interpolation object
    filter_interp[k-1] = new TACSBVecInterp(filter_maps[k],
                                            filter_maps[k-1], vars_per_node);
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

  // Get the rank
  int mpi_rank;
  MPI_Comm_rank(_tacs[0]->getMPIComm(), &mpi_rank);

  // Set the max. number of local variables
  max_local_vars = 0;

  // Compute the max filter size
  for ( int k = 0; k < nlevels; k++ ){
    // Set the maximum local size
    const int *range;
    filter_maps[k]->getOwnerRange(&range);

    // Get the global indices and add its contribution to the
    // max. size of the local design variable vector
    TACSBVecIndices *indices = filter_dist[k]->getIndices();
    int size = vars_per_node*((range[mpi_rank+1] - range[mpi_rank]) +
                              indices->getNumIndices());

    // Compute the total number of referenced nodes on this filter object
    int size2 = 0;
    if (oct_filter){
      size2 = vars_per_node*oct_filter[k]->getNodeNumbers(NULL);
    }
    else {
      size2 = vars_per_node*quad_filter[k]->getNodeNumbers(NULL);
    }

    // Update the maximum local size
    if (size > max_local_vars){
      max_local_vars = size;
    }
    if (size2 > max_local_vars){
      max_local_vars = size2;
    }
  }
}

/*
  Free all the data for the filter object
*/
TMRLagrangeFilter::~TMRLagrangeFilter(){
  // Decrease the reference counts
  for ( int k = 0; k < nlevels; k++ ){
    tacs[k]->decref();
    if (quad_filter && quad_filter[k]){
      quad_filter[k]->decref();
    }
    if (oct_filter && oct_filter[k]){
      oct_filter[k]->decref();
    }
    filter_maps[k]->decref();
    filter_dist[k]->decref();
    x[k]->decref();
  }
  delete [] tacs;
  if (quad_filter){
    delete [] quad_filter;
  }
  if (oct_filter){
    delete [] oct_filter;
  }

  delete [] filter_maps;
  delete [] filter_dist;
  delete [] x;

  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->decref();
  }
  delete [] filter_interp;
}

/*
  Get the MPI communicator
*/
MPI_Comm TMRLagrangeFilter::getMPIComm(){
  return tacs[0]->getMPIComm();
}

/*
  Get the design variable mapping for the finest mesh
*/
TACSVarMap* TMRLagrangeFilter::getDesignVarMap(){
  return filter_maps[0];
}

/*
  Get the root assembler object
*/
TACSAssembler* TMRLagrangeFilter::getAssembler(){
  return tacs[0];
}

/*
  Get problem definitions maximum local size of the design variable values
*/
int TMRLagrangeFilter::getVarsPerNode(){
  return vars_per_node;
}

/*
  Get the number of local owned variables
*/
int TMRLagrangeFilter::getNumLocalVars(){
  return x[0]->getArray(NULL);
}

/*
  Get the number of local variables referenced by this proc
*/
int TMRLagrangeFilter::getMaxNumLocalVars(){
  return max_local_vars;
}

/*
  Create a design vector for this problem
*/
TACSBVec* TMRLagrangeFilter::createVec(){
  return new TACSBVec(filter_maps[0], vars_per_node, filter_dist[0]);
}

/*
  Set the design variables for each level
*/
void TMRLagrangeFilter::setDesignVars( TACSBVec *xvec ){
  // Copy the values to the local design variable vector
  x[0]->copyValues(xvec);

  // Distribute the design variable values
  x[0]->beginDistributeValues();
  x[0]->endDistributeValues();

  // Temporarily allocate an array to store the variables
  TacsScalar *xlocal = new TacsScalar[ getMaxNumLocalVars() ];

  // Copy the values to the local array
  int size = getLocalValuesFromBVec(0, x[0], xlocal);
  tacs[0]->setDesignVars(xlocal, size);

  // Set the design variable values on all processors
  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->multWeightTranspose(x[k], x[k+1]);
    // Distribute the design variable values
    x[k+1]->beginDistributeValues();
    x[k+1]->endDistributeValues();

    // Set the design variable values
    size = getLocalValuesFromBVec(k+1, x[k+1], xlocal);
    tacs[k+1]->setDesignVars(xlocal, size);
  }

  delete [] xlocal;
}

/*
  Add values to the output TACSBVec
*/
void TMRLagrangeFilter::addValues( TacsScalar *xlocal, TACSBVec *vec ){
  setBVecFromLocalValues(0, xlocal, vec, TACS_ADD_VALUES);
  vec->beginSetValues(TACS_ADD_VALUES);
  vec->endSetValues(TACS_ADD_VALUES);
}

/*
  Set values to the output vector
*/
void TMRLagrangeFilter::setValues( TacsScalar *xlocal, TACSBVec *vec ){
  setBVecFromLocalValues(0, xlocal, vec, TACS_INSERT_VALUES);
  vec->beginDistributeValues();
  vec->endDistributeValues();
  vec->beginSetValues(TACS_INSERT_VALUES);
  vec->endSetValues(TACS_INSERT_VALUES);
}

/*
  Get the local values of the design variables from the TACSBVec
  object and set them in xlocal.

  This function requires that the values already be distributed
  between processors so that the external values (owned by other
  processors) are correct.

  input:
  vec:     the design variable vector

  output:
  xlocal:  the local design variable array

  returns: the number of design variables on this processor
*/
int TMRLagrangeFilter::getLocalValuesFromBVec( int level,
                                               TACSBVec *vec,
                                               TacsScalar *xloc ){
  TacsScalar *x_vals, *x_ext_vals;
  int size = vec->getArray(&x_vals);
  int ext_size = vec->getExtArray(&x_ext_vals);

  memcpy(xloc, x_vals, size*sizeof(TacsScalar));
  if (x_ext_vals){
    memcpy(&xloc[size], x_ext_vals, ext_size*sizeof(TacsScalar));
  }
  return size + ext_size;
}

/*
  Set the values into the TACSBVec object using the local
  size
*/
void TMRLagrangeFilter::setBVecFromLocalValues( int level,
                                                const TacsScalar *xloc,
                                                TACSBVec *vec,
                                                TACSBVecOperation op ){
  TacsScalar *x_vals, *x_ext_vals;
  int size = vec->getArray(&x_vals);
  int ext_size = vec->getExtArray(&x_ext_vals);
  memcpy(x_vals, xloc, size*sizeof(TacsScalar));
  if (x_ext_vals){
    memcpy(x_ext_vals, &xloc[size], ext_size*sizeof(TacsScalar));
  }
}
