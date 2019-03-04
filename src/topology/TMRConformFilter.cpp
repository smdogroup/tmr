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

#include "TMRConformFilter.h"

TMRConformFilter::TMRConformFilter( int _nlevels,
                                    TACSAssembler *_tacs[],
                                    TMROctForest *_filter[],
                                    int _vars_per_node ){
  initialize(_nlevels, _tacs, _filter, NULL, _vars_per_node);
}

TMRConformFilter::TMRConformFilter( int _nlevels,
                                    TACSAssembler *_tacs[],
                                    TMRQuadForest *_filter[],
                                    int _vars_per_node ){
  initialize(_nlevels, _tacs, NULL, _filter, _vars_per_node);
}

void TMRConformFilter::initialize( int _nlevels,
                                   TACSAssembler *_tacs[],
                                   TMROctForest *_oct_filter[],
                                   TMRQuadForest *_quad_filter[],
                                   int _vars_per_node ){
  nlevels = _nlevels;

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
  filter_dep_nodes = new TACSBVecDepNodes*[ nlevels ];

  // The design variable vector for each level
  x = new TACSBVec*[ nlevels ];

  // Copy over the assembler objects and filters
  for ( int k = 0; k < nlevels; k++ ){
    TACSBVecIndices *filter_indices;
    createMapIndices(quad_filter[k], oct_filter[k],
                     &filter_maps[k], &filter_indices);
    filter_maps[k]->incref();

    // Create the distribution object for the design variables
    filter_dist[k] = new TACSBVecDistribute(filter_maps[k],
                                            filter_indices);
    filter_dist[k]->incref();

    // Extract the dependent node info from the oct filter
    int *dep_ptr, *dep_conn;
    const int *_dep_ptr, *_dep_conn;
    double *dep_weights;
    const double *_dep_weights;
    int num_dep_nodes = 0;
    if (oct_filter){
      num_dep_nodes =
        oct_filter[k]->getDepNodeConn(&_dep_ptr, &_dep_conn,
                                      &_dep_weights);
    }
    else {
      num_dep_nodes =
        quad_filter[k]->getDepNodeConn(&_dep_ptr, &_dep_conn,
                                       &_dep_weights);
    }

    // Copy over the dependent node data
    dep_ptr = new int[ num_dep_nodes+1 ];
    memcpy(dep_ptr, _dep_ptr, (num_dep_nodes+1)*sizeof(int));

    int dep_size = dep_ptr[num_dep_nodes];
    dep_conn = new int[ dep_size ];
    memcpy(dep_conn, _dep_conn, dep_size*sizeof(int));

    dep_weights = new double[ dep_size ];
    memcpy(dep_weights, _dep_weights, dep_size*sizeof(double));

    // Create the dependent node object
    filter_dep_nodes[k] = new TACSBVecDepNodes(num_dep_nodes,
                                               &dep_ptr, &dep_conn,
                                               &dep_weights);
    filter_dep_nodes[k]->incref();

    // Create the vectors for each mesh level
    x[k] = new TACSBVec(filter_maps[k], vars_per_node, filter_dist[k],
                        filter_dep_nodes[k]);
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

  // Compute the max filter size
  max_local_vars = 0;
  for ( int k = 0; k < nlevels; k++ ){
    // Compute the total number of referenced nodes on this filter object
    int size = 0;
    if (oct_filter){
      size = vars_per_node*oct_filter[k]->getNodeNumbers(NULL);
    }
    else {
      size = vars_per_node*quad_filter[k]->getNodeNumbers(NULL);
    }

    // Update the maximum local size
    if (size > max_local_vars){
      max_local_vars = size;
    }
  }
}

TMRConformFilter::~TMRConformFilter(){
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
    filter_dep_nodes[k]->decref();
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
  delete [] filter_dep_nodes;
  delete [] x;

  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->decref();
  }
  delete [] filter_interp;
}

/*
  Based on the input filter information, create the variable range/
  index data that will be used to distribute the design variable
  values
*/
void TMRConformFilter::createMapIndices( TMRQuadForest *quad_filter,
                                         TMROctForest *oct_filter,
                                         TACSVarMap **filter_map,
                                         TACSBVecIndices **filter_indices ){
  int mpi_rank, mpi_size;
  MPI_Comm comm = getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Get the node range for the filter design variables
  const int *filter_range;
  if (quad_filter){
    quad_filter->getOwnedNodeRange(&filter_range);
  }
  else {
    oct_filter->getOwnedNodeRange(&filter_range);
  }

  // Set up the variable map for the design variable numbers
  int num_filter_local = filter_range[mpi_rank+1] - filter_range[mpi_rank];
  *filter_map = new TACSVarMap(comm, num_filter_local);

  // Get the local node numbers
  const int *node_nums;
  int num_nodes = 0, num_dep_nodes = 0;

  if (quad_filter){
    quad_filter->getNodeNumbers(&node_nums);
    num_dep_nodes = quad_filter->getDepNodeConn();
  }
  else {
    oct_filter->getNodeNumbers(&node_nums);
    num_dep_nodes = oct_filter->getDepNodeConn();
  }

  // Set the number of filter nodes
  int num_filter_nodes = 0;

  // Compute the number of external node numbers
  for ( int i = num_dep_nodes; i < num_nodes; i++ ){
    if (node_nums[i] >= 0 &&
        (node_nums[i] < filter_range[mpi_rank] ||
         node_nums[i] >= filter_range[mpi_rank+1])){
      num_filter_nodes++;
    }
  }

  // Allocate an array for the external node numbers and set their values into
  // the filter_node_nums array
  int *filter_node_nums = new int[ num_filter_nodes ];
  int count = 0;
  for ( int i = num_dep_nodes; i < num_nodes; i++ ){
    if (node_nums[i] >= 0 &&
        (node_nums[i] < filter_range[mpi_rank] ||
         node_nums[i] >= filter_range[mpi_rank+1])){
      filter_node_nums[count] = node_nums[i];
      count++;
    }
  }

  // Set up the external filter indices for this filter.  The indices
  // objects steals the array for the external nodes.
  *filter_indices = new TACSBVecIndices(&filter_node_nums,
                                        num_filter_nodes);
  (*filter_indices)->setUpInverse();
}

/*
  Get the root assembler object
*/
TACSAssembler* TMRConformFilter::getAssembler(){
  return tacs[0];
}

/*
  Get problem definitions maximum local size of the design variable values
*/
int TMRConformFilter::getVarsPerNode(){
  return vars_per_node;
}

/*
  Get the number of local owned variables
*/
int TMRConformFilter::getNumLocalVars(){
  return x[0]->getArray(NULL);
}

/*
  Get the number of local variables referenced by this proc
*/
int TMRConformFilter::getMaxNumLocalVars(){
  return max_local_vars;
}

/*
  Create the design TACSBVec vector
*/
TACSBVec* TMRConformFilter::createVec(){
  return new TACSBVec(filter_maps[0], vars_per_node,
                      filter_dist[0], filter_dep_nodes[0]);
}

/*
  Set the design variables for each level
*/
void TMRConformFilter::setDesignVars( TACSBVec *xvec ){
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
void TMRConformFilter::addValues( TacsScalar *xlocal, TACSBVec *vec ){
  setBVecFromLocalValues(0, xlocal, vec, TACS_ADD_VALUES);
  vec->beginSetValues(TACS_ADD_VALUES);
  vec->endSetValues(TACS_ADD_VALUES);
}

/*
  Set values to the output vector
*/
void TMRConformFilter::setValues( TacsScalar *xlocal, TACSBVec *vec ){
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
int TMRConformFilter::getLocalValuesFromBVec( int level,
                                              TACSBVec *vec,
                                              TacsScalar *xloc ){
  // If it is for a quadtree problem
  if (quad_filter){
    const int *node_numbers;
    int num_nodes = quad_filter[level]->getNodeNumbers(&node_numbers);
    vec->getValues(num_nodes, node_numbers, &xloc[0]);
    return num_nodes;
  }
  else { // if (oct_filter){
    const int *node_numbers;
    int num_nodes = oct_filter[level]->getNodeNumbers(&node_numbers);
    vec->getValues(num_nodes, node_numbers, &xloc[0]);
    return num_nodes;
  }
}

/*
  Set the values into the TACSBVec object using the local
  size
*/
void TMRConformFilter::setBVecFromLocalValues( int level,
                                               const TacsScalar *xloc,
                                               TACSBVec *vec,
                                               TACSBVecOperation op ){
  if (quad_filter){
    const int *node_numbers;
    int num_nodes = quad_filter[level]->getNodeNumbers(&node_numbers);
    vec->setValues(num_nodes, node_numbers, xloc, op);
  }
  else if (oct_filter){
    const int *node_numbers;
    int num_nodes = oct_filter[level]->getNodeNumbers(&node_numbers);
    vec->setValues(num_nodes, node_numbers, xloc, op);
  }
}
