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

#include "TMRBase.h"
#include "TMRQuadrant.h"
#include "TMROctant.h"
#include <stddef.h>
#include <string.h>

// Static flag to test if TMR is initialized or not
static int TMR_is_initialized = 0;

// The TMR data type for MPI useage
MPI_Datatype TMROctant_MPI_type;
MPI_Datatype TMRQuadrant_MPI_type;
MPI_Datatype TMRPoint_MPI_type;
MPI_Datatype TMRIndexWeight_MPI_type;

/*
  Initialize TMR data type
*/
void TMRInitialize(){
  if (!TMR_is_initialized){
    int counts[2];
    MPI_Aint offset[2];
    MPI_Datatype types[2];
    types[0] = MPI_INT32_T;
    types[1] = MPI_INT16_T;

    // Create the TMRQudrant data type
    counts[0] = 4;
    counts[1] = 2;
    offset[0] = offsetof(TMRQuadrant, face);
    offset[1] = offsetof(TMRQuadrant, level);
    MPI_Type_create_struct(2, counts, offset, types, 
                           &TMRQuadrant_MPI_type);
    MPI_Type_commit(&TMRQuadrant_MPI_type);

    // Create the TMROctant data type
    counts[0] = 5;
    counts[1] = 2;
    offset[0] = offsetof(TMROctant, block);
    offset[1] = offsetof(TMROctant, level);
    MPI_Type_create_struct(2, counts, offset, types, 
                           &TMROctant_MPI_type);
    MPI_Type_commit(&TMROctant_MPI_type);

    // Create the TMRPoint data type
    counts[0] = 3;
    offset[0] = 0;
    types[0] = MPI_DOUBLE;
    MPI_Type_create_struct(1, counts, offset, types, 
                           &TMRPoint_MPI_type);
    MPI_Type_commit(&TMRPoint_MPI_type);

    // Create the index/weight pair data
    int len[2] = {1, 1};
    MPI_Aint disp[2];
    disp[0] = offsetof(TMRIndexWeight, index);
    disp[1] = offsetof(TMRIndexWeight, weight);
    types[0] = MPI_INT;
    types[1] = MPI_DOUBLE;
    MPI_Type_create_struct(2, len, disp, types, 
                           &TMRIndexWeight_MPI_type);
    MPI_Type_commit(&TMRIndexWeight_MPI_type);
    
    // Set the TMR initialization flag
    TMR_is_initialized = 1;
  }
}

/*
  Check whether the TMR data types have been initialized or not
*/
int TMRIsInitialized(){
  return TMR_is_initialized;
}

/*
  Finalize the TMR data type
*/
void TMRFinalize(){
  MPI_Type_free(&TMROctant_MPI_type);
  MPI_Type_free(&TMRQuadrant_MPI_type);
  MPI_Type_free(&TMRPoint_MPI_type);
  MPI_Type_free(&TMRIndexWeight_MPI_type);
}

TMREntity::TMREntity(): entity_id(entity_id_count){
  entity_id_count++;
  attr = NULL;
  ref_count = 0;
}

TMREntity::~TMREntity(){
  if (attr){ delete [] attr; }
}

int TMREntity::entity_id_count = 0;

/*
  Set the attribute/name associate with this object
*/
void TMREntity::setAttribute( const char *_attr ){
  if (attr){
    delete [] attr;
  }
  if (_attr){
    attr = new char[ strlen(_attr)+1 ];
    strcpy(attr, _attr);
  }
}

/*
  Retrieve the attribute associated with this object
*/
const char* TMREntity::getAttribute() const {
  return attr;
}

/*
  Increment the reference count
*/
void TMREntity::incref(){ 
  ref_count++; 
}

/*
  Decrease the reference count
*/
void TMREntity::decref(){ 
  ref_count--;
  if (ref_count == 0){
    delete this;
  }
}

// Set the default convergence criteria
double TMREntity::eps_dist = 1e-6;
double TMREntity::eps_cosine = 1e-6;

/*
  Set the geometric tolerances
*/
void TMREntity::setTolerances( double _eps_dist, double _eps_cosine ){
  eps_dist = _eps_dist;
  eps_cosine = _eps_cosine;
}

/*
  Get the geometric tolerances
*/
void TMREntity::getTolerances( double *_eps_dist, double *_eps_cosine ){
  *_eps_dist = eps_dist;
  *_eps_cosine = eps_cosine;
}
