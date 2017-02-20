#include "TMRBase.h"
#include <string.h>

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

// Static flag to test if TMR is initialized or not
static int TMR_is_initialized = 0;

// The TMR data type for MPI useage
MPI_Datatype TMROctant_MPI_type;
MPI_Datatype TMRQuadrant_MPI_type;
MPI_Datatype TMRPoint_MPI_type;

/*
  Initialize TMR data type
*/
void TMRInitialize(){
  if (!TMR_is_initialized){
    MPI_Aint offset = 0;
    MPI_Datatype type = MPI_INT32_T;

    // Create the TMROctant data type
    int counts = 6;
    MPI_Type_create_struct(1, &counts, &offset, &type, 
                           &TMROctant_MPI_type);
    MPI_Type_commit(&TMROctant_MPI_type);

    // Create the TMRQudrant data type
    counts = 5;
    MPI_Type_create_struct(1, &counts, &offset, &type, 
                           &TMRQuadrant_MPI_type);
    MPI_Type_commit(&TMRQuadrant_MPI_type);

    // Create the TMRPoint data type
    counts = 3;
    type = MPI_DOUBLE;
    MPI_Type_create_struct(1, &counts, &offset, &type, 
                           &TMRPoint_MPI_type);
    MPI_Type_commit(&TMRPoint_MPI_type);
    
    // Set the TMR initialization flag
    TMR_is_initialized = 1;
  }
}

/*
  Check whether the TMR data types have been initialized or not
*/
int TMRIsInitialize(){
  return TMR_is_initialized;
}

/*
  Finalize the TMR data type
*/
void TMRFinalize(){
  MPI_Type_free(&TMROctant_MPI_type);
  MPI_Type_free(&TMRQuadrant_MPI_type);
  MPI_Type_free(&TMRPoint_MPI_type);
}

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
