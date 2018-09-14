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

#ifndef TMR_BASE_H
#define TMR_BASE_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#define TMR_EXTERN_C_BEGIN extern "C" {
#define TMR_EXTERN_C_END }

/*
  The following constants define the maximum octant depth and maximum
  element order within the code.
*/
static const int TMR_MAX_LEVEL = 30;

/*
  Set the type of interpolation to use (only makes a difference
  for order >= 4)
*/
enum TMRInterpolationType { TMR_UNIFORM_POINTS, 
                            TMR_GAUSS_LOBATTO_POINTS,
                            TMR_BERNSTEIN_POINTS };

/*
  Base class for all point-evaluation algorithms
*/
class TMRPoint {
 public:
  TMRPoint(){}
  TMRPoint( const TMRPoint& p ){
    x = p.x;  y = p.y;  z = p.z;
  }
  TMRPoint( double _x, double _y, double _z ){
    x = _x;  y = _y;  z = _z;
  }
  inline void zero(){
    x = y = z = 0.0;
  }
  inline double dot( const TMRPoint& p ) const {
    return x*p.x + y*p.y + z*p.z;
  }
  inline double dot( const TMRPoint* p ) const {
    return x*p->x + y*p->y + z*p->z;
  }
  
  double x, y, z;
};

/*
  The MPI TMROctant data type
*/
extern MPI_Datatype TMROctant_MPI_type;
extern MPI_Datatype TMRQuadrant_MPI_type;
extern MPI_Datatype TMRPoint_MPI_type;
extern MPI_Datatype TMRIndexWeight_MPI_type;

// Initialize and finalize the data type
void TMRInitialize();
int TMRIsInitialized();
void TMRFinalize();

/*
  The following class is used to help create the interpolation and
  restriction operators. It stores both the node index and
  interpolation/restriction weight values in the same structure.
*/
class TMRIndexWeight {
 public:
  // TMRIndexWeight( const TMRIndexWeight& in ){
  //   index = in.index;  weight = in.weight;
  // }
  int index;
  double weight;

  // Sort and uniquify a list of indices and weights
  static int uniqueSort( TMRIndexWeight *array, int size ){
    qsort(array, size, sizeof(TMRIndexWeight), 
	  TMRIndexWeight::compare_weights);
    
    // Now that the weights are sorted, remove duplicates by adding up
    // the weights
    int j = 0; // Location to place entries
  
    for ( int i = 0; i < size; i++, j++ ){
      // Copy over the newest weight/index data
      if (i != j){
        array[j] = array[i];
      }
      
      // While the index is the same as the next, keep
      // adding the weight
      while ((i < size-1) && 
             (array[i].index == array[i+1].index)){
        array[j].weight += array[i+1].weight;
        i++;
      }
    }
    
    // Return the new size of the array
    return j;
  }

 private:
  // Compare two TMRIndexWeight objects - compare them based on
  // their index value only.
  static int compare_weights( const void *a, const void *b ){
    const TMRIndexWeight *A = static_cast<const TMRIndexWeight*>(a);
    const TMRIndexWeight *B = static_cast<const TMRIndexWeight*>(b);
    
    return A->index - B->index;
  }
};

/*
  Reference counted TMR entity
*/
class TMREntity {
 public:
  TMREntity();
  virtual ~TMREntity();
  
  // Set an attribute
  void setName( const char *name );
  const char* getName() const;

  // Reference count the geometric entity objects
  // --------------------------------------------
  void incref();
  void decref();

  // Set/get the tolerances used within the geometric search algorithms
  // ------------------------------------------------------------------
  static void setTolerances( double _eps_dist, double _eps_cosine );
  static void getTolerances( double *_eps_dist, double *_eps_cosine );

  // Retrieve the entity identification number
  // -----------------------------------------
  int getEntityId() const { return entity_id; }

 protected:
  static double eps_dist;
  static double eps_cosine;

 private:
  // Total reference count for this object
  int ref_count;

  // Named attribute
  char *name;

  // The entity identification value
  const int entity_id;
  static int entity_id_count;
};

#endif // TMR_BASE_H
