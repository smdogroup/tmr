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

#ifndef TMR_TOPO_FILTER_H
#define TMR_TOPO_FILTER_H

#include "TMRBase.h"
#include "TMRQuadForest.h"
#include "TMROctForest.h"
#include "TACSAssembler.h"

/*
  Abstract base class for the filter problem
*/
class TMRTopoFilter : public TMREntity {
 public:
  // Get the TACSAssembler instance (on the finest mesh level)
  virtual TACSAssembler* getAssembler() = 0;

  // Get the Quad or OctForest on the finest mesh level
  virtual TMRQuadForest* getFilterQuadForest() = 0;
  virtual TMROctForest* getFilterOctForest() = 0;

  // Set the design variables
  virtual void setDesignVars( TACSBVec *vec ) = 0;

  // Get the unfiltered design variables
  virtual void getDesignVars( TACSBVec **vec ){
    fprintf(stderr, "TMR Filter Error: getDesignVars() method is not implemented!\n");
  };

  // Set values/add values to the vector
  virtual void addValues( TACSBVec *vec ) = 0;
  virtual void setValues( TACSBVec *vec ) = 0;

  // Write the STL file
  virtual void writeSTLFile( int k, double cutoff, const char *filename ){}

  // Apply filter/filter transpose to some vector that has same size as design variable
  virtual void applyFilter( TACSBVec *in, TACSBVec *out ) {
    fprintf(stderr, "TMR Filter Error: applyFilter() method is not implemented!\n");
  }
  virtual void applyTranspose( TACSBVec *in, TACSBVec *out ) {
    fprintf(stderr, "TMR Filter Error: applyTranspose() method is not implemented!\n");
  }
};

#endif // TMR_TOPO_FILTER_H
