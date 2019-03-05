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
#include "TACSAssembler.h"

/*
  Abstract base class for the filter problem
*/
class TMRTopoFilter : public TMREntity {
 public:
  // Get the MPI communicator
  virtual MPI_Comm getMPIComm() = 0;
  
  // Get the TACSAssembler instance (on the finest mesh level)
  virtual TACSAssembler *getAssembler() = 0;
  
  // Get problem definitions maximum local size of the design variable values
  virtual int getVarsPerNode() = 0;
  virtual int getNumLocalVars() = 0;  
  virtual int getMaxNumLocalVars() = 0;
  
  // Create a design vector on the finest mesh level
  virtual TACSBVec *createVec() = 0;

  // Set the design variable values (including all local values) 
  virtual void setDesignVars( TACSBVec *x ) = 0;

  // Set values/add values to the vector
  virtual void addValues( TacsScalar *in, TACSBVec *out ) = 0;
  virtual void setValues( TacsScalar *in, TACSBVec *out ) = 0;
};

#endif // TMR_TOPO_FILTER_H
