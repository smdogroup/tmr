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

#ifndef TACS_MATRIX_FILTER_CREATOR_H
#define TACS_MATRIX_FILTER_CREATOR_H

#include "TMR_TACSCreator.h"
#include "TACSQuadBasis.h"
#include "TACSHexaBasis.h"
#include "TACSQuadBernsteinBasis.h"
#include "TACSHexaBernsteinBasis.h"
#include "TACSElement2D.h"
#include "TACSElement3D.h"

/*
  Matrix filter creator classes
*/
class TMRQuadTACSMatrixCreator : public TMRQuadTACSCreator {
 public:
  TMRQuadTACSMatrixCreator( TACSElementModel *_model ):
    TMRQuadTACSCreator(NULL){
    model = _model;
    model->incref();
  }
  ~TMRQuadTACSMatrixCreator(){
    model->decref();
  }
  void createElements( int order,
                       TMRQuadForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    TACSElementBasis *basis = NULL;
    if (order == 2){
      basis = new TACSLinearQuadBasis();
    }
    else if (order == 3){
      basis = new TACSQuadraticQuadBernsteinBasis();
    }
    else if (order == 4){
      basis = new TACSCubicQuadBernsteinBasis();
    }
    else if (order == 5){
      basis = new TACSQuarticQuadBernsteinBasis();
    }
    else if (order == 6){
      basis = new TACSQuinticQuadBernsteinBasis();
    }

    TACSElement *elem = new TACSElement2D(model, basis);

    for ( int i = 0; i < num_elements; i++ ){
      elements[i] = elem;
    }
  }
 private:
  TACSElementModel *model;
};

class TMROctTACSMatrixCreator : public TMROctTACSCreator {
 public:
  TMROctTACSMatrixCreator( TACSElementModel *_model ):
    TMROctTACSCreator(NULL){
    model = _model;
    model->incref();
  }
  ~TMROctTACSMatrixCreator(){
    model->decref();
  }
  void createElements( int order,
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    TACSElementBasis *basis = NULL;
    if (order == 2){
      basis = new TACSLinearHexaBasis();
    }
    else if (order == 3){
      basis = new TACSQuadraticHexaBernsteinBasis();
    }
    else if (order == 4){
      basis = new TACSCubicHexaBernsteinBasis();
    }

    TACSElement *elem = new TACSElement3D(model, basis);

    for ( int i = 0; i < num_elements; i++ ){
      elements[i] = elem;
    }
  }
 private:
  TACSElementModel *model;
};

#endif // TACS_MATRIX_FILTER_CREATOR_H
