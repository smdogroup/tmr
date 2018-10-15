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

#ifndef TMR_CY_CREATOR_H
#define TMR_CY_CREATOR_H

#include "TMR_TACSCreator.h"
#include "TMR_TACSTopoCreator.h"

/*
  This is a light-weight wrapper for creating elements in python
*/
class TMRCyQuadCreator : public TMRQuadTACSCreator {
 public:
  TMRCyQuadCreator( TMRBoundaryConditions *_bcs ):
    TMRQuadTACSCreator(_bcs){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateQuadElement( TACSElement* (*func)(void*, int, TMRQuadrant*) ){
    createquadelement = func;
  }

  void createElements( int order,
                       TMRQuadForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    // Get the array of quadrants
    int size;
    TMRQuadrant *array;
    TMRQuadrantArray *quadrants;
    forest->getQuadrants(&quadrants);
    quadrants->getArray(&array, &size);

    // Set the element types into the matrix
    memset(elements, 0, num_elements*sizeof(TACSElement*));
    for ( int i = 0; i < num_elements; i++ ){
      TACSElement *elem =
        createquadelement(self, order, &array[i]);
      if (!elem){
        fprintf(stderr, "TMRCyQuadCreator error: Element not created\n");
      }
      else {
        elements[i] = elem;
      }
    }
  }

 private:
  void *self;
  TACSElement* (*createquadelement)( void*, int, TMRQuadrant* );
};

/*
  This is a light-weight wrapper for creating elements in python
*/
class TMRCyOctCreator : public TMROctTACSCreator {
 public:
  TMRCyOctCreator( TMRBoundaryConditions *_bcs ):
    TMROctTACSCreator(_bcs){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateOctElement( TACSElement* (*func)(void*, int, TMROctant*) ){
    createoctelement = func;
  }

  void createElements( int order,
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    // Get the array of octants
    int size;
    TMROctant *array;
    TMROctantArray *octants;
    forest->getOctants(&octants);
    octants->getArray(&array, &size);

    // Set the element types into the matrix
    memset(elements, 0, num_elements*sizeof(TACSElement*));
    for ( int i = 0; i < num_elements; i++ ){
      TACSElement *elem =
        createoctelement(self, order, &array[i]);
      if (!elem){
        fprintf(stderr, "TMRCyOctCreator error: Element not created\n");
      }
      else {
        elements[i] = elem;
      }
    }
  }

 private:
  void *self;
  TACSElement* (*createoctelement)( void*, int, TMROctant* );
};

/*
  Create a wrapper for topology optimization with a filter
*/
class TMRCyTopoQuadCreator : public TMRQuadTACSTopoCreator {
 public:
  TMRCyTopoQuadCreator( TMRBoundaryConditions *_bcs,
                        TMRQuadForest *_filter ):
  TMRQuadTACSTopoCreator(_bcs, _filter){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateQuadTopoElement( 
    TACSElement* (*func)(void*, int, TMRQuadrant*, TMRIndexWeight*, int) ){
    createquadtopoelement = func;
  }

  // Create the element
  TACSElement *createElement( int order, 
                              TMRQuadrant *quad,
                              TMRIndexWeight *weights, 
                              int nweights ){
    TACSElement *elem =
      createquadtopoelement(self, order, quad, weights, nweights);
    return elem;
  }

 private:
  void *self; // Pointer to the python-level object
  TACSElement* (*createquadtopoelement)( 
    void*, int, TMRQuadrant*, TMRIndexWeight *weights, int nweights );
};

/*
  Create a wrapper for topology optimization with a filter
*/
class TMRCyTopoOctCreator : public TMROctTACSTopoCreator {
 public:
  TMRCyTopoOctCreator( TMRBoundaryConditions *_bcs,
                       TMROctForest *_filter ):
  TMROctTACSTopoCreator(_bcs, _filter){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateOctTopoElement( 
    TACSElement* (*func)(void*, int, TMROctant*, TMRIndexWeight*, int) ){
    createocttopoelement = func;
  }

  // Create the element
  TACSElement *createElement( int order, 
                              TMROctant *oct,
                              TMRIndexWeight *weights, 
                              int nweights ){
    TACSElement *elem =
      createocttopoelement(self, order, oct, weights, nweights);
    return elem;
  }

 private:
  void *self; // Pointer to the python-level object
  TACSElement* (*createocttopoelement)( 
    void*, int, TMROctant*, TMRIndexWeight *weights, int nweights );
};

/*
  Create a wrapper for topology optimization with a filter
*/
class TMRCyTopoQuadBernsteinCreator : public TMRQuadBernsteinTACSTopoCreator {
 public:
  TMRCyTopoQuadBernsteinCreator( TMRBoundaryConditions *_bcs,
                                 TMRQuadForest *_forest ):
  TMRQuadBernsteinTACSTopoCreator(_bcs, _forest){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateQuadTopoElement( 
                                TACSElement* (*func)(void*, int, TMRQuadrant*,
                                                     int*, int, TMRQuadForest*) ){
    createquadtopoelement = func;
  }

  // Create the element
  TACSElement *createElement( int order, 
                              TMRQuadrant *quad,
                              int *index, 
                              int nweights,
                              TMRQuadForest *filter){
    TACSElement *elem =
      createquadtopoelement(self, order, quad, index, nweights,
                            filter);
    return elem;
  }

 private:
  void *self; // Pointer to the python-level object
  TACSElement* (*createquadtopoelement)( void*, int, TMRQuadrant*, int *index,
                                         int nweights, TMRQuadForest* );
};

#endif // TMR_CY_CREATOR_H
