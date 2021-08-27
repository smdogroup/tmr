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
  TMRCyQuadCreator( TMRBoundaryConditions *_bcs, int _design_vars_per_node=1,
                    TMRQuadForest *_filter=NULL ):
    TMRQuadTACSCreator(_bcs, _design_vars_per_node, _filter){}

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
  TMRCyOctCreator( TMRBoundaryConditions *_bcs, int _design_vars_per_node=1,
                   TMROctForest *_filter=NULL ):
    TMROctTACSCreator(_bcs, _design_vars_per_node, _filter){}

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
                        int _design_vars_per_node,
                        TMRQuadForest *_filter ):
  TMRQuadTACSTopoCreator(_bcs, _design_vars_per_node, _filter){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateQuadTopoElement( 
    TACSElement* (*func)(void*, int, TMRQuadrant*, int, TMRIndexWeight*) ){
    createquadtopoelement = func;
  }

  // Create the element
  TACSElement *createElement( int order, 
                              TMRQuadrant *quad,
                              int nweights,
                              TMRIndexWeight *weights ){
    TACSElement *elem =
      createquadtopoelement(self, order, quad, nweights, weights);
    return elem;
  }

 private:
  void *self; // Pointer to the python-level object
  TACSElement* (*createquadtopoelement)( 
    void*, int, TMRQuadrant*, int nweights, TMRIndexWeight *weights );
};

/*
  Create a wrapper for topology optimization with a filter
*/
class TMRCyTopoOctCreator : public TMROctTACSTopoCreator {
 public:
  TMRCyTopoOctCreator( TMRBoundaryConditions *_bcs,
                        int _design_vars_per_node,
                       TMROctForest *_filter ):
  TMROctTACSTopoCreator(_bcs, _design_vars_per_node, _filter){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateOctTopoElement( 
    TACSElement* (*func)(void*, int, TMROctant*, int, TMRIndexWeight*) ){
    createocttopoelement = func;
  }

  // Create the element
  TACSElement *createElement( int order, 
                              TMROctant *oct,
                              int nweights,
                              TMRIndexWeight *weights ){
    TACSElement *elem =
      createocttopoelement(self, order, oct, nweights, weights);
    return elem;
  }

 private:
  void *self; // Pointer to the python-level object
  TACSElement* (*createocttopoelement)( 
    void*, int, TMROctant*, int nweights, TMRIndexWeight *weights );
};

/*
  Create a wrapper for topology optimization with a filter
*/
class TMRCyTopoQuadConformCreator : public TMRQuadConformTACSTopoCreator {
 public:
  TMRCyTopoQuadConformCreator( TMRBoundaryConditions *_bcs,
                               int _design_vars_per_node,
                               TMRQuadForest *_forest,
                               int order=-1,
                               TMRInterpolationType interp_type=TMR_UNIFORM_POINTS ):
  TMRQuadConformTACSTopoCreator(_bcs, _design_vars_per_node, _forest, order, interp_type){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateQuadTopoElement( TACSElement* (*func)(void*, int, TMRQuadrant*,
                                                      int, const int*, TMRQuadForest*) ){
    createquadtopoelement = func;
  }

  // Create the element
  TACSElement *createElement( int order, 
                              TMRQuadrant *quad,
                              int nweights,
                              const int *index,
                              TMRQuadForest *fltr ){
    TACSElement *elem =
      createquadtopoelement(self, order, quad, nweights, index, fltr);
    return elem;
  }

 private:
  void *self; // Pointer to the python-level object
  TACSElement* (*createquadtopoelement)( void*, int, TMRQuadrant*, int nweights,
                                         const int *index, TMRQuadForest* );
};

/*
  Create a wrapper for topology optimization with a filter
*/
class TMRCyTopoOctConformCreator : public TMROctConformTACSTopoCreator {
 public:
  TMRCyTopoOctConformCreator( TMRBoundaryConditions *_bcs,
                              int _design_vars_per_node,
                              TMROctForest *_forest,
                              int order=-1,
                              TMRInterpolationType interp_type=TMR_UNIFORM_POINTS ):
  TMROctConformTACSTopoCreator(_bcs, _design_vars_per_node, _forest, order, interp_type){}

  void setSelfPointer( void *_self ){
    self = _self;
  }
  void setCreateOctTopoElement( 
                               TACSElement* (*func)(void*, int, TMROctant*,
                                                    int, const int*, TMROctForest*) ){
    createocttopoelement = func;
  }

  // Create the element
  TACSElement *createElement( int order, 
                              TMROctant *oct,
                              int nweights,
                              const int *index,
                              TMROctForest *fltr ){
    TACSElement *elem =
      createocttopoelement(self, order, oct, nweights, index, fltr);
    return elem;
  }

 private:
  void *self; // Pointer to the python-level object
  TACSElement* (*createocttopoelement)( void*, int, TMROctant*, int nweights,
                                        const int *index, TMROctForest* );
};

#endif // TMR_CY_CREATOR_H
