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

#endif // TMR_CY_CREATOR_H