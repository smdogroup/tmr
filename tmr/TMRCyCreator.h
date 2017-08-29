#ifndef TMR_CY_CREATOR_H
#define TMR_CY_CREATOR_H

#include "TMR_TACSCreator.h"

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

#endif // TMR_CY_CREATOR_H