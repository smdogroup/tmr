#include <string.h>
#include "CyTMR_TACSTopoCreator.h"

/*
  Copyright (c) 2017 Graeme Kennedy. All rights reserved
*/

/*
  The constructor for the TMR_TACSTopoCreator wrapper
*/ 

CyTMROctTACSTopoCreator:: CyTMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                                                   TMROctStiffnessProperties _properties,
                                                   TMROctForest *_filter,
                                                   const char *shell_attr, 
                                                   SolidShellWrapper *_shell ):
TMROctTACSTopoCreator(_bcs,_properties,_filter,shell_attr,_shell){
  // Set the initial values for the callbacks
  self = NULL; 
  ocreateelements = NULL;
}

CyTMROctTACSTopoCreator::~CyTMROctTACSTopoCreator(){}

/*
  Set the member callback functions that are required
*/
void CyTMROctTACSTopoCreator::setSelfPointer( void *_self ){
  self = _self;
}

void CyTMROctTACSTopoCreator::setCreateElements( void (*func)(void*, int, TMROctForest*, 
                                                              int, TACSElement**)){
  ocreateelements = func;
}

void CyTMROctTACSTopoCreator::createElements( int order,
                                              TMROctForest *forest,
                                              int num_elements,
                                              TACSElement **elements ){
  if (!ocreateelements){
    fprintf(stderr, "createelements callback not defined\n");
    return;
  }
  ocreateelements(self, order, forest, 
                  num_elements, elements);
}

CyTMRQuadTACSTopoCreator:: CyTMRQuadTACSTopoCreator( TMRBoundaryConditions *_bcs,
                                                     TMRQuadStiffnessProperties _properties,
                                                     TMRQuadForest *_filter):
TMRQuadTACSTopoCreator(_bcs,_properties,_filter){
  // Set the initial values for the callbacks
  self = NULL; 
  qcreateelements = NULL;
}

CyTMRQuadTACSTopoCreator::~CyTMRQuadTACSTopoCreator(){}

/*
  Set the member callback functions that are required
*/
void CyTMRQuadTACSTopoCreator::setSelfPointer( void *_self ){
  self = _self;
}

void CyTMRQuadTACSTopoCreator::setCreateElements( void (*func)(void*, int, TMRQuadForest*, 
                                                               int, TACSElement**)){
  qcreateelements = func;
}

void CyTMRQuadTACSTopoCreator::createElements( int order,
                                               TMRQuadForest *forest,
                                               int num_elements,
                                               TACSElement **elements ){
  if (!qcreateelements){
    fprintf(stderr, "createelements callback not defined\n");
    return;
  }
  qcreateelements(self, order, forest, 
                  num_elements, elements);
}
