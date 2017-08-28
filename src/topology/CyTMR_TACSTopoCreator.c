#include <string.h>
#include "CyTMR_TACSTopoCreator.h"

/*
  Copyright (c) 2017 Graeme Kennedy. All rights reserved
*/

/*
  The constructor for the TMR_TACSTopoCreator wrapper
*/ 

CyTMROctTACSTopoCreator:: CyTMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                                                   TMRStiffnessProperties _properties,
                                                   TMROctForest *_filter,
                                                   const char *shell_attr, 
                                                   SolidShellWrapper *_shell ):
TMROctTACSTopoCreator(_bcs,_properties,_filter,shell_attr,_shell){
  // Set the initial values for the callbacks
  self = NULL; 
  createelements = NULL;
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
  createelements = func;
}

void CyTMROctTACSTopoCreator::createElements( int order,
                                              TMROctForest *forest,
                                              int num_elements,
                                              TACSElement **elements ){
  if (!createelements){
    fprintf(stderr, "createelements callback not defined\n");
    return;
  }
  createelements(self, order, forest, 
                 num_elements, elements);
}
