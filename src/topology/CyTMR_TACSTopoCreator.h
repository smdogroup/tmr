#ifndef CYTHON_TMR_TACS_TOPO_CREATOR_H
#define CYTHON_TMR_TACS_TOPO_CREATOR_H
/*
  Copyright (c) 2017 Graeme Kennedy. All rights reserved
*/

#include "TMR_TACSTopoCreator.h"

/*
  This code implements a simplifed interface for the TMROctTACSTopoCreator
  where we enable the specification of the element type
*/
class CyTMROctTACSTopoCreator : public TMROctTACSTopoCreator {
 public:
  CyTMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                           TMRStiffnessProperties _properties,
                           TMROctForest *_filter,
                           const char *shell_attr=NULL, 
                           SolidShellWrapper *_shell=NULL );
  ~CyTMROctTACSTopoCreator();

  // Set the member callback functions that are required
  // ---------------------------------------------------
  void setSelfPointer( void *_self );
  void setCreateElements( void (*func)(void*, int, TMROctForest*, 
                                       int, TACSElement**) );

  void createElements( int order,
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements );
 private:
  // Public member function pointers to the callbacks that are
  // required before class can be used
  // ---------------------------------------------------------
  void *self;
  void (*createelements)( void *self, int order, TMROctForest *forest,
                          int num_elements, TACSElement** elements);  
};
/* class CyTMRQuadTACSTopoCreator : public TMRQuadTACSTopoCreator { */
/*   CyTMROctTACSTopoCreator( TMRBoundaryConditions *_bcs, */
/*                            TMRStiffnessProperties _properties, */
/*                            TMROctForest *_filter, */
/*                            const char *shell_attr=NULL,  */
/*                            SolidShellWrapper *_shell=NULL ); */
/*   ~CyTMROctTACSTopoCreator(); */

/*   // Set the member callback functions that are required */
/*   // --------------------------------------------------- */
/*   void setSelfPointer( void *_self ); */
/*   void setCreateElements( void (*func)(void*, int, TMROctForest*,  */
/*                                        int, TACSElement**) ); */

/*   void createElements( int order, */
/*                        TMROctForest *forest, */
/*                        int num_elements, */
/*                        TACSElement **elements ); */
/*  private: */
/*   // Public member function pointers to the callbacks that are */
/*   // required before class can be used */
/*   // --------------------------------------------------------- */
/*   void *self; */
/*   void (*createelements)( void *self, int order, TMROctForest *forest, */
/*                           int num_elements, TACSElement** elements);   */
/* }; */
#endif // CYTHON_TMR_TACS_TOPO_CREATOR_H
