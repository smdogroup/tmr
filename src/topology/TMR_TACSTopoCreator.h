#ifndef TMR_TACS_TOPO_CREATOR_H
#define TMR_TACS_TOPO_CREATOR_H

#include "TMR_TACSCreator.h"
#include "TACSAssembler.h"
#include "TMROctStiffness.h"

/*
  This is an abstract base class used to create octforests specialized
  for topology optimization. This class can be overriden with the 
  createElement class to create different types of topolgy optimization
  problems for 3D structures
*/
class TMROctTACSTopoCreator : public TMROctTACSCreator {
 public:
  TMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                         TMROctForest *_filter );
  ~TMROctTACSTopoCreator();

  // Create the elements
  void createElements( int order, 
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements );

  // Create the element
  virtual TACSElement *createElement( int order, 
                                      TMROctant *oct,
                                      TMRIndexWeight *weights, 
                                      int nweights ) = 0;

  // Get the underlying objects that define the filter
  void getFilter( TMROctForest **filter );
  void getMap( TACSVarMap **_map );
  void getIndices( TACSBVecIndices **_indices );  

 private:
  // Compute the weights for a given point
  void computeWeights( const int mesh_order, const double *knots,
                       TMROctant *node, TMROctant *oct,
                       TMRIndexWeight *weights, double *tmp );

  // The forest that defines the filter
  TMROctForest *filter;

  // The filter map for this object. This defines how the design
  // variables are distributed across all of the processors.
  TACSVarMap *filter_map;

  // The filter indices. This defines the relationship between the
  // local design variable numbers and the global design variable
  // numbers.
  TACSBVecIndices *filter_indices;
};

/*
  This class is an abstract base class for setting up 
  topology optimzation problems using TMRQuadForest objects
*/
class TMRQuadTACSTopoCreator : public TMRQuadTACSCreator {
 public:
  TMRQuadTACSTopoCreator( TMRBoundaryConditions *_bcs,
                          TMRQuadForest *_filter );
  ~TMRQuadTACSTopoCreator();

  // Create the elements
  void createElements( int order,
                       TMRQuadForest *forest,
                       int num_elements,
                       TACSElement **elements );

  // Create the element
  virtual TACSElement *createElement( int order, 
                                      TMRQuadrant *oct,
                                      TMRIndexWeight *weights, 
                                      int nweights ) = 0;

  // Get the underlying objects that define the filter
  void getFilter( TMRQuadForest **filter );
  void getMap( TACSVarMap **_map );
  void getIndices( TACSBVecIndices **_indices );  

 private:
  // Compute the weights for a given point
  void computeWeights( const int mesh_order, const double *knots,
                       TMRQuadrant *node, TMRQuadrant *quad,
                       TMRIndexWeight *weights, double *tmp );

  // The forest that defines the filter
  TMRQuadForest *filter;

  // The filter map for this object. This defines how the design
  // variables are distributed across all of the processors.
  TACSVarMap *filter_map;

  // The filter indices. This defines the relationship between the
  // local design variable numbers and the global design variable
  // numbers.
  TACSBVecIndices *filter_indices;
};

#endif // TMR_TACS_TOPO_CREATOR_H

