#ifndef TMR_TACS_TOPO_CREATOR_H
#define TMR_TACS_TOPO_CREATOR_H

#include "TMR_TACSCreator.h"
#include "TACSAssembler.h"
#include "TMROctStiffness.h"
#include "SolidShellWrapper.h"

class TMROctTACSTopoCreator : public TMROctTACSCreator {
 public:
  TMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                         TMRStiffnessProperties _properties,
                         TMROctForest *_filter,
                         const char *shell_attr=NULL, 
                         SolidShellWrapper *_shell=NULL );
  ~TMROctTACSTopoCreator();

  // Create the connectivity
  void createConnectivity( int order,
                           TMROctForest *forest,
                           int **_conn, int **_ptr,
                           int *_num_elements );

  // Create the elements
  void createElements( int order,
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements );

  // Get the underlying objects that define the filter
  void getForest( TMROctForest **filter );
  void getMap( TACSVarMap **_map );
  void getIndices( TACSBVecIndices **_indices );  

 private:
  // Compute the weights for a given point
  void computeWeights( TMROctant *oct, TMROctant *node,
                       TMRIndexWeight *welem );

  // The stiffness properties
  TMRStiffnessProperties properties;

  // The forest that defines the filter
  TMROctForest *filter;

  // The filter map for this object. This defines how the design
  // variables are distributed across all of the processors.
  TACSVarMap *filter_map;

  // The filter indices. This defines the relationship between the
  // local design variable numbers and the global design variable
  // numbers.
  TACSBVecIndices *filter_indices;

  // Set the top/bottom attributes
  char *shell_attr;
  SolidShellWrapper *shell;
};

#endif // TMR_TACS_TOPO_CREATOR_H
