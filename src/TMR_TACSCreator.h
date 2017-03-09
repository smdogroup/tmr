#ifndef TMR_TACS_CREATOR_H
#define TMR_TACS_CREATOR_H

#include "TMRQuadForest.h"
#include "TACSAssembler.h"

/*
  Get the element and the auxiliary elements
*/
class TMRQuadTACSCreator : public TMREntity {
 public:
  TMRQuadTACSCreator();
  virtual ~TMRQuadTACSCreator();

  virtual TACSElement *createElement( int order,
                                      TMRQuadForest *forest,
                                      TMRQuadrant quad ) = 0;
  virtual TACSElement *createAuxElement( int order,
                                         TMRQuadForest *forest,
                                         TMRQuadrant quad ){
    return NULL;
  }

  // Add a boundary condition associated with the specified attribute  
  // ----------------------------------------------------------------
  void addBoundaryCondition( const char *attr, 
                             int num_bcs, const int bc_nums[] );

  // Create the TACSAssembler object with the given order for this forest
  // --------------------------------------------------------------------
  TACSAssembler *createTACS( int order, TMRQuadForest *forest );

 private:
  // Set the boundary conditions
  void setBoundaryConditions( TMRQuadForest *forest,
                              TACSAssembler *tacs );

  // Set the node locations
  void setNodeLocations( TMRQuadForest *forest, 
                         TACSAssembler *tacs );

  // Linked list sub-clas for the boundary conditions -- there are
  // typically only a handful of these objects
  class BCNode {
  public:
    BCNode( const char *_attr, int _num_bcs, const int *_bc_nums );
    ~BCNode();
    BCNode *next;
    char *attr;
    int num_bcs;
    int *bc_nums;
  } *bc_root, *bc_current;
};

#endif // TMR_TACS_CREATOR
