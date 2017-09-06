#ifndef TMR_TACS_CREATOR_H
#define TMR_TACS_CREATOR_H

/*
  The following file contains objects that facilitate the creation of
  TACSAssembler ojbects. These objects automate setting connectivity,
  node locations and boundary conditions in parallel for both forests
  of quadtrees and forests of octrees.

  These are virtual base classes. It is necessary to override the
  createElement() and optionally the createAuxElement()
  functions. These functions create the appropriate TACS element
  objects needed for TACSAssembler.
*/

#include "TMRQuadForest.h"
#include "TMROctForest.h"
#include "TACSAssembler.h"

/*
  Specify a list of boundary conditions through a list of attributes
  and boundary condition information for each node
*/
class TMRBoundaryConditions : public TMREntity {
 public:
  TMRBoundaryConditions();
  ~TMRBoundaryConditions();

  // Add a boundary condition associated with the specified attribute  
  void addBoundaryCondition( const char *attr, 
                             int num_bcs, const int bc_nums[],
                             const TacsScalar *_bc_vals=NULL );

  // Get the number of boundary conditions
  int getNumBoundaryConditions();
  void getBoundaryCondition( int bc, const char **_attr, 
                             int *_num_bcs,
                             const int **_bcs_nums,
                             const TacsScalar **_bc_vals );

 public:
  // The number of boundary conditions
  int num_bcs;

  // Linked list sub-clas for the boundary conditions -- there are
  // typically only a handful of these objects
  class BCNode {
  public:
    BCNode( const char *_attr, int _num_bcs, const int *_bc_nums,
            const TacsScalar *_bc_vals );
    ~BCNode();
    BCNode *next;
    char *attr;
    int num_bcs;
    int *bc_nums;
    TacsScalar *bc_vals;
  } *bc_root, *bc_current;
};

/*
  The creator object for quadrilateral meshes

  This sets the quads/elements into a new TACSAssembler object
*/
class TMRQuadTACSCreator : public TMREntity {
 public:
  TMRQuadTACSCreator( TMRBoundaryConditions *_bcs );
  virtual ~TMRQuadTACSCreator();

  // Vritual function to create the connectivity (default is usually fine)
  virtual void createConnectivity( int order,
                                   TMRQuadForest *forest,
                                   int **_conn, int **_ptr,
                                   int *_num_elements );

  // Create an array of elements for the given forest
  virtual void createElements( int order,
                               TMRQuadForest *forest,
                               int num_elements,
                               TACSElement **elements ) = 0;

  // Create any auxiliary element for the given quadrant
  virtual TACSAuxElements *createAuxElements( int order,
                                              TMRQuadForest *forest ){
    return NULL;
  }

  // Create the TACSAssembler object with the given order for this forest
  TACSAssembler *createTACS( int order, TMRQuadForest *forest,
                             TacsScalar _scale=1.0 );

 private:
  // Set the boundary conditions
  void setBoundaryConditions( TMRQuadForest *forest,
                              TACSAssembler *tacs );

  // Set the node locations
  void setNodeLocations( TMRQuadForest *forest, 
                         TACSAssembler *tacs,
                         TacsScalar _scale=1.0);
  
  TMRBoundaryConditions *bcs;
};

/*
  The creator object for octant meshes

  This sets octants/elements into a new TACSAssembler object
*/
class TMROctTACSCreator : public TMREntity {
 public:
  TMROctTACSCreator( TMRBoundaryConditions *_bcs );
  virtual ~TMROctTACSCreator();

  // Vritual function to create the connectivity (default is usually fine)
  virtual void createConnectivity( int order,
                                   TMROctForest *forest,
                                   int **_conn, int **_ptr,
                                   int *_num_elements );

  // Create an array of elements for the given forest
  virtual void createElements( int order,
                               TMROctForest *forest,
                               int num_elements,
                               TACSElement **elements ) = 0;

  // Create any auxiliary element for the given quadrant
  virtual TACSAuxElements *createAuxElements( int order,
                                              TMROctForest *forest ){
    return NULL;
  }

  // Add a boundary condition associated with the specified attribute  
  void addBoundaryCondition( const char *attr, 
                             int num_bcs, const int bc_nums[] );

  // Create the TACSAssembler object with the given order for this forest
  TACSAssembler *createTACS( int order, TMROctForest *forest,
                             TacsScalar _scale=1.0);

 private:
  // Set the boundary conditions
  void setBoundaryConditions( TMROctForest *forest,
                              TACSAssembler *tacs );

  // Set the node locations
  void setNodeLocations( TMROctForest *forest, 
                         TACSAssembler *tacs,
                         TacsScalar _scale=1.0);

  TMRBoundaryConditions *bcs;
};

#endif // TMR_TACS_CREATOR
