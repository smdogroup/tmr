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
  Specify a list of boundary conditions through a list of names
  and boundary condition information for each node
*/
class TMRBoundaryConditions : public TMREntity {
 public:
  TMRBoundaryConditions();
  ~TMRBoundaryConditions();

  // Add a boundary condition associated with the specified name  
  void addBoundaryCondition( const char *name, 
                             int num_bcs, const int bc_nums[],
                             const TacsScalar *_bc_vals=NULL );

  // Get the number of boundary conditions
  int getNumBoundaryConditions();
  void getBoundaryCondition( int bc, const char **_name, 
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
    BCNode( const char *_name, int _num_bcs, const int *_bc_nums,
            const TacsScalar *_bc_vals );
    ~BCNode();
    BCNode *next;
    char *name;
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
  TMRQuadTACSCreator( TMRBoundaryConditions *_bcs,
                      int _design_vars_per_node=1,
                      TMRQuadForest *_filter=NULL );
  TMRQuadTACSCreator();
  virtual ~TMRQuadTACSCreator();

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
  TACSAssembler *createTACS( TMRQuadForest *forest,
                             TACSAssembler::OrderingType 
                               ordering=TACSAssembler::NATURAL_ORDER );

  TMRQuadForest* getFilter(){
    return filter;
  }

 protected:
  // Initialize the data
  void initialize( TMRBoundaryConditions *_bcs,
                   int _design_vars_per_node,
                   TMRQuadForest *_filter );

  // Set the boundary conditions
  void setBoundaryConditions( TMRQuadForest *forest,
                              TACSAssembler *tacs );

  // Set the node locations
  void setNodeLocations( TMRQuadForest *forest, 
                         TACSAssembler *tacs );
  
  TMRBoundaryConditions *bcs;
  int design_vars_per_node;
  TMRQuadForest *filter;
};

/*
  The creator object for octant meshes

  This sets octants/elements into a new TACSAssembler object
*/
class TMROctTACSCreator : public TMREntity {
 public:
  TMROctTACSCreator( TMRBoundaryConditions *_bcs,
                     int _design_vars_per_node=1,
                     TMROctForest *_filter=NULL  );
  TMROctTACSCreator();
  virtual ~TMROctTACSCreator();

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

  // Add a boundary condition associated with the specified name  
  void addBoundaryCondition( const char *name, 
                             int num_bcs, const int bc_nums[] );

  // Create the TACSAssembler object with the given order for this forest
  TACSAssembler *createTACS( TMROctForest *forest,
                             TACSAssembler::OrderingType 
                               ordering=TACSAssembler::NATURAL_ORDER );

  TMROctForest* getFilter(){
    return filter;
  }

 protected:
  // Initialize the data
  void initialize( TMRBoundaryConditions *_bcs,
                   int _design_vars_per_node,
                   TMROctForest *_filter );

  // Set the boundary conditions
  void setBoundaryConditions( TMROctForest *forest,
                              TACSAssembler *tacs );

  // Set the node locations
  void setNodeLocations( TMROctForest *forest, 
                         TACSAssembler *tacs );

  TMRBoundaryConditions *bcs;
  int design_vars_per_node;
  TMROctForest *filter;
};

#endif // TMR_TACS_CREATOR
