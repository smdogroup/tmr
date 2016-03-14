#ifndef TMR_OCTANT_TREE_H
#define TMR_OCTANT_TREE_H

#include <stdio.h>
#include "TMROctant.h"

/*
  Assemble a single octree
*/
class TMROctree {
 public:
  TMROctree( int refine_level );
  TMROctree( int nrand, int min_level, int max_level );
  TMROctree( TMROctantArray *_list, int _is_sorted );
  ~TMROctree();

  // Balance the tree to ensure 2-1 balancing
  // ----------------------------------------
  void balance( int balance_type = 3 );  

  // Coarsen the tree uniformly
  // --------------------------
  TMROctree *coarsen();

  // Create the connectivity for the mesh
  // ------------------------------------
  void createNodes( int order );

  // Get the connectivity from the mesh
  // ----------------------------------
  void createMesh( int );

  // Print a representation of the tree to a file
  // --------------------------------------------
  void printTree( const char * filename );

 private:
  // The list octants representing the elements
  TMROctantArray *elements; 
  int is_sorted; // Is the element list sorted?

  // The nodes within the element mesh
  TMROctantArray *nodes; 

  // The number of elements and independent and dependent nodes
  int num_elements;
  int num_nodes, num_dependent_nodes;
};

#endif // TMR_OCTANT_TREE_H
