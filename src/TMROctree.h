#ifndef TMR_OCTANT_TREE_H
#define TMR_OCTANT_TREE_H

#include <stdio.h>
#include "TMROctant.h"

/*
  The main octree class. 

  This is used to create, balance, and coarsen an octree as well as
  create a mesh from the octree data. This data can be used to
  construct meshes in TACS (or other finite-element codes). The octree
  is well-suited for creating adaptive meshes for topology
  optimization problems on regular domains.

  The octree object provides a means to construct inter-grid
  operators. The inter-grid operations are constructed from the nodes
  on the coarse mesh. For each node, we store the largest adjoining
  element level. Taking the next-nearest neighbours for the largest
  adjoining element length gives an intermesh operator that is
  guaranteed to be the independent node set. The regular full
  weighting multigrid scheme can be used. 

  Note that these meshes are 2-1 balanced, otherwise this type of
  construction scheme could not be used. There will be at most a
  single level of refinement difference between adjacent cells.
*/
class TMROctree {
 public:
  TMROctree( int refine_level );
  TMROctree( int nrand, int min_level, int max_level );
  TMROctree( TMROctantArray *_list );
  ~TMROctree();

  // Refine the octree
  // -----------------
  void refine( const int refinement[],
	       int min_level=0, int max_level=TMR_MAX_LEVEL );

  // Balance the tree to ensure 2-1 balancing
  // ----------------------------------------
  void balance( int balance_corner=0 );

  // Coarsen the tree uniformly
  // --------------------------
  TMROctree *coarsen();

  // Find an octant that completely encloses the provided octant
  // -----------------------------------------------------------
  TMROctant *findEnclosing( TMROctant *oct );
  void findEnclosingRange( TMROctant *oct,
			   int *low, int *high );

  // Allocate the local nodes (but do not order them)
  // ------------------------------------------------
  void createNodes( int _order );

  // Print a representation of the tree to a file
  // --------------------------------------------
  void printOctree( const char *filename );

  // Retrieve the octree elements
  // ----------------------------
  void getElements( TMROctantArray **_elements ){
    if (_elements){ *_elements = elements; }
  }
  void setElements( TMROctantArray *_elements ){
    if (elements){ delete elements; }
    elements = _elements;
  }

  // Retrieve the octant nodes
  // -------------------------
  void getNodes( TMROctantArray **_nodes ){
    if (_nodes){ *_nodes = nodes; }
  }
  void setNodes( TMROctantArray *_nodes ){
    if (nodes){ delete nodes; }
    nodes = _nodes;
  }

  // Quickly retrieve the number of nodes/elements
  // ---------------------------------------------
  int getNumNodes(){
    int size = 0;
    if (nodes){ nodes->getArray(NULL, &size); }
    return size; 
  }
  int getNumElements(){
    int size = 0;
    if (elements){ elements->getArray(NULL, &size); }
    return size; 
  }

  // Retrieve the array of points - created concurrently with the nodes
  // ------------------------------------------------------------------
  int getPoints( TMRPoint **_X ){
    int size = 0;
    if (_X){
      *_X = X;
      if (nodes){ nodes->getArray(NULL, &size); }
    }
    return size;
  }

 private:
  // The list octants representing the elements
  TMROctantArray *elements; 

  // The nodes within the element mesh
  TMROctantArray *nodes; 

  // The node locations
  TMRPoint *X;

  // Store the order of the mesh
  int order;
};

#endif // TMR_OCTANT_TREE_H
