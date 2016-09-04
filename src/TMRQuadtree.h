#ifndef TMR_QUADTREE_H
#define TMR_QUADTREE_H

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

#include "TMRQuadrant.h"

/*
  The main quadtree class. 

  This is used to create, balance, and coarsen an quadtree as well as
  create a mesh from the quadtree data. This data can be used to
  construct meshes in TACS (or other finite-element codes). The quadtree
  is well-suited for creating adaptive meshes for topology
  optimization problems on regular domains.

  The quadtree object provides a means to construct inter-grid
  operators. The inter-grid operations are constructed from the nodes
  on the coarse mesh. For each node, we store the largest adjoining
  element level. Taking the next-nearest neighbours for the largest
  adjoining element length gives an intermesh operator that is
  guaranteed to be the independent node set. The regular full
  weighting multigrid scheme can be used. The fine and coarse meshes
  are given as follows:

  Fine:                         Coarse:

  o -- o -- o -- o -- o         o ------- o ------- o
  |    |    |    |    |         |         |         |
  o -- o -- o -- o -- o         |         |         |
  |    |    |    |    |         |         |         |
  o -- o -- o -- o -- o         o ------- o ------- o
  |         |    |    |         |         |         |
  |         o -- o -- o         |         |         |
  |         |    |    |         |         |         |
  o ------- o -- o -- o         o ------- o ------- o
  |         |         |         |                   |
  |         |         |         |                   |
  |         |         |         |                   |
  o ------- o ------- o         |                   |
  |         |         |         |                   |
  |         |         |         |                   |
  |         |         |         |                   |
  o ------- o ------- o         o ----------------- o

  Note that these meshes are 2-1 balanced, otherwise this type of
  construction scheme could not be used. There will be at most a
  single level of refinement difference between adjacent cells.
*/
class TMRQuadtree {
 public:
  TMRQuadtree( int refine_level );
  TMRQuadtree( int nrand, int min_level, int max_level );
  TMRQuadtree( TMRQuadrantArray *_list );
  ~TMRQuadtree();

  // Refine the quadtree
  // -----------------
  void refine( const int refinement[]=NULL,
	       int min_level=0, int max_level=TMR_MAX_LEVEL );

  // Duplicate or coarsen the tree
  // -----------------------------
  TMRQuadtree *coarsen();

  // Find an quadrant that completely encloses the provided quadrant
  // ---------------------------------------------------------------
  TMRQuadrant *findEnclosing( TMRQuadrant *quad );
  void findEnclosingRange( TMRQuadrant *quad,
			   int *low, int *high );

  // Order the nodes but do not create the connectivity
  // --------------------------------------------------
  void createNodes( int _order );

  // Print a representation of the tree to a file
  // --------------------------------------------
  void printQuadtree( const char * filename );

  // Retrieve the quadtree elements
  // ------------------------------
  void getElements( TMRQuadrantArray **_elements ){
    if (_elements){ *_elements = elements; }
  }
  void setElements( TMRQuadrantArray *_elements ){
    if (elements){ delete elements; }
    elements = _elements;
  }

  // Retrieve the quadtree nodes
  // ---------------------------
  void getNodes( TMRQuadrantArray **_nodes ){
    if (_nodes){ *_nodes = nodes; }
  }
  void setNodes( TMRQuadrantArray *_nodes ){
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
  
 private:
  // The list quadrants representing the elements
  TMRQuadrantArray *elements; 

  // The nodes within the element mesh
  TMRQuadrantArray *nodes; 

  // Store the order of the mesh
  int order;
};

#endif // TMR_QUADTREE_H
