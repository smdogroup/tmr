#ifndef TMR_QUADTREE_H
#define TMR_QUADTREE_H

#include <stdio.h>
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
  TMRQuadtree( int refine_level, 
               TMRQuadrant *_domain=NULL, int ndomain=0 );
  TMRQuadtree( int nrand, int min_level, int max_level );
  TMRQuadtree( TMRQuadrantArray *_list,
               TMRQuadrant *_domain=NULL, int ndomain=0 );
  ~TMRQuadtree();

  // Check if the provided quadrant is within the domain
  // ---------------------------------------------------
  int inDomain( TMRQuadrant *p );
  int onBoundary( TMRQuadrant *p );

  // Refine the quadtree
  // -----------------
  void refine( int refinement[],
	       int min_level=0, int max_level=TMR_MAX_LEVEL );

  // Balance the tree to ensure 2-1 balancing
  // ----------------------------------------
  void balance( int balance_corner=0 );

  // Coarsen the tree uniformly
  // --------------------------
  TMRQuadtree *coarsen();

  // Find an quadrant that completely encloses the provided quadrant
  // ---------------------------------------------------------------
  TMRQuadrant* findEnclosing( TMRQuadrant *oct );
  void findEnclosingRange( TMRQuadrant *oct,
			   int *low, int *high );

  // Create the connectivity information for the mesh
  // ------------------------------------------------
  void createMesh( int _order );

  // Order the nodes but do not create the connectivity
  // --------------------------------------------------
  void createNodes( int _order );

  // Retrieve the mesh information
  // -----------------------------
  void getMesh( int *_num_nodes, 
                int *_num_elements, 
                const int **_elem_ptr, 
                const int **_elem_conn );

  // Retrieve the depednent node information
  // ---------------------------------------
  void getDependentMesh( int *_num_dep_nodes, 
                         const int **_dep_ptr,
                         const int **_dep_conn,
                         const double **_dep_weights );

  // Create the interpolation from a coarser mesh
  // --------------------------------------------
  void createInterpolation( TMRQuadtree *coarse,
                            int **_interp_ptr,
                            int **_interp_conn,
                            double **_interp_weights );
  void createRestriction( TMRQuadtree *coarse,
                          int **_interp_ptr,
                          int **_interp_conn,
                          double **_interp_weights );

  // Print a representation of the tree to a file
  // --------------------------------------------
  void printQuadtree( const char * filename );

  // Retrieve the quadtree elements
  // ------------------------------
  void getElements( TMRQuadrantArray **_elements ){
    if (_elements){ *_elements = elements; }
  }

  // Retrieve the quadtree nodes
  // ---------------------------
  void getNodes( TMRQuadrantArray **_nodes ){
    if (_nodes){ *_nodes = nodes; }
  }

 private:
  // The list quadrants representing the elements
  TMRQuadrantArray *elements; 

  // The nodes within the element mesh
  TMRQuadrantArray *nodes; 

  // A list of quadrants that define the domain
  int ndomain;
  TMRQuadrant *domain;

  // Store the order of the mesh
  int order;

  // The number of elements and independent and dependent nodes
  int num_elements;
  int num_nodes, num_dependent_nodes;

  // The mesh connectivity
  int *elem_ptr, *elem_conn;
  
  // The dependent mesh connectivity
  int *dep_ptr, *dep_conn;
  double *dep_weights;
};

#endif // TMR_QUADTREE_H
