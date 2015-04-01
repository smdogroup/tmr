#ifndef TMR_OCTANT_TREE_H
#define TMR_OCTANT_TREE_H

#include <stdio.h>
#include "TMROctant.h"

class TMROctree {
 public:
  static const int MAX_OCTANT_LEVELS = 30;

  TMROctree( int refine_level );
  TMROctree( int nrand, int min_level, int max_level );
  TMROctree( TMROctantArray *_list, int _is_sorted );
  ~TMROctree();

  void printTree( const char * filename );
  void balance( int balance_type = 3 );  
  TMROctree *coarsen();

 private:
  int size, max_size;
  int is_sorted;

  TMROctantArray *list;
};

#endif // TMR_OCTANT_TREE_H
