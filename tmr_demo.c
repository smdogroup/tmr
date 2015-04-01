
#include "TMROctree.h"

int main( int argc, char * argv[] ){
  
  TMROctree *tree = new TMROctree(1000, 2, 10);

  tree->balance(7);
  tree->printTree("octree.dat");

  TMROctree *coarse = tree->coarsen();
  coarse->balance(1);
  coarse->printTree("coarse_octree.dat");

  TMROctree *coarsest = coarse->coarsen();
  coarsest->balance(1);
  coarsest->printTree("coarsest_octree.dat");

  return (0);
}
