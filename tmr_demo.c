#include "TMROctree.h"
#include "mpi.h"

int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);

  TMROctree *next = NULL;
  for ( int k = 0; k < 4; k++ ){
    TMROctree *tree = NULL;
    if (k == 0){
      tree = new TMROctree(5000, 3, 10);
      double t1 = MPI_Wtime();
      tree->balance(7);
      t1 = MPI_Wtime() - t1;
      next = tree;

      printf("Balance time: %15.8f\n", t1);
    }
    else {
      double t1 = MPI_Wtime();
      tree = next->coarsen();
      double t2 = MPI_Wtime();
      tree->balance(1);
      double t3 = MPI_Wtime();      
      next = tree;

      printf("Coarsen time: %15.8f\n", t2-t1);
      printf("Balance time: %15.8f\n", t3-t2);
    }
    
    double t1 = MPI_Wtime();
    tree->createNodes(2);
    t1 = MPI_Wtime() - t1;
    printf("Nodes time:   %15.8f\n", t1);
  }

  MPI_Finalize();
  return (0);
}
