#include "TMRForest.h"

const double Xpts[] = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.3, 0.7, 0.0,
                       0.8, 0.25, 0.0,
                       0.25, 0.2, 0.0,
                       0.75, 0.6, 0.0,
                       0.0, 1.0, 0.0,
                       1.0, 1.0, 0.0};

const int conn[] = {0, 1, 4, 3,
                    2, 4, 5, 3, 
                    6, 0, 2, 4, 
                    2, 5, 6, 7,
                    3, 1, 5, 7};


int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  MPI_Comm comm = MPI_COMM_WORLD;
  TMRQuadForest *forest = new TMRQuadForest(comm);
    
  // Set the connectivity
  int num_nodes = 8;
  int num_faces = 5;
  forest->setConnectivity(num_nodes, conn, num_faces);

  forest->createRandomTrees(5000, 0, TMR_MAX_LEVEL);
  forest->balance(1);

  delete forest;

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
