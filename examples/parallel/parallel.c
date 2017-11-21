#include "TMRQuadForest.h"

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

void getLocation( int face, int x, int y, double *X ){
  // Compute the u/v locations on the face
  const double dh = 1.0/(1 << TMR_MAX_LEVEL);
  double u = dh*x, v = dh*y;

  // Evaluate the shape functions
  double N[4];
  N[0] = (1.0 - u)*(1.0 - v);
  N[1] = u*(1.0 - v);
  N[2] = (1.0 - u)*v;
  N[3] = u*v;
  
  // Compute the node location
  X[0] = (N[0]*Xpts[3*conn[4*face]] +
          N[1]*Xpts[3*conn[4*face+1]] +
          N[2]*Xpts[3*conn[4*face+2]] +
          N[3]*Xpts[3*conn[4*face+3]]);
  X[1] = (N[0]*Xpts[3*conn[4*face]+1] +
          N[1]*Xpts[3*conn[4*face+1]+1] +
          N[2]*Xpts[3*conn[4*face+2]+1] +
          N[3]*Xpts[3*conn[4*face+3]+1]);
  X[2] = (N[0]*Xpts[3*conn[4*face]+2] +
          N[1]*Xpts[3*conn[4*face+1]+2] +
          N[2]*Xpts[3*conn[4*face+2]+2] +
          N[3]*Xpts[3*conn[4*face+3]+2]);
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  MPI_Comm comm = MPI_COMM_WORLD;
  TMRQuadForest *forest = new TMRQuadForest(comm);
    
  // Set the connectivity
  int num_nodes = 8;
  int num_faces = 5;
  forest->setConnectivity(num_nodes, conn, num_faces);

  // Create the random trees
  forest->createRandomTrees(175, 0, 10);
  forest->repartition();
  forest->balance(1);

  // Create the nodes
  int order = 2;
  forest->createNodes(order);

  // Print out the forrest using the quadrants
  int rank;
  MPI_Comm_rank(comm, &rank);
  char filename[128];
  sprintf(filename, "parallel%d.dat", rank);
  FILE *fp = fopen(filename, "w");
 
  fprintf(fp, "Variables = X, Y\n");
  
  TMRQuadrantArray *quadrants;
  forest->getQuadrants(&quadrants);

  int size;
  TMRQuadrant *array;
  quadrants->getArray(&array, &size);

  fprintf(fp, "ZONE T=Quadrants N=%d E=%d ", 4*size, size);
  fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEQUADRILATERAL\n");

  // Write out this portion of the forrest
  for ( int k = 0; k < size; k++ ){
    int32_t h = 1 << (TMR_MAX_LEVEL - array[k].level);
    double X[3];
    
    for ( int jj = 0; jj < 2; jj++ ){
      for ( int ii = 0; ii < 2; ii++ ){
        int x = array[k].x + h*ii;
        int y = array[k].y + h*jj;
        getLocation(array[k].face, x, y, X);
        fprintf(fp, "%e %e\n", X[0], X[1]);
      }
    }
  }
  
  for ( int k = 0; k < size; k++ ){
    fprintf(fp, "%d %d %d %d\n", 4*k+1, 4*k+2, 4*k+4, 4*k+3);
  }

  delete forest;

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
