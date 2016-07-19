#include "TMROctForest.h"

void getLocation( int i, int nx, int ny, int nz,
                  int x, int y, int z, double X[] ){
  int iz = i/(nx*ny);
  int iy = (i - iz*nx*ny)/nx;
  int ix = i % nx;

  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  double u = (1.0*x)/hmax;
  double v = (1.0*y)/hmax;
  double w = (1.0*z)/hmax;

  X[0] = ix + u;
  X[1] = iy + v;
  X[2] = iz + w;
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  MPI_Comm comm = MPI_COMM_WORLD;
  TMROctForest *forest = new TMROctForest(comm);
  
  // Create a regular connectivity
  int nx = 3, ny = 2, nz = 4;
  int num_nodes = (nx+1)*(ny+1)*(nz+1);
  int num_blocks = nx*ny*nz;

  int *conn = new int[ 8*num_blocks ];

  for ( int iz = 0; iz < nz; iz++ ){
    for ( int iy = 0; iy < ny; iy++ ){
      for ( int ix = 0; ix < nx; ix++ ){
        int block = ix + nx*iy + nx*ny*iz;

        for ( int kz = 0; kz < 2; kz++ ){
          for ( int ky = 0; ky < 2; ky++ ){
            for ( int kx = 0; kx < 2; kx++ ){
              conn[8*block + kx + 2*ky + 4*kz] = 
                (ix + kx) + (nx+1)*(iy + ky) + (nx+1)*(ny+1)*(iz + kz);
            }
          }
        }
      }
    }
  }
  
  forest->setConnectivity(num_nodes, conn, num_blocks);
  delete [] conn;
  
  // Create the random trees
  forest->createRandomTrees(175, 0, 10);
  forest->balance(1);

  forest->createNodes(2);

  // Get the octrees within the forest
  TMROctree **trees;
  int ntrees = forest->getOctrees(&trees);

  // Get the rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Write out a file for each processor - bad practice!
  char filename[128];
  sprintf(filename, "parallel%d.dat", rank);
  FILE *fp = fopen(filename, "w");

  // Write the tecplot header
  fprintf(fp, "Variables = X, Y, Z\n");

  for ( int i = 0; i < ntrees; i++ ){
    if (trees[i]){
      TMROctantArray *elements;
      trees[i]->getElements(&elements);

      // Get the elements
      int size;
      TMROctant *array;
      elements->getArray(&array, &size);

      fprintf(fp, "ZONE T=TMR%d N=%d E=%d ", i, 8*size, size);
      fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEBRICK\n");

      // Write out this portion of the forrest
      for ( int k = 0; k < size; k++ ){
        int h = 1 << (TMR_MAX_LEVEL - array[k].level);
        int x = array[k].x;
        int y = array[k].y;
        int z = array[k].z;

        double X[3];
        for ( int kz = 0; kz < 2; kz++ ){
          for ( int ky = 0; ky < 2; ky++ ){
            for ( int kx = 0; kx < 2; kx++ ){
              getLocation(i, nx, ny, nz,
                          x + kx*h, y + ky*h, z + kz*h, X);
              fprintf(fp, "%e %e %e\n", X[0], X[1], X[2]);
            }
          }
        }
      }

      for ( int k = 0; k < size; k++ ){
	fprintf(fp, "%d %d %d %d %d %d %d %d\n",
                8*k+1, 8*k+2, 8*k+4, 8*k+3,
                8*k+5, 8*k+6, 8*k+8, 8*k+7);
      }
    }
  }

  fclose(fp);

  delete forest;

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
