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
  forest->balance(1);

  // Create the nodes
  int order = 2;
  forest->createNodes(order);

  /*
  // Get the nodes/elements
  int nnodes, ndep_nodes, nelems;
  int *elem_conn, *dep_conn;
  double *dep_weights;

  // Get the mesh
  forest->getMesh(&nnodes, &ndep_nodes, &nelems, 
                  &elem_conn, &dep_conn, &dep_weights);
  */
  /*
  // Create the position array
  double *X = new double[ 3*(nnodes + ndep_nodes) ];

  // Get the quadtrees - some may be NULL!
  TMRQuadtree **trees;
  int ntrees = forest->getQuadtrees(&trees);

  for ( int i = 0; i < ntrees; i++ ){
    if (trees[i]){
      TMRQuadrantArray *nodes;
      trees[i]->getNodes(&nodes);

      // Get the elements
      int size;
      TMRQuadrant *array;
      nodes->getArray(&array, &size);

      // Write out this portion of the forrest
      for ( int k = 0; k < size; k++ ){
        int32_t x = array[k].x;
        int32_t y = array[k].y;

        if (array[k].tag >= 0){
          getLocation(i, x, y, &X[3*array[k].tag]);
        }
      }
    }
  }

  // Set the dependent node locations
  for ( int i = 0, n = nnodes; i < ndep_nodes; i++, n++ ){
    X[3*n] = X[3*n+1] = X[3*n+2] = 0.0;

    for ( int j = 0; j < order; j++ ){
      int k = dep_conn[order*i + j];
      if (k >= 0){
        X[3*n] += dep_weights[order*i + j]*X[3*k];
        X[3*n+1] += dep_weights[order*i + j]*X[3*k+1];
        X[3*n+2] += dep_weights[order*i + j]*X[3*k+2];
      }
      else {
        X[3*n] = 1.0;
        X[3*n+1] = 2.0;
      }
    }
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  char filename[128];
  sprintf(filename, "parallel%d.dat", rank);
  FILE *fp = fopen(filename, "w");

  // Write the tecplot header
  fprintf(fp, "Variables = X, Y\n");
  fprintf(fp, "ZONE T=TMR%d N=%d E=%d ", rank, 
          nnodes + ndep_nodes, (order-1)*(order-1)*nelems);
  fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEQUADRILATERAL\n");

  // Write out the node locations
  for ( int i = 0; i < nnodes + ndep_nodes; i++ ){
    fprintf(fp, "%e %e\n", X[3*i], X[3*i+1]);
  }

  for ( int i = 0; i < nelems; i++ ){
    if (order == 2){
      int n[4];
      for ( int k = 0; k < 4; k++ ){
        n[k] = elem_conn[order*order*i + k];
        if (n[k] < 0){
          n[k] = nnodes - n[k]-1;
        }
      }
    
      fprintf(fp, "%d %d %d %d\n",
              n[0]+1, n[1]+1, n[3]+1, n[2]+1);
    }
  }

  fclose(fp);
  */

  // Print out the forrest using the quadrants
  int rank;
  MPI_Comm_rank(comm, &rank);
  char filename[128];
  sprintf(filename, "parallel%d.dat", rank);
  FILE *fp = fopen(filename, "w");
  
  fprintf(fp, "Variables = X, Y, var\n");
  
  // Get the quadtrees - some may be NULL!
  TMRQuadtree **trees;
  int ntrees = forest->getQuadtrees(&trees);
  
  for ( int i = 0; i < ntrees; i++ ){
    if (trees[i]){
      TMRQuadrantArray *elements;
      trees[i]->getElements(&elements);

      // Get the elements
      int size;
      TMRQuadrant *array;
      elements->getArray(&array, &size);

      TMRQuadrantArray *nodes;
      trees[i]->getNodes(&nodes);

      fprintf(fp, "ZONE T=TMR%d N=%d E=%d ", i, 4*size, size);
      fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEQUADRILATERAL\n");

      // Write out this portion of the forrest
      for ( int k = 0; k < size; k++ ){
        int32_t h = 1 << (TMR_MAX_LEVEL - array[k].level);
        int x = array[k].x;
        int y = array[k].y;
        TMRQuadrant q = array[k], *t;
        
        double X[3];
        t = nodes->contains(&q, 1);
        getLocation(i, x, y, X);
        fprintf(fp, "%e %e %d\n", X[0], X[1], t->tag);

        q.x = array[k].x + h;
        t = nodes->contains(&q, 1);
        getLocation(i, x+h, y, X);
        fprintf(fp, "%e %e %d\n", X[0], X[1], t->tag);

        q.x = array[k].x + h;
        q.y = array[k].y + h;
        getLocation(i, x+h, y+h, X);
        t = nodes->contains(&q, 1);
        fprintf(fp, "%e %e %d\n", X[0], X[1], t->tag);

        q.x = array[k].x;
        q.y = array[k].y + h;
        getLocation(i, x, y+h, X);
        t = nodes->contains(&q, 1);
        fprintf(fp, "%e %e %d\n", X[0], X[1], t->tag);
      }

      for ( int k = 0; k < size; k++ ){
	fprintf(fp, "%d %d %d %d\n", 
                4*k+1, 4*k+2, 4*k+3, 4*k+4);
      }
    }
  }

  delete forest;

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
