#include "TMRForrest.h"

// Create the nodal locations for the mesh
const double Xpts[] = {0.0, 0.0,
                       1.0, 1.0,
                       1.0, 0.0, 
                       0.5, 0.0,
                       0.0, 0.5,
                       0.35, 0.35,
                       0.0, 1.0};
const int test_conn[] = {0, 3, 4, 5,
                         5, 3, 1, 2,
                         4, 5, 6, 1};

/*
  Retrieve the x/y location on the face based on the u/v coordinates
*/
void get_location( int face, double u, double v, 
                   double *x, double *y ){
  double N[4];
  N[0] = (1.0 - u)*(1.0 - v);
  N[1] = u*(1.0 - v);
  N[2] = (1.0 - u)*v;
  N[3] = u*v;

  *x = (N[0]*Xpts[2*test_conn[4*face]] +
        N[1]*Xpts[2*test_conn[4*face+1]] +
        N[2]*Xpts[2*test_conn[4*face+2]] +
        N[3]*Xpts[2*test_conn[4*face+3]]);
  *y = (N[0]*Xpts[2*test_conn[4*face]+1] +
        N[1]*Xpts[2*test_conn[4*face+1]+1] +
        N[2]*Xpts[2*test_conn[4*face+2]+1] +
        N[3]*Xpts[2*test_conn[4*face+3]+1]);
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  srand(time(NULL));

  TMRQuadForrest *forrest = new TMRQuadForrest(MPI_COMM_WORLD);

  // Set the connectivity
  int num_nodes = 7;
  int num_faces = 3;
  forrest->setConnectivity(num_nodes, 
                           test_conn, num_faces);
  
  // Allocate the trees (random trees for now)
  int refine_level = 3;
  forrest->createTrees(refine_level);

  // Balance the forrest
  forrest->balance(1);

  TMRQuadtree **quad;
  forrest->getQuadtrees(&quad);

  int num_elements = 0;
  for ( int face = 0; face < num_faces; face++ ){
    TMRQuadrantArray *elements;
    quad[face]->getElements(&elements);

    int size;
    TMRQuadrant *array;
    elements->getArray(&array, &size);
    num_elements += size;
  }

  FILE * fp = fopen("forrest_test.dat", "w");
  fprintf(fp, "Variables = X, Y\n");
  fprintf(fp, "ZONE T=TMR N=%d E=%d ", 4*num_elements, num_elements);
  fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEQUADRILATERAL\n");

  double dh = 1.0/(1 << TMR_MAX_LEVEL);

  for ( int face = 0; face < num_faces; face++ ){
    TMRQuadrantArray *elements;
    quad[face]->getElements(&elements);

    int size;
    TMRQuadrant *array;
    elements->getArray(&array, &size);

    for ( int i = 0; i < size; i++ ){
      int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      
      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          double u = dh*(array[i].x + ii*h);
          double v = dh*(array[i].y + jj*h);
          
          double x, y;
          get_location(face, u, v, &x, &y);
          fprintf(fp, "%e %e\n", x, y);
        }
      }
    }
  }
  
  for ( int i = 0; i < num_elements; i++ ){
    fprintf(fp, "%d %d %d %d\n", 4*i+1, 4*i+2, 4*i+4, 4*i+3);
  }
  fclose(fp);

  MPI_Finalize();
  return (0);
}
