#include "TMRTriangleInterface.h"
#include <stdio.h>

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Define the points
  int npts = 16;
  const double pts[] =
    {0.0, 0.0,
     0.25, 0.0,
     0.5, 0.0,
     0.75, 0.0,
     1.0, 0.0,
     1.0, 0.25,
     1.0, 0.5,
     1.0, 0.75,
     1.0, 1.0,
     0.75, 1.0,
     0.5, 1.0,
     0.25, 1.0,
     0.0, 1.0,
     0.0, 0.75,
     0.0, 0.5,
     0.0, 0.25};

  int nsegs = 16;
  const int seg[] =
    {0, 1,
     1, 2,
     2, 3,
     3, 4,
     4, 5,
     5, 6,
     6, 7,
     7, 8,
     8, 9,
     9, 10,
     10, 11,
     11, 12,
     12, 13,
     13, 14,
     14, 15,
     15, 0};

  const int segmarks[] =
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  
  double area = 0.5*0.25*0.25;

  int nholes = 0;
  double *holes = NULL;

  // Triangulate the region
  TMRTriangulation *tri = new TMRTriangulation(npts, pts);
  tri->setSegments(nsegs, seg, segmarks);
  tri->create();

  // Get the number of triangles
  int ntris;
  tri->getTriangulation(&ntris, NULL);

  double *areas = new double[ ntris ];
  for ( int i = 0; i < ntris; i++ ){
    areas[i] = area;
  }
  
  // Refine the areas
  printf("Refine\n");
  tri->refine(areas);
  tri->writeToVTK("triangle.vtk");
  
  TMRFinalize();
  MPI_Finalize();
  return (0);
}
