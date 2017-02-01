#include "TMRTriangleInterface.h"
#include <stdio.h>

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  int pts_per_side = 101;
  int npts = 4*(pts_per_side-1);
  double *pts = new double[ 2*(npts+1) ];

  for ( int i = 0; i < pts_per_side; i++ ){
    int node = i;
    pts[2*node] = 1.0*i/(pts_per_side-1);
    pts[2*node+1] = 0.0;

    node = i + pts_per_side-1;
    pts[2*node] = 1.0;
    pts[2*node+1] = 1.0*i/(pts_per_side-1);;

    node = i + 2*(pts_per_side-1);
    pts[2*node] = 1.0 - 1.0*i/(pts_per_side-1);
    pts[2*node+1] = 1.0;

    node = i + 3*(pts_per_side-1);
    pts[2*node] = 0.0;
    pts[2*node+1] = 1.0 - 1.0*i/(pts_per_side-1);
  }

  int nsegs = npts;
  int *seg = new int[ 2*nsegs ];
  int *segmarks = new int[ nsegs ];
  for ( int i = 0; i < nsegs; i++ ){
    seg[2*i] = i;
    seg[2*i+1] = i+1;
    segmarks[i] = 1;
  }
  seg[2*nsegs-2] = 0;

  double length = 1.0/(pts_per_side-1);
  double area = 0.5*length*length;

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
  delete [] areas;

  // Get the number of triangles
  tri->getTriangulation(&ntris, NULL);
  
  areas = new double[ ntris ];
  for ( int i = 0; i < ntris; i++ ){
    areas[i] = area;
  }
  tri->refine(areas);
  delete [] areas;
  tri->writeToVTK("triangle.vtk");

  tri->smooth();

  int new_npts;
  const double *new_pts;
  tri->getPoints(&new_npts, &new_pts);

  TMRPoint *p = new TMRPoint[ new_npts ];
  for ( int i = 0; i < new_npts; i++ ){
    p[i].x = new_pts[2*i];
    p[i].y = new_pts[2*i+1];
    p[i].z = 0.0;
  }
  tri->recombine(p);
  
  TMRFinalize();
  MPI_Finalize();
  return (0);
}
