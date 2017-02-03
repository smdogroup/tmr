#include "TMRTriangleInterface.h"
#include "TMRBspline.h"
#include <stdio.h>
#include <math.h>

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  int nu = 2, ku = 2;
  int nv = 2, kv = 2;
  TMRPoint pts[4];

  pts[0].x = -10.0;
  pts[0].y = -10.0;
  pts[0].z = 0.0;

  pts[1].x = 10.0;
  pts[1].y = -10.0;
  pts[1].z = 0.0;

  pts[2].x = -10.0;
  pts[2].y = 10.0;
  pts[2].z = 0.0;

  pts[3].x = 10.0;
  pts[3].y = 10.0;
  pts[3].z = 0.0;

  TMRBsplineSurface *surf = new TMRBsplineSurface(nu, nv, ku, kv, pts);
  surf->incref();

  int npts = 350;
  double *prms = new double[ 2*npts ];

  for ( int i = 0; i < npts; i++ ){
    TMRPoint P;
    double u = 2.0*M_PI*1.0*i/npts;
    P.x = cos(u);
    P.y = sin(u);
    P.z = 0.0;
    surf->invEvalPoint(P, &prms[2*i], &prms[2*i+1]);
  }

  int nsegs = npts;
  int *seg = new int[ 2*nsegs ];
  for ( int i = 0; i < nsegs; i++ ){
    seg[2*i] = i;
    seg[2*i+1] = i+1;
  }
  seg[2*nsegs-2] = 0;

  double length = 2.0*M_PI/npts;
  double area = 0.5*length*length;

  int nholes = 0;
  double *holes = NULL;

  // Triangulate the region
  TMRTriangulation *tri = new TMRTriangulation(npts, prms, NULL, surf);
  tri->setSegments(nsegs, seg);
  tri->create();

  // Get the number of triangles
  int ntris;
  tri->getTriangulation(&ntris, NULL);

  double *areas = new double[ ntris ];
  for ( int i = 0; i < ntris; i++ ){
    areas[i] = 1.1*area;
  }
  
  // Refine the areas
  tri->refine(areas);
  delete [] areas;

  // Get the number of triangles
  tri->getTriangulation(&ntris, NULL);
  
  area = area/(20*20);

  areas = new double[ ntris ];
  for ( int i = 0; i < ntris; i++ ){
    areas[i] = area;
  }
  tri->refine(areas);
  delete [] areas;

 for ( int k = 0; k < 5; k++ ){
    tri->laplacianSmoothing(100);
    tri->remesh();
  }

  tri->writeToVTK("triangle.vtk");
  tri->recombine();
  tri->printQuadQuality();

  tri->writeQuadToVTK("match.vtk");
  tri->laplacianQuadSmoothing(100);
  tri->printQuadQuality();
  tri->writeQuadToVTK("smoothed.vtk");
  
  TMRFinalize();
  MPI_Finalize();
  return (0);
}
