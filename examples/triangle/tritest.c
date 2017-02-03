#include "TMRTriangleInterface.h"
#include "TMRBspline.h"
#include <stdio.h>
#include <math.h>

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  const int nu = 30, ku = 5;
  const int nv = 30, kv = 5;
  TMRPoint pts[nu*nv];

  for ( int j = 0; j < nu; j++ ){
    for ( int i = 0; i < nv; i++ ){
      double u = -10.0 + 20.0*i/(nu-1);
      double v = -10.0 + 20.0*j/(nv-1);
      pts[nu*j+i].x = 1.0*u;
      pts[nu*j+i].y = 1.0*v;
      pts[nu*j+i].z = 0.1/(0.1 + u*u + v*v);
    }
  }

  TMRBsplineSurface *surf = new TMRBsplineSurface(nu, nv, ku, kv, pts);
  surf->incref();

  int npts = 200;
  double *prms = new double[ 2*npts ];

  for ( int i = 0; i < npts; i++ ){
    TMRPoint P;
    double u = 2.0*M_PI*1.0*i/npts;
    P.x = cos(u);
    P.y = sin(u);
    P.z = 0.1/(0.1 + P.x*P.x + P.y*P.y);
    surf->invEvalPoint(P, &prms[2*i], &prms[2*i+1]);
    surf->evalPoint(prms[2*i], prms[2*i+1], &P);
    printf("P(%f, %f) = %f %f %f\n", prms[2*i], prms[2*i+1], P.x, P.y, P.z);
  }

  int nsegs = npts;
  int *seg = new int[ 2*nsegs ];
  for ( int i = 0; i < nsegs; i++ ){
    seg[2*i] = i;
    seg[2*i+1] = i+1;
  }
  seg[2*nsegs-2] = 0;

  double length = 2.0*M_PI/(npts-1);
  double area = 0.5*length*length;

  int nholes = 0;
  double *holes = NULL;

  // Triangulate the region
  TMRTriangulation *tri = new TMRTriangulation(npts, prms, NULL, surf);
  tri->setSegments(nsegs, seg);
  tri->create();
  tri->refine(length); 

  for ( int k = 0; k < 5; k++ ){
    // tri->laplacianSmoothing(50);
    tri->springSmoothing(50);
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
