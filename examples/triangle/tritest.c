#include "TMRTriangleInterface.h"
#include "TMRBspline.h"
#include <stdio.h>
#include <math.h>

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  double R = 2.0;
  const int nu = 2, ku = 2;
  const int nv = 2, kv = 2;
  TMRPoint pts[nu*nv];

  for ( int j = 0; j < nu; j++ ){
    for ( int i = 0; i < nv; i++ ){
      double u = -10.0 + 20.0*i/(nu-1);
      double v = -10.0 + 20.0*j/(nv-1);
      pts[nu*j+i].x = 1.0*u;
      pts[nu*j+i].y = 1.0*v;
      pts[nu*j+i].z = 0.0;
    }
  }

  TMRBsplineSurface *surf = new TMRBsplineSurface(nu, nv, ku, kv, pts);
  surf->incref();

  int npts = 100;
  double *prms = new double[ 2*npts ];

  for ( int i = 0; i < npts; i++ ){
    TMRPoint P;
    double u = (2.0*M_PI*i)/npts;
    P.x = R*cos(u);
    P.y = R*sin(u);
    P.z = 0.0;
    surf->invEvalPoint(P, &prms[2*i], &prms[2*i+1]);
    surf->evalPoint(prms[2*i], prms[2*i+1], &P);
  }

  int nsegs = npts;
  int *seg = new int[ 2*nsegs ];
  for ( int i = 0; i < nsegs; i++ ){
    if (i == nsegs-1){
      seg[2*i] = i;
      seg[2*i+1] = 0;
    }
    else {
      seg[2*i] = i;
      seg[2*i+1] = i+1;
    }
  }
  
  double length = 2.0*M_PI*R/(npts-1);

  int nholes = 0;
  double *holes = NULL;

  // Triangulate the region
  TMRTriangulation *tri = new TMRTriangulation(npts, prms, NULL, surf);
  tri->setSegments(nsegs, seg);
  tri->create();
  tri->refine(length);

  for ( int k = 0; k < 10; k++ ){
    tri->springSmoothing(100);
    tri->remesh();
  }
  tri->springSmoothing(50);
  tri->printTriQuality();
  tri->writeToVTK("triangle.vtk");
  
  tri->recombine();
  tri->printQuadQuality();
  tri->writeQuadToVTK("match.vtk");
  tri->springQuadSmoothing(100);
  tri->printQuadQuality();
  tri->writeQuadToVTK("smoothed.vtk");

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
