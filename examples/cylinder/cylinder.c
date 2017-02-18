#include "TMRGeometry.h"
#include "TMRBspline.h"
#include "TMRMesh.h"
#include <math.h>

/*
  Create a circle centered at the point c with radius r in the (x,y)
  plane.
*/
TMRBsplineCurve* createSemiCircle( TMRPoint c, double r ){
  // Set the points and weights for the B-spline circle
  const int nctl = 5;
  const int ku = 3;
  TMRPoint p[nctl];
  double wts[nctl];
  memset(p, 0, nctl*sizeof(TMRPoint));

  // Set the knot locations
  double Tu[] = {
    0.0, 0.0, 0.0,
    0.5, 0.5,
    1.0, 1.0, 1.0};
  
  for ( int k = 0; k < nctl; k++ ){
    p[k] = c;
  }

  // Set the weights
  double sqrt2 = 1.0/sqrt(2.0);

  // c + (r,0)
  p[0].x += r;
  wts[0] = 1.0;

  // c + (r,r)
  p[1].x += r;
  p[1].y += r;
  wts[1] = sqrt2;

  // c + (0,r)
  p[2].y += r;
  wts[2] = 1.0;
  
  // c + (-r,r)
  p[3].x -= r;
  p[3].y += r;
  wts[3] = sqrt2;

  // c + (-r,0)
  p[4].x -= r;
  wts[4] = 1.0;

  // Create the circle
  TMRBsplineCurve *curve =
    new TMRBsplineCurve(nctl, ku, Tu, wts, p);

  // Return the curve
  return curve;
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  double htarget = 0.1;
  for ( int i = 0; i < argc; i++ ){
    if (sscanf(argv[i], "h=%lf", &htarget) == 1){
      if (htarget < 0.01){ htarget = 0.01; }
      if (htarget > 0.5){ htarget = 0.5; }
    }
  }
  
  // Set the radius and height of the circle
  double r1 = 1.0;
  double r2 = 0.15;
  double h = 2.0;

  TMRPoint p;
  p.x = p.y = p.z = 0.0;

  TMRBsplineCurve *loft[2];
  loft[0] = createSemiCircle(p, r1);
  loft[0]->incref();

  p.z = h;
  loft[1] = createSemiCircle(p, r2);
  loft[1]->incref();

  // Loft the circles together to form a closed surface
  TMRCurveLofter *lofter = new TMRCurveLofter(loft, 2);
  lofter->incref();

  // Create the surface object from the lofted surface
  TMRBsplineSurface *surface = lofter->createSurface(2);
  surface->incref();
  lofter->decref();

  // Create the boundary curves for the surface
  TMRPoint pts[2];
  pts[0].x = pts[0].y = pts[0].z = 0.0;
  pts[0].x = r1;
  pts[1].x = pts[1].y = pts[1].z = 0.0;
  pts[1].x = r2;
  pts[1].z = h;
  TMRVertexFromPoint *v1 = new TMRVertexFromPoint(pts[0]);
  TMRVertexFromPoint *v2 = new TMRVertexFromPoint(pts[1]);
  TMRBsplineCurve *line1 = new TMRBsplineCurve(2, 2, pts);
  line1->incref();

  // Set the points for the second line connecting the semi-circular
  // segments.
  pts[0].x = pts[0].y = pts[0].z = 0.0;
  pts[0].x = -r1;
  pts[1].x = pts[1].y = pts[1].z = 0.0;
  pts[1].x = -r2;
  pts[1].z = h;
  TMRVertexFromPoint *v3 = new TMRVertexFromPoint(pts[0]);
  TMRVertexFromPoint *v4 = new TMRVertexFromPoint(pts[1]);
  TMRBsplineCurve *line2 = new TMRBsplineCurve(2, 2, pts);
  line2->incref();

  // Create the curves and add them to the surface
  int num_curves = 4;
  TMRCurve *curves[4];
  curves[0] = loft[0];
  curves[1] = line2;
  curves[2] = loft[1];
  curves[3] = line1;

  // Set the vertices
  curves[0]->setVertices(v1, v3);
  curves[1]->setVertices(v3, v4);
  curves[2]->setVertices(v2, v4);
  curves[3]->setVertices(v1, v2);
  
  // Set the directions of the curves
  int dir[4];
  dir[0] = 1;
  dir[1] = 1;
  dir[2] = -1;
  dir[3] = -1;

  surface->addCurveSegment(num_curves, curves, dir);

  // Write out the curves and cylindrical surface
  for ( int k = 0; k < num_curves; k++ ){
    char filename[128];
    sprintf(filename, "curve%d.vtk", k);
    curves[k]->writeToVTK(filename);
  }
  surface->writeToVTK("cylinder.vtk");

  // Create the mesh
  TMRSurfaceMesh *mesh = new TMRSurfaceMesh(surface);
  mesh->incref();
  mesh->mesh(htarget);
  mesh->writeToVTK("quads.vtk");
  mesh->decref();  


  TMRFinalize();
  MPI_Finalize();
  return (0);
}
