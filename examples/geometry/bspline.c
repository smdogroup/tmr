#include "TMRBspline.h"
#include <math.h>

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Set the control point locations
  TMRPoint pts[3];
  pts[0].zero();  pts[1].zero();  pts[2].zero();
  pts[0].x = 1.0;
  pts[1].x = 1.0;
  pts[1].y = 1.0;
  pts[2].y = 1.0;

  double wts[3];
  wts[0] = 1.0;
  wts[1] = 1.0/sqrt(2.0);
  wts[2] = 1.0;

  int nu = 3;
  int ku = 3;

  double Tu[6] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  double Tv[6] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};

  TMRBsplineCurve *bspline = new TMRBsplineCurve(nu, ku, pts);
  TMRBsplineCurve *nurbs = new TMRBsplineCurve(nu, ku, Tu, wts, pts);
  bspline->incref();
  nurbs->incref();

  TMRPoint pt(1.0, 2.0, 2.0);

  TMRPoint point;
  bspline->evalPoint(0.5, &point);
  printf("bspline = %f %f %f\n", point.x, point.y, point.z);

  nurbs->evalPoint(0.5, &point);
  printf("nurbs = %f %f %f\n", point.x, point.y, point.z);

  // Perform the inverse fit
  double t;
  bspline->invEvalPoint(pt, &t);
  printf("inv bspline = %20.15e\n", t);

  nurbs->invEvalPoint(pt, &t);
  printf("inv nubrs   = %20.15e\n", t);

  // Set the points and weights for the cylindrical surface
  TMRPoint spts[12];

  const double sq2 = 1.0/sqrt(2.0);
  const double x[] = 
    {0.0, 0.0, 1.0,
     1.0, 0.0, 1.0,
     1.0, 0.0, 0.0,
     0.0, 0.0, 1.0,
     sq2, sq2, 1.0,
     sq2, sq2, 0.0,
     0.0, 0.0, 1.0,
     0.0, 1.0, 1.0,
     0.0, 1.0, 0.0,
     0.0, 0.0, 1.0,
     -sq2, sq2, 1.0,
     -sq2, sq2, 0.0};

  const double swts[] = {
    1.0, sq2, 1.0, 1.0, sq2, 1.0, 
    1.0, sq2, 1.0, 1.0, sq2, 1.0};

  for ( int k = 0; k < 12; k++ ){
    // Set the point locations
    spts[k].x = x[3*k];
    spts[k].y = x[3*k+1];
    spts[k].z = x[3*k+2];
  }
  
  double sTv[] = {0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};

  nu = 3;
  int nv = 4;
  TMRBsplineSurface *sbspline = new TMRBsplineSurface(nu, nv, 
                                                      ku, ku, spts);
  TMRBsplineSurface *snurbs = new TMRBsplineSurface(nu, nv, ku, ku, 
                                                    Tu, sTv, swts, spts);

  sbspline->evalPoint(0.5, 0.5, &point);
  printf("bspline = %f %f %f\n", point.x, point.y, point.z);

  snurbs->evalPoint(0.5, 0.5, &point);
  printf("nurbs = %f %f %f\n", point.x, point.y, point.z);

  // Perform the inverse fit
  double u, v;
  sbspline->invEvalPoint(pt, &u, &v);
  printf("inv bspline = %20.15e %20.15e\n", u, v);

  snurbs->invEvalPoint(pt, &u, &v);
  printf("inv nubrs   = %20.15e %20.15e\n", u,v);

  int n = 50;
  FILE *fp = fopen("surface.dat", "w");
  fprintf(fp, "Variables = X, Y, Z\n");
  fprintf(fp, "Zone I=%d J=%d datapacking=point\n", n, n);

  // Get the parameter range
  double umin, vmin, umax, vmax;
  snurbs->getRange(&umin, &vmin, &umax, &vmax);

  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < n; i++ ){
      TMRPoint p;
      double u = umin + (umax-umin)*i/(n-1.0);
      double v = vmin + (vmax-vmin)*j/(n-1.0);
      snurbs->evalPoint(u, v, &p);
      fprintf(fp, "%e %e %e\n", p.x, p.y, p.z);
    }
  }

  fclose(fp);

  bspline->decref();
  nurbs->decref();

  // Create the control points
  int nctl = 150;
  TMRPoint line_pts[150];
  memset(line_pts, 0, nctl*sizeof(TMRPoint));
  
  for ( int i = 0; i < nctl; i++ ){
    line_pts[i].x = -5.0 + 10.0*i/(nctl-1);
    line_pts[i].y = 0.1*line_pts[i].x*line_pts[i].x + cos(M_PI*line_pts[i].x);
  }

  // Create the interpolation object
  TMRCurveInterpolation *interper = 
    new TMRCurveInterpolation(line_pts, nctl);
  interper->incref();

  interper->setNumControlPoints(13);
  
  // Create an interpolation
  fp = fopen("line.dat", "w");
  fprintf(fp, "Variables = X, Y, Z\n");

  TMRBsplineCurve *curve = interper->createCurve(4);

  const double Tnew[] = { 0.1, 0.15, 0.2, 0.25, 0.3 };
  int nnew = 5;
  curve = curve->refineKnots(Tnew, nnew);

  // Get the parameter range
  double tmin, tmax;
  curve->getRange(&tmin, &tmax);

  int nvals = 5*nctl;
  fprintf(fp, "Zone T = curve I=%d\n", nvals);
  for ( int i = 0; i < nvals; i++ ){
    TMRPoint p;
    double t = tmin + (tmax-tmin)*i/(nvals-1.0);
    curve->evalPoint(t, &p);
    fprintf(fp, "%e %e %e\n", p.x, p.y, p.z);
  }

  int npts;
  const TMRPoint *ctrl_pts;
  curve->getData(&npts, NULL, NULL, NULL, &ctrl_pts);

  fprintf(fp, "Zone T = pts I=%d\n", npts);
  for ( int i = 0; i < npts; i++ ){
    fprintf(fp, "%e %e %e\n", 
            ctrl_pts[i].x, ctrl_pts[i].y, ctrl_pts[i].z); 
  }

  fprintf(fp, "Zone T = interpolation I=%d\n", nctl);
  for ( int i = 0; i < nctl; i++ ){
    fprintf(fp, "%e %e %e\n", 
            line_pts[i].x, line_pts[i].y, line_pts[i].z); 
  }

  fclose(fp);

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
