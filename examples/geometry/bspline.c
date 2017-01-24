#include "TMRBspline.h"
#include <math.h>

const int rae2822_npts = 128;
const double rae2822_pts[] =
  {1.0000, 0.0000,
   0.9994, 0.0003,
   0.99759, 0.00069,
   0.99459, 0.00132,
   0.99039, 0.00218,
   0.98502, 0.00326,
   0.97847, 0.00455,
   0.97077, 0.00606,
   0.96194, 0.00775,
   0.952, 0.00964,
   0.94096, 0.0117,
   0.92886, 0.01393,
   0.91574, 0.01627,
   0.9016, 0.01874,
   0.88651, 0.02131,
   0.87048, 0.02397,
   0.85355, 0.0267,
   0.83578, 0.02948,
   0.8172, 0.03231,
   0.79785, 0.03514,
   0.77778, 0.03795,
   0.75705, 0.04075,
   0.7357, 0.04338,
   0.71378, 0.04612,
   0.69134, 0.04857,
   0.66845, 0.05112,
   0.64514, 0.05339,
   0.62149, 0.05547,
   0.59754, 0.05733,
   0.57336, 0.05895,
   0.54901, 0.0603,
   0.52453, 0.06135,
   0.5, 0.06212,
   0.47547, 0.06261,
   0.45099, 0.06286,
   0.42663, 0.06285,
   0.40245, 0.06263,
   0.37851, 0.0622,
   0.35486, 0.06155,
   0.33156, 0.0607,
   0.30866, 0.05967,
   0.28622, 0.05848,
   0.2643, 0.05713,
   0.24295, 0.05556,
   0.22221, 0.05377,
   0.20215, 0.05187,
   0.1828, 0.04987,
   0.16422, 0.04778,
   0.14645, 0.04558,
   0.12952, 0.04321,
   0.11349, 0.04073,
   0.0984, 0.03817,
   0.08427, 0.03552,
   0.07114, 0.0328,
   0.05904, 0.03004,
   0.04801, 0.02726,
   0.03806, 0.02445,
   0.02923, 0.02163,
   0.02153, 0.01875,
   0.01498, 0.01579,
   0.00961, 0.01269,
   0.00541, 0.00945,
   0.00241, 0.00642,
   0.0006, 0.00323,
   0.0000, 0.0000,
   0.0006, -0.00317,
   0.00241, -0.00658,
   0.00541, -0.00957,
   0.00961, -0.01273,
   0.01498, -0.0158,
   0.02153, -0.0188,
   0.02923, -0.0218,
   0.03806, -0.02472,
   0.04801, -0.02761,
   0.05904, -0.03042,
   0.07114, -0.03315,
   0.08427, -0.03584,
   0.0984, -0.03844,
   0.11349, -0.04094,
   0.12952, -0.04333,
   0.14645, -0.04561,
   0.16422, -0.04775,
   0.1828, -0.04977,
   0.20215, -0.05167,
   0.22221, -0.0534,
   0.24295, -0.05498,
   0.2643, -0.05638,
   0.28622, -0.05753,
   0.30866, -0.05843,
   0.33156, -0.059,
   0.35486, -0.05919,
   0.37851, -0.05893,
   0.40245, -0.05817,
   0.42663, -0.05689,
   0.45099, -0.05515,
   0.47547, -0.05297,
   0.5, -0.05044,
   0.52453, -0.04761,
   0.54901, -0.04452,
   0.57336, -0.04127,
   0.59754, -0.03791,
   0.62149, -0.03463,
   0.64514, -0.0311,
   0.66845, -0.0277,
   0.69134, -0.02438,
   0.71378, -0.02118,
   0.7357, -0.01812,
   0.75705, -0.01524,
   0.77778, -0.01256,
   0.79785, -0.01013,
   0.8172, -0.00792,
   0.83578, -0.00594,
   0.85355, -0.00422,
   0.87048, -0.00273,
   0.88651, -0.00149,
   0.9016, -0.00049,
   0.91574, 0.00027,
   0.92886, 0.00081,
   0.94096, 0.00113,
   0.952, 0.00125,
   0.96194, 0.00125,
   0.97077, 0.00113,
   0.97847, 0.00094,
   0.98502, 0.00071,
   0.99039, 0.00048,
   0.99459, 0.00026,
   0.99759, 0.00009,
   0.9994, -0.00001,
   1.000, 0.0000 };

void test_surface_lofter(){
  // Create the control points
  int nctl = rae2822_npts;
  TMRPoint line_pts[rae2822_npts];
  memset(line_pts, 0, nctl*sizeof(TMRPoint));
  
  // Create a series of curves
  int num_curves = 5;
  TMRBsplineCurve *curves[5];
  double chord[5] = {6.0, 4.0, 3.0,  2.0, 1.0};
  double twist[5] = {0, -1.0, -2.0, -4.0, -10.0};
  for ( int k = 0; k < num_curves; k++ ){
    twist[k] *= M_PI/180.0;
  }

  for ( int k = 0; k < num_curves; k++ ){
    for ( int i = 0; i < nctl; i++ ){
      double c = cos(twist[k]);
      double s = sin(twist[k]);
      double x = rae2822_pts[2*i];
      double y = rae2822_pts[2*i+1];
      line_pts[i].x = chord[k]*(c*x + y*s) + 3*k;
      line_pts[i].y = chord[k]*(-s*x + y*c);
      line_pts[i].z = 5.0*k;
    }

    // Create the interpolation object
    TMRCurveInterpolation *interper = 
      new TMRCurveInterpolation(line_pts, nctl);
    interper->incref();

    interper->setNumControlPoints(13 + 2*k);
    curves[k] = interper->createCurve(4);
    curves[k]->incref();

    // Free the interpolation object
    interper->decref();
  }  

  // Create the lofter object
  TMRCurveLofter *lofter = new TMRCurveLofter(curves, num_curves);
  lofter->incref();

  // Create the surface object from the lofted surface
  TMRBsplineSurface *surf = lofter->createSurface(4);
  surf->incref();
  lofter->decref();

  // Print out the lofted surface
  int dim = 1;
  int n = 150;
  FILE *fp = fopen("surface.xyz", "wb");
  fwrite(&dim, sizeof(int), dim, fp);
  fwrite(&n, sizeof(int), dim, fp);
  fwrite(&n, sizeof(int), dim, fp);
  fwrite(&dim, sizeof(int), dim, fp);

  // Get the parameter range
  double umin, vmin, umax, vmax;
  surf->getRange(&umin, &vmin, &umax, &vmax);
  
  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < n; i++ ){
      TMRPoint p;
      double u = umin + (umax-umin)*i/(n-1.0);
      double v = vmin + (vmax-vmin)*j/(n-1.0);
      surf->evalPoint(u, v, &p);
      fwrite(&p.x, sizeof(double), dim, fp);
    }
  }

  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < n; i++ ){
      TMRPoint p;
      double u = umin + (umax-umin)*i/(n-1.0);
      double v = vmin + (vmax-vmin)*j/(n-1.0);
      surf->evalPoint(u, v, &p);
      fwrite(&p.y, sizeof(double), dim, fp);
    }
  }

  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < n; i++ ){
      TMRPoint p;
      double u = umin + (umax-umin)*i/(n-1.0);
      double v = vmin + (vmax-vmin)*j/(n-1.0);
      surf->evalPoint(u, v, &p);
      fwrite(&p.z, sizeof(double), dim, fp);
    }
  }

  fclose(fp);

  // Create an interpolation
  fp = fopen("line.dat", "w");
  fprintf(fp, "Variables = X, Y, Z\n");

  // Get the parameter range
  double tmin, tmax;
  curves[0]->getRange(&tmin, &tmax);

  int nvals = 5*nctl;
  fprintf(fp, "Zone T = curve I=%d\n", nvals);
  for ( int i = 0; i < nvals; i++ ){
    TMRPoint p;
    double t = tmin + (tmax-tmin)*i/(nvals-1.0);
    // curves[0]->evalPoint(t, &p);
    surf->evalPoint(t, 0.5, &p);
    fprintf(fp, "%e %e %e\n", p.x, p.y, p.z);
  }

  /* int npts; */
  /* const TMRPoint *ctrl_pts; */
  /* curves[0]->getData(&npts, NULL, NULL, NULL, &ctrl_pts); */

  /* fprintf(fp, "Zone T = pts I=%d\n", npts); */
  /* for ( int i = 0; i < npts; i++ ){ */
  /*   fprintf(fp, "%e %e %e\n",  */
  /*           ctrl_pts[i].x, ctrl_pts[i].y, ctrl_pts[i].z);  */
  /* } */

  /* fprintf(fp, "Zone T = interpolation I=%d\n", nctl); */
  /* for ( int i = 0; i < nctl; i++ ){ */
  /*   fprintf(fp, "%e %e %e\n",  */
  /*           line_pts[i].x, line_pts[i].y, line_pts[i].z);  */
  /* } */

  fclose(fp);

  // Free the objects
  surf->decref();
  for ( int k = 0; k < num_curves; k++ ){
    curves[k]->decref();
  }
}

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

  test_surface_lofter();

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
