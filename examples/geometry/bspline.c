#include "TMRGeometry.h"
#include "TMRBspline.h"
#include "TMRMesh.h"
#include "TMRQuadForest.h"
#include <stdio.h>
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

void test_surface_lofter( double htarget ){
  // Create the control points
  int nctl = rae2822_npts;
  TMRPoint line_pts[rae2822_npts];
  memset(line_pts, 0, nctl*sizeof(TMRPoint));
  
  // Create a series of curves
  int num_curves = 5;
  TMRBsplineCurve *lofts[5];
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
    lofts[k] = interper->createCurve(4);
    lofts[k]->incref();

    // Free the interpolation object
    interper->decref();
  }  

  // Create the lofter object
  TMRCurveLofter *lofter = new TMRCurveLofter(lofts, num_curves);
  lofter->incref();

  // Create the surface object from the lofted surface
  TMRBsplineSurface *surface = lofter->createSurface(4);
  surface->incref();
  lofter->decref();

  surface->writeToVTK("bspline_surface.vtk");

  // (0,1) -- (1,1)
  //   |        |
  // (0,0) -- (1,0)
  double pts1[] = {0.1, 0.0, 0.4, 0.0};
  double pts2[] = {0.4, 0.0, 0.4, 1.0};
  double pts3[] = {0.4, 1.0, 0.1, 1.0};
  double pts4[] = {0.1, 1.0, 0.1, 0.0};
  TMRBsplinePcurve *p1 = new TMRBsplinePcurve(2, 2, pts1);
  TMRBsplinePcurve *p2 = new TMRBsplinePcurve(2, 2, pts2);
  TMRBsplinePcurve *p3 = new TMRBsplinePcurve(2, 2, pts3);
  TMRBsplinePcurve *p4 = new TMRBsplinePcurve(2, 2, pts4);

  // Create the curves and add them to the surface
  int ncurves = 4;
  TMRCurve *curves[4];
  curves[0] = new TMRCurveFromSurface(surface, p1);
  curves[1] = new TMRCurveFromSurface(surface, p2);
  curves[2] = new TMRCurveFromSurface(surface, p3);
  curves[3] = new TMRCurveFromSurface(surface, p4);

  // Create the boundary curves for the surface
  TMRVertexFromCurve *v1 = new TMRVertexFromCurve(curves[0], 0.0);
  TMRVertexFromCurve *v2 = new TMRVertexFromCurve(curves[1], 0.0);
  TMRVertexFromCurve *v3 = new TMRVertexFromCurve(curves[2], 0.0);
  TMRVertexFromCurve *v4 = new TMRVertexFromCurve(curves[3], 0.0);

  int num_verts = 4;
  TMRVertex *verts[] = {v1, v2, v3, v4};

  // Set the vertices
  curves[0]->setVertices(v1, v1);
  curves[1]->setVertices(v2, v3);
  curves[2]->setVertices(v3, v4);
  curves[3]->setVertices(v4, v1);
  
  // Set the directions of the curves
  int dir[4];
  dir[0] = 1;
  dir[1] = 1;
  dir[2] = 1;
  dir[3] = 1;

  surface->addCurveSegment(ncurves, curves, dir);

  // Create the TMRGeometry
  TMRSurface *surf = surface;
  TMRGeometry *geo = new TMRGeometry(num_verts, verts,
                                     num_curves, curves,
                                     1, &surf);

  // Allocate the new mesh
  TMRMesh *mesh = new TMRMesh(geo);

  // Mesh the geometry
  mesh->mesh(htarget);

  // Create the mesh
  TMRGeometry *geo_mesh = mesh->createMeshGeometry();
  geo_mesh->incref();

  TMRTopology *topo = new TMRTopology(geo_mesh);
  topo->incref();

  TMRQuadForest *forest = new TMRQuadForest(MPI_COMM_WORLD);
  forest->incref();

  // Set up the forest 
  forest->setTopology(topo);


  forest->decref();
  topo->incref();
  geo_mesh->decref();

  // Free the objects
  surface->decref();
  for ( int k = 0; k < num_curves; k++ ){
    lofts[k]->decref();
  }
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  double htarget = 0.1;
  for ( int i = 0; i < argc; i++ ){
    if (sscanf(argv[i], "h=%lf", &htarget) == 1){
      if (htarget < 0.01){ htarget = 0.01; }
      if (htarget > 1.0){ htarget = 1.0; }
    }
  }

  test_surface_lofter(htarget);

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
