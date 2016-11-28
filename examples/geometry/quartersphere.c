#include "TMRGeometry.h"
#include "TMROctForest.h"

/*
  Create the primitive geometry classes required for this example 
*/
class QuarterSphere : public TMRGeoSurface {
 public:
  QuarterSphere( double _R ){ R = _R; }
  
  // Get the parameter range for this surface
  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax ){
    *umin = 0.0;  *umax = 0.5*M_PI;
    *vmin = 0.0;  *vmax = 0.5*M_PI;
  }

  // Given the parametric point, compute the x,y,z location
  int evalPoint( double u, double v, TMRPoint *X ){
    X->x = R*cos(u)*cos(v);
    X->y = R*sin(u)*cos(v);
    X->z = R*sin(v);

    return 0;
  }

  // Given the parametric point, evaluate the first derivative 
  int evalDeriv( double u, double v, 
                 TMRPoint *Xu, TMRPoint *Xv ){
    Xu->x = -R*sin(u)*cos(v);
    Xu->y = R*cos(u)*cos(v);
    Xu->z = 0.0;
 
    Xv->x = -R*cos(u)*sin(v);
    Xv->y = -R*sin(u)*sin(v);
    Xv->z = R*cos(v);

    return 0;
  }

 private:
  double R;
};

/*
  Object for the parametrization of the side of the quarter sphere
*/
class QuarterSphereSide : public TMRGeoSurface {
 public:
  QuarterSphereSide( double _R, TMRPoint *_x1, TMRPoint *_x2 ){ 
    R = _R; 
    x1 = *_x1;
    x2 = *_x2;
  }

  // Get the parameter range for this surface
  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax ){
    *umin = 0.0;  *umax = R;
    *vmin = 0.0;  *vmax = 0.5*M_PI;
  }

  // Given the parametric point, compute the x,y,z location
  int evalPoint( double u, double v, TMRPoint *X ){
    X->x = u*(x1.x*cos(v) + x2.x*sin(v));
    X->y = u*(x1.y*cos(v) + x2.y*sin(v));
    X->z = u*(x1.z*cos(v) + x2.z*sin(v));

    return 0;
  }

  // Given the parametric point, evaluate the first derivative 
  int evalDeriv( double u, double v, 
                 TMRPoint *Xu, TMRPoint *Xv ){
    Xu->x = (x1.x*cos(v) + x2.x*sin(v));
    Xu->y = (x1.y*cos(v) + x2.y*sin(v));
    Xu->z = (x1.z*cos(v) + x2.z*sin(v));

    Xv->x = u*(-x1.x*sin(v) + x2.x*cos(v));
    Xv->y = u*(-x1.y*sin(v) + x2.y*cos(v));
    Xv->z = u*(-x1.z*sin(v) + x2.z*cos(v));

    return 0;
  }

 private:
  double R;
  TMRPoint x1, x2;
};

/*
  Object for the parametrization of the edge
*/
class LinearParamEdge : public TMRGeoEdge {
 public:
  LinearParamEdge( TMRGeoSurface *_surf,
             double _u1, double _v1,
             double _u2, double _v2 ){
    surf = _surf;
    u1 = _u1;
    v1 = _v1;
    u2 = _u2;
    v2 = _v2;
  }

  // Get the parameter range for this edge
  void getRange( double *tmin, double *tmax ){
    *tmin = 0.0; *tmax = 1.0;
  }
  
  // Given the parametric point, evaluate the x,y,z location
  int evalPoint( double t, TMRPoint *X ){
    // Evaluate the coordinates
    double u = u1 + (u2 - u1)*t;
    double v = v1 + (v2 - v1)*t;
    
    // Evaluate point
    surf->evalPoint(u, v, X);
    return 0;
  }

  // Given the parametric point, evaluate the derivative 
  int evalDeriv( double t, TMRPoint *Xt ){
    // The derivatives of the surface along the u/v directions
    TMRPoint Xu, Xv;

    // Evaluate the coordinates
    double du = u2 - u1;
    double dv = v2 - v1;
    double u = u1 + du*t;
    double v = v1 + dv*t;
    surf->evalDeriv(u, v, &Xu, &Xv);
    
    // Evaluate the derivative Xt = Xu*u,t + X,v*v,t
    Xt->x = Xu.x*du + Xv.x*dv;
    Xt->y = Xu.y*du + Xv.y*dv;
    Xt->z = Xu.z*du + Xv.z*dv;
    return 0;
  }

  // Find the surface u,v coordinates of the curve p(t)
  int reparamOnSurface( TMRGeoSurface *surface, 
                        double t, int dir,
                        double *u, double *v ){
    *u = *v = 0.0;
    if (surface == surf){
      if (dir >= 0){
        *u = u1 + (u2 - u1)*t;
        *v = v1 + (v2 - v1)*t;
      }
      else {
        *u = u2 + (u1 - u2)*t;
        *v = v2 + (v1 - v2)*t;
      }
      return 0;
    }
    return 1;
  }

 private:
  // The initial and final points on the surface
  double u1, v1, u2, v2;

  // The underlying surface object
  TMRGeoSurface *surf;
};

/*
  A linear edge between two points
*/
class LinearEdge : public TMRGeoEdge {
 public:
  LinearEdge( TMRPoint *_p1, TMRPoint *_p2 ){
    p1 = *_p1;  p2 = *_p2;
  }

  // Get the parameter range for this edge
  void getRange( double *tmin, double *tmax ){
    *tmin = 0.0; *tmax = 1.0;
  }
  
  // Given the parametric point, evaluate the x,y,z location
  int evalPoint( double t, TMRPoint *X ){
    X->x = (1.0 - t)*p1.x + t*p2.x;
    X->y = (1.0 - t)*p1.y + t*p2.y;
    X->z = (1.0 - t)*p1.z + t*p2.z;
    return 0;
  }

  // Given the parametric point, evaluate the derivative 
  int evalDeriv( double t, TMRPoint *Xt ){    
    Xt->x = p2.x - p1.x;
    Xt->y = p2.y - p1.y;
    Xt->z = p2.z - p1.z;
    return 0;
  }

  // Find the surface u,v coordinates of the curve p(t)
  int reparamOnSurface( TMRGeoSurface *surface, 
                        double t, int dir,
                        double *u, double *v ){
    // Find u,v such that psurf(u,v) = pedge(t)
    return 1;
  }

 private:
  TMRPoint p1, p2;
};

/*
  Generate the volume mesh from 
*/
int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Create the quarter sphere
  double R = 1.0;
  QuarterSphere *sphere = new QuarterSphere(R);

  // Create the quarter sphere sides
  TMRPoint e0(0.0, 0.0, 0.0);
  TMRPoint e1(1.0, 0.0, 0.0);
  TMRPoint e2(1.0, 0.0, 0.0);
  TMRPoint e3(1.0, 0.0, 0.0);
  QuarterSphereSide *side1 = new QuarterSphereSide(R, &e1, &e3);
  QuarterSphereSide *side2 = new QuarterSphereSide(R, &e1, &e2);
  QuarterSphereSide *side3 = new QuarterSphereSide(R, &e3, &e2);

  // Create the edges for the quarter sphere
  double end = 0.5*M_PI;
  LinearParamEdge *edge1 = new LinearParamEdge(sphere, 0.0, 0.0, end, 0.0);
  LinearParamEdge *edge2 = new LinearParamEdge(sphere, 0.0, 0.0, 0.0, end);
  LinearParamEdge *edge3 = new LinearParamEdge(sphere, end, end, end, 0.0);
    
  // Create the edges for the quarter sphere
  LinearEdge *edge4 = new LinearEdge(&e0, &e1);
  LinearEdge *edge5 = new LinearEdge(&e0, &e2);
  LinearEdge *edge6 = new LinearEdge(&e0, &e3);

  // Create the nodes that are needed for the model
  TMRGeoVertex vert0(&e0);
  TMRGeoVertex vert1(&e1);
  TMRGeoVertex vert2(&e2);
  TMRGeoVertex vert3(&e3);

  // Set the node range
  const double tol = 1e-6;
  double t1, t2;
  edge3->getRange(&t1, &t2);

  // Evaluate the integral and the distance function
  int nvals;
  double *tvals, *dist;
  edge3->integrate(t1, t2, tol, &tvals, &dist, &nvals);
  
  for ( int i = 0; i < nvals; i++ ){
    printf("[%3d] t: %15.5e   dist: %15.5e\n",
           i, tvals[i], dist[i]);
  }

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
