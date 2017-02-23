#include "TMRGeometry.h"
#include "TMRMesh.h"
#include <math.h>
#include <stdio.h>

/*
  Perform an inverse evaluation by obtaining the underlying 
  parametrization on the specified curve 
*/
int TMRVertex::getParamsOnCurve( TMRCurve *curve, double *t ){
  TMRPoint p;
  int fail = evalPoint(&p);
  fail = fail || curve->invEvalPoint(p, t);
  return fail;
}

/*
  Same thing, except on the specified surface
*/
int TMRVertex::getParamsOnSurface( TMRSurface *surface,
                                   double *u, double *v ){
  TMRPoint p;
  int fail = evalPoint(&p);
  fail = fail || surface->invEvalPoint(p, u, v);
  return fail;
}

/*
  Set/retrieve vertex numbers
*/
int TMRVertex::setNodeNum( int *num ){
  if (var == -1){
    var = *num;
    (*num)++;
    return 1;
  }
  return 0;
}

/*
  Retrieve the vertex number
*/
int TMRVertex::getNodeNum( int *num ){
  if (var != -1){
    *num = var;
    return 1;
  }
  return 0;
}


/*
  Build the curve without specifying the start/end vertices
*/
TMRCurve::TMRCurve(){
  v1 = v2 = NULL;
  mesh = NULL;
}

/*
  Build the curve the specified starting/ending vertices
*/
TMRCurve::TMRCurve( TMRVertex *_v1, TMRVertex *_v2 ){
  v1 = _v1;  v1->incref(); 
  v2 = _v2;  v2->incref();
  mesh = NULL;
}

/*
  Free the curve object
*/
TMRCurve::~TMRCurve(){
  if (v1){ v1->decref(); }
  if (v2){ v2->decref(); }
  if (mesh){ mesh->decref(); }
}

/*
  Find the point on the surface closest to the point C(t)
*/
int TMRCurve::getParamsOnSurface( TMRSurface *surface, double t, 
                                  int dir, double *u, double *v ){
  TMRPoint p;
  int fail = evalPoint(t, &p);
  fail = fail || surface->invEvalPoint(p, u, v);
  return fail;
}

/*
  Compute the inverse: This is not always required. By default it is
  not implemented. Derived classes can implement it if needed.
*/
int TMRCurve::invEvalPoint( TMRPoint X, double *t ){
  int fail = 1;
  *t = 0.0;

  return fail;
}

/*
  Set the step size for the derivative
*/
double TMRCurve::deriv_step_size = 1e-6;

/*
  Evaluate the derivative using a finite-difference step size
*/
int TMRCurve::evalDeriv( double t, TMRPoint *Xt ){
  int fail = 1;

  // Retrieve the parameter bounds for the curve
  double tmin, tmax; 
  getRange(&tmin, &tmax);

  if (t >= tmin && t <= tmax){
    // Evaluate the point at the original 
    TMRPoint p;
    fail = evalPoint(t, &p);
    if (fail){ return fail; }
    
    // Compute the approximate derivative using a forward
    // difference
    if (t + deriv_step_size <= tmax){
      TMRPoint p2;
      fail = evalPoint(t + deriv_step_size, &p2);
      if (fail){ return fail; }

      Xt->x = (p2.x - p.x)/deriv_step_size;
      Xt->y = (p2.y - p.y)/deriv_step_size;
      Xt->z = (p2.z - p.z)/deriv_step_size;
    }
    else if (t >= tmin + deriv_step_size){
      TMRPoint p2;
      fail = evalPoint(t - deriv_step_size, &p2);
      if (fail){ return fail; }

      Xt->x = (p.x - p2.x)/deriv_step_size;
      Xt->y = (p.y - p2.y)/deriv_step_size;
      Xt->z = (p.z - p2.z)/deriv_step_size;
    }
  }

  return fail;
}

/*
  Set the adjacent vertices
*/
void TMRCurve::setVertices( TMRVertex *_v1, TMRVertex *_v2 ){
  _v1->incref();
  _v2->incref();
  if (v1){ v1->decref(); }
  if (v2){ v2->decref(); }
  v1 = _v1;
  v2 = _v2;
}

/*
  Retrieve the adjacent vertices
*/
void TMRCurve::getVertices( TMRVertex **_v1, TMRVertex **_v2 ){
  if (_v1){ *_v1 = v1; }
  if (_v2){ *_v2 = v2; }
}

/*
  An integral entry for the linked list
*/
class IntegralPt {
 public:
  double t;
  double dist;
  IntegralPt *next;
};

/*
  Evaluate the distance between two points
*/
double pointDist( TMRPoint *a, TMRPoint *b ){
  return sqrt((a->x - b->x)*(a->x - b->x) + 
              (a->y - b->y)*(a->y - b->y) + 
              (a->z - b->z)*(a->z - b->z));
}

/*
  Recursive integration on an edge with an adaptive error control to
  ensure that the integral is computed with sufficient accuracy.

  input:
  t1, t2:  the limits of integration for this interval
  tol:     the absolute error measure
  ncalls:  the recursion depth
  pt:      pointer into the linked list
*/
void integrateEdge( TMRCurve *edge,
                    double t1, TMRPoint p1, double t2, double tol,
                    int ncalls, IntegralPt **_pt ){
  // Dereference the pointer to the integral point
  IntegralPt *pt = *_pt;

  // Find the mid point of the interval
  TMRPoint pmid;
  double tmid = 0.5*(t1 + t2);
  edge->evalPoint(tmid, &pmid);

  // Evaluate the point at the end of the interval
  TMRPoint p2; 
  edge->evalPoint(t2, &p2);
  
  // Evaluate the approximate integral contributions
  double int1 = pointDist(&p1, &pmid);
  double int2 = pointDist(&pmid, &p2);
  double int3 = pointDist(&p1, &p2);

  // Compute the integration error
  double error = fabs(int3 - int1 - int2);

  if (((ncalls > 5) && (error < tol)) || (ncalls > 20)){
    // Add the mid point
    pt->next = new IntegralPt;
    pt->next->dist = pt->dist + int1;
    pt->next->t = tmid;
    pt->next->next = NULL;
    pt = pt->next;

    // Add the final point p2
    pt->next = new IntegralPt;
    pt->next->dist = pt->dist + int2;
    pt->next->t = t2;
    pt->next->next = NULL;
    pt = pt->next;

    // Set the pointer to the end of the linked list
    *_pt = pt;
  }
  else {
    // Continue the recursive integration
    integrateEdge(edge, t1, p1, tmid, tol, ncalls+1, _pt);
    integrateEdge(edge, tmid, pmid, t2, tol, ncalls+1, _pt);
  }
}

/*
  Integrate along the edge adaptively, creating a list 
*/
double TMRCurve::integrate( double t1, double t2, double tol,
                            double **_tvals, double **_dist, 
                            int *_nvals ){
  *_tvals = NULL;
  *_dist = NULL;
  *_nvals = 0;

  // Allocate the entry in the linked list
  IntegralPt *root = new IntegralPt;
  root->next = NULL;
  root->dist = 0.0;
  root->t = t1;

  // Evaluate the first point
  TMRPoint p1;
  evalPoint(t1, &p1);

  // Integrate over the edge
  IntegralPt *pt = root;
  integrateEdge(this, t1, p1, t2, tol, 0, &pt);

  // Count up and allocate the num
  int count = 1;
  IntegralPt *curr = root;
  while (curr->next){
    curr = curr->next;
    count++;
  }

  // Allocate arrays to store the parametric location/distance data
  double *tvals = new double[ count ];
  double *dist = new double[ count ];

  // Scan through the linked list, read out the values of the
  // parameter and its integral and delete the entries as we go...
  count = 0;
  curr = root;
  tvals[count] = curr->t;
  dist[count] = curr->dist;
  count++;

  while (curr->next){
    IntegralPt *tmp = curr;
    curr = curr->next;
    tvals[count] = curr->t;
    dist[count] = curr->dist;
    count++;
    delete tmp;
  }

  double len = curr->dist;
  delete curr;

  // Set the pointers for the output
  *_nvals = count;
  *_tvals = tvals;
  *_dist = dist;

  return len;
}

/*
  Set the mesh into the array
*/
void TMRCurve::setMesh( TMRCurveMesh *_mesh ){
  _mesh->incref();
  if (mesh){ mesh->decref(); }
  mesh = _mesh;
}

/*
  Retrieve the mesh pointer
*/
void TMRCurve::getMesh( TMRCurveMesh **_mesh ){
  *_mesh = mesh;
}

/*
  Write out a representation of the curve to a VTK file
*/
void TMRCurve::writeToVTK( const char *filename ){
  double t1, t2;
  getRange(&t1, &t2);

  const int npts = 100;

  // Write out the vtk file
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    
    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts);
    for ( int k = 0; k < npts; k++ ){
      double u = 1.0*k/(npts-1);
      double t = (1.0-u)*t1 + u*t2;

      // Evaluate the point
      TMRPoint p;
      evalPoint(t, &p);
      
      // Write out the point
      fprintf(fp, "%e %e %e\n", p.x, p.y, p.z);
    }
    
    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", npts-1, 3*(npts-1));
    for ( int k = 0; k < npts-1; k++ ){
      fprintf(fp, "2 %d %d\n", k, k+1);
    }
    
    // Write out the cell types
    fprintf(fp, "\nCELL_TYPES %d\n", npts-1);
    for ( int k = 0; k < npts-1; k++ ){
      fprintf(fp, "%d\n", 3);
    }
    
    fclose(fp);
  } 
}

/*
  Initialize data within the TMRSurface object
*/
TMRSurface::TMRSurface(){
  max_num_segments = 0;
  num_segments = 0;
  segments = NULL;
  mesh = NULL;
}

/*
  Deallocate the curve segments (if any)
*/
TMRSurface::~TMRSurface(){
  if (segments){
    for ( int i = 0; i < num_segments; i++ ){
      for ( int j = 0; j < segments[i]->num_curves; j++ ){
        segments[i]->curves[j]->decref();
      }
      delete [] segments[i]->curves;
      delete [] segments[i]->dir;
    }
    delete [] segments;
  }
  if (mesh){
    mesh->decref();
  }
}

/*
  Set the step size for the derivative
*/
double TMRSurface::deriv_step_size = 1e-6;

/*
  Evaluate the derivative using a finite-difference step size
*/
int TMRSurface::evalDeriv( double u, double v, 
                           TMRPoint *Xu, TMRPoint *Xv ){
  int fail = 0;

  // Retrieve the parameter bounds for the curve
  double umin, vmin, umax, vmax;
  getRange(&umin, &vmin, &umax, &vmax);

  if (u >= umin && u <= umax &&
      v >= vmin && v <= vmax){
    // Evaluate the point at the original 
    TMRPoint p;
    fail = evalPoint(u, v, &p);

    // Compute the approximate derivative using a forward
    // difference or backward difference, depending on whether
    // the step is within the domain
    if (u + deriv_step_size <= umax){
      TMRPoint p2;
      fail = fail || evalPoint(u + deriv_step_size, v, &p2);

      Xu->x = (p2.x - p.x)/deriv_step_size;
      Xu->y = (p2.y - p.y)/deriv_step_size;
      Xu->z = (p2.z - p.z)/deriv_step_size;
    }
    else if (u >= umin + deriv_step_size){
      TMRPoint p2;
      fail = fail || evalPoint(u - deriv_step_size, v, &p2);

      Xu->x = (p.x - p2.x)/deriv_step_size;
      Xu->y = (p.y - p2.y)/deriv_step_size;
      Xu->z = (p.z - p2.z)/deriv_step_size;
    }
    else {
      fail = 1;
    }

    // Compute the approximate derivative using a forward
    // difference
    if (v + deriv_step_size <= vmax){
      TMRPoint p2;
      fail = fail || evalPoint(u, v + deriv_step_size, &p2);

      Xv->x = (p2.x - p.x)/deriv_step_size;
      Xv->y = (p2.y - p.y)/deriv_step_size;
      Xv->z = (p2.z - p.z)/deriv_step_size;
    }
    else if (v >= vmin + deriv_step_size){
      TMRPoint p2;
      fail = fail || evalPoint(u, v - deriv_step_size, &p2);

      Xv->x = (p.x - p2.x)/deriv_step_size;
      Xv->y = (p.y - p2.y)/deriv_step_size;
      Xv->z = (p.z - p2.z)/deriv_step_size;
    }
    else {
      fail = 1;
    }
  }

  return fail;
}

/*
  Add the curves that bound the surface
*/
int TMRSurface::addCurveSegment( int ncurves, TMRCurve **curves, 
                                 const int dir[] ){
  int fail = 0;
  if (ncurves == 0){
    fail = 1;
    fprintf(stderr, "TMRSurface::addCurveSegment: Zero length segment\n");
  }

  // First, check whether the loop is closed
  TMRVertex *vinit = NULL;
  TMRVertex *vnext;
  for ( int i = 0; i < ncurves; i++ ){
    TMRVertex *v1, *v2;
    if (dir[i] > 0){
      curves[i]->getVertices(&v1, &v2);
    }
    else {
      curves[i]->getVertices(&v2, &v1);
    }
    if (i == 0){ vinit = v1; }
    vnext = v2;
    if (i == ncurves-1){
      if (vinit != vnext){
        fprintf(stderr, 
                "TMRSurface::addCurveSegment: Curve segment must be closed\n");
        fail = 1;
      }
    }
  }

  // Return if the segment is not closed
  if (fail){
    return fail;
  }
  
  // Add the curves
  for ( int i = 0; i < ncurves; i++ ){
    curves[i]->incref();
  }

  if (num_segments >= max_num_segments){
    // max_num_segments = max(2*num_segments, 10);
    max_num_segments = (2*num_segments > 10 ? 2*num_segments : 10);

    // Allocate the new segment array
    TMRSegment **segs = new TMRSegment*[ max_num_segments ];

    // Copy over any existing segments
    if (num_segments > 0){
      memcpy(segs, segments, num_segments*sizeof(TMRSegment*));
      delete [] segments;
    }

    // Set the new segments array with the newly allocated/copied data
    // (if any)
    segments = segs;
  }

  // Set the new segment array
  segments[num_segments] = new TMRSegment;
  segments[num_segments]->num_curves = ncurves;
  segments[num_segments]->curves = new TMRCurve*[ ncurves ];
  segments[num_segments]->dir = new int[ ncurves ];
  memcpy(segments[num_segments]->curves, curves, ncurves*sizeof(TMRCurve*));
  memcpy(segments[num_segments]->dir, dir, ncurves*sizeof(int));
  num_segments++;
}

/*
  Get the number of closed segments
*/
int TMRSurface::getNumSegments(){
  return num_segments;
}

/*
  Retrieve the information from the given segment number
*/
int TMRSurface::getCurveSegment( int k, int *ncurves, 
                                 TMRCurve ***curves, 
                                 const int **dir ){
  int fail = 0;

  if (k >= 0 && k < num_segments){
    if (ncurves){ *ncurves = segments[k]->num_curves; }
    if (curves){ *curves = segments[k]->curves; }
    if (dir){ *dir = segments[k]->dir; }
    return fail;
  }

  if (ncurves){ *ncurves = 0; }
  if (curves){ *curves = NULL; }
  if (dir){ *dir = NULL; }

  // No curve found
  fail = 1;
  return fail;
}

/*
  Set the mesh into the array
*/
void TMRSurface::setMesh( TMRSurfaceMesh *_mesh ){
  _mesh->incref();
  if (mesh){ mesh->decref(); }
  mesh = _mesh;
}

/*
  Retrieve the mesh pointer
*/
void TMRSurface::getMesh( TMRSurfaceMesh **_mesh ){
  *_mesh = mesh;
}

/*
  Write out a representation of the surface to a VTK file
*/
void TMRSurface::writeToVTK( const char *filename ){
  double umin, vmin, umax, vmax;
  getRange(&umin, &vmin, &umax, &vmax);

  const int npts = 100;

  // Write out the vtk file
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    
    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts*npts);
    for ( int j = 0; j < npts; j++ ){
      for ( int i = 0; i < npts; i++ ){
        double u = 1.0*i/(npts-1);
        double v = 1.0*j/(npts-1);
        u = (1.0 - u)*umin + u*umax;
        v = (1.0 - v)*vmin + v*vmax;

        // Evaluate the point
        TMRPoint p;
        evalPoint(u, v, &p);
        
        // Write out the point
        fprintf(fp, "%e %e %e\n", p.x, p.y, p.z);
      }
    } 
    
    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", (npts-1)*(npts-1), 5*(npts-1)*(npts-1));
    for ( int j = 0; j < npts-1; j++ ){
      for ( int i = 0; i < npts-1; i++ ){
        fprintf(fp, "4 %d %d %d %d\n", 
                i + j*npts, i+1 + j*npts, 
                i+1 + (j+1)*npts, i + (j+1)*npts);
      }
    }
    
    // Write out the cell types
    fprintf(fp, "\nCELL_TYPES %d\n", (npts-1)*(npts-1));
    for ( int k = 0; k < (npts-1)*(npts-1); k++ ){
      fprintf(fp, "%d\n", 9);
    }
    
    fclose(fp);
  } 
}

/*
  Create a vertex from a point
*/
TMRVertexFromPoint::TMRVertexFromPoint( TMRPoint p ){
  pt = p;
}

/*
  Read out the point
*/
int TMRVertexFromPoint::evalPoint( TMRPoint *p ){
  *p = pt;
  return 0;
}

/*
  Create a vertex from curve
*/
TMRVertexFromCurve::TMRVertexFromCurve( TMRCurve *_curve, 
                                        double _t ){
  t = _t;
  curve = _curve;
  curve->incref();
  setAttribute(curve->getAttribute());
}

/*
  Evaluate the vertex based on a node location
*/
TMRVertexFromCurve::TMRVertexFromCurve( TMRCurve *_curve,
                                        TMRPoint p ){
  curve = _curve;
  curve->incref();
  setAttribute(curve->getAttribute());
  
  // Determine the parametric location of p using the initial
  // position
  curve->invEvalPoint(p, &t);
}

/*
  Free the object
*/
TMRVertexFromCurve::~TMRVertexFromCurve(){
  curve->decref();
}

/*
  Evaluate the point
*/
int TMRVertexFromCurve::evalPoint( TMRPoint *p ){
  return curve->evalPoint(t, p);
}

/*
  Retrieve the underlying curve object
*/
TMRCurve* TMRVertexFromCurve::getCurve(){
  return curve;
}

/*
  Get the underlying parametric point
*/
int TMRVertexFromCurve::getParamsOnCurve( TMRCurve *_curve, 
                                          double *_t ){
  int fail = 0;
  if (curve == _curve){
    *_t = t;
    return fail;
  }
  return TMRVertex::getParamsOnCurve(_curve, _t);
}

/*
  Get the underlying parametric point (if any)
*/
int TMRVertexFromCurve::getParamsOnSurface( TMRSurface *surface,
                                            double *u, double *v ){
  curve->getParamsOnSurface(surface, t, 1, u, v);
}

/*
  Determine the vertex location based on a surface location
*/
TMRVertexFromSurface::TMRVertexFromSurface( TMRSurface *_surface, 
                                            double _u, double _v ){
  surface = _surface;
  surface->incref();
  setAttribute(surface->getAttribute());
  u = _u;
  v = _v;
}

/*
  First determine the parametric vertex locations by projecting
  the point onto the surface.
*/
TMRVertexFromSurface::TMRVertexFromSurface( TMRSurface *_surface, 
                                            TMRPoint p ){
  surface = _surface;
  surface->incref();
  setAttribute(surface->getAttribute());
  surface->invEvalPoint(p, &u, &v);
}

/*
  Free the data
*/
TMRVertexFromSurface::~TMRVertexFromSurface(){
  surface->decref();
}

/*
  Evaluate the point on the surface
*/
int TMRVertexFromSurface::evalPoint( TMRPoint *p ){
  return surface->evalPoint(u, v, p);
}

/*
  Get the underlying parametric point (if any)
*/
int TMRVertexFromSurface::getParamsOnSurface( TMRSurface *_surface,
                                              double *_u, double *_v ){
  if (surface == _surface){
    *_u = u;
    *_v = v;
    return 1;
  }
  return TMRVertex::getParamsOnSurface(surface, _u, _v);
}

/*
  Create the curve parametrized on the surface
*/
TMRCurveFromSurface::TMRCurveFromSurface( TMRSurface *_surface, 
                                          TMRPcurve *_pcurve ){
  surface = _surface;
  surface->incref();
  setAttribute(surface->getAttribute());
  pcurve = _pcurve;
  pcurve->incref();
}

/*
  Destroy the curve
*/
TMRCurveFromSurface::~TMRCurveFromSurface(){
  surface->decref();
  setAttribute(surface->getAttribute());
  pcurve->decref();
}

/*
  Get the parameter range for this curve
*/
void TMRCurveFromSurface::getRange( double *tmin, double *tmax ){
  pcurve->getRange(tmin, tmax);
}  
  
/*
  Given the parametric point, evaluate the x,y,z location
*/
int TMRCurveFromSurface::evalPoint( double t, TMRPoint *X ){
  double u, v;
  int fail = pcurve->evalPoint(t, &u, &v);
  fail = fail || surface->evalPoint(u, v, X);
  return fail;
}

/*
  Parametrize the curve on the given surface
*/
int TMRCurveFromSurface::getParamsOnSurface( TMRSurface *surf, 
                                             double t, int dir, 
                                             double *u, double *v ){
  if (surf == surface){
    return pcurve->evalPoint(t, u, v);
  }

  TMRPoint p;
  int fail = evalPoint(t, &p);
  fail = fail || surface->invEvalPoint(p, u, v);
  return fail;
}

/*
  Given the point, find the parametric location
*/
int TMRCurveFromSurface::invEvalPoint( TMRPoint X, double *t ){
  *t = 0.0;
  int fail = 1;
  return fail;
}

/*
  Given the parametric point, evaluate the derivative 
*/
int TMRCurveFromSurface::evalDeriv( double t, TMRPoint *Xt ){
  double u, v, ut, vt;
  pcurve->evalPoint(t, &u, &v);
  pcurve->evalDeriv(t, &ut, &vt);
  TMRPoint Xu, Xv;
  surface->evalDeriv(u, v, &Xu, &Xv);
  Xt->x = ut*Xu.x + vt*Xv.x;
  Xt->y = ut*Xu.y + vt*Xv.y;
  Xt->z = ut*Xu.z + vt*Xv.z;
}

/*
  Split/segment the curve 
*/
TMRSplitCurve::TMRSplitCurve( TMRCurve *_curve, 
                              double _t1, double _t2 ){
  curve = _curve;
  curve->incref();
  setAttribute(curve->getAttribute());

  // Set the parameter values
  t1 = _t1;
  t2 = _t2;

  // Check the range
  double tmin, tmax;
  curve->getRange(&tmin, &tmax);
  if (t1 < tmin){ t1 = tmin; }
  else if (t1 > tmax){ t1 = tmax; }
  if (t2 > tmax){ t2 = tmax; }
  else if (t2 < tmin){ t2 = tmin; }
}

/*
  Split the curve to the nearest points
*/
TMRSplitCurve::TMRSplitCurve( TMRCurve *_curve, 
                              TMRPoint *p1, TMRPoint *p2 ){
  curve = _curve;
  curve->incref();
  setAttribute(curve->getAttribute());

  // Perform the inverse evaluation
  curve->invEvalPoint(*p1, &t1);
  curve->invEvalPoint(*p2, &t2);

  // Check the range
  double tmin, tmax;
  curve->getRange(&tmin, &tmax);
  if (t1 < tmin){ t1 = tmin; }
  else if (t1 > tmax){ t1 = tmax; }
  if (t2 > tmax){ t2 = tmax; }
  else if (t2 < tmin){ t2 = tmin; }
}

/*
  Split the curve and check if the two vertices are evalauted from
  this curve using a parametric location. Otherwise do the same as
  before.
*/
TMRSplitCurve::TMRSplitCurve( TMRCurve *_curve, 
                              TMRVertex *v1, TMRVertex *v2 ){
  curve = _curve;
  curve->incref();
  setAttribute(curve->getAttribute());

  // Get the parameters for this curve for the point v1/v2
  v1->getParamsOnCurve(curve, &t1);
  v2->getParamsOnCurve(curve, &t2);

  // Set the vertices for this curve
  setVertices(v1, v2);

  // Check the range
  double tmin, tmax;
  curve->getRange(&tmin, &tmax);
  if (t1 < tmin){ t1 = tmin; }
  else if (t1 > tmax){ t1 = tmax; }
  if (t2 > tmax){ t2 = tmax; }
  else if (t2 < tmin){ t2 = tmin; } 
}

/*
  Decrease the reference count
*/
TMRSplitCurve::TMRSplitCurve(){
  curve->decref();
}

/*
  Get the parameter range
*/
void TMRSplitCurve::getRange( double *tmin, double *tmax ){
  *tmin = 0.0;
  *tmax = 1.0;
}

/*
  Evaluate the point
*/
int TMRSplitCurve::evalPoint( double t, TMRPoint *X ){
  int fail = 1;
  if (t < 0.0){ return fail; }
  if (t > 1.0){ return fail; }
  t = (1.0 - t)*t1 + t*t2;
  return curve->evalPoint(t, X);
}

/*
  Get the parameter on the split curve
*/
int TMRSplitCurve::getParamsOnSurface( TMRSurface *surface, double t, 
                                       int dir, double *u, double *v ){
  int fail = 1;
  if (t < 0.0){ return fail; }
  if (t > 1.0){ return fail; }
  t = (1.0 - t)*t1 + t*t2;
  return curve->getParamsOnSurface(surface, t, dir, u, v);
}

/*
  Create a parametric TFI

  The transfinite interpolation is performed in the parameter space
  and all points are obtained directly from the surface object
  itself.
*/
TMRParametricTFISurface::TMRParametricTFISurface( TMRSurface *_surf, 
                                                  TMRCurve *_curves[], 
                                                  const int _dir[],
                                                  TMRVertex *verts[] ){
  surf = _surf;
  surf->incref();
  setAttribute(surf->getAttribute());

  for ( int k = 0; k < 4; k++ ){
    // Retrieve the parametric curves on the surface
    curves[k] = _curves[k];
    curves[k]->incref();
    dir[k] = _dir[k];
  
    double tmin = 0.0, tmax = 0.0;
    curves[k]->getRange(&tmin, &tmax);
    if (tmin != 0.0 || tmax != 1.0){
      fprintf(stderr, 
              "TMRParametricTFISurface error: All curves must have t in [0, 1]\n");
    }

    // Reparametrize the vertices on the surface
    verts[k]->getParamsOnSurface(surf, &vupt[k], &vvpt[k]);
  }

  // Set the vertices for this surface
  if (dir[0] > 0){
    curves[0]->setVertices(verts[0], verts[2]); 
  }
  else {
    curves[0]->setVertices(verts[2], verts[0]);
  }

  if (dir[1] > 0){  
    curves[1]->setVertices(verts[1], verts[3]);
  }
  else {
    curves[1]->setVertices(verts[3], verts[1]);
  }

  if (dir[2] > 0){
    curves[2]->setVertices(verts[0], verts[1]);
  }
  else {
    curves[2]->setVertices(verts[1], verts[0]);
  }

  if (dir[3] > 0){
    curves[3]->setVertices(verts[2], verts[3]);
  }
  else {
    curves[3]->setVertices(verts[3], verts[2]);
  }


  // Set the curve segment into the curve
  TMRCurve *c[4];
  int d[4];
  c[0] = curves[2];  d[0] = dir[2];
  c[1] = curves[1];  d[1] = dir[1];
  c[2] = curves[3];  d[2] = -dir[3];
  c[3] = curves[0];  d[3] = -dir[0];
  addCurveSegment(4, c, d);
}

/*
  Destroy the parametric TFI object
*/  
TMRParametricTFISurface::~TMRParametricTFISurface(){
  surf->decref();
  for ( int k = 0; k < 4; k++ ){
    curves[k]->decref();
  }
}

/*
  The range must always be between [0,1] for all curves
*/
void TMRParametricTFISurface::getRange( double *umin, double *vmin,
                                        double *umax, double *vmax ){
  *umin = 0.0;
  *vmin = 0.0;
  *umax = 1.0;
  *vmax = 1.0;
}

/*
  Evaluate the surface at the specified parametric point (u,v)

  This code uses a transfinite interpolation to obtain the 
  parametric surface coordinates (us(u,v), vs(u,v)) in terms 
  of the TFI parameter coordinates (u,v)
*/
int TMRParametricTFISurface::evalPoint( double u, double v, 
                                        TMRPoint *X ){
  // Evaluate the curves along the v-direction
  int fail = 0;
  double cupt[4], cvpt[4];
  double params[4] = {v, v, u, u};

  for ( int k = 0; k < 4; k++ ){
    if (dir[k] > 0){
      fail = fail || 
        curves[k]->getParamsOnSurface(surf, params[k], dir[k],
                                      &cupt[k], &cvpt[k]);
    }
    else {
      fail = fail || 
        curves[k]->getParamsOnSurface(surf, 1.0 - params[k], dir[k],
                                      &cupt[k], &cvpt[k]);
    }
  }
  
  // Compute the parametric coordinates
  double us, vs;
  us = (1.0-u)*cupt[0] + u*cupt[1] + (1.0-v)*cupt[2] + v*cupt[3]
    - ((1.0-u)*(1.0-v)*vupt[0] + u*(1.0-v)*vupt[1] + 
       v*(1.0-u)*vupt[2] + u*v*vupt[3]);

  vs = (1.0-u)*cvpt[0] + u*cvpt[1] + (1.0-v)*cvpt[2] + v*cvpt[3]
    - ((1.0-u)*(1.0-v)*vvpt[0] + u*(1.0-v)*vvpt[1] + 
       v*(1.0-u)*vvpt[2] + u*v*vvpt[3]);

  fail = fail || surf->evalPoint(us, vs, X);
  return fail;
} 

/*  
  Inverse evaluation: This is not yet implemented
*/
int TMRParametricTFISurface::invEvalPoint( TMRPoint p, 
                                           double *u, double *v ){
  *u = 0.0; 
  *v = 0.0;
  int fail = 1;
  return fail;
}

/*
  Derivative evaluation: This is not yet implemented
*/
int TMRParametricTFISurface::evalDeriv( double u, double v, 
                                        TMRPoint *Xu, TMRPoint *Xv ){

  // Evaluate the curves along the v-direction
  int fail = 1;
  Xu->zero();
  Xv->zero();

  return fail;
}

/*
  The TMRGeometry class containing all of the required geometry
  objects.
*/
TMRGeometry::TMRGeometry( int _num_vertices, TMRVertex **_vertices, 
                          int _num_curves, TMRCurve **_curves,
                          int _num_surfaces, TMRSurface **_surfaces ){
  num_vertices = _num_vertices;
  num_curves = _num_curves;
  num_surfaces = _num_surfaces;

  vertices = new TMRVertex*[ num_vertices ];
  curves = new TMRCurve*[ num_curves ];
  surfaces = new TMRSurface*[ num_surfaces ];

  for ( int i = 0; i < num_vertices; i++ ){
    vertices[i] = _vertices[i];
    vertices[i]->incref();
  }

  for ( int i = 0; i < num_curves; i++ ){
    curves[i] = _curves[i];
    curves[i]->incref();
  }

  for ( int i = 0; i < num_surfaces; i++ ){
    surfaces[i] = _surfaces[i];
    surfaces[i]->incref();
  }

  ordered_verts = new OrderedPair<TMRVertex>[ num_vertices ];
  ordered_curves = new OrderedPair<TMRCurve>[ num_curves ];
  ordered_surfaces = new OrderedPair<TMRSurface>[ num_surfaces ];

  for ( int i = 0; i < num_vertices; i++ ){
    ordered_verts[i].num = i;
    ordered_verts[i].obj = vertices[i];
  }

  for ( int i = 0; i < num_curves; i++ ){
    ordered_curves[i].num = i;
    ordered_curves[i].obj = curves[i];
  }

  for ( int i = 0; i < num_surfaces; i++ ){
    ordered_surfaces[i].num = i;
    ordered_surfaces[i].obj = surfaces[i];
  }

  // Sort the vertices, curves and surfaces
  qsort(ordered_verts, num_vertices, sizeof(OrderedPair<TMRVertex>),
        compare_ordered_pairs<TMRVertex>);
  qsort(ordered_curves, num_curves, sizeof(OrderedPair<TMRCurve>),
        compare_ordered_pairs<TMRCurve>);
  qsort(ordered_surfaces, num_surfaces, sizeof(OrderedPair<TMRSurface>),
        compare_ordered_pairs<TMRSurface>);

  verify();
}

/*
  Free the geometry objects
*/
TMRGeometry::~TMRGeometry(){
  for ( int i = 0; i < num_vertices; i++ ){
    vertices[i]->decref();
  }
  for ( int i = 0; i < num_curves; i++ ){
    curves[i]->decref();
  }
  for ( int i = 0; i < num_surfaces; i++ ){
    surfaces[i]->decref();
  }
  delete [] vertices;
  delete [] curves;
  delete [] surfaces;
  delete [] ordered_verts;
  delete [] ordered_curves;
  delete [] ordered_surfaces;
}

/*
  Verify that all objects that are referenced in the geometry have a
  verifiable list name and that objects are referred to one or more
  times (but not zero times!)
*/
int TMRGeometry::verify(){
  int fail = 0;
  int *verts = new int[ num_vertices ];
  int *crvs = new int[ num_curves ];
  memset(verts, 0, num_vertices*sizeof(int));
  memset(crvs, 0, num_curves*sizeof(int));

  for ( int face = 0; face < num_surfaces; face++ ){
    int nseg = surfaces[face]->getNumSegments();
    
    for ( int k = 0; k < nseg; k++ ){
      int ncurves;
      TMRCurve **curves;
      surfaces[face]->getCurveSegment(k, &ncurves, &curves, NULL);

      // Loop over all of the curves and check whether the data exists
      // or not
      for ( int j = 0; j < ncurves; j++ ){
        int cindex = getCurveIndex(curves[j]);
        if (cindex < 0){
          fail = 1;
          fprintf(stderr, 
                  "TMRGeometry error: Curve does not exist within curve list\n");
        }
        else {
          crvs[cindex]++;
        }

        TMRVertex *v1, *v2;      
        curves[j]->getVertices(&v1, &v2);
        if (!v1 || !v2){
          fail = 1;
          fprintf(stderr, 
                  "TMRGeometry error: Vertices not set for curve %d\n",
                  cindex);
        }

        int v1index = getVertexIndex(v1);
        int v2index = getVertexIndex(v2);
        if (v1index < 0 || v2index < 0){
          fail = 1;
          fprintf(stderr, 
                  "TMRGeometry error: Vertices do not exist within vertex list\n");
        }
        else {
          verts[v1index]++;
          verts[v2index]++;
        }
      }
    }
  }

  // Check if any of the counts are zero
  for ( int i = 0; i < num_vertices; i++ ){
    if (verts[i] == 0){
      fprintf(stderr,
              "TMRGeometry error: Vertex %d unreferenced\n", i);
      fail = 1;
    }
  }
  for ( int i = 0; i < num_curves; i++ ){
    if (crvs[i] == 0){
      fprintf(stderr,
              "TMRGeometry error: Curve %d unreferenced\n", i);
      fail = 1;
    }
  }

  delete [] verts;
  delete [] crvs;

  return fail;
}

/*
  Retrieve the vertices
*/
void TMRGeometry::getVertices( int *_num_vertices, 
                               TMRVertex ***_vertices ){
  if (_num_vertices){ *_num_vertices = num_vertices; }
  if (_vertices){ *_vertices = vertices; }
}

/*
  Retrieve the curves
*/
void TMRGeometry::getCurves( int *_num_curves, 
                             TMRCurve ***_curves ){
  if (_num_curves){ *_num_curves = num_curves; }
  if (_curves){ *_curves = curves; }
}

/*
  Retrieve the surfaces
*/
void TMRGeometry::getSurfaces( int *_num_surfaces, 
                               TMRSurface ***_surfaces ){
  if (_num_surfaces){ *_num_surfaces = num_surfaces; }
  if (_surfaces){ *_surfaces = surfaces; }
}

/*
  Static member function for sorting the ordered pairs
*/
template <class ctype>
int TMRGeometry::compare_ordered_pairs( const void *avoid, const void *bvoid ){
  const OrderedPair<ctype> *a = static_cast<const OrderedPair<ctype>*>(avoid);
  const OrderedPair<ctype> *b = static_cast<const OrderedPair<ctype>*>(bvoid);
  return a->obj - b->obj; // Use pointer arithmetic to determine relative positions
}

/*
  Retrieve the index given the vertex point
*/
int TMRGeometry::getVertexIndex( TMRVertex *vertex ){
  OrderedPair<TMRVertex> pair;
  pair.num = -1;
  pair.obj = vertex;

  // Search for the ordered pair
  OrderedPair<TMRVertex> *item = 
    (OrderedPair<TMRVertex>*)bsearch(&pair, ordered_verts, num_vertices, 
                                     sizeof(OrderedPair<TMRVertex>),
                                     compare_ordered_pairs<TMRVertex>);
  if (item){
    return item->num;
  }
  return -1;
}

/*
  Retrieve the index given the pointer to the curve object
*/
int TMRGeometry::getCurveIndex( TMRCurve *curve ){
  OrderedPair<TMRCurve> pair;
  pair.num = -1;
  pair.obj = curve;

  // Search for the ordered pair
  OrderedPair<TMRCurve> *item = 
    (OrderedPair<TMRCurve>*)bsearch(&pair, ordered_curves, num_curves,
                                    sizeof(OrderedPair<TMRCurve>),
                                    compare_ordered_pairs<TMRCurve>);
  if (item){
    return item->num;
  }
  return -1;
}

/*
  Retrieve the index given the pointer to the surface object
*/
int TMRGeometry::getSurfaceIndex( TMRSurface *surf ){
  OrderedPair<TMRSurface> pair;
  pair.num = -1;
  pair.obj = surf;

  // Search for the ordered pair
  OrderedPair<TMRSurface> *item = 
    (OrderedPair<TMRSurface>*)bsearch(&pair, ordered_surfaces, num_surfaces,
                                      sizeof(OrderedPair<TMRSurface>),
                                      compare_ordered_pairs<TMRSurface>);
  if (item){
    return item->num;
  }
  return -1;
}
