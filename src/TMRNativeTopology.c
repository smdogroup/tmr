#include "TMRNativeTopology.h"
#include <stdio.h>

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
TMRVertexFromEdge::TMRVertexFromEdge( TMREdge *_edge, 
                                      double _t ){
  t = _t;
  edge = _edge;
  edge->incref();
  setAttribute(edge->getAttribute());
}

/*
  Evaluate the vertex based on a node location
*/
TMRVertexFromEdge::TMRVertexFromEdge( TMREdge *_edge,
                                      TMRPoint p ){
  edge = _edge;
  edge->incref();
  setAttribute(edge->getAttribute());
  
  // Determine the parametric location of p using the initial
  // position
  edge->invEvalPoint(p, &t);
}

/*
  Free the object
*/
TMRVertexFromEdge::~TMRVertexFromEdge(){
  edge->decref();
}

/*
  Evaluate the point
*/
int TMRVertexFromEdge::evalPoint( TMRPoint *p ){
  return edge->evalPoint(t, p);
}

/*
  Retrieve the underlying curve object
*/
TMREdge* TMRVertexFromEdge::getEdge(){
  return edge;
}

/*
  Get the underlying parametric point
*/
int TMRVertexFromEdge::getParamOnEdge( TMREdge *_edge, 
                                       double *_t ){
  int fail = 0;
  if (edge == _edge){
    *_t = t;
    return fail;
  }
  return TMRVertex::getParamOnEdge(_edge, _t);
}

/*
  Get the underlying parametric point (if any)
*/
int TMRVertexFromEdge::getParamsOnFace( TMRFace *face,
                                        double *u, double *v ){
  return edge->getParamsOnFace(face, t, 1, u, v);
}

/*
  Determine the vertex location based on a surface location
*/
TMRVertexFromFace::TMRVertexFromFace( TMRFace *_face, 
                                      double _u, double _v ){
  face = _face;
  face->incref();
  setAttribute(face->getAttribute());
  u = _u;
  v = _v;
}

/*
  First determine the parametric vertex locations by projecting
  the point onto the surface.
*/
TMRVertexFromFace::TMRVertexFromFace( TMRFace *_face, 
                                      TMRPoint p ){
  face = _face;
  face->incref();
  setAttribute(face->getAttribute());
  face->invEvalPoint(p, &u, &v);
}

/*
  Free the data
*/
TMRVertexFromFace::~TMRVertexFromFace(){
  face->decref();
}

/*
  Evaluate the point on the surface
*/
int TMRVertexFromFace::evalPoint( TMRPoint *p ){
  return face->evalPoint(u, v, p);
}

/*
  Get the underlying parametric point (if any)
*/
int TMRVertexFromFace::getParamsOnFace( TMRFace *_face,
                                        double *_u, double *_v ){
  if (face == _face){
    *_u = u;
    *_v = v;
    return 1;
  }
  return TMRVertex::getParamsOnFace(face, _u, _v);
}

/*
  Create the curve parametrized on the surface
*/
TMREdgeFromFace::TMREdgeFromFace( TMRFace *_face, 
                                  TMRPcurve *_pcurve ){
  face = _face;
  face->incref();
  setAttribute(face->getAttribute());
  pcurve = _pcurve;
  pcurve->incref();
}

/*
  Destroy the curve
*/
TMREdgeFromFace::~TMREdgeFromFace(){
  face->decref();
  pcurve->decref();
}

/*
  Get the parameter range for this curve
*/
void TMREdgeFromFace::getRange( double *tmin, double *tmax ){
  pcurve->getRange(tmin, tmax);
}  
  
/*
  Given the parametric point, evaluate the x,y,z location
*/
int TMREdgeFromFace::evalPoint( double t, TMRPoint *X ){
  double u, v;
  int fail = pcurve->evalPoint(t, &u, &v);
  fail = fail || face->evalPoint(u, v, X);
  return fail;
}

/*
  Parametrize the curve on the given surface
*/
int TMREdgeFromFace::getParamsOnFace( TMRFace *surf, 
                                      double t, int dir, 
                                      double *u, double *v ){
  if (surf == face){
    return pcurve->evalPoint(t, u, v);
  }

  TMRPoint p;
  int fail = evalPoint(t, &p);
  fail = fail || surf->invEvalPoint(p, u, v);
  return fail;
}

/*
  Given the point, find the parametric location
*/
int TMREdgeFromFace::invEvalPoint( TMRPoint X, double *t ){
  *t = 0.0;
  int fail = 1;
  return fail;
}

/*
  Given the parametric point, evaluate the derivative 
*/
int TMREdgeFromFace::evalDeriv( double t, TMRPoint *Xt ){
  double u, v, ut, vt;
  pcurve->evalPoint(t, &u, &v);
  pcurve->evalDeriv(t, &ut, &vt);
  TMRPoint Xu, Xv;
  face->evalDeriv(u, v, &Xu, &Xv);
  Xt->x = ut*Xu.x + vt*Xv.x;
  Xt->y = ut*Xu.y + vt*Xv.y;
  Xt->z = ut*Xu.z + vt*Xv.z;
}

/*
  Split/segment the curve 
*/
TMRSplitEdge::TMRSplitEdge( TMREdge *_edge, 
                            double _t1, double _t2 ){
  edge = _edge;
  edge->incref();
  setAttribute(edge->getAttribute());

  // Set the parameter values
  t1 = _t1;
  t2 = _t2;

  // Check the range
  double tmin, tmax;
  edge->getRange(&tmin, &tmax);
  if (t1 < tmin){ t1 = tmin; }
  else if (t1 > tmax){ t1 = tmax; }
  if (t2 > tmax){ t2 = tmax; }
  else if (t2 < tmin){ t2 = tmin; }
}

/*
  Split the curve to the nearest points
*/
TMRSplitEdge::TMRSplitEdge( TMREdge *_edge, 
                            TMRPoint *p1, TMRPoint *p2 ){
  edge = _edge;
  edge->incref();
  setAttribute(edge->getAttribute());

  // Perform the inverse evaluation
  edge->invEvalPoint(*p1, &t1);
  edge->invEvalPoint(*p2, &t2);

  // Check the range
  double tmin, tmax;
  edge->getRange(&tmin, &tmax);
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
TMRSplitEdge::TMRSplitEdge( TMREdge *_edge, 
                            TMRVertex *v1, TMRVertex *v2 ){
  edge = _edge;
  edge->incref();
  setAttribute(edge->getAttribute());

  // Get the parameters for this curve for the point v1/v2
  v1->getParamOnEdge(edge, &t1);
  v2->getParamOnEdge(edge, &t2);

  // Set the vertices for this curve
  setVertices(v1, v2);

  // Check the range
  double tmin, tmax;
  edge->getRange(&tmin, &tmax);
  if (t1 < tmin){ t1 = tmin; }
  else if (t1 > tmax){ t1 = tmax; }
  if (t2 > tmax){ t2 = tmax; }
  else if (t2 < tmin){ t2 = tmin; } 
}

/*
  Decrease the reference count
*/
TMRSplitEdge::TMRSplitEdge(){
  edge->decref();
}

/*
  Get the parameter range
*/
void TMRSplitEdge::getRange( double *tmin, double *tmax ){
  *tmin = 0.0;
  *tmax = 1.0;
}

/*
  Evaluate the point
*/
int TMRSplitEdge::evalPoint( double t, TMRPoint *X ){
  int fail = 1;
  if (t < 0.0){ return fail; }
  if (t > 1.0){ return fail; }
  t = (1.0 - t)*t1 + t*t2;
  return edge->evalPoint(t, X);
}

/*
  Get the parameter on the split curve
*/
int TMRSplitEdge::getParamsOnFace( TMRFace *face, double t, 
                                   int dir, double *u, double *v ){
  int fail = 1;
  if (t < 0.0){ return fail; }
  if (t > 1.0){ return fail; }
  t = (1.0 - t)*t1 + t*t2;
  return edge->getParamsOnFace(face, t, dir, u, v);
}

/*
  Create a parametric TFI

  The transfinite interpolation is performed in the parameter space
  and all points are obtained directly from the surface object
  itself.
*/
TMRParametricTFIFace::TMRParametricTFIFace( TMRFace *_face, 
                                            TMREdge *_edges[], 
                                            const int _dir[],
                                            TMRVertex *verts[] ){
  face = _face;
  face->incref();
  setAttribute(face->getAttribute());

  for ( int k = 0; k < 4; k++ ){
    // Retrieve the parametric curves on the surface
    edges[k] = _edges[k];
    edges[k]->incref();
    dir[k] = _dir[k];
  
    double tmin = 0.0, tmax = 0.0;
    edges[k]->getRange(&tmin, &tmax);
    if (tmin != 0.0 || tmax != 1.0){
      fprintf(stderr, 
              "TMRParametricTFIFace error: All edges must have t in [0, 1]\n");
    }

    // Reparametrize the vertices on the surface
    verts[k]->getParamsOnFace(face, &vupt[k], &vvpt[k]);
  }

  // Set the vertices for this surface
  if (dir[0] > 0){
    edges[0]->setVertices(verts[0], verts[2]); 
  }
  else {
    edges[0]->setVertices(verts[2], verts[0]);
  }

  if (dir[1] > 0){  
    edges[1]->setVertices(verts[1], verts[3]);
  }
  else {
    edges[1]->setVertices(verts[3], verts[1]);
  }

  if (dir[2] > 0){
    edges[2]->setVertices(verts[0], verts[1]);
  }
  else {
    edges[2]->setVertices(verts[1], verts[0]);
  }

  if (dir[3] > 0){
    edges[3]->setVertices(verts[2], verts[3]);
  }
  else {
    edges[3]->setVertices(verts[3], verts[2]);
  }

  // Set the curve segment into the curve
  TMREdge *c[4];
  int d[4];
  c[0] = edges[2];  d[0] = dir[2];
  c[1] = edges[1];  d[1] = dir[1];
  c[2] = edges[3];  d[2] = -dir[3];
  c[3] = edges[0];  d[3] = -dir[0];
  TMREdgeLoop *loop = new TMREdgeLoop(4, c, d);
  addEdgeLoop(loop);
}

/*
  Destroy the parametric TFI object
*/  
TMRParametricTFIFace::~TMRParametricTFIFace(){
  face->decref();
  for ( int k = 0; k < 4; k++ ){
    edges[k]->decref();
  }
}

/*
  The range must always be between [0,1] for all curves
*/
void TMRParametricTFIFace::getRange( double *umin, double *vmin,
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
int TMRParametricTFIFace::evalPoint( double u, double v, 
                                     TMRPoint *X ){
  // Evaluate the curves along the v-direction
  int fail = 0;
  double cupt[4], cvpt[4];
  double params[4] = {v, v, u, u};

  for ( int k = 0; k < 4; k++ ){
    if (dir[k] > 0){
      fail = fail || 
        edges[k]->getParamsOnFace(face, params[k], dir[k],
                                  &cupt[k], &cvpt[k]);
    }
    else {
      fail = fail || 
        edges[k]->getParamsOnFace(face, 1.0 - params[k], dir[k],
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

  fail = fail || face->evalPoint(us, vs, X);
  return fail;
} 

/*  
  Inverse evaluation: This is not yet implemented
*/
int TMRParametricTFIFace::invEvalPoint( TMRPoint p, 
                                        double *u, double *v ){
  *u = 0.0; 
  *v = 0.0;
  int fail = 1;
  return fail;
}

/*
  Derivative evaluation: This is not yet implemented
*/
int TMRParametricTFIFace::evalDeriv( double u, double v, 
                                     TMRPoint *Xu, TMRPoint *Xv ){

  // Evaluate the curves along the v-direction
  int fail = 1;
  Xu->zero();
  Xv->zero();

  return fail;
}
