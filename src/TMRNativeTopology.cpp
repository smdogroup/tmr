/*
  This file is part of the package TMR for adaptive mesh refinement.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#include "TMRNativeTopology.h"
#include <stdio.h>
#include <math.h>

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
  setName(edge->getName());
}

/*
  Evaluate the vertex based on a node location
*/
TMRVertexFromEdge::TMRVertexFromEdge( TMREdge *_edge,
                                      TMRPoint p ){
  edge = _edge;
  edge->incref();
  setName(edge->getName());

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
  setName(face->getName());
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
  setName(face->getName());
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
                                  TMRPcurve *_pcurve,
                                  int _is_degen ){
  nfaces = 1;
  faces = new TMRFace*[ nfaces ];
  pcurves = new TMRPcurve*[ nfaces ];
  faces[0] = _face;
  faces[0]->incref();
  pcurves[0] = _pcurve;
  pcurves[0]->incref();
  setName(faces[0]->getName());
  is_degen = _is_degen;
}

/*
  Destroy the curve
*/
TMREdgeFromFace::~TMREdgeFromFace(){
  for ( int i = 0; i < nfaces; i++ ){
    faces[i]->decref();
    pcurves[i]->decref();
  }
  delete [] faces;
  delete [] pcurves;
}

/*
  Get the parameter range for this curve
*/
void TMREdgeFromFace::getRange( double *tmin, double *tmax ){
  pcurves[0]->getRange(tmin, tmax);
}

/*
  Given the parametric point, evaluate the x,y,z location
*/
int TMREdgeFromFace::evalPoint( double t, TMRPoint *X ){
  double u, v;
  int fail = pcurves[0]->evalPoint(t, &u, &v);
  fail = fail || faces[0]->evalPoint(u, v, X);
  return fail;
}

/*
  Parametrize the curve on the given surface
*/
int TMREdgeFromFace::getParamsOnFace( TMRFace *surf,
                                      double t, int dir,
                                      double *u, double *v ){
  for ( int i = 0; i < nfaces; i++ ){
    if (surf == faces[i]){
      return pcurves[i]->evalPoint(t, u, v);
    }
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
int TMREdgeFromFace::evalDeriv( double t, TMRPoint *X, TMRPoint *Xt ){
  int fail = 0;
  double u, v, ut, vt;
  pcurves[0]->evalDeriv(t, &u, &v, &ut, &vt);
  TMRPoint Xu, Xv;
  faces[0]->evalDeriv(u, v, X, &Xu, &Xv);
  Xt->x = ut*Xu.x + vt*Xv.x;
  Xt->y = ut*Xu.y + vt*Xv.y;
  Xt->z = ut*Xu.z + vt*Xv.z;
  return fail;
}

/*
  Add an additional parametrization of the edge along the face
*/
void TMREdgeFromFace::addEdgeFromFace( TMRFace *_face,
                                       TMRPcurve *_pcurve ){
  nfaces++;
  TMRFace **_faces = new TMRFace*[ nfaces ];
  TMRPcurve **_pcurves = new TMRPcurve*[ nfaces ];
  for ( int i = 0; i < nfaces-1; i++ ){
    _faces[i] = faces[i];
    _pcurves[i] = pcurves[i];
  }
  _face->incref();
  _pcurve->incref();
  _faces[nfaces-1] = _face;
  _pcurves[nfaces-1] = _pcurve;
  delete [] faces;
  delete [] pcurves;
  faces = _faces;
  pcurves = _pcurves;
}

/*
  Split/segment the curve
*/
TMRSplitEdge::TMRSplitEdge( TMREdge *_edge,
                            double _t1, double _t2 ){
  edge = _edge;
  edge->incref();
  setName(edge->getName());

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
  setName(edge->getName());

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
  setName(edge->getName());

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
  t = (1.0 - t)*t1 + t*t2;
  return edge->evalPoint(t, X);
}

/*
  Get the parameter on the split curve
*/
int TMRSplitEdge::getParamsOnFace( TMRFace *face, double t,
                                   int dir, double *u, double *v ){
  t = (1.0 - t)*t1 + t*t2;
  return edge->getParamsOnFace(face, t, dir, u, v);
}

/*
  Create a transfinite interpolation edge
*/
TMRTFIEdge::TMRTFIEdge( TMRVertex *_v1, TMRVertex *_v2 ){
  setVertices(_v1, _v2);
}

TMRTFIEdge::~TMRTFIEdge(){}

/*
  Set the parameter range
*/
void TMRTFIEdge::getRange( double *tmin, double *tmax ){
  *tmin = 0.0;
  *tmax = 1.0;
}

/*
  Evaluate the point based on the vertex locations
*/
int TMRTFIEdge::evalPoint( double t, TMRPoint *X ){
  TMRVertex *_v1, *_v2;
  getVertices(&_v1, &_v2);

  // Evaluate the points
  TMRPoint p1, p2;;
  int f1 = _v1->evalPoint(&p1);
  int f2 = _v2->evalPoint(&p2);

  // Interpolate between the start/end locations
  X->x = (1.0 - t)*p1.x + t*p2.x;
  X->y = (1.0 - t)*p1.y + t*p2.y;
  X->z = (1.0 - t)*p1.z + t*p2.z;

  return f1 || f2;
}

/*
  Create a transfinite-interpolation (TFI) face from the given set of
  edges and vertices
*/
TMRTFIFace::TMRTFIFace( TMREdge *_edges[],
                        const int _dir[],
                        TMRVertex *verts[] ){
  for ( int k = 0; k < 4; k++ ){
    // Set the edge
    edges[k] = _edges[k];
    edges[k]->incref();
    dir[k] = _dir[k];

    // Evaluate the vertex points
    verts[k]->evalPoint(&c[k]);
  }

  // Set the vertices for this surface
  if (dir[0] > 0){
    edges[0]->setVertices(verts[0], verts[1]);
    edges[0]->getRange(&tmin[0], &tmax[0]);
  }
  else {
    edges[0]->setVertices(verts[1], verts[0]);
    edges[0]->getRange(&tmax[0], &tmin[0]);
  }

  if (dir[1] > 0){
    edges[1]->setVertices(verts[1], verts[2]);
    edges[1]->getRange(&tmin[1], &tmax[1]);
  }
  else {
    edges[1]->setVertices(verts[2], verts[1]);
    edges[1]->getRange(&tmax[1], &tmin[1]);
  }

  if (dir[2] > 0){
    edges[2]->setVertices(verts[2], verts[3]);
    edges[2]->getRange(&tmin[2], &tmax[2]);
  }
  else {
    edges[2]->setVertices(verts[3], verts[2]);
    edges[2]->getRange(&tmax[2], &tmin[2]);
  }

  if (dir[3] > 0){
    edges[3]->setVertices(verts[3], verts[0]);
    edges[3]->getRange(&tmin[3], &tmax[3]);
  }
  else {
    edges[3]->setVertices(verts[0], verts[3]);
    edges[3]->getRange(&tmax[3], &tmin[3]);
  }

  // Set the curve segment into the curve
  TMREdgeLoop *loop = new TMREdgeLoop(4, edges, dir);
  addEdgeLoop(1, loop);
}

/*
  Destroy the parametric TFI object
*/
TMRTFIFace::~TMRTFIFace(){
  for ( int k = 0; k < 4; k++ ){
    edges[k]->decref();
  }
}

/*
  Set the maximum number of Newton iterations
*/
int TMRTFIFace::max_newton_iters = 25;

/*
  The range must always be between [0,1] for all curves
*/
void TMRTFIFace::getRange( double *umin, double *vmin,
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
int TMRTFIFace::evalPoint( double u, double v,
                           TMRPoint *X ){
  int fail = 0;

  // Evaluate the points on the edges
  TMRPoint e[4];
  double params[4] = {u, v, 1.0-u, 1.0-v};

  for ( int k = 0; k < 4; k++ ){
    double p = (1.0 - params[k])*tmin[k] + params[k]*tmax[k];
    fail = fail || edges[k]->evalPoint(p, &e[k]);
  }

  // Evaluate the point on the surface
  X->x = (1.0-u)*e[3].x + u*e[1].x + (1.0-v)*e[0].x + v*e[2].x
    - ((1.0-u)*(1.0-v)*c[0].x + u*(1.0-v)*c[1].x +
       u*v*c[2].x + v*(1.0-u)*c[3].x);

  X->y = (1.0-u)*e[3].y + u*e[1].y + (1.0-v)*e[0].y + v*e[2].y
    - ((1.0-u)*(1.0-v)*c[0].y + u*(1.0-v)*c[1].y +
       u*v*c[2].y + v*(1.0-u)*c[3].y);

  X->z = (1.0-u)*e[3].z + u*e[1].z + (1.0-v)*e[0].z + v*e[2].z
    - ((1.0-u)*(1.0-v)*c[0].z + u*(1.0-v)*c[1].z +
       u*v*c[2].z + v*(1.0-u)*c[3].z);

  return fail;
}

/*
  Inverse evaluation: This is not yet implemented
*/
int TMRTFIFace::invEvalPoint( TMRPoint point,
                              double *uf, double *vf ){
  int fail = 0;
  *uf = 0.0;
  *vf = 0.0;

  // Get the bounds
  double umin, vmin, umax, vmax;
  getRange(&umin, &vmin, &umax, &vmax);

  // Set the initial guess
  double u = 0.5, v = 0.5;

  // Perform a newton iteration until convergence
  for ( int k = 0; k < max_newton_iters; k++ ){
    // Set the point and their derivatives
    TMRPoint X, Xu, Xv;
    fail = fail || evalDeriv(u, v, &X, &Xu, &Xv);

   // Compute the difference between the position on the surface
    // and the point
    TMRPoint r;
    r.x = (X.x - point.x);
    r.y = (X.y - point.y);
    r.z = (X.z - point.z);

    // Compute the residual
    double ru = Xu.dot(r);
    double rv = Xv.dot(r);

    // Compute the elements of the Jacobian matrix
    double Juu = Xu.dot(Xu);
    double Juv = Xu.dot(Xv);
    double Jvv = Xv.dot(Xv);

    double du = 0.0, dv = 0.0;

    // Check for the bounds on u
    if (u <= umin && ru >= 0.0){
      ru = 0.0;
      Juu = 1.0;
      Juv = 0.0;
    }
    else if (u >= umax && ru <= 0.0){
      ru = 0.0;
      Juu = 1.0;
      Juv = 0.0;
    }

    // Check for the bounds on v
    if (v <= vmin && rv >= 0.0){
      rv = 0.0;
      Jvv = 1.0;
      Juv = 0.0;
    }
    else if (v >= vmax && rv <= 0.0){
      rv = 0.0;
      Jvv = 1.0;
      Juv = 0.0;
    }

    // Solve the 2x2 system
    double det = Juu*Jvv - Juv*Juv;
    if (det != 0.0){
      du = (Jvv*ru - Juv*rv)/det;
      dv = (Juu*rv - Juv*ru)/det;
    }

    // Compute the updates
    double unew = u - du;
    double vnew = v - dv;

    // Truncate the u/v values to the bounds
    if (unew < umin){ unew = umin; }
    else if (unew > umax){ unew = umax; }
    if (vnew < vmin){ vnew = vmin; }
    else if (vnew > vmax){ vnew = vmax; }

    // Check if the convergence test is satisfied
    if (fabs(r.x) < eps_dist &&
        fabs(r.y) < eps_dist &&
        fabs(r.z) < eps_dist){
      *uf = u;
      *vf = v;
      return 0;
    }

    // Perform the cosine check
    double dotr = r.dot(r);
    double dotu = Xu.dot(Xu);
    double dotv = Xv.dot(Xv);
    if (ru*ru < eps_cosine*eps_cosine*dotu*dotr &&
        rv*rv < eps_cosine*eps_cosine*dotv*dotr){
      *uf = u;
      *vf = v;
      return 0;
    }

    // Update the new parameter values
    u = unew;
    v = vnew;
  }

  return fail;
}

/*
  Derivative evaluation: This is not yet implemented
*/
int TMRTFIFace::evalDeriv( double u, double v,
                           TMRPoint *X,
                           TMRPoint *Xu, TMRPoint *Xv ){
  int fail = 0;

  // Evaluate the points on the edges
  TMRPoint e[4], ed[4];
  double params[4] = {u, v, 1.0-u, 1.0-v};
  int d[4] = {1, 1, -1, -1};

  // Evaluate the curves along the v-direction
  for ( int k = 0; k < 4; k++ ){
    double p = (1.0 - params[k])*tmin[k] + params[k]*tmax[k];
    fail = fail || edges[k]->evalDeriv(p, &e[k], &ed[k]);

    // Evaluate the derivative w.r.t. u/v
    ed[k].x *= d[k]*(tmax[k] - tmin[k]);
    ed[k].y *= d[k]*(tmax[k] - tmin[k]);
    ed[k].z *= d[k]*(tmax[k] - tmin[k]);
  }

  // Evaluate the point on the surface
  X->x = (1.0-u)*e[3].x + u*e[1].x + (1.0-v)*e[0].x + v*e[2].x
    - ((1.0-u)*(1.0-v)*c[0].x + u*(1.0-v)*c[1].x +
       u*v*c[2].x + v*(1.0-u)*c[3].x);

  X->y = (1.0-u)*e[3].y + u*e[1].y + (1.0-v)*e[0].y + v*e[2].y
    - ((1.0-u)*(1.0-v)*c[0].y + u*(1.0-v)*c[1].y +
       u*v*c[2].y + v*(1.0-u)*c[3].y);

  X->z = (1.0-u)*e[3].z + u*e[1].z + (1.0-v)*e[0].z + v*e[2].z
    - ((1.0-u)*(1.0-v)*c[0].z + u*(1.0-v)*c[1].z +
       u*v*c[2].z + v*(1.0-u)*c[3].z);

  // Evaluate the point on the surface
  Xu->x = e[1].x - e[3].x + (1.0-v)*ed[0].x + v*ed[2].x
    - ((1.0-v)*(c[1].x - c[0].x) + v*(c[2].x - c[3].x));

  Xu->y = e[1].y - e[3].y + (1.0-v)*ed[0].y + v*ed[2].y
    - ((1.0-v)*(c[1].y - c[0].y) + v*(c[2].y - c[3].y));

  Xu->z = e[1].z - e[3].z + (1.0-v)*ed[0].z + v*ed[2].z
    - ((1.0-v)*(c[1].z - c[0].z) + v*(c[2].z - c[3].z));

  // Evaluate the point on the surface
  Xv->x = (1.0-u)*ed[3].x + u*ed[1].x + e[2].x - e[0].x +
    - ((1.0-u)*(c[3].x - c[0].x) + u*(c[2].x - c[1].x));

  Xv->y = (1.0-u)*ed[3].y + u*ed[1].y + e[2].y - e[0].y +
    - ((1.0-u)*(c[3].y - c[0].y) + u*(c[2].y - c[1].y));

  Xv->z = (1.0-u)*ed[3].z + u*ed[1].z + e[2].z - e[0].z +
    - ((1.0-u)*(c[3].z - c[0].z) + u*(c[2].z - c[1].z));

  return fail;
}

/*
  Create a parametric TFI

  The transfinite interpolation is performed in the parameter space
  and all points are obtained directly from the surface object itself.
*/
TMRParametricTFIFace::TMRParametricTFIFace( TMRFace *_face,
                                            TMREdge *_edges[],
                                            const int _dir[],
                                            TMRVertex *verts[] ){
  face = _face;
  face->incref();
  setName(face->getName());

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
    if (dir[k] > 0){
      edges[k]->getParamsOnFace(face, 0.0, dir[k],
                                &vupt[k], &vvpt[k]);
    }
    else {
      edges[k]->getParamsOnFace(face, 1.0, dir[k],
                                &vupt[k], &vvpt[k]);
    }
  }

  // Set the vertices for this surface
  if (dir[0] > 0){
    edges[0]->setVertices(verts[0], verts[1]);
  }
  else {
    edges[0]->setVertices(verts[1], verts[0]);
  }

  if (dir[1] > 0){
    edges[1]->setVertices(verts[1], verts[2]);
  }
  else {
    edges[1]->setVertices(verts[2], verts[1]);
  }

  if (dir[2] > 0){
    edges[2]->setVertices(verts[2], verts[3]);
  }
  else {
    edges[2]->setVertices(verts[3], verts[2]);
  }

  if (dir[3] > 0){
    edges[3]->setVertices(verts[3], verts[0]);
  }
  else {
    edges[3]->setVertices(verts[0], verts[3]);
  }

  // Set the curve segment into the curve
  TMREdgeLoop *loop = new TMREdgeLoop(4, edges, dir);
  addEdgeLoop(1, loop);
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
  int fail = 0;
  double cupt[4], cvpt[4];
  double params[4] = {u, v, 1.0-u, 1.0-v};

  for ( int k = 0; k < 4; k++ ){
    if (dir[k] > 0){
      fail = fail ||
        edges[k]->getParamsOnFace(face, params[k], dir[k],
                                  &cupt[k], &cvpt[k]);
    }
    else {
      fail = fail ||
        edges[k]->getParamsOnFace(face, 1.0-params[k], dir[k],
                                  &cupt[k], &cvpt[k]);
    }
  }

  // Compute the parametric coordinates
  double us, vs;
  us = (1.0-u)*cupt[3] + u*cupt[1] + (1.0-v)*cupt[0] + v*cupt[2]
    - ((1.0-u)*(1.0-v)*vupt[0] + u*(1.0-v)*vupt[1] +
       u*v*vupt[2] + v*(1.0-u)*vupt[3]);

  vs = (1.0-u)*cvpt[3] + u*cvpt[1] + (1.0-v)*cvpt[0] + v*cvpt[2]
    - ((1.0-u)*(1.0-v)*vvpt[0] + u*(1.0-v)*vvpt[1] +
       u*v*vvpt[2] + v*(1.0-u)*vvpt[3]);

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
                                     TMRPoint *X,
                                     TMRPoint *Xu, TMRPoint *Xv ){
  // Evaluate the curves along the v-direction
  int fail = 1;
  X->zero();
  Xu->zero();
  Xv->zero();

  return fail;
}

/*
  Create the TFI volume
*/
TMRTFIVolume::TMRTFIVolume( TMRFace *_faces[],
                            TMREdge *_edges[], const int _dir[],
                            TMRVertex *_verts[] ):
TMRVolume(6, _faces){
  for ( int i = 0; i < 6; i++ ){
    faces[i] = _faces[i];
    faces[i]->incref();
  }
  for ( int i = 0; i < 12; i++ ){
    edges[i] = _edges[i];
    edges[i]->incref();
    edge_dir[i] = _dir[i];
  }
  for ( int i = 0; i < 8; i++ ){
    verts[i] = _verts[i];
    verts[i]->incref();
    verts[i]->evalPoint(&c[i]);
  }
}

/*
  Destroy the TFI volume object
*/
TMRTFIVolume::~TMRTFIVolume(){
  for ( int i = 0; i < 6; i++ ){
    faces[i]->decref();
  }
  for ( int i = 0; i < 12; i++ ){
    edges[i]->decref();
  }
  for ( int i = 0; i < 8; i++ ){
    verts[i]->decref();
  }
}

/*
  Get the parameter range for this volume
*/
void TMRTFIVolume::getRange( double *umin, double *vmin, double *wmin,
                             double *umax, double *vmax, double *wmax ){
  *umin = *vmin = *wmin = 0.0;
  *umax = *vmax = *wmax = 1.0;
}

/*
  Evaluate the point within the volume
*/
int TMRTFIVolume::evalPoint( double u, double v, double w,
                             TMRPoint *X ){
  int fail = 0;

  // Evaluate/retrieve the points on the edges
  TMRPoint e[12];
  for ( int k = 0; k < 12; k++ ){
    double t;
    if (k < 4){
      t = u;
    }
    else if (k < 8){
      t = v;
    }
    else {
      t = w;
    }

    if (edge_dir[k] > 0){
      fail = fail || edges[k]->evalPoint(t, &e[k]);
    }
    else {
      fail = fail || edges[k]->evalPoint(1.0-t, &e[k]);
    }
  }

  X->x = ((1.0-v)*(1.0-w)*e[0].x + v*(1.0-w)*e[1].x +
       (1.0-v)*w*e[2].x + v*w*e[3].x +
       (1.0-u)*(1.0-w)*e[4].x + u*(1.0-w)*e[5].x +
       (1.0-u)*w*e[6].x + u*w*e[7].x +
       (1.0-u)*(1.0-v)*e[8].x + u*(1.0-v)*e[9].x +
       (1.0-u)*v*e[10].x + u*v*e[11].x)
    - 2.0*((1.0-u)*(1.0-v)*(1.0-w)*c[0].x + u*(1.0-v)*(1.0-w)*c[1].x +
           (1.0-u)*v*(1.0-w)*c[2].x + u*v*(1.0-w)*c[3].x +
           (1.0-u)*(1.0-v)*w*c[4].x + u*(1.0-v)*w*c[5].x +
           (1.0-u)*v*w*c[6].x + u*v*w*c[7].x);

  X->y = ((1.0-v)*(1.0-w)*e[0].y + v*(1.0-w)*e[1].y +
       (1.0-v)*w*e[2].y + v*w*e[3].y +
       (1.0-u)*(1.0-w)*e[4].y + u*(1.0-w)*e[5].y +
       (1.0-u)*w*e[6].y + u*w*e[7].y +
       (1.0-u)*(1.0-v)*e[8].y + u*(1.0-v)*e[9].y +
       (1.0-u)*v*e[10].y + u*v*e[11].y)
    - 2.0*((1.0-u)*(1.0-v)*(1.0-w)*c[0].y + u*(1.0-v)*(1.0-w)*c[1].y +
           (1.0-u)*v*(1.0-w)*c[2].y + u*v*(1.0-w)*c[3].y +
           (1.0-u)*(1.0-v)*w*c[4].y + u*(1.0-v)*w*c[5].y +
           (1.0-u)*v*w*c[6].y + u*v*w*c[7].y);

  X->z = ((1.0-v)*(1.0-w)*e[0].z + v*(1.0-w)*e[1].z +
       (1.0-v)*w*e[2].z + v*w*e[3].z +
       (1.0-u)*(1.0-w)*e[4].z + u*(1.0-w)*e[5].z +
       (1.0-u)*w*e[6].z + u*w*e[7].z +
       (1.0-u)*(1.0-v)*e[8].z + u*(1.0-v)*e[9].z +
       (1.0-u)*v*e[10].z + u*v*e[11].z)
    - 2.0*((1.0-u)*(1.0-v)*(1.0-w)*c[0].z + u*(1.0-v)*(1.0-w)*c[1].z +
           (1.0-u)*v*(1.0-w)*c[2].z + u*v*(1.0-w)*c[3].z +
           (1.0-u)*(1.0-v)*w*c[4].z + u*(1.0-v)*w*c[5].z +
           (1.0-u)*v*w*c[6].z + u*v*w*c[7].z);

  return 1;
}

/*
  Get the underlying face, edge and volume entities
*/
void TMRTFIVolume::getEntities( TMRFace ***_faces,
                                TMREdge ***_edges,
                                TMRVertex ***_verts ){
  if (_faces){ *_faces = faces; }
  if (_edges){ *_edges = edges; }
  if (_verts){ *_verts = verts; }
}
