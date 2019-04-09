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

#ifndef TMR_EGADS_H
#define TMR_EGADS_H

#ifdef TMR_HAS_EGADS

/*
  Include the TMR files required
*/
#include "TMRGeometry.h"
#include "TMRTopology.h"

namespace TMR_EgadsInterface {

#include "egads.h"

class TMR_EgadsContext : public TMREntity {
  public:
  TMR_EgadsContext();
  ~TMR_EgadsContext();
  ego getContext();
 private:
  ego ctx;
};

/*
  Topoloy objects for OpenCascade:

  TMR_EgadsNode
  TMR_EgadsEdge
  TMR_EgadsFace
*/
class TMR_EgadsNode : public TMRVertex {
 public:
  TMR_EgadsNode( TMR_EgadsContext* _ctx, ego _node );
  ~TMR_EgadsNode();
  int evalPoint( TMRPoint *p );
  int getParamOnEdge( TMREdge *edge, double *t );
  int getParamsOnFace( TMRFace *face, double *u, double *v );
  void getNodeObject( ego *n );
 private:
  TMR_EgadsContext *ctx;
  ego node;
};

class TMR_EgadsEdge : public TMREdge {
 public:
  TMR_EgadsEdge( TMR_EgadsContext* _ctx, ego _edge, int _is_degenerate=0 );
  ~TMR_EgadsEdge();
  void getRange( double *tmin, double *tmax );
  int evalPoint( double t, TMRPoint *X );
  int invEvalPoint( TMRPoint X, double *t );
  int evalDeriv( double t, TMRPoint *Xt );
  int eval2ndDeriv( double t, TMRPoint *Xtt );
  int getParamsOnFace( TMRFace *face, double t,
                       int dir, double *u, double *v );
  int isDegenerate();
  void getEdgeObject( ego *e );
 private:
  TMR_EgadsContext *ctx;
  ego edge;
  int is_degenerate;
};

class TMR_EgadsFace : public TMRFace {
 public:
  TMR_EgadsFace( TMR_EgadsContext* _ctx, int _normal_dir, ego _face );
  ~TMR_EgadsFace();
  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax );
  int evalPoint( double u, double v, TMRPoint *X );
  int invEvalPoint( TMRPoint p, double *u, double *v );
  int evalDeriv( double u, double v, TMRPoint *Xu, TMRPoint *Xv );
  int eval2ndDeriv( double u, double v,
                    TMRPoint *Xuu, TMRPoint *Xuv, TMRPoint *Xvv );
  void getFaceObject( ego *f );
 private:
  TMR_EgadsContext *ctx;
  ego face;
};

}

/*
  Initialization of the OpenCascade geometry from an IGES/STEP files
*/
TMRModel* TMR_LoadModelFromEGADSFile( const char *filename,
                                      int print_level=0 );

#endif // TMR_HAS_EGADS
#endif // TMR_EGADS_H
