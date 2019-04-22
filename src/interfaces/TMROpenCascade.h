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

#ifndef TMR_OPENCASCADE_H
#define TMR_OPENCASCADE_H

#ifdef TMR_HAS_OPENCASCADE

/*
  Include the TMR files required
*/
#include "TMRGeometry.h"
#include "TMRTopology.h"

/*
  Include all of the files required from OpenCASCADE
*/
#include <Standard.hxx>
#include <Standard_Type.hxx>
#include <Standard_Integer.hxx>
#include <Standard_Real.hxx>
#include <Standard_ErrorHandler.hxx>
#include <Standard_Version.hxx>
#include <Precision.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Compound.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools.hxx>

#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <GeomLib.hxx>
#include <GeomProjLib.hxx>
#include <Geom2dAdaptor_HCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_HSurface.hxx>
#include <Geom2dAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Adaptor3d_HSurface.hxx>
#include <Adaptor2d_HCurve2d.hxx>
#include <Adaptor3d_CurveOnSurface.hxx>

#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepLib.hxx>
#include <BRepBndLib.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepLProp_CLProps.hxx>
#include <BRepBuilderAPI.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_FindPlane.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepFeat_SplitShape.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <BRep_Builder.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeUpgrade_ShapeDivideClosed.hxx>
#include <ShapeFix.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeConstruct.hxx>
#include <ShapeConstruct_Curve.hxx>

#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>

/*
  Topoloy objects for OpenCascade:

  TMR_OCCVertex
  TMR_OCCEdge
  TMR_OCCFace
*/
class TMR_OCCVertex : public TMRVertex {
 public:
  TMR_OCCVertex( TopoDS_Vertex &v );
  ~TMR_OCCVertex();
  int evalPoint( TMRPoint *p );
  int getParamOnEdge( TMREdge *edge, double *t );
  int getParamsOnFace( TMRFace *face, double *u, double *v );
  int isSame( TMRVertex *v );
  void getVertexObject( TopoDS_Vertex &v );
 private:
  TopoDS_Vertex vert;
};

class TMR_OCCEdge : public TMREdge {
 public:
  TMR_OCCEdge( TopoDS_Edge &e );
  ~TMR_OCCEdge();
  void getRange( double *tmin, double *tmax );
  int evalPoint( double t, TMRPoint *X );
  int invEvalPoint( TMRPoint X, double *t );
  int evalDeriv( double t, TMRPoint *X, TMRPoint *Xt );
  int eval2ndDeriv( double t, TMRPoint *X, TMRPoint *Xt, TMRPoint *Xtt );
  int getParamsOnFace( TMRFace *face, double t,
                       int dir, double *u, double *v );
  int isSame( TMREdge *e );
  int isDegenerate();
  void getEdgeObject( TopoDS_Edge &e );
 private:
  TopoDS_Edge edge;
  TopoDS_Edge reverse_edge;
};

class TMR_OCCFace : public TMRFace {
 public:
  TMR_OCCFace( int _normal_dir, TopoDS_Face &f );
  ~TMR_OCCFace();
  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax );
  int evalPoint( double u, double v, TMRPoint *X );
  int invEvalPoint( TMRPoint p, double *u, double *v );
  int evalDeriv( double u, double v, TMRPoint *X,
                 TMRPoint *Xu, TMRPoint *Xv );
  int eval2ndDeriv( double u, double v,
                    TMRPoint *X,
                    TMRPoint *Xu, TMRPoint *Xv,
                    TMRPoint *Xuu, TMRPoint *Xuv, TMRPoint *Xvv );
  int isSame( TMRFace *f );
  void getFaceObject( TopoDS_Face &f );
 private:
  TopoDS_Face face;
};

/*
  Initialization of the OpenCascade geometry from an IGES/STEP files
*/
TMRModel* TMR_LoadModelFromIGESFile( const char *filename,
                                     int print_level=0 );
TMRModel* TMR_LoadModelFromSTEPFile( const char *filename,
                                     int print_level=0 );
TMRModel* TMR_LoadModelFromCompound( TopoDS_Compound &compound,
                                     int print_level=0 );

#endif // TMR_HAS_OPENCASCADE
#endif // TMR_OPENCASCADE_H
