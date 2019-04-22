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

#include "TMROpenCascade.h"

#ifdef TMR_HAS_OPENCASCADE

/*
  The TMR interface to the underlying OpenCascade vertex object
*/
TMR_OCCVertex::TMR_OCCVertex( TopoDS_Vertex &v ){
  vert = v;
}

TMR_OCCVertex::~TMR_OCCVertex(){}

int TMR_OCCVertex::evalPoint( TMRPoint *X ){
  gp_Pnt p = BRep_Tool::Pnt(vert);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();
  return 0;
}

int TMR_OCCVertex::getParamOnEdge( TMREdge *edge, double *t ){
  TMR_OCCEdge *e = dynamic_cast<TMR_OCCEdge*>(edge);
  if (e){
    TopoDS_Edge occ_edge;
    e->getEdgeObject(occ_edge);
    *t = BRep_Tool::Parameter(vert, occ_edge);
    return 0;
  }
  return TMRVertex::getParamOnEdge(edge, t);
}

int TMR_OCCVertex::getParamsOnFace( TMRFace *face,
                                    double *u, double *v ){
  TMR_OCCFace *f = dynamic_cast<TMR_OCCFace*>(face);
  if (f){
    TopoDS_Face occ_face;
    f->getFaceObject(occ_face);
    gp_Pnt2d pt = BRep_Tool::Parameters(vert, occ_face);
    *u = pt.X();
    *v = pt.Y();
    return 0;
  }
  return TMRVertex::getParamsOnFace(face, u, v);
}

int TMR_OCCVertex::isSame( TMRVertex *vt ){
  TMR_OCCVertex *v = dynamic_cast<TMR_OCCVertex*>(vt);
  if (v){
    return vert.IsSame(v->vert);
  }
  return 0;
}

void TMR_OCCVertex::getVertexObject( TopoDS_Vertex &v ){
  v = vert;
}

/*
  TMR interface to the underlying OpenCascade edge object
*/
TMR_OCCEdge::TMR_OCCEdge( TopoDS_Edge &e ){
  edge = e;
  reverse_edge = e;
  reverse_edge.Reverse();
}

TMR_OCCEdge::~TMR_OCCEdge(){}

void TMR_OCCEdge::getRange( double *tmin, double *tmax ){
  BRep_Tool::Range(edge, *tmin, *tmax);
}

int TMR_OCCEdge::getParamsOnFace( TMRFace *surface, double t,
                                  int dir, double *u, double *v ){
  TMR_OCCFace *f = dynamic_cast<TMR_OCCFace*>(surface);
  if (f){
    // Get the face topology
    TopoDS_Face occ_face;
    f->getFaceObject(occ_face);

    // Get the pcurve on the surface
    Handle(Geom2d_Curve) pcurve;
    if (dir > 0){
      double t0, t1;
      pcurve = BRep_Tool::CurveOnSurface(edge, occ_face, t0, t1);
    }
    else{
      double t0, t1;
      pcurve = BRep_Tool::CurveOnSurface(reverse_edge, occ_face, t0, t1);
    }

    if (pcurve.IsNull()){
      return 1;
    }

    // Evaluate the point on the curve
    gp_Pnt2d p;
    pcurve->D0(t, p);
    *u = p.X();
    *v = p.Y();
    return 0;
  }

  // Fall through to the underlying implementation
  return TMREdge::getParamsOnFace(surface, t, dir, u, v);
}

int TMR_OCCEdge::evalPoint( double t, TMRPoint *X ){
  BRepAdaptor_Curve curve(edge);
  gp_Pnt p;
  curve.D0(t, p);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();
  return 0;
}

int TMR_OCCEdge::invEvalPoint( TMRPoint X, double *t ){
  gp_Pnt pt(X.x, X.y, X.z);
  double tmin, tmax;
  getRange(&tmin, &tmax);
  const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, tmin, tmax);
  GeomAPI_ProjectPointOnCurve projection(pt, curve);
  if (projection.NbPoints() == 0){
    return 1;
  }
  else {
    *t = projection.LowerDistanceParameter();
    return 0;
  }
}

int TMR_OCCEdge::evalDeriv( double t, TMRPoint *X, TMRPoint *Xt ){
  int fail = 0;
  gp_Pnt p;
  gp_Vec pt;
  BRepAdaptor_Curve curve(edge);
  curve.D1(t, p, pt);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();

  Xt->x = pt.X();
  Xt->y = pt.Y();
  Xt->z = pt.Z();
  return fail;
}

int TMR_OCCEdge::eval2ndDeriv( double t, TMRPoint *X,
                               TMRPoint *Xt, TMRPoint *Xtt ){
  int fail = 0;
  gp_Pnt p;
  gp_Vec pt, ptt;
  BRepAdaptor_Curve curve(edge);
  curve.D2(t, p, pt, ptt);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();

  Xt->x = pt.X();
  Xt->y = pt.Y();
  Xt->z = pt.Z();

  Xtt->x = ptt.X();
  Xtt->y = ptt.Y();
  Xtt->z = ptt.Z();
  return fail;
}

int TMR_OCCEdge::isSame( TMREdge *et ){
  TMR_OCCEdge *e = dynamic_cast<TMR_OCCEdge*>(et);
  if (e){
    return edge.IsSame(e->edge);
  }
  return 0;
}

int TMR_OCCEdge::isDegenerate(){
  return BRep_Tool::Degenerated(edge);
}

void TMR_OCCEdge::getEdgeObject( TopoDS_Edge &e ){
  e = edge;
}

/*
  TMR interface to the underlying OpenCascade TopoDS_Face object
*/
TMR_OCCFace::TMR_OCCFace( int _normal_dir, TopoDS_Face &f ):
TMRFace(_normal_dir){
  face = f;
}

TMR_OCCFace::~TMR_OCCFace(){}

void TMR_OCCFace::getRange( double *umin, double *vmin,
                            double *umax, double *vmax ){
  BRepTools::UVBounds(face, *umin, *umax, *vmin, *vmax);
}

int TMR_OCCFace::evalPoint( double u, double v, TMRPoint *X ){
  const Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
  gp_Pnt p;
  surf->D0(u, v, p);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();
  return 0;
}

int TMR_OCCFace::invEvalPoint( TMRPoint X, double *u, double *v ){
  gp_Pnt pt(X.x, X.y, X.z);
  const Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
  GeomAPI_ProjectPointOnSurf projection(pt, surf);
  if (projection.NbPoints() == 0){
    *u = *v = 0.0;
    return 1;
  }
  else {
    projection.LowerDistanceParameters(*u, *v);
    return 0;
  }
}

int TMR_OCCFace::evalDeriv( double u, double v,
                            TMRPoint *X,
                            TMRPoint *Xu, TMRPoint *Xv ){
  const Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
  gp_Pnt p;
  gp_Vec pu, pv;
  surf->D1(u, v, p, pu, pv);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();

  Xu->x = pu.X();
  Xu->y = pu.Y();
  Xu->z = pu.Z();

  Xv->x = pv.X();
  Xv->y = pv.Y();
  Xv->z = pv.Z();
  return 0;
}

int TMR_OCCFace::eval2ndDeriv( double u, double v,
                               TMRPoint *X,
                               TMRPoint *Xu,
                               TMRPoint *Xv,
                               TMRPoint *Xuu,
                               TMRPoint *Xuv,
                               TMRPoint *Xvv ){
  const Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
  gp_Pnt p;
  gp_Vec pu, pv;
  gp_Vec puu, puv, pvv;
  surf->D2(u, v, p, pu, pv, puu, pvv, puv);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();

  Xu->x = pu.X();
  Xu->y = pu.Y();
  Xu->z = pu.Z();

  Xv->x = pv.X();
  Xv->y = pv.Y();
  Xv->z = pv.Z();

  Xuu->x = puu.X();
  Xuu->y = puu.Y();
  Xuu->z = puu.Z();

  Xvv->x = pvv.X();
  Xvv->y = pvv.Y();
  Xvv->z = pvv.Z();

  Xuv->x = puv.X();
  Xuv->y = puv.Y();
  Xuv->z = puv.Z();
  return 0;
}

int TMR_OCCFace::isSame( TMRFace *ft ){
  TMR_OCCFace *f = dynamic_cast<TMR_OCCFace*>(ft);
  if (f){
    return face.IsSame(f->face);
  }
  return 0;
}

void TMR_OCCFace::getFaceObject( TopoDS_Face &f ){
  f = face;
}

/*
  Remove from a compound, those objects that are not referenced by the
  given TopoAbs type
*/
void TMR_RemoveFloatingShapes( TopoDS_Compound &compound,
                               TopAbs_ShapeEnum shape_type ){
  // Keep track of what we have already added to the new compound
  TopTools_IndexedMapOfShape shapes;

  TopExp_Explorer solidExp(compound, TopAbs_SOLID);
  for ( ; solidExp.More(); solidExp.Next() ){
    TopoDS_Shape solid = solidExp.Current();
    if (!shapes.Contains(solid)){
      shapes.Add(solid);

      TopExp_Explorer shellExp(solid, TopAbs_SHELL);
      for ( ; shellExp.More(); shellExp.Next() ){
        TopoDS_Shape shell = shellExp.Current();
        if (!shapes.Contains(shell)){
          shapes.Add(shell);

          TopExp_Explorer faceExp(shell, TopAbs_FACE);
          for ( ; faceExp.More(); faceExp.Next() ){
            TopoDS_Shape face = faceExp.Current();
            if (!shapes.Contains(face)){
              shapes.Add(face);

              TopExp_Explorer wireExp(face, TopAbs_WIRE);
              for ( ; wireExp.More(); wireExp.Next() ){
                TopoDS_Shape wire = wireExp.Current();
                if (!shapes.Contains(wire)){
                  shapes.Add(wire);

                  TopExp_Explorer edgeExp(wire, TopAbs_EDGE);
                  for ( ; edgeExp.More(); edgeExp.Next() ){
                    TopoDS_Shape edge = edgeExp.Current();
                    if (!shapes.Contains(edge)){
                      shapes.Add(edge);

                      TopExp_Explorer vertexExp(edge, TopAbs_VERTEX);
                      for ( ; vertexExp.More(); vertexExp.Next() ){
                        TopoDS_Shape vertex = vertexExp.Current();
                        if (!shapes.Contains(vertex)){
                          shapes.Add(vertex);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if (shape_type != TopAbs_SHELL){
    TopExp_Explorer shellExp(compound, TopAbs_SHELL);
    for ( ; shellExp.More(); shellExp.Next() ){
      TopoDS_Shape shell = shellExp.Current();
      if (!shapes.Contains(shell)){
        shapes.Add(shell);

        TopExp_Explorer faceExp(shell, TopAbs_FACE);
        for ( ; faceExp.More(); faceExp.Next() ){
          TopoDS_Shape face = faceExp.Current();
          if (!shapes.Contains(face)){
            shapes.Add(face);

            TopExp_Explorer wireExp(face, TopAbs_WIRE);
            for ( ; wireExp.More(); wireExp.Next() ){
              TopoDS_Shape wire = wireExp.Current();
              if (!shapes.Contains(wire)){
                shapes.Add(wire);

                TopExp_Explorer edgeExp(wire, TopAbs_EDGE);
                for ( ; edgeExp.More(); edgeExp.Next() ){
                  TopoDS_Shape edge = edgeExp.Current();
                  if (!shapes.Contains(edge)){
                    shapes.Add(edge);

                    TopExp_Explorer vertexExp(edge, TopAbs_VERTEX);
                    for ( ; vertexExp.More(); vertexExp.Next() ){
                      TopoDS_Shape vertex = vertexExp.Current();
                      if (!shapes.Contains(vertex)){
                        shapes.Add(vertex);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (shape_type != TopAbs_FACE){
      TopExp_Explorer faceExp(compound, TopAbs_FACE);
      for ( ; faceExp.More(); faceExp.Next() ){
        TopoDS_Shape face = faceExp.Current();
        if (!shapes.Contains(face)){
          shapes.Add(face);

          TopExp_Explorer wireExp(face, TopAbs_WIRE);
          for ( ; wireExp.More(); wireExp.Next() ){
            TopoDS_Shape wire = wireExp.Current();
            if (!shapes.Contains(wire)){
              shapes.Add(wire);

              TopExp_Explorer edgeExp(wire, TopAbs_EDGE);
              for ( ; edgeExp.More(); edgeExp.Next() ){
                TopoDS_Shape edge = edgeExp.Current();
                if (!shapes.Contains(edge)){
                  shapes.Add(edge);

                  TopExp_Explorer vertexExp(edge, TopAbs_VERTEX);
                  for ( ; vertexExp.More(); vertexExp.Next() ){
                    TopoDS_Shape vertex = vertexExp.Current();
                    if (!shapes.Contains(vertex)){
                      shapes.Add(vertex);
                    }
                  }
                }
              }
            }
          }
        }
      }

      if (shape_type != TopAbs_EDGE){
        TopExp_Explorer wireExp(compound, TopAbs_WIRE);
        for ( ; wireExp.More(); wireExp.Next() ){
          TopoDS_Shape wire = wireExp.Current();
          if (!shapes.Contains(wire)){
            shapes.Add(wire);

            TopExp_Explorer edgeExp(wire, TopAbs_EDGE);
            for ( ; edgeExp.More(); edgeExp.Next() ){
              TopoDS_Shape edge = edgeExp.Current();
              if (!shapes.Contains(edge)){
                shapes.Add(edge);

                TopExp_Explorer vertexExp(edge, TopAbs_VERTEX);
                for ( ; vertexExp.More(); vertexExp.Next() ){
                  TopoDS_Shape vertex = vertexExp.Current();
                  if (!shapes.Contains(vertex)){
                    shapes.Add(vertex);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Build the new compound
  TopoDS_Compound new_compound;
  BRep_Builder builder;
  builder.MakeCompound(new_compound);
  for ( int i = 1; i <= shapes.Extent(); i++ ){
    TopoDS_Shape shape = shapes(i);
    builder.Add(new_compound, shape);
  }

  compound = new_compound;
}

/*
  Create the TMRModel based on the IGES input file
*/
TMRModel* TMR_LoadModelFromIGESFile( const char *filename,
                                     int print_level ){
  FILE *fp = fopen(filename, "r");
  if (!fp){
    return NULL;
  }
  fclose(fp);

  IGESControl_Reader reader;
  IFSelect_ReturnStatus status = reader.ReadFile(filename);
  if (status != IFSelect_RetDone){
    fprintf(stderr, "TMR Warning: IGES file reader failed\n");
  }

  // inspect the root transfers
  reader.PrintCheckLoad(Standard_False, IFSelect_ItemsByEntity);

  // The number of root objects for transfer
  int nroots = reader.NbRootsForTransfer();
  for ( int k = 1; k <= nroots; k++ ){
    Standard_Boolean ok = reader.TransferOneRoot(k);
    if (!ok){
      fprintf(stderr, "TMR Warning: Transfer %d not OK!\n", k);
    }
  }

  // Load the different shapes
  int nbs = reader.NbShapes();

  // Build the shape
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  for ( int i = 1; i <= nbs; i++ ){
    TopoDS_Shape shape = reader.Shape(i);
    builder.Add(compound, shape);
  }

  // TMR_RemoveFloatingShapes(compound, TopAbs_FACE);

  return TMR_LoadModelFromCompound(compound, print_level);
}

/*
  Create the TMRModel based on the STEP input file
*/
TMRModel* TMR_LoadModelFromSTEPFile( const char *filename,
                                     int print_level ){
  FILE *fp = fopen(filename, "r");
  if (!fp){
    return NULL;
  }
  fclose(fp);

  STEPControl_Reader reader;
  IFSelect_ReturnStatus status = reader.ReadFile(filename);
  if (status != IFSelect_RetDone){
    fprintf(stderr, "TMR Warning: STEP file reader failed\n");
  }

  // inspect the root transfers
  reader.PrintCheckLoad(Standard_False, IFSelect_ItemsByEntity);

  // The number of root objects for transfer
  int nroots = reader.NbRootsForTransfer();
  for ( int k = 1; k <= nroots; k++ ){
    Standard_Boolean ok = reader.TransferRoot(k);
    if (!ok){
      fprintf(stderr, "TMR Warning: Transfer %d not OK!\n", k);
    }
  }

  // Load the different shapes
  int nbs = reader.NbShapes();

  // Build the shape
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  for ( int i = 1; i <= nbs; i++ ){
    TopoDS_Shape shape = reader.Shape(i);
    builder.Add(compound, shape);
  }

  TMR_RemoveFloatingShapes(compound, TopAbs_EDGE);

  return TMR_LoadModelFromCompound(compound, print_level);
}

/*
  Create the TMRModel based on the TopoDS_Compound object
*/
TMRModel* TMR_LoadModelFromCompound( TopoDS_Compound &compound,
                                     int print_level ){
  // Create the index <--> geometry object
  TopTools_IndexedMapOfShape verts, edges, faces, wires, shells, solids;

  // Extract the vertices, edges, faces, wires etc. from the
  // compound
  TopExp_Explorer Exp;
  Exp.Init(compound, TopAbs_VERTEX);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    verts.Add(shape);
  }

  Exp.Init(compound, TopAbs_EDGE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    edges.Add(shape);
  }

  Exp.Init(compound, TopAbs_WIRE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    wires.Add(shape);
  }

  Exp.Init(compound, TopAbs_FACE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    faces.Add(shape);
  }

  Exp.Init(compound, TopAbs_SHELL);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    shells.Add(shape);
  }

  Exp.Init(compound, TopAbs_SOLID);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    solids.Add(shape);
  }

  int nverts = verts.Extent();
  int nedges = edges.Extent();
  int nwires = wires.Extent();
  int nfaces = faces.Extent();
  int nshells = shells.Extent();
  int nsolids = solids.Extent();

  if (print_level > 0){
    printf("Compound loaded with:\nnverts = %d nedges = %d nfaces = %d "
           "nwires = %d nshells = %d nsolids = %d\n",
           nverts, nedges, nfaces, nwires, nshells, nsolids);
  }

  // Re-iterate through the list and create the objects needed to
  // define the geometry in TMR
  TMRVertex **all_vertices = new TMRVertex*[ nverts ];
  memset(all_vertices, 0, nverts*sizeof(TMRVertex*));
  for ( int index = 1; index <= verts.Extent(); index++ ){
    TopoDS_Vertex v = TopoDS::Vertex(verts(index));
    all_vertices[index-1] = new TMR_OCCVertex(v);
  }

  TMREdge **all_edges = new TMREdge*[ nedges ];
  memset(all_edges, 0, nedges*sizeof(TMREdge*));
  for ( int index = 1; index <= edges.Extent(); index++ ){
    TopoDS_Edge e = TopoDS::Edge(edges(index).Oriented(TopAbs_FORWARD));
    all_edges[index-1] = new TMR_OCCEdge(e);

    // Set the vertices belonging to this edge
    TopoDS_Vertex v1, v2;

    // Note that this call does not consider orientation
    TopExp::Vertices(e, v1, v2);
    int idx1 = verts.FindIndex(v1);
    int idx2 = verts.FindIndex(v2);
    all_edges[index-1]->setVertices(all_vertices[idx1-1],
                                    all_vertices[idx2-1]);
  }

  TMRFace **all_faces = new TMRFace*[ nfaces ];
  memset(all_faces, 0, nfaces*sizeof(TMRFace*));
  for ( int index = 1; index <= faces.Extent(); index++ ){
    TopoDS_Face face = TopoDS::Face(faces(index));

    // Check if the orientation of the face is flipped relative to the
    // orientation of the surface
    int orient = 1;
    if (face.Orientation() == TopAbs_REVERSED){
      orient = -1;
    }
    all_faces[index-1] = new TMR_OCCFace(orient, face);

    // Find the wires connected to the face
    TopTools_IndexedMapOfShape wire_map;
    TopExp::MapShapes(face, TopAbs_WIRE, wire_map);
    for ( int i = 1; i <= wire_map.Extent(); i++ ){
      TopoDS_Wire wire = TopoDS::Wire(wire_map(i));

      // Count up the enumber of edges
      BRepTools_WireExplorer wExp;
      int ne = 0;
      for ( wExp.Init(wire); wExp.More(); wExp.Next()){
        ne++;
      }

      // Create the edge list
      TMREdge **edgs = new TMREdge*[ ne ];
      int *dir = new int[ ne ];
      int k = 0;

      for ( wExp.Init(wire); wExp.More(); wExp.Next(), k++ ){
        TopoDS_Edge edge = wExp.Current();
        dir[k] = 1;
        if (edge.Orientation() == TopAbs_REVERSED){
          dir[k] = -1;
        }
        int edge_index = edges.FindIndex(edge);
        edgs[k] = all_edges[edge_index-1];
      }

      // This is not useful in this context
      int loop_orient = -1;
      all_faces[index-1]->addEdgeLoop(loop_orient,
                                      new TMREdgeLoop(ne, edgs, dir));

      // Free the allocated data
      delete [] edgs;
      delete [] dir;
    }
  }

  // Create the volumes
  TMRVolume **all_vols = new TMRVolume*[ nsolids ];
  memset(all_vols, 0, nsolids*sizeof(TMRVolume*));
  for ( int index = 1; index <= solids.Extent(); index++ ){
    TopoDS_Solid solid = TopoDS::Solid(solids(index));
    int vol_index = solids.FindIndex(solid);

    // Count up the number of faces
    int nvol_faces = 0;
    TopTools_IndexedMapOfShape shell_map;
    TopExp::MapShapes(solid, TopAbs_SHELL, shell_map);
    for ( int i = 1; i <= shell_map.Extent(); i++ ){
      TopoDS_Shell shell = TopoDS::Shell(shell_map(i));

      // Find the number of faces associated with this shell
      TopTools_IndexedMapOfShape face_map;
      TopExp::MapShapes(shell, TopAbs_FACE, face_map);
      nvol_faces += face_map.Extent();
    }

    // Allocate the faces arrays/directions
    TMRFace **vol_faces = new TMRFace*[ nvol_faces ];

    // Now, extract the faces from the underlying shell(s)
    nvol_faces = 0;
    TopExp::MapShapes(solid, TopAbs_SHELL, shell_map);
    for ( int i = 1; i <= shell_map.Extent(); i++ ){
      TopoDS_Shell shell = TopoDS::Shell(shell_map(i));

      TopTools_IndexedMapOfShape face_map;
      TopExp::MapShapes(shell, TopAbs_FACE, face_map);
      for ( int j = 1; j <= face_map.Extent(); j++ ){
        TopoDS_Face face = TopoDS::Face(face_map(j));

        // Find the index of the face and set its direction relative
        // to the face stored in the OCC file. If the direction is
        // reversed, store an orientation of -1.
        int index = faces.FindIndex(face);

        // Assign the face pointer
        vol_faces[nvol_faces] = all_faces[index-1];
        nvol_faces++;
      }
    }

    // Create the volume object
    all_vols[vol_index-1] = new TMRVolume(nvol_faces, vol_faces);
    delete [] vol_faces;
  }

  TMRModel *geo = new TMRModel(nverts, all_vertices,
                               nedges, all_edges,
                               nfaces, all_faces,
                               nsolids, all_vols);

  // Free the arrays
  delete [] all_vertices;
  delete [] all_edges;
  delete [] all_faces;
  delete [] all_vols;

  return geo;
}

#endif // TMR_HAS_OPENCASCADE
