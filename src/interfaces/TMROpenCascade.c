#include "TMROpenCascade.h"

#ifdef TMR_HAS_OPENCASCADE

/*
  TMR interface to OpenCascade for the curve
*/
TMR_OCCCurve::TMR_OCCCurve( Handle(Geom_Curve) &c ){
  curve = c;
}

TMR_OCCCurve::~TMR_OCCCurve(){}

void TMR_OCCCurve::getRange( double *tmin, double *tmax ){
  *tmin = curve->FirstParameter();
  *tmax = curve->LastParameter();
} 
 
int TMR_OCCCurve::evalPoint( double t, TMRPoint *X ){
  gp_Pnt p;
  curve->D0(t, p);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();
  return 0;
}

int TMR_OCCCurve::invEvalPoint( TMRPoint X, double *t ){
  gp_Pnt pt(X.x, X.y, X.z);
  GeomAPI_ProjectPointOnCurve projection(pt, curve);
  if (projection.NbPoints() == 0){
    return 1;
  } 
  else {
    *t = projection.LowerDistanceParameter();
    return 0;
  }
}

int TMR_OCCCurve::evalDeriv( double t, TMRPoint *Xt ){
  gp_Pnt p;
  gp_Vec v1;
  curve->D1(t, p, v1);
  Xt->x = v1.X();
  Xt->y = v1.Y();
  Xt->z = v1.Z();
  return 0;
}

/*
  TMR interface to the OpenCascade geometry
*/
TMR_OCCSurface::TMR_OCCSurface( Handle(Geom_Surface) &s ){
  surf = s;
}

TMR_OCCSurface::~TMR_OCCSurface(){}

void TMR_OCCSurface::getRange( double *umin, double *vmin,
                               double *umax, double *vmax ){
  surf->Bounds(*umin, *umax, *vmin, *vmax);
}
 
int TMR_OCCSurface::evalPoint( double u, double v, TMRPoint *X ){
  gp_Pnt p;
  surf->D0(u, v, p);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();
  return 0;
}
  
int TMR_OCCSurface::invEvalPoint( TMRPoint X, double *u, double *v ){
  gp_Pnt pt(X.x, X.y, X.z);
  GeomAPI_ProjectPointOnSurf projection(pt, surf);
  if (projection.NbPoints() == 0){
    return 1;
  } 
  else {
    projection.LowerDistanceParameters(*u, *v);
    return 0;
  }
}

int TMR_OCCSurface::evalDeriv( double u, double v, 
                               TMRPoint *Xu, TMRPoint *Xv ){
  gp_Pnt p;
  gp_Vec pu, pv;
  surf->D1(u, v, p, pu, pv);
  Xu->x = pu.X();
  Xu->y = pu.Y();
  Xu->z = pu.Z();
  Xv->x = pv.X();
  Xv->y = pv.Y();
  Xv->z = pv.Z();
  return 0;
}

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

int TMR_OCCEdge::evalDeriv( double t, TMRPoint *Xt ){
  gp_Pnt p;
  gp_Vec pt;
  BRepAdaptor_Curve curve(edge);
  curve.D1(t, p, pt);
  Xt->x = pt.X();
  Xt->y = pt.Y();
  Xt->z = pt.Z();
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
TMR_OCCFace::TMR_OCCFace( TopoDS_Face &f ){
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
                            TMRPoint *Xu, TMRPoint *Xv ){
  const Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
  gp_Pnt p;
  gp_Vec pu, pv;
  surf->D1(u, v, p, pu, pv);
  Xu->x = pu.X();
  Xu->y = pu.Y();
  Xu->z = pu.Z();
  Xv->x = pv.X();
  Xv->y = pv.Y();
  Xv->z = pv.Z();
  return 0;
}

void TMR_OCCFace::getFaceObject( TopoDS_Face &f ){
  f = face;
}

/*
  Create the TMRModel based on the STEP input file
*/
TMRModel* TMR_LoadModelFromSTEPFile( const char *filename ){
  FILE *fp = fopen(filename, "r");
  if (!fp){
    return NULL;
  }
  fclose(fp);
  
  STEPControl_Reader reader;
  IFSelect_ReturnStatus status = reader.ReadFile(filename);

  // inspect the root transfers
  reader.PrintCheckLoad(Standard_False, IFSelect_ItemsByEntity);

  // The number of root objects for transfer
  int nroots = reader.NbRootsForTransfer();
  for ( int k = 0; k < nroots; k++ ){
    Standard_Boolean ok = reader.TransferRoot(k+1);
    if (!ok){
      fprintf(stderr, "TMR Warning: Transfer %d not OK!\n", k+1);
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

  return TMR_LoadModelFromCompound(compound);
}

/*
  Create the TMRModel based on the TopoDS_Compound object
*/
TMRModel* TMR_LoadModelFromCompound( TopoDS_Compound &compound ){
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
  printf("Compound loaded with:\nnverts = %d nedges = %d \
nfaces = %d nwires = %d nshells = %d nsolids = %d\n",
         nverts, nedges, nfaces, nwires, nshells, nsolids);

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
    all_faces[index-1] = new TMR_OCCFace(face);

    // Get the face with the normal orientation
    TopoDS_Face face_norm = TopoDS::Face(faces(index).Oriented(TopAbs_FORWARD));

    TopTools_IndexedMapOfShape map;
    TopExp::MapShapes(face_norm, TopAbs_WIRE, map);
    for ( int i = 1; i <= map.Extent(); i++ ){
      TopoDS_Wire wire = TopoDS::Wire(map(i));

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
      for ( wExp.Init(wire, face_norm); wExp.More(); wExp.Next(), k++ ){
        TopoDS_Edge edge = wExp.Current();
        dir[k] = 1;
        if (edge.Orientation() == TopAbs_REVERSED){
          dir[k] *= -1;
        }
        int edge_index = edges.FindIndex(edge);
        edgs[k] = all_edges[edge_index-1];
      }

      if (face.IsEqual(face_norm)){
        // Allocate the loop with the given edges/directions
        all_faces[index-1]->addEdgeLoop(new TMREdgeLoop(ne, edgs, dir));
      }
      else {
        // The face is reversed
        for ( int j = 0; j < ne/2; j++ ){
          // swap edges 
          TMREdge *etmp = edgs[j];
          edgs[j] = edgs[ne-1-j];
          edgs[ne-1-j] = etmp;
          
          // swap the edge directions
          int tmp = dir[j];
          dir[j] = dir[ne-1-j];
          dir[ne-1-j] = tmp;
        }

        // Multiply all the directions by -1
        for ( int j = 0; j < ne; j++ ){
          dir[j] *= -1;
        }

        // Allocate the reverse loop
        all_faces[index-1]->addEdgeLoop(new TMREdgeLoop(ne, edgs, dir));
      }

      // Free the allocated data
      delete [] edgs;
      delete [] dir;
    }
  }

  // Create the volumes
  TMRVolume **all_vols = new TMRVolume*[ nsolids ];
  memset(all_vols, 0, nsolids*sizeof(TMRVolume*));
  Exp.Init(compound, TopAbs_SOLID);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Solid solid = TopoDS::Solid(shape);
    int vol_index = solids.FindIndex(shape);

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
    int *dir = new int[ nvol_faces ];

    // Now, extract the faces from the underlying shell(s)
    nvol_faces = 0;
    TopExp::MapShapes(solid, TopAbs_SHELL, shell_map);
    for ( int i = 1; i <= shell_map.Extent(); i++ ){
      TopoDS_Shell shell = TopoDS::Shell(shell_map(i));

      TopTools_IndexedMapOfShape face_map;
      TopExp::MapShapes(shell, TopAbs_FACE, face_map);
      for ( int j = 1; j <= face_map.Extent(); j++ ){
        TopoDS_Face face = TopoDS::Face(face_map(j));

        dir[nvol_faces] = 1;
        int index = faces.FindIndex(face);
        if (faces(index).IsEqual(face)){
          dir[nvol_faces] = 1;
        }
        else {
          dir[nvol_faces] = -1;
        }

        // Assign the face pointer
        vol_faces[nvol_faces] = all_faces[index-1];
        nvol_faces++;
      }
    }

    // Create the volume object
    all_vols[vol_index-1] = new TMRVolume(nvol_faces, vol_faces, dir);
    delete [] vol_faces;
    delete [] dir;
  }

  TMRModel *geo = new TMRModel(nverts, all_vertices,
                               nedges, all_edges,
                               nfaces, all_faces,
                               nsolids, all_vols);
  return geo;
}

#endif // TMR_HAS_OPENCASCADE
