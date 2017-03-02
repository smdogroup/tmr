#include "TMROpenCascade.h"

/*
  TMR interface to OpenCascade for the curve
*/
TMR_OCCCurve::TMR_OCCCurve( Geom_Curve &c ){
  curve = c;
}

TMR_OCCCurve::~TMR_OCCCurve(){}

void TMR_OCCCurve::getRange( double *tmin, double *tmax ){
  *tmin = curve.FirstParameter();
  *tmax = curve.LastParameter();
} 
 
int TMR_OCCCurve::evalPoint( double t, TMRPoint *X ){
  gp_Pnt p;
  curve.D0(t, p);
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
  curve.D1(t, p, v1);
  Xt->x = v1.X();
  Xt->y = v1.Y();
  Xt->z = v1.Z();
  return 0;
}

/*
  TMR interface to the OpenCascade geometry
*/
TMR_OCCSurface::TMR_OCCSurface( Geom_Surface &s ){
  surf = s;
}

TMR_OCCSurface::~TMR_OCCSurface(){}

void TMR_OCCSurface::getRange( double *umin, double *vmin,
                               double *umax, double *vmax ){
  *umin = surf.FirstUParameter();
  *vmin = surf.FirstVParameter();
  *umax = surf.LastUParameter();
  *vmax = surf.LastVParameter();
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
  surf.D1(u, v, p, pu, pv);
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
  return TMRVertex::getParamsOnCurve(edge, t);
}

int TMR_OCCVertex::getParamsOnFace( TMRFace *face,
                                    double *u, double *v ){
  TMR_OCCFace *f = dynamic_cast<TMR_OCCFace*>(face);  
  if (f){
    TopoDS_Face occ_face;
    f->getSurfaceObject(occ_face);
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
}

TMR_OCCEdge::~TMR_OCCEdge(){}

void TMR_OCCEdge::getRange( double *tmin, double *tmax ){
  BRepAdaptor_Curve curve(edge);
  *tmin = curve.FirstParameter();
  *tmax = curve.LastParameter();
}

int TMR_OCCEdge::getParamsOnSurface( TMRFace *face, double t, 
                                     int dir, double *u, double *v ){
  TMR_OCCFace *f = dynamic_cast<TMR_OCCFace*>(face);
  if (f){
    // Get the face topology
    TopoDS_Face occ_face;
    f->getSurfaceObject(occ_face);

    // Get the min/max range
    double *tmin, *tmax;
    getRange(&tmin, &tmax);

    // Get the pcurve on the surface
    Handle(Geom2d_Curve) pcurve;
    if (dir > 0){
      pcurve = BRep_Tool::CurveOnSurface(edge, occ_face, tmin, tmax);
    }
    else{
      pcurve = BRep_Tool::CurveOnSurface(edge.Reverse(), occ_face, tmin, tmax);
    }

    if (pcurve.IsNull()){
      return 1;
    }

    // Evaluate the point on the curve
    gp_Pnt2d p;
    pcurve.D0(t, p);
    *u = p.X();
    *v = p.Y();
    return 0;
  }

  // Fall through to the underlying implementation
  return TMREdge::getParamsOnFace(face, t, dir, u, v);
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
  BRepAdaptor_Curve curve(edge);
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
  BRepAdaptor_Surface surf(face);
  *umin = surf.FirstUParameter();
  *vmin = surf.FirstVParameter();
  *umax = surf.LastUParameter();
  *vmax = surf.LastVParameter();
}

int TMR_OCCFace::evalPoint( double u, double v, TMRPoint *X ){
  BRepAdaptor_Surface surf(face);
  gp_Pnt p;
  surf.D0(u, v, p);
  X->x = p.X();
  X->y = p.Y();
  X->z = p.Z();
  return 0;
}

int TMR_OCCFace::invEvalPoint( TMRPoint p, double *u, double *v ){
  gp_Pnt pt(X.x, X.y, X.z);
  BRepAdaptor_Surface surf(face);
  GeomAPI_ProjectPointOnSurf projection(pt, surf);
  if (projection.NbPoints() == 0){
    return 1;
  } 
  else {
    projection.LowerDistanceParameters(*u, *v);
    return 0;
  }
}

int TMR_OCCFace::evalDeriv( double u, double v, 
                            TMRPoint *Xu, TMRPoint *Xv ){
  BRepAdaptor_Surface surf(face);
  gp_Pnt p;
  gp_Vec pu, pv;
  surf.D1(u, v, p, pu, pv);
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

TMRGeometry* TMR_LoadGeometryFromSTEPFile( const char *filename ){
  FILE *fp = fopen(filename, "r");
  if (!fp){
    int fail = 1;
    return fail;
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
      printf("TMR Warning: Transfer %d not OK!\n", k+1);
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
  
  // Create the index <--> geometry object
  int nverts = 0, nedges = 0, nfaces = 0, nwires = 0;
  TopTools_IndexedMapOfShape verts, edges, faces, wires;

  // Get the faces, edges and  up the number of faces and a
  TopExp_Explorer Exp;
  Exp.Init(compound, TopAbs_VERTEX);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    verts.Add(shape);
    nverts++;
  }

  Exp.Init(compound, TopAbs_EDGE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    edges.Add(shape);
    nedges++;
  }

  Exp.Init(compound, TopAbs_FACE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    faces.Add(shape);
    nfaces++;
  }

  Exp.Init(compound, TopAbs_WIRE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    wires.Add(shape);
    nwires++;
  }

  // Re-iterate through the list and create the objects needed to
  // define the geometry in TMR
  TMRVertex **all_vertices = new TMRVertex*[ nverts ];
  Exp.Init(compound, TopAbs_VERTEX);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Vertex v = TopoDS::Vertex(shape);
    int index = verts.FindIndex(shape);
    all_vertices[index-1] = new TMR_OCCVertex(v);
  }

  TMREdge **all_edges = new TMREdge*[ nedges ];
  Exp.Init(compound, TopAbs_EDGE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Edge e = TopoDS::Edge(shape);
    int index = edges.FindIndex(shape);
    all_edges[index-1] = new TMR_OCCEdge(e);

    // Set the vertices belonging to this edge
    TopoDS_Vertex v1 = TopExp::FirstVertex(e);
    TopoDS_Vertex v2 = TopExp::LastVertex(e);
    int idx1 = verts.FindIndex(v1);
    int idx2 = verts.FindIndex(v2);
    all_edges[index-1]->setVertices(all_vertices[idx1-1], all_vertices[idx2-1]);
  }

  TMRFace **all_faces = new TMRFace*[ nfaces ];
  Exp.Init(compound, TopAbs_FACE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Edge f = TopoDS::Face(shape);
    int index = faces.FindIndex(shape);
    all_faces[index-1] = new TMR_OCCFace(f);
  }

  // Loop over all of the wire objects in the model
  TMREdgeLoop **all_loops = new TMREdgeLoop[ nwires ];
  Exp.Init(compound, TopAbs_WIRE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Wire wire = TopoDS::Wire(shape);
    int index = wires.FindIndex(shape);

    // Check if the wire is closed
    TopExp::Vertices(wire, v1, v2);
    int closed = 0;
    if (!v1.IsNull() && !v2.IsNull()){
      if (v1.IsSame(v2)){
        closed = 1;
      }
    }

    // Count up th enumber of wires
    BRepTools_WireExplorer wExp;
    int ne = 0;
    for ( wExp.Init(wire); wExp.More(); wExp.Next()){
      ne++;
    }

    // Create the edge list
    TMREdge **edgs = new TMREdge[ ne ];
    int *dir = new int[ ne ];
    wExp.Init(wire);
    for ( int k = 0; wExp.More; wExp.Next(), k++ ){
      TopoDS_Shape shape = wExp.Current()
      dir[k] = 1;
      if (shape.Orientation() == TopAbs_REVERSED){
        dir[k] = -1;
      } 
      int index = edges.FindIndex(shape);
      edgs[k] = all_edges[index-1];
    }

    // Allocate the loop with the given edges/directions
    all_loops[index-1] = new TMREdgeLoop(ne, edgs, dir);

    // Free the allocated data
    delete [] edgs;
    delete [] dir;
  }

  Exp.Init(compound, TopAbs_FACE);
  for ( ; Exp.More(); Exp.Next() ){
    TopoDS_Shape shape = Exp.Current();    
    TopoDS_Shape face = TopoDS::Face(shape);
    int index = faces.FindIndex(shape);

    // Get the wires attached to this face
    TopTools_IndexedMapOfShape map;
    TopExp::MapShapes(face, TopAbs_WIRE, map);
    for ( int i = 1; i <= map.Extent(); i++ ){
      TopoDS_Wire wire_shape = map(i);
      int loop_index = wires.FindIndex(wire_shape);
      all_faces[index-1]->addEdgeLoop(all_loops[loop_index-1]);
    }
  }
}
