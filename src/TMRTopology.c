#include "TMRTopology.h"
#include "TMRMesh.h"
#include <stdio.h>
#include <math.h>

/*
  Perform an inverse evaluation by obtaining the underlying 
  parametrization on the specified curve 
*/
int TMRVertex::getParamOnEdge( TMREdge *edge, double *t ){
  TMRPoint p;
  int fail = evalPoint(&p);
  fail = fail || edge->invEvalPoint(p, t);
  return fail;
}

/*
  Same thing, except on the specified surface
*/
int TMRVertex::getParamsOnFace( TMRFace *face,
                                double *u, double *v ){
  TMRPoint p;
  int fail = evalPoint(&p);
  fail = fail || face->invEvalPoint(p, u, v);
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
  Set the original vertices
*/
TMREdge::TMREdge(){
  v1 = v2 = NULL;
  mesh = NULL;
}

/*
  Free the edge
*/
TMREdge::~TMREdge(){
  if (v1){ v1->decref(); }
  if (v2){ v2->decref(); }
}

/*
  Perform the inverse evaluation
*/
int TMREdge::invEvalPoint( TMRPoint p, double *t ){
  *t = 0.0;
  int fail = 1;
  return fail;
}

/*
  Set the step size for the derivative
*/
double TMREdge::deriv_step_size = 1e-6;

/*
  Evaluate the derivative using a finite-difference step size
*/
int TMREdge::evalDeriv( double t, TMRPoint *Xt ){
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
  Find the point on the surface closest to the point C(t)
*/
int TMREdge::getParamsOnFace( TMRFace *face, double t, 
                              int dir, double *u, double *v ){
  TMRPoint p;
  int fail = evalPoint(t, &p);
  fail = fail || face->invEvalPoint(p, u, v);
  return fail;
}

/*
  Set the adjacent vertices
*/
void TMREdge::setVertices( TMRVertex *_v1, TMRVertex *_v2 ){
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
void TMREdge::getVertices( TMRVertex **_v1, TMRVertex **_v2 ){
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
void integrateEdge( TMREdge *edge,
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
double TMREdge::integrate( double t1, double t2, double tol,
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
void TMREdge::setMesh( TMREdgeMesh *_mesh ){
  _mesh->incref();
  if (mesh){ mesh->decref(); }
  mesh = _mesh;
}

/*
  Retrieve the mesh pointer
*/
void TMREdge::getMesh( TMREdgeMesh **_mesh ){
  *_mesh = mesh;
}

/*
  Create the edge loop object
*/
TMREdgeLoop::TMREdgeLoop( int _nedges, TMREdge *_edges[], 
                          const int _dir[] ){
  int fail = 0;
  if (_nedges == 0){
    fail = 1;
    fprintf(stderr, "TMREdgeLoop error: Zero length segment\n");
  }

  // First, check whether the loop is closed
  TMRVertex *vinit = NULL;
  TMRVertex *vnext;
  for ( int i = 0; i < _nedges; i++ ){
    TMRVertex *v1, *v2;
    if (_dir[i] > 0){
      _edges[i]->getVertices(&v1, &v2);
    }
    else {
      _edges[i]->getVertices(&v2, &v1);
    }
    if (i == 0){ vinit = v1; }
    vnext = v2;
    if (i == _nedges-1){
      if (vinit != vnext){
        fprintf(stderr, 
                "TMREdgeLoop error: Edge loop must be closed\n");
        fail = 1;
      }
    }
  }

  // Return if the segment is not closed
  if (fail){
    nedges = 0;
    edges = NULL;
    dir = NULL;
    return;
  }

  nedges = _nedges;
  edges = new TMREdge*[ nedges ];
  dir = new int[ nedges ];
  for ( int k = 0; k < nedges; k++ ){
    edges[k] = _edges[k];
    edges[k]->incref();
    dir[k] = _dir[k];
  }
}

/*
  Free the data from the edge loop object
*/
TMREdgeLoop::~TMREdgeLoop(){
  for ( int k = 0; k < nedges; k++ ){
    edges[k]->decref();
  }
  delete [] edges;
  delete [] dir;
}

/*
  Retrieve the edge loop edges
*/
void TMREdgeLoop::getEdgeLoop( int *_nedges, TMREdge **_edges[], 
                               const int *_dir[] ){
  if (_nedges){ *_nedges = nedges; }
  if (_edges){ *_edges = edges; }
  if (_dir){ *_dir = dir; }
}

/*
  Initialize data within the TMRSurface object
*/
TMRFace::TMRFace(){
  max_num_loops = 0;
  num_loops = 0;
  loops = NULL;
  mesh = NULL;
}

/*
  Deallocate the curve segments (if any)
*/
TMRFace::~TMRFace(){
  if (loops){
    for ( int i = 0; i < num_loops; i++ ){
      loops[i]->decref();
    }
    delete [] loops;
  }
  if (mesh){
    mesh->decref();
  }
}

/*
  Set the step size for the derivative
*/
double TMRFace::deriv_step_size = 1e-6;

/*
  Evaluate the derivative using a finite-difference step size
*/
int TMRFace::evalDeriv( double u, double v, 
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
  Perform the inverse evaluation
*/
int TMRFace::invEvalPoint( TMRPoint p, double *u, double *v ){
  *u = *v = 0.0;
  int fail = 1;
  return fail;
}

/*
  Add the curves that bound the surface
*/
void TMRFace::addEdgeLoop( TMREdgeLoop *loop ){
  int fail = 0;
  
  // Increase the reference count
  loop->incref();

  if (num_loops >= max_num_loops){
    // Extend the loops array
    max_num_loops = (2*num_loops > 10 ? 2*num_loops : 10);

    // Allocate the new loops array
    TMREdgeLoop **lps = new TMREdgeLoop*[ max_num_loops ];

    // Copy over any existing segments
    if (num_loops > 0){
      memcpy(lps, loops, num_loops*sizeof(TMREdgeLoop*));
      delete [] loops;
    }
    loops = lps;
  }

  // Set the new segment array
  loops[num_loops] = loop;
  num_loops++;
}

/*
  Get the number of closed segments
*/
int TMRFace::getNumEdgeLoops(){
  return num_loops;
}

/*
  Retrieve the information from the given segment number
*/
void TMRFace::getEdgeLoop( int k, TMREdgeLoop **_loop ){
  *_loop = NULL;
  if (k >= 0 && k < num_loops){
    if (_loop){ *_loop = loops[k]; }
  }
}

/*
  Set the mesh into the array
*/
void TMRFace::setMesh( TMRFaceMesh *_mesh ){
  _mesh->incref();
  if (mesh){ mesh->decref(); }
  mesh = _mesh;
}

/*
  Retrieve the mesh pointer
*/
void TMRFace::getMesh( TMRFaceMesh **_mesh ){
  *_mesh = mesh;
}

/*
  The TMRModel class containing all of the required geometry
  objects.
*/
TMRModel::TMRModel( int _num_vertices, TMRVertex **_vertices, 
                    int _num_edges, TMREdge **_edges,
                    int _num_faces, TMRFace **_faces ){
  num_vertices = _num_vertices;
  num_edges = _num_edges;
  num_faces = _num_faces;

  vertices = new TMRVertex*[ num_vertices ];
  edges = new TMREdge*[ num_edges ];
  faces = new TMRFace*[ num_faces ];

  for ( int i = 0; i < num_vertices; i++ ){
    vertices[i] = _vertices[i];
    vertices[i]->incref();
  }

  for ( int i = 0; i < num_edges; i++ ){
    edges[i] = _edges[i];
    edges[i]->incref();
  }

  for ( int i = 0; i < num_faces; i++ ){
    faces[i] = _faces[i];
    faces[i]->incref();
  }

  ordered_verts = new OrderedPair<TMRVertex>[ num_vertices ];
  ordered_edges = new OrderedPair<TMREdge>[ num_edges ];
  ordered_faces = new OrderedPair<TMRFace>[ num_faces ];

  for ( int i = 0; i < num_vertices; i++ ){
    ordered_verts[i].num = i;
    ordered_verts[i].obj = vertices[i];
  }

  for ( int i = 0; i < num_edges; i++ ){
    ordered_edges[i].num = i;
    ordered_edges[i].obj = edges[i];
  }

  for ( int i = 0; i < num_faces; i++ ){
    ordered_faces[i].num = i;
    ordered_faces[i].obj = faces[i];
  }

  // Sort the vertices, curves and surfaces
  qsort(ordered_verts, num_vertices, sizeof(OrderedPair<TMRVertex>),
        compare_ordered_pairs<TMRVertex>);
  qsort(ordered_edges, num_edges, sizeof(OrderedPair<TMREdge>),
        compare_ordered_pairs<TMREdge>);
  qsort(ordered_faces, num_faces, sizeof(OrderedPair<TMRFace>),
        compare_ordered_pairs<TMRFace>);

  verify();
}

/*
  Free the geometry objects
*/
TMRModel::~TMRModel(){
  for ( int i = 0; i < num_vertices; i++ ){
    vertices[i]->decref();
  }
  for ( int i = 0; i < num_edges; i++ ){
    edges[i]->decref();
  }
  for ( int i = 0; i < num_faces; i++ ){
    faces[i]->decref();
  }
  delete [] vertices;
  delete [] edges;
  delete [] faces;
  delete [] ordered_verts;
  delete [] ordered_edges;
  delete [] ordered_faces;
}

/*
  Verify that all objects that are referenced in the geometry have a
  verifiable list name and that objects are referred to one or more
  times (but not zero times!)
*/
int TMRModel::verify(){
  int fail = 0;
  int *verts = new int[ num_vertices ];
  int *crvs = new int[ num_edges ];
  memset(verts, 0, num_vertices*sizeof(int));
  memset(crvs, 0, num_edges*sizeof(int));

  for ( int face = 0; face < num_faces; face++ ){
    int nloops = faces[face]->getNumEdgeLoops();
    
    for ( int k = 0; k < nloops; k++ ){
      TMREdgeLoop *loop;
      faces[face]->getEdgeLoop(k, &loop);

      int ncurves;
      TMREdge **curves;
      loop->getEdgeLoop(&ncurves, &curves, NULL);

      // Loop over all of the curves and check whether the data exists
      // or not
      for ( int j = 0; j < ncurves; j++ ){
        int cindex = getEdgeIndex(curves[j]);
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
  for ( int i = 0; i < num_edges; i++ ){
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
void TMRModel::getVertices( int *_num_vertices, 
                            TMRVertex ***_vertices ){
  if (_num_vertices){ *_num_vertices = num_vertices; }
  if (_vertices){ *_vertices = vertices; }
}

/*
  Retrieve the curves
*/
void TMRModel::getEdges( int *_num_edges, 
                         TMREdge ***_edges ){
  if (_num_edges){ *_num_edges = num_edges; }
  if (_edges){ *_edges = edges; }
}

/*
  Retrieve the surfaces
*/
void TMRModel::getFaces( int *_num_faces, 
                         TMRFace ***_faces ){
  if (_num_faces){ *_num_faces = num_faces; }
  if (_faces){ *_faces = faces; }
}

/*
  Static member function for sorting the ordered pairs
*/
template <class ctype>
int TMRModel::compare_ordered_pairs( const void *avoid, const void *bvoid ){
  const OrderedPair<ctype> *a = static_cast<const OrderedPair<ctype>*>(avoid);
  const OrderedPair<ctype> *b = static_cast<const OrderedPair<ctype>*>(bvoid);
  return a->obj - b->obj; // Use pointer arithmetic to determine relative positions
}

/*
  Retrieve the index given the vertex point
*/
int TMRModel::getVertexIndex( TMRVertex *vertex ){
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
int TMRModel::getEdgeIndex( TMREdge *edge ){
  OrderedPair<TMREdge> pair;
  pair.num = -1;
  pair.obj = edge;

  // Search for the ordered pair
  OrderedPair<TMREdge> *item = 
    (OrderedPair<TMREdge>*)bsearch(&pair, ordered_edges, num_edges,
                                   sizeof(OrderedPair<TMREdge>),
                                   compare_ordered_pairs<TMREdge>);
  if (item){
    return item->num;
  }
  return -1;
}

/*
  Retrieve the index given the pointer to the surface object
*/
int TMRModel::getFaceIndex( TMRFace *face ){
  OrderedPair<TMRFace> pair;
  pair.num = -1;
  pair.obj = face;

  // Search for the ordered pair
  OrderedPair<TMRFace> *item = 
    (OrderedPair<TMRFace>*)bsearch(&pair, ordered_faces, num_faces,
                                   sizeof(OrderedPair<TMRFace>),
                                   compare_ordered_pairs<TMRFace>);
  if (item){
    return item->num;
  }
  return -1;
}

/*
  The main topology class that contains the objects used to build the
  underlying mesh.  
*/
TMRTopology::TMRTopology( TMRModel *_geo ){
  // Increase the ref. count to the geometry object
  geo = _geo;
  geo->incref();

  // Get the geometry objects
  int num_vertices, num_edges, num_faces;
  TMRVertex **vertices;
  TMREdge **edges;
  TMRFace **faces;
  geo->getVertices(&num_vertices, &vertices);
  geo->getEdges(&num_edges, &edges);
  geo->getFaces(&num_faces, &faces);

  // Build the topology information based on the inputs
  face_to_vertices = new int[ 4*num_faces ];
  face_to_edges = new int[ 4*num_faces ];
  edge_to_vertices = new int[ 2*num_edges ];

  // Sort the addresses to make the array
  const int coordinate_to_edge[] = {3, 1, 0, 2};

  // Create the face -> edge information
  for ( int i = 0; i < num_faces; i++ ){
    TMREdgeLoop *loop;
    faces[i]->getEdgeLoop(0, &loop);
    TMREdge **e;
    loop->getEdgeLoop(NULL, &e, NULL);

    // Search for the face indices
    for ( int jp = 0; jp < 4; jp++ ){
      int j = coordinate_to_edge[jp];
      face_to_edges[4*i+jp] = geo->getEdgeIndex(e[j]);
    }
  }

  // Create the edge -> vertex information
  for ( int i = 0; i < num_edges; i++ ){
    TMRVertex *v1, *v2;
    edges[i]->getVertices(&v1, &v2);
    
    edge_to_vertices[2*i] = geo->getVertexIndex(v1);
    edge_to_vertices[2*i+1] = geo->getVertexIndex(v2);
  }

  // Create the face -> vertex information. Within the TMR
  // geometry routines, the ordering is as shown on the left.
  // The quadrant/octree code uses the coordinate ordering shown
  // on the right.
  // 
  //  From TMRSurface         Coordinate-ordering
  //  v3---e2---v2            v2---e3---v3
  //  |         |             |         |
  //  e3        e1            e0        e1
  //  |         |             |         |
  //  v0---e0---v1            v0---e2---v1
  //  C.C.W. direction        Coordinate direction
  
  for ( int i = 0; i < num_faces; i++ ){
    TMREdgeLoop *loop;
    faces[i]->getEdgeLoop(0, &loop);

    const int *edge_dir;
    loop->getEdgeLoop(NULL, NULL, &edge_dir);

    // Coordinate-ordered edge e0
    int edge = face_to_edges[4*i];
    if (edge_dir[3] > 0){
      face_to_vertices[4*i] = edge_to_vertices[2*edge+1];
      face_to_vertices[4*i+2] = edge_to_vertices[2*edge];
    }
    else {
      face_to_vertices[4*i] = edge_to_vertices[2*edge];
      face_to_vertices[4*i+2] = edge_to_vertices[2*edge+1];
    }

    // Coordinate-ordered edge e1
    edge = face_to_edges[4*i+1];
    if (edge_dir[1] > 0){
      face_to_vertices[4*i+1] = edge_to_vertices[2*edge];
      face_to_vertices[4*i+3] = edge_to_vertices[2*edge+1];
    }
    else {
      face_to_vertices[4*i+1] = edge_to_vertices[2*edge+1];
      face_to_vertices[4*i+3] = edge_to_vertices[2*edge];
    }
  }
}

/*
  Free the topology data
*/
TMRTopology::~TMRTopology(){
  // Free the geometry object
  geo->decref();
  delete [] face_to_vertices;
  delete [] face_to_edges;
  delete [] edge_to_vertices;
}

/*
  Retrieve the face object
*/
void TMRTopology::getSurface( int face_num, TMRFace **face ){
  int num_faces;
  TMRFace **faces;
  geo->getFaces(&num_faces, &faces);
  *face = NULL;
  if (faces && (face_num >= 0 && face_num < num_faces)){
    *face = faces[face_num];
  }
}

/*
  Retrieve the curve object associated with the given face/edge index
*/
void TMRTopology::getFaceCurve( int face_num, int edge_index, 
                                TMREdge **edge ){
  int num_faces;
  TMRFace **faces;
  geo->getFaces(&num_faces, &faces);

  // Coordinate edge number to actual edge number
  const int coordinate_to_edge[] = {3, 1, 0, 2};

  *edge = NULL;
  if (faces && (face_num >= 0 && face_num < num_faces)){
    TMREdgeLoop *loop;
    faces[face_num]->getEdgeLoop(0, &loop);

    // Get the edge loop
    TMREdge **edgs;
    loop->getEdgeLoop(NULL, &edgs, NULL);
    if (edge_index >= 0 && edge_index < 4){
      edge_index = coordinate_to_edge[edge_index];
      *edge = edgs[edge_index];
    }
  }
}

/*
  Retrieve the connectivity information
*/
void TMRTopology::getConnectivity( int *nnodes, 
                                   int *nedges, int *nfaces,
                                   const int **face_nodes,
                                   const int **face_edges ){
  int num_vertices, num_edges, num_faces;
  geo->getVertices(&num_vertices, NULL);
  geo->getEdges(&num_edges, NULL);
  geo->getFaces(&num_faces, NULL);
  *nnodes = num_vertices;
  *nedges = num_edges;
  *nfaces = num_faces;
  *face_nodes = face_to_vertices;
  *face_edges = face_to_edges;
}
