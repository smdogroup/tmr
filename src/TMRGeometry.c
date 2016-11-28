#include "TMRGeometry.h"
#include <math.h>

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
void integrateEdge( TMRGeoEdge *edge,
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
  double int1 = (tmid - t1)*pointDist(&p1, &pmid);
  double int2 = (t2 - tmid)*pointDist(&pmid, &p2);
  double int3 = (tmid - t1)*pointDist(&p1, &p2);

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
double TMRGeoEdge::integrate( double t1, double t2, double tol,
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
  Compare two edge entries
*/
int compare_edge_entries( const void *aobj, const void *bobj ){
  const EdgeEntry *a = static_cast<const EdgeEntry*>(aobj);
  const EdgeEntry *b = static_cast<const EdgeEntry*>(bobj);

  return a->t - b->t;
}

/*
  Compare two face entries based on their parametric location on the
  surface. (Note that the parameteric locations are integeters, enabling
  exact comparisons. 
*/
int compare_face_entries( const void *aobj, const void *bobj ){
  const FaceEntry *a = static_cast<const FaceEntry*>(aobj);
  const FaceEntry *b = static_cast<const FaceEntry*>(bobj);
  
  // Take an exclusive or between the x/y coordinates
  uint32_t xxor = a->u ^ b->u;
  uint32_t yxor = a->v ^ b->v;

  // This is now a mask for the most significant bit
  uint32_t sor = xxor | yxor;

  // Check for the most-significant bit
  int discrim = 0;
  if (xxor > (sor ^ xxor)){
    discrim = a->u - b->u;
  }
  else {
    discrim = a->v - b->v;
  }

  if (discrim > 0){
    return 1;
  }
  else if (discrim < 0){
    return -1;
  }
  return 0;
}

/*
  Create the edge look up table
*/
TMR_EdgeLookup::TMR_EdgeLookup( int32_t *t, TMRPoint *pts, int _npts ){
  // Create the points
  npts = _npts;
  points = new EdgeEntry[ npts ];

  // Copy over the values
  for ( int k = 0; k < npts; k++ ){
    points[k].t = t[k];
    points[k].x = pts[k].x;
    points[k].y = pts[k].y;
    points[k].z = pts[k].z;
  }

  // Sort the points by edge index
  qsort(points, npts, sizeof(EdgeEntry), compare_edge_entries);
}

/*
  Free the edge entries
*/
TMR_EdgeLookup::~TMR_EdgeLookup(){
  delete [] points;
}

/*
  Edge look up based on an edge integer parameter
*/
int TMR_EdgeLookup::evalPoint( int32_t t_int, TMRPoint *pt ){
  // Set the fail flag
  int fail = 0;

  // Search the list of edge points to find this one
  EdgeEntry dummy;
  dummy.t = t_int;

  // Search for this entry
  void *ptr = bsearch(&dummy, points, npts, sizeof(EdgeEntry),
                      compare_edge_entries);
  
  // Copy over the entries
  if (ptr){
    EdgeEntry *entry = static_cast<EdgeEntry*>(ptr);
    if (entry){
      // copy over the values
      pt->x = entry->x;
      pt->y = entry->y;
      pt->z = entry->z;
    }
    else {
      fail = 1;
    }
  }
  else {
    fail = 1;
  }

  return fail;
}

/*
  Create the edge look up table
*/
TMR_SurfaceLookup::TMR_SurfaceLookup( int32_t *u, int32_t *v, 
                                      TMRPoint *pts, int _npts ){
  // Create the points
  npts = _npts;
  points = new FaceEntry[ npts ];

  // Copy over the values
  for ( int k = 0; k < npts; k++ ){
    points[k].u = u[k];
    points[k].v = v[k];
    points[k].x = pts[k].x;
    points[k].y = pts[k].y;
    points[k].z = pts[k].z;
  }

  // Sort the points by surface index
  qsort(points, npts, sizeof(FaceEntry), compare_face_entries);
}

/*
  Free the edge entries
*/
TMR_SurfaceLookup::~TMR_SurfaceLookup(){
  delete [] points;
}

/*
  Face look up based on u/v surface integer parameters.
  
  Evaluate the point by checking for values from the table stored in
  this surface object.
*/
int TMR_SurfaceLookup::evalPoint( int32_t u_int, int32_t v_int,
                                  TMRPoint *pt ){
  // Set the fail flag
  int fail = 0;

  // Search the list of edge points to find this one
  FaceEntry dummy;
  dummy.u = u_int;
  dummy.v = v_int;

  // Search for this entry
  void *ptr = bsearch(&dummy, points, npts, sizeof(FaceEntry),
                      compare_face_entries);
  
  // Copy over the entries
  if (ptr){
    FaceEntry *entry = static_cast<FaceEntry*>(ptr);
    if (entry){
      // copy over the values
      pt->x = entry->x;
      pt->y = entry->y;
      pt->z = entry->z;
    }
    else {
      fail = 1;
    }
  }
  else {
    fail = 1;
  }

  return fail;
}

/*
  Set up the transfinite-interpolation edge class.

  This class just interpolates between the two end points that are
  provided. 
*/
TMR_TFIEdge::TMR_TFIEdge( TMRPoint *_a, TMRPoint *_b ){
  // Copy over the point locations
  a = *_a;
  b = *_b;
}

/*
  There is no allocated data that needs to be freed
*/
TMR_TFIEdge::~TMR_TFIEdge(){}

/*
  Trans-finite interpolation 
*/
int TMR_TFIEdge::evalPoint( int32_t t_int, TMRPoint *pt ){
  // Convert the input to a parameter
  const double t = 1.0*t_int/(1 << TMR_MAX_LEVEL);

  // Interpolate between the two edges
  pt->x = (1.0 - t)*a.x + t*b.x;
  pt->y = (1.0 - t)*a.y + t*b.y;
  pt->z = (1.0 - t)*a.z + t*b.z;

  return 0;
}

/*
  Store the data for the transfinite interpolation
*/
TMR_TFISurface::TMR_TFISurface( TMRPoint *corners, 
                                TMR_EdgeLookup *_edges[] ){
  for ( int k = 0; k < 4; k++ ){
    c[k] = corners[k];
    edges[k] = _edges[k];
  }  
}

/*
  Transfinite interpolation along a face
*/
int TMR_TFISurface::evalPoint( int32_t u_int, int32_t v_int,
                               TMRPoint *pt ){
  // Set the fail flag to false
  int fail = 0;

  // Evaluate the points along the edges
  TMRPoint e0, e1, e2, e3;
  fail = fail || edges[0]->evalPoint(v_int, &e0);
  fail = fail || edges[1]->evalPoint(v_int, &e1);
  fail = fail || edges[2]->evalPoint(u_int, &e2);
  fail = fail || edges[3]->evalPoint(u_int, &e3);

  // Convert the integer coordinates into real coordinates
  const double u = 1.0*u_int/(1 << TMR_MAX_LEVEL);
  const double v = 1.0*v_int/(1 << TMR_MAX_LEVEL);

  // Perform the transfinite interpolation
  pt->x = (1.0-u)*e0.x + u*e1.x + (1.0-v)*e2.x + v*e3.x
    - ((1.0-u)*(1.0-v)*c[0].x + u*(1.0-v)*c[1].x + u*(1.0-v)*c[2].x + u*v*c[3].x);
  pt->y = (1.0-u)*e0.y + u*e1.y + (1.0-v)*e2.y + v*e3.y
    - ((1.0-u)*(1.0-v)*c[0].y + u*(1.0-v)*c[1].y + u*(1.0-v)*c[2].y + u*v*c[3].y);
  pt->z = (1.0-u)*e0.z + u*e1.z + (1.0-v)*e2.z + v*e3.z
    - ((1.0-u)*(1.0-v)*c[0].z + u*(1.0-v)*c[1].z + u*(1.0-v)*c[2].z + u*v*c[3].z);
}


/*
  Create the transfinite volume
*/
TMR_TFIVolume::TMR_TFIVolume( TMRPoint *corners, 
                              TMR_EdgeLookup *_edges[], 
                              TMR_SurfaceLookup *_faces[] ){
  // Copy the corners, edges and faces
  for ( int k = 0; k < 8; k++ ){
    c[k] = corners[k];
  }
  for ( int k = 0; k < 12; k++ ){
    edges[k] = _edges[k];
  }
  for ( int k = 0; k < 6; k++ ){
    faces[k] = _faces[k];
  }
}

/*
  Transfinite interpolation within the volume
*/
int TMR_TFIVolume::evalPoint( int32_t u_int, int32_t v_int, int32_t w_int,
                              TMRPoint *pt ){
  int fail = 0;
  // Evaluate/retrieve the points on the surfaces
  TMRPoint f[6];
  for ( int k = 0; k < 6; k++ ){
    if (k < 2){
      fail = fail || faces[k]->evalPoint(v_int, w_int, &f[k]);
    }
    else if (k < 4){
      fail = fail || faces[k]->evalPoint(u_int, w_int, &f[k]);
    }
    else {
      fail = fail || faces[k]->evalPoint(u_int, v_int, &f[k]);
    }
  }

  // Evaluate/retrieve the points on the edges
  TMRPoint e[12];
  for ( int k = 0; k < 12; k++ ){
    if (k < 4){
      fail = fail || edges[k]->evalPoint(u_int, &e[k]);
    }
    else if (k < 8){
      fail = fail || edges[k]->evalPoint(v_int, &e[k]);
    }
    else {
      fail = fail || edges[k]->evalPoint(w_int, &e[k]);
    }
  }

  // If any of the point evaluations failed, then everything else must
  // fail as well.
  if (fail){
    return fail;
  }

  // Convert the integer coordinates into real coordinates
  const double u = 1.0*u_int/(1 << TMR_MAX_LEVEL);
  const double v = 1.0*v_int/(1 << TMR_MAX_LEVEL);
  const double w = 1.0*w_int/(1 << TMR_MAX_LEVEL);

  // Evaluate the point based on the values along the corners, edges and faces
  pt->x = (1.0-u)*f[0].x + u*f[1].x + (1.0-v)*f[2].x + v*f[3].x + (1.0-w)*f[4].x + w*f[5].x 
    - ((1.0-v)*(1.0-w)*e[0].x + v*(1.0-w)*e[1].x + (1.0-v)*w*e[2].x + v*w*e[3].x +
       (1.0-u)*(1.0-w)*e[4].x + u*(1.0-w)*e[5].x + (1.0-u)*w*e[6].x + u*w*e[7].x +
       (1.0-u)*(1.0-v)*e[8].x + u*(1.0-v)*e[9].x + (1.0-u)*v*e[10].x + u*v*e[11].x)
    + ((1.0-u)*(1.0-v)*(1.0-w)*c[0].x + u*(1.0-v)*(1.0-w)*c[1].x + 
       (1.0-u)*v*(1.0-w)*c[2].x + u*v*(1.0-w)*c[3].x + 
       (1.0-u)*(1.0-v)*w*c[4].x + u*(1.0-v)*w*c[5].x + 
       (1.0-u)*v*w*c[6].x + u*v*w*c[7].x);

  pt->y = (1.0-u)*f[0].y + u*f[1].y + (1.0-v)*f[2].y + v*f[3].y + (1.0-w)*f[4].y + w*f[5].y 
    - ((1.0-v)*(1.0-w)*e[0].y + v*(1.0-w)*e[1].y + (1.0-v)*w*e[2].y + v*w*e[3].y +
       (1.0-u)*(1.0-w)*e[4].y + u*(1.0-w)*e[5].y + (1.0-u)*w*e[6].y + u*w*e[7].y +
       (1.0-u)*(1.0-v)*e[8].y + u*(1.0-v)*e[9].y + (1.0-u)*v*e[10].y + u*v*e[11].y)
    + ((1.0-u)*(1.0-v)*(1.0-w)*c[0].y + u*(1.0-v)*(1.0-w)*c[1].y + 
       (1.0-u)*v*(1.0-w)*c[2].y + u*v*(1.0-w)*c[3].y + 
       (1.0-u)*(1.0-v)*w*c[4].y + u*(1.0-v)*w*c[5].y + 
       (1.0-u)*v*w*c[6].y + u*v*w*c[7].y);

  pt->z = (1.0-u)*f[0].z + u*f[1].z + (1.0-v)*f[2].z + v*f[3].z + (1.0-w)*f[4].z + w*f[5].z 
    - ((1.0-v)*(1.0-w)*e[0].z + v*(1.0-w)*e[1].z + (1.0-v)*w*e[2].z + v*w*e[3].z +
       (1.0-u)*(1.0-w)*e[4].z + u*(1.0-w)*e[5].z + (1.0-u)*w*e[6].z + u*w*e[7].z +
       (1.0-u)*(1.0-v)*e[8].z + u*(1.0-v)*e[9].z + (1.0-u)*v*e[10].z + u*v*e[11].z)
    + ((1.0-u)*(1.0-v)*(1.0-w)*c[0].z + u*(1.0-v)*(1.0-w)*c[1].z + 
       (1.0-u)*v*(1.0-w)*c[2].z + u*v*(1.0-w)*c[3].z + 
       (1.0-u)*(1.0-v)*w*c[4].z + u*(1.0-v)*w*c[5].z + 
       (1.0-u)*v*w*c[6].z + u*v*w*c[7].z);

  return 0;
}
