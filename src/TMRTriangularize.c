#include "TMRTriangularize.h"
#include "TMRPerfectMatchInterface.h"
#include "TMRHashFunction.h"
#include "predicates.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>

// Define the maximum distance
const double MAX_QUAD_DISTANCE = 1e40;

/*
  Compare coordinate pairs of points. This uses a Morton ordering
  comparison to facilitate sorting/searching the list of edges.

  This function can be used by the stdlib functions qsort and
  bsearch to sort/search the edge pairs.
*/
static int compare_edges( const void *avoid, const void *bvoid ){
  // Cast the input to uint32_t types
  const uint32_t *a = static_cast<const uint32_t*>(avoid);
  const uint32_t *b = static_cast<const uint32_t*>(bvoid);
  
  // Extract the x/y locations for the a and b points
  uint32_t ax = a[0], ay = a[1];
  uint32_t bx = b[0], by = b[1];

  uint32_t xxor = ax ^ bx;
  uint32_t yxor = ay ^ by;
  uint32_t sor = xxor | yxor;

  // Note that here we do not distinguish between levels
  // Check for the most-significant bit
  int discrim = 0;
  if (xxor > (sor ^ xxor)){
    discrim = ax - bx;
  }
  else {
    discrim = ay - by;
  }

  return discrim;
}

/*
  A light-weight queue for keeping track of groups of triangles
*/
class TriQueueNode {
 public:
  TMRTriangle *tri;
  TriQueueNode *next;
};

/*
  The triangle priority queue

  This enables the storage of a list of triangles that can
  be accessed 
*/
class TriQueue {
 public:
  // Initialize the priority queue without and data
  TriQueue(){
    start = end = NULL;
    size = 0;
  }

  // Free all of the queue entries
  ~TriQueue(){
    while (start){
      TriQueueNode *tmp = start->next;
      delete start;
      start = tmp;      
    }
  }

  // Pop the first triangle pointer off the queue and delete its
  // corresponding node. This is used to extract the first node in the
  // list.
  TMRTriangle *pop(){
    TMRTriangle *t = NULL;
    if (start){
      t = start->tri;
      TriQueueNode *tmp = start;
      start = start->next;
      delete tmp;
      size--;
    }
    return t;
  }

  // Append the given triangle pointer to the end of the priority
  // queue. If the size of the queue is zero, then readjust the
  // start/end pointers and the size accordingly.
  void append( TMRTriangle *t ){
    if (size == 0){
      end = new TriQueueNode();
      end->tri = t;
      end->next = NULL;
      start = end;
      size++;
    }
    else {
      end->next = new TriQueueNode();
      end = end->next;
      end->tri = t;
      end->next = NULL;
      size++;
    }
  }
  
  // Delete the *next* entry in the list after the provided
  // pointer. Note that this is required since we only have a pointer
  // to the next entry in the list (and not the previous one). This
  // function returns 1 for success and 0 for failure.
  int deleteNext( TriQueueNode *node ){
    if (node && node->next){
      // Pointer to the entry to be deleted
      TriQueueNode *next = node->next;
      
      // Adjust the 'next' pointer to point past the entry that will
      // be deleted
      node->next = next->next;

      // Decrease the length of the queue
      size--;

      // Delete the entry
      delete next;
      return 1;
    }

    // The node cannot be deleted
    return 0;
  }

  // Note that we grant public access to these members. Use this
  // access wisely!
  int size;
  TriQueueNode *start;
  TriQueueNode *end;
};

/*
  The following class implments a simple quadtree data structure for fast
  geometric searching queries.
*/
TMRQuadNode::TMRQuadNode( TMRQuadDomain *_domain ){
  initialize(_domain, 0, 0, 0);
}

/*
  Create a child node
*/
TMRQuadNode::TMRQuadNode( TMRQuadDomain *_domain,
                          uint32_t _u, uint32_t _v, int _level ){
  initialize(_domain, _u, _v, _level);
}

/*
  Initialize the data for the quadtree root/node/leaf
*/
void TMRQuadNode::initialize( TMRQuadDomain *_domain,
                              uint32_t _u, uint32_t _v, int _level ){
  // Set the domain
  domain = _domain;

  // Set the level of this quadndoe
  level = _level;

  // Set the domain level
  u = _u;
  v = _v;

  // Set the maximum length in parameter space for
  // the domain
  const uint32_t hmax = 1 << MAX_DEPTH;
  const uint32_t h = 1 << (MAX_DEPTH - level - 1);

  // Compute the distance along each edge for this node
  // within the quadtree
  double ax = 1.0*(u + h)/hmax;
  double ay = 1.0*(v + h)/hmax;

  // Set the new planes that spatially divide this node
  x = (1.0 - ax)*domain->xlow + ax*domain->xhigh;
  y = (1.0 - ay)*domain->ylow + ay*domain->yhigh;

  // Set the pointers to the different children (all to NULL for now)
  low_left = NULL;
  low_right = NULL;
  up_left = NULL;
  up_right = NULL;

  // Allocate space for the points that will be added
  num_points = 0;
  pts = new double[ 2*NODES_PER_LEVEL ];
  pt_nums = new uint32_t[ NODES_PER_LEVEL ];
}

/*
  Free the quadtree node and its data
*/
TMRQuadNode::~TMRQuadNode(){
  if (pts){ delete [] pts; }
  if (pt_nums){ delete [] pt_nums; }
  if (low_left){
    delete low_left;
    delete low_right;
    delete up_left;
    delete up_right;
  }
}

/*
  Add a node to the quadtree.

  This code does not check for duplicated geometric entities at the
  same point or duplicated indices. It is the user's responsibility
  not to add duplicates.
*/
void TMRQuadNode::addNode( uint32_t num, const double pt[] ){
  // If any of the children have been allocated, searh them for
  // the place where the node should be added
  if (low_left){    
    if (pt[0] < x && pt[1] < y){
      low_left->addNode(num, pt);
    }
    else if (pt[0] < x){
      up_left->addNode(num, pt);
    }
    else if (pt[1] < y){
      low_right->addNode(num, pt);
    }
    else {
      up_right->addNode(num, pt);
    }
    return;
  }

  // Check if adding this node would exceed the number of nodes that
  // are allowed per level
  if (num_points < NODES_PER_LEVEL){
    pts[2*num_points] = pt[0];
    pts[2*num_points+1] = pt[1];
    pt_nums[num_points] = num;
    num_points++;
    return;
  }
  else {
    // Allocate new children for the new nodes
    const uint32_t h = 1 << (MAX_DEPTH - level - 1);

    low_left = new TMRQuadNode(domain, u, v, level+1);
    low_right = new TMRQuadNode(domain, u+h, v, level+1);
    up_left = new TMRQuadNode(domain, u, v+h, level+1);
    up_right = new TMRQuadNode(domain, u+h, v+h, level+1);

    // Add the points to myself
    for ( int k = 0; k < NODES_PER_LEVEL; k++ ){
      addNode(pt_nums[k], &pts[2*k]);
    }

    // Free the space that was allocated
    delete [] pt_nums;
    delete [] pts;
    pt_nums = NULL;
    pts = NULL;
    num_points = 0;
  }
}

/*
  Delete a point from the quadtree.

  Note that the coordinate provided is used to match the point, but
  needs to be used to scan through the quadtree hierarchy to find the
  correct quadrant location.
*/
int TMRQuadNode::deleteNode( uint32_t num, const double pt[] ){
  if (low_left){
    if (pt[0] < x && pt[1] < y){
      return low_left->deleteNode(num, pt);
    }
    else if (pt[0] < x){
      return up_left->deleteNode(num, pt);
    }
    else if (pt[1] < y){
      return low_right->deleteNode(num, pt);
    }
    else {
      return up_right->deleteNode(num, pt);
    }
  }

  // The node could not be deleted because no nodes exist on this
  // level.
  if (num_points == 0){ return 0; }

  // Scan through the list of nodes, and check for the one that
  // needs to be deleted.
  int i = 0;
  for ( ; i < num_points; i++ ){
    if (pt_nums[i] == num){
      break;
    }
  }

  // We ended up past the last entry, therefore the point was not found
  if (i == num_points){ return 0; }

  // Increment the index one past the node that will be eliminated
  i++;
  for ( ; i < num_points; i++ ){
    // Move all the data back to cover the deleted entry
    pt_nums[i-1] = pt_nums[i];
    pts[2*(i-1)] = pts[2*i];
    pts[2*(i-1)+1] = pts[2*i+1];
  }
  // Decrease the number of points to account for the deleted node
  num_points--;
    
  // The point was found
  return 1;
}

/*
  Find the closest indexed point to the provided (x,y) location

  This function is the top-level function that calls the recursive
  version after initializing the data.
*/
uint32_t TMRQuadNode::findClosest( const double pt[], double *_dist ){
  double dist = MAX_QUAD_DISTANCE;
  uint32_t index = 0;

  // Make the recursive call
  findClosest(pt, &index, &dist);
  if (_dist){ *_dist = dist; }
  return index;
}

/*
  The recursive call to find the closest point

  This scans each quadrant, testing whether it will be necessary to
  scan multiple quadrants to ensure that we have, in fact, located the
  closest point to the provided point.
*/
void TMRQuadNode::findClosest( const double pt[], uint32_t *index, double *dist ){
  // Scan through the quadtree
  if (low_left){
    if (pt[0] < x && pt[1] < y){
      low_left->findClosest(pt, index, dist);
      if (x - pt[0] <= *dist){
        low_right->findClosest(pt, index, dist);
      }
      if (y - pt[1] <= *dist){
        up_left->findClosest(pt, index, dist);
      }
      if (x - pt[0] <= *dist || y - pt[1] <= *dist){
        up_right->findClosest(pt, index, dist);
      }
    }
    else if (pt[0] < x){
      up_left->findClosest(pt, index, dist);
      if (x - pt[0] <= *dist){
        up_right->findClosest(pt, index, dist);
      }
      if (pt[1] - y <= *dist){
        low_left->findClosest(pt, index, dist);
      }
      if (x - pt[0] <= *dist || pt[1] - y <= *dist){
        low_right->findClosest(pt, index, dist);
      }
    }
    else if (pt[1] < y){
      low_right->findClosest(pt, index, dist);
      if (pt[0] - x <= *dist){
        low_left->findClosest(pt, index, dist);
      }
      if (y - pt[1] <= *dist){
        up_right->findClosest(pt, index, dist);
      }
      if (pt[0] - x <= *dist || y - pt[1] <= *dist){
        up_left->findClosest(pt, index, dist);
      }
    }
    else {
      up_right->findClosest(pt, index, dist);
      if (pt[0] - x <= *dist){
        up_left->findClosest(pt, index, dist);
      }
      if (pt[1] - y <= *dist){
        low_right->findClosest(pt, index, dist);
      }
      if (pt[0] - x <= *dist || pt[1] - y <= *dist){
        low_left->findClosest(pt, index, dist);
      }
    }
    return;
  }

  // This is a leaf, search the points
  for ( int i = 0; i < num_points; i++ ){
    double dx = pt[0] - pts[2*i];
    double dy = pt[1] - pts[2*i+1];
    double d = sqrt(dx*dx + dy*dy);
    if (d < *dist){
      *dist = d;
      *index = pt_nums[i];
    }
  }
}

/*
  Create the triangularization object. 

  This code uses Schewchuk's geometric predicates for testing whether
  points lie within a triangle's circumcircle or whether a point forms
  a properly oriented triangle (orient2d).

  This object takes in a series of points, a nuber of holes, and a
  series of segments. The points are in parametric space while the
  face defines the mapping between parameter space and physical
  space. Meshing is performed by using a local metric within each
  triangle to define distance. The incircle operations that test
  whether a point lies within the circumcircle of a triangle are
  modified to reflect this local metric. (This is not completely
  justified theoretically, but it works in practice.) The orient2d
  operations still take place in parameter space, since their purpose
  is to test whether the point lies within the triangle.

  The holes in the domain must be specified so that the triangulation
  does not cover them. The boundary of the hole itself must be
  represented by the segments and must be closed. The points
  indicating the location of the holes must be numbered last in the
  point list and cannot be included in the segment edges.

  If no face object is provided, then the meshing takes place in
  the parameter space directly.

  input:
  npts:   the number of points (including holes) in the domain
  inpts:  the parametric input points
  nholes: the number of holes that are specified
  nsegs:  the number of segments
  segs:   the segments: consecutive point numbers indicating edges
  surf:   the optional face
*/
TMRTriangularize::TMRTriangularize( int npts, const double inpts[],
                                    int nsegs, const int segs[],
                                    TMRFace *surf ){
  initialize(npts, inpts, 0, nsegs, segs, surf);
}

/*
  Alternate constructor: Triangularize with holes
*/
TMRTriangularize::TMRTriangularize( int npts, const double inpts[], int nholes,
                                    int nsegs, const int segs[],
                                    TMRFace *surf ){
  initialize(npts, inpts, nholes, nsegs, segs, surf);
}


/*
  The 'true' constructor. This is called, depending on which arguments
  are provided.
*/
void TMRTriangularize::initialize( int npts, const double inpts[], int nholes,
                                   int nsegs, const int segs[],
                                   TMRFace *surf ){
  // Initialize the predicates code
  exactinit();

  // Set the surface
  face = surf;
  if (face){
    face->incref();
  }

  // Allocate and initialize/zero the hash table data for the edges
  num_hash_nodes = 0;
  num_buckets = 100;
  buckets = new EdgeHashNode*[ num_buckets ];
  memset(buckets, 0, num_buckets*sizeof(EdgeHashNode*));

  // Allocate and initialize/zero the hash table data for the edges
  num_active_triangles = 0;
  num_active_buckets = 100;
  active_buckets = new ActiveHashNode*[ num_active_buckets ];
  memset(active_buckets, 0, num_active_buckets*sizeof(ActiveHashNode*));

  // Set the initial root/current location of the doubly-linked
  // triangle list structure. These are allocated as we add new triangles.
  list_start = NULL;
  list_end = NULL;
  list_marker = NULL;

  // Keep track of the total number of triangles
  num_triangles = 0;

  // Set the initial number of points
  init_boundary_points = npts - nholes;
  num_points = FIXED_POINT_OFFSET + npts;

  // Set the initial for the maximum number of points
  max_num_points = 1024;
  if (max_num_points < num_points){
    max_num_points = num_points;
  }

  // Allocate the initial set of points
  pts = new double[ 2*max_num_points ];
  pts_to_tris = new TMRTriangle*[ max_num_points ];

  // If we have a face object
  if (face){
    X = new TMRPoint[ max_num_points ];
  }
  else {
    X = NULL;
  }

  // Find the maximum domain size
  domain.xlow = domain.xhigh = inpts[0];
  domain.ylow = domain.yhigh = inpts[1];
  for ( int i = 1; i < npts; i++ ){
    if (inpts[2*i] < domain.xlow){
      domain.xlow = inpts[2*i];
    }
    if (inpts[2*i+1] < domain.ylow){
      domain.ylow = inpts[2*i+1];
    }
    if (inpts[2*i] > domain.xhigh){
      domain.xhigh = inpts[2*i];
    }
    if (inpts[2*i+1] > domain.yhigh){
      domain.yhigh = inpts[2*i+1];
    }    
  }

  // Re-adjust the domain boundary to ensure that it is sufficiently
  // large
  domain.xhigh += 1.0;
  domain.xlow -= 1.0;
  domain.yhigh += 1.0;
  domain.ylow -= 1.0;

  // Allocate the new root mode
  root = new TMRQuadNode(&domain);
  search_tag = 0;

  // Set up the PSLG edges
  setUpPSLGEdges(nsegs, segs);

  // Set the initial points
  num_points = FIXED_POINT_OFFSET;

  // Set the point (xlow, ylow)
  pts[0] = domain.xlow;
  pts[1] = domain.ylow;
  // Set the point (xhigh, ylow)
  pts[2] = domain.xhigh;
  pts[3] = domain.ylow;
  // Set the point (xlow, yhigh)
  pts[4] = domain.xlow;
  pts[5] = domain.yhigh;
  // Set the point (xhigh, yhigh)
  pts[6] = domain.xhigh;
  pts[7] = domain.yhigh;
  
  // Evaluate the points on the face
  if (face){
    for ( int i = 0; i < num_points; i++ ){
      face->evalPoint(pts[2*i], pts[2*i+1], &X[i]);
    }
  }

  // Add the extreme points to the quadtree
  for ( int i = 0; i < FIXED_POINT_OFFSET; i++ ){
    root->addNode(i, &pts[2*i]);
  }

  // Add the initial triangles
  addTriangle(TMRTriangle(0, 1, 2));
  addTriangle(TMRTriangle(2, 1, 3));

  // Add the points to the triangle. This creates a CDT of the
  // original set of points.
  for ( int i = 0; i < npts; i++ ){
    addPointToMesh(&inpts[2*i]);
  }

  // Set the triangle tags to zero
  setTriangleTags(0);

  // Mark all the triangles in the list that contain or touch nodes
  // that are in the FIXED_POINT_OFFSET list that are not separated by
  // a PSLG edge. These triangles will be deleted.
  TriListNode *node = list_start;
  while (node){
    if (node->tri.tag == 0 && 
        ((node->tri.u < FIXED_POINT_OFFSET ||
          node->tri.v < FIXED_POINT_OFFSET ||
          node->tri.w < FIXED_POINT_OFFSET) ||
         (node->tri.u >= num_points - nholes ||
          node->tri.v >= num_points - nholes ||
          node->tri.w >= num_points - nholes))){
      tagTriangles(&node->tri);
    }
    node = node->next;
  }

  // Free the triangles that we've found
  list_marker = NULL;
  node = list_start;
  while (node){
    // Keep track of what node will be next, since we may be deleting
    // 'node' itself, we cannot access this afterwards
    TriListNode *tmp = node->next;

    // Delete the triangle from the triangle list
    if (node->tri.tag == 1){
      deleteTriangle(node->tri);
    }

    // Set the next node
    node = tmp;
  }

  // Free the points and holes from the quadtree
  for ( uint32_t num = 0; num < FIXED_POINT_OFFSET; num++ ){
    root->deleteNode(num, &pts[2*num]);
  }

  // Free the points associated with the number of holes
  for ( uint32_t num = num_points - nholes; num < num_points; num++ ){
    root->deleteNode(num, &pts[2*num]);
  }
  num_points -= nholes;

  // Reset the node->traingle pointers to avoid referring to a
  // triangle that belonged to a hole and was deleted.
  node = list_start;
  while (node){
    pts_to_tris[node->tri.u] = &(node->tri);
    pts_to_tris[node->tri.v] = &(node->tri);
    pts_to_tris[node->tri.w] = &(node->tri);
    node = node->next;
  }
}

/*
  Free the triangularization object
*/
TMRTriangularize::~TMRTriangularize(){
  delete root;
  delete [] pts;
  if (X){ delete [] X; }
  delete [] pts_to_tris;

  // Dereference the face
  if (face){
    face->decref();
  }

  // Free the PSLG edges
  if (pslg_edges){
    delete [] pslg_edges;
  }

  // Free the data for the edge hash table
  for ( int i = 0; i < num_buckets; i++ ){
    EdgeHashNode *node = buckets[i];
    while (node){
      EdgeHashNode *tmp = node;
      node = node->next;
      delete tmp;
    }
  }

  // Free the data for the active triangles
  for ( int i = 0; i < num_active_buckets; i++ ){
    ActiveHashNode *node = active_buckets[i];
    while (node){
      ActiveHashNode *tmp = node;
      node = node->next;
      delete tmp;
    }
  }

  // Free the doubly linked edge list
  while (list_start){
    TriListNode *tmp = list_start;
    list_start = list_start->next;
    delete tmp;
  }
}

/*
  Retrieve the underlying mesh
*/
void TMRTriangularize::getMesh( int *_num_points,
                                int *_num_triangles, 
                                int **_conn, double **_pts,
                                TMRPoint **_X ){
  int npts = (num_points - FIXED_POINT_OFFSET);
  *_num_points = npts;
  *_num_triangles = num_triangles;
  *_conn = new int[ 3*num_triangles ];
  *_pts = new double[ 2*npts ];
  if (face){
    *_X = new TMRPoint[ npts ];
  }

  // Set the points
  memcpy(*_pts, &pts[2*FIXED_POINT_OFFSET], 2*npts*sizeof(double));
  if (face){
    memcpy(*_X, &X[FIXED_POINT_OFFSET], npts*sizeof(TMRPoint));
  }

  // Set the pointer into the connectivity array
  int *t = *_conn;

  // Determine the connectivity
  TriListNode *node = list_start;
  while (node){
    t[0] = node->tri.u - FIXED_POINT_OFFSET;
    t[1] = node->tri.v - FIXED_POINT_OFFSET;
    t[2] = node->tri.w - FIXED_POINT_OFFSET;
    t += 3;    
    node = node->next;
  }
}

/*
  Reset the tags of all the triangles within the list
*/
void TMRTriangularize::setTriangleTags( uint32_t tag ){
  TriListNode *node = list_start;
  while (node){
    TriListNode *next = node->next;
    node->tri.tag = tag;
    node = node->next;
  }
}

/*
  Mark triangles that should be deleted.

  This is used to mark triangles that are in holes, or have points that
  lie outside the
*/
void TMRTriangularize::tagTriangles( TMRTriangle *tri ){
  // Set the combinations of edge pairs that will be added
  // to the hash table
  uint32_t edge_pairs[][2] = {{tri->u, tri->v},
                              {tri->v, tri->w},
                              {tri->w, tri->u}};
  for ( int k = 0; k < 3; k++ ){
    if (!edgeInPSLG(edge_pairs[k][0], edge_pairs[k][1])){
      TMRTriangle *t;
      completeMe(edge_pairs[k][1], edge_pairs[k][0], &t);
      if (t && t->tag == 0){
        t->tag = 1;
        tagTriangles(t);
      }
    }
  }                           
}

/*
  Write out the triangularization to a file
*/
void TMRTriangularize::writeToVTK( const char *filename ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
    // Write out the points
    fprintf(fp, "POINTS %d float\n", num_points);
    if (face){
      for ( int k = 0; k < num_points; k++ ){
        fprintf(fp, "%e %e %e\n", X[k].x, X[k].y, X[k].z);
      }
    }
    else {
      for ( int k = 0; k < num_points; k++ ){
        fprintf(fp, "%e %e %e\n", pts[2*k], pts[2*k+1], 0.0);
      }
    }
    
    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", num_triangles, 4*num_triangles);

    TriListNode *node = list_start;
    while (node){
      fprintf(fp, "3 %d %d %d\n", node->tri.u, node->tri.v, node->tri.w);
      node = node->next;
    }

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", num_triangles);
    for ( int k = 0; k < num_triangles; k++ ){
      fprintf(fp, "%d\n", 5);
    }

    // Print out the rest as fields one-by-one
    fprintf(fp, "CELL_DATA %d\n", num_triangles);
    fprintf(fp, "SCALARS status float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    node = list_start;
    while (node){
      fprintf(fp, "%d\n", node->tri.status);
      node = node->next;
    }

    fclose(fp);
  }
}

/*
  Retrieve the edge hash
*/
inline uint32_t TMRTriangularize::getEdgeHash( uint32_t x, uint32_t y ){
  return TMRIntegerPairHash(x, y);
}

/*
  Add a triangle to the mesh
*/
int TMRTriangularize::addTriangle( TMRTriangle tri ){
  int success = 1;

  // Add the triangle to the list of nodes
  if (!list_start){
    // Create the root triangle node
    list_start = new TriListNode();
    list_start->tri = tri;
    list_start->next = NULL;
    list_start->prev = NULL;

    // Update the current node
    list_end = list_start;
  }
  else {
    // Create the new member in the triangle list
    TriListNode *next = new TriListNode();
    next->tri = tri;
    next->next = NULL;
    next->prev = list_end;

    // Update the current member in the triangle list
    list_end->next = next;
    list_end = list_end->next;
  }

  list_end->tri.tag = 0;
  list_end->tri.status = NO_STATUS;

  // Set the pointer to the list of triangles
  pts_to_tris[tri.u] = &(list_end->tri);
  pts_to_tris[tri.v] = &(list_end->tri);
  pts_to_tris[tri.w] = &(list_end->tri);

  // Add the triangle to the triangle count
  num_triangles++;

  // Redistribute the members in the hash table if required
  if (num_hash_nodes > 10*num_buckets){
    // Create a new array of buckets, twice the size of the old array
    // of buckets and re-insert the entries back into the new array
    int num_old_buckets = num_buckets;
    num_buckets *= 2;

    // Create a new array pointing to the new buckets and the
    // end of the bucket array
    EdgeHashNode **new_buckets = new EdgeHashNode*[ num_buckets ];
    EdgeHashNode **end_buckets = new EdgeHashNode*[ num_buckets ];
    memset(new_buckets, 0, num_buckets*sizeof(EdgeHashNode*));
    memset(end_buckets, 0, num_buckets*sizeof(EdgeHashNode*));

    // Assign the new buckets
    for ( int i = 0; i < num_old_buckets; i++ ){
      EdgeHashNode *node = buckets[i];

      while (node){
        // Get the new hash values
        EdgeHashNode *tmp = node->next;
        uint32_t value = getEdgeHash(node->u, node->v);
        uint32_t bucket = value % num_buckets;
        
        // If the new bucket linked list does not exist
        if (!new_buckets[bucket]){
          new_buckets[bucket] = new EdgeHashNode();
          new_buckets[bucket]->next = NULL;
          new_buckets[bucket]->u = node->u;
          new_buckets[bucket]->v = node->v;
          new_buckets[bucket]->tri_node = node->tri_node;
          delete node;

          end_buckets[bucket] = new_buckets[bucket];

          // new_buckets[bucket] = node;
          // node->next = NULL;
          // end_buckets[bucket] = new_buckets[bucket];
        }
        else {
          end_buckets[bucket]->next = new EdgeHashNode();
          end_buckets[bucket] = end_buckets[bucket]->next;
          end_buckets[bucket]->next = NULL;
          end_buckets[bucket]->u = node->u;
          end_buckets[bucket]->v = node->v;
          end_buckets[bucket]->tri_node = node->tri_node;
          delete node;

          // end_buckets[bucket]->next = node;
          // end_buckets[bucket] = end_buckets[bucket]->next;
          // node->next = NULL;
        }
        
        node = tmp;
      }
    }

    delete [] end_buckets;
    delete [] buckets;
    buckets = new_buckets;
  }

  // Set the combinations of edge pairs that will be added
  // to the hash table
  uint32_t edge_pairs[][2] = {{tri.u, tri.v},
                              {tri.v, tri.w},
                              {tri.w, tri.u}};
                              
  // Add the triangle to the hash table
  for ( int k = 0; k < 3; k++ ){
    // Add a hash for each pair of edges around the triangle
    uint32_t u = edge_pairs[k][0];
    uint32_t v = edge_pairs[k][1];
    uint32_t value = getEdgeHash(u, v);
    uint32_t bucket = value % num_buckets;
    if (!buckets[bucket]){
      // Create the new buckets node and assign the values
      buckets[bucket] = new EdgeHashNode();
      buckets[bucket]->u = u;
      buckets[bucket]->v = v;
      buckets[bucket]->tri_node = list_end;
      buckets[bucket]->next = NULL;
      num_hash_nodes++;
    }
    else {
      // Scan through to the end of the array and set the pointer to
      // the triangle node that contains the 
      EdgeHashNode *node = buckets[bucket];
      while (node){
        // Check whether the edge is already in the hash table
        if (node->u == u && node->v == v){
          // The edge already exists, it was overwritten, but we'll
          // call this a failure...
          success = 0;

          // Overwite the triangle
          node->tri_node = list_end;

          // Set the pointer to NULL so that the new node will not be
          // created
          node = NULL;
          break;
        }

        // If the next node does not exist, then break 
        if (!node->next){
          break;
        }

        // Increment the node to the next position
        node = node->next;
      }

      if (node){
        // Create the new hash node
        node->next = new EdgeHashNode();
        node = node->next;

        // Set the data into the object
        node->u = u;
        node->v = v;
        node->tri_node = list_end;
        node->next = NULL;
        num_hash_nodes++;
      }
      else {
        success = 0;
      }
    }
  }

  return success;
}

/*
  Delete a triangle from the mesh
*/
int TMRTriangularize::deleteTriangle( TMRTriangle tri ){
  // Keep track of whether we successfully delete all of the edges, or
  // just some of the edges (1 or 2 out of 3 is bad!)
  int success = 1;

  // Delete the triangle if the status is active
  if (tri.status == ACTIVE){
    deleteActiveTriangle(&tri);
  }

  // Set the combinations of edge pairs that will be added
  // to the hash table
  uint32_t edge_pairs[][2] = {{tri.u, tri.v},
                              {tri.v, tri.w},
                              {tri.w, tri.u}};

  // Keep track if this is the first edge we find
  int first = 1;
  
  // Add the triangle to the hash table
  for ( int k = 0; k < 3; k++ ){
    // Add a hash for each pair of edges around the triangle
    uint32_t u = edge_pairs[k][0];
    uint32_t v = edge_pairs[k][1];

    // Keep track of whether we've deleted this edge
    int edge_success = 0;

    // Get the hash values and access the first entry of the bucket
    // it's listed under
    uint32_t value = getEdgeHash(u, v);
    uint32_t bucket = value % num_buckets;
    EdgeHashNode *node = buckets[bucket];
    EdgeHashNode *prev = node;
    while (node){
      // The edge matches, we have to delete this triangle
      if (u == node->u && v == node->v){
        // If this is the first edge we've found, delete the
        // corresponding list entry as well.
        if (first){
          // This triangle will be successfully deleted, adjust the
          // triangle count
          num_triangles--;

          // Adjust the double linking so that when the node is
          // deleted, the pointers will still work
          TriListNode *ptr = node->tri_node;
          if (ptr == list_start){
            list_start = list_start->next;
            if (list_start){
              list_start->prev = NULL;
            }
          }
          else if (ptr == list_end){
            list_end = list_end->prev;
            if (list_end){
              list_end->next = NULL;
            }
          }
          else {
            // Set the pointers from the previous and next objects
            // in the list so that they point past the ptr member
            if (ptr->prev){
              ptr->prev->next = ptr->next;
            }
            if (ptr->next){
              ptr->next->prev = ptr->prev;
            }
          }

          // Adjust the list marker if we're about to delete the
          // triangle that the marker points to. Move the marker one
          // triangle back (or set it to NULL if it's the beginning of
          // the list)
          if (list_marker && ptr == list_marker){
            list_marker = ptr->prev;
          }
          delete ptr;
          first = 0;
        }

        // Delete the edge from the hash table
        if (node == prev){
          buckets[bucket] = node->next;
          delete node;
          num_hash_nodes--;
        }
        else {
          prev->next = node->next;
          delete node;
          num_hash_nodes--;
        }

        edge_success = 1;
        break;
      }

      prev = node;
      node = node->next;
    }

    // Keep track of each individual edge
    success = success && edge_success;
  }

  return success;
}

/*
  Retrieve the triangle hash
*/
inline uint32_t TMRTriangularize::getTriangleHash( TMRTriangle *tri ){
  return TMRIntegerTripletHash(tri->u, tri->v, tri->w);
}

/*
  Add a triangle to the mesh
*/
int TMRTriangularize::addActiveTriangle( TMRTriangle *tri ){
  int success = 1;

  // Redistribute the members in the hash table if required
  if (num_active_triangles > 10*num_active_buckets){
    // Create a new array of buckets, twice the size of the old array
    // of buckets and re-insert the entries back into the new array
    int num_old_buckets = num_active_buckets;
    num_active_buckets *= 2;

    // Create a new array pointing to the new buckets and the
    // end of the bucket array
    ActiveHashNode **new_buckets = new ActiveHashNode*[ num_active_buckets ];
    ActiveHashNode **end_buckets = new ActiveHashNode*[ num_active_buckets ];
    memset(new_buckets, 0, num_active_buckets*sizeof(ActiveHashNode*));
    memset(end_buckets, 0, num_active_buckets*sizeof(ActiveHashNode*));

    // Assign the new buckets
    for ( int i = 0; i < num_old_buckets; i++ ){
      ActiveHashNode *node = active_buckets[i];

      while (node){
        // Get the new hash values
        ActiveHashNode *tmp = node->next;
        uint32_t value = getTriangleHash(node->tri);
        uint32_t bucket = value % num_active_buckets;
        
        // If the new bucket linked list does not exist
        if (!new_buckets[bucket]){
          new_buckets[bucket] = new ActiveHashNode();
          new_buckets[bucket]->next = NULL;
          new_buckets[bucket]->tri = node->tri;
          delete node;
          end_buckets[bucket] = new_buckets[bucket];
        }
        else {
          end_buckets[bucket]->next = new ActiveHashNode();
          end_buckets[bucket] = end_buckets[bucket]->next;
          end_buckets[bucket]->next = NULL;
          end_buckets[bucket]->tri = node->tri;
          delete node;
        }
        
        node = tmp;
      }
    }

    delete [] end_buckets;
    delete [] active_buckets;
    active_buckets = new_buckets;
  }

  // Add the triangle to the hash table
  uint32_t value = getTriangleHash(tri);
  uint32_t bucket = value % num_active_buckets;
  if (!active_buckets[bucket]){
    // Create the node and assign the values
    active_buckets[bucket] = new ActiveHashNode();
    active_buckets[bucket]->tri = tri;
    active_buckets[bucket]->next = NULL;

    // Increase the number of active triangles
    num_active_triangles++;
  }
  else {
    // Scan through to the end of the array and set the pointer to
    // the triangle node that contains the 
    ActiveHashNode *node = active_buckets[bucket];
    while (node){
      // Check whether the edge is already in the hash table
      if (node->tri->u == tri->u && 
          node->tri->v == tri->v && 
          node->tri->w == tri->w){
        // The triangle already exists. Replace the pointer
        // and quit. Call this a failure...
        success = 0;

        // Overwite the triangle
        node->tri = tri;

        // Set the pointer to NULL so that the new node will not be
        // created
        node = NULL;
        break;
      }

      // If the next node does not exist, then break 
      if (!node->next){
        break;
      }

      // Increment the node to the next position
      node = node->next;
    }

    if (node){
      // Create the new hash node
      node->next = new ActiveHashNode();
      node = node->next;

      // Set the data into the object
      node->tri = tri;
      node->next = NULL;

      // Increase the number of active triangles
      num_active_triangles++;
    }
    else {
      success = 0;
    }
  }

  return success;
}

/*
  Delete a triangle from the mesh
*/
int TMRTriangularize::deleteActiveTriangle( TMRTriangle *tri ){
  int success = 1;

  // Get the hash values and access the first entry of the bucket
  // it's listed under
  uint32_t value = getTriangleHash(tri);
  uint32_t bucket = value % num_active_buckets;
  ActiveHashNode *node = active_buckets[bucket];
  ActiveHashNode *prev = node;
  
  while (node){
    // The edge matches, we have to delete this triangle
    if (tri->u == node->tri->u && 
        tri->v == node->tri->v && 
        tri->w == node->tri->w){
      // Delete the edge from the hash table
      if (node == prev){
        active_buckets[bucket] = node->next;
        delete node;
        num_active_triangles--;
      }
      else {
        prev->next = node->next;
        delete node;
        num_active_triangles--;
      }

      success = 1;
      break;
    }

    prev = node;
    node = node->next;
  }

  return success;
}

/*
  Complete the adjacent triangle

  You complete me: Find the triangle that completes the specified
  edge. This can be used to find the triangle that is adjacent to
  another triangle.
*/
void TMRTriangularize::completeMe( uint32_t u, uint32_t v,
                                   TMRTriangle **tri ){
  *tri = NULL;

  // Retrieve the hash value/bucket for this edge
  uint32_t value = getEdgeHash(u, v);
  uint32_t bucket = value % num_buckets;
  EdgeHashNode *node = buckets[bucket];

  // Loop over the edges until we find the matching triangle
  while (node){
    if (node->u == u && node->v == v){
      *tri = &(node->tri_node->tri);
      break;
    }

    // Increment the pointer to the next member of this bucket
    node = node->next;
  }
}

/*
  Create the list of PSLG edges
*/
void TMRTriangularize::setUpPSLGEdges( int nsegs, const int segs[] ){
  num_pslg_edges = 2*nsegs;
  pslg_edges = new uint32_t[ 2*num_pslg_edges ];

  for ( int i = 0; i < nsegs; i++ ){
    uint32_t u = 0, v = 0;
    if (segs[2*i] >= 0){
      u = segs[2*i] + FIXED_POINT_OFFSET;
    }
    if (segs[2*i+1] >= 0){
      v = segs[2*i+1] + FIXED_POINT_OFFSET;
    }

    pslg_edges[4*i] = u;
    pslg_edges[4*i+1] = v;
    pslg_edges[4*i+2] = v;
    pslg_edges[4*i+3] = u;
  }

  qsort(pslg_edges, num_pslg_edges, 2*sizeof(uint32_t), compare_edges);
}

/*
  Search the list of sorted edges
*/
int TMRTriangularize::edgeInPSLG( uint32_t u, uint32_t v ){
  uint32_t edge[2];
  edge[0] = u;  edge[1] = v;

  int result = (bsearch(edge, pslg_edges, num_pslg_edges, 2*sizeof(uint32_t),
                        compare_edges) != NULL);

  return result;
}

/*
  Does the given triangle enclose the point?
*/
inline int TMRTriangularize::enclosed( const double p[],
                                       uint32_t u, uint32_t v, uint32_t w ){
  double pt[2] = {p[0], p[1]};
  if (orient2d(&pts[2*u], &pts[2*v], pt) >= 0.0 &&
      orient2d(&pts[2*v], &pts[2*w], pt) >= 0.0 &&
      orient2d(&pts[2*w], &pts[2*u], pt) >= 0.0){
    return 1;
  }

  return 0;
}

/*
  Does the final given point lie within the circumcircle of the remaining
  points?
*/
inline double TMRTriangularize::inCircle( uint32_t u, uint32_t v,
                                          uint32_t w, uint32_t x ){

  if (face){
    // Compute the mid-point of the triangle
    double mpt[2];
    const double frac = 1.0/3.0;
    mpt[0] = frac*(pts[2*u] + pts[2*v] + pts[2*w]);
    mpt[1] = frac*(pts[2*u+1] + pts[2*v+1] + pts[2*w+1]);

    // Compute the metric components at the center of the triangle (u,v,w)
    TMRPoint Xu, Xv;
    face->evalDeriv(mpt[0], mpt[1], &Xu, &Xv);
    double g11 = Xu.dot(Xu);
    double g12 = Xu.dot(Xv);
    double g22 = Xv.dot(Xv);

    // Compute a multiplicative decomposition such that G = L*L^{T}
    // [l11    ][l11 l21] = [g11  g12]
    // [l21 l22][    l22]   [g12  g22]
    // l11 = sqrt(g11)
    // l11*l21 = g12 => l21 = g12/l11
    // l21*l21 + l22^2 = g22 => l22 = sqrt(g22 - l21*l21);
    double l11 = sqrt(g11);
    double inv11 = 1.0/l11;
    double l21 = inv11*g12;
    double l22 = sqrt(g22 - l21*l21);

    // Compute p' = L^{T}*p to transform into the local coordinates
    double pu[2];
    pu[0] = l11*pts[2*u] + l21*pts[2*u+1];
    pu[1] = l22*pts[2*u+1];

    double pv[2];
    pv[0] = l11*pts[2*v] + l21*pts[2*v+1];
    pv[1] = l22*pts[2*v+1];

    double pw[2];
    pw[0] = l11*pts[2*w] + l21*pts[2*w+1];
    pw[1] = l22*pts[2*w+1];

    double px[2];
    px[0] = l11*pts[2*x] + l21*pts[2*x+1];
    px[1] = l22*pts[2*x+1];

    return incircle(pu, pv, pw, px);
  }

  return incircle(&pts[2*u], &pts[2*v], &pts[2*w], &pts[2*x]);
}

/*
  Add a point to the point list
*/
uint32_t TMRTriangularize::addPoint( const double pt[] ){
  // If the size of the array is exceeded, multiply it by 2
  if (num_points >= max_num_points){
    max_num_points *= 2;

    // Allocate a new array for the parametric locations of all the
    // points
    double *new_pts = new double[ 2*max_num_points ];
    memcpy(new_pts, pts, 2*num_points*sizeof(double));
    delete [] pts;
    pts = new_pts;

    // Allocate a new array for the pointer from the triangle vertices
    // to an attaching triangle
    TMRTriangle **new_pts_to_tris = new TMRTriangle*[ 2*max_num_points ];
    memcpy(new_pts_to_tris, pts_to_tris, num_points*sizeof(TMRTriangle*));
    delete [] pts_to_tris;
    pts_to_tris = new_pts_to_tris;
    
    if (face){
      // Allocate the space for the triangle vertices in physical space
      TMRPoint *new_X = new TMRPoint[ max_num_points ];
      memcpy(new_X, X, num_points*sizeof(TMRPoint));
      delete [] X;
      X = new_X;
    }
  }

  // Add the point to the quadtree
  root->addNode(num_points, pt);

  // Set the new point location
  pts[2*num_points] = pt[0];
  pts[2*num_points+1] = pt[1];

  // No new triangle has been assigned yet
  pts_to_tris[num_points] = NULL;

  // Evaluate the face location
  if (face){
    face->evalPoint(pt[0], pt[1], &X[num_points]);
  }
  
  // Increase the number of points by one
  num_points++;

  // Return the new point index
  return num_points-1;
}

/*
  Add the vertex to the underlying Delaunay triangularization.
*/
void TMRTriangularize::addPointToMesh( const double pt[] ){
  // Find the enclosing triangle
  TMRTriangle *tri;
  findEnclosing(pt, &tri);

  // Add the point to the quadtree
  uint32_t u = addPoint(pt);

  if (tri){
    uint32_t v = tri->u;
    uint32_t w = tri->v;
    uint32_t x = tri->w;
    deleteTriangle(*tri);
    digCavity(u, v, w);
    digCavity(u, w, x);
    digCavity(u, x, v);
  }

  return;
}

/*
  Add the point to the mesh, given that we've already found the 
  triangle that encloses the given point. 
*/
void TMRTriangularize::addPointToMesh( const double pt[], 
                                       TMRTriangle *tri ){
  uint32_t u = addPoint(pt);

  if (tri){
    uint32_t v = tri->u;
    uint32_t w = tri->v;
    uint32_t x = tri->w;
    deleteTriangle(*tri);
    digCavity(u, v, w);
    digCavity(u, w, x);
    digCavity(u, x, v);
  }
}

/*
  The following code tests whether the triangle formed from the point
  (u, v, w) is constrained Delaunay. 

  If the edge (w, v) is in the PSLG, then the triangle is added immediately. 
  If not, and if the point lies within the circumcircle of the triangle 
  (u, w, v) with the directed edge (v, w).
*/
void TMRTriangularize::digCavity( uint32_t u, uint32_t v, uint32_t w ){ 
  // If the edge is along the polynomial straight line graph, then we
  // add the triangle as it exists and we're done, even though it may
  // not be Delaunay (it will be constrained Delaunay). We cannot
  // split a PSLG edge.
  if (edgeInPSLG(w, v)){
    addTriangle(TMRTriangle(u, v, w));
    return; 
  }

  // Complete the triangle
  TMRTriangle *tri;
  completeMe(w, v, &tri);

  if (tri){
    // Get the index of the final vertex
    uint32_t x = 0;
    if (tri->u == w && tri->v == v){
      x = tri->w;
    }
    else if (tri->v == w && tri->w == v){
      x = tri->u;
    }
    else if (tri->w == w && tri->u == v){
      x = tri->v;
    }

    // Check whether the point lies within the circumcircle or
    // exactly on its boundary
    if (inCircle(u, v, w, x) > 0.0){
      deleteTriangle(*tri);
      digCavity(u, v, x);
      digCavity(u, x, w);
      return;
    }
  }

  addTriangle(TMRTriangle(u, v, w));
}

/*
  Find the enclosing triangle within the mesh. 

  This code uses a quadtree for geometric searching. First, we find
  the node that is closest to the query point. This node is not
  necessarily connected with the enclosing triangle that we
  want. Next, we find one triangle associated with this node. If this
  triangle does not contain the point, we march over the mesh, marking
  the triangles that we have visited.
*/
void TMRTriangularize::findEnclosing( const double pt[],
                                      TMRTriangle **ptr ){
  *ptr = NULL;
  if (search_tag == UINT_MAX){
    search_tag = 0;
    setTriangleTags(0);
  }
  search_tag++;

  // Find the closest point to the given 
  uint32_t u = root->findClosest(pt, NULL);

  // Obtain the triangle associated with the node u. This triangle may
  // not contain the node, but will hopefully be close to the node.
  // We'll walk the mesh to nearby elements until we find the proper
  // enclosing triangle.
  TMRTriangle *tri = pts_to_tris[u];

  if (enclosed(pt, tri->u, tri->v, tri->w)){
    *ptr = tri;
    return;
  }
  else {
    tri->tag = search_tag;
  }

  // Members of the list do not enclose the point and have
  // been labeled that they are searched
  TriQueue queue;
  queue.append(tri);

  while (queue.size > 0){
    // Pop the top member from the queue
    TMRTriangle *t = queue.pop();

    // Search the adjacent triangles across edges
    uint32_t edge_pairs[][2] = {{t->u, t->v},
                                {t->v, t->w},
                                {t->w, t->u}};

    // Search the adjacent triangles and determine whether they have
    // been tagged
    for ( int k = 0; k < 3; k++ ){
      TMRTriangle *t2;
      completeMe(edge_pairs[k][1], edge_pairs[k][0], &t2);
      if (t2 && t2->tag != search_tag){
        // Check whether the point is enclosed by t2
        if (enclosed(pt, t2->u, t2->v, t2->w)){
          *ptr = t2;
          return;
        }

        // If the triangle does not enclose the point, label this guy
        // as having been searched
        t2->tag = search_tag;

        // Append this guy to the queue
        queue.append(t2);
      }
    }
  }

  return;
}

/*
  Compute the distance to the intersection of the line m + alpha*e
  with the two lines from u to w and v to w

  This code is used to ensure that a point will lie within the current
  triangle if we've failed to pick a point that lies within the
  domain.
*/
double TMRTriangularize::computeIntersection( const double m[], 
                                              const double e[], 
                                              uint32_t u, uint32_t v, uint32_t w ){
  // m + alpha*e = pts[u] + beta*(pts[w] - pts[u])
  // => alpha*e + beta*(pts[u] - pts[w]) = pts[u] - m
  double a11 = e[0], a21 = e[1];
  double a12 = pts[2*u] - pts[2*w];
  double a22 = pts[2*u+1] - pts[2*w+1];
  double det = a11*a22 - a12*a21;

  // Check if the line is orthogonal to this direction
  if (fabs(det) > 1e-12*fabs(a11*a22)){
    double b1 = pts[2*u] - m[0];
    double b2 = pts[2*u+1] - m[1];
    double beta = (a11*b2 - a21*b1)/det;

    if (beta >= 0.0 && beta <= 1.0){
      double alpha = (a22*b1 - a12*b2)/det;
      return alpha;
    }
  }

  // If the line from u to w was orthogonal to the e direction in the
  // previous case, it cannot be orthogonal in this case and still
  // form a triangle. No check is performed here, but perhaps it
  // should be.
  a12 = pts[2*v] - pts[2*w];
  a22 = pts[2*v+1] - pts[2*w+1];
  det = a11*a22 - a12*a21;

  if (fabs(det) > 1e-12*fabs(a11*a22)){
    double b1 = pts[2*v] - m[0];
    double b2 = pts[2*v+1] - m[1];
    double beta = (a11*b2 - a21*b1)/det;

    if (beta >= 0.0 && beta <= 1.0){
      double alpha = (a22*b1 - a12*b2)/det;
      return alpha;
    }
  }
 
  // If we get here, something really bad has happened. (u, v, w) do
  // not form a triangle?
  return -1.0;
}

/*
  Compute the maximum edge length between any two points.

  This edge length is used as a metric for quality, although having
  the 'right' edge length doesn't imply that we've found a
  high-quality triangle.
*/
double TMRTriangularize::computeMaxEdgeLength( TMRTriangle *tri ){
  double D1, D2, D3;
  if (face){
    // Compute the edge lengths in physical space
    TMRPoint d;
    d.x = X[tri->v].x - X[tri->u].x;
    d.y = X[tri->v].y - X[tri->u].y;
    d.z = X[tri->v].z - X[tri->u].z;
    D1 = d.dot(d);

    d.x = X[tri->w].x - X[tri->v].x;
    d.y = X[tri->w].y - X[tri->v].y;
    d.z = X[tri->w].z - X[tri->v].z;
    D2 = d.dot(d);

    d.x = X[tri->u].x - X[tri->w].x;
    d.y = X[tri->u].y - X[tri->w].y;
    d.z = X[tri->u].z - X[tri->w].z;
    D3 = d.dot(d);
  }
  else {
    // Compute the edge lengths in parameter space
    double d1[2];
    d1[0] = pts[2*tri->v] - pts[2*tri->u];
    d1[1] = pts[2*tri->v+1] - pts[2*tri->u+1];
    D1 = d1[0]*d1[0] + d1[1]*d1[1];
    
    double d2[2];
    d2[0] = pts[2*tri->w] - pts[2*tri->v];
    d2[1] = pts[2*tri->w+1] - pts[2*tri->v+1];
    D2 = d2[0]*d2[0] + d2[1]*d2[1];
    
    double d3[2];
    d3[0] = pts[2*tri->u] - pts[2*tri->w];
    d3[1] = pts[2*tri->u+1] - pts[2*tri->w+1];
    D3 = d3[0]*d3[0] + d3[1]*d3[1];
  }

  double Dmax = D1;
  if (D2 > Dmax){ Dmax = D2; }
  if (D3 > Dmax){ Dmax = D3; }

  return sqrt(Dmax);
}

/*
  Compute the circumcircle for the given triangle.

  This can be used to evaluate an effective 'h' value (based on the 
  equilateral radius length r_eq = h/sqrt(3)) which  is used as a
  metric to determine whether to retain the triangle, or search for
  a better one.
*/
double TMRTriangularize::computeCircumcircle( TMRTriangle *tri ){
  double R = 0.0;
  if (face){
    // Compute the edge lengths in physical space
    TMRPoint d1;
    d1.x = X[tri->v].x - X[tri->u].x;
    d1.y = X[tri->v].y - X[tri->u].y;
    d1.z = X[tri->v].z - X[tri->u].z;
    
    TMRPoint d2;
    d2.x = X[tri->w].x - X[tri->u].x;
    d2.y = X[tri->w].y - X[tri->u].y;
    d2.z = X[tri->w].z - X[tri->u].z;
 
    double dot = d1.dot(d2)/d1.dot(d1);

    // Compute the perpendicular component along the second
    // direction that we'll use to determine the center point
    TMRPoint n1;
    n1.x = d2.x - dot*d1.x;
    n1.y = d2.y - dot*d1.y;
    n1.z = d2.z - dot*d1.z;

    // Compute alpha = 0.5*(d2, d2 - d1)/(d2, n1)
    double alpha = 0.5*(d2.x*(d2.x - d1.x) + 
                        d2.y*(d2.y - d1.y) + 
                        d2.z*(d2.z - d1.z));
    alpha = alpha/d2.dot(n1);

    // Compute the distance from the point
    d1.x = 0.5*d1.x + alpha*n1.x;
    d1.y = 0.5*d1.y + alpha*n1.y;
    d1.z = 0.5*d1.z + alpha*n1.z;

    // Compute the radius
    R = sqrt(d1.dot(d1));
  }
  else {
    double d1[2];
    d1[0] = pts[2*tri->v] - pts[2*tri->u];
    d1[1] = pts[2*tri->v+1] - pts[2*tri->u+1];
    
    double d2[2];
    d2[0] = pts[2*tri->w] - pts[2*tri->v];
    d2[1] = pts[2*tri->w+1] - pts[2*tri->v+1];

    double b[2];
    b[0] = 0.5*(pts[2*tri->w] - pts[2*tri->v]);
    b[1] = 0.5*(pts[2*tri->w+1] - pts[2*tri->v+1]);

    // Evaluate alpha
    double invdet = 1.0/(d1[1]*d2[0] - d2[1]*d1[0]);
    double alpha = -invdet*(d2[0]*b[0] + d2[1]*b[1]);

    // Compute the point locations
    double pt[2];
    pt[0] = 0.5*(pts[2*tri->u] + pts[2*tri->v]) - alpha*d1[1];
    pt[1] = 0.5*(pts[2*tri->u+1] + pts[2*tri->v+1]) + alpha*d1[0];

    // Determine the radius
    d1[0] = pts[2*tri->u] - pt[0];
    d1[1] = pts[2*tri->u+1] - pt[1];
    R = sqrt(d1[0]*d1[0] + d1[1]*d1[1]);
  }

  return R;
}

/*
  Perform a frontal mesh generation algorithm and simultaneously
  create a constrained Delaunay triangularization of the generated
  mesh.

  The Delaunay triangularization based on the Bowyer-Watson mesh
  generation algorithm. The frontal mesh generation technique is based
  on Rebay's 1993 paper in JCP.
*/
void TMRTriangularize::frontal( double h ){
  // The length of the edge of the triangle
  const double sqrt3 = sqrt(3.0);
  const double de = 0.5*sqrt3*h;

  // Add the triangles to the active set that 
  TriListNode *node = list_start;
  while (node){
    // Set the status by default as waiting
    node->tri.status = WAITING;
    
    // Compute the 'quality' indicator for this triangle
    double hval = sqrt3*computeCircumcircle(&node->tri);
    node->tri.quality = hval;
    if (hval < 1.25*h){
      node->tri.status = ACCEPTED;
    }
    else {
      // If any of the triangles touches an edge in the planar
      // straight line graph, change it to a waiting triangle
      uint32_t edge_pairs[][2] = {{node->tri.u, node->tri.v},
                                  {node->tri.v, node->tri.w},
                                  {node->tri.w, node->tri.u}};
      for ( int k = 0; k < 3; k++ ){
        if (edgeInPSLG(edge_pairs[k][0], edge_pairs[k][1])){
          node->tri.status = ACTIVE;
          addActiveTriangle(&node->tri);
          break;
        }
      }
    }

    // Go to the next triangle in the list
    node = node->next;
  }

  // Iterate over the list again and add any triangles that are
  // adjacent to an ACCEPTED triangle to the ACTIVE set of triangles
  node = list_start;
  while (node){
    if (node->tri.status == ACCEPTED){
      // Check if any of the adjacent triangles are WAITING.  If so,
      // change their status to ACTIVE
      uint32_t edge_pairs[][2] = {{node->tri.u, node->tri.v},
                                  {node->tri.v, node->tri.w},
                                  {node->tri.w, node->tri.u}};

      for ( int k = 0; k < 3; k++ ){
        TMRTriangle *adjacent;
        completeMe(edge_pairs[k][1], edge_pairs[k][0], &adjacent);
        if (adjacent && adjacent->status == WAITING){
          node->tri.status = ACTIVE;
          addActiveTriangle(&node->tri);
          break;
        }
      }
    }

    // Increment the pointer to the next member of the list
    node = node->next;
  }

  printf("%10s %10s %10s\n", "Iteration", "Triangles", "Active");
 
  int iter = 0;
  while (1){
    if (iter % 1000 == 0){
      printf("%10d %10d %10d\n",
             iter, num_triangles, num_active_triangles);
    }
    iter++;

    // The pointer to the triangle that we're going to use next
    TMRTriangle *tri = NULL;

    // Search through all active nodes
    for ( int i = 0; i < num_active_buckets; i++ ){
      ActiveHashNode *act = active_buckets[i];
      while (act){
        if (!tri){
          tri = act->tri;
        }
        else if (act->tri->quality > tri->quality){
          tri = act->tri;
        }

        // Increment the active node pointer
        act = act->next;
      }
    }

    // We've failed to find any triangle in the active set of
    // triangles
    if (!tri){
      break;
    }   

    // Now insert a new point based on the Voronoi criteria.
    int found = 0;
    uint32_t u = 0, v = 0;
    uint32_t edge_pairs[][2] = {{tri->u, tri->v},
                                {tri->v, tri->w},
                                {tri->w, tri->u}};
    
    // Check if the edge is on the PSLG
    for ( int k = 0; k < 3; k++ ){
      u = edge_pairs[k][0];
      v = edge_pairs[k][1];
      if (edgeInPSLG(u, v)){
        found = 1;
        break;
      }
    }

    // Check if an adjacent triangle is accepted
    if (!found){
      for ( int k = 0; k < 3; k++ ){
        u = edge_pairs[k][0];
        v = edge_pairs[k][1];

        // Compute the completed triangle
        TMRTriangle *t;
        completeMe(v, u, &t);
        if (t && t->status == ACCEPTED){
          found = 1;
          break;
        }
      }
    }

    // Compute the location of the new point
    // | i   j   k |
    // | 0   0   1 | = - i*dy + j*dx
    // | dx  dy  0 |

    // Compute the parametric mid-point
    double m[2], e[2];
    m[0] = 0.5*(pts[2*u] + pts[2*v]);
    m[1] = 0.5*(pts[2*u+1] + pts[2*v+1]);

    // Half the parametric distance between the two edge points
    double p = 0.0;

    double pt[2];
    if (face){
      TMRPoint Xu, Xv;
      face->evalDeriv(m[0], m[1], &Xu, &Xv);

      // Compute the metric tensor components
      double g11 = Xu.dot(Xu);
      double g12 = Xu.dot(Xv);
      double g22 = Xv.dot(Xv);

      // Compute the inverse metric components
      double invdet = 1.0/(g11*g22 - g12*g12);
      double G11 = invdet*g22;
      double G12 = -invdet*g12;
      double G22 = invdet*g11;

      // Compute the parametric direction along the curvilinear line
      // connecting from u to v
      double d[2];
      d[0] = (pts[2*v] - pts[2*u]);
      d[1] = (pts[2*v+1] - pts[2*u+1]);

      // Compute the orthogonal coordinate contributions to the vector
      e[0] = G12*d[0] - G11*d[1];
      e[1] = G22*d[0] - G12*d[1];

      // Compute the direction in physical space
      TMRPoint dir;
      dir.x = e[0]*Xu.x + e[1]*Xv.x;
      dir.y = e[0]*Xu.y + e[1]*Xv.y;
      dir.z = e[0]*Xu.z + e[1]*Xv.z;

      // Compute the ratio between the desired distance in physical
      // space and the length of the direction dir in physical space
      // for a unit change in parameter space.
      double f = de/sqrt(dir.dot(dir));
      pt[0] = m[0] + f*e[0];
      pt[1] = m[1] + f*e[1];
    }
    else {
      // Compute the direction perpendicular to the line from u to v
      e[0] = (pts[2*u+1] - pts[2*v+1]);
      e[1] = (pts[2*v] - pts[2*u]);
      
      // Compute half the distance between the u and v points
      p = 0.5*sqrt(e[0]*e[0] + e[1]*e[1]);
      e[0] = 0.5*e[0]/p;
      e[1] = 0.5*e[1]/p;
      
      // Compute the new optimal point
      pt[0] = m[0] + de*e[0];
      pt[1] = m[1] + de*e[1];
    }

    // Find the enclosing triangle for the 
    TMRTriangle *pt_tri = tri;
    if (!enclosed(pt, pt_tri->u, pt_tri->v, pt_tri->w)){
      findEnclosing(pt, &pt_tri);
      
      // We've tried a new point and it was outside the domain.  That
      // is not allowed, so next we pick a point that is actually in
      // the current triangle and therefore within the allowed domain
      if (!pt_tri){
        pt_tri = tri;

        if (face){
          // Pick a point at the middle of the current triangle
          const double frac = 1.0/3.0;
          pt[0] = frac*(pts[2*pt_tri->u] + pts[2*pt_tri->v] + pts[2*pt_tri->w]);
          pt[1] = frac*(pts[2*pt_tri->u+1] + pts[2*pt_tri->v+1] + pts[2*pt_tri->w+1]);
        }
        else {
          // Pick a point in the triangle
          uint32_t w = 0;
          if (u == pt_tri->u && v == pt_tri->v){
            w = pt_tri->w;
          }
          else if (u == pt_tri->v && v == pt_tri->w){
            w = pt_tri->u;
          }
          else if (u == pt_tri->w && v == pt_tri->u){
            w = pt_tri->v;
          }
          
          double q = computeIntersection(m, e, u, v, w);
          double rho = (p*p + q*q)/q;
          pt[0] = m[0] + rho*e[0];
          pt[1] = m[1] + rho*e[1];
        }
      }
    }
  
    // Set the pointer to the last member in the list
    list_marker = list_end; 
    addPointToMesh(pt, pt_tri);
    pt_tri = NULL;

    // Compute the circumcircles of the new triangles and check
    // whether they belong in the done category or not...
    TriListNode *ptr = list_start;
    if (list_marker){
      ptr = list_marker->next;
    }
    while (ptr){
      double hval = sqrt3*computeCircumcircle(&ptr->tri);
      ptr->tri.quality = hval;
      if (hval < 1.25*h){
        ptr->tri.status = ACCEPTED;
      }
      else {
        ptr->tri.status = WAITING;
      }
      ptr = ptr->next;
    }

    // Complete me with the newly created triangle with the
    completeMe(u, v, &pt_tri);
    if (pt_tri){
      if (pt_tri->status == ACTIVE){
        deleteActiveTriangle(pt_tri);
      }
      pt_tri->status = ACCEPTED;
    }

    // Scan through the list of the added triangles and mark which
    // ones are active/working/accepted.
    ptr = list_start;
    if (list_marker){ 
      ptr = list_marker->next;
    }

    while (ptr){
      if (ptr->tri.status != ACCEPTED){
        // If any of the triangles touches an edge in the planar
        // straight line graph, change it to a waiting triangle
        int flag = 0;
        uint32_t edge_pairs[][2] = {{ptr->tri.u, ptr->tri.v},
                                    {ptr->tri.v, ptr->tri.w},
                                    {ptr->tri.w, ptr->tri.u}};

        // Loop over all of the edges in the triangle and check
        // whether they're in the PSLG 
        for ( int k = 0; k < 3; k++ ){
          if (edgeInPSLG(edge_pairs[k][0], edge_pairs[k][1])){
            ptr->tri.status = ACTIVE;
            addActiveTriangle(&ptr->tri);
            flag = 1;
            break;
          }
        }

        // If any triangle is adjacent to a triangle that is accepted,
        // then change the status of the new triangle to be active
        if (!flag){
          for ( int k = 0; k < 3; k++ ){
            TMRTriangle *adjacent;
            completeMe(edge_pairs[k][1], edge_pairs[k][0], &adjacent);
            if (adjacent && adjacent->status == ACCEPTED){
              ptr->tri.status = ACTIVE;
              addActiveTriangle(&ptr->tri);
              break;
            }
          }
        }
      }

      // Increment the pointer to the next member of the list
      ptr = ptr->next;
    }
  }

  printf("%10d %10d %10d\n",
         iter, num_triangles, num_active_triangles);
}
