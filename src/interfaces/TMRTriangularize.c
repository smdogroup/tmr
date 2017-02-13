#include "TMRTriangularize.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>

#include "predicates.h"

// Define the maximum distance
const double MAX_QUAD_DISTANCE = 1e40;

/*
  Compare coordinate pairs of points. This uses a Morton ordering
  comparison to facilitate sorting/searching the list of edges.

  This function can be used by the stdlib functions qsort and
  bsearch to sort/search the edge pairs.
*/
int compare_edges( const void *avoid, const void *bvoid ){
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
  
*/
TMRTriangularize::TMRTriangularize( int npts, const double inpts[], 
                                    int nsegs, const int segs[] ){
  // Initialize the predicates code
  exactinit();

  // Allocate and initialize/zero the hash table data for the edges
  num_hash_nodes = 0;
  num_buckets = 100;
  buckets = new EdgeHashNode*[ num_buckets ];
  memset(buckets, 0, num_buckets*sizeof(EdgeHashNode*));

  // Set the initial root/current location of the doubly-linked
  // triangle list structure. These are allocated as we add new triangles.
  list_start = NULL;
  list_end = NULL;
  list_marker = NULL;

  // Keep track of the total number of triangles
  num_triangles = 0;

  // Set the initial number of points
  fixed_point_offset = 4;
  num_points = fixed_point_offset + npts;

  // Set the initial for the maximum number of points
  max_num_points = 1024;
  if (max_num_points < num_points){
    max_num_points = num_points;
  }

  // Allocate the initial set of points
  pts = new double[ 2*max_num_points ];
  pts_to_tris = new TMRTriangle*[ max_num_points ];

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

  // Re-adjust the domain boundary to ensure that it is
  // sufficiently large
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
  num_points = 4;

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

  // Add the extreme points to the quadtree
  for ( int i = 0; i < 4; i++ ){
    root->addNode(i, &pts[2*i]);
  }

  // Add the initial triangles
  addTriangle(TMRTriangle(0, 1, 2));
  addTriangle(TMRTriangle(2, 1, 3));

  // Add the points to the triangle
  for ( int i = 0; i < npts; i++ ){
    addPointToMesh(&inpts[2*i]);
  }

  // Set the triangle tags to zero
  setTriangleTags(0);

  // Mark all the triangles in the list that contain or touch nodes
  // that are in the fixed_point_offset list that are not separated by
  // a PSLG edge. These triangles will be deleted.
  int flag = 1;
  while (flag){
    flag = 0;

    TriListNode *node = list_start;
    while (node){
      if (node->tri.tag == 0 && 
          (node->tri.u < fixed_point_offset ||
           node->tri.v < fixed_point_offset ||
           node->tri.w < fixed_point_offset)){
        tagTriangles(&node->tri);
      }
      node = node->next;
    }
  }

  // Free the triangles that we've found
  list_marker = NULL;
  TriListNode *node = list_start;
  while (node){
    // Keep track of what node will be next, since we may be deleting
    // 'node' itself, we cannot access this afterwards
    TriListNode *tmp = node->next;

    // Delete the triangle from the triangle list
    if (node->tri.tag){
      deleteTriangle(node->tri);
    }

    // Set the next node
    node = tmp;
  }

  // Free the points and holes from the quadtree
  for ( uint32_t num = 0; num < fixed_point_offset; num++ ){
    root->deleteNode(num, &pts[2*num]);
  }
}

/*
  Free the triangularization object
*/
TMRTriangularize::~TMRTriangularize(){
  delete root;
  delete [] pts;

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
void TMRTriangularize::getMesh( int *_num_triangles, 
                                int **_conn, double **_pts ){
  *_num_triangles = num_triangles;
  *_conn = new int[ 3*num_triangles ];

  int npts = (num_points - fixed_point_offset);
  *_pts = new double[ 2*npts ];

  // Set the points
  memcpy(*_pts, &pts[2*fixed_point_offset], 2*npts*sizeof(double));

  // Set the pointer into the connectivity array
  int *t = *_conn;

  // Determine the connectivity
  TriListNode *node = list_start;
  while (node){
    t[0] = node->tri.u - fixed_point_offset;
    t[1] = node->tri.v - fixed_point_offset;
    t[2] = node->tri.w - fixed_point_offset;
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
  Write out the
*/
void TMRTriangularize::writeToVTK( const char *filename ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
    // Write out the points
    fprintf(fp, "POINTS %d float\n", num_points);
    for ( int k = 0; k < num_points; k++ ){
      fprintf(fp, "%e %e %e\n", pts[2*k], pts[2*k+1], 0.0);
    }
    
    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", 
            num_triangles, 4*num_triangles);

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
  return (x * 0x1f1f1f1f) ^ y;
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

  // Set the search tag to zero
  list_end->tri.status = 0;
  list_end->tri.tag = 0;

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
      u = segs[2*i] + fixed_point_offset;
    }
    if (segs[2*i+1] >= 0){
      v = segs[2*i+1] + fixed_point_offset;
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

  return (bsearch(edge, pslg_edges, num_pslg_edges, 2*sizeof(uint32_t),
                  compare_edges) != NULL);
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
  return incircle(&pts[2*u], &pts[2*v], &pts[2*w], &pts[2*x]);
}

/*
  Add a point to the point list
*/
uint32_t TMRTriangularize::addPoint( const double pt[] ){
  // If the size of the array is exceeded, multiply it by 2
  if (num_points >= max_num_points){
    max_num_points *= 2;
    double *new_pts = new double[ 2*max_num_points ];
    memcpy(new_pts, pts, 2*num_points*sizeof(double));
    delete [] pts;
    pts = new_pts;

    TMRTriangle **new_pts_to_tris = new TMRTriangle*[ 2*max_num_points ];
    memcpy(new_pts_to_tris, pts_to_tris, num_points*sizeof(TMRTriangle*));
    delete [] pts_to_tris;
    pts_to_tris = new_pts_to_tris;
  }

  // Add the point to the quadtree
  root->addNode(num_points, pt);

  // Set the new point location
  pts[2*num_points] = pt[0];
  pts[2*num_points+1] = pt[1];
  pts_to_tris[num_points] = NULL;
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

  // Scan through the mesh, using a linear search
  *ptr = NULL;
  TriListNode *node = list_start;
  while (node){
    if (enclosed(pt, node->tri.u, node->tri.v, node->tri.w)){
      *ptr = &(node->tri);
      return;
    }
    node = node->next;
  }
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

  return -1.0;
}

/*
  Compute the maximum edge length between any two points.

  This edge length is used as a metric for quality, although having
  the 'right' edge length doesn't imply that we've found a
  high-quality triangle.
*/
double TMRTriangularize::computeMaxEdgeLength( TMRTriangle *tri ){
  double d1[2], D1;
  d1[0] = pts[2*tri->v] - pts[2*tri->u];
  d1[1] = pts[2*tri->v+1] - pts[2*tri->u+1];
  D1 = d1[0]*d1[0] + d1[1]*d1[1];

  double d2[2], D2;
  d2[0] = pts[2*tri->w] - pts[2*tri->v];
  d2[1] = pts[2*tri->w+1] - pts[2*tri->v+1];
  D2 = d2[0]*d2[0] + d2[1]*d2[1];

  double d3[2], D3;
  d3[0] = pts[2*tri->u] - pts[2*tri->w];
  d3[1] = pts[2*tri->u+1] - pts[2*tri->w+1];
  D3 = d3[0]*d3[0] + d3[1]*d3[1];

  double Dmax = D1;
  if (D2 > Dmax){ Dmax = D2; }
  if (D3 > Dmax){ Dmax = D3; }

  return sqrt(Dmax);
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
  double de = 0.5*sqrt(3.0)*h;

  // Set the status of the different 
  const uint32_t WAITING = 1;
  const uint32_t ACTIVE = 2;
  const uint32_t ACCEPTED = 3;

  // Add the triangles to the active set that 
  TriListNode *node = list_start;
  while (node){
    node->tri.status = WAITING;
  
    // If any of the triangles touches an edge in the planar
    // straight line graph, change it to a waiting triangle
    uint32_t edge_pairs[][2] = {{node->tri.u, node->tri.v},
                                {node->tri.v, node->tri.w},
                                {node->tri.w, node->tri.u}};
    for ( int k = 0; k < 3; k++ ){
      if (edgeInPSLG(edge_pairs[k][0], edge_pairs[k][1])){
        node->tri.status = ACTIVE;
        break;
      }
    }

    // Go to the next triangle in the list
    node = node->next;
  }
  
  int iter = 0;
  while (1){
    if (iter % 250 == 0){
      printf("Iteration %d  num_triangles %d\n", iter, num_triangles);
    }
    iter++;

    // Find the first triangle in the list that is tagged as being
    // active. This search could be made more efficient by storing the
    // active triangles separately.
    node = NULL;
    TriListNode *tmp = list_start;
    double hmax = 0.0;
    while (tmp){
      if (tmp->tri.status == ACTIVE){
        if (tmp->tri.quality > hmax){
          hmax = tmp->tri.quality;
          node = tmp;
        }
      }
      tmp = tmp->next;
    }
    if (!node){
      break;
    }

    // Now insert a new point based on the Voronoi criteria.
    int found = 0;
    uint32_t u = 0, v = 0;
    uint32_t edge_pairs[][2] = {{node->tri.u, node->tri.v},
                                {node->tri.v, node->tri.w},
                                {node->tri.w, node->tri.u}};
    
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
        TMRTriangle *tri;
        completeMe(v, u, &tri);
        if (tri && tri->status == ACCEPTED){
          found = 1;
          break;
        }
      }
    }

    // Compute the location of the new point
    // | i   j   k |
    // | 0   0   1 | = - i*dy + j*dx
    // | dx  dy  0 |

    // Compute
    double m[2];
    m[0] = 0.5*(pts[2*u] + pts[2*v]);
    m[1] = 0.5*(pts[2*u+1] + pts[2*v+1]);

    // Compute the direction
    double e[2];
    e[0] = (pts[2*u+1] - pts[2*v+1]);
    e[1] = (pts[2*v] - pts[2*u]);

    // Compute half the distance between the u and v points
    double p = 0.5*sqrt(e[0]*e[0] + e[1]*e[1]);
    e[0] = 0.5*e[0]/p;
    e[1] = 0.5*e[1]/p;

    // Compute the new optimal point
    double pt[2];
    pt[0] = m[0] + de*e[0];
    pt[1] = m[1] + de*e[1];

    // Find the enclosing triangle for the 
    TMRTriangle *pt_tri = &(node->tri);
    if (!enclosed(pt, pt_tri->u, pt_tri->v, pt_tri->w)){
      findEnclosing(pt, &pt_tri);

      // We've tried a new point and it was outside the domain.
      // That is not allowed, so next we pick a point that is 
      // actually in the current triangle and therefore within the
      // allowed domain
      if (!pt_tri){
        pt_tri = &(node->tri);
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
        double rho = 0.5*(p*p + q*q)/q;
        pt[0] = m[0] + rho*e[0];
        pt[1] = m[1] + rho*e[1];
      }
    }

    // Set the pointer to the last member in the list
    list_marker = list_end; 
    addPointToMesh(pt, pt_tri);
    pt_tri = NULL;

    // Complete me with the newly created triangle with the
    completeMe(u, v, &pt_tri);
    if (pt_tri){
      pt_tri->status = ACCEPTED;
    }

    // Compute the circumcircles of the new triangles and check
    // whether they belong in the done category or not...
    TriListNode *ptr = list_start;
    int count = 0;
    if (list_marker){
      ptr = list_marker->next;
    }
    while (ptr){
      double hval = computeMaxEdgeLength(&ptr->tri);
      ptr->tri.quality = hval;
      if (hval < 1.5*h){
        ptr->tri.status = ACCEPTED;
      }
      else {
        ptr->tri.status = WAITING;
      }
      ptr = ptr->next;
      count++;
    }

    // Scan through the list of the 
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
              break;
            }
          }
        }
      }

      // Increment the pointer to the next member of the list
      ptr = ptr->next;
    }
  }
}
