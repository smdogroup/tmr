#include "TMRTriangularize.h"
#include <math.h>
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
  needs to be used to scan through the quadtree hierarchy to find
  the correct quadrant location.
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

  This scans each quadrant, testing whether it will be necessary to scan
  multiple quadrants to ensure that we have, in fact, located the closest
  point to the provided point.
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
  list_root = NULL;
  list_current = NULL;

  // Keep track of the number of triangles
  num_triangles = 0;

  // Set the initial number of points
  num_points = 4+npts;

  // Set the initial for the maximum number of points
  max_num_points = 1024;
  if (max_num_points < num_points){
    max_num_points = num_points;
  }

  // Allocate the initial set of points
  pts = new double[ 2*max_num_points ];

  // Find the maximum domain size
  domain.xlow = domain.xhigh = pts[0];
  domain.ylow = domain.yhigh = pts[1];
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
    printf("Adding point %d\n", i);
    addPointToMesh(&inpts[2*i]);
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
  while (list_root){
    TriListNode *tmp = list_root;
    list_root = list_root->next;
    delete tmp;
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

    TriListNode *node = list_root;
    int count = 0;
    while (node){
      count++;
      fprintf(fp, "3 %d %d %d\n", node->tri.u, node->tri.v, node->tri.w);
      node = node->next;
    }
    printf("count = %d  num_triangles = %d\n", count, num_triangles);

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", num_triangles);
    for ( int k = 0; k < num_triangles; k++ ){
      fprintf(fp, "%d\n", 5);
    }

    //     // Print out the rest as fields one-by-one
    // fprintf(fp, "CELL_DATA %d\n", num_triangles);
    // fprintf(fp, "SCALARS quality float 1\n");
    // fprintf(fp, "LOOKUP_TABLE default\n");
    // for ( int k = 0; k < ntris; k++ ){
    //   fprintf(fp, "%e\n", computeTriQuality(&tris[3*k], pts));
    // }

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
  if (!list_root){
    // Create the root triangle node
    list_root = new TriListNode();
    list_root->tri = tri;
    list_root->next = NULL;
    list_root->prev = NULL;

    // Update the current node
    list_current = list_root;
  }
  else {
    // Create the new member in the triangle list
    TriListNode *next = new TriListNode();
    next->tri = tri;
    next->next = NULL;
    next->prev = list_current;

    // Update the current member in the triangle list
    list_current->next = next;
    list_current = list_current->next;
  }

  // Redistribute the members in the hash table if required
  if (num_hash_nodes > 10*num_buckets){
    printf("Rebucket hash tables\n");
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
        int value = getEdgeHash(node->u, node->v);
        int bucket = value % num_buckets;
        
        // If the new bucket linked list does not exist
        if (!new_buckets[bucket]){
          new_buckets[bucket] = node;
          end_buckets[bucket] = new_buckets[bucket];
        }
        else {
          end_buckets[bucket]->next = node;
          end_buckets[bucket] = end_buckets[bucket]->next;
          node->next = NULL;
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
  uint32_t edge_pairs[] = {tri.u, tri.v,
                           tri.v, tri.w,
                           tri.w, tri.u};
  
  // Add the triangle to the hash table
  for ( int k = 0; k < 3; k++ ){
    // Add a hash for each pair of edges around the triangle
    uint32_t u = edge_pairs[2*k];
    uint32_t v = edge_pairs[2*k+1];
    int value = getEdgeHash(u, v);
    int bucket = value % num_buckets;
    if (!buckets[bucket]){
      // Create the new buckets node and assign the values
      buckets[bucket] = new EdgeHashNode();
      buckets[bucket]->u = u;
      buckets[bucket]->v = v;
      buckets[bucket]->tri_node = list_current;
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
          node->tri_node = list_current;

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
        node->tri_node = list_current;
        node->next = NULL;
        num_hash_nodes++;
      }
      else {
        success = 0;
      }
    }
  }

  if (success){
    num_triangles++;
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

  printf("Deleting triangle (%d, %d, %d)\n", tri.u, tri.v, tri.w);

  // Set the combinations of edge pairs that will be added
  // to the hash table
  uint32_t edge_pairs[] = {tri.u, tri.v,
                          tri.v, tri.w,
                          tri.w, tri.u};

  // Keep track if this is the first edge we find
  int first = 1;
  
  // Add the triangle to the hash table
  for ( int k = 0; k < 3; k++ ){
    // Add a hash for each pair of edges around the triangle
    uint32_t u = edge_pairs[2*k];
    uint32_t v = edge_pairs[2*k+1];

    // Keep track of whether we've deleted this edge
    int edge_success = 0;

    // Get the hash values and access the first entry of the bucket
    // it's listed under
    int value = getEdgeHash(u, v);
    int bucket = value % num_buckets;
    EdgeHashNode *node = buckets[bucket];
    EdgeHashNode *prev = node;
    while (node){
      // The edge matches, we have to delete this triangle
      if (u == node->u && v == node->v){
        // If this is the first edge we've found, delete the
        // corresponding list entry as well.
        if (first){
          // Adjust the double linking so that when the node is
          // deleted, the pointers will still work
          TriListNode *ptr = node->tri_node;
          if (ptr == list_root){
            list_root = list_root->next;
            if (list_root){
              list_root->prev = NULL;
            }
          }
          else if (ptr == list_current){
            list_current = list_current->prev;
            if (list_current){
              list_current->next = NULL;
            }
          }
          else {
            if (ptr->prev){
              ptr->prev->next = ptr->next;
            }
            if (ptr->next){
              ptr->next->prev = ptr->prev;
            }
          }
          delete ptr;
          first = 0;
        }

        printf("Deleting edge (%d, %d)\n", u, v);

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

  // If this triangle was successfully deleted, adjust
  // the triangle count
  if (success){
    num_triangles--;
  }

  printf("Done deleting triangle\n");

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
  int value = getEdgeHash(u, v);
  int bucket = value % num_buckets;
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
      u = segs[2*i];
    }
    if (segs[2*i+1] >= 0){
      v = segs[2*i+1];
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
double TMRTriangularize::inCircle( uint32_t u, uint32_t v,\
                                   uint32_t w, uint32_t x ){
  return incircle(&pts[2*u], &pts[2*v], &pts[2*w], &pts[2*x]);
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
  printf("Check if edge is PSLG\n"); fflush(stdout);
  if (edgeInPSLG(w, v)){
    TMRTriangle tri(u, v, w);
    addTriangle(tri);
    return; 
  }

  printf("Complete me\n"); fflush(stdout);
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
    printf("In circle test\n"); fflush(stdout);
    if (inCircle(w, v, x, u) >= 0.0){
      printf("delete triangle\n"); fflush(stdout);
      deleteTriangle(TMRTriangle(w, v, x));
      printf("digCavity()\n"); fflush(stdout);
      digCavity(u, v, x);
      digCavity(u, x, w);
      return;
    }
  }

  addTriangle(TMRTriangle(u, v, w));
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
  }

  // Add the point to the quadtree
  root->addNode(num_points, pt);

  // Set the new point location
  pts[2*num_points] = pt[0];
  pts[2*num_points+1] = pt[1];
  num_points++;

  // Return the new point index
  return num_points-1;
}

/*
  Add the vertex to the underlying Delaunay triangularization
*/
void TMRTriangularize::addPointToMesh( const double pt[] ){
  // Find the enclosing triangle
  TMRTriangle *tri;
  printf("Find enclosing\n");
  findEnclosing(pt, &tri);

  // Add the point to the quadtree
  printf("Adding point\n");
  uint32_t u = addPoint(pt);
  printf("Adding node to quadtree\n");
  root->addNode(u, pt);

  if (tri){
    uint32_t v = tri->u;
    uint32_t w = tri->v;
    uint32_t x = tri->w;
    printf("Triangle (%d, %d, %d)\n", v, w, x);
    deleteTriangle(*tri);
    digCavity(u, v, w);
    digCavity(u, w, x);
    digCavity(u, x, v);
  }

  return;
}

/*
  Add the point to the mesh
*/
void TMRTriangularize::addPointToMesh( const double pt[], 
                                       TMRTriangle *tri ){
  uint32_t u = addPoint(pt);
  root->addNode(u, pt);

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
  Find the enclosing triangle within the mesh. 

  This code uses the quadtree for geometric searching. First,
  we find the node that is closest to the query point. This 
  node is not necessarily connected with the enclosing triangle 
  that we want. Next, we find one triangle associated with this
  node. If this triangle does not contain the point, we march over
  the mesh, marking the triangles that we have visited. 
*/
void TMRTriangularize::findEnclosing( const double pt[],
                                      TMRTriangle **tri ){
  // Scan through the mesh, using a linear search
  *tri = NULL;
  TriListNode *node = list_root;
  while (node){
    if (enclosed(pt, node->tri.u, node->tri.v, node->tri.w)){
      *tri = &(node->tri);
      return;
    }
    node = node->next;
  }
}


/*
  Compute the distance to the intersection of the line m + alpha*e
  with the two lines from u to w and v to w
*/
void TMRTriangularize::computeIntersection( const double m[], 
                                            const double e[], 
                                            uint32_t u, uint32_t v, uint32_t w ){
  // m + alpha*e = pts[u] + beta*(pts[w] - pts[u])
  // => alpha*e + beta*(pts[u] - pts[w]) = pts[u] - m
  // a2 = [self.pts[u][i] - self.pts[w][i] for i in range(2)]

  // Check if the line is orthogonal to this direction
  /* if (e[0]*a2[1] - e[1]*a2[0] != 0.0){ */
  /*   b = [self.pts[u][i] - m[i] for i in range(2)] */
  /*     A = np.vstack((e, a2)).T */
  /*           ans = np.linalg.solve(A, b) */

  /*           # If beta is on the interval [0,1], alpha is the distance */
  /*           if ans[1] >= 0.0 and ans[1] <= 1.0: */
  /*               return ans[0] */

  /*       # If we get here and there's no solution, then we have a degenerate triangle */
  /*       b = [self.pts[v][i] - m[i] for i in range(2)] */
  /*       a2 = [self.pts[v][i] - self.pts[w][i] for i in range(2)] */
  /*       A = np.vstack((e, a2)).T */
  /*       ans = np.linalg.solve(A, b) */
            
  /*       # If beta is on the interval [0,1], alpha is the distance */
  /*       if ans[1] >= 0.0 and ans[1] <= 1.0: */
  /*           return ans[0] */

  /*       return -1.0 */
}

/*
void TMRTriangularize::frontal( double h ){
  // The length of the edge of the triangle
  double de = 0.5*sqrt(3.0)*h;

  // Add the triangles to the active set that 
  // status = {}
  // circumcircles = {}
  for key in self.tris:
    status[key] = 'waiting'

    // Check whether the triangle should be labeled as active
    t = self.tris[key]
    u, v, w = t
    for e in [(u,v), (v,w), (w,u)]:
      if e in self.pslg_edges:
        status[key] = 'active'

  // Compute the circumcircles
  for key in self.tris:
    circumcircles[key] = self.circumcircle(self.tris[key])

  // Compute the circumcircle radii of all triangles
  int itr = 0;
  
  while (num_active > 0){
            
            # Find the wors triangle    
            rmax = 0.0
            tri = -1
            # If we maintained a sorted list of the worst offenders, this
            # could be faster
            for key in circumcircles:
                if status[key] == 'active' and circumcircles[key] > rmax:
                    rmax = circumcircles[key]
                    tri = key

            # Now insert a new point based on the Voronoi criteria
            # Determine the triangle's accepted or boundary edge
            t = self.tris[tri]

            # Check if we are along an edge of the PSLG
            edge = None
            u, v, w = t
            for e in [(u,v), (v,w), (w,u)]:
                if e in self.pslg_edges:
                    edge = e

            # Check whether one of the adjacent edges is done
            if edge is None:
                for e in [(u,v), (v,w), (w,u)]:
                    index = self.get_triangle(e[1], e[0])
                    if index is not None and status[index] == 'done':
                        edge = e
                        break

            # Compute the location of the new point
            # | i   j   k |
            # | dx  dy  0 | = i*dy - j*dx
            # | 0   0   1 |
            u, v = edge
            w = self.adjacent(u, v)
            m = 0.5*np.array(
                [self.pts[u][0] + self.pts[v][0],
                 self.pts[u][1] + self.pts[v][1]])
            e = -np.array(
                [self.pts[v][1] - self.pts[u][1],
                 self.pts[u][0] - self.pts[v][0]])

            # Compute half the distance between the u and v points
            p = 0.5*np.sqrt(np.sum(e**2))
            e = 0.5*e/p

            # Try the size-optimal point first
            pt = m + de*e

            # Find the enclosing triangle for the 
            if not self.enclosed(t[0], t[1], t[2], pt):
                t = self.find_enclosing(pt)

                if t is None:
                    # Pick a point that is actually in the current triangle
                    q = self.computeIntersection(m, e, u, v, w)
                    rho = 0.5*q
                    pt = m + rho*e
                    t = self.tris[tri]

            # Use the edge for this triangle to determine where to insert
            # the new point
            old_tri = 1*self.tri_key
            self.deleted_tris = []
            self.add_vertex_frontal(pt, t)
            new_tri = self.tri_key

            # Free the deleted triangles
            for key in self.deleted_tris:
                circumcircles.pop(key)
                status.pop(key)

            # Compute the circumcircles of the new triangles and check
            # whether they belong in the done category or not...
            for key in range(old_tri, new_tri):
                circumcircles[key] = self.circumcircle(self.tris[key])
                if circumcircles[key] <= 1.5*h:
                    status[key] = 'done'
                else:
                    status[key] = 'waiting'

            # Mark the triangle with the edge that was just computed 
            # as being done, regardless of whether it is or not... 
            # Otherwise the same extrapolation point will be used which 
            # will duplicate points within the domain leading to chaos!
            if (u,v) in self.edge_to_tris:
                tri = self.edge_to_tris[(u,v)]
                status[tri] = 'done'

            # Label the new triangles with active status
            for key in range(old_tri, new_tri):
                if status[key] == 'done':
                    continue

                # Extract the triangle
                t = self.tris[key]
                u, v, w = t

                # Check the adjacent triangles
                for edge in [(u,v), (v,w), (w,u)]:
                    index = self.get_triangle(edge[1], edge[0])
                    if index is not None and status[index] == 'done':
                        status[key] = 'active'
                        break
                    if edge in self.pslg_edges:
                        status[key] = 'active'
                        break

        self.status = status

        return
*/
