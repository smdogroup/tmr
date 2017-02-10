#include "TMRTriangularize.h"

/*
  The following class implments a simple quadtree data structure for fast
  geometric searching queries.
*/
TMRQuadNode::TMRQuadNode( QuadNodeDomain *_domain,
                          int32_t _u, int32_t _v, int level ){
  // Set the level of this quadndoe
  level = _level;

  // Set the domain level
  u = _u;
  v = _v;

  // Set the maximum length in parameter space for
  // the domain
  const int32_t hmax = 1 << MAX_DEPTH;
  const int32_t h = 1 << (MAX_DEPTH - level - 1);

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
  pt_nums = new int[ NODES_PER_LEVEL ];
}
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
  not to put in duplicated entries.
*/
void TMRQuadNode::addNode( int num, const double pt[] ){
  // If any of the children have been allocated, searh them for
  // the place where the node should be added
  if (low_left){    
    if (pt[0] < x && pt[1] < y){
      low_left->add_point(num, pt);
    }
    else if (pt[0] < x){
      up_left->add_point(num, pt);
    }
    else if (pt[1] < y){
      low_right->add_point(num, pt);
    }
    else {
      up_right->add_point(num, pt);
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
    const int32_t h = 1 << (QUAD_NODE_MAX_DEPTH - level - 1);

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

  Note that the coordinate provided is not matched exactly, but
  needs to be used to scan through the quadtree hierarchy to find
  the correct quadrant location.
*/
int TMRQuadNode::deleteNode( int num, const double pt[] ){
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
    
*/
int TMRQuadNode::findClosest( const double pt[], double *_dist ){
  double dist = MAX_QUAD_DISTANCE;
  int index = -1;

  // Make the recursive call
  findClosest(pt, &index, &dist);
  if (_dist){ *_dist = dist; }
  return index;
}

  
// The recursive call to the quadtree data structure
void TMRQuadNode::findClosest( const double pt[], int *index, double *dist ){
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
  Retrieve the edge hash
*/
inline int TMRTriangularize::getEdgeHash( int32_t x, int32_t y ){
  return (x * 0x1f1f1f1f) ^ y;
}

/*
  Add a triangle to the mesh
*/
void TMRTriangularize::addTriangle( TMRTriangle tri ){
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
    // Create a new array of buckets, twice the size of the old array
    // of buckets and re-insert the entries back into the new array
    int num_old_buckets = num_buckets;
    num_buckets *= 2;
    TriHashNode **new_buckets = new TriHashNode*[ num_buckets ];
    memset(new_buckets, 0, num_buckets*sizeof(TriHashNode*));

    // Assign the new buckets
    for ( int i = 0; i < num_old_buckets; i++ ){
      TriHashNode *node = buckets[i];

      while (node){
        // Get the new hash values
        TriHashNode *tmp = node->next;
        int value = getEdgeHash(node->u, node->v);
        int bucket = value % num_buckets;
        
        if (!new_buckets[bucket]){
          new_buckets[bucket] = node;
          end_buckets[bucket] = new_buckets[bucket];
        }
        else {
          end_bucket[bucket]->next = node;
          end_bucket[bucket] = end_bucket[bucket]->next;
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
  int32_t edge_pairs[] = {tri.u, tri.v,
                          tri.v, tri.w,
                          tri.w, tri.u};
  
  // Add the triangle to the hash table
  for ( int k = 0; k < 3; k++ ){
    // Add a hash for each pair of edges around the triangle
    int32_t u = egde_pairs[2*k];
    int32_t v = edge_pairs[2*k+1];
    int hash = getEdgeHash(u, v);
    int bucket = hash % num_buckets;
    if (!hash[bucket]){
      // Create the new hash node and assign the values
      hash[bucket] = new TriHashNode();
      hash[bucket]->u = u;
      hash[bucket]->v = v;
      hash[bucket]->node = list_current;
      hash[bucket]->next = NULL;
    }
    else {
      // Scan through to the end of the array and set the pointer to
      // the triangle node that contains the 
      TriHashNode *node = hash[bucket];
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
        node->next = new TriHashNode();
        node = node->next;
        node->u = u;
        node->v = v;
        node->tri_node = list_current;
        node->next = NULL;
      }
    }
  }

  return success;
}

/*
  Delete a triangle from the mesh
*/
void TMRTriangularize::deleteTriangle( TMRTriangle tri ){
  // Keep track of whether we successfully delete all of the edges, or
  // just some of the edges (1 or 2 out of 3 is bad!)
  int success = 1;

  // Set the combinations of edge pairs that will be added
  // to the hash table
  int32_t edge_pairs[] = {tri.u, tri.v,
                          tri.v, tri.w,
                          tri.w, tri.u};

  // Keep track if this is the first edge we find
  int first = 1;
  
  // Add the triangle to the hash table
  for ( int k = 0; k < 3; k++ ){
    // Add a hash for each pair of edges around the triangle
    int32_t u = egde_pairs[2*k];
    int32_t v = edge_pairs[2*k+1];

    // Keep track of whether we've deleted this edge
    int edge_success = 0;

    // Get the hash values and access the first entry of the bucket
    // it's listed under
    int value = getEdgeHash(u, v);
    int bucket = value % num_buckets;
    TriHashNode *node = buckets[bucket];
    TriHashNode *prev = node;
    while (node){
      // The edge matches, we have to delete this triangle
      if (u == node->u && v == node->v){
        // If this is the first edge we've found, delete the
        // corresponding list entry as well.
        if (first){
          // Adjust the double linking so that when the node is
          // deleted, the pointers will still work
          TriListNode *ptr = node->tri_node;
          if (ptr->prev){
            ptr->prev->next = ptr->next;
          }
          if (ptr->next){
            ptr->next->prev = ptr->prev;
          }
          delete ptr;
          first = 0;
        }

        // Delete the edge from the hash table
        if (node == prev){
          buckets[bucket] = node->next;
          delete node;
        }
        else {
          prev->next = node->next;
          delete node;
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
void TMRTriangularize::completeMe( int32_t u, int32_t v,
                                   TMRTriangle **tri ){
  *tri = NULL;

  // Retrieve the hash value/bucket for this edge
  int value = getEdgeHash(u, v);
  int bucket = value % num_buckets;
  TriHashNode *node = buckets[bucket];

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

    def enclosed(self, u, v, w, pt):
        '''Does this triangle enclose this point?'''
        if (self.orient2d(u, v, pt) and
            self.orient2d(v, w, pt) and
            self.orient2d(w, u, pt)):
            return True
        return False


    def orient2d(self, a, b, pt):
        '''Check the relative orientation of the points a, b, and pt'''
        A = np.array([
            [self.pts[a][0] - pt[0], self.pts[a][1] - pt[1]],
            [self.pts[b][0] - pt[0], self.pts[b][1] - pt[1]]])
        # This is NOT robust to precision errors
        return  A[0,0]*A[1,1] - A[0,1]*A[1,0] >= 0.0
/*
  

 */
void Triangularize::digCavity( int32_t u, int32_t v, int32_t w ){ 
   // If the edge is along the polynomial straight line graph, then we
   // add the triangle as it exists and we're done, even though it may
   // not be Delaunay (it will be constrained Delaunay). We cannot
   // split a PSLG edge.
   if (edgeInPSLG(w, v)){
     TMRTriangle tri(u, v, w);
     addTriangle(tri);
     return; 
   }

   // Complete the triangle
   TMRTriangle *tri;
   completeMe(w, v, &tri);

   if (tri){
     if (inCircle(u, w, v, x)){
       deleteTriangle(TMRTriangle(w, v, x));
       digCavity(u, v, x);
       digCavity(u, x, w);
       return;
     }
   }

   addTriangle(u, v, w);
 }

/*
  Add the vertex to the underlying Delaunay triangularization
*/
void Triangularize::addVertex( const double pt[] ){
  // Find the enclosing triangle
  TMRTriangle *tri;
  findEnclosing(pt, &tri);

  // Add the point to the quadtree
  root->addNode(pt);

  int32_t u = num_points;
  num_points++;

  if (tri){
    int32_t v = tri->u;
    int32_t w = tri->v;
    int32_t x = tri->x;
    deleteTriangle(*tri);
    
    digCavity(u, v, w);
    digCavity(u, w, x);
    digCavity(u, x, v);
  }

  return;
}

void TMRTriangularize::addVertex( const double pt[], 
                                  TMRTriangle *tri ){
  root->addNode(pt);
  int32_t u = num_points;
  num_points++;

  if (tri){
    int32_t v = tri->u;
    int32_t w = tri->v;
    int32_t x = tri->x;
    deleteTriangle(*tri);
    
    digCavity(u, v, w);
    digCavity(u, w, x);
    digCavity(u, x, v);
  }
}


void TMRTriangularize::findEnclosing( const double pt[] ){

}


/*
  Compute the distance to the intersection of the line m + alpha*e
  with the two lines from u to w and v to w
*/
def computeIntersection( const double m[], const double e[], 
                         int32_t u, int32_t v, int32_t w ){
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
  void frontal( 
    def frontal(self, h, plot=False, freq=50):
        '''Refine the triangles using the frontal method'''

        # The length of the edge of the triangle
        de = 0.5*np.sqrt(3.0)*h

        # Add all of the triangles to the active set
        status = {}
        circumcircles = {}
        for key in self.tris:
            status[key] = 'waiting'

            # Check whether the triangle should be labeled as active
            t = self.tris[key]
            u, v, w = t
            for e in [(u,v), (v,w), (w,u)]:
                if e in self.pslg_edges:
                    status[key] = 'active'

        # Compute the circumcircles
        for key in self.tris:
            circumcircles[key] = self.circumcircle(self.tris[key])

        # Compute the circumcircle radii of all triangles
        itr = 0
        while 'active' in status.values():
            if itr % freq == 0:
                print 'iteration = %d'%(itr)
                if plot:
                    self.status = status
                    self.plot()
                    plt.show()
            itr += 1
            
            # Find the wors ttriangle    
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
