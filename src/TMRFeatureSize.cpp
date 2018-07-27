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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TMRFeatureSize.h"
#include "tmrlapack.h"

/*
  Create the element feature size.

  This corresponds to a constant feature size across the entire
  domain. Element grading can be implemented by overriding this class
  to modify the feature size as a function of position.
*/
TMRElementFeatureSize::TMRElementFeatureSize( double _hmin ){
  hmin = _hmin;
}

TMRElementFeatureSize::~TMRElementFeatureSize(){}

/*
  Return the feature size - hmin
*/
double TMRElementFeatureSize::getFeatureSize( TMRPoint pt ){
  return hmin;
}

/*
  Create a feature size dependency that is linear but does not
  exceed hmin or hmax anywhere in the domain
*/
TMRLinearElementSize::TMRLinearElementSize( double _hmin, double _hmax,
                                            double _c,
                                            double _ax,
                                            double _ay,
                                            double _az ):
TMRElementFeatureSize(_hmin){
  hmax = _hmax;
  c = _c;
  ax = _ax;
  ay = _ay;
  az = _az;
}

TMRLinearElementSize::~TMRLinearElementSize(){}

/*
  Get the feature size based on the input coefficients
*/
double TMRLinearElementSize::getFeatureSize( TMRPoint pt ){
  double h = c + ax*pt.x + ay*pt.y + az*pt.z;
  if (h < hmin){ h = hmin; }
  if (h > hmax){ h = hmax; }
  return h;
}

/*
  Create the feature size within a box
*/
TMRBoxFeatureSize::TMRBoxFeatureSize( TMRPoint p1, TMRPoint p2,
                                      double _hmin, double _hmax ):
TMRElementFeatureSize(_hmin){
  hmax = _hmax;

  TMRPoint m1, d1;
  m1.x = 0.5*(p1.x + p2.x);
  m1.y = 0.5*(p1.y + p2.y);
  m1.z = 0.5*(p1.z + p2.z);
  d1.x = 0.5*fabs(p1.x - p2.x);
  d1.y = 0.5*fabs(p1.y - p2.y);
  d1.z = 0.5*fabs(p1.z - p2.z);

  // Allocate the root for the boxes
  list_current = new BoxList();
  list_root = list_current;

  // Set the initial number of boxes
  num_boxes = 0;

  // Allocate a box
  list_current->boxes[num_boxes].m = m1;
  list_current->boxes[num_boxes].d = d1;
  list_current->boxes[num_boxes].h = hmax;

  // Add a box to the data structure
  root = new BoxNode(&(list_current->boxes[num_boxes]), m1, d1);
  num_boxes++;
}

/*
  Free the data allocated by the feature-size object
*/
TMRBoxFeatureSize::~TMRBoxFeatureSize(){
  delete root;
  while (list_root){
    BoxList *tmp = list_root;
    list_root = list_root->next;
    delete tmp;
  }
}

/*
  Add a box that is designed to constrain the feature size
*/
void TMRBoxFeatureSize::addBox( TMRPoint p1, TMRPoint p2, double hval ){
  TMRPoint m1, d1;
  m1.x = 0.5*(p1.x + p2.x);
  m1.y = 0.5*(p1.y + p2.y);
  m1.z = 0.5*(p1.z + p2.z);
  d1.x = 0.5*fabs(p1.x - p2.x);
  d1.y = 0.5*fabs(p1.y - p2.y);
  d1.z = 0.5*fabs(p1.z - p2.z);

  // Check if we need to expand the array of boxes
  if (num_boxes >= MAX_LIST_BOXES){
    list_current->next = new BoxList();
    list_current = list_current->next;
    num_boxes = 0;
  }

  // Allocate a box
  list_current->boxes[num_boxes].m = m1;
  list_current->boxes[num_boxes].d = d1;
  list_current->boxes[num_boxes].h = hval;

  // Add a box to the data structure
  root->addBox(&(list_current->boxes[num_boxes]));
  num_boxes++;
}

/*
  Retrieve the feature size
*/
double TMRBoxFeatureSize::getFeatureSize( TMRPoint pt ){
  double h = hmax;
  root->getSize(pt, &h);
  if (h < hmin){ h = hmin; }
  if (h > hmax){ h = hmax; }
  return h;
}

/*
  Check if the box contains the point
*/
int TMRBoxFeatureSize::BoxSize::contains( TMRPoint p ){
  // Compute the lower/upper bounds
  double xl = m.x - d.x;
  double xu = m.x + d.x;
  double yl = m.y - d.y;
  double yu = m.y + d.y;
  double zl = m.z - d.z;
  double zu = m.z + d.z;

  // Return true if the box contains the point
  if ((p.x >= xl && p.x <= xu) &&
      (p.y >= yl && p.y <= yu) &&
      (p.z >= zl && p.z <= zu)){
    return 1;
  }
  return 0;
}

/*
  Create the box node and allocate a single box (for now)
*/
TMRBoxFeatureSize::BoxNode::BoxNode( BoxSize *cover,
                                     TMRPoint m1,
                                     TMRPoint d1 ){
  m = m1;
  d = d1;
  num_boxes = 1;
  boxes = new BoxSize*[ MAX_NUM_BOXES ];
  boxes[0] = cover;
  memset(c, 0, sizeof(c));
}

/*
  Deallocate all of the data
*/
TMRBoxFeatureSize::BoxNode::~BoxNode(){
  if (boxes){
    delete [] boxes;
  }
  for ( int i = 0; i < 8; i++ ){
    if (c[i]){
      delete c[i];
    }
  }
}

/*
  Add a box to the octree data structure
*/
void TMRBoxFeatureSize::BoxNode::addBox( BoxSize *ptr ){
  // Get the dimensions for the box
  double xl = ptr->m.x - ptr->d.x;
  double xu = ptr->m.x + ptr->d.x;
  double yl = ptr->m.y - ptr->d.y;
  double yu = ptr->m.y + ptr->d.y;
  double zl = ptr->m.z - ptr->d.z;
  double zu = ptr->m.z + ptr->d.z;

  // Get the dimensions for this node
  double nxl = m.x - d.x;
  double nxu = m.x + d.x;
  double nyl = m.y - d.y;
  double nyu = m.y + d.y;
  double nzl = m.z - d.z;
  double nzu = m.z + d.z;

  // Check if this box coverss any part of this node, if not we're
  // done
  if (!((xu >= nxl && xl <= nxu) &&
        (yu >= nyl && yl <= nyu) &&
        (zu >= nzl && zl <= nzu))){
    return;
  }

  // Check if the boxes at this level are currently in use
  if (boxes){
    // Check if this box covers the entire node or not
    if ((xl <= nxl && xu >= nxu) &&
        (yl <= nyl && yu >= nyu) &&
        (zl <= nzl && zu >= nzu)){
      // This box fully covers the node. Check if the h-dimension is
      // more strict than the current covering box
      if (boxes[0]->h > ptr->h){
        boxes[0] = ptr;
      }
      return;
    }

    // Otherwise, the box only covers part of the node and we have
    // to add it to the list of local boxes
    if (num_boxes < MAX_NUM_BOXES){
      boxes[num_boxes] = ptr;
      num_boxes++;
      return;
    }
    else {
      // Allocate new children
      for ( int k = 0; k < 2; k++ ){
        double zl = m.z + (k-1)*d.z;
        double zu = m.z + k*d.z;

        for ( int j = 0; j < 2; j++ ){
          double yl = m.y + (j-1)*d.y;
          double yu = m.y + j*d.y;

          for ( int i = 0; i < 2; i++ ){
            double xl = m.x + (i-1)*d.x;
            double xu = m.x + i*d.x;

            // Set the mid-points and distances for each child
            TMRPoint m1, d1;
            m1.x = 0.5*(xu + xl);
            m1.y = 0.5*(yu + yl);
            m1.z = 0.5*(zu + zl);
            d1.x = 0.5*(xu - xl);
            d1.y = 0.5*(yu - yl);
            d1.z = 0.5*(zu - zl);

            // Create the child nodes
            c[i + 2*j + 4*k] = new BoxNode(boxes[0], m1, d1);
          }
        }
      }

      // Add the boxes to the rest of the tree
      while (num_boxes > 0){
        for ( int k = 0; k < 2; k++ ){
          for ( int j = 0; j < 2; j++ ){
            for ( int i = 0; i < 2; i++ ){
              c[i + 2*j + 4*k]->addBox(boxes[num_boxes-1]);
            }
          }
        }
        num_boxes--;
      }

      // Free the boxes
      delete [] boxes;
      boxes = NULL;
    }
  }

  // Add the boxes to the child - if applicable
  for ( int k = 0; k < 2; k++ ){
    for ( int j = 0; j < 2; j++ ){
      for ( int i = 0; i < 2; i++ ){
        c[i + 2*j + 4*k]->addBox(ptr);
      }
    }
  }
}

/*
  Get the element size based on the extent of the box
*/
void TMRBoxFeatureSize::BoxNode::getSize( TMRPoint p, double *hval ){
  // Scan through the tree and find the most-constraining box size
  for ( int k = 0; k < 2; k++ ){
    double zl = m.z + (k-1)*d.z;
    double zu = m.z + k*d.z;

    for ( int j = 0; j < 2; j++ ){
      double yl = m.y + (j-1)*d.y;
      double yu = m.y + j*d.y;

      for ( int i = 0; i < 2; i++ ){
        double xl = m.x + (i-1)*d.x;
        double xu = m.x + i*d.x;

        // Check if the point lies in the box and the box exists
        if ((p.x >= xl && p.x <= xu) &&
            (p.y >= yl && p.y <= yu) &&
            (p.z >= zl && p.z <= zu)){
          if (c[i + 2*j + 4*k]){
            return c[i + 2*j + 4*k]->getSize(p, hval);
          }
        }
      }
    }
  }

  for ( int i = 0; i < num_boxes; i++ ){
    if (boxes[i]->contains(p) && *hval > boxes[i]->h){
      *hval = boxes[i]->h;
    }
  }
}

/*
  Set the feature size based on a point cloud and the mesh size at
  those points
*/
TMRPointFeatureSize::TMRPointFeatureSize( int _npts, TMRPoint *_pts,
                                          double *_hvals,
                                          double _hmin, double _hmax ):
  TMRElementFeatureSize(_hmin){
  npts = _npts;
  hmin = _hmin;
  hmax = _hmax;

  pts = new TMRPoint[ npts ];
  hvals = new double[ npts ];
  memcpy(pts, _pts, npts*sizeof(TMRPoint));
  memcpy(hvals, _hvals, npts*sizeof(double));

  // Calculate approximately how many nodes there should be
  max_num_nodes = int(2.0*npts/MAX_BIN_SIZE) + 1;
  num_nodes = 0;

  // The point indicies
  indices = new int[ npts ];
  for ( int i = 0; i < npts; i++ ){
    indices[i] = i;
  }

  // Set up the data structure that represents the splitting planes
  nodes = new int[ 2*max_num_nodes ];
  index_offset = new int[ max_num_nodes ];
  index_count = new int[ max_num_nodes ];

  // The base point and normal direction for the splitting planes
  node_loc = new TMRPoint[ max_num_nodes ];
  node_normal = new TMRPoint[ max_num_nodes ];
  memset(node_loc, 0, max_num_nodes*sizeof(TMRPoint));
  memset(node_normal, 0, max_num_nodes*sizeof(TMRPoint));

  // Recursively split the points
  split(0, npts);
}

TMRPointFeatureSize::~TMRPointFeatureSize(){
  delete [] pts;
  delete [] hvals;
  delete [] indices;
  delete [] nodes;
  delete [] index_offset;
  delete [] index_count;
  delete [] node_loc;
  delete [] node_normal;
}

double TMRPointFeatureSize::getFeatureSize( TMRPoint pt ){
  // Get the closest points and use them to compute a set of weights
  int indx[MAX_CLOSEST_POINTS];
  double dist[MAX_CLOSEST_POINTS];
  double weights[MAX_CLOSEST_POINTS];

  // Find the closest points
  int n = 0;
  int root = 0;
  locateClosest(MAX_CLOSEST_POINTS, root, pt, &n, indx, dist);

  // The maximum h-size squared
  double d = 10*hmax;
  double dinv = 1.0/d;

  // Set the first weight (on the closest point) to 1
  double wsum = 1.0;
  weights[0] = 1.0;

  // Compute the weights
  for ( int i = 0; i < n; i++ ){
    dist[i] = sqrt(dist[i]);
    if (dist[i] < d){
      weights[i] = 1.0 - dinv*dist[i];
    }
    else {
      weights[i] = 0.0;
    }
    if (i == 0){
      wsum = weights[i];
    }
    else {
      wsum += weights[i];
    }
  }

  // Set the weights so that they satisfy a partition of unity
  // property
  wsum = 1.0/wsum;

  // Compute the new value of h
  double h = 0.0;
  for ( int i = 0; i < n; i++ ){
    h += wsum*weights[i]*hvals[indx[i]];
  }

  if (h < hmin){ h = hmin; }
  if (h > hmax){ h = hmax; }

  return h;
}

/*
  Split the list of indices into approximately two. Those on one half
  of a plane and those on the other.
*/
int TMRPointFeatureSize::split( int start, int end ){
  int root = num_nodes;

  num_nodes++;

  // Need to extend the arrays to make them fit
  if (num_nodes >= max_num_nodes){
    int max_num_nodes = 2*(num_nodes+1);

    // Allocate and set a new nodes pointer
    int *temp_nodes = new int[ 2*max_num_nodes ];
    memcpy(temp_nodes, nodes, 2*num_nodes*sizeof(int));
    delete [] nodes;
    nodes = temp_nodes;

    // Allocate and set a new offset array
    int *temp_index_offset = new int[ max_num_nodes ];
    memcpy(temp_index_offset, index_offset, num_nodes*sizeof(int));
    delete [] index_offset;
    index_offset = temp_index_offset;

    // Allocate more space for the index counts
    int *temp_index_count = new int[ max_num_nodes ];
    memcpy(temp_index_count, index_count, num_nodes*sizeof(int));
    delete [] index_count;
    index_count = temp_index_count;

    // Allocate more space for the average plane locations
    TMRPoint *temp_node_loc = new TMRPoint[ max_num_nodes ];
    memcpy(temp_node_loc, node_loc, num_nodes*sizeof(TMRPoint));
    delete [] node_loc;
    node_loc = temp_node_loc;

    // Allocate more space for the plane normals
    TMRPoint *temp_node_normal = new TMRPoint[ max_num_nodes ];
    memcpy(temp_node_normal, node_normal, num_nodes*sizeof(TMRPoint));
    delete [] node_normal;
    node_normal = temp_node_normal;
  }

  // If there are fewer than the max number of points in a bin, end
  // the recursion here and store the result of the last split
  if (end - start <= MAX_BIN_SIZE){
    nodes[2*root] = -1;
    nodes[2*root+1] = -1;

    // Set the offset into the node
    index_offset[root] = start;
    index_count[root] = end - start;

    return root;
  }

  // If this is not a leaf, set the number of indices to zero and the
  // index pointer to a negative value to indicate that this is not a
  // leaf. Sort the indices so that the points on one side of the
  // plane are below the mid point and the others are above.
  index_offset[root] = -1;
  index_count[root] = 0;

  // Split the list
  int mid = partitionPoints(&node_loc[root], &node_normal[root],
                            &indices[start], end - start);

  if (mid == 0 || mid == end-start){
    fprintf(stderr, "TMRPointFeatureSize: Error, splitting points did nothing \
-- problem with your nodes?\n");
    return root;
  }

  // Now, split the right and left hand sides of the list to continue
  // the recursion.
  int left_node = split(start, start + mid);
  int right_node = split(start + mid, end);
  nodes[2*root] = left_node;
  nodes[2*root+1] = right_node;

  return root;
}

/*
  Split the array of indices into two sets: those indices that
  correspond to points on either side of a plane in three-space.
*/
int TMRPointFeatureSize::partitionPoints( TMRPoint *loc, TMRPoint *normal,
                                          int *indx, int np ){
  // The inertia about the average location of the point cloud
  double I[9];
  memset(I, 0, 9*sizeof(double));

  // Zero the average location
  double x = 0.0, y = 0.0, z = 0.0;

  // Find the average point and the moment of inertia about the
  // average point
  for ( int i = 0; i < np; i++ ){
    int n = indx[i];

    // Keep track of the average location
    x += pts[n].x;
    y += pts[n].y;
    z += pts[n].z;

    // The moments of inertia
    I[0] += (pts[n].y*pts[n].y + pts[n].z*pts[n].z); // y^2 + z^2
    I[4] += (pts[n].x*pts[n].x + pts[n].z*pts[n].z); // x^2 + z^2
    I[8] += (pts[n].x*pts[n].x + pts[n].y*pts[n].y); // x^2 + y^2

    // The products of inertia
    I[1] -= pts[n].x*pts[n].y; // Ixy = - xy
    I[2] -= pts[n].x*pts[n].z; // Ixz = - xz
    I[5] -= pts[n].y*pts[n].z; // Ixz = - yz
  }

  // Set the average location
  x = x/np;
  y = y/np;
  z = z/np;

  // Convert the inertia so that is about the average location using
  // the parallel axis theorem
  I[0] = I[0] - np*(y*y + z*z);
  I[4] = I[4] - np*(x*x + z*z);
  I[8] = I[8] - np*(x*x + y*y);

  I[1] = I[1] + np*x*y;
  I[2] = I[2] + np*x*z;
  I[5] = I[5] + np*y*z;

  // Copy over the symmetric part of the products of inertia
  I[3] = I[1];
  I[6] = I[2];
  I[7] = I[5];

  // Find the eigenvalues/eigenvectors
  int info = 0, N = 3;
  int lwork = 64, liwork = 64;
  double eigs[3], work[64];
  int iwork[64];
  TmrLAPACKsyevd("V", "U", &N, I, &N,
                 eigs, work, &lwork, iwork, &liwork, &info);

  // Extract the normal
  double nx = I[0];
  double ny = I[1];
  double nz = I[2];

  // Loop over all the indices and decide where they should go
  int low = 0;
  int high = np-1;

  // Now, split the index array such that the indices below
  // the lower index
  while (high > low){
    while (high > low &&
           ((pts[indx[low]].x - x)*nx +
            (pts[indx[low]].y - y)*ny +
            (pts[indx[low]].z - z)*nz) < 0.0){
      low++;
    }
    while (high > low &&
           ((pts[indx[high]].x - x)*nx +
            (pts[indx[high]].y - y)*ny +
            (pts[indx[high]].z - z)*nz) >= 0.0){
      high--;
    }

    if (high > low){
      // Switch the two indices that don't match
      int temp = indx[high];
      indx[high] = indx[low];
      indx[low] = temp;
    }
  }

  // Set the average point location
  loc->x = x;
  loc->y = y;
  loc->z = z;

  // Set the normal point location
  normal->x = nx;
  normal->y = ny;
  normal->z = nz;

  // If this didn't work, issue an error message
  if (low == 0 || low == np){
    fprintf(stderr, "TMRPointFeatureSize: Error split points\n");
  }

  return low;
}

/*!
  Locate the closest points to a given point

  input:
  K:      The number of closest points to find
  root:   The root node index to search
  pt:     The point

  input/output:
  dist:   A sorted list of the K-closest distances
  indx:   The indices of the K-closest values
  nk:     The actual number of points in the list nk <= K
*/
void TMRPointFeatureSize::locateClosest( const int K, const int root,
                                         const TMRPoint pt,
                                         int *nk, int *indx, double *dist ){
  int start = index_offset[root];
  int left_node = nodes[2*root];
  int right_node = nodes[2*root+1];

  if (start != -1){
    // This node is a leaf. Do an exhaustive search of the points
    // to find the ones that are closest to the given point
    int end = start + index_count[root];

    // Loop over the indices in the leaf
    for ( int k = start; k < end; k++ ){
      int n = indices[k];

      // Compute the square of the distances
      double t = ((pts[n].x - pt.x)*(pts[n].x - pt.x) +
                  (pts[n].y - pt.y)*(pts[n].y - pt.y) +
                  (pts[n].z - pt.z)*(pts[n].z - pt.z));

      // Insert the point if needed
      if ((*nk < K) || (t < dist[K-1])){
        insertIndex(K, n, t, nk, indx, dist);
      }
    }
  }
  else {
    // This is not a leaf node. Figure out which side we should search
    // and then perform the recursive search on that side.
    double x = node_loc[root].x;
    double y = node_loc[root].y;
    double z = node_loc[root].z;

    // Find the normal
    double nx = node_normal[root].x;
    double ny = node_normal[root].y;
    double nz = node_normal[root].z;

    // The normal distance
    double ndist = ((pt.x - x)*nx +
                    (pt.y - y)*ny +
                    (pt.z - z)*nz);

    if (ndist < 0.0){ // The point lies to the 'left' of the plane
      locateClosest(K, left_node, pt, nk, indx, dist);

      // If the minimum distance to the plane is less than the minimum
      // distance then search the other branch too - there could be a
      // point on that branch that lies closer than *dist
      if (*nk < K || ndist*ndist < dist[*nk-1]){
        locateClosest(K, right_node, pt, nk, indx, dist);
      }
    }
    else { // The point lies to the 'right' of the plane
      locateClosest(K, right_node, pt, nk, indx, dist);

      // If the minimum distance to the plane is less than the minimum
      // distance then search the other branch too - there could be a
      // point on that branch that lies closer than *dist
      if (*nk < K || ndist*ndist < dist[*nk-1]){
        locateClosest(K, left_node, pt, nk, indx, dist);
      }
    }
  }
}

/*
  Insert a point into a sorted list based upon the distance from the
  given point

  input:
  K:      the maximum length of the array
  dindx:  the index of the new point to insert
  d:      the square of the distance to the new point

  input/output:
  nk:     the length of the sorted index/distance arrays nk <= K
  indx:   the sorted list of index values
  dist;   the distances
*/
void TMRPointFeatureSize::insertIndex( const int K, int dindx,
                                       double d, int *nk,
                                       int *indx, double *dist ){
  if (*nk == 0){
    dist[*nk] = d;
    indx[*nk] = dindx;
    *nk += 1;
    return;
  }
  else if (*nk < K && dist[*nk-1] <= d){
    dist[*nk] = d;
    indx[*nk] = dindx;
    *nk += 1;
    return;
  }

  // Place it into the list
  int i = 0;
  while (i < *nk && (d >= dist[i])){
    i++;
  }

  for ( ; i < *nk; i++ ){
    int tindx = indx[i];
    double t = dist[i];
    indx[i] = dindx;
    dist[i] = d;
    dindx = tindx;
    d = t;
  }

  if (*nk < K){
    indx[*nk] = dindx;
    dist[*nk] = d;
    *nk += 1;
  }
}
