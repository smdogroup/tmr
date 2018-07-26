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

#ifndef TMR_FEATURE_SIZE_H
#define TMR_FEATURE_SIZE_H

#include "TMRBase.h"

/*
  The element feature size class
*/
class TMRElementFeatureSize : public TMREntity {
 public:
  TMRElementFeatureSize( double _hmin );
  virtual ~TMRElementFeatureSize();
  virtual double getFeatureSize( TMRPoint pt );

 protected:
  // The min local feature size
  double hmin;
};

/*
  Set a min/max element feature size
*/
class TMRLinearElementSize : public TMRElementFeatureSize {
 public:
  TMRLinearElementSize( double _hmin, double _hmax,
                        double c, double _ax, double _ay, double _az );
  ~TMRLinearElementSize();
  double getFeatureSize( TMRPoint pt );

 private:
  double hmax;
  double c, ax, ay, az;
};

/*
  Set a min/max feature size and feature sizes dictated within boxes
*/
class TMRBoxFeatureSize : public TMRElementFeatureSize {
 public:
  TMRBoxFeatureSize( TMRPoint p1, TMRPoint p2,
                     double _hmin, double _hmax );
  ~TMRBoxFeatureSize();
  void addBox( TMRPoint p1, TMRPoint p2, double h );
  double getFeatureSize( TMRPoint pt );

 private:
  // Maximum feature size
  double hmax;

  // Store the information about the points
  class BoxSize {
   public:
    // Check whether the box contains a point
    int contains( TMRPoint p );

    // Data for the box and its location
    TMRPoint m; // Center of the box
    TMRPoint d; // Half-edge length of each box
    double h; // Mesh size within the box
  };

  // The list of boxes that are stored
  int num_boxes;

  // A list object
  static const int MAX_LIST_BOXES = 256;
  class BoxList {
   public:
    BoxList(){ next = NULL; }
    BoxSize boxes[MAX_LIST_BOXES];
    BoxList *next;
  } *list_root, *list_current;

  // Store the information about the size
  class BoxNode {
   public:
    static const int MAX_NUM_BOXES = 10;

    BoxNode( BoxSize *cover, TMRPoint m, TMRPoint d );
    ~BoxNode();

    // Add a box
    void addBox( BoxSize *ptr );

    // Add a box to the data structure
    void getSize( TMRPoint pt, double *h );

    // The mid-point of the node and the distance from the mid-point
    // to the box sides
    TMRPoint m, d;

    // The children from this node
    BoxNode *c[8];

    // The number of boxes and the pointers to them
    int num_boxes;
    BoxSize **boxes;
  } *root;
};

/*
  Create a feature size on a finite-element mesh
*/
class TMRPointFeatureSize : public TMRElementFeatureSize {
 public:
  TMRPointFeatureSize( int npts, TMRPoint *pts, double *hvals,
                      double _hmin, double _hmax );
  ~TMRPointFeatureSize();
  double getFeatureSize( TMRPoint pt );

 private:
  // Maximum number of points used
  static const int MAX_CLOSEST_POINTS = 16;
  static const int MAX_BIN_SIZE = 8;

  // Private functions
  int split( int start, int end );
  int partitionPoints( TMRPoint *loc, TMRPoint *normal,
                       int *indx, int np );

  // Find the points that are closest to the provided point
  void locateClosest( const int K, const int root, const TMRPoint pt,
                      int *nk, int *indx, double *dist );

  // Insert the point into the list
  void insertIndex( const int K, int dindx, double d,
                    int *nk, int *indx, double *dist );

  // Max/min feature sizes
  double hmax, hmin;

  // Point data
  int npts;
  TMRPoint *pts;
  double *hvals;

  // Private data for a sorted spatial data structure
  int num_nodes, max_num_nodes;
  int *indices, *nodes;
  int *index_offset, *index_count;
  TMRPoint *node_loc, *node_normal;
};

#endif // TMR_FEATURE_SIZE_H
