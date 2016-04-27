#ifndef LOCATE_POINT_H
#define LOCATE_POINT_H

#define LAPACKsyevd dsyevd_

extern "C" {
  // Compute the eigenvalues of a symmetric matrix
  extern void LAPACKsyevd( const char *jobz, const char *uplo, int *N, 
                           double *A, int *lda, double *w, 
                           double *work, int *lwork, int *iwork, int *liwork, 
                           int *info );
}

/*! 
  Given a set of points in R^3, locate the closest one to a given
  point in O(log(N)) time -- after an initial O(N) setup time.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

class LocatePoint {
 public:
  LocatePoint( double *_Xpts, int _npts, int _max_num_points );
  ~LocatePoint();

  // Return the index of the point in the array
  // ------------------------------------------
  int locateClosest( double xpt[] ); 

  // Locate the K-closest points (note that dist/indices must of length K)
  // ---------------------------------------------------------------------
  void locateKClosest( int K, int indices[], 
                       double dist[], double xpt[] );

 private:
  // The recursive versions of the above functions
  void locateClosest( int root, double xpt[], 
                      double * dist, int * index );
  void locateKClosest( int K, int root, double xpt[], 
		       double * dist, int * indices, int * nk );

  // Insert the index into the sorted list of indices
  void insertIndex( double * dist, int * indices, int *nk, 
		    double d, int dindex, int K );

  // Sort the list of initial indices into the tree data structure
  int split( int start, int end );
  int splitList( double xav[], double normal[], 
                 int * indices, int npts );
  
  // Functions for array management
  void extendArrays( int old_len, int new_len );
  int * newIntArray( int * array, int old_len, int new_len );
  double * newDoubleArray( double * array, int old_len, int new_len );

  // The cloud of points to match
  double *Xpts;
  int npts; 

  int max_num_points; // Maximum number of points stored at a leaf

  // Keep track of the nodes that have been created
  int max_nodes;
  int num_nodes;

  int *indices; // Indices into the array of points
  int *nodes;  // Indices from the current node to the two child nodes
  int *indices_ptr; // Pointer into the global indices array
  int *num_indices; // Number of indices associated with this node
  double *node_xav; // Origin point for the array
  double *node_normal; // Normal direction of the plane
};

#endif
