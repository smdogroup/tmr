#ifndef TMR_GEOMETRY_H
#define TMR_GEOMETRY_H

#include "TMRBase.h"

/*
  The following header file contains the interface for the geometry/
  topology for the TMR objects. These vertex/edge/surface and volume
  objects are used to map the 

  The vertex, edge, face and volume classes are used in conjunction
  with the TMROct(Quad)Forest class to evaluate nodal locations with
  the mesh. These interfaces are designed to be overriden with an
  external geometry engine.
*/

// Declaration of the basic geometric types
class TMRGeoVertex;
class TMRGeoEdge;
class TMRGeoSurface;

/*
  The vertex class
*/
class TMRGeoVertex {
 public:
  TMRGeoVertex( TMRPoint *_x ){
    x = x;
  }

  // The point associated with the object
  TMRPoint *x;
};

/*
  The parametrization for the edge
*/
class TMRGeoEdge {
 public:
  // Get the parameter range for this edge
  virtual void getRange( double *tmin, double *tmax ) = 0;
  
  // Given the parametric point, evaluate the x,y,z location
  virtual int evalPoint( double t, TMRPoint *X ) = 0;

  // Given the parametric point, evaluate the derivative 
  virtual int evalDeriv( double t, TMRPoint *Xt ) = 0;

  // Find the surface u,v coordinates of the curve p(t)
  virtual int reparamOnSurface( TMRGeoSurface *surface, 
                                double t, int dir,
                                double *u, double *v ) = 0;

  // Integrate along the edge and return an array containing
  // the parametric locations to provide an even spacing
  double integrate( double t1, double t2, double tol,
                    double **tvals, double **dist, int *nvals );

 private:
  // The start and end vertices for this list
  TMRGeoVertex *v1, *v2;
  
  // The list of faces that refer to this edge
  class SurfaceList {
  public:
    TMRGeoSurface *surf;
    SurfaceList *next;
  } *faces;
};

/*
  The parametrization for the surface
*/
class TMRGeoSurface {
 public:
  // Get the parameter range for this surface
  virtual void getRange( double *umin, double *vmin,
                         double *umax, double *vmax ) = 0;
 
  // Given the parametric point, compute the x,y,z location
  virtual int evalPoint( double u, double v, TMRPoint *X ) = 0;

  // Given the parametric point, evaluate the first derivative 
  virtual int evalDeriv( double u, double v, 
                         TMRPoint *Xu, TMRPoint *Xv ) = 0;

 private:

  // The list of edges associated with the surface
  class EdgeList {
  public:
    TMRGeoEdge *edge;
    EdgeList *next;
  } *edges;
};

/*
  The geometric model container class
*/
class TMRGeoModel {
 public:
  TMRGeoModel(){}

  // ....
};






/*
  A searchable/sortable list of edge points
*/
class EdgeEntry {
 public:
  int32_t t;
  double x, y, z;
};

/*
  Class that stores e a searchable/sortable list of the face points
*/
class FaceEntry {
 public:
  int32_t u, v;
  double x, y, z;
};

/*
  The following class implements a look up table for either edges or
  faces. 

  The index along the face is an integer (int32_t) that is on the same
  interval used for the octants. The reason is that comparison between
  integers is exact and the edge locations can be generated once using
  either a TFI or by evaluating the points within the domain.
*/
class TMR_EdgeLookup {
 public:
  TMR_EdgeLookup( int32_t *t, TMRPoint *pts, int npts );
  ~TMR_EdgeLookup();

  // Evaluate a point 
  int evalPoint( int32_t t_int, TMRPoint *X );

 private:
  // Set 
  EdgeEntry *points;

  // The number of points along this edge
  int npts;
};

/*
  The following class implements a look up table for a face
*/
class TMR_SurfaceLookup {
 public:
  TMR_SurfaceLookup( int32_t *u, int32_t *v, TMRPoint *pts, int npts );
  ~TMR_SurfaceLookup();

  // Evaluate a point on the surface
  int evalPoint( int32_t u_int, int32_t v_int, TMRPoint *X );
  
 private:  
  // Store the points on this face
  FaceEntry *points;
  
  // The number of points
  int npts;
};

/*
  Transfinite interpolation edge class.

  This class is used to create a straight line interpolated between
  the starting and end locations.
*/
class TMR_TFIEdge {
 public:
  TMR_TFIEdge( TMRPoint *a, TMRPoint *b ); 
  ~TMR_TFIEdge();
  
  // Evaluate the point on the edge
  int evalPoint( int32_t t_int, TMRPoint *pt );

 private:
  // The start and end locations of the edge
  TMRPoint a, b;
};

/*
  Transfinite surface interpolation class

  Perform a transfinite interpolation between the edges that form a
  face. Note that the surface points are only defined at the look up
  points that are provided along each edge. 

  Note that this is not a continuous representation of a surface.
*/
class TMR_TFISurface {
 public:
  TMR_TFISurface( TMRPoint *corners, 
                  TMR_EdgeLookup *edges[] );

  // Evaluate a point on the surface determined using a transfinite
  // interpolation from the edges
  int evalPoint( int32_t u_int, int32_t v_int, TMRPoint *pt );

 private:
  // Store the corner points
  TMRPoint c[4];

  // Store pointer to the edges
  TMR_EdgeLookup *edges[4];
};

/*
  Transfinite interpolation class for the volumes/blocks

  Perform a transfinite interpolation of the points within a volume.
  In a similar manner to the surface, the interpolation is only
  defined along points.
*/
class TMR_TFIVolume {
 public:
  TMR_TFIVolume( TMRPoint *corners, 
                 TMR_EdgeLookup *edges[], 
                 TMR_SurfaceLookup *faces[] );

  // Evaluate the point on this surface
  int evalPoint( int32_t u_int, int32_t v_int, int32_t w_int,
                 TMRPoint *pt );

 private:
  // The corner points for this volume
  TMRPoint c[8];

  // The edge objects for this volume
  TMR_EdgeLookup *edges[12];

  // The face objects 
  TMR_SurfaceLookup *faces[6];
};

#endif // TMR_GEOMETRY_H
