#ifndef TMR_GEOMETRY_H
#define TMR_GEOMETRY_H

/*
  The following header file contains the interface for the geometry/
  topology for the TMR objects. These vertex/edge/surface and volume
  objects are used to map the 

  The vertex, edge, face and volume classes are used in conjunction
  with the TMROct(Quad)Forest class to evaluate nodal locations with
  the mesh. These interfaces are designed to be overriden with an
  external geometry engine.

*/

// Abstract base classes
class TMRPoint {
 public:
  virtual void eval( double *X ) = 0;
};

class TMREdge {
 public:
  virtual void getRange( double *lb, double *ub ) = 0;
  virtual void eval( double t, double *X ) = 0;
  virtual void invEval( const double *X, double *t ) = 0;
};

class TMRSurface {
 public:
  virtual void getRange( double *lb, double *ub ) = 0;
  virtual void eval( const double *u, double *X ) = 0;
  virtual void invEval( const double *X, double *u ) = 0;
};

class TMRVolume {
 public:
  virtual void getRange( double *lb, double *ub ) = 0;
  virtual void eval( const double *u, double *X ) = 0;
  virtual void invEval( const double *X, double *u ) = 0;
};







class TMR_TFIEdge : TMREdge {

};

class TMR_TFISurface : TMRSurface {
  
};

class TMR_TFIVolume {
 public:
  TMR_TFIVolume( );

};




#endif // TMR_GEOMETRY_H
