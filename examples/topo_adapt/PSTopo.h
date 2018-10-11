#ifndef PLANE_STRESS_TOPO_H
#define PLANE_STRESS_TOPO_H

#include "TACSAssembler.h"
#include "PlaneStressStiffness.h"
#include "TACS2DElement.h"

class PSTopo : public PlaneStressStiffness {
 public:
  static const int MAX_NUM_MATERIALS = 12;
  static const int MAX_NUM_WEIGHTS = 15;

  PSTopo( TacsScalar _rho, TacsScalar E, TacsScalar nu, TacsScalar ys,
          double q, double eps,
          int _nodes[], double _weights[], int _nweights );
  ~PSTopo();

  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lb[], TacsScalar ub[],
                          int numDVs );
  void calculateStress( const double pt[], const TacsScalar e[],
                        TacsScalar s[] );
  void addStressDVSens( const double pt[], const TacsScalar e[],
                        TacsScalar alpha, const TacsScalar psi[],
                        TacsScalar fdvSens[], int dvLen );
  void getPointwiseMass( const double pt[], TacsScalar mass[] );
  void addPointwiseMassDVSens( const double pt[],
                               const TacsScalar alpha[],
                               TacsScalar fdvSens[], int dvLen );
  void failure( const double pt[], const TacsScalar strain[],
                TacsScalar * fail );
  void failureStrainSens( const double pt[], const TacsScalar strain[],
                          TacsScalar sens[] );
  void addFailureDVSens( const double pt[], const TacsScalar strain[],
                         TacsScalar alpha, TacsScalar dvSens[], int dvLen );
  TacsScalar getDVOutputValue( int dvIndex, const double pt[] ){
    return xw;
  }
  TacsScalar getDensity(){ return xw; }

 public:
  // Material property information
  TacsScalar rho, D, G, nu, ys;
  double q, eps;

  // Maximum number of weights/nodes
  int nweights;

  // weights/index values
  class Weight {
  public:
    int index;
    double weight;
  };

  // Weights used in the computation
  Weight *w;

  // The design variable values
  TacsScalar xw;
};

class PSTopo4 : public PlaneStressStiffness {
 public:
  static const int MAX_NUM_MATERIALS = 12;
  static const int MAX_NUM_WEIGHTS = 15;

  PSTopo4( TacsScalar _rho, TacsScalar E, TacsScalar nu, TacsScalar ys,
           double q, double eps,
           int *_nodes[], double *_weights[], int _nweights[] );
  ~PSTopo4();

  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lb[], TacsScalar ub[],
                          int numDVs );
  void calculateStress( const double pt[], const TacsScalar e[],
                        TacsScalar s[] );
  void addStressDVSens( const double pt[], const TacsScalar e[],
                        TacsScalar alpha, const TacsScalar psi[],
                        TacsScalar fdvSens[], int dvLen );
  void getPointwiseMass( const double pt[], TacsScalar mass[] );
  void addPointwiseMassDVSens( const double pt[],
                               const TacsScalar alpha[],
                               TacsScalar fdvSens[], int dvLen );
  void failure( const double pt[], const TacsScalar strain[],
                TacsScalar * fail );
  void failureStrainSens( const double pt[], const TacsScalar strain[],
                          TacsScalar sens[] );
  void addFailureDVSens( const double pt[], const TacsScalar strain[],
                         TacsScalar alpha, TacsScalar dvSens[], int dvLen );
  TacsScalar getDVOutputValue( int dvIndex, const double pt[] );

 public:
  // Material property information
  TacsScalar rho, D, G, nu, ys;
  double q, eps;

  // Maximum number of weights/nodes
  int nweights[4];

  // weights/index values
  class Weight {
  public:
    int index;
    double weight;
  };

  // Weights used in the computation
  Weight *w[4];

  // The design variable values
  TacsScalar xw[4];
};

#endif // PLANE_STRESS_TOPO_H
