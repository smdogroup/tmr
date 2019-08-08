#ifndef PLANE_STRESS_TOPO_H
#define PLANE_STRESS_TOPO_H

#include "TACSAssembler.h"
#include "PlaneStressStiffness.h"
#include "PlaneStressQuad.h"
#include "FElibrary.h"

class LobattoQuad2 : public PlaneStressQuad<2> {
 public:
 LobattoQuad2( PlaneStressStiffness *_stiff ):
  PlaneStressQuad<2>(_stiff){
    this->numGauss = 3;
    this->gaussWts = FElibrary::lobattoWts3;
    this->gaussPts = FElibrary::lobattoPts3;
  }  
};

class LobattoQuad3 : public PlaneStressQuad<3> {
 public:
 LobattoQuad3( PlaneStressStiffness *_stiff ):
  PlaneStressQuad<3>(_stiff){
    this->numGauss = 4;
    this->gaussWts = FElibrary::lobattoWts4;
    this->gaussPts = FElibrary::lobattoPts4;
  }  
};

class LobattoQuad4 : public PlaneStressQuad<4> {
 public:
 LobattoQuad4( PlaneStressStiffness *_stiff ):
  PlaneStressQuad<4>(_stiff){
    this->numGauss = 5;
    this->gaussWts = FElibrary::lobattoWts5;
    this->gaussPts = FElibrary::lobattoPts5;
  }  
};

class LobattoQuad5 : public PlaneStressQuad<5> {
 public:
 LobattoQuad5( PlaneStressStiffness *_stiff ):
  PlaneStressQuad<5>(_stiff){
    this->numGauss = 6;
    this->gaussWts = FElibrary::lobattoWts6;
    this->gaussPts = FElibrary::lobattoPts6;
  }  
};

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
