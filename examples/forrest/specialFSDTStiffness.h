#ifndef SPECIAL_FSDT_STIFFNESS_H
#define SPECIAL_FSDT_STIFFNESS_H

#include "FSDTStiffness.h"
#include "MaterialProperties.h"

/*
  An FSDTStiffness class of it's own! 

  This class is specifically for testing the cylinder - it's required
  so that the failure criterion is evaluated only at the upper surface
  - and does not return the max of the upper/lower values.
*/
class specialFSDTStiffness : public FSDTStiffness {
 public:
  specialFSDTStiffness( OrthoPly * _ply, 
                        int _orthotropic_flag, 
                        TacsScalar _t, TacsScalar _kcorr ){
    ply = _ply; 
    ply->incref();
    
    orthotropic_flag = _orthotropic_flag;
    t = _t;
    kcorr = _kcorr;
  }
  ~specialFSDTStiffness(){
    ply->decref();
  }

  void getPointwiseMass( const double pt[], 
                         TacsScalar mass[] ){
    mass[0] = ply->getRho()*t;
    mass[1] = ply->getRho()*(t*t*t)/12.0;
  } 

  // Compute the stiffness
  TacsScalar getStiffness( const double pt[],
                           TacsScalar A[], TacsScalar B[],
                           TacsScalar D[], TacsScalar As[] ){
    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = 0.0;
    B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;
    D[0] = D[1] = D[2] = D[3] = D[4] = D[5] = 0.0;
    As[0] = As[1] = As[2] = 0.0;

    TacsScalar Qbar[6], Abar[3];
    ply->calculateQbar(Qbar, 0.0);
    ply->calculateAbar(Abar, 0.0);

    TacsScalar t0 = -0.5*t;
    TacsScalar t1 = 0.5*t;

    TacsScalar a = (t1 - t0);
    TacsScalar b = 0.5*(t1*t1 - t0*t0);
    TacsScalar d = 1.0/3.0*(t1*t1*t1 - t0*t0*t0);
    
    for ( int i = 0; i < 6; i++ ){
      A[i] += a*Qbar[i];
      B[i] += b*Qbar[i];
      D[i] += d*Qbar[i];
    }
    
    for ( int i = 0; i < 3; i++ ){
      As[i] += kcorr*a*Abar[i];
    }

    return 0.5*DRILLING_REGULARIZATION*(As[0] + As[2]);
  }

  // Compute the failure function at the given point
  void failure( const double pt[], const TacsScalar strain[], 
                TacsScalar * fail ){
    TacsScalar z = 0.5*t;
    
    TacsScalar e[3];
    e[0] = strain[0] + z*strain[3];
    e[1] = strain[1] + z*strain[4];
    e[2] = strain[2] + z*strain[5];

    if (orthotropic_flag){
      *fail = ply->failure(0.0, e);
    }
    else {
      *fail = sqrt(ply->failure(0.0, e));
    }
  }

  void failureStrainSens( const double pt[], const TacsScalar strain[],
                          TacsScalar sens[] ){
    TacsScalar z = 0.5*t;
    
    TacsScalar e[3], eSens[3];
    e[0] = strain[0] + z*strain[3];
    e[1] = strain[1] + z*strain[4];
    e[2] = strain[2] + z*strain[5];
    
    ply->failureStrainSens(eSens, 0.0, e);
    
    TacsScalar scale = 1.0;
    if (!orthotropic_flag){
      TacsScalar fail = sqrt(ply->failure(0.0, e));
      if (fail == 0.0){
        scale = 0.0;
      }
      else {
        scale = 0.5/fail;
      }
    }
    
    // Set the output sensitivity
    sens[0] = scale*eSens[0];
    sens[1] = scale*eSens[1];
    sens[2] = scale*eSens[2];
    sens[3] = scale*z*eSens[0];
    sens[4] = scale*z*eSens[1];
    sens[5] = scale*z*eSens[2];
    sens[6] = 0.0;
    sens[7] = 0.0;
  }


 private:
  OrthoPly * ply;
  TacsScalar t, kcorr;
  int orthotropic_flag;
};

#endif // SPECIAL_FSDT_STIFFNESS_H
