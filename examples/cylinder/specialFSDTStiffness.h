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
  specialFSDTStiffness(OrthoPly* _ply, int _orthotropic_flag, TacsScalar _t,
                       TacsScalar _kcorr, int _tNum = -1) {
    ply = _ply;
    ply->incref();

    orthotropic_flag = _orthotropic_flag;
    t = _t;
    kcorr = _kcorr;
    tNum = _tNum;
  }
  ~specialFSDTStiffness() { ply->decref(); }

  void setDesignVars(const TacsScalar dvs[], int numDVs) {
    if (tNum >= 0 && tNum < numDVs) {
      t = dvs[tNum];
    }
  }
  void getDesignVars(TacsScalar dvs[], int numDVs) {
    if (tNum >= 0 && tNum < numDVs) {
      dvs[tNum] = t;
    }
  }

  void getPointwiseMass(const double pt[], TacsScalar mass[]) {
    mass[0] = ply->getRho() * t;
    mass[1] = ply->getRho() * (t * t * t) / 12.0;
  }

  // Compute the stiffness
  TacsScalar getStiffness(const double pt[], TacsScalar A[], TacsScalar B[],
                          TacsScalar D[], TacsScalar As[]) {
    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = 0.0;
    B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;
    D[0] = D[1] = D[2] = D[3] = D[4] = D[5] = 0.0;
    As[0] = As[1] = As[2] = 0.0;

    TacsScalar Qbar[6], Abar[3];
    ply->calculateQbar(Qbar, 0.0);
    ply->calculateAbar(Abar, 0.0);

    TacsScalar t0 = -0.5 * t;
    TacsScalar t1 = 0.5 * t;

    TacsScalar a = (t1 - t0);
    TacsScalar d = 1.0 / 3.0 * (t1 * t1 * t1 - t0 * t0 * t0);

    for (int i = 0; i < 6; i++) {
      A[i] += a * Qbar[i];
      D[i] += d * Qbar[i];
    }

    for (int i = 0; i < 3; i++) {
      As[i] += kcorr * a * Abar[i];
    }

    return 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);
  }

  // Compute the failure function at the given point
  void failure(const double pt[], const TacsScalar strain[], TacsScalar* fail) {
    // The failure criteria evaluated at the top/bottom surface
    TacsScalar f1, f2;

    // Compute the strain at the top surface
    TacsScalar z1 = 0.5 * t;
    TacsScalar e1[3];
    e1[0] = strain[0] + z1 * strain[3];
    e1[1] = strain[1] + z1 * strain[4];
    e1[2] = strain[2] + z1 * strain[5];

    if (orthotropic_flag) {
      f1 = ply->failure(0.0, e1);
    } else {
      f1 = sqrt(ply->failure(0.0, e1));
    }

    // Compute the strain at the bottom surface
    TacsScalar z2 = -0.5 * t;
    TacsScalar e2[3];
    e2[0] = strain[0] + z2 * strain[3];
    e2[1] = strain[1] + z2 * strain[4];
    e2[2] = strain[2] + z2 * strain[5];

    if (orthotropic_flag) {
      f2 = ply->failure(0.0, e2);
    } else {
      f2 = sqrt(ply->failure(0.0, e2));
    }

    // Compute the weighted combination of the two
    double P = 10.0;

    if (f1 > f2) {
      *fail = f1 + log(1.0 + exp(P * (f2 - f1))) / P;
    } else {
      *fail = f2 + log(1.0 + exp(P * (f1 - f2))) / P;
    }
  }

  void failureStrainSens(const double pt[], const TacsScalar strain[],
                         TacsScalar sens[]) {
    // The failure criteria evaluated at the top/bottom surface
    TacsScalar f1, f2;

    // Compute the strain at the top surface
    TacsScalar z1 = 0.5 * t;
    TacsScalar e1[3], de1[3];
    e1[0] = strain[0] + z1 * strain[3];
    e1[1] = strain[1] + z1 * strain[4];
    e1[2] = strain[2] + z1 * strain[5];

    if (orthotropic_flag) {
      f1 = ply->failure(0.0, e1);
      ply->failureStrainSens(de1, 0.0, e1);
    } else {
      f1 = sqrt(ply->failure(0.0, e1));
      ply->failureStrainSens(de1, 0.0, e1);
      TacsScalar scale = 0.5 / f1;
      de1[0] *= scale;
      de1[1] *= scale;
      de1[2] *= scale;
    }

    // Compute the strain at the bottom surface
    TacsScalar z2 = -0.5 * t;
    TacsScalar e2[3], de2[3];
    e2[0] = strain[0] + z2 * strain[3];
    e2[1] = strain[1] + z2 * strain[4];
    e2[2] = strain[2] + z2 * strain[5];

    if (orthotropic_flag) {
      f2 = ply->failure(0.0, e2);
      ply->failureStrainSens(de2, 0.0, e2);
    } else {
      f2 = sqrt(ply->failure(0.0, e2));
      ply->failureStrainSens(de2, 0.0, e2);
      TacsScalar scale = 0.5 / f2;
      de2[0] *= scale;
      de2[1] *= scale;
      de2[2] *= scale;
    }

    // Compute the weighted combination of the two
    double P = 10.0;

    // Compute the factors on the scaled sensitivities
    TacsScalar d1 = 0.0, d2 = 0.0;
    if (f1 > f2) {
      TacsScalar expf = exp(P * (f2 - f1));
      d1 = 1.0 / (1.0 + expf);
      d2 = expf / (1.0 + expf);
    } else {
      TacsScalar expf = exp(P * (f1 - f2));
      d1 = expf / (1.0 + expf);
      d2 = 1.0 / (1.0 + expf);
    }

    // Set the output sensitivity
    sens[0] = d1 * de1[0] + d2 * de2[0];
    sens[1] = d1 * de1[1] + d2 * de2[1];
    sens[2] = d1 * de1[2] + d2 * de2[2];
    sens[3] = z1 * d1 * de1[0] + z2 * d2 * de2[0];
    sens[4] = z1 * d1 * de1[1] + z2 * d2 * de2[1];
    sens[5] = z1 * d1 * de1[2] + z2 * d2 * de2[2];
    sens[6] = 0.0;
    sens[7] = 0.0;
  }

  TacsScalar getDVOutputValue(int dv_index, const double pt[]) {
    if (dv_index == 0) {
      return t;
    }
    return 0.0;
  }

 private:
  OrthoPly* ply;
  int tNum;
  TacsScalar t, kcorr;
  int orthotropic_flag;
};

#endif  // SPECIAL_FSDT_STIFFNESS_H
