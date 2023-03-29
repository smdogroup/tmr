#ifndef TMR_KS_FSDT_H
#define TMR_KS_FSDT_H

#include "FSDTStiffness.h"
#include "YSlibrary.h"

class ksFSDTStiffness : public FSDTStiffness {
 public:
  static const int NUM_STRESSES = FSDTStiffness::NUM_STRESSES;

  ksFSDTStiffness(TacsScalar _ksweight, TacsScalar _rho, TacsScalar _E,
                  TacsScalar _nu, TacsScalar _kcorr, TacsScalar _yieldStress,
                  TacsScalar _thickness, int _tNum = -1,
                  TacsScalar _minThickness = 1e-3,
                  TacsScalar _maxThickness = 50.0) {
    // Set the material properties
    ksweight = _ksweight;
    rho = _rho;
    E = _E;
    nu = _nu;
    G = 0.5 * E / (1.0 + nu);
    kcorr = _kcorr;
    yieldStress = _yieldStress;

    // Copy over the thickness data
    t = _thickness;
    tNum = _tNum;
    minThickness = _minThickness;
    maxThickness = _maxThickness;

    // Set the locator (if any is specified)
  }
  ~ksFSDTStiffness() {}

  // Functions for design variable control
  // -------------------------------------
  void setDesignVars(const TacsScalar dvs[], int dvLen) {
    if (tNum >= 0 && tNum < dvLen) {
      t = dvs[tNum];
    }
  }
  void getDesignVars(TacsScalar dvs[], int dvLen) {
    if (tNum >= 0 && tNum < dvLen) {
      dvs[tNum] = t;
    }
  }
  void getDesignVarRange(TacsScalar lowerBound[], TacsScalar upperBound[],
                         int dvLen) {
    if (tNum >= 0 && tNum < dvLen) {
      lowerBound[tNum] = minThickness;
      upperBound[tNum] = maxThickness;
    }
  }

  // Functions required by TACSConstitutive
  // --------------------------------------
  void getPointwiseMass(const double pt[], TacsScalar mass[]) {
    mass[0] = rho * t;
    mass[1] = (rho * t * t * t) / 12.0;
  }
  void addPointwiseMassDVSens(const double pt[], const TacsScalar alpha[],
                              TacsScalar fdvSens[], int dvLen) {
    if (tNum >= 0 && tNum < dvLen) {
      fdvSens[tNum] += rho * (alpha[0] + 0.25 * alpha[1] * t * t);
    }
  }

  // Functions required by FSDTStiffness
  // -----------------------------------
  TacsScalar getStiffness(const double pt[], TacsScalar A[], TacsScalar B[],
                          TacsScalar D[], TacsScalar As[]) {
    TacsScalar a = t * E / (1.0 - nu * nu);
    TacsScalar d = t * t * a / 12.0;

    // Compute the in-plane stiffness
    A[2] = A[4] = 0.0;
    A[0] = A[3] = a;
    A[1] = nu * a;
    A[5] = G * t;

    // Set the bending-stretching coupling to zero
    B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;

    // Compute the bending moments
    D[2] = D[4] = 0.0;
    D[0] = D[3] = d;
    D[1] = nu * d;
    D[5] = 0.5 * d * (1.0 - nu);

    // Compute the shear resultants
    As[0] = kcorr * G * t;
    As[1] = 0.0;
    As[2] = kcorr * G * t;

    return DRILLING_REGULARIZATION * G * t;
  }
  void addStiffnessDVSens(const double pt[], const TacsScalar e[],
                          const TacsScalar psi[], TacsScalar rotPsi,
                          TacsScalar fdvSens[], int dvLen) {
    if (tNum >= 0 && tNum < dvLen) {
      // Compute the derivative of the stiffness coefficients
      TacsScalar A = E / (1.0 - nu * nu);
      TacsScalar D = t * t * A / 4.0;

      // Store the derivative of the stress values
      TacsScalar s[8];

      // Compute the in-plane resultants
      s[0] = A * (e[0] + nu * e[1]);
      s[1] = A * (e[1] + nu * e[0]);
      s[2] = G * e[2];

      // Compute the bending moments
      s[3] = D * (e[3] + nu * e[4]);
      s[4] = D * (e[4] + nu * e[3]);
      s[5] = 0.5 * D * (1.0 - nu) * e[5];

      // Compute the shear resultants
      s[6] = kcorr * G * e[6];
      s[7] = kcorr * G * e[7];

      TacsScalar ksens = DRILLING_REGULARIZATION * G;

      // Add the result to the design variable vector
      fdvSens[tNum] += (s[0] * psi[0] + s[1] * psi[1] + s[2] * psi[2] +
                        s[3] * psi[3] + s[4] * psi[4] + s[5] * psi[5] +
                        s[6] * psi[6] + s[7] * psi[7] + rotPsi * ksens);
    }
  }

  // Functions to compute the failure properties
  // -------------------------------------------
  void failure(const double pt[], const TacsScalar strain[], TacsScalar *fail) {
    TacsScalar stress[3];
    TacsScalar ht = 0.5 * t;

    // Determine whether the failure will occur on the top or the bottom
    // Test the top of the plate
    calculatePlaneStress(stress, ht, strain);
    TacsScalar failTop = VonMisesFailurePlaneStress(stress, yieldStress);

    // Test the bottom of the plate
    calculatePlaneStress(stress, -ht, strain);
    TacsScalar failBot = VonMisesFailurePlaneStress(stress, yieldStress);

    // Compute the max of the top/bottom failure values
    TacsScalar fmax =
        (TacsRealPart(failTop) > TacsRealPart(failBot) ? failTop : failBot);

    // Compute the aggregate of the top/bottom
    *fail = (fmax + log(exp(ksweight * (failTop - fmax)) +
                        exp(ksweight * (failBot - fmax))) /
                        ksweight);
  }
  void failureStrainSens(const double pt[], const TacsScalar strain[],
                         TacsScalar sens[]) {
    TacsScalar stress[3];
    TacsScalar ht = 0.5 * t;

    // Determine whether the failure will occur on the top or the bottom
    // Test the top of the plate
    calculatePlaneStress(stress, ht, strain);
    TacsScalar failTop = VonMisesFailurePlaneStress(stress, yieldStress);

    // Test the bottom of the plate
    calculatePlaneStress(stress, -ht, strain);
    TacsScalar failBot = VonMisesFailurePlaneStress(stress, yieldStress);

    // Compute the max value
    TacsScalar fmax =
        (TacsRealPart(failTop) > TacsRealPart(failBot) ? failTop : failBot);

    // Top/bottom exponential values
    TacsScalar topExp = exp(ksweight * (failTop - fmax));
    TacsScalar botExp = exp(ksweight * (failBot - fmax));
    TacsScalar ksSum = topExp + botExp;

    // Compute the weights
    topExp = topExp / ksSum;
    botExp = botExp / ksSum;

    // The stress sensitivity of the top/bottom
    TacsScalar st[8], sb[8];

    // Test the top of the plate
    TacsScalar stressSens[3];
    calculatePlaneStress(stress, ht, strain);
    VonMisesFailurePlaneStressSens(stressSens, stress, yieldStress);
    calculatePlaneStressTranspose(st, ht, stressSens);

    // Test the bottom of the plate
    calculatePlaneStress(stress, -ht, strain);
    VonMisesFailurePlaneStressSens(stressSens, stress, yieldStress);
    calculatePlaneStressTranspose(sb, -ht, stressSens);

    for (int i = 0; i < 8; i++) {
      sens[i] = (topExp * st[i] + botExp * sb[i]);
    }
  }
  void addFailureDVSens(const double pt[], const TacsScalar strain[],
                        TacsScalar alpha, TacsScalar fdvSens[], int dvLen) {
    if (tNum >= 0 && tNum < dvLen) {
      TacsScalar stress[3];
      TacsScalar ht = 0.5 * t;

      // Determine whether the failure will occur on the top or the bottom
      // Test the top of the plate
      calculatePlaneStress(stress, ht, strain);
      TacsScalar failTop = VonMisesFailurePlaneStress(stress, yieldStress);

      // Test the bottom of the plate
      calculatePlaneStress(stress, -ht, strain);
      TacsScalar failBot = VonMisesFailurePlaneStress(stress, yieldStress);

      // Compute the max value
      TacsScalar fmax =
          (TacsRealPart(failTop) > TacsRealPart(failBot) ? failTop : failBot);

      // Top/bottom exponential values
      TacsScalar topExp = exp(ksweight * (failTop - fmax));
      TacsScalar botExp = exp(ksweight * (failBot - fmax));
      TacsScalar ksSum = topExp + botExp;

      // Compute the weights
      topExp = alpha * topExp / ksSum;
      botExp = alpha * botExp / ksSum;

      // Test the top of the plate
      TacsScalar stressSens[3];
      calculatePlaneStress(stress, ht, strain);
      VonMisesFailurePlaneStressSens(stressSens, stress, yieldStress);
      fdvSens[tNum] +=
          0.5 * topExp * (calculatePlaneStressTSensProduct(stressSens, strain));

      calculatePlaneStress(stress, -ht, strain);
      VonMisesFailurePlaneStressSens(stressSens, stress, yieldStress);
      fdvSens[tNum] += -0.5 * botExp *
                       (calculatePlaneStressTSensProduct(stressSens, strain));
    }
  }

  // Retrieve the design variable for plotting purposes
  // --------------------------------------------------
  TacsScalar getDVOutputValue(int dv_index, const double pt[]) {
    if (dv_index == 0) {
      return t;
    }
    if (dv_index == 1) {
      return tNum;
    }
    return 0.0;
  }

  // Return the name of the constitutive object
  // ------------------------------------------
  const char *constitutiveName() { return "ksFSDT"; }

 private:
  // Calculate the state of plane stress within the element
  // ------------------------------------------------------
  inline void calculatePlaneStress(TacsScalar pstr[], const TacsScalar ht,
                                   const TacsScalar strain[]) {
    const TacsScalar D = E / (1.0 - nu * nu);
    pstr[0] =
        D * ((strain[0] + ht * strain[3]) + nu * (strain[1] + ht * strain[4]));
    pstr[1] =
        D * (nu * (strain[0] + ht * strain[3]) + (strain[1] + ht * strain[4]));
    pstr[2] = G * (strain[2] + ht * strain[5]);
  }

  inline void calculatePlaneStressTranspose(TacsScalar strain[],
                                            const TacsScalar ht,
                                            const TacsScalar pstr[]) {
    const TacsScalar D = E / (1.0 - nu * nu);
    strain[0] = D * (pstr[0] + nu * pstr[1]);
    strain[1] = D * (nu * pstr[0] + pstr[1]);
    strain[2] = G * pstr[2];

    strain[3] = ht * D * (pstr[0] + nu * pstr[1]);
    strain[4] = ht * D * (nu * pstr[0] + pstr[1]);
    strain[5] = G * ht * pstr[2];

    strain[6] = 0.0;
    strain[7] = 0.0;
  }

  TacsScalar calculatePlaneStressTSensProduct(const TacsScalar pstr[],
                                              const TacsScalar strain[]) {
    return E / (1.0 - nu * nu) *
               (pstr[0] * (strain[3] + nu * strain[4]) +
                pstr[1] * (nu * strain[3] + strain[4])) +
           G * pstr[2] * strain[5];
  }

  // The material properties required for this object
  TacsScalar ksweight;
  TacsScalar kcorr;        // The shear correction factor
  TacsScalar rho;          // The material density
  TacsScalar E, nu, G;     // The stiffness properties
  TacsScalar yieldStress;  // The material yield stress

  // The thickness information required by this object
  int tNum;                               // The design variable number
  TacsScalar t;                           // The thickness
  TacsScalar minThickness, maxThickness;  // The min/max thickness values
};

#endif  // TMR_KS_FSDT_H
