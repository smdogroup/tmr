#include "PSTopo.h"

#include <math.h>

#include "tacslapack.h"

PSTopo::PSTopo(TacsScalar _rho, TacsScalar _E, TacsScalar _nu, TacsScalar _ys,
               double _q, double _eps, int _nodes[], double _weights[],
               int _nweights) {
  // Set the material properties
  rho = _rho;
  D = _E / (1.0 - _nu * _nu);
  G = 0.5 * _E / (1.0 + _nu);
  nu = _nu;
  ys = _ys;
  q = _q;
  eps = _eps;

  // Set the weights
  nweights = _nweights;
  w = new Weight[nweights];
  for (int i = 0; i < nweights; i++) {
    w[i].index = _nodes[i];
    w[i].weight = _weights[i];
  }

  // Set the design variable value
  xw = 0.5;
}

PSTopo::~PSTopo() { delete[] w; }

/*
  Set the design variable values
*/
void PSTopo::setDesignVars(const TacsScalar dvs[], int numDVs) {
  // Record the design variable values
  xw = 0.0;
  for (int i = 0; i < nweights; i++) {
    xw += w[i].weight * dvs[w[i].index];
  }
}

/*
  Retrieve the design variable values. This call has no effect...
*/
void PSTopo::getDesignVars(TacsScalar dvs[], int numDVs) {
  for (int i = 0; i < nweights; i++) {
    if (w[i].index >= 0 && w[i].index < numDVs) {
      dvs[w[i].index] = 0.95;
    }
  }
}

/*
  Get the design variable range
*/
void PSTopo::getDesignVarRange(TacsScalar lb[], TacsScalar ub[], int numDVs) {
  // number of variables per node
  for (int i = 0; i < nweights; i++) {
    if (w[i].index >= 0 && w[i].index < numDVs) {
      lb[w[i].index] = 0.0;
      ub[w[i].index] = 1.0;
    }
  }
}

/*
  Compute the stress at a parametric point in the element based on
  the local strain value
*/
void PSTopo::calculateStress(const double pt[], const TacsScalar e[],
                             TacsScalar s[]) {
  // Compute the penalty
  // double p = xw/(1.0 + q*(1.0 - xw)) + k0;
  // double p = xw*xw*xw;
  double p = pow(xw, q);
  s[0] = p * D * (e[0] + nu * e[1]);
  s[1] = p * D * (e[1] + nu * e[0]);
  s[2] = p * G * e[2];
}

/*
  Add the derivative of the product of the strain with an input vector
  psi and add it to the array fdvSens
*/
void PSTopo::addStressDVSens(const double pt[], const TacsScalar e[],
                             TacsScalar alpha, const TacsScalar psi[],
                             TacsScalar fdvSens[], int dvLen) {
  // Compute the contribution to the derivative
  // double p = (q + 1.0)/((1.0 + q*(1.0 - xw))*(1.0 + q*(1.0 - xw)));
  double p = q * pow(xw, q - 1.0);
  TacsScalar s[3];
  s[0] = p * D * (e[0] + nu * e[1]);
  s[1] = p * D * (e[1] + nu * e[0]);
  s[2] = p * G * e[2];

  TacsScalar scale = alpha * (psi[0] * s[0] + psi[1] * s[1] + psi[2] * s[2]);

  // Add the dependency on the filter
  for (int i = 0; i < nweights; i++) {
    fdvSens[w[i].index] += scale * w[i].weight;
  }
}

/*
  Compute the mass at this point
*/
void PSTopo::getPointwiseMass(const double pt[], TacsScalar mass[]) {
  mass[0] = rho * xw;
}

/*
  Add the derivative of the mass at this point
*/
void PSTopo::addPointwiseMassDVSens(const double pt[], const TacsScalar alpha[],
                                    TacsScalar fdvSens[], int dvLen) {
  for (int i = 0; i < nweights; i++) {
    fdvSens[w[i].index] += w[i].weight * alpha[0] * rho;
  }
}

void PSTopo::failure(const double pt[], const TacsScalar e[],
                     TacsScalar *fail) {
  TacsScalar r = xw / (eps * (1.0 - xw) + xw);
  TacsScalar s[3];
  s[0] = D * (e[0] + nu * e[1]);
  s[1] = D * (e[1] + nu * e[0]);
  s[2] = G * e[2];
  *fail = r *
          sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3.0 * s[2] * s[2]) /
          ys;
}

void PSTopo::failureStrainSens(const double pt[], const TacsScalar e[],
                               TacsScalar sens[]) {
  TacsScalar r = xw / (eps * (1.0 - xw) + xw);
  TacsScalar s[3];
  s[0] = D * (e[0] + nu * e[1]);
  s[1] = D * (e[1] + nu * e[0]);
  s[2] = G * e[2];

  TacsScalar fact =
      sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3.0 * s[2] * s[2]);
  if (fact != 0.0) {
    fact = 1.0 / (ys * fact);
  }

  TacsScalar ds[3];
  ds[0] = r * fact * (s[0] - 0.5 * s[1]);
  ds[1] = r * fact * (s[1] - 0.5 * s[0]);
  ds[2] = 3.0 * r * fact * s[2];

  sens[0] = D * (ds[0] + nu * ds[1]);
  sens[1] = D * (ds[1] + nu * ds[0]);
  sens[2] = G * ds[2];
}

void PSTopo::addFailureDVSens(const double pt[], const TacsScalar e[],
                              TacsScalar alpha, TacsScalar fdvSens[],
                              int dvLen) {
  TacsScalar d = 1.0 / (eps * (1.0 - xw) + xw);
  TacsScalar r = eps * d * d;

  TacsScalar s[3];
  s[0] = D * (e[0] + nu * e[1]);
  s[1] = D * (e[1] + nu * e[0]);
  s[2] = G * e[2];
  TacsScalar fail =
      sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3.0 * s[2] * s[2]) / ys;

  TacsScalar scale = alpha * r * fail;
  for (int i = 0; i < nweights; i++) {
    fdvSens[w[i].index] += scale * w[i].weight;
  }
}

PSTopo4::PSTopo4(TacsScalar _rho, TacsScalar _E, TacsScalar _nu, TacsScalar _ys,
                 double _q, double _eps, int *_nodes[], double *_weights[],
                 int _nweights[]) {
  // Set the material properties
  rho = _rho;
  D = _E / (1.0 - _nu * _nu);
  G = 0.5 * _E / (1.0 + _nu);
  nu = _nu;
  ys = _ys;
  q = _q;
  eps = _eps;

  // Set the weights
  for (int k = 0; k < 4; k++) {
    nweights[k] = _nweights[k];
    w[k] = new Weight[nweights[k]];
    for (int i = 0; i < nweights[k]; i++) {
      w[k][i].index = _nodes[k][i];
      w[k][i].weight = _weights[k][i];
    }

    // Set the design variable value
    xw[k] = 0.5;
  }
}

PSTopo4::~PSTopo4() {
  for (int k = 0; k < 4; k++) {
    delete[] w[k];
  }
}

/*
  Set the design variable values
*/
void PSTopo4::setDesignVars(const TacsScalar dvs[], int numDVs) {
  // Record the design variable values
  for (int k = 0; k < 4; k++) {
    xw[k] = 0.0;
    for (int i = 0; i < nweights[k]; i++) {
      xw[k] += w[k][i].weight * dvs[w[k][i].index];
    }
  }
}

/*
  Retrieve the design variable values. This call has no effect...
*/
void PSTopo4::getDesignVars(TacsScalar dvs[], int numDVs) {
  for (int k = 0; k < 4; k++) {
    for (int i = 0; i < nweights[k]; i++) {
      if (w[k][i].index >= 0 && w[k][i].index < numDVs) {
        dvs[w[k][i].index] = 0.95;
      }
    }
  }
}

/*
  Get the design variable range
*/
void PSTopo4::getDesignVarRange(TacsScalar lb[], TacsScalar ub[], int numDVs) {
  // number of variables per node
  for (int k = 0; k < 4; k++) {
    for (int i = 0; i < nweights[k]; i++) {
      if (w[k][i].index >= 0 && w[k][i].index < numDVs) {
        lb[w[k][i].index] = 0.0;
        ub[w[k][i].index] = 1.0;
      }
    }
  }
}

/*
  Compute the stress at a parametric point in the element based on
  the local strain value
*/
void PSTopo4::calculateStress(const double pt[], const TacsScalar e[],
                              TacsScalar s[]) {
  // Compute the penalty
  // double p = xw/(1.0 + q*(1.0 - xw)) + k0;
  // double p = xw*xw*xw;
  double N[4];
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  // Compute the local value of the design variable
  TacsScalar xp = (N[0] * xw[0] + N[1] * xw[1] + N[2] * xw[2] + N[3] * xw[3]);

  // Compute the stress
  double p = pow(xp, q);
  s[0] = p * D * (e[0] + nu * e[1]);
  s[1] = p * D * (e[1] + nu * e[0]);
  s[2] = p * G * e[2];
}

/*
  Add the derivative of the product of the strain with an input vector
  psi and add it to the array fdvSens
*/
void PSTopo4::addStressDVSens(const double pt[], const TacsScalar e[],
                              TacsScalar alpha, const TacsScalar psi[],
                              TacsScalar fdvSens[], int dvLen) {
  // Compute the contribution to the derivative
  // double p = (q + 1.0)/((1.0 + q*(1.0 - xw))*(1.0 + q*(1.0 - xw)));
  double N[4];
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  // Compute the local value of the design variable
  TacsScalar xp = (N[0] * xw[0] + N[1] * xw[1] + N[2] * xw[2] + N[3] * xw[3]);

  // Compute the stress contribution
  double p = q * pow(xp, q - 1.0);
  TacsScalar s[3];
  s[0] = p * D * (e[0] + nu * e[1]);
  s[1] = p * D * (e[1] + nu * e[0]);
  s[2] = p * G * e[2];

  TacsScalar scale = alpha * (psi[0] * s[0] + psi[1] * s[1] + psi[2] * s[2]);

  // Add the dependency on the filter
  for (int k = 0; k < 4; k++) {
    for (int i = 0; i < nweights[k]; i++) {
      fdvSens[w[k][i].index] += scale * N[k] * w[k][i].weight;
    }
  }
}

/*
  Compute the mass at this point
*/
void PSTopo4::getPointwiseMass(const double pt[], TacsScalar mass[]) {
  double N[4];
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  // Compute the local
  TacsScalar xp = (N[0] * xw[0] + N[1] * xw[1] + N[2] * xw[2] + N[3] * xw[3]);

  mass[0] = rho * xp;
}

/*
  Add the derivative of the mass at this point
*/
void PSTopo4::addPointwiseMassDVSens(const double pt[],
                                     const TacsScalar alpha[],
                                     TacsScalar fdvSens[], int dvLen) {
  double N[4];
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  for (int k = 0; k < 4; k++) {
    for (int i = 0; i < nweights[k]; i++) {
      fdvSens[w[k][i].index] += w[k][i].weight * alpha[0] * rho * N[k];
    }
  }
}

void PSTopo4::failure(const double pt[], const TacsScalar e[],
                      TacsScalar *fail) {
  double N[4];
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  // Compute the local value of the design variable
  TacsScalar xp = (N[0] * xw[0] + N[1] * xw[1] + N[2] * xw[2] + N[3] * xw[3]);

  TacsScalar r = xp / (eps * (1.0 - xp) + xp);
  TacsScalar s[3];
  s[0] = D * (e[0] + nu * e[1]);
  s[1] = D * (e[1] + nu * e[0]);
  s[2] = G * e[2];
  *fail = r *
          sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3.0 * s[2] * s[2]) /
          ys;
}

void PSTopo4::failureStrainSens(const double pt[], const TacsScalar e[],
                                TacsScalar sens[]) {
  double N[4];
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  // Compute the local value of the design variable
  TacsScalar xp = (N[0] * xw[0] + N[1] * xw[1] + N[2] * xw[2] + N[3] * xw[3]);

  TacsScalar r = xp / (eps * (1.0 - xp) + xp);

  TacsScalar s[3];
  s[0] = D * (e[0] + nu * e[1]);
  s[1] = D * (e[1] + nu * e[0]);
  s[2] = G * e[2];

  TacsScalar fact =
      sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3.0 * s[2] * s[2]);
  if (fact != 0.0) {
    fact = 1.0 / (ys * fact);
  }

  TacsScalar ds[3];
  ds[0] = r * fact * (s[0] - 0.5 * s[1]);
  ds[1] = r * fact * (s[1] - 0.5 * s[0]);
  ds[2] = 3.0 * r * fact * s[2];

  sens[0] = D * (ds[0] + nu * ds[1]);
  sens[1] = D * (ds[1] + nu * ds[0]);
  sens[2] = G * ds[2];
}

void PSTopo4::addFailureDVSens(const double pt[], const TacsScalar e[],
                               TacsScalar alpha, TacsScalar fdvSens[],
                               int dvLen) {
  double N[4];
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  // Compute the local value of the design variable
  TacsScalar xp = (N[0] * xw[0] + N[1] * xw[1] + N[2] * xw[2] + N[3] * xw[3]);
  TacsScalar d = 1.0 / (eps * (1.0 - xp) + xp);
  TacsScalar r = eps * d * d;

  TacsScalar s[3];
  s[0] = D * (e[0] + nu * e[1]);
  s[1] = D * (e[1] + nu * e[0]);
  s[2] = G * e[2];
  TacsScalar fail =
      sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3.0 * s[2] * s[2]) / ys;

  TacsScalar scale = alpha * r * fail;
  for (int k = 0; k < 4; k++) {
    for (int i = 0; i < nweights[k]; i++) {
      fdvSens[w[k][i].index] += scale * w[k][i].weight * N[k];
    }
  }
}

TacsScalar PSTopo4::getDVOutputValue(int dvIndex, const double pt[]) {
  double N[4];
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  TacsScalar xp = (N[0] * xw[0] + N[1] * xw[1] + N[2] * xw[2] + N[3] * xw[3]);

  return xp;
}
