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

#ifndef TMR_QUADRANT_STIFFNESS_H
#define TMR_QUADRANT_STIFFNESS_H

/*
  The following class defines a SolidStiffness base class for topology
  optimization using TACS.
*/

#include "TACSMaterialProperties.h"
#include "TACSPlaneStressConstitutive.h"
#include "TMROctConstitutive.h"
#include "TMRQuadForest.h"

/*
  The TMROctStiffness class

  This defines the TMROctStiffness class which takes the weights from
  up to 8 adjacent vertices. This class uses the RAMP method for
  penalization.
*/
class TMRQuadConstitutive : public TACSPlaneStressConstitutive {
 public:
  TMRQuadConstitutive(TMRStiffnessProperties *_props, TMRQuadForest *_forest);
  ~TMRQuadConstitutive();

  // Get the number of design variables at each "design node"
  int getDesignVarsPerNode() { return nvars; }

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Evaluate the material density
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]);

  // Add the derivative of the density
  void addDensityDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the material density
  TacsScalar evalMassMatrixDensity(int elemIndex, const double pt[],
                                   const TacsScalar X[]);

  // Add the derivative of the density
  void addMassMatrixDensityDVSens(int elemIndex, TacsScalar scale,
                                  const double pt[], const TacsScalar X[],
                                  int dvLen, TacsScalar dfdx[]);

  // Evaluate the specific heat
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]);

  // Add the derivative of the density
  void addSpecificHeatDVSens(int elemIndex, TacsScalar scale, const double pt[],
                             const TacsScalar X[], int dvLen,
                             TacsScalar dfdx[]);

  // Evaluate the stresss
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar strain[], TacsScalar stress[]);

  // Evaluate the tangent stiffness
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]);

  // Add the contribution
  void addStressDVSens(int elemIndex, TacsScalar scale, const double pt[],
                       const TacsScalar X[], const TacsScalar strain[],
                       const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the geometric stiffness constitutive matrix
  void evalGeometricTangentStiffness(int elemIndex, const double pt[],
                                     const TacsScalar X[], TacsScalar C[]);

  // Add the derivative of the geometric constitutive matrix
  void addGeometricTangentStressDVSens(int elemIndex, TacsScalar scale,
                                       const double pt[], const TacsScalar X[],
                                       const TacsScalar e[],
                                       const TacsScalar psi[], int dvLen,
                                       TacsScalar dfdx[]);

  // Evaluate the thermal strain
  void evalThermalStrain(int elemIndex, const double pt[], const TacsScalar X[],
                         TacsScalar theta, TacsScalar strain[]);

  // Evaluate the thermal strain sensitivity
  void addThermalStrainDVSens(int elemIndex, const double pt[],
                              const TacsScalar X[], TacsScalar theta,
                              const TacsScalar psi[], int dvLen,
                              TacsScalar dfdx[]);

  // Evaluate the heat flux, given the thermal gradient
  void evalHeatFlux(int elemIndex, const double pt[], const TacsScalar X[],
                    const TacsScalar grad[], TacsScalar flux[]);

  // Evaluate the tangent of the heat flux
  void evalTangentHeatFlux(int elemIndex, const double pt[],
                           const TacsScalar X[], TacsScalar C[]);

  // Add the derivative of the heat flux
  void addHeatFluxDVSens(int elemIndex, TacsScalar scale, const double pt[],
                         const TacsScalar X[], const TacsScalar grad[],
                         const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the material failure index
  TacsScalar evalFailure(int elemIndex, const double pt[], const TacsScalar X[],
                         const TacsScalar e[]);

  // Add the derivative of the function w.r.t. design variables
  void addFailureDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], const TacsScalar e[], int dvLen,
                        TacsScalar dfdx[]);

  // Evaluate the derivative of the failure criteria w.r.t. strain
  TacsScalar evalFailureStrainSens(int elemIndex, const double pt[],
                                   const TacsScalar X[], const TacsScalar e[],
                                   TacsScalar sens[]);

  // Evaluate the design field value
  TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                  const TacsScalar X[], int index);

  // Extra info about the constitutive class
  const char *getObjectName();

 private:
  // The stiffness properties
  TMRStiffnessProperties *props;
  TMRQuadForest *forest;

  // Information about the design variable values
  int nmats, nvars;
  TacsScalar *x;           // All the design variable values
  double *N;               // Space for the shape functions
  TacsScalar *temp_array;  // Temporary array
};

#endif  // TMR_QUADRANT_STIFFNESS_H
