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

#include "TMRMatrixFilterModel.h"

TMRQuadMatrixModel::TMRQuadMatrixModel(){}

int TMRQuadMatrixModel::getSpatialDim(){
  return 2;
}

int TMRQuadMatrixModel::getVarsPerNode(){
  return 1;
}

void TMRQuadMatrixModel::evalWeakIntegrand( int elemIndex, const double time,
                                            int n, const double pt[],
                                            const TacsScalar X[], const TacsScalar Xd[],
                                            const TacsScalar Ut[], const TacsScalar Ux[],
                                            TacsScalar DUt[], TacsScalar DUx[] ){
  DUt[0] = DUt[1] = DUt[2] = 0.0;
  DUx[0] = DUx[1] = 0.0;
  DUt[0] = Ut[0];
}

void TMRQuadMatrixModel::evalWeakJacobian( int elemIndex, const double time,
                                           int n, const double pt[],
                                           const TacsScalar X[], const TacsScalar Xd[],
                                           const TacsScalar Ut[], const TacsScalar Ux[],
                                           TacsScalar DUt[], TacsScalar DUx[],
                                           int *Jac_nnz, const int *Jac_pairs[],
                                           TacsScalar Jac[] ){
  DUt[0] = DUt[1] = DUt[2] = 0.0;
  DUx[0] = DUx[1] = 0.0;
  DUt[0] = Ut[0];
  *Jac_nnz = 1;
  *Jac_pairs = elem_Jac_pairs;
  Jac[0] = 1.0;
}

int TMRQuadMatrixModel::elem_Jac_pairs[] = {0, 0};

TMRHexaMatrixModel::TMRHexaMatrixModel(){}

int TMRHexaMatrixModel::getSpatialDim(){
  return 3;
}

int TMRHexaMatrixModel::getVarsPerNode(){
  return 1;
}

void TMRHexaMatrixModel::evalWeakIntegrand( int elemIndex, const double time,
                                            int n, const double pt[],
                                            const TacsScalar X[], const TacsScalar Xd[],
                                            const TacsScalar Ut[], const TacsScalar Ux[],
                                            TacsScalar DUt[], TacsScalar DUx[] ){
  DUt[0] = DUt[1] = DUt[2] = 0.0;
  DUx[0] = DUx[1] = DUx[2] = 0.0;
  DUt[0] = Ut[0];
}

void TMRHexaMatrixModel::evalWeakJacobian( int elemIndex, const double time,
                                           int n, const double pt[],
                                           const TacsScalar X[], const TacsScalar Xd[],
                                           const TacsScalar Ut[], const TacsScalar Ux[],
                                           TacsScalar DUt[], TacsScalar DUx[],
                                           int *Jac_nnz, const int *Jac_pairs[],
                                           TacsScalar Jac[] ){
  DUt[0] = DUt[1] = DUt[2] = 0.0;
  DUx[0] = DUx[1] = DUx[2] = 0.0;
  DUt[0] = Ut[0];
  *Jac_nnz = 1;
  *Jac_pairs = elem_Jac_pairs;
  Jac[0] = 1.0;
}

int TMRHexaMatrixModel::elem_Jac_pairs[] = {0, 0};
