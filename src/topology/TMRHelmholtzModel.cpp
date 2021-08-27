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

#include "TMRHelmholtzModel.h"

TMRQuadHelmholtzModel::TMRQuadHelmholtzModel( double _r ){
  r = _r;
}

int TMRQuadHelmholtzModel::getNumParameters(){
  return 2;
}

int TMRQuadHelmholtzModel::getVarsPerNode(){
  return 1;
}

void TMRQuadHelmholtzModel::evalWeakIntegrand( int elemIndex, const double time,
                                               int n, const double pt[],
                                               const TacsScalar X[], const TacsScalar Xd[],
                                               const TacsScalar Ut[], const TacsScalar Ux[],
                                               TacsScalar DUt[], TacsScalar DUx[] ){
  DUt[0] = Ut[0];
  DUt[1] = DUt[2] = 0.0;
  DUx[0] = r*r*Ux[0];
  DUx[1] = r*r*Ux[1];
}

void TMRQuadHelmholtzModel::getWeakMatrixNonzeros( ElementMatrixType matType,
                                                   int elemIndex,
                                                   int *Jac_nnz,
                                                   const int *Jac_pairs[] ){
  *Jac_nnz = 3;
  *Jac_pairs = elem_Jac_pairs;
}

void TMRQuadHelmholtzModel::evalWeakMatrix( ElementMatrixType matType,
                                            int elemIndex, const double time,
                                            int n, const double pt[],
                                            const TacsScalar X[], const TacsScalar Xd[],
                                            const TacsScalar Ut[], const TacsScalar Ux[],
                                            TacsScalar DUt[], TacsScalar DUx[],
                                            TacsScalar Jac[] ){
  DUt[0] = Ut[0];
  DUt[1] = DUt[2] = 0.0;
  DUx[0] = r*r*Ux[0];
  DUx[1] = r*r*Ux[1];
  Jac[0] = 1.0;
  Jac[1] = Jac[2] = r*r;
}

int TMRQuadHelmholtzModel::elem_Jac_pairs[] = {0, 0, 3, 3, 4, 4};

TMRHexaHelmholtzModel::TMRHexaHelmholtzModel( double _r ){
  r = _r;
}

int TMRHexaHelmholtzModel::getNumParameters(){
  return 3;
}

int TMRHexaHelmholtzModel::getVarsPerNode(){
  return 1;
}

void TMRHexaHelmholtzModel::evalWeakIntegrand( int elemIndex, const double time,
                                               int n, const double pt[],
                                               const TacsScalar X[], const TacsScalar Xd[],
                                               const TacsScalar Ut[], const TacsScalar Ux[],
                                               TacsScalar DUt[], TacsScalar DUx[] ){
  DUt[0] = Ut[0];
  DUt[1] = DUt[2] = 0.0;
  DUx[0] = r*r*Ux[0];
  DUx[1] = r*r*Ux[1];
  DUx[2] = r*r*Ux[2];
}

void TMRHexaHelmholtzModel::getWeakMatrixNonzeros( ElementMatrixType matType,
                                                   int elemIndex,
                                                   int *Jac_nnz,
                                                   const int *Jac_pairs[] ){
  *Jac_nnz = 4;
  *Jac_pairs = elem_Jac_pairs;
}

void TMRHexaHelmholtzModel::evalWeakMatrix( ElementMatrixType matType,
                                            int elemIndex, const double time,
                                            int n, const double pt[],
                                            const TacsScalar X[], const TacsScalar Xd[],
                                            const TacsScalar Ut[], const TacsScalar Ux[],
                                            TacsScalar DUt[], TacsScalar DUx[],
                                            TacsScalar Jac[] ){
  DUt[0] = Ut[0];
  DUt[1] = DUt[2] = 0.0;
  DUx[0] = r*r*Ux[0];
  DUx[1] = r*r*Ux[1];
  DUx[2] = r*r*Ux[2];
  Jac[0] = 1.0;
  Jac[1] = Jac[2] = Jac[3] = r*r;
}

int TMRHexaHelmholtzModel::elem_Jac_pairs[] = {0, 0, 3, 3, 4, 4, 5, 5};
