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

#ifndef TMR_MATRIX_FILTER_MODEL_H
#define TMR_MATRIX_FILTER_MODEL_H

#include "TACSElementModel.h"

class TMRQuadMatrixModel : public TACSElementModel {
 public:
  TMRQuadMatrixModel();
  int getSpatialDim();
  int getVarsPerNode();
  void evalWeakIntegrand( int elemIndex, const double time,
                          int n, const double pt[],
                          const TacsScalar X[], const TacsScalar Xd[],
                          const TacsScalar Ut[], const TacsScalar Ux[],
                          TacsScalar DUt[], TacsScalar DUx[] );
  void evalWeakJacobian( int elemIndex, const double time,
                         int n, const double pt[],
                         const TacsScalar X[], const TacsScalar Xd[],
                         const TacsScalar Ut[], const TacsScalar Ux[],
                         TacsScalar DUt[], TacsScalar DUx[],
                         int *Jac_nnz, const int *Jac_pairs[],
                         TacsScalar Jac[] );
 private:
  static int elem_Jac_pairs[2];
};

class TMRHexaMatrixModel : public TACSElementModel {
 public:
  TMRHexaMatrixModel();
  int getSpatialDim();
  int getVarsPerNode();
  void evalWeakIntegrand( int elemIndex, const double time,
                          int n, const double pt[],
                          const TacsScalar X[], const TacsScalar Xd[],
                          const TacsScalar Ut[], const TacsScalar Ux[],
                          TacsScalar DUt[], TacsScalar DUx[] );
  void evalWeakJacobian( int elemIndex, const double time,
                         int n, const double pt[],
                         const TacsScalar X[], const TacsScalar Xd[],
                         const TacsScalar Ut[], const TacsScalar Ux[],
                         TacsScalar DUt[], TacsScalar DUx[],
                         int *Jac_nnz, const int *Jac_pairs[],
                         TacsScalar Jac[] );
 private:
  static int elem_Jac_pairs[2];
};

#endif // TMR_MATRIX_FILTER_MODEL_H
