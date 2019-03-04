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

#ifndef TMR_MATRIX_FILTER_H
#define TMR_MATRIX_FILTER_H

#include "TMRQuadForest.h"
#include "TMROctForest.h"

/*
  The following class creates and stores approximate M-filters for
  topology optimization.  These filters can be constructed via
  different strategies that ensure that the matrices satisfy the
  positivity and sum-to-unity properties.

  Mass-matrix approach:

  Given a matrix M with non-negative entries, and a scalar s > 1, a
  filter matrix F can be formed as follows. First construct a diagonal
  matrix D such that

  D_{i} = 1/(s-1)*sum_{j} M_{ij}

  Then the matrix A = (sI - D^{-1}M) is an M matrix, and F = A^{-1}
  satisfies the properties of a partition-of-unity filter.

  Using a Neumann series, the inverse of A can be expressed as

  A^{-1} ~ \sum_{n=0}^{infty} 1/s**(n+1)*(D^{-1}*M)**n

  However, this series has to be truncated. Since all entries of D and
  M are positive, this filter maintains positivity. However, the
  property Fe = e will be lost. This can be restored by introducing a
  truncation value N, and a scaling diagonal matrix

  T = Diag{\sum_{n=0}^{N} 1/s**(n+1)(D^{-1}*M)**n e}

  So the final filter can be written as follows:

  F = T^{-1}*[ \sum_{n=0}^{infty} 1/s**(n+1)*(D^{-1}*M)**n ]
*/

class TMRMatrixFilter : public TMREntity {
 public:
  TMRMatrixFilter( TACSMat *_M, double _s, int _N );
  ~TMRMatrixFilter();

  // Apply the filter to get the density values
  void applyFilter( TACSVec *in, TACSVec *out );
  
  // Apply the transpose of the filter for sensitivities
  void applyTranspose( TACSVec *in, TACSVec *out );

 private:
  // The non-negative matrix M
  TACSMat *M;

  // The number of terms to include in the approximate inverse
  int N;

  // The scalar > 1.0
  double s;

  // Store the inverse of the diagonal matrices
  TACSVec *Dinv, *Tinv;

  // Temporary vectors required for the matrix computation
  TACSVec *t1, *t2;

  // Compute the Kronecker product
  void kronecker( TACSVec *c, TACSVec *x, TACSVec *y=NULL );
};

#endif // TMR_MATRIX_FILTER_H
