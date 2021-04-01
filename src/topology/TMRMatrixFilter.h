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

#include "TMRConformFilter.h"

/*
  The following class creates and stores approximate M-filters for
  topology optimization.  These filters can be constructed via
  different strategies that ensure that the matrices satisfy the
  positivity and sum-to-unity properties.

  Mass-matrix approach:

  Given a matrix M with non-negative entries, and a scalar r >= 0, a
  filter matrix F can be formed as follows. First construct a diagonal
  matrix D such that

  D_{i} = sum_{j} M_{ij}

  The filter matrix F, can be formed implicitly as

  F = (D + r^2*D^{-2/d}*(D - M))^{-1}*D

  or

  F = (1 + r^2*D^{-2/d-1}*(D - M))^{-1}

  Note that the filter matrix satisfies the property that F*e = e, where
  e is a vector of all unit entries.

  Introducing the definitions:

  Ainv = (1 + r^2*D^{-2/d})^{-1}

  and

  B = r^2*D^{-2/d-1}*Ainv

  Then the approximate action of the filter matrix can be written using a
  Neumann series as

  F ~ T^{T} * \sum_{n=0}^{N} (B*M)^{n} A^{-1}

  Where T is defined as

  T = \sum_{n=0}^{N} (B*M)^{n} A^{-1} e
*/

class TMRMatrixFilter : public TMRConformFilter {
 public:
  TMRMatrixFilter( double _r, int _N, int _nlevels,
                   TACSAssembler *_tacs[],
                   TMROctForest *_filter[] );
  TMRMatrixFilter( double _r, int _N, int _nlevels,
                   TACSAssembler *_tacs[],
                   TMRQuadForest *_filter[] );
  ~TMRMatrixFilter();

  // Set the design variable values (including all local values)
  void setDesignVars( TACSBVec *x );

  // Set values/add values to the vector
  void addValues( TACSBVec *vec );

  // Apply the filter to get the density values
  void applyFilter( TACSBVec *in, TACSBVec *out );

  // Apply the transpose of the filter for sensitivities
  void applyTranspose( TACSBVec *in, TACSBVec *out );

 private:
  void initialize_matrix( double _r, int _N,
                          TMROctForest *oct_filter,
                          TMRQuadForest *quad_filter );

  // The non-negative matrix M
  TACSMat *M;

  // The number of terms to include in the approximate inverse
  int N;

  // The scalar for the approximate Helmholtz
  double r;

  // Store the inverse of the diagonal matrices
  TACSBVec *Ainv, *B, *Tinv;

  // Temporary vectors required for the matrix computation
  TACSBVec *t1, *t2;

  // Another set of temporary vectors
  TACSBVec *y1, *y2;

  // Temporary design variable vector
  TACSBVec *temp;

  // Compute the Kronecker product
  void kronecker( TACSBVec *c, TACSBVec *x, TACSBVec *y=NULL );
};

#endif // TMR_MATRIX_FILTER_H
