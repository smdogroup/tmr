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

#ifndef TMR_HELMHOLTZ_PARTITION_UNITY_FILTER_H
#define TMR_HELMHOLTZ_PARTITION_UNITY_FILTER_H

#include "TMRConformFilter.h"

/*
  Create a partition of unity filter
*/
class TMRHelmholtzPUFilter : public TMRConformFilter {
 public:
  TMRHelmholtzPUFilter(int _N, int _nlevels, TACSAssembler *_assembler[],
                       TMROctForest *_filter[]);
  TMRHelmholtzPUFilter(int _N, int _nlevels, TACSAssembler *_assembler[],
                       TMRQuadForest *_filter[]);
  ~TMRHelmholtzPUFilter();

  // Compute the stencil at an interior node
  virtual int getInteriorStencil(int diagonal_index, int npts,
                                 const TacsScalar Xpts[], double alpha[]) = 0;

  // Compute the stencil at a boundary node with the normal n
  virtual int getBoundaryStencil(int diagonal_index, const TacsScalar n[],
                                 int npts, const TacsScalar Xpts[],
                                 double alpha[]) = 0;

  // Set the design variable values (including all local values)
  void setDesignVars(TACSBVec *x);

  // Set values/add values to the vector
  void addValues(TACSBVec *vec);

  void initialize();

 private:
  // Apply the filter to get the density values
  void applyFilter(TACSBVec *in, TACSBVec *out);

  // Apply the transpose of the filter for sensitivities
  void applyTranspose(TACSBVec *in, TACSBVec *out);

  // The non-negative matrix M
  TACSMat *B;

  // The number of terms to include in the approximate inverse
  int N;

  // Store the inverse of the diagonal matrices
  TACSBVec *Dinv, *Tinv;

  // Temporary vectors required for the matrix computation
  TACSBVec *t1, *t2;

  // Another set of temporary vectors
  TACSBVec *y1, *y2;

  // Temporary design variable vector
  TACSBVec *temp;

  // Compute the Kronecker product
  void kronecker(TACSBVec *c, TACSBVec *x, TACSBVec *y = NULL);
};

/*
  Create a partition of unity filter with callbacks
*/
class TMRCallbackHelmholtzPUFilter : public TMRHelmholtzPUFilter {
 public:
  TMRCallbackHelmholtzPUFilter(int _N, int _nlevels,
                               TACSAssembler *_assembler[],
                               TMROctForest *_filter[])
      : TMRHelmholtzPUFilter(_N, _nlevels, _assembler, _filter) {
    self = NULL;
    getinteriorstencil = NULL;
    getboundarystencil = NULL;
  }
  TMRCallbackHelmholtzPUFilter(int _N, int _nlevels,
                               TACSAssembler *_assembler[],
                               TMRQuadForest *_filter[])
      : TMRHelmholtzPUFilter(_N, _nlevels, _assembler, _filter) {
    self = NULL;
    getinteriorstencil = NULL;
    getboundarystencil = NULL;
  }
  ~TMRCallbackHelmholtzPUFilter() {}

  // Set the pointers/function pointers
  void setSelfPointer(void *_self) { self = _self; }
  void setGetInteriorStencil(int (*func)(void *, int, int, const TacsScalar *,
                                         double *)) {
    getinteriorstencil = func;
  }
  void setGetBoundaryStencil(int (*func)(void *, int, const TacsScalar *, int,
                                         const TacsScalar *, double *)) {
    getboundarystencil = func;
  }

  // Compute the stencil at an interior node
  int getInteriorStencil(int diagonal_index, int npts, const TacsScalar Xpts[],
                         double alpha[]) {
    if (self && getinteriorstencil) {
      return getinteriorstencil(self, diagonal_index, npts, Xpts, alpha);
    }
    return 1;
  }

  // Compute the stencil at a boundary node with the normal n
  int getBoundaryStencil(int diagonal_index, const TacsScalar n[], int npts,
                         const TacsScalar Xpts[], double alpha[]) {
    if (self && getboundarystencil) {
      return getboundarystencil(self, diagonal_index, n, npts, Xpts, alpha);
    }
    return 1;
  }

 private:
  void *self;
  int (*getinteriorstencil)(void *, int, int, const TacsScalar *, double *);
  int (*getboundarystencil)(void *, int, const TacsScalar *, int,
                            const TacsScalar *, double *);
};

#endif  // TMR_HELMHOLTZ_PARTITION_UNITY_FILTER_H
