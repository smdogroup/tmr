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

#ifndef TMR_LAPACK_H
#define TMR_LAPACK_H

/*
  LAPACK definitions
*/
#define TmrLAPACKdgbsv dgbsv_
#define TmrLAPACKdgels dgels_
#define TmrLAPACKdgetrf dgetrf_
#define TmrLAPACKdgetrs dgetrs_
#define TmrLAPACKsyevd dsyevd_

extern "C" {
extern void TmrLAPACKdgbsv(int *n, int *ku, int *kv, int *nrhs, double *ab,
                           int *ldab, int *ipiv, double *b, int *ldb,
                           int *info);

extern void TmrLAPACKdgels(char *trans, int *m, int *n, int *nrhs, double *A,
                           int *lda, double *rhs, int *ldrhs, double *work,
                           int *lwork, int *info);

extern void TmrLAPACKdgetrf(int *m, int *n, double *a, int *lda, int *ipiv,
                            int *info);

extern void TmrLAPACKdgetrs(const char *c, int *n, int *nrhs, double *a,
                            int *lda, int *ipiv, double *b, int *ldb,
                            int *info);

extern void TmrLAPACKsyevd(const char *jobz, const char *uplo, int *N,
                           double *A, int *lda, double *w, double *work,
                           int *lwork, int *iwork, int *liwork, int *info);
}

#endif
