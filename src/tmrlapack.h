#ifndef TMR_LAPACK_H
#define TMR_LAPACK_H

/*
  LAPACK definitions
*/
#define TmrLAPACKdgbsv dgbsv_
#define TmrLAPACKdgels dgels_
#define TmrLAPACKdgetrf dgetrf_
#define TmrLAPACKdgetrs dgetrs_

extern "C" {
  extern void TmrLAPACKdgbsv( int *n, int *ku, int *kv, int *nrhs, 
                              double *ab, int *ldab, int *ipiv, double *b, int *ldb, 
                              int *info );

  extern void TmrLAPACKdgels( char *trans, int *m, int *n, int *nrhs,
                              double *A, int *lda, double *rhs, int *ldrhs, 
                              double *work, int *lwork, int *info );

  extern void TmrLAPACKdgetrf( int *m, int *n, 
                               double *a, int *lda, int *ipiv, int *info );

  extern void TmrLAPACKdgetrs( const char *c, int *n, int *nrhs, 
                               double *a, int *lda, int *ipiv, 
                               double *b, int *ldb, int *info );
}

#endif
