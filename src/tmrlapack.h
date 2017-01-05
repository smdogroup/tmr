#ifndef TMR_LAPACK_H
#define TMR_LAPACK_H

/*
  LAPACK definitions
*/
#define LAPACKdgbsv dgbsv_
#define LAPACKdgels dgels_

extern "C" {
  extern void LAPACKdgbsv( int *n, int *ku, int *kv, int *nrhs, 
                           double *ab, int *ldab, int *ipiv, double *b, int *ldb, 
                           int *info );

  extern void LAPACKdgels( char *trans, int *m, int *n, int *nrhs,
                           double *A, int *lda, double *rhs, int *ldrhs, 
                           double *work, int *lwork, int *info );
}

#endif
