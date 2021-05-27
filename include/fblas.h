#ifndef _FBLAS_H
#define _FBLAS_H

#include <stddef.h>
#include "config_sk.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  typedef INT blas_int;
  
#define dcopy FWRAPPER(dcopy,DCOPY)
  extern void dcopy( const blas_int *n,
                     const double *dx, const blas_int *incx,
                     const double *dy, const blas_int *incy );

#define dgemv FWRAPPER(dgemv,DGEMV)
  extern void dgemv( const char *trans,
                     const blas_int *m, const blas_int *n,
                     const double *alpha,
                     const double *a, const blas_int *lda,
                     const double *x, const blas_int *incx,
                     const double *beta,
                     const double *y, const blas_int *incy );

#define dgemm FWRAPPER(dgemm,DGEMM)
  extern void dgemm( const char *trans_a, const char *trans_b,
                     const blas_int *m, const blas_int *n, const blas_int *k,
                     const double *alpha,
                     const double *a, const blas_int *lda,
                     const double *b, const blas_int *ldb,
                     const double *beta,
                     const double *c, const blas_int *ldc );

#define dnrm2 FWRAPPER(dnrm2,DNRM2)
  extern double dnrm2( const blas_int *n,
                             const double *x, const blas_int *incx );
  
#define dscal FWRAPPER(dscal,DSCAL)
  extern void dscal( const blas_int *n,
                     const double *da,
                     const double *dx, const blas_int *incx );

#define daxpy FWRAPPER(daxpy,DAXPY)
  extern void daxpy( const blas_int *n,
                     const double *da,
                     const double *dx, const blas_int *incx,
                     const double *dy, const blas_int *incy );

#define dtrsv FWRAPPER(dtrsv,DTRSV)
  extern void dtrsv( char* uplo, char* trans, char* diag,
                     const blas_int *n, const double* A, const blas_int* lda,
                     const double* x, const blas_int* incx );

#define dswap FWRAPPER(dswap,DSWAP)
  extern void dswap( blas_int *, double *, blas_int *, double *, blas_int *);

#define dtrmm FWRAPPER(dtrmm,DTRMM)
  extern void dtrmm(char *, char *, char *, char *, blas_int *, blas_int *, const double *, const double *, blas_int *, double *, blas_int *);

#define xerbla FWRAPPER(xerbla,XERBLA)
  extern void xerbla(char *, void *);

#ifdef __cplusplus
}
#endif
  
#endif  /* _FBLAS_H */
