#ifndef _FLAPACK_H
#define _FLAPACK_H

#include <stddef.h>
#include "config_sk.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef INT lapack_int;
  
#define dgesdd FWRAPPER(dgesdd,DGESDD)
  extern void dgesdd( const char* jobz, const lapack_int *m, const lapack_int *n,
                      const double* a, const lapack_int *lda,
                      const double* s,
                      const double* u, const lapack_int *ldu,
                      const double* vt, const lapack_int *ldvt,
                      const double* work, const lapack_int *lwork, const lapack_int *iwork,
                      const lapack_int *info );
  
  
#define dgelsd FWRAPPER(dgelsd,DGELSD)
  extern void dgelsd( const lapack_int *m, const lapack_int *n, const lapack_int *nrhs,
                      const double *a, const lapack_int *lda,
                      const double *b, const lapack_int *ldb,
                      const double *s, const double *rcond, const lapack_int *rank,
                      const double *work, const lapack_int *lwork, const lapack_int *iwork,
                      const lapack_int *info );
  
#define dpotrf FWRAPPER(dpotrf,DPOTRF)
  extern void dpotrf( const char *uplo, const lapack_int *n,
                      const double *a, const lapack_int *lda, const lapack_int *info );

#define dpotrs FWRAPPER(dpotrs,DPOTRS)
  extern void dpotrs( const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                      const double *a, const lapack_int *lda,
                      const double *b, const lapack_int *ldb, const lapack_int *info );
  
#define dgeqp3 FWRAPPER(dgeqp3,DGEQP3)
  extern void dgeqp3( const lapack_int *m, const lapack_int *n,
                      const double *a, const lapack_int *lda, const lapack_int *E, const double*tau,
                      const double *work, const lapack_int *lwork, const lapack_int *info);

#define dgeqrf FWRAPPER(dgeqrf,DGEQRF)
extern void dgeqrf( const lapack_int *m, const lapack_int *n,
                    const double *a, const lapack_int *lda, const double*tau,
                    const double *work, const lapack_int *lwork, const lapack_int *info);


#define dtrcon FWRAPPER(dtrcon,DTRCON)
  extern void dtrcon( const char* norm, const char* uplo, const char* diag, 
                      const lapack_int *n, const double* a, const lapack_int * ld,
                      const double* rcond, const double* work, const lapack_int *lwork, 
                      const lapack_int *info);

#define dlapmt FWRAPPER(dlapmt,DLAPMT)
  extern void dlapmt( lapack_int *forward, const long *m, const long *n, const double* x, 
                      const long *lda, long* E);
#define dormqr FWRAPPER(dormqr,DORMQR)
  extern void dormqr( const char *side, const char *trans, const lapack_int *m, 
                      const lapack_int *n, const lapack_int *k, 
                      const double *a, const lapack_int *lda, const double *t, 
                      const double *c, const lapack_int *ldc, 
                      const double* work, const lapack_int *lwork, const lapack_int *info);
#define dgelsy FWRAPPER(dgelsy,DGELSY)
  extern void dgelsy( const lapack_int *m, const lapack_int *n, const lapack_int *nrhs, 
                      const double* a, const lapack_int *lda, const double* b,
                      const lapack_int *ldb, const double* jpvt, const double* rcond,
                      const lapack_int *rank, const double* work, const lapack_int *lwork,
                      const lapack_int *info);


#define dlacpy FWRAPPER(dlacpy,DLACPY)
extern void dlacpy(
    char const* uplo,
    lapack_int const* m, lapack_int const* n,
    double const* A, lapack_int const* lda,
    double* B, lapack_int const* ldb );

#define dorgqr FWRAPPER(dorgqr,DORGQR)
extern void dorgqr(
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double* A, lapack_int const* lda,
    double const* tau,
    double* work, lapack_int const* lwork,
    lapack_int* info );

#define dlarf FWRAPPER(dlarf,DLARF)
extern void dlarf(
    char const* side,
    lapack_int const* m, lapack_int const* n,
    double const* V, lapack_int const* incv,
    double const* tau,
    double* C, lapack_int const* ldc,
    double* work );

#define dlarfb FWRAPPER(dlarfb,DLARFB)
extern void dlarfb(
    char const* side, char const* trans, char const* direct, char const* storev,
    lapack_int const* m, lapack_int const* n, lapack_int const* k,
    double const* V, lapack_int const* ldv,
    double const* T, lapack_int const* ldt,
    double* C, lapack_int const* ldc,
    double* work, lapack_int const* ldwork );

#define dlarfg FWRAPPER(dlarfg,DLARFG)
extern void dlarfg(
    lapack_int const* n,
    double* alpha,
    double* X, lapack_int const* incx,
    double* tau );

#define dlarft FWRAPPER(dlarft,DLARFT)
extern void dlarft(
    char const* direct, char const* storev,
    lapack_int const* n, lapack_int const* k,
    double const* V, lapack_int const* ldv,
    double const* tau,
    double* T, lapack_int const* ldt );



#ifdef __cplusplus
}
#endif

#endif  /* _FLAPACK_H */
