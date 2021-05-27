/**
 * @file   LinAlg.cpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Mon Oct  3 15:51:58 2011
 * 
 * @brief  
 * 
 * 
 */

#include <cmath>

#include <iostream>
#include <algorithm>

#include "Config.hpp"
#include "Vec.hpp"
#include "LinOp.hpp"
#include "LinAlg.hpp"

template <>
void svd<double>( Mat_d A, Vec_d sgm, Mat_d U_VT )
{
  long    m = A.m();
  long    n = A.n();

  char    jobz = 'O';
  double  lwork_q;
  long    lwork = -1;
  long   *iwork;
  double *work;
  long    info;
  
  if( m >= n )
  {
    iwork = (long *) malloc( 8*n * sizeof(long) );
    assert( iwork != NULL );
    
    /* query */
    dgesdd( &jobz, &m, &n, A.data(), &m, sgm.data(), NULL, &m, U_VT.data(), &n,
            &lwork_q, &lwork, iwork, &info );

    lwork = (long) lwork_q;
    work  = (double *) malloc( lwork * sizeof(double) );
    assert( work != NULL ); 
    
    /* svd */
    dgesdd( &jobz, &m, &n, A.data(), &m, sgm.data(), NULL, &m, U_VT.data(), &n,
            work,     &lwork, iwork, &info );

    free(work);
    free(iwork);
  }
  else
  {
    iwork = (long *) malloc( 8*m * sizeof(long) );
    assert( iwork != NULL );

    /* query */
    dgesdd( &jobz, &m, &n, A.data(), &m, sgm.data(), U_VT.data(), &m, NULL, &m, 
            &lwork_q, &lwork, iwork, &info );

    lwork = (long) lwork_q;
    work  = (double *) malloc( lwork * sizeof(double) );
    assert( work != NULL );

    /* svd */
    dgesdd( &jobz, &m, &n, A.data(), &m, sgm.data(), U_VT.data(), &m, NULL, &m, 
            work,     &lwork, iwork, &info );

    free(work);
    free(iwork);
  }
}

template <>
void chol<double>( Mat_d A )
{
  assert( A.m() == A.n() );

  char uplo = 'L';
  lapack_int info;

  Mat_d const& Ac = A;
  
  dpotrf( &uplo, &Ac.m(), Ac.data(), &Ac.ld(), &info );

  assert( info == 0 );
}

template <>
void chol_sol<double>( const Mat_d A, Mat_d B )
{
  assert( A.m() == A.n() );
  assert( A.m() == B.m() );

  char uplo = 'L';
  lapack_int info;

  Mat_d const& Ac = A;
  Mat_d const& Bc = B;
  
  dpotrs( &uplo, &Ac.m(), &Bc.n(), Ac.data(), &Ac.ld(), Bc.data(), &Bc.ld(), &info );

  assert( info == 0 );
}

template <>
void chol_sol<double>( const Mat_d A, Vec_d b )
{
  assert( A.m() == A.n() );
  assert( A.m() == b.n() );

  char uplo = 'L';
  lapack_int info;
  lapack_int one = 1;

  Mat_d const& Ac = A;
  
  if( b.inc() == 1 )
  {
    dpotrs( &uplo, &Ac.m(), &one, Ac.data(), &Ac.ld(), b.data(), &b.n(), &info );
    assert( info == 0 );
  }
  else
  {
    Vec_d x = b.copy();
    dpotrs( &uplo, &Ac.m(), &one, Ac.data(), &Ac.ld(), b.data(), &b.n(), &info );
    assert( info == 0 );
    copy( x, b );
  }
}

template <>
void lstsq<double>( Mat_d A, Mat_d& BX, Vec<double> s,
                    const double rcond, long& rank )
{
  assert( A.m() == BX.m() );
  assert( A.n() <= BX.ld() );
  
  long m      = A.m();
  long n      = A.n();
  long nrhs   = BX.n();
  
  long min_mn = std::min(m,n);

  assert( s.n() == min_mn );
  if( s.inc() != 1 )
    s = s.copy();
  
  double lwork_q;
  long   lwork = -1;          /* query */
  long   info;

  /* dgelsd doesn't return best liwork in iwork(1) */
  long nlvl    = (long) ceil( log2( min_mn/2.0 ) + 1 );
  nlvl         = std::max( nlvl, 0L );
  long liwork  = 3*min_mn*nlvl + 11*min_mn;

  long  *iwork = new long [liwork];

  // query
  dgelsd( &m, &n, &nrhs,
          A.data(), &A.ld(),
          BX.data(), &BX.ld(),
          s.data(), &rcond, &rank,
          &lwork_q, &lwork, iwork,
          &info );

  lwork = (long) lwork_q;

  double *work = new double [lwork];
  
  dgelsd( &m, &n, &nrhs,
          A.data(), &A.ld(),
          BX.data(), &BX.ld(),
          s.data(), &rcond, &rank,
          work, &lwork, iwork,
          &info );

  delete [] iwork;
  delete [] work;
  
  assert( info == 0 );

  BX = Mat_d(n,nrhs,BX);
}
