/**
 * @file   SpMat.cpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Sat Oct  8 15:50:37 2011
 * 
 * @brief  
 * 
 * 
 */

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <pthread.h>

#include "Config.hpp"
#include "Utils.hpp"
#include "SpMat.hpp"
#include "fblas.h"

class CSC_MV_Task 
{
 public:
  char       trans;
  double     alpha;
  CSC_Mat_d  A;
  Vec_d      x;
  double     beta;
  Vec_d      y;
};

void csc_mv( const char trans, const long m, const long n, const double alpha,
             long * const ir, long * const jc, double * const a,
             double * const x, const long incx,
             const double beta, double * const y, const long incy )
{
  long i, j, idx;

  if( trans == 'n' )
  {
    /* y = beta*y */
    if( beta == 0.0 )
      for( i=0; i<m; ++i )
        y[i*incy] = 0.0;
    else
      dscal( &m, &beta, y, &incy );
    /* y += alpha*A*x */
    for( j=0; j<n; ++j )
      for( idx = jc[j]; idx < jc[j+1]; ++idx )
        y[ir[idx]*incy] += alpha*a[idx]*x[j*incx];
  }
  else if( trans == 't' )
  {
    /* y = beta*y */
    if( beta == 0.0 )
      for( i=0; i<n; ++i )
        y[i*incy] = 0.0;
    else
      dscal( &n, &beta, y, &incy );
    /* y += alpha*A'*x */
    for( j=0; j<n; ++j )
      for( idx = jc[j]; idx < jc[j+1]; ++idx )
        y[j*incy] += alpha*a[idx]*x[ir[idx]*incx];
  }
  else
  {
    throw std::runtime_error( "trans must be either 'n' or 't'." );
  }
}

void csc_mv( const char trans, const double alpha, const CSC_Mat_d A,
             const Vec_d x, const double beta, Vec_d y )
{
  csc_mv( trans, A.m(), A.n(), alpha, A.ir(), A.jc(), A.data(),
          x.data(), x.inc(), beta, y.data(), y.inc() );  
}

void *csc_mv_thread( void *task )
{
  CSC_MV_Task *t = (CSC_MV_Task *) task;
  csc_mv( t->trans, t->alpha, t->A, t->x, t->beta, t->y );
  pthread_exit(NULL);
}

void csc_mv_mt( const char trans, const double alpha, const CSC_Mat_d A,
                const Vec<double> x, const double beta, Vec<double> y )
{
  long n_threads = get_num_cores();

  pthread_t *threads = new pthread_t [n_threads];
  CSC_MV_Task *tasks = new CSC_MV_Task [n_threads];

  long n = A.n();
  
  if( trans == 'n' )
  {
    Mat_d Y( A.m(), n_threads );
    
    for( long i=0; i<n_threads; ++i )
    {
      long begin = ceil(1.0*n* i   /n_threads);
      long end   = ceil(1.0*n*(i+1)/n_threads);

      CSC_MV_Task *t = &tasks[i];
      
      t->trans = trans;
      t->alpha = alpha;
      t->A     = A.cols(begin,end);
      t->x     = x.subvec(begin,end);
      t->beta  = 0.0;
      t->y     = Y.col(i);

      pthread_create( &threads[i], NULL, csc_mv_thread, (void *) t );
    }

    for( int i=0; i<n_threads; ++i )
      pthread_join( threads[i], NULL );
    
    Vec_d ones(n_threads);
    ones = 1.0;
    mv( 'n', 1.0, Y, ones, beta, y );
  }
  else if( trans == 't' )
  {
    for( long i=0; i<n_threads; ++i )
    {
      long begin = ceil(1.0*n* i   /n_threads);
      long end   = ceil(1.0*n*(i+1)/n_threads);

      CSC_MV_Task *t = &tasks[i];

      t->trans = trans;
      t->alpha = alpha;
      t->A     = A.cols(begin,end);
      t->x     = x;
      t->beta  = beta;
      t->y     = y.subvec(begin,end);

      pthread_create( &threads[i], NULL, csc_mv_thread, (void *) t );
    }

    for( long i=0; i<n_threads; ++i )
      pthread_join( threads[i], NULL );
  }
  else
  {
    throw std::runtime_error( "trans must be either 'n' or 't'." );
  }

  delete [] threads;
  delete [] tasks;
}

template <>
void mv<double>( const char trans, const double alpha, const CSC_Mat_d A,
                 const Vec_d x, const double beta, Vec_d y )
{
  csc_mv_mt( trans, alpha, A, x, beta, y );
}

void csc_mm( const char trans_a, const char trans_b, const double alpha,
             const CSC_Mat_d A, const Mat_d B, const double beta, Mat_d C )
{
  if( trans_b == 'n' )
  {
    for( long j=0; j<B.n(); ++j )
      csc_mv( trans_a, alpha, A, B.col(j), beta, C.col(j) );
  }
  else if( trans_b == 't' )
  {
    for( long j=0; j<B.m(); ++j )
      csc_mv( trans_a, alpha, A, B.row(j), beta, C.col(j) );
  }
  else
  {
    throw std::runtime_error( "trans must be either 'n' or 't'." );
  }
}

class CSC_MM_Task
{
 public:
  char      trans_a;
  char      trans_b;
  double    alpha;
  CSC_Mat_d A;
  Mat_d     B;
  double    beta;
  Mat_d     C;
};

void *csc_mm_thread( void *task )
{
  CSC_MM_Task *t = (CSC_MM_Task *) task;
  csc_mm( t->trans_a, t->trans_b, t->alpha, t->A, t->B, t->beta, t->C );
  pthread_exit(NULL);
}

void csc_mm_mt( const char trans_a, const char trans_b, const double alpha,
                const CSC_Mat_d A, const Mat_d B, const double beta, Mat_d C )
{
  long n_threads = get_num_cores();

  pthread_t *threads = new pthread_t [n_threads];
  CSC_MM_Task *tasks = new CSC_MM_Task [n_threads];

  if( trans_a == 'n' )
  {
    long n = A.n();
    long n_elem = C.m()*C.n();
    Mat_d W( n_elem, n_threads );

    for( long i=0; i<n_threads; ++i )
    {
      long begin = ceil(1.0*n* i   /n_threads);
      long end   = ceil(1.0*n*(i+1)/n_threads);
      
      CSC_MM_Task *t = &tasks[i];

      t->trans_a = trans_a;
      t->trans_b = trans_b;
      t->alpha   = alpha;
      t->A       = A.cols(begin,end);
      if( trans_b == 'n' )
        t->B     = B.submat(begin,end,0,B.n());
      else if( trans_b == 't' )
        t->B     = B.submat(0,B.m(),begin,end);
      else
        throw std::runtime_error( "trans must be either 'n' or 't'." );
      t->beta    = 0.0;
      t->C       = Mat<double>( C.m(), C.n(), W.data()+i*n_elem );
      
      pthread_create( &threads[i], NULL, csc_mm_thread, (void *) t ); 
   }

    for( long i=0; i<n_threads; ++i )
      pthread_join( threads[i], NULL );
    
    Vec_d ones(n_threads);
    ones = 1.0;

    for( long j=0; j<C.n(); ++j )
    {
      Mat_d Wj = W.submat(j*C.m(),(j+1)*C.m(),0,n_threads);
      mv( 'n', 1.0, Wj, ones, beta, C.col(j) );
    }
  }
  else if( trans_a == 't' )
  {
    long n = A.n();
    for( long i=0; i<n_threads; ++i )
    {
      long begin = ceil(1.0*n* i   /n_threads);
      long end   = ceil(1.0*n*(i+1)/n_threads);
      
      CSC_MM_Task *t = &tasks[i];
      
      t->trans_a = trans_a;
      t->trans_b = trans_b;
      t->alpha   = alpha;
      t->A       = A.cols(begin,end);
      t->B       = B;
      t->beta    = beta;
      t->C       = C.submat(begin,end,0,C.n());

      pthread_create( &threads[i], NULL, csc_mm_thread, (void *) t );
    }

    for( long i=0; i<n_threads; ++i )
      pthread_join( threads[i], NULL );
  }

  delete [] threads;
  delete [] tasks;
}

template <>
void mm<double>( const char trans_a, const char trans_b, const double alpha,
                 const CSC_Mat_d A, const Mat_d B, const double beta, Mat_d C )
{
  csc_mm_mt( trans_a, trans_b, alpha, A, B, beta, C );
}

void csc_to_csr( long m, long n,
                 long *csc_ia, long *csc_ja, double *csc_a,
                 long *csr_ia, long *csr_ja, double *csr_a )
{
  long *row_nnz = (long *) malloc( m * sizeof(long) );
  assert( row_nnz != NULL );
  
  /* count number of non-zeros in each row */
  memset( row_nnz, 0, m * sizeof(long) );
  for( long k=csc_ja[0]; k<csc_ja[n]; ++k )
    row_nnz[csc_ia[k]]++;

  /* build csr_ia */
  csr_ia[0] = 0;
  for( long i=1; i<m+1; ++i )
    csr_ia[i] = csr_ia[i-1] + row_nnz[i-1];

  /* build csr_ja & csr_a */
  memset( row_nnz, 0, m * sizeof(long) );
  for( long j=0; j<n; ++j )
    for( long k=csc_ja[j]; k<csc_ja[j+1]; ++k )
    {
      long i      = csc_ia[k];
      long idx    = csr_ia[i] + row_nnz[i];
      csr_ja[idx] = j;
      csr_a[idx]  = csc_a[k];
      row_nnz[i]++;
    }
  free( row_nnz );
}

template <>
void csc_to_csr<double>( const CSC_Mat_d C, CSR_Mat_d R )
{
  assert( C.m() == R.m() && C.n() == R.n() );
  assert( C.nnz() <= R.nzmax() );
  csc_to_csr( C.m(), C.n(),
              C.ir(), C.jc(), C.data(),
              R.ir(), R.jc(), R.data() );
}

template <>
void mv<double>( const char trans, const double alpha, const CSR_Mat_d A,
                 const Vec_d x, const double beta, Vec_d y )
{
  CSC_Mat_d At( A.n(), A.m(), A.nzmax(), A.jc(), A.ir(), A.data(), A.jc_blk(), A.ir_blk(), A.data_blk() );
  char trans_t = trans == 'n' ? 't' : 'n';
  mv( trans_t, alpha, At, x, beta, y );
}


template <>
void mm<double>( const char trans_a, const char trans_b, const double alpha,
                 const CSR_Mat_d A, const Mat_d B, const double beta, Mat_d C )
{
  CSC_Mat_d At( A.n(), A.m(), A.nzmax(), A.jc(), A.ir(), A.data(), A.jc_blk(), A.ir_blk(), A.data_blk() );
  char trans_t = trans_a == 'n' ? 't' : 'n';
  mm( trans_t, trans_b, alpha, At, B, beta, C );
}

