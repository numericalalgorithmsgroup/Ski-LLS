/**
 * @file   Mat.cpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Tue Oct  4 23:46:26 2011
 * 
 * @brief  
 * 
 * 
 */

#include <cassert>

#include "fblas.h"
#include "Mat.hpp"

template <>
void mv<double>( const char trans, const double alpha, const Mat_d A,
                 const Vec_d x, const double beta, Vec_d y )
{
  if( trans == 'n' )
    assert( A.m() == y.n() && A.n() == x.n() );
  else
    assert( A.m() == x.n() && A.n() == y.n() );

  dgemv( &trans, &A.m(), &A.n(), &alpha, A.data(), &A.ld(),
         x.data(), &x.inc(), &beta, y.data(), &y.inc() );
}

template <>
void mm<double>( const char trans_a, const char trans_b, const double alpha,
                 const Mat_d A, const Mat_d B,
                 const double beta, Mat_d C )
{
  if( trans_a == 'n' )
    assert( A.m() == C.m() );
  else
    assert( A.n() == C.m() );
  
  Mat_d::idx_t k = trans_a == 'n' ? A.n() : A.m();

  if( trans_b == 'n' )
    assert( B.m() == k && B.n() == C.n() );
  else
    assert( B.n() == k && B.m() == C.n() );

  dgemm( &trans_a, &trans_b, &C.m(), &C.n(), &k,
         &alpha, A.data(), &A.ld(), B.data(), &B.ld(),
         &beta, C.data(), &C.ld() );
}

