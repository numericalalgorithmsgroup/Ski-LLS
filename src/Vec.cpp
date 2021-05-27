/**
 * @file   Vec.cpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Thu Jul 28 15:46:52 2011
 * 
 * @brief  
 * 
 * 
 */

#include <cassert>

#include "Vec.hpp"
#include "fblas.h"

template <>
void scal<double>( const double alpha, Vec_d v )
{
  return dscal( &v.n(), &alpha, v.data(), &v.inc() );
}

template <>
double nrm2<double>( const Vec_d v )
{
  return dnrm2( &v.n(), v.data(), &v.inc() );
}

template <>
void copy<double>( const Vec_d x, Vec_d y )
{
  assert( x.n() == y.n() );
  dcopy( &x.n(), x.data(), &x.inc(), y.data(), &y.inc() );
}

template <>
void axpy<double>( const double alpha, const Vec_d x, Vec_d y )
{
  assert( x.n() == y.n() );
  daxpy( &x.n(), &alpha, x.data(), &x.inc(), y.data(), &y.inc() );
}
