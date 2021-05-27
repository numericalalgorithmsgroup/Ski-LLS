/**
 * @file   Vec.hpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Thu Jul 28 15:09:13 2011
 * 
 * @brief  vector class
 * 
 * 
 */

#ifndef _VEC_HPP
#define _VEC_HPP

#include <iostream>

#include "Blk.hpp"

template <typename T>
class Vec;

typedef Vec<double> Vec_d;
typedef Vec<long>   Vec_l;

template <typename T>
T nrm2( const Vec<T> v );

template <typename T>
void copy( const Vec<T> src, Vec<T> dst )
{
  assert( src.n() == dst.n() );
  for( typename Vec<T>::idx_t i=0; i<src.n(); ++i )
    dst(i) = src(i);
}

template <typename T>
void axpy( const T alpha, const Vec<T> x, Vec<T> y );

template <typename T>
void scal( const T alpha, Vec<T> x );

/**
 * A vector class.
 * 
 */
template <typename T>
class Vec
{
 public:

  typedef T         val_t;
  typedef ptrdiff_t idx_t;
  typedef Blk<T>    blk_t;
  typedef Vec<T>    self_t;
  
  /** 
   * create a vector of size sz
   * 
   * @param sz size
   */
  Vec( const idx_t sz = 0 )
      : _blk(sz)
  {
    _n    = sz;
    _inc  = 1;
    _data = _blk.begin();
  }
  
  /** 
   * create a vector from allocated ram
   * 
   * @param sz size
   * @param data pointer to the beginning of the data
   * @param inc increments
   * @param blk (optional)
   */
  Vec( const idx_t sz, val_t * const data, const idx_t inc = 1, const blk_t blk = blk_t() )
  {
    _n    = sz;
    _data = data;
    _inc  = inc;
    _blk  = blk;

    _assert_consistency();
  }
  
  Vec( const idx_t sz, const self_t v, const idx_t offset = 0, const idx_t inc = 1 )
      : _blk(v._blk)
  {
    _n    = sz;
    _data = v._data + offset*v._inc;
    _inc  = v._inc*inc;
    
    _assert_consistency();
  }

  self_t& operator=( const val_t x )
  {
    for( idx_t i=0; i<_n; ++i )
      (*this)(i) = x;
    return *this;
  }

  self_t copy() const
  {
    self_t v(_n);
    ::copy( *this, v );
    return v;
  }

  idx_t   const& n()    const { return _n; }
  val_t * const& data() const { return _data; }
  idx_t   const& inc()  const { return _inc; }

  val_t const& operator()( const idx_t i ) const
  {
    return _data[i*_inc];
  }
  
  val_t& operator()( const idx_t i )
  {
    return _data[i*_inc];
  }

  self_t subvec( idx_t begin, idx_t end ) const
  {
    return self_t( end-begin, *this, begin );
  }
  
 private:

  idx_t  _n;
  val_t *_data;
  idx_t  _inc;
  blk_t  _blk;
  
  void _assert_consistency() const
  {
    if( !_blk.empty() )
      assert( ( _data >= _blk.begin() ) && ( _data+(_n-1)*_inc < _blk.end() ) );
  }
};

template <typename T>
std::ostream& operator<<( std::ostream& out, Vec<T> v )
{
  out << "[";
  for( typename Vec<T>::idx_t i=0; i<v.n(); ++i )
    out << v(i) << ",";
  out << "]";
  return out;
}

#endif  // _VEC_HPP
