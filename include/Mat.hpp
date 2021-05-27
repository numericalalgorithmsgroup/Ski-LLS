/**
 * @file   Mat.hpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Sun Jul 31 16:28:43 2011
 * 
 * @brief  matrix class
 * 
 * 
 */

#ifndef _MAT_HPP
#define _MAT_HPP

#include <cstddef>

#include <iostream>
#include <algorithm>

#include "Blk.hpp"
#include "Vec.hpp"
#include "LinOp.hpp"

template <typename T>
class LinOp;

template <typename T>
class Mat;

typedef Mat<double> Mat_d;
typedef Mat<long>   Mat_l;

template <typename T>
void mv( const char trans, const T alpha, const Mat<T> A, const Vec<T> x, const T beta, Vec<T> y );

template <typename T>
void mm( const char trans_a, const char trans_b, const T alpha, const Mat<T> A, const Mat<T> B, const T beta, Mat<T> C );

template <typename T>
class Mat
    : public LinOp<T>
{
 public:
  
  typedef T         val_t;
  typedef ptrdiff_t idx_t;
  typedef Mat<T>    self_t;
  typedef Blk<T>    blk_t;
  
  /** 
   * create a matrix
   * 
   * @param sz_m number of rows
   * @param sz_n number of columns
   */
  Mat( const idx_t sz_m = 0, const idx_t sz_n = 0, const idx_t ld = -1 )
      : _blk(sz_m*sz_n)
  {
    _m     = sz_m;
    _n     = sz_n;
    _ld    = std::max( sz_m, ld );
    _data  = _blk.begin();
  }

  /** 
   * create a matrix from allocated ram
   * 
   * @param sz_m number of rows
   * @param sz_n number of columns
   * @param data pointer to the beginning of the data
   * @param ld leading dimension (distance in memory between two consecutive columns, e.g. (1,1) and (1,2)
   * @param blk (optional)
   */
  Mat( const idx_t sz_m, const idx_t sz_n, val_t * const data,
       const idx_t ld = -1, const blk_t blk = blk_t() )
  {
    _m     = sz_m;
    _n     = sz_n;
    _data  = data;
    _ld    = std::max( _m, ld );
    _blk   = blk;
    
    _assert_consistency();
  }

  /** 
   * create a matrix from an existing matrix
   * 
   * @param sz_m number of rows
   * @param sz_n number of columns
   * @param mat matrix
   * @param offset_m row offset
   * @param offset_n column offset
   * @param inc_n inc of column
   */
  Mat( const idx_t sz_m, const idx_t sz_n, const self_t mat,
       const idx_t offset_m = 0, const idx_t offset_n = 0, const idx_t inc_n = 1 )
      : _blk(mat._blk)
  {
    _m    = sz_m;
    _n    = sz_n;
    _ld   = mat._ld*inc_n;
    _data = mat._data + offset_m + offset_n*mat._ld;

    _assert_consistency();
  }

  self_t& operator=( const val_t x )
  {
    for( idx_t i=0; i<_m; ++i )
      for( idx_t j=0; j<_n; ++j )
        (*this)(i,j) = x;
    return *this;
  }

  self_t copy() const
  {
    self_t A( _m, _n );
    for( idx_t i=0; i<_m; ++i )
      for( idx_t j=0; j<_n; ++j )
        A(i,j) = (*this)(i,j);
    return A;
  }

  virtual idx_t const& m() const { return _m; }
  virtual idx_t const& n() const { return _n; }

  Vec<idx_t> shape() const
  {
    Vec<idx_t> dim(2);
    dim(0) = _m;
    dim(1) = _n;
    return dim;
  }

  bool empty() const
  {
    return (_m*_n==0L);
  }
  
  val_t * const& data() const { return _data; }
  idx_t const& ld() const { return _ld; }
  blk_t const& blk() const { return _blk; }
  
  val_t const& operator()( idx_t i, idx_t j ) const
  {
    return _data[i+j*_ld];
  }

  val_t& operator()( const idx_t i, const idx_t j )
  {
    return _data[i+j*_ld];
  }
  
  Vec<T> row( const idx_t i ) const
  {
    return Vec<T>( _n, _data+i, _ld, _blk );
  }

  Vec<T> col( const idx_t j ) const
  {
    return Vec<T>( _m, _data+j*_ld, 1, _blk );
  }

  self_t submat( const idx_t i_begin, const idx_t i_end,
                 const idx_t j_begin, const idx_t j_end ) const 
  {
    return self_t( i_end-i_begin, j_end-j_begin, *this, i_begin, j_begin );
  }
  
  Vec<T> diag() const
  {
    idx_t min_mn = _m < _n ? _m : _n;
    Vec<T> v( min_mn, _data, _ld+1, _blk );
    return v;
  }

  Mat<T> trans() const
  {
    Mat<T> A = *this;
    A._trans = !A._trans;
    return A;
  }
  
  virtual void mv( const char trans, const T alpha, const Vec<T> x, const T beta, Vec<T> y )
  {
    ::mv<double>( trans, alpha, *this, x, beta, y );
  }

  virtual void mm( const char trans, const char trans_b, const T alpha, const Mat<T> B, const T beta, Mat<T> C )
  {
    ::mm<double>( trans, trans_b, alpha, *this, B, beta, C );
  }
      
 private:

  idx_t  _m;
  idx_t  _n;
  idx_t  _ld;
  val_t *_data;
  blk_t  _blk;

  bool   _trans;
  
  void _assert_consistency() const
  {
    if( !_blk.empty() )
      assert( ( _data >= _blk.begin() ) && ( _data + (_n-1)*_ld + _m-1 < _blk.end() ) );
  }
};

template <typename T>
std::ostream& operator<<( std::ostream& out, const Mat<T> A )
{
  out << "[";
  for( typename Mat<T>::idx_t i=0; i<A.n(); ++i )
    out << A.col(i) << "," << std::endl;
  out << "]";
  return out;
}

template <typename T>
Vec<T> diag( const Mat<T> A )
{
  return A.diag();
}

#endif  // _MAT_HPP
