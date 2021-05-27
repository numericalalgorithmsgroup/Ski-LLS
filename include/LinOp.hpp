// @author Xiangrui Meng

#ifndef _LINOP_HPP
#define _LINOP_HPP

#include <iostream>

#include "Config.hpp"
#include "Vec.hpp"
#include "Mat.hpp"
#include "SpMat.hpp"

template <typename T>
class Mat;

template <typename T>
class LinOp // mv: matrix vector multiply, mm: matrix matrix multiply
{
 public:
  
  typedef T    val_t;
  typedef long idx_t; 
  
  virtual idx_t const & m() const = 0; 
  virtual idx_t const & n() const = 0; 
  virtual void mv( const char trans, const T alpha, const Vec<T> x, const T beta, Vec<T> y ) = 0; 
  virtual void mm( const char trans, const char trans_b, const T alpha, const Mat<T> B, const T beta, Mat<T> C )
  {
    if( trans_b == 'n' )
    {
      for( idx_t j=0; j<C.n(); ++j )
        mv( trans, alpha, B.col(j), beta, C.col(j) );
    }
    else
    {
      for( idx_t j=0; j<C.n(); ++j )
        mv( trans, alpha, B.row(j).copy(), beta, C.col(j) );
    } 
  }
  
  virtual void mv( const char trans, const Vec<T> x, Vec<T> y )
  {
    mv( trans, 1, x, 0, y );
  }

  virtual void mm( const char trans, const char trans_b, const Mat<T> B, Mat<T> C )
  {
    mm( trans, trans_b, (T) 1, B, (T) 0, C );
  }
  
  virtual ~LinOp() {}
};

template <typename T>
class Composite_LinOP
    : public LinOp<T>
{
 public:

  typedef long idx_t;
  
  Composite_LinOP( const char op,
                   const char trans_a, LinOp<T>& A,
                   const char trans_b, LinOp<T>& B )
      : _op(op)
      , _trans_a(trans_a), _A(A)
      , _trans_b(trans_b), _B(B)
  {
    if( _op == '*' )
    {
      if( trans_a == 'n' && trans_b == 'n' )
      {
        assert( _A.n() == _B.m() );
        tmp = Vec<T>(_A.n());
      }
      else if( trans_a == 'n' && trans_b == 't' )
      {
        assert( _A.n() == _B.n() );
        tmp = Vec<T>(_A.n());
      }
      else if( trans_a == 't' && trans_b == 'n' )
      {
        assert( _A.m() == _B.m() );
        tmp = Vec<T>(_A.m());
      }
      else
      {
        assert( _A.m() == _B.n() );
        tmp = Vec<T>(_A.m());
      }
    }
  }

  virtual idx_t const& m() const
  {
    if( _trans_a == 'n' )
      return _A.m();
    else
      return _A.n();
  }
  
  virtual idx_t const& n() const
  {
    if( _trans_b == 'n' )
      return _B.n();
    else
      return _B.m();
  }
                         
  virtual void mv( const char trans, const T alpha, const Vec<T> x, const T beta, Vec<T> y )
  {
    if( _op == '*' )
    {
      if( trans == 'n' )
      {
        _B.mv( _trans_add( trans, _trans_b ), 1.0, x, 0.0, tmp );
        _A.mv( _trans_add( trans, _trans_a ), alpha, tmp, beta, y );
      }
      else
      {
        _A.mv( _trans_add( trans, _trans_a ), 1.0, x, 0.0, tmp );
        _B.mv( _trans_add( trans, _trans_b ), alpha, tmp, beta, y );
      }
    }
  }
                         
 private:

  char _trans_add( const char trans_1, const char trans_2 )
  {
    if( trans_1 == trans_2 )
      return 'n';
    else
      return 't';
  }
                         
  char      _op;
  char      _trans_a;
  LinOp<T>& _A;
  char      _trans_b;
  LinOp<T>& _B;

  Vec_d     tmp;
};

#endif  // _LINOP_HPP
