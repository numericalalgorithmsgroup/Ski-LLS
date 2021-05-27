/**
 * @file   SpMat.hpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Mon Oct  3 10:28:27 2011
 * 
 * @brief  
 * 
 * 
 */

#ifndef _SPMAT_HPP
#define _SPMAT_HPP

#include "Config.hpp"
#include "LinOp.hpp"

template <typename T>
class LinOp;

template <typename T>
class Mat;

template <typename T>
class CSC_Mat;

typedef CSC_Mat<double> CSC_Mat_d; 
typedef CSC_Mat<long>   CSC_Mat_l;

template <typename T>
class CSR_Mat;

typedef CSR_Mat<double> CSR_Mat_d;
typedef CSR_Mat<long>   CSR_Mat_l;

template <typename T>
void csc_to_csr( const CSC_Mat<T> C, CSR_Mat<T> R );

template <typename T>
class SpMat_Op;

template <typename T>
void mv( const char trans, const T alpha, const CSC_Mat<T> A,
         const Vec<T> x, const T beta, Vec<T> y );

template <typename T>
void mm( const char trans_a, const char trans_b, const T alpha,
         const CSC_Mat<T> A, const Mat<T> B, const T beta, Mat<T> C );

template <typename T>
void mv( const char trans, const T alpha, const CSR_Mat<T> A,
         const Vec<T> x, const T beta, Vec<T> y );

template <typename T>
void mm( const char trans_a, const char trans_b, const T alpha,
         const CSR_Mat<T> A, const Mat<T> B, const T beta, Mat<T> C );

template <typename T>
void csc_to_csr( const CSC_Mat<T> C, CSR_Mat<T> R );

template <typename T>
class CSC_Mat
    : public LinOp<T> 
{
 public:

  typedef T          val_t;
  typedef long       idx_t;
  typedef CSC_Mat<T> self_t; 

  CSC_Mat( const idx_t m = 0, const idx_t n = 0, const idx_t nzmax = 0 )
      : _ir_blk(nzmax),
        _jc_blk(n+1),
        _data_blk(nzmax) 
  {
    _m        = m; 
    _n        = n;
    _nzmax    = nzmax;

    _ir   = _ir_blk.begin();
    _jc   = _jc_blk.begin();
    _data = _data_blk.begin(); 

    for( idx_t j=0; j<n+1; ++j )
      _jc[j] = 0;  
  }

  /** 
   * Contruct a CSC_Mat from existing data
   * 
   * @param m 
   * @param n
   * @param nzmax
   * @param ir 
   * @param jc 
   * @param data 
   * @param data_blk 
   * 
   */
  CSC_Mat( const idx_t m, const idx_t n, const idx_t nzmax, 
           idx_t * const ir, idx_t * const jc, T * const data,
           const Blk<idx_t> ir_blk = Blk<idx_t>(),
           const Blk<idx_t> jc_blk = Blk<idx_t>(),
           const Blk<val_t> data_blk = Blk<val_t>() )
      : _ir_blk(ir_blk),
        _jc_blk(jc_blk),
        _data_blk(data_blk)   
  {
    _m     = m;
    _n     = n;
    _nzmax = nzmax;
    
    _ir   = ir;
    _jc   = jc;
    _data = data;
  }  

  virtual idx_t const& m() const { return _m; }
  virtual idx_t const& n() const { return _n; }

  idx_t const& nnz() const { return _jc[_n]; } 
  idx_t const& nzmax() const { return _nzmax; }

  idx_t * const& ir() const { return _ir; }
  Blk<idx_t> const& ir_blk() const { return _ir_blk; }

  idx_t * const& jc() const { return _jc; }
  Blk<idx_t> const& jc_blk() const { return _jc_blk; }

  T* const& data() const { return _data; }
  Blk<val_t> const& data_blk() const { return _data_blk; } 
  
  virtual void mv( const char trans, const T alpha,
                   const Vec<T> x, const T beta, Vec<T> y )
  {
    return ::mv( trans, alpha, *this, x, beta, y ); 
  }

  virtual void mm( const char trans, const char trans_b, const T alpha,
                   const Mat<T> B, const T beta, Mat<T> C )
  {
     return ::mm( trans, trans_b, alpha, *this, B, beta, C );
  } 

  /** 
   * Get submatrix formed by columns [begin,end)
   * 
   * @param begin 
   * @param end 
   * 
   * @return 
   */
  self_t cols( const idx_t begin, const idx_t end ) const
  {
    self_t A( _m, end-begin, _jc[end]-_jc[begin], 
              _ir, _jc+begin, _data, 
              _ir_blk, _jc_blk, _data_blk ); 
    return A;
  }

  CSR_Mat<T> to_csr() const
  {
    CSR_Mat<T> R( _m, _n, nnz() );
    csc_to_csr<T>( *this, R );
    return R;
  }
  
 private:

  idx_t _m; 
  idx_t _n;
  idx_t _nzmax; 

  idx_t  *_ir;
  idx_t  *_jc;
  val_t  *_data;

  Blk<idx_t> _ir_blk; 
  Blk<idx_t> _jc_blk; 
  Blk<val_t> _data_blk; 

};

template <typename T>
class CSR_Mat 
    : public LinOp<T>
{
 public:

  typedef T          val_t;
  typedef long       idx_t;
  typedef CSR_Mat<T> self_t;

  CSR_Mat( const idx_t m = 0, const idx_t n = 0, const idx_t nzmax = 0 )
      : _ir_blk(m+1),
        _jc_blk(nzmax),
        _data_blk(nzmax) 
  {
    _m        = m;
    _n        = n;
    _nzmax    = nzmax;

    _ir   = _ir_blk.begin();
    _jc   = _jc_blk.begin();
    _data = _data_blk.begin();

    for( idx_t i=0; i<m+1; ++i )
      _ir[i] = 0;
  }

  /** 
   * Contruct a CSR_Mat from existing data
   * 
   * @param m 
   * @param n
   * @param nzmax
   * @param ir 
   * @param jc 
   * @param data 
   * @param data_blk 
   * 
   */
  CSR_Mat( const idx_t m, const idx_t n, const idx_t nzmax, 
           idx_t * const ir, idx_t * const jc, T * const data,
           const Blk<idx_t> ir_blk = Blk<idx_t>(),
           const Blk<idx_t> jc_blk = Blk<idx_t>(),
           const Blk<val_t> data_blk = Blk<val_t>() )
      : _ir_blk(ir_blk),
        _jc_blk(jc_blk),
        _data_blk(data_blk)
  {
    _m     = m;
    _n     = n;
    _nzmax = nzmax;
    
    _ir   = ir;
    _jc   = jc;
    _data = data;
  }

  virtual idx_t const& m() const { return _m; }
  virtual idx_t const& n() const { return _n; }

  idx_t const& nnz() const { return _ir[_m]; }
  idx_t const& nzmax() const { return _nzmax; }

  idx_t * const& ir() const { return _ir; }
  Blk<idx_t> const& ir_blk() const { return _ir_blk; }

  idx_t * const& jc() const { return _jc; }
  Blk<idx_t> const& jc_blk() const { return _jc_blk; }

  T* const& data() const { return _data; }
  Blk<val_t> const& data_blk() const { return _data_blk; }
  
  virtual void mv( const char trans, const T alpha,
                   const Vec<T> x, const T beta, Vec<T> y )
  {
    return ::mv( trans, alpha, *this, x, beta, y );
  }

  virtual void mm( const char trans, const char trans_b, const T alpha,
                   const Mat<T> B, const T beta, Mat<T> C )
  {
     return ::mm( trans, trans_b, alpha, *this, B, beta, C );
  }

  /** 
   * Get submatrix formed by columns [begin,end) // formed by rows!!
   * 
   * @param begin 
   * @param end 
   * 
   * @return 
   */
  self_t rows( const idx_t begin, const idx_t end ) const
  {
    self_t A( end-begin, _n, _ir[end]-_ir[begin],
              _ir+begin, _jc, _data,
              _ir_blk, _jc_blk, _data_blk );
    return A;
  }
  
 private:

  idx_t _m;
  idx_t _n;
  idx_t _nzmax;

  idx_t  *_ir;
  idx_t  *_jc;
  val_t  *_data;

  Blk<idx_t> _ir_blk;
  Blk<idx_t> _jc_blk;
  Blk<val_t> _data_blk;
};

template <typename T>
class SpMat_Op
    : public LinOp<T>
{
 public:

  typedef typename LinOp<T>::idx_t idx_t;
  
  SpMat_Op( const CSC_Mat<T> C )
  {
    _m = C.m();
    _n = C.n();
    _C = C;
    _R = _C.to_csr();
  }

  virtual idx_t const& m() const { return _m; }
  virtual idx_t const& n() const { return _n; }

  virtual void mv( const char trans, const T alpha, const Vec<T> x, const T beta, Vec<T> y )
  {
    if( _m >= _n )                     
      ::mv( trans, alpha, _R, x, beta, y );
    else
      ::mv( trans, alpha, _C, x, beta, y );
  }
  
  virtual void mm( const char trans, const char trans_b, const T alpha, const Mat<T> B, const T beta, Mat<T> C )
  {
    if( _m >= _n )
      ::mm( trans, trans_b, alpha, _R, B, beta, C );
    else
      ::mm( trans, trans_b, alpha, _C, B, beta, C );
  }

 private:

  long _m, _n;
  CSC_Mat<T> _C;
  CSR_Mat<T> _R;
}; 

#endif  // _SPMAT_HPP
