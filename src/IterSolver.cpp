// @author Xiangrui Meng
// Extended by Zhen Shao
// Fixed a small bug, allow sparse and incomplete Cholesky preconditioner to be used

#include <cassert>
#include <cmath>
#include "LinOp.hpp"
#include "LinAlg.hpp"
#include "IterSolver.hpp"
#include "cs.h"
#include "SuiteSparseQR.hpp"
#include "cholmod.h"
#include "fblas.h"
#include "fhsl_mi35.h"
#include "lsex_all_in_c.h"

using namespace std;
int PRINTING=0;

// ----------- SPARSE PRE --------------//
template <typename T>
class LinOP_sparse_preconditioner
    : public LinOp<T>
{
 public:
  typedef long idx_t;
  LinOP_sparse_preconditioner( LinOp<T>& A, cholmod_sparse* R, 
  int success=1, double perturb=1e-6)
      : _A(A) , _R(R), _n(_R->ncol), _m(_A.m()), 
        _success(success), _perturb(perturb)
  {
    assert( A.n() == (long) R->nrow );
    tmp = Vec<T>(A.n());
  }
  virtual idx_t const& m() const { return _m;}
  virtual idx_t const& n() const { return _n;}
  
  // overwrite matrix-vector multiplication by using preconditioner               
  virtual void mv( const char trans, const T alpha, const Vec<T> x, const T beta, Vec<T> y )
  {
    if( trans == 'n' )
    {
      tmp = x.copy();
      cs_usolve_cholmod_structure(_R, tmp.data(), _success, _perturb);
      _A.mv( 'n', alpha, tmp, beta, y );
    }
    else
    {
      _A.mv( 't', 1.0, x, 0.0, tmp );
      cs_utsolve_cholmod_structure(_R, tmp.data(), _success, _perturb);
      scal(beta, y);
      axpy(alpha, tmp, y);
    }
  }                
 private:
  Vec_d     tmp;
  LinOp<T>&  _A;
  cholmod_sparse* _R;
  idx_t _n;
  idx_t _m;
  int _success;
  double _perturb;
};

// ---------------------------- DENSE PRE -----------------//
template <typename T>
class LinOP_dense_preconditioner
    : public LinOp<T>
{
 public:
  typedef long idx_t;
  LinOP_dense_preconditioner( LinOp<T>& A, Mat_d& R)
      : _A(A) , _R(R), 
      _up('u'), _trans('t'), _no_trans('n'), _diag('n')
  {
    assert( A.n() == R.m());
    tmp = Vec<T>(A.n());
  }
  virtual idx_t const& m() const { return _A.m();}
  virtual idx_t const& n() const { return _R.n();}
  
  // overwrite matrix-vector multiplication by using preconditioner               
  virtual void mv( const char trans, const T alpha, const Vec<T> x, const T beta, Vec<T> y )
  {

    if( trans == 'n' )
    {
      tmp = x.copy();
      dtrsv(&_up, &_no_trans, &_diag, &_R.n(), _R.data(), &_R.ld(), tmp.data(), &tmp.inc());
      _A.mv( 'n', alpha, tmp, beta, y );
    }
    else
    {
      _A.mv( 't', 1.0, x, 0.0, tmp );
      dtrsv(&_up, &_trans, &_diag, &_R.n(), _R.data(), &_R.ld(), tmp.data(), &tmp.inc());
      scal(beta, y);
      axpy(1.0, tmp, y);
    }
  }                
 private:
  Vec_d     tmp;
  LinOp<T>&  _A;
  Mat_d& _R;
  char _up;
  char _trans; 
  char _no_trans; 
  char _diag;
};

// ---------- IC PRE ----------------------------//
#ifdef HAVE_HSL

template <typename T>
class LinOP_ic_preconditioner
    : public LinOp<T>
{
 public:
  typedef long idx_t;
  LinOP_ic_preconditioner( LinOp<T>& A, void* pkeep)
      : _A(A) , _pkeep(pkeep),
      _trans(1), _no_trans(0), _ifail(0)
  {
    // assert( A.n() == R.m()); cannot do the assertion because pkeep is difficult to handle
    tmp = Vec<T>(A.n());
  }
  virtual idx_t const& m() const { return _A.m();}
  virtual idx_t const& n() const { return _A.n();} // assuming preconditioning size match
  
  // overwrite matrix-vector multiplication by using preconditioner               
  virtual void mv( const char trans, const T alpha, const Vec<T> x, const T beta, Vec<T> y )
  {

    if( trans == 'n' )
    {
      // Note IC returns Lower trianglar, so transpose or no transpose is reversed
      hsl_mi35_solve(&_trans, &_A.n(), &_pkeep, x.data() , tmp.data() , &_ifail); 
      _A.mv( 'n', alpha, tmp, beta, y );
    }
    else
    {
      _A.mv( 't', 1.0, x, 0.0, tmp );
      hsl_mi35_solve(&_no_trans, &_A.n(), &_pkeep, tmp.data(), tmp.data(), &_ifail);
      scal(beta, y);
      axpy(1.0, tmp, y);
    }
  }                
 private:
  Vec_d     tmp;
  LinOp<T>&  _A;
  void* _pkeep;
  fint _trans;
  fint _no_trans;
  fint _ifail;
};

#endif

// -------------------------- Vanilla LSQR --------------------------//

void lsqr( LinOp<double>& A, const Vec_d b,
           const double tol, const long maxit,
           Vec_d& x, int& flag, long& it, int debug ) 
{
  assert( A.m() == b.n() );
  assert( A.n() == x.n() );

  // if ( debug==1 && (A.m()>100.0 || A.n()>100.0) ){
  //   std::cout << "Matrix too large, turning off debug mode..." << std::endl;
  //   debug=0;
  // }


  Vec_d v( A.n() );
  Vec_d w(v.n());
  double alpha=0;

  // explicitly make sure it is all zero...
  for (long i=0; i<A.n(); i++){
    v(i) = 0;
  }

  Vec_d u = b.copy();
  double beta = nrm2(u);
  if(beta > 0){
    scal(1.0/beta,u);    
    if (debug==1)
    {
      cout << "A.n(): " << A.n() << endl;
      cout << "u: " << u.subvec(0, min(5, (int)u.n())) << endl;
      cout << "v: " << v.subvec(0, min(5,(int)v.n())) << endl;
      A.mv( 't', 1.0, u, 0.0, v );
      cout << "u: " << u.subvec(0, min(5, (int)u.n())) << endl;
      cout << "v: " << v.subvec(0, min(5,(int)v.n())) << endl << endl;
    } else{
      A.mv( 't', 1.0, u, 0.0, v );
    }
    alpha = nrm2(v); // norm of v might be zero!!!
  }

  if (alpha >0){
    scal(1.0/alpha,v);    
    w= v.copy();
  }

  x = 0.0;
  double nrm_ar = alpha*beta;
  if (nrm_ar ==0){
    it=0;
    flag=0; // converge in 0 iteration
    return;
  } 
  
  double phi, rho;
  double cs, sn;
  double phibar = beta;
  double rhobar = alpha;
  double theta;
  
  double nrm_a    = 0.0;
  double nrm_r;
  
  // double nrm_ar_0 = alpha*beta;
   flag =1;
   it = 0;
  for( long k = 0; k < maxit; ++k )
  {
    it++;

    if (debug==1)
    {
      cout << "u: " << u.subvec(0, min(5,(int)u.n())) << endl;
      cout << "v: " << v.subvec(0, min(5,(int)v.n())) << endl;
      A.mv( 'n', 1.0, v, -alpha, u );
      cout << "u: " << u.subvec(0, min(5,(int)u.n())) << endl;
      cout << "v: " << v.subvec(0, min(5,(int)v.n())) << endl << endl;
    } else{
      A.mv( 'n', 1.0, v, -alpha, u );
    }

    beta = nrm2(u);

    if (beta>0){
      scal(1.0/beta,u);      
      nrm_a = sqrt( nrm_a*nrm_a + alpha*alpha + beta*beta );

      if (debug==1){
        cout << "u: " << u.subvec(0, min(5,(int)u.n())) << endl;
        cout << "v: " << v.subvec(0, min(5,(int)v.n())) << endl;
        A.mv( 't', 1.0, u, -beta, v );
        cout << "u: " << u.subvec(0, min(5,(int)u.n())) << endl;
        cout << "v: " << v.subvec(0, min(5,(int)v.n())) << endl << endl;
      } else{
        A.mv( 't', 1.0, u, -beta, v );
      }

      alpha = nrm2(v);
      if (alpha >0){
        scal(1.0/alpha,v);    
      }
    }

    rho    = sqrt( rhobar*rhobar + beta*beta );
    cs     = rhobar/rho;
    sn     = beta/rho;
    theta  = sn*alpha;
    rhobar = -cs*alpha;
    phi    = cs*phibar;
    phibar = sn*phibar;

    axpy( phi/rho, w, x );

    scal( -theta/rho, w );
    axpy( 1.0, v,  w );

    nrm_r  = phibar;
    nrm_ar = phibar*alpha*fabs(cs);

    if( nrm_ar < tol*nrm_a*nrm_r ){
        flag = 0;
        break;
    }

  }
}

// ---------------------- Preconditioned version, sparse preconditioner ------------------------------//
void lsqr( LinOp<double>& A, const Vec_d b,
           const double tol, const long maxit,
           Vec_d& x, int& flag, long& it, cholmod_sparse* R_11, 
           int success, double perturb, int debug) // LSQR with preconditioner
{
  LinOP_sparse_preconditioner<double> Apre(A, R_11, success, perturb);
  lsqr( Apre, b, tol, maxit, x, flag, it, debug );
}


// ----------------------- Preconditioned lsqr, dense preconditioner ----------------------//

void lsqr_dense_pre( LinOp<double>& A, const Vec_d b,
           const double tol, const long maxit,
           Vec_d& x, int& flag, long& it, Mat_d& R, int debug) // LSQR with preconditioner
{
  LinOP_dense_preconditioner<double> Apre(A, R);
  lsqr( Apre, b, tol, maxit, x, flag, it, debug);
}


// ------------------------ Preconditioned lsqr, using incomplete cholesky ----------//
#ifdef HAVE_HSL
void lsqr_ic( LinOp<double>& A, const Vec_d b,
           const double tol, const long maxit,
           Vec_d& x, int& flag, long& it, void* pkeep, int debug) // LSQR with ic preconditioner
{
  LinOP_ic_preconditioner<double> Apre(A, pkeep);
  lsqr( Apre, b, tol, maxit, x, flag, it, debug);
}
#endif




