/*
// Author: Xingrui Meng
*/
#include <cstddef>
#include <cassert>
#include <cmath>

#include <stdexcept>
#include <algorithm>

#include "Utils.hpp"
#include "LinOp.hpp"
#include "LinAlg.hpp"
#include "Random.hpp"
#include "IterSolver.hpp"
#include "RandSolver.hpp"

void lsrn( LinOp<double>& A, Vec_d& b, double rcond, Vec_d& x, long& rank, int& flag, long& it, 
  double gamma, double it_tol, long max_it) //where lsrn is defined
{
  assert( A.m() == b.n() ); // input: A,b, parameter rcond, return by reference, x and rank
  assert( A.n() == x.n() );

  typedef ptrdiff_t idx_t;
  
  idx_t m = A.m();
  idx_t n = A.n();
  
  if( m > gamma*n )
  {
    idx_t s = ceil(gamma*n);
    Mat_d N = rnpre( A, s, rcond );
    rank = N.n();
    Composite_LinOP<double> AN( '*', 'n', A, 'n', N );
    gamma = 1.0*s/rank;
    //long maxit = ceil( 2.0 * ( log(2.0) - log(tol) ) / log(gamma) );
    Vec_d y(rank);
    lsqr( AN, b, it_tol, max_it, y, flag, it);
    mv( 'n', 1.0, N, y, 0.0, x );
  }
  else if( n > gamma*m )
  {
    idx_t s = ceil(gamma*m); // so s is the final sketched dimension
    Mat_d M = rnpre( A, s, rcond );
    rank = M.n();
    Composite_LinOP<double> MtA( '*', 't', M, 'n', A );
    Vec_d Mtb(rank);
    M.mv( 't', 1.0, b, 0.0, Mtb );
    gamma = 1.0*s/rank;
    lsqr( MtA, Mtb, it_tol, max_it, x, flag, it, 0);
  }
  else
  {
    throw std::runtime_error( "Matrix is not thin enough or fat enough." ); // only accepts sketchable matrices...
  }
}

Mat_d rnpj( LinOp<double>& A, ptrdiff_t s ) // called in rnpre for m>n, input: A, pass by reference, s, probably desired size
{
  typedef ptrdiff_t idx_t;

  idx_t m = A.m();
  idx_t n = A.n();
  
  idx_t blk_sz = 512;

  if( m >= n )
  {
    blk_sz = std::min( blk_sz, s );
    Mat_d G;
    while( G.empty() )
    {
      try
      {
        G = Mat_d( m, blk_sz );
      }
      catch( std::bad_alloc& )
      {
        blk_sz /= 2;
        if( blk_sz == 0 )
          throw;
      }
    }
    Mat_d As_T( n, s );
    for( int j=0; j<s; j+=blk_sz )
    {
      idx_t len = std::min( j+blk_sz, s ) - j; // length to be generated
      randn( m*len, G.data() ); //using the parallel random_normal generation to generate the sketching matrix
      // generating m*len random normals
      Mat_d sub_As_T( n, len, As_T, 0, j );
      Mat_d sub_G( m, len, G, 0, 0 );
      
      A.mm( 't', 'n', 1.0, sub_G, 0.0, sub_As_T ); //doing multiplication on the fly
    }

    return As_T;
  } // underdetermined case
  else if( m < n ) //over determined case, what we are interested!
  {
    blk_sz = std::min( blk_sz, s );
    Mat_d G;
    while( G.empty() )
    {
      try
      {
        G = Mat_d( n, blk_sz );
      }
      catch( std::bad_alloc& )
      {
        blk_sz /= 2;
        if( blk_sz == 0 )
          throw;
      }
    }
    Mat_d As( m, s );

    for( int j=0; j<s; j+=blk_sz )
    {
      idx_t len = std::min( j+blk_sz, s ) - j;
      randn( n*len, G.data() );

      Mat_d sub_As( m, len, As, 0, j );
      Mat_d sub_G( n, len, G, 0, 0 );

      A.mm( 'n', 'n', 1.0, sub_G, 0.0, sub_As ); // a simple matrix matrix multiplication
    }

    return As; //As, the sketched matrix
  }

  return Mat_d();
}

Mat_d rnpre( LinOp<double>& A, ptrdiff_t s, double rcond ) // generate preconditioner
{ // input: A,
  typedef ptrdiff_t idx_t; //ptrdiff: type of difference between two pointers...

  idx_t m = A.m();
  idx_t n = A.n();

  if( m >= n )
  {
    Mat_d As_T = rnpj( A, s );
    Vec_d sgm(n);
    Mat_d V(n,n);

    svd( As_T, sgm, V );

    double r_tol = rcond*sgm(0);
    idx_t rank = 0;
    for( idx_t i=0; i<n; ++i )
    {
      if( sgm(i) > r_tol )
        rank++;
      else
        break;
    }

    for( idx_t j=0; j<rank; ++j )
      scal(1.0/sgm(j), V.col(j));

    Mat_d N(n,rank,V);

    return N;
  }
  else
  {
    Mat_d As = rnpj( A, s );
    Vec_d sgm(m);
    Mat_d U(m,m);

    svd( As, sgm, U ); // so As is the sketched matrix

    double r_tol = rcond*sgm(0);
    idx_t rank = 0;
    for( idx_t i=0; i<m; ++i )
    {
      if( sgm(i) > r_tol )
        rank++;
      else
        break;
    }

    for( idx_t j=0; j<rank; ++j )
      scal(1.0/sgm(j), U.col(j));
    
    Mat_d M(m,rank,U);

    return M;
  }
}