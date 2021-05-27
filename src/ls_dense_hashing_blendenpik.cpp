// Author: Zhen Shao
#include "lsex_all_in_c.h"
#include <stdlib.h>
#include "Config.hpp"
#include <cmath>
#include "RandSolver.hpp"
#include "LinAlg.hpp"
#include "IterSolver.hpp"
#include "SpMat.hpp"
#include "cs.h"         
#include <complex>
#include "cholmod.h"
#include "SuiteSparseQR.hpp"
#include "cholmod_core.h"
#include "fblas.h"
#include "flapack.h"
#include "blendenpik.hpp"

#include "ski-lls.h"

// C++ solver for dense linear least squares using hashing
// Solving the linear least sqaures min_x |Ax-b|_2
void ls_dense_hashing_blendenpik( 
    Mat_d & A, /* m by n */
    Vec_d& b, /* m by 1 */
    Vec_d& x, /* solution, n by 1*/
    long& rank, /* detected rank (in case CPQR is used, otherwise equals to n*/
    int& flag, /* LSQR convergence flag (0=not convergent, 1=convergent) */
    long& it, /* LSQR iteration count */
    double gamma, /* over-sampling ratio */
    long k, /* nnz per column in the hashing matrix */
    double abs_tol, /* absolute tolerance for the residual */
    double rcond, /* control minimal diagonal entries of R_11 */
    double it_tol, /* LSQR relative tolerance */
    long max_it, /* LSQR max iteration */
    int debug, /* no use */
    int wisdom /* flag that if fftw wisdom is to be used */
    )
{
    assert( A.m() == b.n() ); 
    assert( A.n() == x.n() );
    if (A.m()<A.n()) {
        throw std::runtime_error( "Matrix is not overdetermined." );
    }
    long n = A.n();

    // sketch
    Vec_d c(ceil(A.n()*gamma)); // RHS to be sketched
    Mat_d As = rand_hashing_dct(A, b, c, gamma, k, wisdom);

    // build preconditioner, test explicit sketching solution
          long* E;
          E = (long*) calloc(As.n(), sizeof(long));
          double *tau = (double*) calloc(n, sizeof(double));
          // --------------- Compute CPQR factorization----------------//
          column_pivoted_qr(As, E, tau);

          // --------------- Compute Q^T c, store the result in c------//
          char left= 'L';
          char trans = 'T';
          char no_trans = 'N';
          long one = 1;
          long workspace_size, info;
          double *workspace;
          double wsize_d;

          /* Query workspace size */
          workspace_size = -1;
          workspace = &wsize_d;
          dormqr(&left, &trans, &c.n(), &one, &As.n(), As.data(), &As.ld(), tau, c.data(), &c.n(), \
            workspace, &workspace_size, &info);  
          /* Compute */
          workspace_size = (long)wsize_d;
          workspace = (double *)malloc(sizeof(double) * workspace_size);
          dormqr(&left, &trans, &c.n(), &one, &As.n(), As.data(), &As.ld(), tau, c.data(), &c.n(), \
            workspace, &workspace_size, &info);  

          // ----------------------- get the R and the rank --------//
            rank=0;
            for (long i=0; i< As.n(); i++){
              if (fabs( As(i,i) ) > rcond)
              {
                rank++;
              }
              else
              {
                break;
              }
            }
            Mat_d R;
            R = As.submat(0, rank, 0, rank);  

          // ----------------- Calculate the least square solution by back substitution ---------//
          char up = 'u';
          char diag = 'n';
          dtrsv(&up, &no_trans, &diag, &R.n(), R.data(), &R.ld(), c.data(), &c.inc());  

          // ---------- Get the basic solution ---------------------------- //

          for (long i=0; i < rank; i++){
              x.data()[E[i]-1] = c.data()[i]; 
          }

          for (long i = rank; i<As.n(); i++){
              x.data()[E[i]-1] = 0;
          }  
          // ----------- Test residual ----------------------------- //
          Vec_d rs = b.copy();
          A.mv('n', -1, x, 1, rs);
          std::cout << "Explicit Sketching Result is: "<<  nrm2(rs) << std::endl;

          if (nrm2(rs) < abs_tol){
            std::cout << "Return solution found by explicit sketching with residual: "<< nrm2(rs) << std::endl;
            it = 0;
            flag = 0;
            return;
          }
    Vec_d xs = x.copy();

    // subselect columns of A
    long forward = 1;
    long backward = 0;
    dlapmt(&forward, &A.m(), &A.n(), A.data(),&A.ld(), E);
    Mat_d Areduced;
    Areduced = A.submat(0, A.m(), 0, rank);
    Vec_d y(rank);

    // preconditioned lsqr (dense preconditioner)
    lsqr_dense_pre(Areduced, rs, it_tol , max_it, y, flag, it, R, debug);
    dtrsv(&up, &no_trans, &diag, &R.n(), R.data(), &R.ld(), y.data(), &y.inc());
    // recover original solution by doing the permutation

    for (long i=0; i < rank; i++){
        x.data()[E[i]-1] = y.data()[i]; 
    }

    for (long i = rank; i<A.n(); i++){
        x.data()[E[i]-1] = 0;
    }

    // add the explicit sketching solution xs (initial guess)
    for (long i=0; i< x.n(); i++){
      x(i) = xs(i) + x(i); 
    }

    // bring A back to its original form
    dlapmt(&backward, &A.m(), &A.n(), A.data(),&A.ld(), E);
    free(E);
}
