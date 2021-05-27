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

// C++ solver for dense linear least squares using LAPACK QR
// Solving the linear least sqaures min_x |Ax-b|_2
void ls_dense_lapack_qr( 
    Mat_d & A, /* m by n */
    Vec_d& b, /* m by 1 */
    Vec_d& x /* solution, n by 1*/
    )
{
    assert( A.m() == b.n() ); 
    assert( A.n() == x.n() );
    if (A.m()<A.n()) {
        throw std::runtime_error( "Matrix is not overdetermined." );
    }
    long m = A.m();
    long n = A.n();

    // Dense QR
      Mat_d Acopy = A.copy();
      x = b.copy();
      double *tau;
      long workspace_size, info;
      double *workspace;
      double wsize_d;
      tau = (double*)malloc(n * sizeof(double));

      /* Query workspace size */
      workspace_size = -1;
      workspace = &wsize_d;
      dgeqrf(&m, &n, Acopy.data(), &m, tau, workspace, &workspace_size, &info);

      /* Compute */

      workspace_size = (long)wsize_d;
      workspace = (double *)malloc(sizeof(double) * workspace_size);
      dgeqrf(&m, &n, Acopy.data(), &m, tau, workspace, &workspace_size, &info);
      free(workspace);


    // Compute least square solution
          // --------------- Compute Q^T c, store the result in c------//
          char left= 'L';
          char trans = 'T';
          char no_trans = 'N';
          long one = 1;

          /* Query workspace size */
          workspace_size = -1;
          workspace = &wsize_d;
          dormqr(&left, &trans, &x.n(), &one, &Acopy.n(), Acopy.data(), &Acopy.ld(), tau, x.data(), &x.n(), \
            workspace, &workspace_size, &info);  
          /* Compute */
          workspace_size = (long)wsize_d;
          workspace = (double *)malloc(sizeof(double) * workspace_size);
          dormqr(&left, &trans, &x.n(), &one, &Acopy.n(), Acopy.data(), &Acopy.ld(), tau, x.data(), &x.n(), \
            workspace, &workspace_size, &info);  

          // ----------------- Calculate the least square solution by back substitution ---------//
          char up = 'u';
          char diag = 'n';
          dtrsv(&up, &no_trans, &diag, &Acopy.n(), Acopy.data(), &Acopy.ld(), x.data(), &x.inc());  

          x = x.subvec(0,n);

    free(workspace);
}