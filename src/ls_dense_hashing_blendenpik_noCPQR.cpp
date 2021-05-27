// Author: Zhen Shao
// Note, this function does not check the explicitly sketched solution
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
#include "blendenpik.hpp"

#include "ski-lls.h"

// C++ solver for dense linear least squares using hashing
// Solving the linear least sqaures min_x |Ax-b|_2
void ls_dense_hashing_blendenpik_noCPQR( 
    Mat_d & A, /* m by n */
    Vec_d& b, /* m by 1 */
    Vec_d& x, /* solution, n by 1*/
    long& rank, /* detected rank (in case CPQR is used, otherwise equals to n*/
    int& flag, /* LSQR convergence flag (0=not convergent, 1=convergent) */
    long& it, /* LSQR iteration count */
    double gamma, /* over-sampling ratio */
    long k, /* nnz per column in the hashing matrix */
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

    // sketch
    Mat_d As = rand_hashing_dct(A, gamma, k, wisdom);

    // build preconditioner,
    dense_qr(As); // overwriting As
    Mat_d R;
    R = As.submat(0, As.n(), 0, As.n()); // truncation, assuming n<m
    rank = R.n();

    // preconditioned lsqr (dense preconditioner)
    lsqr_dense_pre(A, b, it_tol , max_it, x, flag, it, R, debug);

    char up = 'u';
    char diag = 'n';
    char no_trans = 'n';
    dtrsv(&up, &no_trans, &diag, &R.n(), R.data(), &R.ld(), x.data(), &x.inc());
}