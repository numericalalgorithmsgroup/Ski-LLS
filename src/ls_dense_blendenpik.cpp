//Copyright (c) 2009-2016, Haim Avron and Sivan Toledo.
// C++ interface written by Zhen Shao
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

// C++ solver interface for Avron's Blendenpik algorithm
// Solving the linear least sqaures min_x |Ax-b|_2
void ls_dense_blendenpik( 
    Mat_d & A, /* m by n */
    Vec_d& b, /* m by 1 */
    Vec_d& x, /* solution, n by 1*/
    long& rank, /* detected rank (in case CPQR is used, otherwise equals to n*/
    int& flag, /* LSQR convergence flag (0=not convergent, 1=convergent) */
    long& it, /* LSQR iteration count */
    double gamma, /* over-sampling ratio */
    double it_tol, /* LSQR relative tolerance */
    long max_it, /* LSQR max iteration */
    int debug, /* no use*/
    int wisdom /* flag that if fftw wisdom is to be used */
    )
{
    // input check
    assert( A.m() == b.n() ); 
    assert( A.n() == x.n() );
    if (A.m()<A.n()) {
        throw std::runtime_error( "Matrix is not overdetermined." );
    }
    double t1;
    double t2;

    // sketch
    t1 = read_wtime();
    Mat_d As = rand_subsampled_dct(A, gamma, wisdom);
    t2 = read_wtime();
    std::cout << "Sketching time is: " << t2-t1 << std::endl; 
  
    // build preconditioner,
    t1 = read_wtime();
    dense_qr(As); // overwriting As
    Mat_d R;
    R = As.submat(0, As.n(), 0, As.n()); 
    rank = R.n();
    t2 = read_wtime();
    std::cout << "Build preconditioner time is: " << t2-t1 << std::endl;  


    // preconditioned lsqr (dense preconditioner)
    t1 = read_wtime();
    lsqr_dense_pre(A, b, it_tol , max_it, x, flag, it, R, debug);
    t2 = read_wtime();
    std::cout << "LSQR time is: "<<t2-t1 << std::endl;

    char up = 'u';
    char diag = 'n';
    char no_trans = 'n';
    dtrsv(&up, &no_trans, &diag, &R.n(), R.data(), &R.ld(), x.data(), &x.inc());
}