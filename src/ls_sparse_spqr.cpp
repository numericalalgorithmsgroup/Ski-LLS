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

#include "ski-lls.h"

// generate a preconditioner use sparse QR factorization
// return the preconditioner R_11, also return c= Q^T b
cholmod_sparse* rnpre_sparse(
    cs_dl* As /* Sketched A, to be factorized */,
    double* b /*size Am*/, 
    cholmod_dense* &c_cholmod, /* size rank(A), c = Q^T b */
    cholmod_common *cc, /* for SuiteSparse QR */
    SuiteSparse_long* &E, /* array storing the permutation of columns of As */
    int& success, /* success=1 means R is estimated to be well conditioned */
    int ordering, /* ordering used in the sparse QR factorization */
    double rcond_threshold, /* Threshold for 'ill-conditioned' of R */
    int debug /* no use */
    ) 
{ 
    success =0; // assume unsuccessfuly unless proven otherwise

    // fisrt convert As to cholmod format
    double t1= read_wtime();
    cholmod_sparse* As_cholmod = to_cholmod_sparse(As, cc);
    double t2 = read_wtime();
    std::cout << "conversion time is: " << t2-t1 << std::endl;

    // Do sparse, rank-revealing QR
    cholmod_sparse *R, *R_11;
    SuiteSparse_long est_rank;

    cholmod_dense *bcholmod = 
        cholmod_l_allocate_dense(As_cholmod->nrow, 1, As_cholmod->nrow, CHOLMOD_REAL, cc);
    bcholmod->x = b;

    est_rank = SuiteSparseQR <double>(
        ordering, SPQR_DEFAULT_TOL, 
        -1, As_cholmod, bcholmod, &c_cholmod,
        &R, &E, cc);
    if (est_rank ==0){
        throw std::runtime_error("SPQR Factorization fails");
    }

    // Truncate out detected rank-deficiency
    long sub_mat_indices[est_rank];
    for(long i=0; i<est_rank; i++){
        sub_mat_indices[i] = i;
    }
    R_11 = cholmod_l_submatrix(R, sub_mat_indices, 
        est_rank, sub_mat_indices,est_rank, 1,1,cc );

    // test condition number, in case there is undetected rank-deficiency
    t1 = read_wtime();
    double rcond_R_11;
    cholmod_cond_est(R_11, cc, rcond_R_11);
    if (rcond_R_11 > rcond_threshold) // not-ill-conditioned
        success = 1;
    t2=read_wtime();
    std::cout << "condition estimation time is: " << t2-t1 << std::endl;

    // output
    return R_11;
}


// Solvers sparse linear least squares using hashing and sparse QR
void ls_sparse_spqr(
    cs_dl & A, /* input A */
    Vec_d& b, /* input b */
    Vec_d& x, /* output x*/
    long& rank, /* detected rank */
    int& flag, /* LSQR convergence flag, 1=convergent */
    long& it, /* LSQR iteration count */
    double gamma, /* oversampling ratio */
    long k, /* nnz per column of the hashing matrix */
    double abs_tol, /* absolute residual tolerance for the solution */ 
    int ordering, /* ordering used in sparse QR*/
    double it_tol, /* relative tolerance for LSQR termination */
    long max_it, /* LSQR iteration limit */
    double rcond_threshold, /* quantify 'well-conditioned' of the preconditioner R_11 */
    double perturb, /* potential pivot perturbation is R_11 ill-conditioned */
    int debug /* no use*/
    )
{
    assert( A.m == b.n() ); 
    assert( A.n == x.n() );
    if (A.m<A.n) {
        throw std::runtime_error( "Matrix is not overdetermined." );
    }

    cs_dl* As;
    cholmod_sparse *R_11;
    double* b_sketched;
    cholmod_common Common, *cc ;
    cc = &Common;
    cholmod_l_start(cc);
    if (MAX_NUM_THREADS ==1 ) // no thread
    {
        cc->SPQR_grain = 1;
    } else{ // use recommended value, Note this may cause trouble if BLAS is OPENMP
        cc->SPQR_grain = 2 * MAX_NUM_THREADS;
        cc->SPQR_nthreads = MAX_NUM_THREADS;
    }
    SuiteSparse_long* E;
    cholmod_dense *c_cholmod;
    int success = 0;
    Vec_d r_s = b.copy();
    cholmod_sparse* Areduced;
    cholmod_sparse* A_cholmod;

    // ----------------- sketch -------------------------//
    double t1= read_wtime();
    As = sparse_sketch(A, k,  gamma, b.data(), b_sketched);
    double t2 = read_wtime();
    std::cout << "sketching time is: " << t2-t1 << std::endl;

    //  -------------- build preconditioner ------------------//
    t1 = read_wtime();
    R_11 = rnpre_sparse(As, b_sketched, c_cholmod, cc, E, 
        success, ordering, rcond_threshold, debug);
    t2 = read_wtime();
    std::cout << "build preconditioner time is: " << t2-t1 << std::endl;

    if (success==0){
        std::cout <<  "SPQR fails to detect the rank." << std::endl;
        std::cout << " preconditioner is ill conditioned " << std::endl;
    }
    rank = R_11->nrow;

    // ----------- Test the accuracy of explicit sketching solution -----//

    cs_usolve_cholmod_structure(R_11, (double*)c_cholmod->x); // multiply by R_11^{-1}

    for (long i=0; i < rank; i++){
        x.data()[E[i]] = ((double*)c_cholmod->x)[i];
    }
    for (long i = rank; i<A.n; i++){
        x.data()[E[i]] = 0;
    } // multiply by P_1^T (see paper)

    CSC_Mat_d A_CSC(A.m, A.n, A.nzmax, 
        (long*)A.i, (long*)A.p, (double*)A.x);
    A_CSC.mv('n', -1, x, 1, r_s);
    double r = nrm2(r_s);

    if (r < abs_tol){
        std::cout << "return solution found by explicit sketching with residual: "
            << r << std::endl;
        cholmod_l_finish(cc);
        it=0;
        flag=0;
        return;
    }
    else{
        std::cout << "explicit solution residual is: "<< r <<std::endl;
        std::cout << "refine accuracy by CG iterations.... "<< std::endl;
    }

    // ------------------- preconditioned LSQR iteration --------------- //

    A_cholmod = to_cholmod_sparse(&A, cc);
    Areduced = cholmod_l_submatrix(A_cholmod, NULL, -1, 
    E, rank , 1,1,cc ); // transform A by E (Truncate A, to get rid of rank-deficient part)

    t1 = read_wtime();
    CSC_Mat_d Areduced_CSC(Areduced->nrow, Areduced->ncol, Areduced->nzmax, 
        (long*)Areduced->i, (long*)Areduced->p, (double*)Areduced->x);
    Vec_d y(rank);
    lsqr( Areduced_CSC, r_s, it_tol, max_it, y, flag, it, R_11, 
        success, perturb, debug); // solving min_x || AP_1 R_11^{-1} y - (b - AP_1 c_cholmod) ||_2
                                  // Using the explicit sketching solution as the starting point
    t2 = read_wtime();
    std::cout << "lsqr time is: "<< t2-t1 << std::endl;

    // ----- get the solution for Areduced*y = b -------- //
    cs_usolve_cholmod_structure(R_11, y.data());  // get w = R11^{-1} * y

    // add the initial guess c_cholmod
    for (long i=0; i<rank; i++){
        y(i) = y(i) + ((double*)c_cholmod->x)[i];
    }

    // ------- finally multiply by a permutation matrix to get original solution---//

    t1 = read_wtime();
    long permMatrix_p[Areduced->ncol+1];
    double allOnes[rank];
    for (long i = 0; i<(long)(Areduced->ncol+1); i++){
        permMatrix_p[i] = i;
    }  
    for (long i =0; i< rank; i++){
        allOnes[i] = 1;
    }
    CSC_Mat_d permMatrix(A.n, Areduced->ncol, Areduced->ncol,
        (long*)E, permMatrix_p, allOnes);
    permMatrix.mv('n', 1, y, 0, x);
    t2 = read_wtime();
    std::cout << "final post processing time: "<< t2- t1<< std::endl;

    // clean up
    cholmod_l_finish(cc);
}



