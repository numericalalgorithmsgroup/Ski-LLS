#include "lsex_all_in_c.h"
#include <stdlib.h>
#include "Config.hpp"
#include <cmath>
#include "RandSolver.hpp"
#include "LinAlg.hpp"
#include "IterSolver.hpp"
#include "SpMat.hpp"
#include "cs.h"         
#include "cholmod.h"
#include "SuiteSparseQR.hpp"
#include "cholmod_core.h"
#include "fhsl_mi35.h"

void ls_sparse_hsl(cs_dl & A, Vec_d& b, Vec_d& x, 
    int& flag, long& it,double it_tol, long max_it, 
    fint ordering, fint lsize, fint rsize, fint scaling, double shifting,
    int debug)
{
#ifdef HAVE_HSL
    assert( A.m == b.n() ); // input: A,b, parameter rcond, return by reference, x and rank
    assert( A.n == x.n() );

    if (A.m<A.n) {
        throw std::runtime_error( "Matrix is not overdetermined." );
    }

    // Do incomplete Cholesky, 
    long * col = A.p;
    long * row = A.i;
    double *val = A.x;
    fint nnz = A.nzmax;
    fint ifail = 0;
    void * pkeep = NULL;
    double rcontrol[20];
    fint icontrol[20];
    rcontrol[0] = shifting;
    icontrol[0] = ordering;
    icontrol[1] = scaling;
    // fint no_trans = 0;
    fint trans = 1; //various set-up

    long* col_oneBased = (long*) malloc( (A.n+1) * sizeof(long));
    long* row_oneBased = (long*) malloc( (nnz) * sizeof(long));

    for (fint i=0; i<A.n+1; i++){
        col_oneBased[i] = col[i]+1;
    }
    for (fint i =0; i<nnz; i++){
        row_oneBased[i] = row[i]+1;
    } //need one based indices

    hsl_mi35_factorize(&A.m, &A.n, &nnz, 
                (fint*)col_oneBased /*[n+1]*/, (fint*)row_oneBased /*[nnz]*/,
                (double*)val /*nnz*/, 
                &lsize , &rsize, & pkeep, rcontrol, 
                icontrol, &ifail);
    std::cout << "ifail for factorize is: "<<ifail << std::endl;

    // ---------------------- Do LSQR -------------------------- // 
    CSC_Mat_d A_CSC(A.m, A.n, A.nzmax, 
        (long*)A.i, (long*)A.p, (double*)A.x);
    Vec_d y(A.n);
    lsqr_ic( A_CSC, b, it_tol, max_it, y, flag, it, pkeep, debug);

    hsl_mi35_solve(&trans, &A.n, &pkeep, 
      y.data() /*[n]*/, x.data() /*[n]*/, &ifail);
    hsl_mi35_finalise(&pkeep, &ifail);
#else
    throw std::runtime_error( "HSL not available, compile with -DHAVE_HSL." );
#endif
}
