#include <iostream>
#include "ski-lls.h"

int main(int argc, char **argv)
{
	// Create a sparse matrix in compressed column format
    long m = 4;
    long n = 3;
    long nnz = 5;

    long col[4] = {0,  1,  3,  5};
    long row[5] = {0,  0,  2,  1,  3}; 
    double val[5] = {2.0,  3.0,  1.0,  4.0,  5.0};

    cs_dl* A;
    A = cs_dl_spalloc(m, n, nnz, 1,0 );

    A->p = (long*)col;
    A->i = (long*)row;
    A->x = (double*)val;


    double data_b[4] = {0.0649338,0.845946,-0.0164085,0.247119};
    Vec_d b(m, data_b);

    // parameters
    long k = NNZ_PER_COLUMN;
    double gamma = OVER_SAMPLING_RATIO;
    long max_it = MAX_IT;
    double it_tol = IT_TOL;
    double abs_tol = ABS_TOL;
    int ordering = SPQR_ORDERING;
    int debug = 0;

    // storage

    long it;
    int flag;
    long rank;
    double t_start;
    double t_finish;
    double residual;   
    Vec_d x(A->n);

    // solve
    ls_sparse_spqr(*A, b, x, 
    rank, flag, it, gamma, k, abs_tol, ordering, it_tol, max_it, RCOND_THRESHOLD, PERTURB, debug);

    // compute residual
    CSC_Mat_d A_CSC(A->m, A->n, A->nzmax, 
        (long*)A->i, (long*)A->p, (double*)A->x);
    A_CSC.mv('n', 1, x, -1, b);

    std::cout << "A is: " << std::endl;
    cs_dl_print(A,0);
    std::cout << "b is: "<< b << std::endl;
    std::cout << "solution is: "<< x << std::endl;
    // Should get     x = (0.0571 -0.0164 0.1127)

    std::cout << "Residual of ls_sparse_spqr is: "<<  nrm2(b) << std::endl;

}