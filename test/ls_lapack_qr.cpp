#include <iostream>
#include "Config.hpp"
#include "SpMat.hpp"
#include "Random.hpp"
#include <math.h>
#include "RandSolver.hpp"
#include <time.h>
#include "cs.h"
#include <stdio.h> 
#include "SuiteSparseQR.hpp"
#include "cholmod.h"
#include "lsex_all_in_c.h"
#include "IterSolver.hpp"
#include <sys/timeb.h>
#include "bench_config.hpp"
#include "blendenpik.hpp"

/* Command-line interface to use the solver ls_dense_lapack_qr (see ls_dense_blendenpik)
    Solving linear least square with a given dense matrix A and 
        the right-hand side a vector of all ones
     argv[1]: Matrix path, points to a dense matrix in matrix market format
*/
int main(int argc, char **argv){

    int debug=0;
    std::cout.precision(PRINTING_PRECISION);
    std::cout << std::fixed;

    std::cout << std::endl << "==============BEGIN==================" << std::endl;
    std::cout << "Matrix path: " << argv[1] << std::endl;
    std::cout << "Reading input.........." << std::endl;

	// input
    cholmod_common Common, *cc;
    cc = &Common;
    cholmod_l_start(cc);

    cholmod_dense* Acholmod;
    int mtype;
    FILE *fp;
    fp = fopen(argv[1],"r");
    Acholmod = (cholmod_dense *) cholmod_l_read_matrix(fp, 1, &mtype, cc) ;
    fclose(fp);
    std::cout << "number of rows is: " << Acholmod->nrow << std::endl;
    std::cout << "number of columns is: " << Acholmod->ncol << std::endl;
    std::cout << "number of non-zeros is: " << Acholmod->nzmax << std::endl;

    Mat_d A(Acholmod->nrow, Acholmod->ncol, (double*) Acholmod->x);

    // construct b
    double *b_data = (double*)malloc((A.m())*sizeof(double));
    for (long i=0; i<A.m(); i++){
        b_data[i] = 1;
    }
    Vec_d b(A.m(), b_data);

    // parameters
    long k = NNZ_PER_COLUMN;
    double gamma = OVER_SAMPLING_RATIO;
    long max_it = MAX_IT;
    double it_tol = IT_TOL;
    double rcond = 1e-12;
    int wisdom = WISDOM;


    std::cout << "Solver Name: Ls_blendenpik" << std::endl;
    std::cout << "Configuration: " << std::endl;
    std::cout << "NNZ_PER_COLUMN is: "<< k << std::endl;
    std::cout << "OVER_SAMPLING_RATIO is: "<< gamma << std::endl;
    std::cout << "MAX_IT is: "<< max_it << std::endl;
    std::cout << "IT_TOL is: "<< IT_TOL << std::endl;
    std::cout << "rcond for SVD is: "<< rcond << std::endl;
    std::cout << "WISDOM is: " << wisdom << std::endl;


    // storage

    long it;
    int flag;
    long rank;
    double t_start;
    double t_finish;
    double residual;   
    Vec_d x(A.n());

    // solve
    t_start = read_wtime();
    ls_dense_lapack_qr(A, b, x);
    t_finish = read_wtime();

    std::cout << "Time taken by ls_lapack is: "<< (t_finish-t_start) << std::endl;

    // std::cout << A << std::endl;
    // std::cout << x << std::endl;
    // std::cout << b << std::endl;
    A.mv('n', 1, x, -1, b);

    std::cout << "Residual of ls_lapack is: "<<  nrm2(b) << std::endl;

    std::cout << "==============FINISH==================" << std::endl << std::endl;
   
    cholmod_l_finish(cc);

}