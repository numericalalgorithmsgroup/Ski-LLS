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

/* Command-line interface to use the solver LSRN (See lsrn)
    Solving linear least square with a given sparse matrix A and 
        the right-hand side a vector of all ones, 
     argv[1]: Matrix path, points to a sparse matrix in matrix market format
*/

int main(int argc, char **argv){
    std::cout.precision(PRINTING_PRECISION);
    std::cout << std::fixed;
        std::cout << std::endl << "==============BEGIN==================" << std::endl;
    std::cout << "Matrix path: " << argv[1] << std::endl;
    std::cout << "Reading input.........." << std::endl;


	// input
	cs_dl *Araw, *A;
    cholmod_common Common, *cc;
    cc = &Common;
    cholmod_l_start(cc);

    cholmod_sparse* Acholmod;
    int mtype;
    FILE *fp;
    fp = fopen(argv[1],"r");
    Acholmod = (cholmod_sparse *) cholmod_l_read_matrix(fp, 1, &mtype, cc) ;

    fclose(fp);

    Araw = to_cs_dl(Acholmod);

    // if A is fat, transpose it
    if ( (Araw->m) < (Araw->n) ){
        A = cs_dl_transpose(Araw,true);
        cs_dl_free(Araw);
    }
    else{
        A = Araw;
    }


    std::cout << "number of rows is: " << A->m << std::endl;
    std::cout << "number of columns is: " << A->n << std::endl;
    std::cout << "number of non-zeros is: " << A->nzmax << std::endl;


    cholmod_l_finish(cc);

    
    double *b_data = (double*)malloc((A->m)*sizeof(double));
    for (long i=0; i<A->m; i++){
        b_data[i] = 1;
    }
    Vec_d b(A->m, b_data);

    // parameters
    long k = NNZ_PER_COLUMN;
    double gamma = OVER_SAMPLING_RATIO;
    long max_it = MAX_IT;
    double it_tol = IT_TOL;
    double rcond = 1e-12;
    /* If there is more input, overwrite k and gamma*/
    if (argc>=3){
        gamma = atof(argv[2]);
    }
    if (argc>=4){
        k = atoi(argv[3]);
    }

        std::cout << "Solver Name: Ls_lsrn" << std::endl;
    std::cout << "Configuration: " << std::endl;
    std::cout << "NNZ_PER_COLUMN is: "<< k << std::endl;
    std::cout << "OVER_SAMPLING_RATIO is: "<< gamma << std::endl;
    std::cout << "MAX_IT is: "<< max_it << std::endl;
    std::cout << "IT_TOL is: "<< IT_TOL << std::endl;
    std::cout << "rcond for SVD is: "<< rcond << std::endl;

    // storage

    long it;
    int flag;
    long rank;
    double t_start;
    double t_finish;
    double residual;   
    Vec_d x_ls_qr(A->n);

    CSC_Mat_d A_CSC(A->m, A->n, A->nzmax, 
        (long*)A->i, (long*)A->p, (double*)A->x);

    // solve, sparse_spqr
    t_start = read_wtime();
    lsrn(A_CSC, b, rcond, x_ls_qr, 
    rank, flag, it, gamma, it_tol, max_it);
    t_finish = read_wtime();

    std::cout << "Iteration flag is: "<<  flag << std::endl;

    std::cout << "Time taken by lsrn is: "<< (t_finish-t_start) << std::endl;



    A_CSC.mv('n', 1, x_ls_qr, -1, b);

    std::cout << "Residual of lsrn is: "<<  nrm2(b) << std::endl;

    std::cout << "Iteration of lsrn is: " << it << std::endl;

        if (flag==0){
        std::cout << "Status: LSQR Converged, a minimal residual solution found" << std::endl;
    }
    else{
        std::cout << "Status: Maximum iteration reached, a suboptimal solution found" << std::endl;
    }

    std::cout << "==============FINISH==================" << std::endl << std::endl;

}