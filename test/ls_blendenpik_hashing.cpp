#include <iostream>
#include <math.h>
#include <sys/timeb.h>

/* Here we need this for conversion between cholmod to cs_dl, cs_dl to Mat_d and read_wtime functions.
    It is not essential to our solver.
*/
#include "lsex_all_in_c.h"

// Macro for configurations
#include "bench_config.hpp"

#include "ski-lls.h"

/* Command-line interface to use the solver Blendenpik_hashing  (See ls_dense_hashing_blendenpik)
    Solving linear least square with a given dense matrix A and 
        the right-hand side a vector of all ones
     argv[1]: Matrix path, points to a dense matrix in matrix market format
     argv[2]: (Optional) Over-sampling ratio used in sketching
     argv[3]: (Optional) nnz per column in the hashing matrix
*/
int main(int argc, char **argv){
    double t1;
    double t2;
    int debug=0;
    std::cout.precision(PRINTING_PRECISION);
    std::cout << std::fixed;
    std::cout << std::endl << "==============BEGIN==================" << std::endl;
    std::cout << "Matrix path: " << argv[1] << std::endl;
    t1 = read_wtime();
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
    t2 = read_wtime();
    std::cout << "Reading input time is: "<< t2-t1 << std::endl;
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
    long k = 1;
    double gamma = 1.7;
    long max_it = MAX_IT;
    double it_tol = IT_TOL;
    double rcond = 1e-12;
    int wisdom = WISDOM;
    double abs_tol = ABS_TOL;

    /* If there is more input, overwrite k and gamma*/
    if (argc>=3){
        gamma = atof(argv[2]);
    }
    if (argc>=4){
        k = atoi(argv[3]);
    }

    std::cout << "Solver Name: Ls_blendenpik_hashing" << std::endl;
    std::cout << "Configuration: " << std::endl;
    std::cout << "NNZ_PER_COLUMN is: "<< k << std::endl;
    std::cout << "OVER_SAMPLING_RATIO is: "<< gamma << std::endl;
    std::cout << "MAX_IT is: "<< max_it << std::endl;
    std::cout << "IT_TOL is: "<< IT_TOL << std::endl;
    std::cout << "rcond for SVD is: "<< rcond << std::endl;
    std::cout << "absolute tolerance is: "<< abs_tol << std::endl;    
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
    ls_dense_hashing_blendenpik(A, b, x, rank, flag, it, gamma, k, abs_tol, rcond, it_tol, max_it, debug, wisdom);
    t_finish = read_wtime();

    std::cout << "Iteration flag is: "<<  flag << std::endl;

    std::cout << "Time taken by ls_blendenpik_hashing is: "<< (t_finish-t_start) << std::endl;

    if (debug==1){
        std::cout << "A is: "<< A << std::endl;
        std::cout << "b is: "<< b << std::endl;
        std::cout << "solution is: "<< x << std::endl;
    }

    A.mv('n', 1, x, -1, b); //Replace with A after the bug with CPQR is fixed

    std::cout << "Residual of ls_blendenpik_hashing is: "<<  nrm2(b) << std::endl;

    std::cout << "Iteration of ls_blendenpik_hashing is: " << it << std::endl;

    if (flag==0){
        std::cout << "Status: LSQR Converged, a minimal residual solution found" << std::endl;
    }
    else{
        std::cout << "Status: Maximum iteration reached, a suboptimal solution found" << std::endl;
    }

    std::cout << "==============FINISH==================" << std::endl << std::endl;
   
    cholmod_l_finish(cc);

}