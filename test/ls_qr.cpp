#include <iostream>
#include <math.h>
#include <sys/timeb.h>

/* Here we need this for conversion between cholmod to cs_dl and read_wtime functions.
    It is not essential to our solver.
*/
#include "lsex_all_in_c.h"

// Macro for configurations
#include "bench_config.hpp"

#include "ski-lls.h"

/* Command-line interface to use the solver LS_QR (See ls_qr)
    Solving linear least square with a given sparse matrix A and 
        the right-hand side a vector of all ones, 
     argv[1]: Matrix path, points to a sparse matrix in matrix market format
     argv[2]: (Optional) Over-sampling ratio for sketching
     argv[3]: (Optional) nnz per column in the hashing matrix
     argv[4]: Absolute tolerance used in the solver
     argv[5]: Ordering used in sparse QR factorization
*/

int main(int argc, char **argv){
    std::cout.precision(PRINTING_PRECISION);
    std::cout << std::fixed;
    int debug=0;
    // double rcond_thres=1e-10;
    // double perturb=1e-10;   

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

    if(debug==1){
        cc->print=4;
        cholmod_l_print_sparse(Acholmod, "Acholmod_raw",cc);
    }

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
    long k = 2;
    double gamma = 1.4;
    long max_it = MAX_IT;
    double it_tol = IT_TOL;
    double abs_tol = ABS_TOL;
    double rcond = 1e-12;
    int ordering = SPQR_ORDERING;
    if(debug==1){
        max_it=20;
    }

    /* If there is more input, overwrite k and gamma*/
    if (argc>=3){
        gamma = atof(argv[2]);
    }
    if (argc>=4){
        k = atoi(argv[3]);
    }

    if (argc>=5){
        abs_tol = atof(argv[4]);
    }

    if (argc>=6){
        ordering=atoi(argv[5]);
    } 


    std::cout << "Solver Name: Ls_qr" << std::endl;
    std::cout << "Configuration: " << std::endl;
    std::cout << "NNZ_PER_COLUMN is: "<< k << std::endl;
    std::cout << "OVER_SAMPLING_RATIO is: "<< gamma << std::endl;
    std::cout << "MAX_IT is: "<< max_it << std::endl;
    std::cout << "IT_TOL is: "<< IT_TOL << std::endl;
    std::cout << "rcond for SVD is: "<< RCOND << std::endl;
    std::cout << "rcond_threshold is: "<< RCOND_THRESHOLD << std::endl;
    std::cout << "potential perturbation of pivot is: "<< PERTURB << std::endl;
    std::cout << "absolute tolerance is: "<< abs_tol << std::endl;
    std::cout << "ordering is: "<< ordering << std::endl;


    // storage

    long it;
    int flag;
    long rank;
    double t_start;
    double t_finish;
    double residual;   
    Vec_d x_ls_qr(A->n);

    // solve, sparse_spqr
    t_start = read_wtime();
    ls_sparse_spqr(*A, b, x_ls_qr, 
    rank, flag, it, gamma, k, abs_tol, ordering, it_tol, max_it, RCOND_THRESHOLD, PERTURB, debug);
    t_finish = read_wtime();

    std::cout << "Iteration flag is: "<<  flag << std::endl;

    std::cout << "Time taken by ls_sparse_spqr is: "<< (t_finish-t_start) << std::endl;

    CSC_Mat_d A_CSC(A->m, A->n, A->nzmax, 
        (long*)A->i, (long*)A->p, (double*)A->x);

    A_CSC.mv('n', 1, x_ls_qr, -1, b);

    std::cout << "Residual of ls_sparse_spqr is: "<<  nrm2(b) << std::endl;

    std::cout << "Iteration of ls_sparse_spqr is: " << it << std::endl;

        if (flag==0){
        std::cout << "Status: LSQR Converged, a minimal residual solution found" << std::endl;
    }
    else{
        std::cout << "Status: Maximum iteration reached, a suboptimal solution found" << std::endl;
    }

    std::cout << "==============FINISH==================" << std::endl << std::endl;

}
