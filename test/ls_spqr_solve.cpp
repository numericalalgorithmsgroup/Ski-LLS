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

/* Command-line interface to use the solver SPQR (see Tim Davis's paper)
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

    cholmod_common Common, *cc;
    cc = &Common;
    cholmod_l_start(cc);
    cholmod_sparse *Araw, *A;

    cholmod_dense *X, *B, *Residual = NULL ;
    double rnorm, one [2] = {1,0}, minusone [2] = {-1,0} ;


    int mtype;
    FILE *fp;
    fp = fopen(argv[1],"r");
    Araw = (cholmod_sparse *) cholmod_l_read_matrix(fp, 1, &mtype, cc) ;

    fclose(fp);

    // if A is fat, transpose it
    if ( (Araw->nrow) < (Araw->ncol) ){
        A = cholmod_l_transpose(Araw, 1, cc);
    }
    else{
        A = Araw;
    }


    std::cout << "number of rows is: " << A->nrow << std::endl;
    std::cout << "number of columns is: " << A->ncol << std::endl;
    std::cout << "number of non-zeros is: " << A->nzmax << std::endl;   
    B = cholmod_l_ones (A->nrow, 1, A->xtype, cc) ;
    std::cout << "Solver Name: Ls_spqr_solve" << std::endl;


    // storage

    double t_start;
    double t_finish;
    double residual;   

    // solve, sparse_spqr
    t_start = read_wtime();
    X = SuiteSparseQR <double> (A, B, cc) ;
    t_finish = read_wtime();

    std::cout << "Time taken by ls_sparse_spqr is: "<< (t_finish-t_start) << std::endl;

    Residual = cholmod_l_copy_dense (B, cc) ;
    cholmod_l_sdmult (A, 0, minusone, one, X, Residual, cc) ;
    rnorm = cholmod_l_norm_dense (Residual, 2, cc) ;

    std::cout << "Residual of ls_sparse_spqr is: "<<  rnorm << std::endl;


    std::cout << "==============FINISH==================" << std::endl << std::endl;

    cholmod_l_free_dense (&Residual, cc) ;
    cholmod_l_free_sparse (&A, cc) ;
    cholmod_l_free_dense (&X, cc) ;
    cholmod_l_free_dense (&B, cc) ;
    cholmod_l_finish (cc) ;
}