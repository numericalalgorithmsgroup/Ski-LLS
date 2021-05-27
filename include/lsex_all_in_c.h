// @author Zhen Shao
#ifndef LSRN_CPP_0_1_LSEX_ALL_IN_C_H
#define LSRN_CPP_0_1_LSEX_ALL_IN_C_H

#include "Config.hpp"
#include "Mat.hpp"
#include "cs.h"
#include "SuiteSparseQR.hpp"
#include "cholmod.h"
#include "fhsl_mi35.h"
#include "bench_config.hpp"

// helpers
cholmod_sparse* to_cholmod_sparse(cs_dl* A, cholmod_common *cc);
Mat_d toDense(cs_dl *S);
void gen_0_to_n(long n, long* returned_array); 
void select_k_from_n(long n, long k, long* base_array);
void gen_hashing_matrix(long m, long n, long nnz_per_column, 
	long *row_indices, long *col_array,double *values);
cs_dl* to_cs_dl(cholmod_sparse* A);
double read_wtime();
void cs_usolve_cholmod_structure(cholmod_sparse* R ,double *x, 
	int success=1, double perturb=PERTURB);
void cs_utsolve_cholmod_structure(cholmod_sparse* R ,double *x, 
	int success=1, double perturb=PERTURB);

void cholmod_cond_est(cholmod_sparse* A, cholmod_common* cc, double &rcond);

// sketching
Mat_d sketch(Mat_d& A, const long k, const double gamma);
Mat_d sketch(Mat_d& A, const long k, const double gamma, const Mat_d& b, Mat_d& b_hashed);
cs_dl* sparse_sketch(cs_dl &A, const long k, const double gamma);
cs_dl* sparse_sketch(cs_dl &A, const long k, const double gamma, 
    double* b, double* &b_sketched);


// solvers


void ls_sparse_hsl(cs_dl & A, Vec_d& b, Vec_d& x, int& flag, long& it, 
    double it_tol=IT_TOL, long max_it=MAX_IT, fint ordering=ORDERING,
    fint lsize = 20, fint rsize = 20, fint scaling = SCALING, 
    double shifting = SHIFTING, int debug=0);

void column_pivoted_qr(Mat_d& A, long* E);
void column_pivoted_qr(Mat_d& A, long* E, double* tau);
Mat_d reducedQR(Mat_d &A, double rcond, long* &E);


#endif //LSRN_CPP_0_1_LSEX_ALL_IN_C_H
