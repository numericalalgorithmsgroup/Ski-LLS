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
#include "cholmod_core.h"
#include <sys/timeb.h>

double read_wtime()
{
  struct timeb T;
  (void) ftime( &T );
  return ((double)T.time) + (1e-3) * ((double)T.millitm);
}

// Convert a sparse matrix A in cs_dl format to cholmod_sparse format
cholmod_sparse* to_cholmod_sparse(cs_dl* A, 
    cholmod_common *cc) // cc has to be started with!
{
    cholmod_sparse* A_cholmod;
    A_cholmod = cholmod_l_allocate_sparse(
    A->m, A->n, A->nzmax, false, true, 0, CHOLMOD_REAL, cc);

    A_cholmod->p = A->p;
    A_cholmod->i = A->i;
    A_cholmod->x = A->x;
    return A_cholmod;
}

// Convert a sparse matrix A in cholmod_sparse format to cs_dl format
cs_dl* to_cs_dl(cholmod_sparse* A)
{
    cs_dl* A_cs;
    A_cs = cs_dl_spalloc(A->nrow, A->ncol, A->nzmax, 1,0 );

    A_cs->p = (long*)A->p;
    A_cs->i = (long*)A->i;
    A_cs->x = (double*)A->x;
    return A_cs;
}

// Convert a sparse matrix A in cs_dl format to dense format
Mat_d toDense(cs_dl *S){
    double* data_SM = (double *)calloc((S->m)*(S->n), sizeof(double));

        for (long j = 0; j<(S->n); j++)
    {
        for (long k = (S->p)[j]; k<(S->p)[j+1]; k++) 
        {
            data_SM[j* (S->m) + (S->i)[k]] = (S->x)[k];
        }
    }

    Mat_d SDense(S->m, S->n, data_SM);
    return SDense;
}

// Generate an array 1:n
void gen_0_to_n(long n, long *returned_array)
{
    for (long i = 0; i<n; i++){
        returned_array[i] = i;
    }
}

// Randomly shuffle an array so that the first k indices are sampled from 1:n without replacement
void select_k_from_n(
    long n /* length of array */, 
    long k /*number of indices to be shuffled*/, 
    long* array_to_be_shuffled)
{
    // input check
    assert(n >= k);
    for (long i=0; i<k; i++){
        long j = ( rand() % (n-i) ) +i;
        long tmp = array_to_be_shuffled[j];
        array_to_be_shuffled[j] = array_to_be_shuffled[i];
        array_to_be_shuffled[i] = tmp;
    }
}

// Generate a hashing matrix, return in the compressed column sparse format
void gen_hashing_matrix(
    long m /* number of rows */, 
    long n /* number of columns */, 
    long nnz_per_column /* number of non-zeros per column */, 
    // output
    long *row_indices , long *col_array, double *values)
{
    long one_to_m[m];
    gen_0_to_n(m, one_to_m);
    for (int i=0; i<n; i++){
        select_k_from_n(m, nnz_per_column, one_to_m);
        for (int j=0; j<nnz_per_column; j++){
            row_indices[nnz_per_column*i + j] = one_to_m[j];
        }
    } // the above loop generate the row indices

    for (long i=0; i<n+1; i++){
        col_array[i] = nnz_per_column*i;
    } // generate the column array

    double scaling = sqrt(nnz_per_column);
    for (long i=0; i<n; i++){
        for (long j =0; j<nnz_per_column; j++){
            values[i*nnz_per_column +j] = 1/scaling* ((rand()%2)*2-1);
        }
    }
}

// x -> U^{-1} x, the matrix U assumed to be upper triangular
// Overwrite x
// Optionally apply perturbation for pivots if flag = 0
// Adapted from Tim Davis's CSparse
void cs_usolve_cholmod_structure (cholmod_sparse* U ,double *x, 
    int flag, double perturb)
{
    long n = U->ncol;
    long* Up = (long*)U->p;
    long* Ui = (long*)U->i;
    double* Ux = (double*)U->x;
    if (flag==1){
        for (long j = n-1 ; j >= 0 ; j--)
            {
                x [j] /= Ux [Up [j+1]-1] ;
                for (long p = Up [j] ; p < Up [j+1]-1 ; p++)
                {
                    x [Ui [p]] -= Ux [p] * x [j] ;
                }
            }
    } else{
        for (long j = n-1 ; j >= 0 ; j--)
            {
                x [j] /= Ux [Up [j+1]-1] + perturb ;
                for (long p = Up [j] ; p < Up [j+1]-1 ; p++)
                {
                    x [Ui [p]] -= Ux [p] * x [j] ;
                }
            }
    }

}

// x -> U^T^{-1} x, the matrix U assumed to be upper triangular
// Overwrite x
// Optionally apply perturbation for pivots if flag = 0
// Adapted from Tim Davis's CSparse
void cs_utsolve_cholmod_structure (cholmod_sparse* U, double *x, 
    int success, double perturb)
{
    long n = U->ncol;
    long* Up = (long*)U->p;
    long* Ui = (long*)U->i;
    double* Ux = (double*)U->x;
    if (success==1){
        for (long j = 0 ; j < n ; j++)
            {
                for (long p = Up [j] ; p < Up [j+1]-1 ; p++)
                {
                    x [j] -= CS_CONJ (Ux [p]) * x [Ui [p]] ;
                }
                x [j] /= CS_CONJ (Ux [Up [j+1]-1]) ;
            }
    } else{
        for (long j = 0 ; j < n ; j++)
            {
                for (long p = Up [j] ; p < Up [j+1]-1 ; p++)
                {
                    x [j] -= CS_CONJ (Ux [p]) * x [Ui [p]] ;
                }
                x [j] /= CS_CONJ (Ux [Up [j+1]-1]) + perturb ;
            }
    }
}