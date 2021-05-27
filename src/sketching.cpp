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

// Sketch a dense matrix A using a hashing matrix
Mat_d sketch(Mat_d& A, 
    const long k, /* nnz per column in the hashing matrix */
     const double gamma /* over-sampling ratio */
    )
{
    long n = A.m();
    long d = A.n();
    long m = ceil(d*gamma);
    long nnz = n*k;
    long *row = (long*)malloc(nnz*sizeof(long));
    long *col = (long*)malloc((n+1)*sizeof(long));
    double *vals = (double*)malloc(nnz*sizeof(double));
    gen_hashing_matrix(m, n, (long)k, (long*)row, (long*)col, vals);
    CSC_Mat_d S( m, n, (long)nnz, (long *) row, (long *) col, (double *)vals );

    SpMat_Op<double> S_op(S);
    Mat_d As(m,d);
    S_op.mm('n','n',1/sqrt(k), A, 0.0, As);
    return As;
}

// Sketch a dense matrix A using a hashing matrix, also sketch b
Mat_d sketch(Mat_d& A, 
    const long k, /* nnz per column in the hashing matrix */
    const double gamma, /* over-sampling ratio */
    const Mat_d &b, /* changes on output */
    Mat_d &b_hashed
    )
{
    long n = A.m();
    long d = A.n();
    long m = ceil(d*gamma);
    long nnz = n*k;
    long *row = (long*)malloc(nnz*sizeof(long));
    long *col = (long*)malloc((n+1)*sizeof(long));
    double *vals = (double*)malloc(nnz*sizeof(double));
    gen_hashing_matrix(m, n, (long)k, (long*)row, (long*)col, vals);
    CSC_Mat_d S( m, n, (long)nnz, (long *) row, (long *) col, (double *)vals );

    b_hashed = b.copy();
    SpMat_Op<double> S_op(S);
    Mat_d As(m,d);
    S_op.mm('n','n',1/sqrt(k), A, 0.0, As);
    S_op.mm('n','n',1/sqrt(k), b_hashed, 0.0, b_hashed);
    return As;
}

// Sketch a sparse matrix A using a hashing matrix, return a sparse matrix
cs_dl* sparse_sketch(cs_dl &A, 
    const long k, /* nnz per column in the hashing matrix */
    const double gamma /* over-sampling ratio */
    )    
{ 
    long n = A.m;
    long d = A.n;
    long m = ceil(d*gamma);
    long nnz = n*k;
    long *row = (long*)malloc(nnz*sizeof(long));
    long *col = (long*)malloc((n+1)*sizeof(long));
    double *vals = (double*)malloc(nnz*sizeof(double));
    gen_hashing_matrix(m, n, (long)k, (long*)row, (long*)col, vals);

    cs_dl *S;
    S = cs_dl_spalloc(m, n, nnz, true, false);
    S->p = col;
    S->i = row;
    S->x = vals;

    return cs_dl_multiply(S, &A);
}


// Sketch the augmented matrix [A,b], A is a sparse matrix, b is a dense vector
cs_dl* sparse_sketch(cs_dl &A, const long k, const double gamma, 
    double* b, double* &b_sketched)    
{ // return spase matrix for sparse factorization + no need to do transpose
    long n = A.m;
    long d = A.n;
    long m = ceil(d*gamma);
    long nnz = n*k;
    long *row = (long*)malloc(nnz*sizeof(long));
    long *col = (long*)malloc((n+1)*sizeof(long));
    double *vals = (double*)malloc(nnz*sizeof(double));
    gen_hashing_matrix(m, n, (long)k, (long*)row, (long*)col, vals);


    cs_dl *S;
    S = cs_dl_spalloc(m, n, nnz, true, false);
    S->p = col;
    S->i = row;
    S->x = vals;

    b_sketched = (double*) calloc(m,sizeof(double));
    cs_dl_gaxpy(S, b, b_sketched);
    return cs_dl_multiply(S, &A);
}