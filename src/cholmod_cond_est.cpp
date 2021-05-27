// Author Zhen Shao
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <cmath>

#include "cholmod.h"
#include "flapack.h"

// Estimate the condition number of a sparse matrix A,
// The matrix is converted to dense format so that the LAPACK routine could be used
void cholmod_cond_est(cholmod_sparse* A, cholmod_common* cc, double &rcond){
	// input: A, cc
	// output: inverse condition number estimate rcond
	cholmod_dense* Adense;
	Adense = cholmod_l_sparse_to_dense(A,cc);

	// call LAPACK dtrcon
	long info;
	double * workspace1;
	long * workspace2;
	long n, ld;
	n = Adense->ncol;
	ld = Adense->d;
	char norm = 'I';
	char uplo = 'U';
	char diag = 'N';


	workspace1 = (double *)malloc(3 * n * sizeof(double));
  	workspace2 = (long *)malloc(n * sizeof(long));

	dtrcon( &norm, &uplo, &diag, &n, (double*)Adense->x, &ld, &rcond, 
             workspace1, workspace2, &info);
}