// Author Zhen Shao
#include <random>
#include "lsex_all_in_c.h"
#include "blendenpik.hpp"

// Generate an array of random signs, multiplied by a scale
void gen_scaled_rand_signs(
    long m, /* Array length */
    double* returned_array, 
    double scale)
{
	srand(time(NULL));
	for (long i=0; i<m; i++){
		returned_array[i] = scale*((rand()%2)*2-1);
	}
}

// Return a randomly (row) subsampled matrix
Mat_d random_sample(
    const Mat_d &A, /* Original matrix */
    long num_of_rows /* Number of rows to be sampled */
     )
{

    srand(time(NULL));
	long m = A.m();
	long index_array[m];
    Mat_d A_sub(num_of_rows, A.n());

    gen_0_to_n(m, index_array);
    // Shuffle the array so that the first num_of_rows indices point to rows sampled
    select_k_from_n(m, num_of_rows, index_array); 
    for(long j=0; j<A.n(); j++){
    	for (long i=0; i< num_of_rows; i++){
    		A_sub(i,j) = A(index_array[i],j);
    	}
    }
    return A_sub;
}

// Perform the subsampled randomised discrete Cosin transform
Mat_d rand_subsampled_dct(
    const Mat_d &A, /* Input */
    double gamma, /* over-sampling ratio (number of sampled rows/number of columns of A) */
    int wisdom /* fftw wisdom */
    )
{
    srand(time(NULL));
    long sketched_dim = ceil(gamma* A.n());
    double* random_signs = (double*)malloc(sizeof(double) * A.m());

    gen_scaled_rand_signs(A.m(), random_signs, 1.0);
    Mat_d A_transformed = fast_unitary_transform(A, random_signs, wisdom);
    Mat_d A_sub = random_sample(A_transformed, sketched_dim);

    free(random_signs);
    free(A_transformed.data());
    return A_sub;
}

// Perform the hashed randomised discrete cosin transform
Mat_d rand_hashing_dct(
    const Mat_d&A,  /* input */
    double gamma, /* over-sampling ratio (number of hashed rows/number of columns of A) */
    long k,  /* nnz per column of the hashing matrix */
    int wisdom /* fftw wisdom */
    )
{
    srand(time(NULL));
    double* random_signs = (double*)malloc(sizeof(double) * A.m());
    gen_scaled_rand_signs(A.m(), random_signs, 1.0);
    Mat_d A_transformed = fast_unitary_transform(A, random_signs, wisdom);

    Mat_d A_hashed = sketch(A_transformed, k, gamma);
    free(random_signs);
    free(A_transformed.data());
    return A_hashed;
}


// Perform the hashed randomised discrete cosin transform
// Transform b as well
Mat_d rand_hashing_dct(
    const Mat_d&A,  /* input */
    const Vec_d &b, /* on input, a vector b, on output, sketched vector b*/
    Vec_d &b_out, /* output, sketched b*, data allocated at least gamma*A.n()*/
    double gamma, /* over-sampling ratio (number of hashed rows/number of columns of A) */
    long k,  /* nnz per column of the hashing matrix */
    int wisdom /* fftw wisdom */
    )
{
    srand(time(NULL));
    double* random_signs = (double*)malloc(sizeof(double) * A.m());
    gen_scaled_rand_signs(A.m(), random_signs, 1.0);

    Mat_d b_mat(A.m(), 1, b.data());

    Mat_d A_transformed = fast_unitary_transform(A, random_signs, wisdom);
    Mat_d b_transformed = fast_unitary_transform(b_mat, random_signs, wisdom);

    Mat_d b_hashed;
    Mat_d A_hashed = sketch(A_transformed, k, gamma, b_transformed, b_hashed);

    for (long i=0; i< ceil(A.n()*gamma); i++){
        b_out.data()[i] = b_hashed.data()[i];
    }
    free(random_signs);
    free(A_transformed.data());
    free(b_transformed.data());

    return A_hashed;
}