//Copyright (C) 2009-2016, Haim Avron and Sivan Toledo.
// Extended by Zhen Shao

#ifndef BLENDENPIK_H
#define BLENDENPIK_H

#include "Config.hpp"
#include "Mat.hpp"
#include "cs.h"
#include "SuiteSparseQR.hpp"
#include "cholmod.h"
#include "bench_config.hpp"

#define FFTW_R2R_DCT   1
#define FFTW_R2R_DHT   2


// For original blendenpik

void gen_scaled_rand_signs(long m, double* returned_array /* allocated to be size m */, double scale=1.0);

Mat_d fast_unitary_transform(const Mat_d &A /* matrix m_0 by d */, 
	double* D /* array of size m_0 */, int wisdom=0);

Mat_d random_sample(const Mat_d &A, long num_of_rows);

void dense_qr(Mat_d& A);

Mat_d rand_subsampled_dct(const Mat_d &A, double gamma, int wisdom=0);

void ls_dense_blendenpik( Mat_d & A, Vec_d& b, Vec_d& x,
 long& rank, int& flag, long& it, double gamma,
 double it_tol, long max_it, int debug, int wisdom);

// For hashing Blendenpik (new)

Mat_d rand_hashing_dct(const Mat_d&A, double gamma, long k, int wisdom=0);
Mat_d rand_hashing_dct(const Mat_d&A, const Vec_d&b, 
		Vec_d &b_out, double gamma, long k, int wisdom=0);


void ls_dense_lapack_qr( Mat_d & A, Vec_d& b, Vec_d& x);

#endif 
