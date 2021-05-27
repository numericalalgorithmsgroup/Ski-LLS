#ifndef SKILLS_H
#define SKILLS_H

#include "Config.hpp"
#include "Mat.hpp"
#include "cs.h"
#include "SuiteSparseQR.hpp"
#include "cholmod.h"
#include "bench_config.hpp"

#define FFTW_R2R_DCT   1
#define FFTW_R2R_DHT   2

// From blendenpik.hpp
void ls_dense_hashing_blendenpik( Mat_d & A, Vec_d& b, Vec_d& x,
 long& rank, int& flag, long& it, double gamma, long k, double abs_tol,
 double rcond, double it_tol, long max_it, int debug, int wisdom);

void ls_dense_hashing_blendenpik_noCPQR( Mat_d & A, Vec_d& b, Vec_d& x,
 long& rank, int& flag, long& it, double gamma, long k,
 double it_tol, long max_it, int debug, int wisdom);

// From lsex_all_in_c.h
void ls_sparse_spqr(cs_dl& A, Vec_d& b, Vec_d& x, 
	long& rank, int& flag, long& it, double gamma, long k,
	double abs_tol = ABS_TOL, int ordering=SPQR_ORDERING, 
	double it_tol=IT_TOL, long max_it=MAX_IT, 
	double rcond_threshold=RCOND_THRESHOLD, double perturb=PERTURB, int debug=0);

#endif 
