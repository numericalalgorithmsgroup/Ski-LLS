// @author Xiangrui Meng
// Extended by Zhen Shao

#ifndef _RANDSOLVER_HPP
#define _RANDSOLVER_HPP

#include "Config.hpp"
#include "Vec.hpp"
#include "LinOp.hpp"

/** 
 * Random normal projection. Project A's columns/rows to an s-dimensional space.
 * Return projected matrix (transposed) As_T
 */
Mat_d rnpj( LinOp<double>& A, ptrdiff_t s );

/** 
 * Preconditioning using random normal projection
 * Return right pre-conditioner N
 */
Mat_d rnpre( LinOp<double>& A, long s, double rcond = 1e-12 );

/*
	Preconditioning using given projection matrix S
	Return pre-conditioner N
*/
Mat_d rnpre_given_sketch( LinOp<double>& A, long s, Mat_d As_T,double rcond = 1e-12);

/*
	Linear least square solver by dense Gaussian projection and SVD
	Solution x, flag =0 means success. 
	Record rank of matrix and iteration taken by lsqr
*/
void lsrn( LinOp<double>& A, Vec_d& b, double rcond, 
	Vec_d& x, long& rank, int& flag, long& it,
	 double gamma, double it_tol=1e-12, long max_it=1000);

/*
	Linear least square solver given preconditioner, use SVD
	Solution x, flag =0 means success. 
	Record rank of matrix and iteration taken by lsqr
*/
void lsex( LinOp<double>& A, Vec_d& b, double rcond, 
	Vec_d& x, long& rank, int& flag, long& it, Mat_d& S );

#endif  // _RANDSOLVER_HPP
