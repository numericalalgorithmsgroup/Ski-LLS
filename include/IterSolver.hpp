// @author Xiangrui Meng
// Extended by Zhen Shao


#ifndef _ITERSOLVER_HPP
#define _ITERSOLVER_HPP

#include "Vec.hpp"
#include "LinOp.hpp"

#include "cs.h"
#include "SuiteSparseQR.hpp"
#include "cholmod.h"
#include "lsex_all_in_c.h"

void lsqr( LinOp<double>& A, const Vec_d b,
           const double tol, const long maxit,
           Vec_d& x, int& flag, long& it, int debug=0);

void lsqr( LinOp<double>& A, const Vec_d b,
           const double tol, const long maxit,
           Vec_d& x, int& flag, long& it, cholmod_sparse* R_11, 
           int success=1, double perturb=1e-6, int debug=0); // LSQR with sparse preconditioner

void lsqr_dense_pre( LinOp<double>& A, const Vec_d b,
           const double tol, const long maxit,
           Vec_d& x, int& flag, long& it, Mat_d& R, int debug=0); // LSQR with dense preconditioner

void lsqr_ic( LinOp<double>& A, const Vec_d b,
           const double tol, const long maxit,
           Vec_d& x, int& flag, long& it, void* pkeep, int debug=0); 
           // LSQR with incomplete cholesky preconditioner

#endif  // _ITERSOLVER_HPP
