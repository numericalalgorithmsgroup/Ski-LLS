/**
 * @file   LinAlg.hpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Thu Jul 28 14:43:26 2011
 * 
 * @brief  linear algebra
 * 
 * 
 */

#ifndef _LINALG_HPP
#define _LINALG_HPP

#include "Config.hpp"
#include "Vec.hpp"
#include "LinOp.hpp"
#include "flapack.h"

template <typename T>
void svd( Mat<T> A, Vec<T> sgm, Mat<T> U_VT );

template <typename T>
void chol( Mat<T> A );

template <typename T>
void chol_sol( const Mat<T> L, Mat<T> BX );

template <typename T>
void chol_sol( const Mat<T> L, Vec<T> bx );

template <typename T>
void lstsq( Mat<T> A, Mat<T>& BX, Vec<T> s, const T rcond, long& rank );

#endif  // _LINALG_HPP
