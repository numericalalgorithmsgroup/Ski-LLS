/**
 * @file   Random.hpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Tue Oct  4 22:58:43 2011
 * 
 * @brief  
 * 
 * 
 */

#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include <cstddef>

#include "Vec.hpp"
#include "LinOp.hpp"

/** 
 * Multi-threaded pseudorandom number generator for the standard normal distribution.
 * 
 * @param n number of samples
 * @param a array to store random numbers
 */
void randn( ptrdiff_t n, double *a );

void randn( Vec_d v );
void randn( Mat_d A );

#endif  // _RANDOM_HPP
