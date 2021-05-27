/**
 * @file   Utils.cpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Mon Oct 10 21:14:25 2011
 * 
 * @brief  
 * 
 * 
 */

#include <sys/time.h>
#include <unistd.h>

#include "Config.hpp"
#include "Utils.hpp"

long get_num_cores()
{
  long n_threads = sysconf( _SC_NPROCESSORS_ONLN );
  if( n_threads > MAX_NUM_THREADS )
    n_threads = MAX_NUM_THREADS;
  return n_threads;
}

double get_time()
{
  struct timeval tv;
  gettimeofday( &tv, NULL );
  return (double) tv.tv_sec + 1e-6*tv.tv_usec;
}
