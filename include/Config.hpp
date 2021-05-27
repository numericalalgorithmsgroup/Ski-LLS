#ifndef _CONFIG_HPP
#define _CONFIG_HPP

#include <cstddef>
#include "config_sk.h"

#define MAX_NUM_THREADS 1 // this needs to be changed if we want to use more cores
#define FFTW_QUANT 1000
#define FFTW_TIMES 150
#define FFTW_WISDOM_FILE "/fserver/zhens/fftw_wisdom-20-Jan.dat"
// #define FFTW_WISDOM_FILE "/Users/zhen/Dropbox/Codes/MATLAB_codes/NAG_mini_project/Blendenpik-master/dat/wisdom_file.dhcp-163-1-81-22.maths.ox.ac.uk.0.dat"
#define FFTW_WISDOM_FLAG FFTW_ESTIMATE
#define BLAS_UNDERSCORE
#define USE_FFTW

#endif // _CONFIG_HPP
