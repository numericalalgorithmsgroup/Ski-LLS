# External library Locations

# If your environment is set up correctly, MKLROOT should point to the right directory
LIBS_LAPACK =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a \
-Wl,--end-group -lpthread -lm -ldl

SUITESPARSEROOT ?= /fserver/zhens/testInstall/Dependencies/SuiteSparse_intel/SuiteSparse-5.6.0/# Change this
SPARSE_INCLUDE ?= $(SUITESPARSEROOT)/include
LIBS_SPARSE = -L${SUITESPARSEROOT}/lib -lcxsparse -lcholmod -lspqr

FFTWROOT ?= /fserver/zhens/testInstall/Dependencies/fftw-3.3.8/# Change this
FFTW_INCLUDE ?= ${FFTWROOT}/api
LIBS_FFTW = ${FFTWROOT}/.libs/libfftw3.a

BOOSTROOT ?= /fserver/zhens/testInstall/Dependencies/boost_1_72_0/# Change this

# compiler
CXX = icc
CC      = icc
CXXFLAGS =  -DMKL_ILP64  -m64  -I"${MKLROOT}/include" -Wall -O3 -march=native -fPIC -std=c++11
CFLAGS = -O3 -fopenmp

