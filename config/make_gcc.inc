# External library Locations

LAPACKROOT ?= /fserver/zhens/testInstall/Dependencies/lapack-3.9.0/# Change this
LIBS_LAPACK = ${LAPACKROOT}/liblapack.a ${LAPACKROOT}/librefblas.a -lgfortran

SUITESPARSEROOT ?= /fserver/zhens/testInstall/Dependencies/SuiteSparse-5.6.0/# Change this
SPARSE_INCLUDE ?= $(SUITESPARSEROOT)/include
LIBS_SPARSE = -L${SUITESPARSEROOT}/lib -lcxsparse -lcholmod -lspqr

FFTWROOT ?= /fserver/zhens/testInstall/Dependencies/fftw-3.3.8/# Change this
FFTW_INCLUDE ?= ${FFTWROOT}/api
LIBS_FFTW = ${FFTWROOT}/.libs/libfftw3.a

BOOSTROOT ?= /fserver/zhens/testInstall/Dependencies/boost_1_72_0/# Change this

# compiler
CXX = g++
CC      = gcc
CXXFLAGS = -Wall -O3 -march=native -fPIC -std=c++11
CFLAGS = -O3 -fopenmp
