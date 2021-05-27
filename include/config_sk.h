#ifndef __CONFIG_SK_H
#define __CONFIG_SK_H

/* the type of integer across Ski-LLS, must match all the libraries */
typedef long INT;

/* name mangling for Fortran code (BLAS/LAPACK) */
#define FWRAPPER(func,FUNC) func ## _

#endif

