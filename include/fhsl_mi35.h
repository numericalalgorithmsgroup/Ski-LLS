/* C interface for HSL_MI35 incomplete Cholesky solver
 * HSL. A collection of Fortran codes for large scale scientific computation. http://www.hsl.rl.ac.uk/
 */
#ifndef _FHSL_MI35_H
#define _FHSL_MI35_H

#include "config_sk.h"

typedef INT fint; 

#ifdef __cplusplus
extern "C" {
#endif

#define hsl_mi35_factorize FWRAPPER(cmi35_factorize,CMI35_FACTORIZE)
  extern void hsl_mi35_factorize(const fint *m, const fint *n, const fint *nnz,
    const fint *ptr /*[n+1]*/, const fint *row /*[nnz]*/, const double *val /*nnz*/, 
    const fint *lsize, const fint *rsize, void **pkeep, const double rcontrol[20], 
    const fint icontrol[20], fint *ifail);

#define hsl_mi35_finalise FWRAPPER(cmi35_finalise,CMI35_FINALISE)
  extern void hsl_mi35_finalise(void **pkeep, fint *ifail);

#define hsl_mi35_solve FWRAPPER(cmi35_solve,CMI35_SOLVE)
  extern void hsl_mi35_solve(const fint *itrans, const fint *n, void **pkeep, 
    const double *z /*[n]*/, double *y /*[n]*/, fint *ifail);

#ifdef __cplusplus
}
#endif

#endif

