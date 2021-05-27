#ifndef __HQRRP_H
#define __HQRRP_H

#include "config_sk.h"

#ifdef __cplusplus
extern "C" {
#endif

void dgeqp4_INT( INT * m, INT * n, double * A, INT * lda, INT * jpvt, double * tau,
         double * work, INT * lwork, INT * info );

INT NoFLA_HQRRP_WY_blk_var4_INT( INT m_A, INT n_A, double * buff_A, INT ldim_A,
        INT * buff_jpvt, double * buff_tau,
        INT nb_alg, INT pp, INT panel_pivoting );

#ifdef __cplusplus
}
#endif

#endif

