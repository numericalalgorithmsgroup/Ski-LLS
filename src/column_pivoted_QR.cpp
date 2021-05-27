// Author Zhen Shao
#include "Random.hpp"
#include "flapack.h"
#include "hqrrp.h"

typedef ptrdiff_t Long;
void column_pivoted_qr(Mat_d& A, long* E)
{
  // return R as upper part of A, E column permutation of A, use LAPACK
  // if size(A) = (m,n), E has to contain at least min(m,n) elements.

  Long m = A.m();
  Long n = A.n();
  double *tau;
  Long workspace_size, info;
  double *workspace;
  double wsize_d;
  tau = (double*)malloc(n * sizeof(double));

  /* Query workspace size */
  workspace_size = -1;
  workspace = &wsize_d;
  dgeqp3(&m, &n, A.data(), &m, E, tau, workspace, &workspace_size, &info);

  /* Compute */

  workspace_size = (Long)wsize_d;
  workspace = (double *)malloc(sizeof(double) * workspace_size);
  dgeqp3(&m, &n, A.data(), &m, E, tau, workspace, &workspace_size, &info);

  free(workspace);
  free(tau);
}

void column_pivoted_qr(Mat_d& A, long* E, double* tau)
{
  // also return tau to reconstruct the Q matrix

  Long m = A.m();
  Long n = A.n();
  Long workspace_size, info;
  double *workspace;
  double wsize_d;

  /* Query workspace size */
  workspace_size = -1;
  workspace = &wsize_d;
  dgeqp4_INT(&m, &n, A.data(), &m, E, tau, workspace, &workspace_size, &info);

  /* Compute */

  workspace_size = (Long)wsize_d;
  workspace = (double *)malloc(sizeof(double) * workspace_size);
  dgeqp4_INT(&m, &n, A.data(), &m, E, tau, workspace, &workspace_size, &info);


  free(workspace);
}
