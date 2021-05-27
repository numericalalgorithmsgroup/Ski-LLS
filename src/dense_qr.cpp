// Author Zhen Shao
#include "Random.hpp"
#include "flapack.h"

typedef ptrdiff_t Long;

// Dense QR factorization of A using LAPACK dgeqrf
void dense_qr(
  Mat_d& A /* Input matrix, on output, stores R and Q */
  )
{
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
  dgeqrf(&m, &n, A.data(), &m, tau, workspace, &workspace_size, &info);

  /* Compute */

  workspace_size = (Long)wsize_d;
  workspace = (double *)malloc(sizeof(double) * workspace_size);
  dgeqrf(&m, &n, A.data(), &m, tau, workspace, &workspace_size, &info);

  free(workspace);
  free(tau);
}