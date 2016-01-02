#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#define expit(a) (exp(a) / (1 + exp(a)))

void
rsnmm(int *size, int *times, double *beta, double *eta, double *mu, double *xi,
      double *cv, double *b, int *v, int *a, double *p, double *y,
      double *err, double *vc, double *ac, double *sd)
{
  int i, j, k, n = *size, T = *times;
  double q, ym, s = *sd;
  GetRNGstate();
  for (i = 0; i < n; i++) {
    for (j = 1; j < T; j++) {
      q = expit(cv[0]
                + cv[1] * b[i*T + j - 1]
                + cv[2] * v[i*T + j - 1]
                + cv[3] * a[i*T + j - 1]);
      v[i*T + j] = rbinom(1, q) < 1 ? -1 : 1;
      vc[i*T + j] = v[i*T + j] - (q - (1 - q));
      p[i*T + j] = expit(eta[0]
                         + eta[1] * b[i*T + j]
                         + eta[2] * v[i*T + j]
                         + eta[3] * a[i*T + j - 1]
                         + eta[4] * y[i*T + j - 1]);
      a[i*T + j] = (int) rbinom(1, p[i*T + j]);
      ac[i*T + j] = a[i*T + j] - p[i*T + j];
      ym = mu[0] + mu[1] * b[i*T + j]
        + ac[i*T + j] * (beta[0]
                         + beta[1] * b[i*T + j]
                         + beta[2] * v[i*T + j]
                         + beta[3] * a[i*T + j - 1])
        + ac[i*T + j - 1] * (beta[4]
                             + beta[5] * b[i*T + j]
                             + beta[6] * v[i*T + j - 1])
        + xi[0] * vc[i*T + j]
        + xi[1] * err[i*T + j - 1]
        + xi[2] * vc[i*T + j - 1];
      y[i*T + j] = ym + err[i*T + j];
    }
  }
  PutRNGstate();
  return;
}
