#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#define expit(a) (exp(a) / (1 + exp(a)))

void
rsnmm(int *size, int *tmax,
      double *ty, double *tmod, double *tavail, double *tstate,
      double *beta, double *eta, double *mu, double *theta,
      double *coefavail,  double *coefstate, double *coeferr,
      int *avail, double *base, int *state, int *a, double *prob,
      double *y, double *err, double *statec, double *ac, double *availc)
{
  int i, j, n = *size, T = *tmax;
  double r, q, ym;
  GetRNGstate();
  for (i = 0; i < n; i++) {
    for (j = 1; j < T; j++) {
      /* probability of availabilty */
      r = expit(coefavail[0]
                + coefavail[1] * tavail[j]
                + coefavail[2] * a[i*T + j - 1]
                + coefavail[3] * y[i*T + j - 1]);
      /* availability - uncentered and centered */
      avail[i*T + j] = (int) rbinom(1, r);
      availc[i*T + j] = avail[i*T + j] - r;
      /* probability that binary state is +1 */
      q = expit(coefstate[0]
                + coefstate[1] * tstate[j]
                + coefstate[2] * base[i*T + j - 1]
                + coefstate[3] * state[i*T + j - 1]
                + coefstate[4] * a[i*T + j - 1]);
      /* binary state on {-1, 1} - uncentered and centered */
      state[i*T + j] = rbinom(1, q) < 1 ? -1 : 1;
      statec[i*T + j] = state[i*T + j] - (q - (1 - q));
      /* treatment probability */
      prob[i*T + j] = avail[i*T + j] * expit(eta[0]
                                          + eta[1] * base[i*T + j]
                                          + eta[2] * state[i*T + j]
                                          + eta[3] * a[i*T + j - 1]
                                          + eta[4] * y[i*T + j - 1]);
      /* treatment indicator - uncentered and centered */
      a[i*T + j] = (int) rbinom(1, prob[i*T + j]);
      ac[i*T + j] = a[i*T + j] - prob[i*T + j];
      /* conditional mean response */
      ym = mu[0]
        + mu[1] * ty[j]  /* pre-evaluated time function */
        + mu[2] * base[i*T + j]
        + ac[i*T + j] * (beta[0]
                         + beta[1] * tmod[j] /* pre-evaluated time function */
                         + beta[2] * base[i*T + j]
                         + beta[3] * state[i*T + j]
                         + beta[4] * a[i*T + j - 1])
        + ac[i*T + j - 1] * (beta[5]
                             + beta[6] * tmod[j - 1]
                             + beta[7] * base[i*T + j - 1]
                             + beta[8] * state[i*T + j - 1])
        + theta[0] * availc[i*T + j]
        + theta[1] * statec[i*T + j]
        + theta[2] * availc[i*T + j - 1]
        + theta[3] * statec[i*T + j - 1];
      /* error */
      err[i*T + j] += *coeferr * err[i*T + j - 1];
      /* response */
      y[i*T + j] = ym + err[i*T + j];
    }
  }
  PutRNGstate();
  return;
}
