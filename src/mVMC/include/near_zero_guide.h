#ifndef _NEAR_ZERO_GUIDE_H
#define _NEAR_ZERO_GUIDE_H

#include <complex.h>

#ifdef _mpi_use
#include <mpi.h>
typedef MPI_Comm NearZeroGuideComm;
#else
typedef int NearZeroGuideComm;
#endif

/*
 * Near-zero guide for the legacy first-step Lanczos sampler.
 *   P = sum_q Re(w_q) pf_q
 *   A = sum_q |Re(w_q) pf_q|
 *   q_eps / |J|^2 = max(P^2, eps A^2)
 */
typedef struct {
  double abs_ip;
  double abs_sum;
  double log_abs_ip;
  double log_guide;
  int floored;
} NearZeroGuide;

/* pfM is the rank-local slice. pfM[i] is paired with
 * qpFullWeight[i + qpStart]. Returns 0 on success and 1 on invalid input. */
int NearZeroGuideEvaluate(const double complex *qpFullWeight,
                          const double *pfM,
                          int qpStart,
                          int qpEnd,
                          double sqrtEps,
                          NearZeroGuideComm comm,
                          NearZeroGuide *out);

/* |psi|^2 / q_eps from log|P|, log|J|, and stored log q_eps. */
double NearZeroGuideWeight(double logAbsIp,
                           double logProj,
                           double logStoredGuideSq);

enum {
  NZG_SAMPLES = 0,
  NZG_FLOORED,
  NZG_EXACT_ZERO,
  NZG_NUMERIC_SKIPPED,
  NZG_SUM_W,
  NZG_SUM_W2,
  NZG_MIN_W,
  NZG_MAX_RATIO,
  NZG_COUNT
};

#define NZG_STAT_SIZE 8

void NearZeroGuideStatReset(double *stat);

#endif
