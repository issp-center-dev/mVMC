#include <math.h>
#include <stddef.h>

#include "./include/near_zero_guide.h"

int NearZeroGuideEvaluate(const double complex *qpFullWeight,
                          const double *pfM,
                          int qpStart,
                          int qpEnd,
                          double sqrtEps,
                          NearZeroGuideComm comm,
                          NearZeroGuide *out) {
  double local[2] = {0.0, 0.0};
  double total[2];
  int qpidx;

  if (out == NULL || qpFullWeight == NULL || pfM == NULL ||
      qpEnd < qpStart || !(sqrtEps >= 0.0) || !isfinite(sqrtEps)) {
    return 1;
  }

  for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
    const double term =
        creal(qpFullWeight[qpidx + qpStart]) * pfM[qpidx];
    local[0] += term;
    local[1] += fabs(term);
  }

  total[0] = local[0];
  total[1] = local[1];
#ifdef _mpi_use
  {
    int size = 1;
    MPI_Comm_size(comm, &size);
    if (size > 1) {
      MPI_Allreduce(local, total, 2, MPI_DOUBLE, MPI_SUM, comm);
    }
  }
#else
  (void)comm;
#endif

  out->abs_ip = fabs(total[0]);
  out->abs_sum = total[1];
  out->log_abs_ip = log(out->abs_ip);
  out->floored = (out->abs_ip < sqrtEps * out->abs_sum);
  out->log_guide = log(fmax(out->abs_ip, sqrtEps * out->abs_sum));
  return 0;
}

double NearZeroGuideWeight(double logAbsIp,
                           double logProj,
                           double logStoredGuideSq) {
  return exp(2.0 * (logAbsIp + logProj) - logStoredGuideSq);
}

void NearZeroGuideStatReset(double *stat) {
  int i;
  for (i = 0; i < NZG_STAT_SIZE; i++) {
    stat[i] = 0.0;
  }
  stat[NZG_MIN_W] = INFINITY;
}
