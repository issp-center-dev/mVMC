#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../src/mVMC/include/near_zero_guide.h"

static int failures = 0;

#define CHECK(cond, msg)                                                     \
  do {                                                                       \
    if (!(cond)) {                                                           \
      failures++;                                                            \
      fprintf(stderr, "FAIL: %s\n", (msg));                                \
    }                                                                        \
  } while (0)

static void evaluate_all(const double complex *weight,
                         const double *pf,
                         int n,
                         double sqrtEps,
                         NearZeroGuideComm comm,
                         NearZeroGuide *guide) {
  int start = 0;
  int end = n;
#ifdef _mpi_use
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  start = (n * rank) / size;
  end = (n * (rank + 1)) / size;
#endif
  CHECK(NearZeroGuideEvaluate(weight, pf + start, start, end, sqrtEps,
                              comm, guide) == 0,
        "evaluate returns 0");
}

int main(int argc, char **argv) {
  NearZeroGuideComm comm = 0;
  const double complex one[4] = {1.0, 1.0, 1.0, 1.0};
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
#else
  (void)argc;
  (void)argv;
#endif

  {
    const double pf[4] = {0.3, -0.1, 0.25, 0.05};
    NearZeroGuide guide;
    evaluate_all(one, pf, 4, 0.0, comm, &guide);
    CHECK(fabs(guide.abs_ip - 0.5) < 4e-16 &&
              fabs(guide.abs_sum - 0.7) < 4e-16,
          "P and A sums");
    CHECK(guide.log_guide == guide.log_abs_ip && guide.floored == 0,
          "eps=0 keeps the legacy log weight");
#ifdef _mpi_use
    {
      NearZeroGuide serial;
      CHECK(NearZeroGuideEvaluate(one, pf, 0, 4, 0.0, MPI_COMM_SELF,
                                  &serial) == 0,
            "serial evaluate");
      CHECK(fabs(guide.abs_ip - serial.abs_ip) <= 4e-16 &&
                fabs(guide.abs_sum - serial.abs_sum) <= 4e-16,
            "rank split matches the serial sum");
    }
#endif
  }
  {
    const double pf[2] = {3.2066864e-3, -3.2083847e-3};
    const double p = pf[0] + pf[1];
    const double a = fabs(pf[0]) + fabs(pf[1]);
    const double eps = 1e-6;
    const double expected = eps * a * a / (p * p);
    NearZeroGuide guide;
    evaluate_all(one, pf, 2, sqrt(eps), comm, &guide);
    CHECK(guide.floored == 1,
          "stage-one-like cancellation is floored at eps=1e-6");
    CHECK(fabs(exp(2.0 * (guide.log_guide - guide.log_abs_ip)) -
                   expected) <
              1e-9 * expected,
          "q/P^2 = eps A^2/P^2");
    CHECK(fabs(guide.log_guide - log(sqrt(eps) * a)) < 1e-15,
          "floored branch uses sqrt(eps) A");
  }
  {
    const double pf[2] = {0.6, 0.4};
    NearZeroGuide guide;
    evaluate_all(one, pf, 2, 1e-3, comm, &guide);
    CHECK(guide.floored == 0 && guide.log_guide == guide.log_abs_ip,
          "same-sign components are unchanged");
  }
  {
    const double pf[2] = {1.0, -1.0};
    NearZeroGuide guide;
    evaluate_all(one, pf, 2, 1e-3, comm, &guide);
    CHECK(guide.abs_ip == 0.0 && isinf(guide.log_abs_ip) &&
              guide.log_abs_ip < 0.0,
          "P=0 gives log|P|=-inf");
    CHECK(isfinite(guide.log_guide) &&
              fabs(guide.log_guide - log(2e-3)) < 1e-15 &&
              guide.floored == 1,
          "guide stays finite for P=0 and A>0");
  }
  {
    const double pf[2] = {0.0, 0.0};
    NearZeroGuide guide;
    evaluate_all(one, pf, 2, 1e-3, comm, &guide);
    CHECK(isinf(guide.log_guide) && guide.log_guide < 0.0,
          "A=0 keeps log guide=-inf");
  }
  {
    const double pf[2] = {3.2066864e-3, -3.2083847e-3};
    const double pf7[2] = {7.0 * pf[0], 7.0 * pf[1]};
    NearZeroGuide guide;
    NearZeroGuide guide7;
    evaluate_all(one, pf, 2, 1e-3, comm, &guide);
    evaluate_all(one, pf7, 2, 1e-3, comm, &guide7);
    CHECK(fabs((guide.log_guide - guide.log_abs_ip) -
                   (guide7.log_guide - guide7.log_abs_ip)) <
              1e-12,
          "ratio is scale invariant");
  }
  {
    const double pf[2] = {3.2066864e-3, -3.2083847e-3};
    const double x = 0.37;
    NearZeroGuide guide;
    double weight;
    double p2;
    double q;
    evaluate_all(one, pf, 2, 1e-3, comm, &guide);
    weight = NearZeroGuideWeight(guide.log_abs_ip, x,
                                 2.0 * (x + guide.log_guide));
    p2 = guide.abs_ip * guide.abs_ip;
    q = fmax(p2, 1e-6 * guide.abs_sum * guide.abs_sum);
    CHECK(fabs(weight - p2 / q) < 1e-12 * (p2 / q),
          "measurement weight identity");
  }
  {
    const double pf[1] = {1.0};
    NearZeroGuide guide;
    CHECK(NearZeroGuideEvaluate(one, pf, 0, 1, -1.0, comm, &guide) == 1,
          "negative sqrt eps rejected");
    CHECK(NearZeroGuideEvaluate(one, pf, 1, 0, 0.0, comm, &guide) == 1,
          "reversed range rejected");
    CHECK(NearZeroGuideEvaluate(one, pf, 0, 1, 0.0, comm, NULL) == 1,
          "NULL output rejected");
  }
  {
    double stat[NZG_STAT_SIZE];
    NearZeroGuideStatReset(stat);
    CHECK(stat[NZG_SAMPLES] == 0.0 &&
              stat[NZG_NUMERIC_SKIPPED] == 0.0 &&
              stat[NZG_SUM_W] == 0.0 && isinf(stat[NZG_MIN_W]) &&
              stat[NZG_MIN_W] > 0.0 && stat[NZG_MAX_RATIO] == 0.0,
          "statistics reset");
  }

#ifdef _mpi_use
  MPI_Finalize();
#endif
  if (failures != 0) {
    fprintf(stderr, "%d failure(s)\n", failures);
    return EXIT_FAILURE;
  }
  printf("near_zero_guide_unit: all checks passed\n");
  return EXIT_SUCCESS;
}
