#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _mpi_use
#include <mpi.h>
#else
#include "parameter.h"
#endif
extern int omp_get_thread_num(void);
extern int omp_get_max_threads(void);
#include "global.h"
#include "../src/mVMC/workspace.c"
#include "../src/mVMC/gauleg.c"
#include "../src/mVMC/legendrepoly.c"
static int reductions = 0;
#ifdef _mpi_use
static int counted_allreduce(const void *send, void *recv, int count,
                             MPI_Datatype type, MPI_Op op, MPI_Comm comm) {
  ++reductions;
  return MPI_Allreduce(send, recv, count, type, op, comm);
}
#define MPI_Allreduce counted_allreduce
#endif
#include "../src/mVMC/qp.c"
#ifdef _mpi_use
#undef MPI_Allreduce
#endif
#include "../src/mVMC/splitloop.c"


static int failures = 0;
#define CHECK(ok, message) do { if (!(ok)) { \
  fprintf(stderr, "GC projection FAIL: G=%d Wcase=%d N=%d %s\n", \
          NSPGaussLeg, control, n, message); ++failures; } } while (0)

static int same_bits(double complex x, double complex y) {
  return memcmp(&x, &y, sizeof x) == 0;
}

int main(int argc, char **argv) {
  const int nodes[] = {1, 2, 4, 8, 16};
  int rank = 0, size = 1, g, control, n = -1;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc; (void)argv;
#endif
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  NThread = omp_get_max_threads();
  initializeWorkSpaceAll();
  Nsite = 4; Nsite2 = 8; NMPTrans = 4; NSPStot = 0; FlagOptTrans = 0;
  ParaQPTrans = malloc(4 * sizeof *ParaQPTrans);
  for (g = 0; g < 5; ++g) {
    int start, end, q;
    double complex *values;
    NSPGaussLeg = nodes[g];
    NQPFull = NQPFix = 4 * NSPGaussLeg;
    QPFullWeight = calloc(2 * NQPFull + 5 * NSPGaussLeg, sizeof *QPFullWeight);
    QPFixWeight = QPFullWeight + NQPFull;
    SPGLCos = QPFixWeight + NQPFix;
    SPGLSin = SPGLCos + NSPGaussLeg;
    SPGLCosSin = SPGLSin + NSPGaussLeg;
    SPGLCosCos = SPGLCosSin + NSPGaussLeg;
    SPGLSinSin = SPGLCosCos + NSPGaussLeg;
    SplitLoop(&start, &end, NQPFull, rank, size);
    values = calloc(NQPFull, sizeof *values);
    for (control = 0; control < 4; ++control) {
      FlagGrandCanonical = control != 3;
      ParaQPTrans[0] = 0.25; ParaQPTrans[1] = -0.25 * I;
      ParaQPTrans[2] = -0.25; ParaQPTrans[3] = 0.25 * I;
      if (control == 1) ParaQPTrans[2] = nextafter(-0.25, 0.0);
      if (control == 2) {
        ParaQPTrans[0] = 1.0; ParaQPTrans[1] = 0x1p-53;
        ParaQPTrans[2] = -1.0; ParaQPTrans[3] = 0.0;
      }
      InitQPWeight();
      /* Nonidentical components deliberately leave finite residuals. The
       * override contract depends on W and the sector, not on their size. */
      for (q = start; q < end; ++q)
        values[q - start] = (1.0 + 0.031 * q) + I * (0.13 - 0.017 * q);
      for (n = 0; n <= Nsite2; n += 2) {
        double complex old = CalculateIP_fcmp(values, start, end, MPI_COMM_WORLD);
        double complex oldLog = CalculateLogIP_fcmp(values, start, end, MPI_COMM_WORLD);
        double complex now, nowLog;
        reductions = 0;
        now = CalculateIPGC_fcmp(values, n, NULL, start, end, MPI_COMM_WORLD);
        CHECK(reductions == (size > 1), "value collective count changed");
        reductions = 0;
        nowLog = CalculateLogIPGC_fcmp(values, n, NULL, start, end, MPI_COMM_WORLD);
        CHECK(reductions == (size > 1), "log collective count changed");
        if (control == 0 && (n == 0 || n == Nsite2)) {
          CHECK(same_bits(now, 0.0 + 0.0 * I), "boundary is not exact +0");
          CHECK(creal(nowLog) == -INFINITY && isfinite(cimag(nowLog)), "zero log");
          CHECK(cabs(old) > 0.0, "missing nonzero residual control");
        } else {
          CHECK(same_bits(now, old), "value changed outside zero boundary");
          CHECK(same_bits(nowLog, oldLog), "log changed outside zero boundary");
        }
      }
      if (control == 0) {
        const double faults[] = {NAN, INFINITY, -INFINITY};
        int fault, part;
        for (fault = 0; fault < 3; ++fault) for (part = 0; part < 2; ++part) {
          for (q = start; q < end; ++q) values[q - start] = 1.0;
          if (start == 0 && end > 0) {
            const double parts[2] = {part == 0 ? faults[fault] : 1.0,
                                     part == 1 ? faults[fault] : 0.0};
            memcpy(values, parts, sizeof parts);
          }
          for (n = 0; n <= Nsite2; n += Nsite2) {
            const double complex now = CalculateIPGC_fcmp(values, n, NULL, start, end, MPI_COMM_WORLD);
            CHECK(!isfinite(creal(now)) || !isfinite(cimag(now)), "nonfinite hidden");
          }
        }
      }
    }
    free(values); free(QPFullWeight);
  }
  free(ParaQPTrans); FreeWorkSpaceAll();
#ifdef _mpi_use
  {
    int total;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    failures = total;
  }
  MPI_Finalize();
#endif
  if (rank == 0 && !failures) puts("GC projection: PASS (sectors, bitwise, faults, collectives)");
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
