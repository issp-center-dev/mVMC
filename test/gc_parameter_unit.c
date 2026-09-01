#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

#include "global.h"
#include "SFMT.h"

void BFInitParameters(void) {}
void BFRefreshRealLookupTables(void) {}

#include "../src/mVMC/parameter.c"

#define PROJECTION_COUNT 18
#define SLATER_COUNT 2

static int failures = 0;

#define CHECK(condition, ...)                                                  \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "GCParameter_Unit FAIL: ");                             \
      fprintf(stderr, __VA_ARGS__);                                            \
      fprintf(stderr, "\n");                                                 \
      failures++;                                                              \
    }                                                                          \
  } while (0)

static void reset_parameters(double complex *projection,
                             double complex *slater) {
  int index;
  for (index = 0; index < PROJECTION_COUNT; index++) {
    projection[index] = 0.07 * (index + 1) + 0.013 * (index + 2) * I;
  }
  slater[0] = 2.0 + 0.0 * I;
  slater[1] = -0.25 + 0.75 * I;
}

static void test_gc_preserves_sector_weights(void) {
  double complex projection[PROJECTION_COUNT];
  double complex projectionSnapshot[PROJECTION_COUNT];
  double complex slater[SLATER_COUNT];
  double complex slaterSnapshot[SLATER_COUNT];

  NProj = PROJECTION_COUNT;
  NGutzwillerIdx = 1;
  NJastrowIdx = 1;
  NSpinJastrowIdx = 0;
  NDoublonHolon2siteIdx = 1;
  NDoublonHolon4siteIdx = 1;
  NSlater = SLATER_COUNT;
  NOptTrans = 0;
  FlagRBM = 0;
  NProjBF = 0;
  FlagOptTrans = 0;
  FlagShiftGJ = 1;
  FlagShiftDH2 = 1;
  FlagShiftDH4 = 1;
  NPara = 0;
  Para = NULL;
  Proj = projection;
  Slater = slater;
  OptTrans = NULL;

  reset_parameters(projection, slater);
  memcpy(projectionSnapshot, projection, sizeof(projectionSnapshot));
  memcpy(slaterSnapshot, slater, sizeof(slaterSnapshot));
  FlagGrandCanonical = 1;
  SyncModifiedParameter(MPI_COMM_WORLD);
  CHECK(memcmp(projection, projectionSnapshot, sizeof(projection)) == 0,
        "GC correlation parameters were shifted");
  CHECK(memcmp(slater, slaterSnapshot, sizeof(slater)) == 0,
        "GC pairing scale was rescaled");

  reset_parameters(projection, slater);
  memcpy(projectionSnapshot, projection, sizeof(projectionSnapshot));
  FlagGrandCanonical = 0;
  SyncModifiedParameter(MPI_COMM_WORLD);
  CHECK(memcmp(projection, projectionSnapshot, sizeof(projection)) != 0,
        "canonical correlation shifts were unexpectedly disabled");
  CHECK(fabs(cabs(slater[0]) - 4.0) < 1.0e-14,
        "canonical Slater normalization changed: %.17g", cabs(slater[0]));
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  test_gc_preserves_sector_weights();
#ifdef _mpi_use
  MPI_Finalize();
#endif
  if (failures != 0) return EXIT_FAILURE;
  printf("GC parameter checks: PASS\n");
  return EXIT_SUCCESS;
}
