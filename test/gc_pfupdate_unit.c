#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

extern int omp_get_thread_num(void);

#include "global.h"
#include "blas_externs.h"

/* Production implementations share the legacy process-wide globals owned by
 * global.h, so exercise the exact translation-unit chain used by vmc.out. */
#include "../src/mVMC/gc_size.c"
#include "../src/mVMC/workspace.c"
#include "../src/mVMC/matrix_gc.c"
#include "../src/mVMC/gc_config.c"
#include "../src/mVMC/pfupdate_gc.c"

#define MAX_DIMENSION 8
#define MAX_QP 2

typedef struct {
  double complex inv[MAX_QP * MAX_DIMENSION * MAX_DIMENSION];
  double complex pf[MAX_QP];
} PfState;

static int failures = 0;

#define CHECK(condition, ...)                                                  \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "GCPfUpdate_Unit FAIL: ");                             \
      fprintf(stderr, __VA_ARGS__);                                            \
      fprintf(stderr, "\n");                                                 \
      failures++;                                                              \
    }                                                                          \
  } while (0)

static int popcount(unsigned int value) {
  int count = 0;
  while (value != 0U) {
    count += (int)(value & 1U);
    value >>= 1;
  }
  return count;
}

static int close_complex(const double complex actual,
                         const double complex expected) {
  return cabs(actual - expected) <= 2.0e-9 * (1.0 + cabs(expected));
}

static void fill_slater_block(double complex *slater, const int dimension,
                              const double scale, const double phase) {
  int row;
  memset(slater, 0,
         (size_t)dimension * (size_t)dimension * sizeof(*slater));
  for (row = 0; row < dimension; row++) {
    int column;
    for (column = row + 1; column < dimension; column++) {
      const double realPart =
          scale * (0.37 + 0.29 * (row + 1) + 0.11 * (column + 1) +
                   0.017 * (row + 1) * (column + 2));
      const double imagPart =
          scale * (phase + 0.13 * (row + 1) - 0.071 * (column + 1) +
                   0.019 * (row + 2) * (column + 1));
      const double complex value = realPart + imagPart * I;
      slater[(size_t)row * (size_t)dimension + (size_t)column] = value;
      slater[(size_t)column * (size_t)dimension + (size_t)row] = -value;
    }
  }
}

static void setup_dimension(const int dimension) {
  const size_t block = (size_t)dimension * (size_t)dimension;
  Nsite = dimension / 2;
  Nsite2 = dimension;
  Nsize = dimension;
  NsizeMax = dimension;
  NQPFull = MAX_QP;
  LapackLWork = 256;
  SlaterElm = realloc(SlaterElm, MAX_QP * block * sizeof(*SlaterElm));
  InvM = realloc(InvM, MAX_QP * block * sizeof(*InvM));
  PfM = realloc(PfM, MAX_QP * sizeof(*PfM));
  CHECK(SlaterElm != NULL && InvM != NULL && PfM != NULL,
        "allocation failed for dimension=%d", dimension);
  if (SlaterElm == NULL || InvM == NULL || PfM == NULL) exit(EXIT_FAILURE);
  fill_slater_block(SlaterElm, dimension, 1.0, 0.031);
  fill_slater_block(SlaterElm + block, dimension, 1.43, -0.047);
}

static void save_state(PfState *state) {
  const size_t matrixCount =
      (size_t)MAX_QP * (size_t)NsizeMax * (size_t)NsizeMax;
  memcpy(state->inv, InvM, matrixCount * sizeof(*InvM));
  memcpy(state->pf, PfM, MAX_QP * sizeof(*PfM));
}

static void restore_state(const PfState *state) {
  const size_t matrixCount =
      (size_t)MAX_QP * (size_t)NsizeMax * (size_t)NsizeMax;
  memcpy(InvM, state->inv, matrixCount * sizeof(*InvM));
  memcpy(PfM, state->pf, MAX_QP * sizeof(*PfM));
}

static void initialize_old_state(const int ncur, const int *eleIdx,
                                 PfState *oldState, const char *label) {
  const double complex sentinel = 83.0 - 29.0 * I;
  const size_t matrixCount =
      (size_t)MAX_QP * (size_t)NsizeMax * (size_t)NsizeMax;
  size_t i;
  for (i = 0; i < matrixCount; i++) InvM[i] = sentinel;
  PfM[0] = -101.0 + 3.0 * I;
  PfM[1] = -102.0 + 4.0 * I;
  CHECK(CalculateMAllGC_fcmp(ncur, eleIdx, 0, MAX_QP) == GC_MALL_OK,
        "%s old rebuild failed", label);
  save_state(oldState);
}

static void rebuild_candidate(const int ncur, const int *eleIdx,
                              PfState *candidateState, const char *label) {
  CHECK(CalculateMAllGC_fcmp(ncur, eleIdx, 0, MAX_QP) == GC_MALL_OK,
        "%s candidate rebuild failed", label);
  save_state(candidateState);
}

static void check_pf_values(const double complex *actual,
                            const PfState *expected, const char *label) {
  int qp;
  for (qp = 0; qp < MAX_QP; qp++) {
    CHECK(close_complex(actual[qp], expected->pf[qp]),
          "%s qp=%d Pf actual=(%.17g,%.17g) expected=(%.17g,%.17g)",
          label, qp, creal(actual[qp]), cimag(actual[qp]),
          creal(expected->pf[qp]), cimag(expected->pf[qp]));
  }
}

static void check_updated_state(const int ncur, const PfState *expected,
                                const PfState *oldState, const char *label) {
  int qp;
  for (qp = 0; qp < MAX_QP; qp++) {
    int row;
    CHECK(close_complex(PfM[qp], expected->pf[qp]),
          "%s qp=%d updated Pf mismatch", label, qp);
    for (row = 0; row < NsizeMax; row++) {
      int column;
      for (column = 0; column < NsizeMax; column++) {
        const size_t index =
            (size_t)qp * (size_t)NsizeMax * (size_t)NsizeMax +
            (size_t)row * (size_t)NsizeMax + (size_t)column;
        if (row < ncur && column < ncur) {
          CHECK(close_complex(InvM[index], expected->inv[index]),
                "%s qp=%d inverse row=%d col=%d", label, qp, row,
                column);
        } else {
          CHECK(InvM[index] == oldState->inv[index],
                "%s qp=%d inactive tail row=%d col=%d changed", label, qp,
                row, column);
        }
      }
    }
  }
}

static void check_candidate_noncommit(const PfState *oldState,
                                      const char *label) {
  const size_t matrixCount =
      (size_t)MAX_QP * (size_t)NsizeMax * (size_t)NsizeMax;
  CHECK(memcmp(InvM, oldState->inv, matrixCount * sizeof(*InvM)) == 0,
        "%s candidate changed InvM", label);
  CHECK(memcmp(PfM, oldState->pf, MAX_QP * sizeof(*PfM)) == 0,
        "%s candidate changed PfM", label);
}

static void check_nrow_candidate(const int n, const int *msa,
                                 const int *rsa, const int *eleIdx,
                                 const int ncur, const PfState *expected,
                                 const PfState *oldState, const char *label) {
  double complex vec[3 * MAX_DIMENSION];
  double complex w[MAX_DIMENSION];
  double complex smallMat[36];
  int pfIWork[6];
  double complex pfWork[36];
  double pfRWork[36];
  double complex actual[MAX_QP];
  int qp;
  restore_state(oldState);
  for (qp = 0; qp < MAX_QP; qp++) {
    actual[qp] = CalculateNewPfMNGC(
        qp, n, msa, rsa, eleIdx, ncur, vec, w, smallMat, pfIWork,
        pfWork, 36, pfRWork);
  }
  check_candidate_noncommit(oldState, label);
  check_pf_values(actual, expected, label);
}

static void test_hop_and_nrow(const int ncur, const int *configuration,
                              const unsigned int mask) {
  int original[MAX_DIMENSION];
  int ma;
  memcpy(original, configuration, (size_t)ncur * sizeof(*original));
  for (ma = 0; ma < ncur; ma++) {
    int rsa;
    for (rsa = 0; rsa < Nsite2; rsa++) {
      PfState oldState;
      PfState expected;
      double complex candidate[MAX_QP];
      int msa[1] = {ma};
      int destinations[1] = {rsa};
      char label[128];
      if ((mask & (1U << rsa)) != 0U) continue;
      memcpy(original, configuration, (size_t)ncur * sizeof(*original));
      snprintf(label, sizeof(label), "hop M=%d n=%d mask=%u ma=%d rs=%d",
               Nsite2, ncur, mask, ma, rsa);
      initialize_old_state(ncur, original, &oldState, label);
      original[ma] = rsa;
      CalculateNewPfMHopGC(ma, candidate, original, ncur, 0, MAX_QP);
      check_candidate_noncommit(&oldState, label);
      rebuild_candidate(ncur, original, &expected, label);
      check_pf_values(candidate, &expected, label);
      check_nrow_candidate(1, msa, destinations, original, ncur, &expected,
                           &oldState, label);
      restore_state(&oldState);
      UpdateMAllHopGC(ma, original, ncur, 0, MAX_QP);
      check_updated_state(ncur, &expected, &oldState, label);
      original[ma] = configuration[ma];
    }
  }
}

static void test_two_hop_and_nrow(const int ncur, const int *configuration,
                                  const unsigned int mask) {
  int ma;
  if (ncur < 2 || Nsite2 - ncur < 2) return;
  for (ma = 0; ma < ncur; ma++) {
    int mb;
    for (mb = ma + 1; mb < ncur; mb++) {
      int rsa;
      for (rsa = 0; rsa < Nsite2; rsa++) {
        int rsb;
        if ((mask & (1U << rsa)) != 0U) continue;
        for (rsb = 0; rsb < Nsite2; rsb++) {
          int candidateConfiguration[MAX_DIMENSION];
          PfState oldState;
          PfState expected;
          double complex candidate[MAX_QP];
          int msa[2] = {ma, mb};
          int destinations[2] = {rsa, rsb};
          char label[160];
          if (rsb == rsa || (mask & (1U << rsb)) != 0U) continue;
          memcpy(candidateConfiguration, configuration,
                 (size_t)ncur * sizeof(*candidateConfiguration));
          snprintf(label, sizeof(label),
                   "two-hop M=%d n=%d mask=%u ma=%d mb=%d rs=%d,%d",
                   Nsite2, ncur, mask, ma, mb, rsa, rsb);
          initialize_old_state(ncur, candidateConfiguration, &oldState, label);
          candidateConfiguration[ma] = rsa;
          candidateConfiguration[mb] = rsb;
          CalculateNewPfMTwoHopGC(ma, mb, candidate, candidateConfiguration,
                                  ncur, 0, MAX_QP);
          check_candidate_noncommit(&oldState, label);
          rebuild_candidate(ncur, candidateConfiguration, &expected, label);
          check_pf_values(candidate, &expected, label);
          check_nrow_candidate(2, msa, destinations, candidateConfiguration,
                               ncur, &expected, &oldState, label);
        }
      }
    }
  }
}

static void test_m6_exhaustive(void) {
  unsigned int mask;
  setup_dimension(6);
  for (mask = 0; mask < (1U << 6); mask++) {
    const int ncur = popcount(mask);
    int configuration[MAX_DIMENSION];
    int count = 0;
    int rs;
    if (ncur != 2 && ncur != 4) continue;
    for (rs = 0; rs < 6; rs++) {
      if ((mask & (1U << rs)) != 0U) configuration[count++] = rs;
    }
    test_hop_and_nrow(ncur, configuration, mask);
    test_two_hop_and_nrow(ncur, configuration, mask);
  }
}

static void test_m8_noncontracted_n3(void) {
  int configuration[MAX_DIMENSION] = {0, 1, 2, 3};
  const int msa[3] = {0, 1, 2};
  const int destinations[3] = {4, 5, 6};
  PfState oldState;
  PfState expected;
  const char *label = "M=8 ncur=4 noncontracted n=3";
  setup_dimension(8);
  initialize_old_state(4, configuration, &oldState, label);
  configuration[0] = destinations[0];
  configuration[1] = destinations[1];
  configuration[2] = destinations[2];
  rebuild_candidate(4, configuration, &expected, label);
  check_nrow_candidate(3, msa, destinations, configuration, 4, &expected,
                       &oldState, label);
}

static void test_add_case(const int ncurOld, const int *configuration,
                          const unsigned int mask, const int rsa,
                          const int rsb) {
  int candidateConfiguration[MAX_DIMENSION];
  PfState oldState;
  PfState expected;
  double complex candidate[MAX_QP];
  char label[160];
  memcpy(candidateConfiguration, configuration,
         (size_t)ncurOld * sizeof(*candidateConfiguration));
  snprintf(label, sizeof(label),
           "add M=%d n=%d mask=%u rs=%d,%d", Nsite2, ncurOld, mask, rsa,
           rsb);
  initialize_old_state(ncurOld, candidateConfiguration, &oldState, label);
  candidateConfiguration[ncurOld] = rsa;
  candidateConfiguration[ncurOld + 1] = rsb;
  CalculateNewPfMAddGC(rsa, rsb, candidate, candidateConfiguration, ncurOld,
                       0, MAX_QP);
  check_candidate_noncommit(&oldState, label);
  rebuild_candidate(ncurOld + 2, candidateConfiguration, &expected, label);
  check_pf_values(candidate, &expected, label);
  restore_state(&oldState);
  UpdateMAllAddGC(rsa, rsb, candidateConfiguration, ncurOld, 0, MAX_QP);
  check_updated_state(ncurOld + 2, &expected, &oldState, label);
  if (ncurOld == 0) {
    int qp;
    for (qp = 0; qp < MAX_QP; qp++) {
      const double complex *slater =
          SlaterElm + (size_t)qp * (size_t)Nsite2 * (size_t)Nsite2;
      const double complex base =
          -slater[(size_t)rsa * (size_t)Nsite2 + (size_t)rsb];
      CHECK(close_complex(expected.pf[qp], base),
            "%s qp=%d vacuum base Pf != X_ab", label, qp);
    }
  }
}

static void test_remove_case(const int ncurOld, const int *configuration,
                             const unsigned int mask, const int pos0,
                             const int pos1) {
  int candidateConfiguration[MAX_DIMENSION];
  int eleCfg[MAX_DIMENSION];
  int eleNum[MAX_DIMENSION];
  int ncur = ncurOld;
  PfState oldState;
  PfState expected;
  double complex candidate[MAX_QP];
  char label[176];
  int i;
  memcpy(candidateConfiguration, configuration,
         (size_t)ncurOld * sizeof(*candidateConfiguration));
  for (i = 0; i < Nsite2; i++) {
    eleCfg[i] = -1;
    eleNum[i] = 0;
  }
  for (i = 0; i < ncurOld; i++) {
    eleCfg[candidateConfiguration[i]] = i;
    eleNum[candidateConfiguration[i]] = 1;
  }
  snprintf(label, sizeof(label),
           "remove M=%d n=%d mask=%u pos=%d,%d", Nsite2, ncurOld, mask,
           pos0, pos1);
  initialize_old_state(ncurOld, candidateConfiguration, &oldState, label);
  (void)GCRemovePair(pos0, pos1, candidateConfiguration, eleCfg, eleNum,
                     &ncur);
  CHECK(ncur == ncurOld - 2, "%s mutation count", label);
  CalculateNewPfMRemoveGC(pos0, pos1, candidate, candidateConfiguration,
                          ncurOld, 0, MAX_QP);
  check_candidate_noncommit(&oldState, label);
  rebuild_candidate(ncur, candidateConfiguration, &expected, label);
  check_pf_values(candidate, &expected, label);
  restore_state(&oldState);
  UpdateMAllRemoveGC(pos0, pos1, candidateConfiguration, ncurOld, 0,
                     MAX_QP);
  check_updated_state(ncur, &expected, &oldState, label);
}

static void test_add_remove_exhaustive(void) {
  unsigned int mask;
  setup_dimension(6);
  for (mask = 0; mask < (1U << 6); mask++) {
    const int ncur = popcount(mask);
    int configuration[MAX_DIMENSION];
    int count = 0;
    int rs;
    for (rs = 0; rs < 6; rs++) {
      if ((mask & (1U << rs)) != 0U) configuration[count++] = rs;
    }
    if (ncur == 0 || ncur == 2 || ncur == 4) {
      int rsa;
      for (rsa = 0; rsa < 6; rsa++) {
        int rsb;
        if ((mask & (1U << rsa)) != 0U) continue;
        for (rsb = rsa + 1; rsb < 6; rsb++) {
          if ((mask & (1U << rsb)) == 0U) {
            test_add_case(ncur, configuration, mask, rsa, rsb);
          }
        }
      }
    }
    if (ncur == 2 || ncur == 4 || ncur == 6) {
      int pos0;
      for (pos0 = 0; pos0 < ncur; pos0++) {
        int pos1;
        for (pos1 = pos0 + 1; pos1 < ncur; pos1++) {
          test_remove_case(ncur, configuration, mask, pos0, pos1);
        }
      }
    }
  }
}

static void test_add_conjugate_transpose_mutation(void) {
  int configuration[MAX_DIMENSION] = {0, 1, 2, 3};
  PfState oldState;
  PfState expected;
  double complex y0[2];
  double complex y1[2];
  double complex b0[2];
  double complex b1[2];
  double complex d01;
  double complex wrong01;
  const double complex *slater;
  const char *label = "add ordinary-transpose mutation";
  int i;
  int j;
  setup_dimension(6);
  initialize_old_state(2, configuration, &oldState, label);
  configuration[2] = 2;
  configuration[3] = 3;
  rebuild_candidate(4, configuration, &expected, label);
  restore_state(&oldState);
  slater = SlaterElm;
  for (i = 0; i < 2; i++) {
    b0[i] = slater[(size_t)2 * (size_t)Nsite2 +
                    (size_t)configuration[i]];
    b1[i] = slater[(size_t)3 * (size_t)Nsite2 +
                    (size_t)configuration[i]];
    y0[i] = 0.0;
    y1[i] = 0.0;
    for (j = 0; j < 2; j++) {
      y0[i] += InvM[(size_t)i * (size_t)NsizeMax + (size_t)j] * b0[j];
      y1[i] += InvM[(size_t)i * (size_t)NsizeMax + (size_t)j] * b1[j];
    }
  }
  d01 = -slater[(size_t)2 * (size_t)Nsite2 + 3U];
  for (i = 0; i < 2; i++) d01 += b0[i] * y1[i];
  wrong01 = InvM[1] +
            (-y0[0] * conj(y1[1]) + y1[0] * conj(y0[1])) / d01;
  CHECK(cabs(wrong01 - expected.inv[1]) > 1.0e-5,
        "complex fixture does not kill conjugate-transpose mutation");
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  NThread = 1;
  initializeWorkSpaceAll();
  test_m6_exhaustive();
  test_m8_noncontracted_n3();
  test_add_remove_exhaustive();
  test_add_conjugate_transpose_mutation();
  FreeWorkSpaceAll();
  free(SlaterElm);
  free(InvM);
  free(PfM);
  SlaterElm = NULL;
  InvM = NULL;
  PfM = NULL;
#ifdef _mpi_use
  MPI_Finalize();
#endif
  if (failures != 0) {
    fprintf(stderr, "GCPfUpdate_Unit: %d failure(s)\n", failures);
    return EXIT_FAILURE;
  }
  printf("GCPfUpdate_Unit: PASS\n");
  return EXIT_SUCCESS;
}
