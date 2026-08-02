#include "classic_pfaffian_state.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int Nsite;
extern int Ne;
extern int Nsize;
extern int Nsite2;
extern int NQPFull;
extern int NThread;
extern int LapackLWork;
extern double *SlaterElm_real;
extern double *InvM_real;
extern double *PfM_real;
extern double complex *SlaterElm;
extern double complex *InvM;
extern double complex *PfM;

extern int omp_get_max_threads(void);
extern void initializeWorkSpaceAll(void);
extern void FreeWorkSpaceAll(void);
extern int CalculateMAll_real(const int *, int, int);
extern int CalculateMAll_fcmp(const int *, int, int);
extern void CalculateNewPfM2_real(int, int, double *, const int *, int, int);
extern void CalculateNewPfM2(int, int, double complex *, const int *, int, int);
extern void CalculateNewPfMTwo2_real(
    int, int, int, int, double *, const int *, int, int);
extern void CalculateNewPfMTwo2_fcmp(
    int, int, int, int, double complex *, const int *, int, int);

static int failure_count = 0;

#define CHECK(condition, message)                                      \
  do {                                                                 \
    if (!(condition)) {                                                \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__);   \
      ++failure_count;                                                 \
    }                                                                  \
  } while (0)

static int close_complex(double complex actual, double complex expected) {
  return cabs(actual - expected) <= 4.0e-12 * fmax(1.0, cabs(expected));
}

static void set_real_pair(double *matrix, int row, int column, double value) {
  matrix[(size_t)row * (size_t)Nsite2 + (size_t)column] = value;
  matrix[(size_t)column * (size_t)Nsite2 + (size_t)row] = -value;
}

static void set_complex_pair(double complex *matrix, int row, int column,
                             double complex value) {
  matrix[(size_t)row * (size_t)Nsite2 + (size_t)column] = value;
  matrix[(size_t)column * (size_t)Nsite2 + (size_t)row] = -value;
}

static void fill_slater(int qp_total, int nsize) {
  int qp;
  memset(SlaterElm_real, 0,
         (size_t)qp_total * (size_t)Nsite2 * (size_t)Nsite2 *
             sizeof(*SlaterElm_real));
  memset(SlaterElm, 0,
         (size_t)qp_total * (size_t)Nsite2 * (size_t)Nsite2 *
             sizeof(*SlaterElm));
  for (qp = 0; qp < qp_total; ++qp) {
    const double real_value = 0.75 + (double)qp;
    const double complex complex_value =
        (0.5 + 0.25 * I) * (1.0 + (double)qp);
    double *real_matrix =
        SlaterElm_real + (size_t)qp * (size_t)Nsite2 * (size_t)Nsite2;
    double complex *complex_matrix =
        SlaterElm + (size_t)qp * (size_t)Nsite2 * (size_t)Nsite2;

    set_real_pair(real_matrix, 0, 2, real_value);
    set_complex_pair(complex_matrix, 0, 2, complex_value);
    if (nsize == 4) {
      set_real_pair(real_matrix, 0, 1, 1.0 * real_value);
      set_real_pair(real_matrix, 0, 3, 0.5 * real_value);
      set_real_pair(real_matrix, 1, 2, 1.5 * real_value);
      set_real_pair(real_matrix, 1, 3, -0.75 * real_value);
      set_real_pair(real_matrix, 2, 3, 2.25 * real_value);
      set_complex_pair(complex_matrix, 0, 1, 1.0 * complex_value);
      set_complex_pair(complex_matrix, 0, 3, 0.5 * complex_value);
      set_complex_pair(complex_matrix, 1, 2, 1.5 * complex_value);
      set_complex_pair(complex_matrix, 1, 3, -0.75 * complex_value);
      set_complex_pair(complex_matrix, 2, 3, 2.25 * complex_value);
    }
  }
}

static void run_ab(int qp_total, int nsize, const int *ele_idx) {
  const double complex real_weights[4] = {1.0, 1.0, 1.0, 1.0};
  const double complex complex_weights[4] = {1.0, 1.0, 1.0, 1.0};
  MVMCClassicPfaffianRealWorkspace *real_workspace = NULL;
  MVMCClassicPfaffianComplexWorkspace *complex_workspace = NULL;
  double *real_mirror;
  double complex *complex_mirror;
  size_t matrix_elements = (size_t)nsize * (size_t)nsize;
  size_t matrix_count = (size_t)qp_total * matrix_elements;
  size_t mirror_count = (size_t)qp_total * (matrix_elements + 1);
  size_t slater_count =
      (size_t)qp_total * (size_t)Nsite2 * (size_t)Nsite2;
  size_t index;

  NQPFull = qp_total;
  SlaterElm_real = (double *)calloc(slater_count, sizeof(*SlaterElm_real));
  InvM_real = (double *)calloc(matrix_count, sizeof(*InvM_real));
  PfM_real = (double *)calloc((size_t)qp_total, sizeof(*PfM_real));
  SlaterElm = (double complex *)calloc(slater_count, sizeof(*SlaterElm));
  InvM = (double complex *)calloc(matrix_count, sizeof(*InvM));
  PfM = (double complex *)calloc((size_t)qp_total, sizeof(*PfM));
  real_mirror = (double *)calloc(mirror_count, sizeof(*real_mirror));
  complex_mirror = (double complex *)calloc(mirror_count,
                                             sizeof(*complex_mirror));
  CHECK(SlaterElm_real != NULL && InvM_real != NULL && PfM_real != NULL &&
            SlaterElm != NULL && InvM != NULL && PfM != NULL &&
            real_mirror != NULL && complex_mirror != NULL,
        "allocate legacy A/B fixture");
  if (SlaterElm_real == NULL || InvM_real == NULL || PfM_real == NULL ||
      SlaterElm == NULL || InvM == NULL || PfM == NULL ||
      real_mirror == NULL || complex_mirror == NULL) {
    goto cleanup;
  }
  fill_slater(qp_total, nsize);
  CHECK(CalculateMAll_real(ele_idx, 0, qp_total) == 0,
        "execute legacy CalculateMAll_real oracle");
  CHECK(CalculateMAll_fcmp(ele_idx, 0, qp_total) == 0,
        "execute legacy CalculateMAll_fcmp oracle");

  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, nsize, qp_total, 0, qp_total, &real_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create new real A/B workspace");
  CHECK(mvmc_classic_pfaffian_real_prepare(
            real_workspace, SlaterElm_real, ele_idx, real_weights,
            0.0, 0.0) == MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_classic_pfaffian_real_publish(
                real_workspace, real_mirror, mirror_count) ==
                MVMC_PFAFFIAN_STATUS_OK,
        "execute new real prepare/publish");
  CHECK(mvmc_classic_pfaffian_complex_workspace_create(
            2, nsize, qp_total, 0, qp_total, &complex_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create new complex A/B workspace");
  CHECK(mvmc_classic_pfaffian_complex_prepare(
            complex_workspace, SlaterElm, ele_idx, complex_weights,
            0.0, 0.0) == MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_classic_pfaffian_complex_publish(
                complex_workspace, complex_mirror, mirror_count) ==
                MVMC_PFAFFIAN_STATUS_OK,
        "execute new complex prepare/publish");

  for (index = 0; index < matrix_count; ++index) {
    CHECK(fabs(real_mirror[index] - InvM_real[index]) <=
              4.0e-12 * fmax(1.0, fabs(InvM_real[index])),
          "real new inverse mirror matches executed legacy wrapper");
    CHECK(close_complex(complex_mirror[index], InvM[index]),
          "complex new inverse mirror matches executed legacy wrapper");
  }
  for (index = 0; index < (size_t)qp_total; ++index) {
    CHECK(fabs(real_mirror[matrix_count + index] - PfM_real[index]) <=
              4.0e-12 * fmax(1.0, fabs(PfM_real[index])),
          "real new PfM mirror matches executed legacy wrapper");
    CHECK(close_complex(complex_mirror[matrix_count + index], PfM[index]),
          "complex new PfM mirror matches executed legacy wrapper");
  }

cleanup:
  mvmc_classic_pfaffian_real_workspace_destroy(real_workspace);
  mvmc_classic_pfaffian_complex_workspace_destroy(complex_workspace);
  free(complex_mirror);
  free(real_mirror);
  free(PfM);
  free(InvM);
  free(SlaterElm);
  free(PfM_real);
  free(InvM_real);
  free(SlaterElm_real);
  PfM = NULL;
  InvM = NULL;
  SlaterElm = NULL;
  PfM_real = NULL;
  InvM_real = NULL;
  SlaterElm_real = NULL;
}

static void run_measurement_singular_blocker_negative_control(void) {
  const int ele_idx[2] = {0, 0};
  const double complex weights[2] = {1.0, 1.0};
  const int qp_total = 2;
  const size_t matrix_elements = 4;
  const size_t matrix_count = (size_t)qp_total * matrix_elements;
  const size_t mirror_count = (size_t)qp_total * (matrix_elements + 1);
  const size_t slater_count =
      (size_t)qp_total * (size_t)Nsite2 * (size_t)Nsite2;
  MVMCClassicPfaffianRealWorkspace *real_workspace = NULL;
  MVMCClassicPfaffianComplexWorkspace *complex_workspace = NULL;
  const MVMCClassicPfaffianState *real_state;
  const MVMCClassicPfaffianState *complex_state;
  double *real_mirror = NULL;
  double complex *complex_mirror = NULL;

  NQPFull = qp_total;
  SlaterElm_real = (double *)calloc(slater_count, sizeof(*SlaterElm_real));
  InvM_real = (double *)calloc(matrix_count, sizeof(*InvM_real));
  PfM_real = (double *)calloc((size_t)qp_total, sizeof(*PfM_real));
  SlaterElm = (double complex *)calloc(slater_count, sizeof(*SlaterElm));
  InvM = (double complex *)calloc(matrix_count, sizeof(*InvM));
  PfM = (double complex *)calloc((size_t)qp_total, sizeof(*PfM));
  real_mirror = (double *)calloc(mirror_count, sizeof(*real_mirror));
  complex_mirror = (double complex *)calloc(mirror_count,
                                             sizeof(*complex_mirror));
  CHECK(SlaterElm_real != NULL && InvM_real != NULL && PfM_real != NULL &&
            SlaterElm != NULL && InvM != NULL && PfM != NULL &&
            real_mirror != NULL && complex_mirror != NULL,
        "allocate measurement singular blocker fixture");
  if (SlaterElm_real == NULL || InvM_real == NULL || PfM_real == NULL ||
      SlaterElm == NULL || InvM == NULL || PfM == NULL ||
      real_mirror == NULL || complex_mirror == NULL) {
    goto cleanup;
  }

  /* QP 0 is regular and QP 1 is exactly singular, but the total is nonzero. */
  set_real_pair(SlaterElm_real, 0, 2, 1.0);
  set_complex_pair(SlaterElm, 0, 2, 1.0 + 0.25 * I);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, qp_total, 0, qp_total, &real_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_classic_pfaffian_real_prepare(
                real_workspace, SlaterElm_real, ele_idx, weights,
                0.0, 0.0) == MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_classic_pfaffian_real_publish(
                real_workspace, real_mirror, mirror_count) ==
                MVMC_PFAFFIAN_STATUS_OK,
        "absolute real state keeps partial-singular nonzero-total sample");
  real_state = mvmc_classic_pfaffian_real_accepted(real_workspace);
  CHECK(real_state != NULL && real_state->valid &&
            real_state->local_aggregate.singular_count == 1 &&
            cabs(real_state->local_aggregate.total) > 0.0,
        "real partial-singular sample remains a valid projected state");
  CHECK(mvmc_classic_pfaffian_complex_workspace_create(
            2, 2, qp_total, 0, qp_total, &complex_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_classic_pfaffian_complex_prepare(
                complex_workspace, SlaterElm, ele_idx, weights,
                0.0, 0.0) == MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_classic_pfaffian_complex_publish(
                complex_workspace, complex_mirror, mirror_count) ==
                MVMC_PFAFFIAN_STATUS_OK,
        "absolute complex state keeps partial-singular nonzero-total sample");
  complex_state = mvmc_classic_pfaffian_complex_accepted(complex_workspace);
  CHECK(complex_state != NULL && complex_state->valid &&
            complex_state->local_aggregate.singular_count == 1 &&
            cabs(complex_state->local_aggregate.total) > 0.0,
        "complex partial-singular sample remains a valid projected state");

  CHECK(CalculateMAll_real(ele_idx, 0, qp_total) != 0,
        "legacy real measurement rebuild rejects the valid singular sample");
  CHECK(CalculateMAll_fcmp(ele_idx, 0, qp_total) != 0,
        "legacy complex measurement rebuild rejects the valid singular sample");

cleanup:
  mvmc_classic_pfaffian_real_workspace_destroy(real_workspace);
  mvmc_classic_pfaffian_complex_workspace_destroy(complex_workspace);
  free(complex_mirror);
  free(real_mirror);
  free(PfM);
  free(InvM);
  free(SlaterElm);
  free(PfM_real);
  free(InvM_real);
  free(SlaterElm_real);
  PfM = NULL;
  InvM = NULL;
  SlaterElm = NULL;
  PfM_real = NULL;
  InvM_real = NULL;
  SlaterElm_real = NULL;
}

static void run_actual_fast_absolute_observation(void) {
  const int initial_idx[2] = {0, 0};
  const int hopping_idx[2] = {1, 0};
  const int exchange_idx[2] = {0, 1};
  const double complex weights[4] = {1.0, 1.0, 1.0, 1.0};
  const int qp_total = 4;
  const size_t matrix_count = (size_t)qp_total * 4;
  const size_t slater_count =
      (size_t)qp_total * (size_t)Nsite2 * (size_t)Nsite2;
  double fast_rank_one_real[4];
  double fast_rank_two_real[4];
  double complex fast_rank_one_complex[4];
  double complex fast_rank_two_complex[4];
  MVMCClassicPfaffianRealWorkspace *real_workspace = NULL;
  MVMCClassicPfaffianComplexWorkspace *complex_workspace = NULL;
  const MVMCClassicPfaffianState *candidate;
  int qp;

  NQPFull = qp_total;
  SlaterElm_real = (double *)calloc(slater_count, sizeof(*SlaterElm_real));
  InvM_real = (double *)calloc(matrix_count, sizeof(*InvM_real));
  PfM_real = (double *)calloc((size_t)qp_total, sizeof(*PfM_real));
  SlaterElm = (double complex *)calloc(slater_count, sizeof(*SlaterElm));
  InvM = (double complex *)calloc(matrix_count, sizeof(*InvM));
  PfM = (double complex *)calloc((size_t)qp_total, sizeof(*PfM));
  CHECK(SlaterElm_real != NULL && InvM_real != NULL && PfM_real != NULL &&
            SlaterElm != NULL && InvM != NULL && PfM != NULL,
        "allocate actual fast/absolute observation fixture");
  if (SlaterElm_real == NULL || InvM_real == NULL || PfM_real == NULL ||
      SlaterElm == NULL || InvM == NULL || PfM == NULL) {
    goto cleanup;
  }
  fill_slater(qp_total, 2);
  CHECK(CalculateMAll_real(initial_idx, 0, qp_total) == 0 &&
            CalculateMAll_fcmp(initial_idx, 0, qp_total) == 0,
        "initialize actual legacy fast proposal cache");
  CalculateNewPfM2_real(0, 0, fast_rank_one_real, hopping_idx, 0, qp_total);
  CalculateNewPfM2(0, 0, fast_rank_one_complex, hopping_idx, 0, qp_total);
  CalculateNewPfMTwo2_real(
      0, 0, 0, 1, fast_rank_two_real, exchange_idx, 0, qp_total);
  CalculateNewPfMTwo2_fcmp(
      0, 0, 0, 1, fast_rank_two_complex, exchange_idx, 0, qp_total);

  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, qp_total, 0, qp_total, &real_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_classic_pfaffian_complex_workspace_create(
                2, 2, qp_total, 0, qp_total, &complex_workspace) ==
                MVMC_PFAFFIAN_STATUS_OK,
        "create actual fast/absolute observation workspaces");
  CHECK(mvmc_classic_pfaffian_real_prepare(
            real_workspace, SlaterElm_real, hopping_idx, weights,
            0.0, 0.0) == MVMC_PFAFFIAN_STATUS_OK,
        "prepare real rank-one absolute candidate");
  candidate = mvmc_classic_pfaffian_real_candidate(real_workspace);
  for (qp = 0; qp < qp_total; ++qp) {
    CHECK(candidate != NULL &&
              fabs(fast_rank_one_real[qp] -
                   creal(candidate->components[qp].pfaffian)) <=
                  1.0e-11 * fmax(1.0, cabs(candidate->components[qp].pfaffian)),
          "actual real rank-one fast Pfaffian matches absolute candidate");
  }
  mvmc_classic_pfaffian_real_discard_candidate(real_workspace);
  CHECK(mvmc_classic_pfaffian_complex_prepare(
            complex_workspace, SlaterElm, hopping_idx, weights,
            0.0, 0.0) == MVMC_PFAFFIAN_STATUS_OK,
        "prepare complex rank-one absolute candidate");
  candidate = mvmc_classic_pfaffian_complex_candidate(complex_workspace);
  for (qp = 0; qp < qp_total; ++qp) {
    CHECK(candidate != NULL &&
              close_complex(fast_rank_one_complex[qp],
                            candidate->components[qp].pfaffian),
          "actual complex rank-one fast Pfaffian matches absolute candidate");
  }
  mvmc_classic_pfaffian_complex_discard_candidate(complex_workspace);

  CHECK(mvmc_classic_pfaffian_real_prepare(
            real_workspace, SlaterElm_real, exchange_idx, weights,
            0.0, 0.0) == MVMC_PFAFFIAN_STATUS_OK,
        "prepare real rank-two singular absolute candidate");
  candidate = mvmc_classic_pfaffian_real_candidate(real_workspace);
  CHECK(candidate != NULL && candidate->components[0].state ==
            MVMC_PFAFFIAN_SINGULAR,
        "rank-two fast observation includes a singular absolute component");
  for (qp = 0; qp < qp_total; ++qp) {
    CHECK(candidate != NULL &&
              fabs(fast_rank_two_real[qp] -
                   creal(candidate->components[qp].pfaffian)) <=
                  1.0e-11 * fmax(1.0, cabs(candidate->components[qp].pfaffian)),
          "actual real rank-two fast Pfaffian compares singular without skip");
  }
  mvmc_classic_pfaffian_real_discard_candidate(real_workspace);
  CHECK(mvmc_classic_pfaffian_complex_prepare(
            complex_workspace, SlaterElm, exchange_idx, weights,
            0.0, 0.0) == MVMC_PFAFFIAN_STATUS_OK,
        "prepare complex rank-two singular absolute candidate");
  candidate = mvmc_classic_pfaffian_complex_candidate(complex_workspace);
  CHECK(candidate != NULL && candidate->components[0].state ==
            MVMC_PFAFFIAN_SINGULAR,
        "complex rank-two observation includes singular component");
  for (qp = 0; qp < qp_total; ++qp) {
    CHECK(candidate != NULL &&
              close_complex(fast_rank_two_complex[qp],
                            candidate->components[qp].pfaffian),
          "actual complex rank-two fast Pfaffian compares singular without skip");
  }

cleanup:
  mvmc_classic_pfaffian_real_workspace_destroy(real_workspace);
  mvmc_classic_pfaffian_complex_workspace_destroy(complex_workspace);
  free(PfM);
  free(InvM);
  free(SlaterElm);
  free(PfM_real);
  free(InvM_real);
  free(SlaterElm_real);
  PfM = NULL;
  InvM = NULL;
  SlaterElm = NULL;
  PfM_real = NULL;
  InvM_real = NULL;
  SlaterElm_real = NULL;
}

int main(void) {
  const int ele_idx_n2[2] = {0, 0};
  const int ele_idx_n4[4] = {1, 0, 0, 1};

  Nsite = 2;
  Nsite2 = 4;
  NThread = omp_get_max_threads();

  Ne = 1;
  Nsize = 2;
  LapackLWork = 128;
  initializeWorkSpaceAll();
  run_ab(1, Nsize, ele_idx_n2);
  run_ab(4, Nsize, ele_idx_n2);
  run_measurement_singular_blocker_negative_control();
  run_actual_fast_absolute_observation();
  FreeWorkSpaceAll();

  Ne = 2;
  Nsize = 4;
  LapackLWork = 128;
  initializeWorkSpaceAll();
  run_ab(1, Nsize, ele_idx_n4);
  run_ab(4, Nsize, ele_idx_n4);
  FreeWorkSpaceAll();
  if (failure_count == 0) {
    printf("classic Pfaffian legacy execution A/B unit: PASS\n");
  }
  return failure_count == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
