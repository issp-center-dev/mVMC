#include "classic_krylov_amplitude.h"
#include "classic_pfaffian_matrix.h"

#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;

#define CHECK(condition, message)                                    \
  do {                                                               \
    if (!(condition)) {                                              \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__); \
      ++failure_count;                                               \
    }                                                                \
  } while (0)

static int close_complex(double complex actual, double complex expected,
                         double tolerance) {
  return cabs(actual - expected) <= tolerance * fmax(1.0, cabs(expected));
}

static void rank_and_size(int *rank, int *size) {
#ifdef _mpi_use
  MPI_Comm_rank(MPI_COMM_WORLD, rank);
  MPI_Comm_size(MPI_COMM_WORLD, size);
#else
  *rank = 0;
  *size = 1;
#endif
}

static MVMCClassicPfaffianCommunicator world_communicator(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static MVMCClassicKrylovAmplitudeLayout base_layout(int qp_total,
                                                     int rank, int size) {
  MVMCClassicKrylovAmplitudeLayout layout;
  memset(&layout, 0, sizeof(layout));
  layout.site_count = 2;
  layout.up_electron_count = 1;
  layout.down_electron_count = 1;
  layout.qp_total = qp_total;
  layout.qp_start = qp_total * rank / size;
  layout.qp_end = qp_total * (rank + 1) / size;
  layout.communicator = world_communicator();
  return layout;
}

static void set_real_pair(double *slater, int nsite, int up, int down,
                          double pfaffian) {
  const int nsite2 = 2 * nsite;
  const int down_orbital = nsite + down;
  slater[(size_t)down_orbital * (size_t)nsite2 + (size_t)up] = -pfaffian;
  slater[(size_t)up * (size_t)nsite2 + (size_t)down_orbital] = pfaffian;
}

static void set_complex_pair(double complex *slater, int nsite, int up,
                             int down, double complex pfaffian) {
  const int nsite2 = 2 * nsite;
  const int down_orbital = nsite + down;
  slater[(size_t)down_orbital * (size_t)nsite2 + (size_t)up] = -pfaffian;
  slater[(size_t)up * (size_t)nsite2 + (size_t)down_orbital] = pfaffian;
}

static void fill_real_slater(double *slater, int qp_total, int singular_qp,
                             double *pfaffians) {
  int qp, up, down;
  memset(slater, 0, (size_t)qp_total * 16 * sizeof(*slater));
  for (qp = 0; qp < qp_total; ++qp) {
    double *component = slater + (size_t)qp * 16;
    const double multiplier = qp == singular_qp ? 0.0 : (double)(qp + 1);
    for (up = 0; up < 2; ++up) {
      for (down = 0; down < 2; ++down) {
        set_real_pair(component, 2, up, down,
                      multiplier * (1.5 + (double)up + 0.25 * (double)down));
      }
    }
    pfaffians[qp] = multiplier * 1.5;
  }
}

static void fill_complex_slater(double complex *slater, int qp_total,
                                int singular_qp,
                                double complex *pfaffians) {
  int qp, up, down;
  memset(slater, 0, (size_t)qp_total * 16 * sizeof(*slater));
  for (qp = 0; qp < qp_total; ++qp) {
    double complex *component = slater + (size_t)qp * 16;
    const double complex multiplier =
        qp == singular_qp ? 0.0 : (double)(qp + 1) + 0.2 * (double)qp * I;
    for (up = 0; up < 2; ++up) {
      for (down = 0; down < 2; ++down) {
        set_complex_pair(
            component, 2, up, down,
            multiplier * (1.5 + (double)up + 0.25 * (double)down));
      }
    }
    pfaffians[qp] = multiplier * 1.5;
  }
}

static void enable_gutzwiller(MVMCClassicKrylovAmplitudeLayout *layout,
                              const int *indices,
                              const double complex *parameter) {
  layout->nproj = 1;
  layout->ngutzwiller_idx = 1;
  layout->gutzwiller_idx = indices;
  layout->projection_parameters = parameter;
}

static void test_matrix_helper(void) {
  double real_slater[16] = {0.0};
  double real_original[16];
  double real_matrix[4] = {0.0};
  double complex complex_slater[16] = {0.0};
  double complex complex_original[16];
  double complex complex_matrix[4] = {0.0};
  const int ele_idx[2] = {1, 0};
  int invalid_idx[2] = {2, 0};

  set_real_pair(real_slater, 2, 1, 0, 3.25);
  set_complex_pair(complex_slater, 2, 1, 0, -0.5 + 2.0 * I);
  memcpy(real_original, real_slater, sizeof(real_slater));
  memcpy(complex_original, complex_slater, sizeof(complex_slater));
  CHECK(mvmc_classic_pfaffian_build_real_matrix(
            real_slater, 2, 2, ele_idx, real_matrix) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "real matrix helper status");
  CHECK(real_matrix[0] == 0.0 && real_matrix[1] == -3.25 &&
            real_matrix[2] == 3.25 && real_matrix[3] == 0.0,
        "real matrix helper frozen column-major oracle");
  CHECK(memcmp(real_slater, real_original, sizeof(real_slater)) == 0,
        "real matrix helper leaves source unchanged");
  CHECK(mvmc_classic_pfaffian_build_complex_matrix(
            complex_slater, 2, 2, ele_idx, complex_matrix) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "complex matrix helper status");
  CHECK(complex_matrix[0] == 0.0 &&
            complex_matrix[1] == 0.5 - 2.0 * I &&
            complex_matrix[2] == -0.5 + 2.0 * I &&
            complex_matrix[3] == 0.0,
        "complex matrix helper frozen column-major oracle");
  CHECK(memcmp(complex_slater, complex_original,
               sizeof(complex_slater)) == 0,
        "complex matrix helper leaves source unchanged");
  CHECK(mvmc_classic_pfaffian_build_real_matrix(
            real_slater, 0, 2, ele_idx, real_matrix) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "matrix helper rejects zero site count");
  CHECK(mvmc_classic_pfaffian_build_real_matrix(
            real_slater, INT_MAX, 2, ele_idx, real_matrix) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            mvmc_classic_pfaffian_build_complex_matrix(
                complex_slater, INT_MAX, 2, ele_idx, complex_matrix) ==
                MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "matrix helpers reject overflow dimensions before arithmetic");
  CHECK(mvmc_classic_pfaffian_build_real_matrix(
            real_slater, 2, 3, ele_idx, real_matrix) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "matrix helper rejects odd particle count");
  CHECK(mvmc_classic_pfaffian_build_real_matrix(
            real_slater, 2, 2, invalid_idx, real_matrix) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "matrix helper rejects out-of-range electron index");
}

static void test_real_amplitude(int qp_total, int use_gutzwiller,
                                int singular_qp, int mutate_source) {
  double slater[64];
  double pfaffians[4] = {0.0};
  double complex weights[4] = {1.0, -1.0, 1.0, -1.0};
  double complex parameter = -0.43;
  int gutz_indices[2] = {0, 0};
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  MVMCKrylovAmplitudeResult result;
  const uint64_t configuration = UINT64_C(5);
  double complex expected = 0.0;
  int rank, size, qp;

  rank_and_size(&rank, &size);
  layout = base_layout(qp_total, rank, size);
  if (use_gutzwiller) {
    enable_gutzwiller(&layout, gutz_indices, &parameter);
  }
  fill_real_slater(slater, qp_total, singular_qp, pfaffians);
  for (qp = 0; qp < qp_total; ++qp) expected += weights[qp] * pfaffians[qp];
  if (use_gutzwiller) expected *= exp(creal(parameter));
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "real amplitude workspace create");
  if (workspace == NULL) return;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_bytes(workspace) >
            sizeof(*pfaffians),
        "real amplitude workspace reports preallocated bytes");
  if (mutate_source) {
    memset(slater, 0x5a, sizeof(slater));
    for (qp = 0; qp < qp_total; ++qp) weights[qp] = 99.0;
    parameter = 99.0;
    gutz_indices[0] = gutz_indices[1] = 1;
  }
  CHECK(mvmc_classic_krylov_real_amplitude(
            &configuration, 1, workspace, &result) ==
            MVMC_KRYLOV_STATUS_OK,
        "real amplitude evaluation status");
  CHECK(close_complex(result.value, expected, 2.0e-14),
        "real amplitude matches direct QP/projection oracle");
  CHECK(result.total_zero == (expected == 0.0),
        "real amplitude exact zero classification");
  CHECK(result.global_factorization_count == (uint64_t)qp_total &&
            result.local_factorization_count ==
                (uint64_t)(layout.qp_end - layout.qp_start),
        "real amplitude local/global factorization counts");
  CHECK(result.singular_component_count ==
            (uint64_t)(singular_qp >= 0 ? 1 : 0) &&
            result.regular_component_count +
                    result.near_pivot_component_count +
                    result.singular_component_count ==
                (uint64_t)qp_total,
        "real amplitude component classifications");
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);
}

static void test_complex_amplitude(int qp_total, int use_gutzwiller,
                                   int singular_qp) {
  double complex slater[64];
  double complex pfaffians[4] = {0.0};
  double complex weights[4] = {1.0, -1.0 + 0.25 * I,
                               0.5 - 0.75 * I, -I};
  double complex parameter = -0.43;
  int gutz_indices[2] = {0, 0};
  MVMCClassicKrylovComplexAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  MVMCKrylovAmplitudeResult result;
  const uint64_t configuration = UINT64_C(5);
  double complex expected = 0.0;
  int rank, size, qp;

  rank_and_size(&rank, &size);
  layout = base_layout(qp_total, rank, size);
  if (use_gutzwiller) {
    enable_gutzwiller(&layout, gutz_indices, &parameter);
  }
  fill_complex_slater(slater, qp_total, singular_qp, pfaffians);
  for (qp = 0; qp < qp_total; ++qp) expected += weights[qp] * pfaffians[qp];
  if (use_gutzwiller) expected *= exp(creal(parameter));
  CHECK(mvmc_classic_krylov_complex_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "complex amplitude workspace create");
  if (workspace == NULL) return;
  CHECK(mvmc_classic_krylov_complex_amplitude_workspace_bytes(workspace) >
            sizeof(*pfaffians),
        "complex amplitude workspace reports preallocated bytes");
  if (qp_total == 4 && use_gutzwiller) {
    memset(slater, 0x5a, sizeof(slater));
    for (qp = 0; qp < qp_total; ++qp) weights[qp] = 99.0 + 17.0 * I;
    parameter = 99.0 + 13.0 * I;
    gutz_indices[0] = gutz_indices[1] = 1;
  }
  CHECK(mvmc_classic_krylov_complex_amplitude(
            &configuration, 1, workspace, &result) ==
            MVMC_KRYLOV_STATUS_OK,
        "complex amplitude evaluation status");
  CHECK(close_complex(result.value, expected, 3.0e-14),
        "complex amplitude matches direct phase oracle");
  CHECK(result.global_factorization_count == (uint64_t)qp_total &&
            result.local_factorization_count ==
                (uint64_t)(layout.qp_end - layout.qp_start),
        "complex amplitude local/global factorization counts");
  CHECK(result.singular_component_count ==
            (uint64_t)(singular_qp >= 0 ? 1 : 0),
        "complex amplitude singular component remains valid");
  mvmc_classic_krylov_complex_amplitude_workspace_destroy(workspace);
}

static void test_all_zero_and_configuration_guards(void) {
  double slater[16] = {0.0};
  double complex weights[1] = {1.0};
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  MVMCKrylovAmplitudeResult result;
  uint64_t configuration = UINT64_C(5);
  int rank, size;

  rank_and_size(&rank, &size);
  layout = base_layout(1, rank, size);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "all-zero amplitude workspace create");
  if (workspace == NULL) return;
  CHECK(mvmc_classic_krylov_real_amplitude(
            &configuration, 1, workspace, &result) ==
            MVMC_KRYLOV_STATUS_OK &&
            result.value == 0.0 && result.total_zero == 1 &&
            result.singular_component_count == 1,
        "all singular components produce a valid exact-zero total");
  configuration = UINT64_C(7);
  CHECK(mvmc_classic_krylov_real_amplitude(
            &configuration, 1, workspace, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            result.value == 0.0,
        "wrong particle sector is rejected atomically");
  configuration = UINT64_C(21);
  CHECK(mvmc_classic_krylov_real_amplitude(
            &configuration, 1, workspace, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            result.value == 0.0,
        "nonzero configuration padding bit is rejected atomically");
  if (size > 1) {
    configuration = rank == size - 1 ? UINT64_C(6) : UINT64_C(5);
    CHECK(mvmc_classic_krylov_real_amplitude(
              &configuration, 1, workspace, &result) ==
              MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE,
          "rank-mismatched configuration fails exact collective preflight");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);
}

static void test_projection_allowlist(void) {
  double slater[16] = {0.0};
  double complex weights[1] = {1.0};
  double complex parameter[10] = {-0.43, 0.25};
  int gutz_indices[2] = {0, 0};
  int jastrow_row_0[2] = {-1, 0};
  int jastrow_row_1[2] = {0, -1};
  const int *jastrow_idx[2] = {jastrow_row_0, jastrow_row_1};
  int dh2_row[4] = {0, 1, 1, 0};
  const int *dh2_idx[1] = {dh2_row};
  int dh4_row[8] = {0, 1, 0, 1, 1, 0, 1, 0};
  const int *dh4_idx[1] = {dh4_row};
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  int rank, size;

  rank_and_size(&rank, &size);
  layout = base_layout(1, rank, size);
  enable_gutzwiller(&layout, gutz_indices, parameter);
  parameter[0] = -0.43 + I;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "real adapter rejects a finite imaginary projection parameter");
  parameter[0] = -0.43;

  layout = base_layout(1, rank, size);
  layout.nproj = 1;
  layout.njastrow_idx = 1;
  layout.jastrow_idx = jastrow_idx;
  layout.projection_parameters = parameter;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL && workspace == NULL,
        "well-formed Jastrow is outside the initial allowlist");

  jastrow_row_0[1] = jastrow_row_1[0] = 1;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "out-of-range unsupported Jastrow index is malformed input");
  jastrow_row_0[1] = jastrow_row_1[0] = 0;

  layout = base_layout(1, rank, size);
  layout.nproj = 2;
  layout.ngutzwiller_idx = 1;
  layout.njastrow_idx = 1;
  gutz_indices[0] = gutz_indices[1] = 1;
  layout.gutzwiller_idx = gutz_indices;
  layout.jastrow_idx = jastrow_idx;
  layout.projection_parameters = parameter;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "Gutzwiller index cannot address a later projection class");
  gutz_indices[0] = gutz_indices[1] = 0;

  layout = base_layout(1, rank, size);
  layout.nproj = 1;
  layout.nspin_jastrow_idx = 1;
  layout.spin_jastrow_idx = jastrow_idx;
  layout.projection_parameters = parameter;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL && workspace == NULL,
        "well-formed Spin Jastrow is outside the initial allowlist");

  layout = base_layout(1, rank, size);
  layout.nproj = 6;
  layout.ndoublon_holon_2site_idx = 1;
  layout.doublon_holon_2site_idx = dh2_idx;
  layout.projection_parameters = parameter;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL && workspace == NULL,
        "well-formed two-site doublon-holon is outside the allowlist");
  dh2_row[0] = -1;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "out-of-range two-site doublon-holon site is malformed input");
  dh2_row[0] = 0;

  layout = base_layout(1, rank, size);
  layout.nproj = 10;
  layout.ndoublon_holon_4site_idx = 1;
  layout.doublon_holon_4site_idx = dh4_idx;
  layout.projection_parameters = parameter;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL && workspace == NULL,
        "well-formed four-site doublon-holon is outside the allowlist");
  dh4_row[7] = 2;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "out-of-range four-site doublon-holon site is malformed input");
  dh4_row[7] = 0;

  layout = base_layout(1, rank, size);
  layout.nproj = 1;
  layout.ngutzwiller_idx = -1;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "negative projection count is malformed input");

  layout = base_layout(1, rank, size);
  layout.nproj = 2;
  layout.ngutzwiller_idx = 1;
  layout.gutzwiller_idx = gutz_indices;
  layout.projection_parameters = parameter;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "projection count formula mismatch is malformed input");

  layout = base_layout(1, rank, size);
  layout.nproj = 2;
  layout.ngutzwiller_idx = 2;
  layout.gutzwiller_idx = gutz_indices;
  layout.projection_parameters = parameter;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
        "multiple Gutzwiller parameters are outside the initial allowlist");

  {
    const uint32_t unsupported_flags[] = {
        MVMC_CLASSIC_KRYLOV_UNSUPPORTED_RBM,
        MVMC_CLASSIC_KRYLOV_UNSUPPORTED_BACKFLOW,
        MVMC_CLASSIC_KRYLOV_UNSUPPORTED_GENERAL_ORBITAL,
        MVMC_CLASSIC_KRYLOV_UNSUPPORTED_BLOCK_UPDATE};
    size_t flag;
    for (flag = 0; flag < sizeof(unsupported_flags) /
                              sizeof(unsupported_flags[0]); ++flag) {
      layout = base_layout(1, rank, size);
      layout.unsupported_features = unsupported_flags[flag];
      CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
                &layout, slater, weights, &workspace) ==
                MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL && workspace == NULL,
            "unsupported feature is rejected without Pfaffian fallback");
    }
  }
}

static void test_weight_binding_guards(void) {
  double real_slater[16] = {0.0};
  double complex complex_slater[16] = {0.0};
  double complex weights[1] = {1.0 + I};
  MVMCClassicKrylovRealAmplitudeWorkspace *real_workspace = NULL;
  MVMCClassicKrylovComplexAmplitudeWorkspace *complex_workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  int rank, size;

  rank_and_size(&rank, &size);
  layout = base_layout(1, rank, size);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, real_slater, weights, &real_workspace) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            real_workspace == NULL,
        "real adapter rejects a finite non-real QP weight");

  weights[0] = NAN + 0.0 * I;
  CHECK(mvmc_classic_krylov_complex_amplitude_workspace_create(
            &layout, complex_slater, weights, &complex_workspace) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            complex_workspace == NULL,
        "complex adapter rejects NaN QP weight real part");
  weights[0] = 1.0 + INFINITY * I;
  CHECK(mvmc_classic_krylov_complex_amplitude_workspace_create(
            &layout, complex_slater, weights, &complex_workspace) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            complex_workspace == NULL,
        "complex adapter rejects infinite QP weight imaginary part");
}

static void test_rank_invariant_cancellation(void) {
  const double terms[4] = {1.0e16, 1.0, -1.0e16, 1.0};
  double slater[64] = {0.0};
  double complex weights[4] = {1.0, 1.0, 1.0, 1.0};
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  MVMCKrylovAmplitudeResult result;
  const uint64_t configuration = UINT64_C(5);
  int rank, size, qp;

  rank_and_size(&rank, &size);
  for (qp = 0; qp < 4; ++qp) {
    set_real_pair(slater + (size_t)qp * 16, 2, 0, 0, terms[qp]);
  }
  layout = base_layout(4, rank, size);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "rank-invariant cancellation workspace create");
  if (workspace != NULL) {
    CHECK(mvmc_classic_krylov_real_amplitude(
              &configuration, 1, workspace, &result) ==
                  MVMC_KRYLOV_STATUS_OK &&
              result.value == 2.0 &&
              result.global_factorization_count == 4,
          "global QP-order Neumaier sum is exactly rank-count invariant");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);
}

static void test_nonfinite_and_near_pivot_classification(void) {
  double slater[64] = {0.0};
  double complex weights[4] = {2.0, 1.0, 1.0, 1.0};
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  MVMCKrylovAmplitudeResult result;
  const uint64_t two_particle_configuration = UINT64_C(5);
  const uint64_t four_particle_configuration = UINT64_C(15);
  int rank, size, qp;

  rank_and_size(&rank, &size);
  set_real_pair(slater, 2, 0, 0, DBL_MAX);
  layout = base_layout(1, rank, size);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "finite product-overflow workspace create");
  if (workspace != NULL) {
    CHECK(mvmc_classic_krylov_real_amplitude(
              &two_particle_configuration, 1, workspace, &result) ==
                  MVMC_KRYLOV_STATUS_NONFINITE &&
              result.value == 0.0,
          "finite weight times finite Pfaffian overflow is NONFINITE");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);

  memset(slater, 0, sizeof(slater));
  weights[0] = weights[1] = weights[2] = weights[3] = 1.0;
  for (qp = 0; qp < 4; ++qp) {
    set_real_pair(slater + (size_t)qp * 16, 2, 0, 0, DBL_MAX);
  }
  workspace = NULL;
  layout = base_layout(4, rank, size);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "finite sum-overflow workspace create");
  if (workspace != NULL) {
    CHECK(mvmc_classic_krylov_real_amplitude(
              &two_particle_configuration, 1, workspace, &result) ==
                  MVMC_KRYLOV_STATUS_NONFINITE &&
              result.value == 0.0,
          "finite projected sum overflow is NONFINITE");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);

  memset(slater, 0, sizeof(slater));
  set_real_pair(slater, 2, 0, 0, NAN);
  workspace = NULL;
  layout = base_layout(1, rank, size);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "nonfinite Slater workspace create");
  if (workspace != NULL) {
    CHECK(mvmc_classic_krylov_real_amplitude(
              &two_particle_configuration, 1, workspace, &result) ==
                  MVMC_KRYLOV_STATUS_NONFINITE &&
              result.value == 0.0,
          "nonfinite Slater entry is classified NONFINITE atomically");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);

  memset(slater, 0, sizeof(slater));
  set_real_pair(slater, 2, 0, 0, 1.0);
  set_real_pair(slater, 2, 1, 1, 1.0e-15);
  workspace = NULL;
  layout = base_layout(1, rank, size);
  layout.up_electron_count = 2;
  layout.down_electron_count = 2;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "near-pivot adapter workspace create");
  if (workspace != NULL) {
    CHECK(mvmc_classic_krylov_real_amplitude(
              &four_particle_configuration, 1, workspace, &result) ==
                  MVMC_KRYLOV_STATUS_OK &&
              result.near_pivot_component_count == 1 &&
              result.regular_component_count == 0,
          "adapter preserves value-only near-pivot classification");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);
}

static void test_projection_range_failure(void) {
  double slater[16];
  double pfaffians[1];
  double complex weights[1] = {1.0};
  double complex parameter = 1000.0;
  int gutz_indices[2] = {0, 0};
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  MVMCKrylovAmplitudeResult result;
  const uint64_t configuration = UINT64_C(5);
  int rank, size;

  rank_and_size(&rank, &size);
  fill_real_slater(slater, 1, -1, pfaffians);
  layout = base_layout(1, rank, size);
  enable_gutzwiller(&layout, gutz_indices, &parameter);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "overflow projection workspace create");
  if (workspace != NULL) {
    CHECK(mvmc_classic_krylov_real_amplitude(
              &configuration, 1, workspace, &result) ==
                  MVMC_KRYLOV_STATUS_NONFINITE &&
              result.value == 0.0,
          "projection exponential overflow is not a physical zero");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);

  workspace = NULL;
  parameter = -1000.0;
  layout = base_layout(1, rank, size);
  enable_gutzwiller(&layout, gutz_indices, &parameter);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "underflow projection workspace create");
  if (workspace != NULL) {
    CHECK(mvmc_classic_krylov_real_amplitude(
              &configuration, 1, workspace, &result) ==
                  MVMC_KRYLOV_STATUS_NONFINITE &&
              result.value == 0.0,
          "projection exponential underflow is not a physical zero");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);
}

static void test_split_communicator(void) {
#ifdef _mpi_use
  double slater[16];
  double pfaffians[1];
  double complex weights[1] = {1.0};
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  MVMCKrylovAmplitudeResult result;
  uint64_t configuration = UINT64_C(5);
  MPI_Comm split = MPI_COMM_NULL;
  int world_rank, world_size, rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size < 2) return;
  MPI_Comm_split(MPI_COMM_WORLD, world_rank % 2, world_rank, &split);
  MPI_Comm_rank(split, &rank);
  MPI_Comm_size(split, &size);
  layout = base_layout(1, rank, size);
  layout.communicator = split;
  fill_real_slater(slater, 1, -1, pfaffians);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "split-communicator amplitude workspace create");
  if (workspace != NULL) {
    CHECK(mvmc_classic_krylov_real_amplitude(
              &configuration, 1, workspace, &result) ==
                  MVMC_KRYLOV_STATUS_OK &&
              close_complex(result.value, pfaffians[0], 1.0e-14),
          "split-communicator amplitude evaluation");
  }
  mvmc_classic_krylov_real_amplitude_workspace_destroy(workspace);
  MPI_Comm_free(&split);
#endif
}

static void test_rank_mismatched_binding(void) {
#ifdef _mpi_use
  double slater[16];
  double pfaffians[1];
  double complex weights[1] = {1.0};
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  int rank, size;

  rank_and_size(&rank, &size);
  if (size < 2) return;
  fill_real_slater(slater, 1, -1, pfaffians);
  if (rank == size - 1) slater[2 * 4] += 0.125;
  layout = base_layout(1, rank, size);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
                MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE &&
            workspace == NULL,
        "rank-mismatched Slater binding fails bitwise collective audit");

  fill_real_slater(slater, 1, -1, pfaffians);
  layout = base_layout(1, rank, size);
  layout.qp_start = 0;
  layout.qp_end = 1;
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &workspace) ==
                MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE &&
            workspace == NULL,
        "overlapping QP ownership fails create-time partition audit");
#endif
}

static void test_generic_engine_integration(void) {
  double slater[16];
  double pfaffians[1];
  double complex weights[1] = {1.0};
  MVMCClassicKrylovRealAmplitudeWorkspace *amplitude_workspace = NULL;
  MVMCClassicKrylovAmplitudeLayout layout;
  MVMCKrylovLimits limits;
  MVMCKrylovFockModel model;
  MVMCKrylovWorkspace *krylov_workspace = NULL;
  MVMCKrylovResult result;
  MVMCKrylovStatus status;
  const uint64_t configuration = UINT64_C(5);
  int rank, size;

  rank_and_size(&rank, &size);
  fill_real_slater(slater, 1, -1, pfaffians);
  layout = base_layout(1, rank, size);
  CHECK(mvmc_classic_krylov_real_amplitude_workspace_create(
            &layout, slater, weights, &amplitude_workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "generic integration amplitude workspace create");
  if (amplitude_workspace == NULL) return;
  memset(&limits, 0, sizeof(limits));
  limits.max_states = 8;
  limits.max_transitions = 8;
  limits.max_amplitude_evaluations = 8;
  limits.max_bytes = 65536;
  limits.max_order = 0;
  krylov_workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(krylov_workspace != NULL && status == MVMC_KRYLOV_STATUS_OK,
        "generic integration Krylov workspace create");
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.down_electron_count = 1;
  model.hermitian = 1;
  if (krylov_workspace != NULL) {
    CHECK(mvmc_krylov_evaluate(
              krylov_workspace, &model, &configuration, 1,
              mvmc_classic_krylov_real_amplitude, amplitude_workspace,
              &result) == MVMC_KRYLOV_STATUS_OK &&
              result.valid && result.evaluated_order == 0 &&
              close_complex(result.value[0], pfaffians[0], 1.0e-14) &&
              result.statistics.global_factorization_count == 1,
          "classic amplitude callback integrates with generic Krylov engine");
  }
  mvmc_krylov_workspace_destroy(krylov_workspace);
  mvmc_classic_krylov_real_amplitude_workspace_destroy(amplitude_workspace);
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  test_matrix_helper();
  test_real_amplitude(1, 0, -1, 1);
  test_real_amplitude(4, 1, 1, 1);
  test_complex_amplitude(1, 0, -1);
  test_complex_amplitude(4, 1, 2);
  test_all_zero_and_configuration_guards();
  test_projection_allowlist();
  test_weight_binding_guards();
  test_rank_invariant_cancellation();
  test_nonfinite_and_near_pivot_classification();
  test_projection_range_failure();
  test_rank_mismatched_binding();
  test_split_communicator();
  test_generic_engine_integration();

  if (failure_count != 0) {
    fprintf(stderr, "%d classic Krylov amplitude checks failed\n",
            failure_count);
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 1;
  }
  puts("classic Krylov amplitude unit checks passed");
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return 0;
}
