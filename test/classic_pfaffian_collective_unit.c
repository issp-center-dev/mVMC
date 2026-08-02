#include "classic_pfaffian_collective.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;
static int world_rank = 0;
static int world_size = 1;

#define CHECK(condition, message)                                         \
  do {                                                                    \
    if (!(condition)) {                                                   \
      fprintf(stderr, "FAIL rank %d: %s (line %d)\n", world_rank,       \
              (message), __LINE__);                                       \
      ++failure_count;                                                    \
    }                                                                     \
  } while (0)

static int close_complex(double complex actual, double complex expected,
                         double tolerance) {
  return cabs(actual - expected) <= tolerance * fmax(1.0, cabs(expected));
}

static void rank_partition(int qp_total, int rank, int size,
                           int *qp_start, int *qp_end) {
  *qp_start = (qp_total * rank) / size;
  *qp_end = (qp_total * (rank + 1)) / size;
}

static void make_local_aggregate(int qp_total, int qp_start, int qp_end,
                                 MVMCProjectedAmplitudeResult *aggregate) {
  static const double complex terms[2] = {1.0 + 2.0 * I,
                                           -0.25 + 0.5 * I};
  int qp;

  memset(aggregate, 0, sizeof(*aggregate));
  for (qp = qp_start; qp < qp_end; ++qp) {
    aggregate->total += terms[qp];
    aggregate->sum_abs += cabs(terms[qp]);
    if (qp == 0) {
      ++aggregate->regular_count;
    } else {
      ++aggregate->near_singular_count;
    }
  }
  aggregate->cancellation_ratio =
      aggregate->sum_abs == 0.0
          ? 0.0
          : fmin(1.0, cabs(aggregate->total) / aggregate->sum_abs);
  aggregate->valid = qp_total == 2;
}

static void test_global_aggregate(
    MVMCClassicPfaffianCollectiveWorkspace *workspace) {
  static const MVMCPfaffianStatus lower_severity[3] = {
      MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
      MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE,
      MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE};
  static const MVMCPfaffianStatus higher_severity[3] = {
      MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE,
      MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE,
      MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE};
  const int qp_total = 2;
  const double complex expected_total = 0.75 + 2.5 * I;
  const double expected_sum_abs = cabs(1.0 + 2.0 * I) +
                                  cabs(-0.25 + 0.5 * I);
  MVMCProjectedAmplitudeResult local;
  MVMCClassicPfaffianCollectiveResult result;
  MVMCPfaffianStatus status;
  int qp_start, qp_end;
  int target = world_size - 1;
  int severity_case;

  rank_partition(qp_total, world_rank, world_size, &qp_start, &qp_end);
  make_local_aggregate(qp_total, qp_start, qp_end, &local);
  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      qp_total, qp_start, qp_end, -1, -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.aggregate.valid,
        "valid aggregate succeeds on every rank");
  CHECK(close_complex(result.aggregate.total, expected_total, 2.0e-14),
        "global complex total is reduced");
  CHECK(fabs(result.aggregate.sum_abs - expected_sum_abs) <= 2.0e-14,
        "global sum_abs is reduced");
  CHECK(result.aggregate.regular_count == 1 &&
            result.aggregate.near_singular_count == 1 &&
            result.aggregate.singular_count == 0,
        "global state counters are reduced");

  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, NULL,
      qp_total, qp_start, qp_end, -1, -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT && !result.valid &&
            result.failure_rank == 0,
        "OK status with a null local aggregate fails collectively");

  memset(&local, 0, sizeof(local));
  if (qp_start == 0 && qp_end > qp_start) {
    local.total += 1.0;
    local.sum_abs += 1.0;
    ++local.regular_count;
  }
  if (qp_end == 2 && qp_start < qp_end) {
    local.total -= 1.0;
    local.sum_abs += 1.0;
    ++local.regular_count;
  }
  local.cancellation_ratio =
      local.sum_abs == 0.0
          ? 0.0
          : fmin(1.0, cabs(local.total) / local.sum_abs);
  local.valid = 1;
  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      qp_total, qp_start, qp_end, -1, -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK &&
            result.aggregate.total == 0.0 &&
            result.aggregate.sum_abs == 2.0 &&
            result.aggregate.cancellation_ratio == 0.0,
        "global cancellation ratio is recomputed after reduction");

  make_local_aggregate(qp_total, qp_start, qp_end, &local);

  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace,
      world_rank == target ? MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE
                           : MVMC_PFAFFIAN_STATUS_OK,
      world_rank == target ? NULL : &local,
      qp_total, qp_start, qp_end,
      world_rank == target ? 7 : -1,
      world_rank == target ? 11 : -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE &&
            !result.valid && !result.aggregate.valid,
        "one-rank factorization failure is collective");
  CHECK(result.failure_rank == target && result.failure_local_qp == 7 &&
            result.failure_global_qp == 11,
        "collective failure payload identifies source rank and QP");

  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace,
      world_rank == 0 ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                      : (world_rank == target
                             ? MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE
                             : MVMC_PFAFFIAN_STATUS_OK),
      (world_rank == 0 || world_rank == target) ? NULL : &local,
      qp_total, qp_start, qp_end,
      world_rank == target ? 8 : 3,
      world_rank == target ? 12 : 4, &result);
  CHECK(status == (world_size == 1
                       ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                       : MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE),
        "collective status uses fixed severity order");
  if (world_size > 1) {
    CHECK(result.failure_rank == target && result.failure_local_qp == 8 &&
              result.failure_global_qp == 12,
          "higher-severity failure payload wins deterministically");

    for (severity_case = 0; severity_case < 3; ++severity_case) {
      MVMCPfaffianStatus local_severity =
          world_rank == 0
              ? lower_severity[severity_case]
              : (world_rank == target ? higher_severity[severity_case]
                                      : MVMC_PFAFFIAN_STATUS_OK);
      status = mvmc_classic_pfaffian_collective_aggregate(
          workspace, local_severity,
          local_severity == MVMC_PFAFFIAN_STATUS_OK ? &local : NULL,
          qp_total, qp_start, qp_end,
          20 + world_rank, 30 + world_rank, &result);
      CHECK(status == higher_severity[severity_case] &&
                result.failure_rank == target &&
                result.failure_local_qp == 20 + target &&
                result.failure_global_qp == 30 + target,
            "every adjacent severity level has the documented order");
    }
  }

  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT, NULL,
      qp_total, qp_start, qp_end,
      40 + world_rank, 50 + world_rank, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == 0 && result.failure_local_qp == 40 &&
            result.failure_global_qp == 50,
        "equal severity selects the lowest communicator rank");

  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      qp_total, qp_start,
      world_rank == target ? qp_total + 1 : qp_end,
      -1, -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target,
        "malformed QP ownership is rejected collectively");

  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      qp_total,
      world_rank == target ? (world_size == 1 ? 1 : 0) : qp_start,
      qp_end, -1, -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target,
        "QP ownership gap or overlap is rejected collectively");

  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      world_rank == target ? qp_total + 1 : qp_total,
      qp_start, qp_end, -1, -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target,
        "inconsistent qp_total is rejected collectively");

  if (world_rank == target) local.total = NAN;
  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      qp_total, qp_start, qp_end, -1, -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target,
        "nonfinite rank-local aggregate is rejected collectively");

  make_local_aggregate(qp_total, qp_start, qp_end, &local);
  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      qp_total, qp_start, qp_end, -1, -1,
      world_rank == target ? NULL : &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "rank-local null result remains a collective failure");
  if (world_rank != target) {
    CHECK(result.failure_rank == target,
          "null result failure identifies the source rank");
  }
}

static void test_passed_communicator(void) {
#ifdef _mpi_use
  MPI_Comm subcomm = MPI_COMM_NULL;
  MVMCClassicPfaffianCollectiveWorkspace *workspace = NULL;
  MVMCClassicPfaffianCollectiveResult result;
  MVMCProjectedAmplitudeResult local;
  MVMCPfaffianStatus status;
  int color = world_rank % 2;
  int rank, size, qp_start, qp_end;
  double complex term;

  MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &subcomm);
  MPI_Comm_rank(subcomm, &rank);
  MPI_Comm_size(subcomm, &size);
  status = mvmc_classic_pfaffian_collective_workspace_create(
      subcomm, &workspace);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && workspace != NULL,
        "collective workspace accepts passed subcommunicator");
  rank_partition(1, rank, size, &qp_start, &qp_end);
  memset(&local, 0, sizeof(local));
  term = (double)(color + 1) + (double)(color + 2) * I;
  if (qp_end > qp_start) {
    local.total = term;
    local.sum_abs = cabs(term);
    local.cancellation_ratio = 1.0;
    local.regular_count = 1;
  }
  local.valid = 1;
  status = mvmc_classic_pfaffian_collective_aggregate(
      workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      1, qp_start, qp_end, -1, -1, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK &&
            close_complex(result.aggregate.total, term, 1.0e-15),
        "aggregate is confined to passed subcommunicator");
  mvmc_classic_pfaffian_collective_workspace_destroy(workspace);
  MPI_Comm_free(&subcomm);
#endif
}

static void set_real_pair(double *matrix, int n, int row, int column,
                          double value) {
  matrix[(size_t)row * (size_t)n + (size_t)column] = value;
  matrix[(size_t)column * (size_t)n + (size_t)row] = -value;
}

static void set_complex_pair(double complex *matrix, int n, int row,
                             int column, double complex value) {
  matrix[(size_t)row * (size_t)n + (size_t)column] = value;
  matrix[(size_t)column * (size_t)n + (size_t)row] = -value;
}

static void fill_real_slater(double *matrix, double scale) {
  memset(matrix, 0, 16 * sizeof(*matrix));
  set_real_pair(matrix, 4, 0, 1, 2.0 * scale);
  set_real_pair(matrix, 4, 0, 2, -1.0 * scale);
  set_real_pair(matrix, 4, 0, 3, 0.5 * scale);
  set_real_pair(matrix, 4, 1, 2, 3.0 * scale);
  set_real_pair(matrix, 4, 1, 3, 4.0 * scale);
  set_real_pair(matrix, 4, 2, 3, -2.0 * scale);
}

static void fill_complex_slater(double complex *matrix,
                                double complex scale) {
  memset(matrix, 0, 16 * sizeof(*matrix));
  set_complex_pair(matrix, 4, 0, 1, (1.0 + 0.5 * I) * scale);
  set_complex_pair(matrix, 4, 0, 2, (-0.25 + 0.75 * I) * scale);
  set_complex_pair(matrix, 4, 0, 3, (2.0 - 0.5 * I) * scale);
  set_complex_pair(matrix, 4, 1, 2, (0.5 + 1.25 * I) * scale);
  set_complex_pair(matrix, 4, 1, 3, (-1.0 + 0.25 * I) * scale);
  set_complex_pair(matrix, 4, 2, 3, (1.5 - 0.75 * I) * scale);
}

static void test_real_collective_state(
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace) {
  const int qp_total = 2;
  const int ele_idx[4] = {1, 0, 0, 1};
  const double complex weights[2] = {1.0, 0.5};
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveResult result;
  const MVMCClassicPfaffianState *state;
  double slater[32];
  double failed_slater[32];
  double legacy[34];
  double complex local_weights[2];
  MVMCPfaffianStatus status;
  int qp_start, qp_end;
  int target = world_size - 1;

  rank_partition(qp_total, world_rank, world_size, &qp_start, &qp_end);
  fill_real_slater(slater, 1.0);
  fill_real_slater(slater + 16, 0.75);
  memcpy(failed_slater, slater, sizeof(slater));
  memcpy(local_weights, weights, sizeof(weights));
  memset(legacy, 0, sizeof(legacy));
  status = mvmc_classic_pfaffian_real_workspace_create(
      2, 4, qp_total, qp_start, qp_end, &state_workspace);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && state_workspace != NULL,
        "create rank-local real state workspace");

  status = mvmc_classic_pfaffian_real_prepare_collective(
      state_workspace, collective_workspace, slater, ele_idx, local_weights,
      0.0, 0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && result.valid,
        "real absolute proposal succeeds collectively");
  CHECK(result.aggregate.regular_count == 2 &&
            result.aggregate.near_singular_count == 0 &&
            result.aggregate.singular_count == 0,
        "real collective proposal reports regular components");
  state = mvmc_classic_pfaffian_real_candidate(state_workspace);
  CHECK(state != NULL && state->generation == 1,
        "successful collective proposal leaves local candidate ready");
  if (state_workspace != NULL) {
    status = mvmc_classic_pfaffian_real_publish(
        state_workspace, legacy, sizeof(legacy) / sizeof(legacy[0]));
    CHECK(status == MVMC_PFAFFIAN_STATUS_OK &&
              mvmc_classic_pfaffian_real_accepted(state_workspace)->generation ==
                  1,
          "collective proposal can be published on every rank");
  }

  status = mvmc_classic_pfaffian_real_prepare_collective(
      state_workspace, collective_workspace, slater, ele_idx, local_weights,
      0.0, 0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.aggregate.regular_count == 2,
        "regular-to-regular proposal succeeds collectively");
  if (state_workspace != NULL) {
    status = mvmc_classic_pfaffian_real_publish(
        state_workspace, legacy, sizeof(legacy) / sizeof(legacy[0]));
    CHECK(status == MVMC_PFAFFIAN_STATUS_OK,
          "regular-to-regular state publishes on every rank");
  }

  memset(failed_slater + 16, 0, 16 * sizeof(*failed_slater));
  status = mvmc_classic_pfaffian_real_prepare_collective(
      state_workspace, collective_workspace, failed_slater, ele_idx,
      local_weights, 0.0, 0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.aggregate.regular_count == 1 &&
            result.aggregate.singular_count == 1,
        "regular-to-singular proposal remains a valid global proposal");
  if (state_workspace != NULL) {
    status = mvmc_classic_pfaffian_real_publish(
        state_workspace, legacy, sizeof(legacy) / sizeof(legacy[0]));
    CHECK(status == MVMC_PFAFFIAN_STATUS_OK,
          "regular-to-singular state publishes on every rank");
  }

  status = mvmc_classic_pfaffian_real_prepare_collective(
      state_workspace, collective_workspace, failed_slater, ele_idx,
      local_weights, 0.0, 0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.aggregate.regular_count == 1 &&
            result.aggregate.singular_count == 1,
        "singular-to-singular proposal succeeds collectively");
  if (state_workspace != NULL) {
    status = mvmc_classic_pfaffian_real_publish(
        state_workspace, legacy, sizeof(legacy) / sizeof(legacy[0]));
    CHECK(status == MVMC_PFAFFIAN_STATUS_OK,
          "singular-to-singular state publishes on every rank");
  }

  status = mvmc_classic_pfaffian_real_prepare_collective(
      state_workspace, collective_workspace, slater, ele_idx, local_weights,
      0.0, 0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.aggregate.regular_count == 2 &&
            result.aggregate.singular_count == 0,
        "singular-to-regular proposal succeeds collectively");
  if (state_workspace != NULL) {
    status = mvmc_classic_pfaffian_real_publish(
        state_workspace, legacy, sizeof(legacy) / sizeof(legacy[0]));
    CHECK(status == MVMC_PFAFFIAN_STATUS_OK,
          "singular-to-regular state publishes on every rank");
  }

  if (world_rank == target) local_weights[0] = NAN;
  status = mvmc_classic_pfaffian_real_prepare_collective(
      state_workspace, collective_workspace, slater, ele_idx, local_weights,
      0.0, 0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target &&
            result.failure_global_qp == 0,
        "one-rank real prepare failure propagates with global QP");
  CHECK(state_workspace == NULL ||
            mvmc_classic_pfaffian_real_candidate(state_workspace) == NULL,
        "collective prepare failure discards every local candidate");
  CHECK(state_workspace == NULL ||
            mvmc_classic_pfaffian_real_accepted(state_workspace)->generation ==
                5,
        "collective prepare failure preserves accepted generation");

  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

static void test_complex_collective_state(
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace) {
  const int qp_total = 2;
  const int ele_idx[4] = {1, 0, 0, 1};
  const double complex weights[2] = {1.0 + 0.25 * I, -0.5 + 0.75 * I};
  MVMCClassicPfaffianComplexWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveResult result;
  double complex slater[32];
  MVMCPfaffianStatus status;
  int qp_start, qp_end;

  rank_partition(qp_total, world_rank, world_size, &qp_start, &qp_end);
  fill_complex_slater(slater, 1.0);
  fill_complex_slater(slater + 16, 0.5 - 0.25 * I);
  status = mvmc_classic_pfaffian_complex_workspace_create(
      2, 4, qp_total, qp_start, qp_end, &state_workspace);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && state_workspace != NULL,
        "create rank-local complex state workspace");
  status = mvmc_classic_pfaffian_complex_prepare_collective(
      state_workspace, collective_workspace, slater, ele_idx, weights,
      0.0, 0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.aggregate.regular_count == 2 &&
            isfinite(creal(result.aggregate.total)) &&
            isfinite(cimag(result.aggregate.total)),
        "complex absolute proposal reduces a finite global total");
  mvmc_classic_pfaffian_complex_discard_candidate(state_workspace);
  mvmc_classic_pfaffian_complex_workspace_destroy(state_workspace);
}

int main(int argc, char **argv) {
  MVMCClassicPfaffianCollectiveWorkspace *workspace = NULL;
#ifdef _mpi_use
  MVMCClassicPfaffianCollectiveWorkspace *invalid_workspace = NULL;
#endif
  MVMCPfaffianStatus status;
  int total_failures;

#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  status = mvmc_classic_pfaffian_collective_workspace_create(
      MPI_COMM_NULL, &invalid_workspace);
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            invalid_workspace == NULL,
        "MPI_COMM_NULL workspace creation is rejected");
  status = mvmc_classic_pfaffian_collective_workspace_create(
      MPI_COMM_WORLD, &workspace);
#else
  (void)argc;
  (void)argv;
  status = mvmc_classic_pfaffian_collective_workspace_create(0, &workspace);
#endif
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK && workspace != NULL,
        "create world collective workspace");
  if (workspace != NULL) {
    test_global_aggregate(workspace);
    test_passed_communicator();
    test_real_collective_state(workspace);
    test_complex_collective_state(workspace);
  }
  mvmc_classic_pfaffian_collective_workspace_destroy(workspace);

#ifdef _mpi_use
  MPI_Allreduce(&failure_count, &total_failures, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  if (world_rank == 0 && total_failures == 0) {
    printf("classic_pfaffian_collective_unit: PASS (%d ranks)\n", world_size);
  }
  MPI_Finalize();
#else
  total_failures = failure_count;
  if (total_failures == 0) {
    printf("classic_pfaffian_collective_unit: PASS (1 rank)\n");
  }
#endif
  return total_failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
