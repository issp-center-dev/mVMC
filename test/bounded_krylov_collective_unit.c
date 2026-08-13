#include "bounded_krylov_collective.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;
static int world_rank = 0;
static int world_size = 1;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "BoundedKrylovCollective_Unit FAIL rank %d: ",          \
              world_rank);                                                      \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, " (line %d)\n", __LINE__);                              \
      ++failures;                                                               \
    }                                                                           \
  } while (0)

static MVMCScaledComplex finite_value(double complex raw) {
  MVMCScaledComplex value;
  memset(&value, 0, sizeof(value));
  CHECK(creal(raw) != 0.0 || cimag(raw) != 0.0,
        "finite fixture cannot import raw zero");
  CHECK(mvmc_scaled_complex_from_raw_testing(raw, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "finite scaled fixture construction failed");
  return value;
}

static MVMCAbsolutePfaffianScaledValueResult component_value(
    double complex raw) {
  MVMCAbsolutePfaffianScaledValueResult component;
  memset(&component, 0, sizeof(component));
  component.matrix_scale = 1.0;
  component.scaled_min_pivot = 0.5;
  if (creal(raw) == 0.0 && cimag(raw) == 0.0) {
    component.factor_state = MVMC_PFAFFIAN_VALUE_SINGULAR;
    component.factor_info = 1;
    CHECK(mvmc_scaled_complex_make_exact_zero(&component.value) ==
              MVMC_PFAFFIAN_STATUS_OK,
          "exact-zero component construction failed");
  } else {
    component.factor_state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
    component.value = finite_value(raw);
  }
  return component;
}

static int scaled_fields_bitwise_equal(const MVMCScaledComplex *left,
                                       const MVMCScaledComplex *right) {
  return left->state == right->state &&
         memcmp(&left->phase, &right->phase, sizeof(left->phase)) == 0 &&
         memcmp(&left->log_abs, &right->log_abs,
                sizeof(left->log_abs)) == 0 &&
         memcmp(&left->log_abs_error_bound, &right->log_abs_error_bound,
                sizeof(left->log_abs_error_bound)) == 0 &&
         memcmp(&left->max_input_log_abs, &right->max_input_log_abs,
                sizeof(left->max_input_log_abs)) == 0 &&
         memcmp(&left->cancellation_log_abs,
                &right->cancellation_log_abs,
                sizeof(left->cancellation_log_abs)) == 0 &&
         memcmp(&left->cancellation_ratio, &right->cancellation_ratio,
                sizeof(left->cancellation_ratio)) == 0;
}

static int component_fields_bitwise_equal(
    const MVMCAbsolutePfaffianScaledValueResult *left,
    const MVMCAbsolutePfaffianScaledValueResult *right) {
  return left->factor_state == right->factor_state &&
         left->factor_info == right->factor_info &&
         memcmp(&left->matrix_scale, &right->matrix_scale,
                sizeof(left->matrix_scale)) == 0 &&
         memcmp(&left->scaled_min_pivot, &right->scaled_min_pivot,
                sizeof(left->scaled_min_pivot)) == 0 &&
         scaled_fields_bitwise_equal(&left->value, &right->value);
}

static void rank_partition(int total, int rank, int size,
                           int *start, int *end) {
  *start = total * rank / size;
  *end = total * (rank + 1) / size;
}

static void test_workspace_creation(void) {
  MVMCKrylovBoundedCollectiveWorkspace *workspace = NULL;
  MVMCKrylovStatus status;
#ifdef _mpi_use
  const size_t mismatched = world_rank == world_size - 1 ? 5U : 4U;
  status = mvmc_bounded_krylov_collective_workspace_create(
      MPI_COMM_WORLD, mismatched, 3, &workspace);
  if (world_size > 1) {
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
          "rank-mismatched creation capacity was accepted");
  } else {
    CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL,
          "single-rank creation fixture failed");
    mvmc_bounded_krylov_collective_workspace_destroy(workspace);
    workspace = NULL;
  }
  status = mvmc_bounded_krylov_collective_workspace_create(
      MPI_COMM_WORLD, SIZE_MAX, 3, &workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "overflowing QP capacity was accepted");
  status = mvmc_bounded_krylov_collective_workspace_create(
      MPI_COMM_WORLD, 4, SIZE_MAX, &workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "overflowing accumulator capacity was accepted");
#else
  status = mvmc_bounded_krylov_collective_workspace_create(0, 0, 3,
                                                            &workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "zero QP capacity was accepted");
  status = mvmc_bounded_krylov_collective_workspace_create(
      0, SIZE_MAX, 3, &workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "overflowing QP capacity was accepted");
  status = mvmc_bounded_krylov_collective_workspace_create(
      0, 4, SIZE_MAX, &workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "overflowing accumulator capacity was accepted");
#endif
}

static void test_synchronize(
    MVMCKrylovBoundedCollectiveWorkspace *workspace) {
  const int target = world_size - 1;
  const uint64_t hash = UINT64_C(0x123456789abcdef0);
  MVMCKrylovBoundedCollectiveResult result;
  MVMCKrylovStatus status;
  status = mvmc_bounded_krylov_collective_synchronize(
      workspace, MVMC_KRYLOV_STATUS_OK, hash, 1, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            result.failure_rank == -1 && result.plan_hash == hash,
        "equal plan synchronization failed");

  status = mvmc_bounded_krylov_collective_synchronize(
      workspace, MVMC_KRYLOV_STATUS_OK,
      world_rank == target && world_size > 1 ? hash + 1 : hash, 1, &result);
  if (world_size > 1) {
    CHECK(status == MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
              !result.valid && result.failure_rank == target,
          "bad plan hash was not propagated from target rank");
  } else {
    CHECK(status == MVMC_KRYLOV_STATUS_OK,
          "single-rank plan hash synchronization failed");
  }

  status = mvmc_bounded_krylov_collective_synchronize(
      workspace, MVMC_KRYLOV_STATUS_OK, hash,
      world_rank == target ? 0 : 1, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE &&
            result.failure_rank == target,
        "one-rank allocation failure did not propagate");

  status = mvmc_bounded_krylov_collective_synchronize(
      workspace, MVMC_KRYLOV_STATUS_OK, hash,
      world_rank == target ? 2 : 1, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target,
        "one-rank invalid allocation flag did not remain collective");

  status = mvmc_bounded_krylov_collective_synchronize(
      workspace,
      world_rank == target ? MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE
                           : MVMC_KRYLOV_STATUS_OK,
      hash, 1, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            result.failure_rank == target,
        "one-rank callback failure did not propagate");

  status = mvmc_bounded_krylov_collective_synchronize(
      workspace, MVMC_KRYLOV_STATUS_OK, hash, 1,
      world_rank == target ? NULL : &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "rank-local null result did not remain collective");
  if (world_rank != target) {
    CHECK(result.failure_rank == target,
          "null-result failure rank mismatch");
  }
}

static void test_max_u64(
    MVMCKrylovBoundedCollectiveWorkspace *workspace) {
  uint64_t maximum = 0;
  const int target = world_size - 1;
  CHECK(mvmc_bounded_krylov_collective_max_u64(
            workspace, (uint64_t)(world_rank + 3), &maximum) ==
            MVMC_KRYLOV_STATUS_OK &&
            maximum == (uint64_t)(world_size + 2),
        "max_u64 mismatch");
  CHECK(mvmc_bounded_krylov_collective_max_u64(
            workspace, (uint64_t)world_rank,
            world_rank == target ? NULL : &maximum) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "rank-local null max output did not propagate");
}

static void test_component_gather(
    MVMCKrylovBoundedCollectiveWorkspace *workspace) {
  const MVMCAbsolutePfaffianScaledValueResult expected[4] = {
      component_value(1.0 + 2.0 * I),
      component_value(0.0),
      component_value(-0.25 + 0.5 * I),
      component_value(3.0 - I),
  };
  MVMCAbsolutePfaffianScaledValueResult local[4];
  MVMCAbsolutePfaffianScaledValueResult global[4];
  MVMCKrylovBoundedCollectiveResult result;
  MVMCKrylovStatus status;
  int start;
  int end;
  int index;
  int total = 2;
  rank_partition(total, world_rank, world_size, &start, &end);
  memset(local, 0, sizeof(local));
  memset(global, 0x5a, sizeof(global));
  for (index = start; index < end; ++index) {
    local[index - start] = expected[index];
  }
  status = mvmc_bounded_krylov_collective_gather_scaled_components(
      workspace, MVMC_KRYLOV_STATUS_OK, local, (size_t)(end - start),
      total, start, end, global, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            component_fields_bitwise_equal(&global[0], &expected[0]) &&
            component_fields_bitwise_equal(&global[1], &expected[1]),
        "global QP-order gather failed, including empty slices");

  total = 4;
  rank_partition(total, world_rank, world_size, &start, &end);
  for (index = start; index < end; ++index) {
    local[index - start] = expected[index];
  }
  status = mvmc_bounded_krylov_collective_gather_scaled_components(
      workspace, MVMC_KRYLOV_STATUS_OK, local, (size_t)(end - start),
      total, start, end, global, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK &&
            component_fields_bitwise_equal(&global[0], &expected[0]) &&
            component_fields_bitwise_equal(&global[1], &expected[1]) &&
            component_fields_bitwise_equal(&global[2], &expected[2]) &&
            component_fields_bitwise_equal(&global[3], &expected[3]),
        "four-component global gather order mismatch");

  {
    const int target = world_size - 1;
    MVMCAbsolutePfaffianScaledValueResult sentinel[4];
    memcpy(sentinel, global, sizeof(sentinel));
    rank_partition(total, world_rank, world_size, &start, &end);
    status = mvmc_bounded_krylov_collective_gather_scaled_components(
        workspace,
        world_rank == target ? MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE
                             : MVMC_KRYLOV_STATUS_OK,
        local, (size_t)(end - start), total, start, end, global, &result);
    CHECK(status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
              result.failure_rank == target &&
              memcmp(global, sentinel, sizeof(sentinel)) == 0,
          "failed gather changed global output or lost failure rank");
  }

  rank_partition(total, world_rank, world_size, &start, &end);
  if (end > start && world_rank == world_size - 1) {
    local[0].value.state = MVMC_SCALED_COMPLEX_NONFINITE;
  }
  status = mvmc_bounded_krylov_collective_gather_scaled_components(
      workspace, MVMC_KRYLOV_STATUS_OK, local, (size_t)(end - start),
      total, start, end, global, &result);
  if (world_size == 1 || end > start) {
    CHECK(status == MVMC_KRYLOV_STATUS_NONFINITE,
          "nonfinite rank-local component was accepted");
  }

  rank_partition(total, world_rank, world_size, &start, &end);
  for (index = start; index < end; ++index) {
    local[index - start] = expected[index];
  }
  if (end > start && world_rank == world_size - 1) {
    local[0].factor_state = MVMC_PFAFFIAN_VALUE_NONFINITE;
  }
  status = mvmc_bounded_krylov_collective_gather_scaled_components(
      workspace, MVMC_KRYLOV_STATUS_OK, local, (size_t)(end - start),
      total, start, end, global, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_NONFINITE,
        "nonfinite factor state was not propagated");

  for (index = start; index < end; ++index) {
    local[index - start] = expected[index];
  }
  if (end > start && world_rank == world_size - 1) {
    local[0].factor_state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
    CHECK(mvmc_scaled_complex_make_exact_zero(&local[0].value) ==
              MVMC_PFAFFIAN_STATUS_OK,
          "inconsistent factor/value fixture construction failed");
  }
  status = mvmc_bounded_krylov_collective_gather_scaled_components(
      workspace, MVMC_KRYLOV_STATUS_OK, local, (size_t)(end - start),
      total, start, end, global, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "inconsistent factor/value state was accepted");
}

static void test_projected_amplitude(
    MVMCKrylovBoundedCollectiveWorkspace *workspace) {
  const MVMCAbsolutePfaffianScaledValueResult components[4] = {
      component_value(1.0 + 2.0 * I),
      component_value(0.0),
      component_value(-0.25 + 0.5 * I),
      component_value(3.0 - I),
  };
  double complex weights[4] = {1.0, -0.5 + 0.25 * I, 2.0, -I};
  MVMCAbsolutePfaffianScaledValueResult local[4];
  MVMCProjectedScaledAmplitudeResult oracle;
  MVMCKrylovScaledAmplitudeResult amplitude;
  MVMCKrylovScaledAmplitudeResult sentinel;
  MVMCKrylovBoundedCollectiveResult result;
  MVMCKrylovStatus status;
  int start;
  int end;
  int index;
  rank_partition(4, world_rank, world_size, &start, &end);
  for (index = start; index < end; ++index) {
    local[index - start] = components[index];
  }
  CHECK(mvmc_projected_scaled_amplitude_values(components, weights, 4,
                                                &oracle) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "allocating P4-A projection oracle failed");
  memset(&amplitude, 0x5a, sizeof(amplitude));
  status =
      mvmc_bounded_krylov_collective_gather_projected_amplitude(
          workspace, MVMC_KRYLOV_STATUS_OK, local, (size_t)(end - start),
          weights, 4, start, end, &amplitude, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            scaled_fields_bitwise_equal(&amplitude.value, &oracle.total) &&
            amplitude.well_pivoted_component_count == 3 &&
            amplitude.exact_zero_component_count == 1 &&
            amplitude.local_factorization_count ==
                (uint64_t)(unsigned int)(end - start) &&
            amplitude.global_factorization_count == 4,
        "allocation-free projected amplitude disagrees with P4-A oracle");

  rank_partition(1, world_rank, world_size, &start, &end);
  if (end > start) local[0] = components[0];
  status =
      mvmc_bounded_krylov_collective_gather_projected_amplitude(
          workspace, MVMC_KRYLOV_STATUS_OK, local, (size_t)(end - start),
          weights, 1, start, end, &amplitude, &result);
  {
    double complex qp_one_raw = NAN + I * NAN;
    CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
              mvmc_scaled_complex_export_common_scale(
                  &amplitude.value, 0.0, &qp_one_raw) ==
                  MVMC_SCALED_EXPORT_OK &&
              cabs(qp_one_raw - (1.0 + 2.0 * I)) <= 1.0e-14 &&
              amplitude.global_factorization_count == 1,
          "QP=1 projection with empty rank slices failed");
  }

  sentinel = amplitude;
  rank_partition(4, world_rank, world_size, &start, &end);
  for (index = start; index < end; ++index) {
    local[index - start] = components[index];
  }
  if (world_rank == world_size - 1 && world_size > 1) weights[0] += 1.0;
  status =
      mvmc_bounded_krylov_collective_gather_projected_amplitude(
          workspace, MVMC_KRYLOV_STATUS_OK, local, (size_t)(end - start),
          weights, 4, start, end, &amplitude, &result);
  if (world_size > 1) {
    CHECK(status == MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
              memcmp(&amplitude, &sentinel, sizeof(amplitude)) == 0,
          "replicated QP-weight mismatch was accepted or changed output");
  } else {
    CHECK(status == MVMC_KRYLOV_STATUS_OK,
          "single-rank changed-weight projection failed");
  }
}

static void test_ordered_accumulator_merge(
    MVMCKrylovBoundedCollectiveWorkspace *workspace) {
  MVMCScaledComplex local[3];
  MVMCScaledComplex global[3];
  MVMCScaledComplex expected[3];
  MVMCScaledComplex *rank_values =
      (MVMCScaledComplex *)calloc((size_t)world_size, sizeof(*rank_values));
  MVMCKrylovBoundedCollectiveResult result;
  MVMCKrylovStatus status;
  int rank;
  int index;
  CHECK(rank_values != NULL, "test rank scratch allocation failed");
  if (rank_values == NULL) return;
  local[0] = finite_value((double)(world_rank + 1));
  if (world_rank == 0) {
    local[1] = finite_value(1.0e16);
  } else if (world_rank == 1) {
    local[1] = finite_value(1.0);
  } else if (world_rank == 2) {
    local[1] = finite_value(-1.0e16);
  } else {
    local[1] = finite_value(2.0);
  }
  if ((world_rank & 1) == 0) {
    CHECK(mvmc_scaled_complex_make_exact_zero(&local[2]) ==
              MVMC_PFAFFIAN_STATUS_OK,
          "merge exact-zero fixture failed");
  } else {
    local[2] = finite_value((double)world_rank * I);
  }
  for (index = 0; index < 3; ++index) {
    for (rank = 0; rank < world_size; ++rank) {
      if (index == 0) {
        rank_values[rank] = finite_value((double)(rank + 1));
      } else if (index == 1) {
        rank_values[rank] = finite_value(
            rank == 0 ? 1.0e16
                      : (rank == 1 ? 1.0 : (rank == 2 ? -1.0e16 : 2.0)));
      } else if ((rank & 1) == 0) {
        CHECK(mvmc_scaled_complex_make_exact_zero(&rank_values[rank]) ==
                  MVMC_PFAFFIAN_STATUS_OK,
              "expected exact-zero fixture failed");
      } else {
        rank_values[rank] = finite_value((double)rank * I);
      }
    }
    CHECK(mvmc_scaled_complex_sum_ordered(rank_values, (size_t)world_size,
                                           &expected[index]) ==
              MVMC_PFAFFIAN_STATUS_OK,
          "independent ordered merge construction failed");
  }
  status = mvmc_bounded_krylov_collective_merge_scaled_accumulators(
      workspace, MVMC_KRYLOV_STATUS_OK, local, 3, global, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            scaled_fields_bitwise_equal(&global[0], &expected[0]) &&
            scaled_fields_bitwise_equal(&global[1], &expected[1]) &&
            scaled_fields_bitwise_equal(&global[2], &expected[2]),
        "rank-ordered scaled accumulator merge mismatch");

  {
    const int target = world_size - 1;
    MVMCScaledComplex sentinel[3];
    memcpy(sentinel, global, sizeof(sentinel));
    status = mvmc_bounded_krylov_collective_merge_scaled_accumulators(
        workspace,
        world_rank == target ? MVMC_KRYLOV_STATUS_RESOURCE_LIMIT
                             : MVMC_KRYLOV_STATUS_OK,
        local, 3, global, &result);
    CHECK(status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
              result.failure_rank == target &&
              memcmp(global, sentinel, sizeof(sentinel)) == 0,
          "failed ordered merge changed output or failure rank");
  }

  status = mvmc_bounded_krylov_collective_merge_scaled_accumulators(
      workspace, MVMC_KRYLOV_STATUS_OK, local,
      world_rank == world_size - 1 && world_size > 1 ? 2U : 3U,
      global, &result);
  if (world_size > 1) {
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "rank-mismatched accumulator count was accepted");
  } else {
    CHECK(status == MVMC_KRYLOV_STATUS_OK,
          "single-rank shortened accumulator merge failed");
  }
  free(rank_values);
}

int main(int argc, char **argv) {
  MVMCKrylovBoundedCollectiveWorkspace *workspace = NULL;
  MVMCKrylovStatus status;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
  (void)argc;
  (void)argv;
#endif
  test_workspace_creation();
#ifdef _mpi_use
  status = mvmc_bounded_krylov_collective_workspace_create(
      MPI_COMM_WORLD, 4, 3, &workspace);
#else
  status = mvmc_bounded_krylov_collective_workspace_create(0, 4, 3,
                                                            &workspace);
#endif
  CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL &&
            mvmc_bounded_krylov_collective_workspace_bytes(workspace) > 0,
        "collective workspace creation failed");
  if (workspace != NULL) {
    test_synchronize(workspace);
    test_max_u64(workspace);
    test_component_gather(workspace);
    test_projected_amplitude(workspace);
    test_ordered_accumulator_merge(workspace);
  }
  mvmc_bounded_krylov_collective_workspace_destroy(workspace);

#ifdef _mpi_use
  {
    int global_failures = 0;
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failures = global_failures;
  }
#endif
  if (failures != 0) {
    if (world_rank == 0) {
      fprintf(stderr, "BoundedKrylovCollective_Unit: %d failure(s)\n",
              failures);
    }
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return EXIT_FAILURE;
  }
  if (world_rank == 0) {
    printf("BoundedKrylovCollective_Unit: PASS (%d rank%s)\n", world_size,
           world_size == 1 ? "" : "s");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
