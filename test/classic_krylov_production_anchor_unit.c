#include "bounded_krylov_engine.h"
#include "classic_krylov_amplitude.h"
#include "classic_krylov_model.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;
static int world_rank = 0;
static int world_size = 1;

#define CHECK(condition, ...)                                               \
  do {                                                                      \
    if (!(condition)) {                                                     \
      fprintf(stderr, "ClassicKrylovProductionAnchor FAIL rank %d: ",     \
              world_rank);                                                  \
      fprintf(stderr, __VA_ARGS__);                                         \
      fprintf(stderr, " (line %d)\n", __LINE__);                          \
      ++failure_count;                                                       \
    }                                                                       \
  } while (0)

static MVMCClassicPfaffianCommunicator world_communicator(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static void set_real_pair(double *slater, int up, int down,
                          double pfaffian) {
  const int site_count = 2;
  const int orbital_count = 2 * site_count;
  const int down_orbital = site_count + down;
  slater[(size_t)down_orbital * (size_t)orbital_count + (size_t)up] =
      -pfaffian;
  slater[(size_t)up * (size_t)orbital_count +
         (size_t)down_orbital] = pfaffian;
}

static void fill_slater(double *slater, int qp_total, int singular_qp) {
  int qp;
  int up;
  int down;
  memset(slater, 0, (size_t)qp_total * 16u * sizeof(*slater));
  for (qp = 0; qp < qp_total; ++qp) {
    double *component = slater + (size_t)qp * 16u;
    const double multiplier = qp == singular_qp ? 0.0 : (double)(qp + 1);
    for (up = 0; up < 2; ++up) {
      for (down = 0; down < 2; ++down) {
        set_real_pair(component, up, down,
                      multiplier * (1.5 + (double)up +
                                    0.25 * (double)down));
      }
    }
  }
}

static MVMCKrylovBoundedLimits production_limits(uint64_t generation_hash) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash = generation_hash;
  limits.cache_bytes = 4096;
  limits.max_row_transitions = 16;
  limits.max_workspace_bytes = (size_t)1024 * 1024;
  limits.max_node_expansions = 1000;
  limits.max_terminal_amplitude_calls = 1000;
  limits.max_total_row_transitions = 10000;
  limits.max_order = 3;
  return limits;
}

static int scaled_is_finite_nonzero(const MVMCScaledComplex *value) {
  double complex raw = NAN + I * NAN;
  return value->state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
         mvmc_scaled_complex_export_common_scale(value, 0.0, &raw) ==
             MVMC_SCALED_EXPORT_OK &&
         isfinite(creal(raw)) && isfinite(cimag(raw)) && cabs(raw) > 0.0;
}

static void check_rank_invariance(const MVMCKrylovBoundedResult *result,
                                  uint64_t generation_hash,
                                  uint64_t plan_hash) {
#ifdef _mpi_use
  uint64_t local[8];
  uint64_t minimum[8];
  uint64_t maximum[8];
  int order;
  local[0] = generation_hash;
  local[1] = plan_hash;
  local[2] = result->statistics.trace_hash;
  local[3] = result->statistics.node_expansions;
  local[4] = result->statistics.terminal_amplitude_calls;
  local[5] = result->statistics.global_factorization_count;
  local[6] = result->statistics.exact_zero_component_count;
  local[7] = (uint64_t)(unsigned int)result->status;
  MPI_Allreduce(local, minimum, 8, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(local, maximum, 8, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  CHECK(memcmp(minimum, maximum, sizeof(minimum)) == 0,
        "production metadata/statistics differ by rank");
  for (order = 0; order <= 3; ++order) {
    const double local_value[6] = {
        creal(result->value[order].phase),
        cimag(result->value[order].phase),
        result->value[order].log_abs,
        result->value[order].log_abs_error_bound,
        result->value[order].max_input_log_abs,
        result->value[order].cancellation_ratio,
    };
    double minimum_value[6];
    double maximum_value[6];
    const int local_state = (int)result->value[order].state;
    int minimum_state;
    int maximum_state;
    MPI_Allreduce(local_value, minimum_value, 6, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local_value, maximum_value, 6, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_state, &minimum_state, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_state, &maximum_state, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    CHECK(memcmp(minimum_value, maximum_value, sizeof(minimum_value)) == 0 &&
              minimum_state == maximum_state,
          "production v%d differs by rank", order);
  }
#else
  (void)result;
  (void)generation_hash;
  (void)plan_hash;
#endif
}

static void run_production_anchor(void) {
  enum { qp_total = 8, singular_qp = 3, nproj = 20 };
  MVMCClassicKrylovTransfer transfers[4] = {
      {1, 0, 0, 0, 0.75},
      {0, 0, 1, 0, 0.75},
      {1, 1, 0, 1, 0.50},
      {0, 1, 1, 1, 0.50},
  };
  MVMCClassicKrylovSiteCoupling intra[1] = {{0, 0.40}};
  MVMCClassicKrylovRawModel raw_model;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  const MVMCKrylovFockModel *model = NULL;
  double slater[qp_total * 16];
  double slater_snapshot[qp_total * 16];
  double complex weights[qp_total];
  double complex weights_snapshot[qp_total];
  double complex parameters[nproj];
  double complex parameters_snapshot[nproj];
  int gutzwiller[2] = {0, 1};
  int gutzwiller_snapshot[2];
  int jastrow_row_0[2] = {-1, 0};
  int jastrow_row_1[2] = {0, -1};
  const int *jastrow[2] = {jastrow_row_0, jastrow_row_1};
  int spin_row_0[2] = {-1, 0};
  int spin_row_1[2] = {0, -1};
  const int *spin_jastrow[2] = {spin_row_0, spin_row_1};
  int dh2_row[4] = {1, 1, 0, 0};
  const int *dh2[1] = {dh2_row};
  int dh4_row[8] = {1, 1, 1, 1, 0, 0, 0, 0};
  const int *dh4[1] = {dh4_row};
  int binding_snapshot[20];
  MVMCClassicKrylovAmplitudeLayout amplitude_layout;
  MVMCClassicKrylovRealAmplitudeWorkspace *amplitude_workspace = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *bounded_workspace = NULL;
  MVMCKrylovBoundedResult result;
  uint64_t root = UINT64_C(5);
  const uint64_t root_snapshot = root;
  uint64_t generation_hash = 0;
  uint64_t plan_hash = 0;
  MVMCKrylovStatus status;
  int index;

  memset(&raw_model, 0, sizeof(raw_model));
  raw_model.site_count = 2;
  raw_model.up_electron_count = 1;
  raw_model.down_electron_count = 1;
  raw_model.transfer_count = 4;
  raw_model.transfers = transfers;
  raw_model.coulomb_intra_count = 1;
  raw_model.coulomb_intra = intra;

  fill_slater(slater, qp_total, singular_qp);
  for (index = 0; index < qp_total; ++index) {
    weights[index] = 1.0 / (double)(index + 1);
  }
  for (index = 0; index < nproj; ++index) {
    parameters[index] = (double)(index + 1) / 100.0;
  }
  memcpy(slater_snapshot, slater, sizeof(slater));
  memcpy(weights_snapshot, weights, sizeof(weights));
  memcpy(parameters_snapshot, parameters, sizeof(parameters));
  memcpy(gutzwiller_snapshot, gutzwiller, sizeof(gutzwiller));
  memcpy(binding_snapshot, jastrow_row_0, sizeof(jastrow_row_0));
  memcpy(binding_snapshot + 2, jastrow_row_1, sizeof(jastrow_row_1));
  memcpy(binding_snapshot + 4, spin_row_0, sizeof(spin_row_0));
  memcpy(binding_snapshot + 6, spin_row_1, sizeof(spin_row_1));
  memcpy(binding_snapshot + 8, dh2_row, sizeof(dh2_row));
  memcpy(binding_snapshot + 12, dh4_row, sizeof(dh4_row));

  status = mvmc_classic_krylov_model_workspace_create_from_root(
      world_rank == 0 ? &raw_model : NULL, world_communicator(),
      &model_workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && model_workspace != NULL,
        "production model adapter failed with status %d", (int)status);
  if (status != MVMC_KRYLOV_STATUS_OK || model_workspace == NULL) goto cleanup;
  model = mvmc_classic_krylov_model(model_workspace);
  CHECK(model != NULL && model->term_count == 5 &&
            model->operator_count == 12,
        "production model expansion is unexpected");
  if (model == NULL) goto cleanup;

  memset(&amplitude_layout, 0, sizeof(amplitude_layout));
  amplitude_layout.site_count = 2;
  amplitude_layout.up_electron_count = 1;
  amplitude_layout.down_electron_count = 1;
  amplitude_layout.qp_total = qp_total;
  amplitude_layout.qp_start = qp_total * world_rank / world_size;
  amplitude_layout.qp_end = qp_total * (world_rank + 1) / world_size;
  amplitude_layout.nproj = nproj;
  amplitude_layout.ngutzwiller_idx = 2;
  amplitude_layout.njastrow_idx = 1;
  amplitude_layout.nspin_jastrow_idx = 1;
  amplitude_layout.ndoublon_holon_2site_idx = 1;
  amplitude_layout.ndoublon_holon_4site_idx = 1;
  amplitude_layout.gutzwiller_idx = gutzwiller;
  amplitude_layout.jastrow_idx = jastrow;
  amplitude_layout.spin_jastrow_idx = spin_jastrow;
  amplitude_layout.doublon_holon_2site_idx = dh2;
  amplitude_layout.doublon_holon_4site_idx = dh4;
  amplitude_layout.projection_parameters = parameters;
  amplitude_layout.communicator = world_communicator();
  status = mvmc_classic_krylov_real_amplitude_workspace_create(
      &amplitude_layout, slater, weights, &amplitude_workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && amplitude_workspace != NULL,
        "production amplitude adapter failed with status %d", (int)status);
  if (status != MVMC_KRYLOV_STATUS_OK || amplitude_workspace == NULL) {
    goto cleanup;
  }
  generation_hash =
      mvmc_classic_krylov_real_amplitude_generation_hash(amplitude_workspace);
  CHECK(generation_hash != 0,
        "production amplitude generation hash is zero");

  limits = production_limits(generation_hash);
  status = mvmc_bounded_krylov_plan_create(model, &limits, &plan);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && plan != NULL,
        "production bounded plan failed with status %d", (int)status);
  if (status != MVMC_KRYLOV_STATUS_OK || plan == NULL) goto cleanup;
  plan_hash = mvmc_bounded_krylov_plan_hash(plan);
  CHECK(plan_hash != 0 && mvmc_bounded_krylov_plan_max_row_transitions(plan) ==
                              model->term_count,
        "production bounded plan metadata is invalid");
  status = mvmc_bounded_krylov_workspace_create(plan, &bounded_workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && bounded_workspace != NULL,
        "production bounded workspace failed with status %d", (int)status);
  if (status != MVMC_KRYLOV_STATUS_OK || bounded_workspace == NULL) {
    goto cleanup;
  }

  status = mvmc_bounded_krylov_evaluate(
      bounded_workspace, &root, 1,
      mvmc_classic_krylov_real_scaled_amplitude, amplitude_workspace,
      &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            result.evaluated_order == 3 &&
            result.statistics.completed_order == 3,
        "production bounded v0..v3 failed with status %d", (int)status);
  if (status == MVMC_KRYLOV_STATUS_OK && result.valid) {
    for (index = 0; index <= 3; ++index) {
      CHECK(scaled_is_finite_nonzero(&result.value[index]),
            "production v%d is not finite and nonzero", index);
    }
    CHECK(result.statistics.terminal_amplitude_calls > 0 &&
              result.statistics.global_factorization_count ==
                  result.statistics.terminal_amplitude_calls * qp_total,
          "production QP factorization accounting is inconsistent");
    CHECK(result.statistics.exact_zero_component_count ==
              result.statistics.terminal_amplitude_calls &&
              result.statistics.well_pivoted_component_count ==
                  result.statistics.terminal_amplitude_calls *
                      (qp_total - 1),
          "singular QP component was not isolated from projected totals");
    CHECK(result.statistics.computed_exact_zero_values == 0,
          "nonzero projected production amplitudes became exact zero");
    CHECK(result.statistics.engine_heap_allocations == 0,
          "bounded evaluation allocated from the heap");
    check_rank_invariance(&result, generation_hash, plan_hash);
  }

  CHECK(root == root_snapshot &&
            memcmp(slater, slater_snapshot, sizeof(slater)) == 0 &&
            memcmp(weights, weights_snapshot, sizeof(weights)) == 0 &&
            memcmp(parameters, parameters_snapshot, sizeof(parameters)) ==
                0 &&
            memcmp(gutzwiller, gutzwiller_snapshot, sizeof(gutzwiller)) ==
                0 &&
            memcmp(jastrow_row_0, binding_snapshot,
                   sizeof(jastrow_row_0)) == 0 &&
            memcmp(jastrow_row_1, binding_snapshot + 2,
                   sizeof(jastrow_row_1)) == 0 &&
            memcmp(spin_row_0, binding_snapshot + 4,
                   sizeof(spin_row_0)) == 0 &&
            memcmp(spin_row_1, binding_snapshot + 6,
                   sizeof(spin_row_1)) == 0 &&
            memcmp(dh2_row, binding_snapshot + 8, sizeof(dh2_row)) == 0 &&
            memcmp(dh4_row, binding_snapshot + 12, sizeof(dh4_row)) == 0,
        "production model/amplitude evaluation mutated caller state");

cleanup:
  mvmc_bounded_krylov_workspace_destroy(bounded_workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_classic_krylov_real_amplitude_workspace_destroy(amplitude_workspace);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
  (void)argc;
  (void)argv;
#endif
  run_production_anchor();

#ifdef _mpi_use
  {
    int global_failures = 0;
    MPI_Allreduce(&failure_count, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failure_count = global_failures;
  }
#endif
  if (failure_count != 0) {
    if (world_rank == 0) {
      fprintf(stderr, "%d classic production anchor checks failed\n",
              failure_count);
    }
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 1;
  }
  if (world_rank == 0) {
    puts("classic Krylov production anchor checks passed");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return 0;
}
