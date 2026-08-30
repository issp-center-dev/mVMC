#include "bounded_krylov_engine.h"
#include "classic_krylov_amplitude.h"
#include "classic_krylov_model.h"
#include "power_lanczos_classic_bridge.h"
#include "power_lanczos_corrected_dispatch.h"
#include "power_lanczos_production_session.h"

#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;
static int world_rank = 0;
static int world_size = 1;

static int primitive_support_is_valid(uint8_t support) {
  const uint8_t allowed =
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO |
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO |
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
  return support != 0 && (support & (uint8_t)~allowed) == 0;
}

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

static MVMCKrylovStatus create_anchor_session(
    const MVMCPowerLanczosProductionSessionConfig *config,
    const uint64_t *canonical_words, size_t word_count,
    MVMCPowerLanczosProductionSession **session) {
  MVMCPowerLanczosTerminalSnapshotInput terminal;
  const MVMCPowerLanczosTerminalSnapshotInput *root_terminal = NULL;
  memset(&terminal, 0, sizeof(terminal));
  if (config != NULL && config->mpi_world_rank == 0) {
    terminal.canonical_words = canonical_words;
    terminal.word_count = word_count;
    terminal.configuration_bit_count = config->configuration_bit_count;
    terminal.generation = config->base_generation;
    terminal.terminal_proposal_counter = UINT64_C(4096);
    terminal.expected_particle_count = 2;
    terminal.log_sampling_support = 0.0;
    root_terminal = &terminal;
  }
  return mvmc_power_lanczos_production_session_create(
      config, world_communicator(), root_terminal, session);
}

static int build_observable_plan(
    MVMCPowerLanczosObservablePlan *plan) {
  static const char contents[] =
      "====================\n"
      "NCisAjs 4\n"
      "====================\n"
      "observable rows\n"
      "====================\n"
      "0 0 0 0\n"
      "1 0 1 0\n"
      "0 1 0 1\n"
      "1 1 1 1\n";
  MVMCPowerLanczosObservableCensusStatus status =
      MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
  unsigned char *wire = NULL;
  size_t wire_size = 0;
  uint64_t wire_size_u64 = 0;
  int local_ok = 1;
  mvmc_power_lanczos_observable_plan_init(plan);
  if (world_rank == 0) {
    char path[] = "./power_lanczos_observable_anchor_XXXXXX";
    const char *paths[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT] = {
        path, NULL, NULL};
    char diagnostic[256];
    int descriptor = mkstemp(path);
    FILE *stream = descriptor < 0 ? NULL : fdopen(descriptor, "w");
    if (stream == NULL) {
      if (descriptor >= 0) close(descriptor);
      if (descriptor >= 0) (void)remove(path);
      status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR;
    } else {
      if (fputs(contents, stream) == EOF || fclose(stream) != 0) {
        status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR;
      }
      if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
        status = mvmc_power_lanczos_observable_plan_build_from_files(
            2, 2, paths, plan, diagnostic, sizeof(diagnostic));
      }
      if (remove(path) != 0 &&
          status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
        status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR;
      }
    }
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
      status = mvmc_power_lanczos_observable_plan_wire_size(plan,
                                                            &wire_size);
      wire_size_u64 = (uint64_t)wire_size;
    }
  }
#ifdef _mpi_use
  {
    int encoded = (int)status;
    MPI_Bcast(&encoded, 1, MPI_INT, 0, MPI_COMM_WORLD);
    status = (MVMCPowerLanczosObservableCensusStatus)encoded;
    MPI_Bcast(&wire_size_u64, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    wire_size = (size_t)wire_size_u64;
  }
#endif
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK || wire_size == 0 ||
      wire_size_u64 != (uint64_t)wire_size) {
    mvmc_power_lanczos_observable_plan_destroy(plan);
    return 0;
  }
  wire = (unsigned char *)malloc(wire_size);
  local_ok = wire != NULL;
#ifdef _mpi_use
  {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    local_ok = global_ok;
  }
#endif
  if (!local_ok) {
    free(wire);
    mvmc_power_lanczos_observable_plan_destroy(plan);
    return 0;
  }
  if (world_rank == 0) {
    size_t packed_size = 0;
    status = mvmc_power_lanczos_observable_plan_pack(
        plan, wire, wire_size, &packed_size);
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
        packed_size != wire_size) {
      status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW;
    }
  }
#ifdef _mpi_use
  {
    int encoded = (int)status;
    MPI_Bcast(&encoded, 1, MPI_INT, 0, MPI_COMM_WORLD);
    status = (MVMCPowerLanczosObservableCensusStatus)encoded;
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
      if (wire_size > (size_t)INT_MAX ||
          MPI_Bcast(wire, (int)wire_size, MPI_BYTE, 0,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
        status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR;
      }
    }
  }
#endif
  if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK && world_rank != 0) {
    char diagnostic[256];
    status = mvmc_power_lanczos_observable_plan_unpack(
        wire, wire_size, plan, diagnostic, sizeof(diagnostic));
  }
  free(wire);
  local_ok = status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
#ifdef _mpi_use
  {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    local_ok = global_ok;
  }
#endif
  if (!local_ok) {
    mvmc_power_lanczos_observable_plan_destroy(plan);
  }
  return local_ok;
}

static MVMCKrylovBoundedLimits production_limits(
    uint64_t generation_hash);

static MVMCPowerLanczosProductionExecution session_execution(
    int step, const MVMCKrylovFockModel *model,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t generation_hash,
    MVMCKrylovBoundedCommunicator chain_communicator) {
  MVMCPowerLanczosProductionExecution execution;
  const size_t basis_count = (size_t)step + 1;
  size_t index;
  memset(&execution, 0, sizeof(execution));
  execution.world_communicator = world_communicator();
  execution.chain_communicator = chain_communicator;
  execution.proposal_model = model;
  execution.amplitude = amplitude;
  execution.amplitude_context = amplitude_context;
  execution.amplitude_generation_hash = generation_hash;
  execution.bounded_limits = production_limits(generation_hash);
  execution.bounded_limits.max_order = (int)basis_count;
  execution.coefficient_guide_policy.order = (int)basis_count;
  execution.coefficient_guide_policy.eta = 0.25;
  execution.coefficient_guide_policy.policy_hash =
      generation_hash ^ UINT64_C(0x5036444755494445) ^ (uint64_t)step;
  execution.matrix_policy.order = (int)basis_count;
  execution.matrix_policy.eta = 0.25;
  for (index = 0; index <= basis_count; ++index) {
    execution.coefficient_guide_policy.lambda[index] = 1.0;
    execution.matrix_policy.guide_lambda[index] = 1.0;
    execution.matrix_policy.target_weight[index] = 1.0;
  }
  CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
            1, 2, &execution.proposal_policy) == MVMC_KRYLOV_STATUS_OK,
        "production session proposal policy");
  CHECK(mvmc_power_lanczos_gevp_default_policy(
            0x1p-40, &execution.gevp_policy) == MVMC_KRYLOV_GEVP_OK,
        "production session GEVP policy");
  execution.block_count = 4;
  execution.coefficient_warm_up = 16;
  execution.coefficient_sample_count = 512;
  execution.coefficient_interval = 1;
  execution.final_warm_up = 16;
  execution.final_sample_count = 1024;
  execution.final_interval = 1;
  execution.maximum_leave_one_projective_distance = 1.0;
  return execution;
}

static void check_observable_prepare_boundaries(
    const MVMCKrylovFockModel *model,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t generation_hash,
    MVMCKrylovBoundedCommunicator chain_communicator, size_t split_size,
    const MVMCPowerLanczosObservablePlan *observable_plan) {
  const uint64_t root = UINT64_C(5);
  int test_case;
  for (test_case = 0; test_case < 6; ++test_case) {
    const int step = test_case == 3 ? 2 : 1;
    MVMCPowerLanczosProductionSessionConfig config;
    MVMCPowerLanczosProductionExecution execution = session_execution(
        step, model, amplitude, amplitude_context, generation_hash,
        chain_communicator);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary summary;
    MVMCPowerLanczosObservablePlan corrupt_plan = *observable_plan;
    MVMCPowerLanczosObservableRecord corrupt_records[4];
    MVMCKrylovStatus create_status;
    MVMCKrylovStatus prepare_status;
    MVMCKrylovStatus summary_status;
    memcpy(corrupt_records, observable_plan->records,
           sizeof(corrupt_records));
    memset(&config, 0, sizeof(config));
    config.power_step = step;
    config.resolved_base_seed = UINT64_C(0x5036434c41535349);
    config.run_index = UINT64_C(10);
    config.mpi_world_rank = (size_t)world_rank;
    config.mpi_world_size = (size_t)world_size;
    config.split_size = split_size;
    config.base_generation = UINT64_C(22);
    config.configuration_bit_count = 4;
    execution.block_count = 32;
    execution.coefficient_sample_count = 512;
    execution.final_sample_count = 512;
    execution.observable_enabled = 1;
    execution.observable_layout.nsite = 2;
    execution.observable_layout.up_electron_count = 1;
    execution.observable_layout.down_electron_count = 1;
    execution.observable_layout.pure_spin = model->pure_spin;
    execution.observable_plan = observable_plan;
    if (test_case == 0) {
      execution.observable_plan = NULL;
    } else if (test_case == 1) {
      execution.block_count = 31;
      execution.coefficient_sample_count = 496;
      execution.final_sample_count = 496;
    } else if (test_case == 2) {
      execution.final_sample_count = 1024;
    } else if (test_case == 4) {
      corrupt_records[0].raw_indices[0] = 1;
      corrupt_plan.records = corrupt_records;
      execution.observable_plan = &corrupt_plan;
    } else if (test_case == 5) {
      execution.observable_enabled = 0;
      execution.observable_plan = NULL;
      execution.matrix_policy.guide_lambda[0] = 2.0;
    }
    create_status = create_anchor_session(
        &config, &root, 1, &session);
    prepare_status = create_status == MVMC_KRYLOV_STATUS_OK
                         ? mvmc_power_lanczos_production_session_prepare_execution(
                               session, &execution)
                         : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    summary_status =
        mvmc_power_lanczos_production_session_summary(session, &summary);
    CHECK(create_status == MVMC_KRYLOV_STATUS_OK &&
              prepare_status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              summary_status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !summary.valid &&
              summary.state ==
                  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FAILED,
          "observable prepare boundary case %d create=%s prepare=%s "
          "summary=%s",
          test_case, mvmc_krylov_status_string(create_status),
          mvmc_krylov_status_string(prepare_status),
          mvmc_krylov_status_string(summary_status));
    mvmc_power_lanczos_production_session_destroy(session);
  }
}

static void run_session_anchor(
    const MVMCKrylovFockModel *model,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t generation_hash, MVMCKrylovBoundedCommunicator chain_communicator,
    size_t split_size) {
  const uint64_t root = UINT64_C(5);
  MVMCPowerLanczosObservablePlan observable_plan;
  char expected_raw_sha256[MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY];
  char expected_semantic_sha256[
      MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY];
  char expected_operator_id[4][
      MVMC_POWER_LANCZOS_OBSERVABLE_OPERATOR_ID_CAPACITY];
  int step;
  if (!build_observable_plan(&observable_plan)) {
    CHECK(0, "observable plan build/broadcast");
    return;
  }
  memcpy(expected_raw_sha256, observable_plan.raw_observable_census_sha256,
         sizeof(expected_raw_sha256));
  memcpy(expected_semantic_sha256,
         observable_plan.semantic_observable_census_sha256,
         sizeof(expected_semantic_sha256));
  for (step = 0; step < 4; ++step) {
    memcpy(expected_operator_id[step],
           observable_plan.records[step].canonical_operator_id,
           sizeof(expected_operator_id[step]));
  }
  check_observable_prepare_boundaries(
      model, amplitude, amplitude_context, generation_hash, chain_communicator,
      split_size, &observable_plan);
  for (step = 1; step <= 2; ++step) {
    MVMCPowerLanczosProductionSessionConfig config;
    MVMCPowerLanczosProductionExecution execution = session_execution(
        step, model, amplitude, amplitude_context, generation_hash,
        chain_communicator);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary summary;
    MVMCKrylovStatus create_status;
    MVMCKrylovStatus prepare_status;
    MVMCKrylovStatus coefficient_status;
    MVMCKrylovStatus freeze_status;
    MVMCKrylovStatus final_status;
    MVMCKrylovStatus finalize_status;
    MVMCKrylovStatus summary_status;
    const size_t expected_coefficient_samples =
        step == 1 ? 512u : execution.coefficient_sample_count;
    const size_t expected_final_samples =
        step == 1 ? 512u : execution.final_sample_count;
    if (step == 1) {
      execution.block_count = 32;
      execution.coefficient_sample_count = expected_coefficient_samples;
      execution.final_sample_count = expected_final_samples;
      execution.observable_enabled = 1;
      execution.observable_layout.nsite = 2;
      execution.observable_layout.up_electron_count = 1;
      execution.observable_layout.down_electron_count = 1;
      execution.observable_layout.pure_spin = model->pure_spin;
      execution.observable_plan = &observable_plan;
    }
    memset(&config, 0, sizeof(config));
    config.power_step = step;
    config.resolved_base_seed = UINT64_C(0x5036434c41535349);
    config.run_index = UINT64_C(11);
    config.mpi_world_rank = (size_t)world_rank;
    config.mpi_world_size = (size_t)world_size;
    config.split_size = split_size;
    config.base_generation = UINT64_C(23);
    config.configuration_bit_count = 4;
    create_status = create_anchor_session(
        &config, &root, 1, &session);
    prepare_status = create_status == MVMC_KRYLOV_STATUS_OK
                         ? mvmc_power_lanczos_production_session_prepare_execution(
                               session, &execution)
                         : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    if (step == 1 && prepare_status == MVMC_KRYLOV_STATUS_OK) {
      MVMCPowerLanczosProductionObservableResult premature_result;
      MVMCPowerLanczosProductionObservableBlockResult
          premature_block_result;
      MVMCPowerLanczosProductionMatrixResult premature_matrix;
      MVMCPowerLanczosProductionCoefficientBlockResult
          premature_coefficient_block;
      MVMCPowerLanczosProductionFinalBlockResult premature_final_block;
      memset(&premature_result, 0xff, sizeof(premature_result));
      memset(&premature_block_result, 0xff,
             sizeof(premature_block_result));
      memset(&premature_matrix, 0xff, sizeof(premature_matrix));
      memset(&premature_coefficient_block, 0xff,
             sizeof(premature_coefficient_block));
      memset(&premature_final_block, 0xff,
             sizeof(premature_final_block));
      CHECK(mvmc_power_lanczos_production_session_observable_result(
                session, 0, &premature_result) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                !premature_result.valid,
            "classic observable getter before finalization");
      CHECK(mvmc_power_lanczos_production_session_observable_block_result(
                session, 0, 0, &premature_block_result) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                !premature_block_result.valid,
            "classic observable block getter before finalization");
      CHECK(mvmc_power_lanczos_production_session_matrix_result(
                session, &premature_matrix) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                !premature_matrix.valid,
            "classic matrix getter before finalization");
      CHECK(mvmc_power_lanczos_production_session_coefficient_block_result(
                session, 0, &premature_coefficient_block) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                !premature_coefficient_block.valid,
            "classic coefficient block getter before finalization");
      CHECK(mvmc_power_lanczos_production_session_final_block_result(
                session, 0, &premature_final_block) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                !premature_final_block.valid,
            "classic final block getter before finalization");
      observable_plan.records[0].raw_indices[0] = 1;
      observable_plan.records[0].canonical_operator_id[0] = 'X';
    }
    coefficient_status = prepare_status == MVMC_KRYLOV_STATUS_OK
                             ? mvmc_power_lanczos_production_session_run_coefficient(
                                   session)
                             : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    freeze_status = coefficient_status == MVMC_KRYLOV_STATUS_OK
                        ? mvmc_power_lanczos_production_session_freeze_coefficient(
                              session)
                        : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    final_status = freeze_status == MVMC_KRYLOV_STATUS_OK
                       ? mvmc_power_lanczos_production_session_run_final(
                             session)
                       : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    finalize_status = final_status == MVMC_KRYLOV_STATUS_OK
                          ? mvmc_power_lanczos_production_session_finalize(
                                session)
                          : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    summary_status =
        mvmc_power_lanczos_production_session_summary(session, &summary);
    CHECK(create_status == MVMC_KRYLOV_STATUS_OK &&
              prepare_status == MVMC_KRYLOV_STATUS_OK &&
              coefficient_status == MVMC_KRYLOV_STATUS_OK &&
              freeze_status == MVMC_KRYLOV_STATUS_OK &&
              final_status == MVMC_KRYLOV_STATUS_OK &&
              finalize_status == MVMC_KRYLOV_STATUS_OK &&
              summary_status == MVMC_KRYLOV_STATUS_OK &&
              summary.valid &&
              summary.state ==
                  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED &&
              summary.coefficient_gevp.valid &&
              summary.coefficient_gevp.retained_rank == step + 1 &&
              summary.coefficient_sample_count ==
                  (uint64_t)expected_coefficient_samples *
                      (uint64_t)summary.chain_count &&
              summary.final_sample_count ==
                  (uint64_t)expected_final_samples *
                      (uint64_t)summary.chain_count &&
              summary.execution_fingerprint != 0 &&
              summary.coefficient_provenance_hash != 0 &&
              summary.final_policy_hash != 0 &&
              isfinite(summary.final_energy) &&
              isfinite(summary.final_energy_imaginary) &&
              isfinite(summary.corrected_variance) &&
              isfinite(summary.corrected_variance_imaginary),
          "classic production session lifecycle step %d "
          "create=%s prepare=%s coefficient=%s freeze=%s final=%s "
          "finalize=%s summary=%s state=%s rank=%d E=%.17g var=%.17g",
          step, mvmc_krylov_status_string(create_status),
          mvmc_krylov_status_string(prepare_status),
          mvmc_krylov_status_string(coefficient_status),
          mvmc_krylov_status_string(freeze_status),
          mvmc_krylov_status_string(final_status),
          mvmc_krylov_status_string(finalize_status),
          mvmc_krylov_status_string(summary_status),
          mvmc_power_lanczos_production_session_state_string(summary.state),
          summary.coefficient_gevp.retained_rank, summary.final_energy,
          summary.corrected_variance);
    if (finalize_status == MVMC_KRYLOV_STATUS_OK) {
      MVMCPowerLanczosProductionMatrixResult matrix;
      double complex overlap_sum[6] = {0};
      double complex hamiltonian_sum[6] = {0};
      double complex adjoint_sum[6] = {0};
      double complex squared_sum[6] = {0};
      double complex energy_sum = 0.0;
      double complex diagnostic_sum = 0.0;
      uint64_t coefficient_count = 0;
      uint64_t final_count = 0;
      size_t block;
      size_t entry;
      CHECK(mvmc_power_lanczos_production_session_matrix_result(
                session, &matrix) == MVMC_KRYLOV_STATUS_OK &&
                matrix.valid && matrix.basis_count == (size_t)step + 1 &&
                matrix.upper_count ==
                    ((size_t)step + 1) * ((size_t)step + 2) / 2 &&
                matrix.sample_count == summary.coefficient_sample_count,
            "classic immutable matrix evidence step %d", step);
      for (block = 0; block < execution.block_count; ++block) {
        MVMCPowerLanczosProductionCoefficientBlockResult coefficient_block;
        MVMCPowerLanczosProductionFinalBlockResult final_block;
        CHECK(mvmc_power_lanczos_production_session_coefficient_block_result(
                  session, block, &coefficient_block) ==
                      MVMC_KRYLOV_STATUS_OK &&
                  coefficient_block.valid &&
                  coefficient_block.block_index == block &&
                  coefficient_block.basis_count == matrix.basis_count &&
                  coefficient_block.upper_count == matrix.upper_count &&
                  coefficient_block.sample_count > 0 &&
                  coefficient_block.leave_one_sample_count > 0 &&
                  coefficient_block.sample_count +
                          coefficient_block.leave_one_sample_count ==
                      summary.coefficient_sample_count &&
                  isfinite(coefficient_block.leave_one_energy) &&
                  isfinite(creal(
                      coefficient_block.leave_one_second_moment)) &&
                  isfinite(cimag(
                      coefficient_block.leave_one_second_moment)) &&
                  coefficient_block.leave_one_projective_distance >= 0.0 &&
                  coefficient_block.leave_one_projective_distance <=
                      execution.maximum_leave_one_projective_distance,
              "classic immutable coefficient block %zu step %d", block,
              step);
        CHECK(mvmc_power_lanczos_production_session_final_block_result(
                  session, block, &final_block) == MVMC_KRYLOV_STATUS_OK &&
                  final_block.valid && final_block.block_index == block &&
                  final_block.sample_count > 0,
              "classic immutable final block %zu step %d", block, step);
        coefficient_count += coefficient_block.sample_count;
        final_count += final_block.sample_count;
        energy_sum += final_block.energy_sum;
        diagnostic_sum += final_block.local_energy_abs_squared_sum;
        for (entry = 0; entry < matrix.upper_count; ++entry) {
          overlap_sum[entry] += coefficient_block.overlap_sum[entry];
          hamiltonian_sum[entry] +=
              coefficient_block.hamiltonian_sum[entry];
          adjoint_sum[entry] +=
              coefficient_block.hamiltonian_adjoint_sum[entry];
          squared_sum[entry] +=
              coefficient_block.hamiltonian_squared_sum[entry];
        }
      }
      CHECK(coefficient_count == summary.coefficient_sample_count &&
                final_count == summary.final_sample_count,
            "classic immutable evidence sample census step %d", step);
      for (entry = 0; entry < matrix.upper_count; ++entry) {
        CHECK(cabs(matrix.overlap[entry] -
                       overlap_sum[entry] / (double)coefficient_count) <=
                      1.0e-12 &&
                  cabs(matrix.hamiltonian[entry] -
                       hamiltonian_sum[entry] /
                           (double)coefficient_count) <= 1.0e-12 &&
                  cabs(matrix.hamiltonian_adjoint[entry] -
                       adjoint_sum[entry] / (double)coefficient_count) <=
                      1.0e-12 &&
                  cabs(matrix.hamiltonian_squared[entry] -
                       squared_sum[entry] / (double)coefficient_count) <=
                      1.0e-12,
              "classic matrix reconstruction entry %zu step %d", entry,
              step);
      }
      CHECK(cabs(energy_sum / (double)final_count -
                     (summary.final_energy +
                      I * summary.final_energy_imaginary)) <= 1.0e-12 &&
                cabs(diagnostic_sum / (double)final_count -
                     summary.final_local_energy_abs_squared_diagnostic) <=
                    1.0e-12,
            "classic final block reconstruction step %d", step);
      {
        MVMCPowerLanczosProductionCoefficientBlockResult invalid_coefficient;
        MVMCPowerLanczosProductionFinalBlockResult invalid_final;
        memset(&invalid_coefficient, 0xff, sizeof(invalid_coefficient));
        memset(&invalid_final, 0xff, sizeof(invalid_final));
        CHECK(mvmc_power_lanczos_production_session_coefficient_block_result(
                  session, execution.block_count, &invalid_coefficient) ==
                      MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                  !invalid_coefficient.valid,
              "classic coefficient evidence range boundary");
        CHECK(mvmc_power_lanczos_production_session_final_block_result(
                  session, execution.block_count, &invalid_final) ==
                      MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                  !invalid_final.valid,
              "classic final evidence range boundary");
      }
    }
    if (step == 1 && finalize_status == MVMC_KRYLOV_STATUS_OK) {
      enum {
        observable_count = 4,
        observable_matrix_entries = 4,
        coefficient_primitive_count = 29,
        final_primitive_count = 6
      };
      size_t block;
      size_t request;
      double complex coefficient_sum = 0.0;
      double complex final_sum = 0.0;
      CHECK(summary.observable_enabled &&
                summary.observable_request_count == 4 &&
                memcmp(summary.observable_census_sha256,
                       expected_raw_sha256,
                       sizeof(expected_raw_sha256)) == 0 &&
                memcmp(summary.observable_semantic_sha256,
                       expected_semantic_sha256,
                       sizeof(expected_semantic_sha256)) == 0,
            "classic observable summary/deep-copy");
      for (request = 0; request < 4; ++request) {
        MVMCPowerLanczosProductionObservableResult observable;
        CHECK(mvmc_power_lanczos_production_session_observable_result(
                  session, request, &observable) ==
                      MVMC_KRYLOV_STATUS_OK &&
                  observable.valid &&
                  observable.request_index == request &&
                  strcmp(observable.canonical_operator_id,
                         expected_operator_id[request]) == 0 &&
                  isfinite(creal(observable.coefficient_estimate)) &&
                  isfinite(cimag(observable.coefficient_estimate)) &&
                  isfinite(creal(observable.final_estimate)) &&
                  isfinite(cimag(observable.final_estimate)),
              "classic observable request %zu", request);
        coefficient_sum += observable.coefficient_estimate;
        final_sum += observable.final_estimate;
      }
      CHECK(cabs(coefficient_sum - 2.0) <= 1.0e-10 &&
                cabs(final_sum - 2.0) <= 1.0e-10,
            "classic particle-census identity coefficient=(%.17g,%.17g) "
            "final=(%.17g,%.17g)",
            creal(coefficient_sum), cimag(coefficient_sum),
            creal(final_sum), cimag(final_sum));
      {
        MVMCPowerLanczosPrimitiveTraceSummary coefficient_trace;
        MVMCPowerLanczosPrimitiveTraceSummary final_trace;
        const MVMCKrylovStatus coefficient_trace_status =
            mvmc_power_lanczos_production_session_coefficient_trace_summary(
                session, &coefficient_trace);
        const MVMCKrylovStatus final_trace_status =
            mvmc_power_lanczos_production_session_final_trace_summary(
                session, &final_trace);
        if (world_rank == 0) {
          MVMCPowerLanczosPrimitiveTraceGroup coefficient_group;
          MVMCPowerLanczosPrimitiveTraceGroup final_group;
          double complex coefficient_values[coefficient_primitive_count];
          double coefficient_bounds[coefficient_primitive_count];
          uint8_t coefficient_support[coefficient_primitive_count];
          double complex final_values[final_primitive_count];
          double final_bounds[final_primitive_count];
          uint8_t final_support[final_primitive_count];
          const size_t last_group = summary.chain_count - 1;
          double coefficient_tail = NAN;
          double final_tail = NAN;
          size_t primitive;
          CHECK(coefficient_trace_status == MVMC_KRYLOV_STATUS_OK &&
                    coefficient_trace.valid && coefficient_trace.frozen &&
                    coefficient_trace.stage ==
                        MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_COEFFICIENT &&
                    coefficient_trace.group_count == summary.chain_count &&
                    coefficient_trace.samples_per_group ==
                        expected_coefficient_samples &&
                    coefficient_trace.appended_sample_count ==
                        expected_coefficient_samples * summary.chain_count &&
                    coefficient_trace.primitive_count ==
                        coefficient_primitive_count &&
                    coefficient_trace.scalar_count ==
                        2 * coefficient_primitive_count,
                "classic observable coefficient trace summary");
          CHECK(final_trace_status == MVMC_KRYLOV_STATUS_OK &&
                    final_trace.valid && final_trace.frozen &&
                    final_trace.stage ==
                        MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_FINAL &&
                    final_trace.group_count == summary.chain_count &&
                    final_trace.samples_per_group == expected_final_samples &&
                    final_trace.appended_sample_count ==
                        expected_final_samples * summary.chain_count &&
                    final_trace.primitive_count == final_primitive_count &&
                    final_trace.scalar_count ==
                        2 * final_primitive_count,
                "classic observable final trace summary");
          CHECK(mvmc_power_lanczos_production_session_coefficient_trace_group(
                    session, last_group, &coefficient_group) ==
                        MVMC_KRYLOV_STATUS_OK &&
                    coefficient_group.valid &&
                    coefficient_group.group_ordinal == last_group &&
                    coefficient_group.leader_world_rank ==
                        last_group * split_size &&
                    coefficient_group.chain_size ==
                        (size_t)world_size - last_group * split_size &&
                    coefficient_group.sample_count ==
                        expected_coefficient_samples,
                "classic observable coefficient trace group identity");
          CHECK(mvmc_power_lanczos_production_session_final_trace_group(
                    session, last_group, &final_group) ==
                        MVMC_KRYLOV_STATUS_OK &&
                    final_group.valid &&
                    final_group.group_ordinal == last_group &&
                    final_group.leader_world_rank == last_group * split_size &&
                    final_group.chain_size ==
                        (size_t)world_size - last_group * split_size &&
                    final_group.sample_count == expected_final_samples,
                "classic observable final trace group identity");
          CHECK(mvmc_power_lanczos_production_session_coefficient_trace_sample(
                    session, last_group, expected_coefficient_samples - 1,
                    coefficient_values, coefficient_primitive_count,
                    coefficient_bounds, coefficient_primitive_count,
                    coefficient_support, coefficient_primitive_count,
                    &coefficient_tail) == MVMC_KRYLOV_STATUS_OK &&
                    coefficient_tail == 0.0,
                "classic observable coefficient trace sample");
          CHECK(mvmc_power_lanczos_production_session_final_trace_sample(
                    session, last_group, expected_final_samples - 1,
                    final_values, final_primitive_count, final_bounds,
                    final_primitive_count, final_support,
                    final_primitive_count, &final_tail) ==
                        MVMC_KRYLOV_STATUS_OK &&
                    final_values[0] == 1.0 && final_bounds[0] == 0.0 &&
                    final_support[0] ==
                        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO &&
                    final_tail == hypot(creal(final_values[1]),
                                        cimag(final_values[1])),
                "classic observable final trace sample");
          for (primitive = 0; primitive < coefficient_primitive_count;
               ++primitive) {
            CHECK(isfinite(creal(coefficient_values[primitive])) &&
                      isfinite(cimag(coefficient_values[primitive])) &&
                      isfinite(coefficient_bounds[primitive]) &&
                      coefficient_bounds[primitive] >= 0.0 &&
                      primitive_support_is_valid(
                          coefficient_support[primitive]),
                  "classic coefficient primitive %zu", primitive);
          }
          for (primitive = 0; primitive < final_primitive_count;
               ++primitive) {
            CHECK(isfinite(creal(final_values[primitive])) &&
                      isfinite(cimag(final_values[primitive])) &&
                      isfinite(final_bounds[primitive]) &&
                      final_bounds[primitive] >= 0.0 &&
                      primitive_support_is_valid(final_support[primitive]),
                  "classic final primitive %zu", primitive);
          }
          for (primitive = 0; primitive < observable_matrix_entries;
               ++primitive) {
            double complex observable_sum = 0.0;
            const size_t observable_offset = 1 + 4 * summary.upper_count;
            const double complex overlap =
                primitive == 0
                    ? coefficient_values[1]
                    : (primitive == 1
                           ? coefficient_values[2]
                           : (primitive == 2
                                  ? conj(coefficient_values[2])
                                  : coefficient_values[3]));
            for (request = 0; request < observable_count; ++request) {
              observable_sum +=
                  coefficient_values[observable_offset +
                                     request * observable_matrix_entries +
                                     primitive];
            }
            CHECK(cabs(observable_sum - 2.0 * overlap) <= 1.0e-10,
                  "classic coefficient trace observable census entry %zu",
                  primitive);
          }
          {
            double complex observable_sum = 0.0;
            for (request = 0; request < observable_count; ++request) {
              observable_sum += final_values[2 + request];
            }
            CHECK(cabs(observable_sum - 2.0) <= 1.0e-10,
                  "classic final trace observable census");
          }
        } else {
          CHECK(coefficient_trace_status ==
                        MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                    final_trace_status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "classic observable trace must be world-root-owned");
        }
      }
      {
        MVMCPowerLanczosProductionObservableResult invalid_result;
        MVMCPowerLanczosProductionObservableBlockResult
            invalid_block_result;
        memset(&invalid_result, 0xff, sizeof(invalid_result));
        memset(&invalid_block_result, 0xff, sizeof(invalid_block_result));
        CHECK(mvmc_power_lanczos_production_session_observable_result(
                  session, 4, &invalid_result) ==
                      MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                  !invalid_result.valid,
              "classic observable request range boundary");
        CHECK(mvmc_power_lanczos_production_session_observable_block_result(
                  session, execution.block_count, 0,
                  &invalid_block_result) ==
                      MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                  !invalid_block_result.valid,
              "classic observable block range boundary");
      }
      for (block = 0; block < execution.block_count; ++block) {
        double complex leave_one_sum = 0.0;
        double complex final_block_sum = 0.0;
        for (request = 0; request < 4; ++request) {
          MVMCPowerLanczosProductionObservableBlockResult block_result;
          CHECK(mvmc_power_lanczos_production_session_observable_block_result(
                    session, block, request, &block_result) ==
                        MVMC_KRYLOV_STATUS_OK &&
                    block_result.valid &&
                    block_result.block_index == block &&
                    block_result.request_index == request &&
                    isfinite(creal(
                        block_result.coefficient_leave_one_estimate)) &&
                    isfinite(cimag(
                        block_result.coefficient_leave_one_estimate)) &&
                    isfinite(creal(block_result.final_block_mean)) &&
                    isfinite(cimag(block_result.final_block_mean)) &&
                    block_result.leave_one_projective_distance >= 0.0 &&
                    block_result.leave_one_projective_distance <= 1.0,
                "classic observable block %zu request %zu", block,
                request);
          leave_one_sum += block_result.coefficient_leave_one_estimate;
          final_block_sum += block_result.final_block_mean;
        }
        CHECK(cabs(leave_one_sum - 2.0) <= 1.0e-10 &&
                  cabs(final_block_sum - 2.0) <= 1.0e-10,
              "classic block particle-census identity block=%zu "
              "coefficient=(%.17g,%.17g) final=(%.17g,%.17g)",
              block, creal(leave_one_sum), cimag(leave_one_sum),
              creal(final_block_sum), cimag(final_block_sum));
      }
    }
    if (step == 2 && world_size == 1) {
      MVMCPowerLanczosProductionObservableResult disabled_result;
      memset(&disabled_result, 0xff, sizeof(disabled_result));
      CHECK(mvmc_power_lanczos_production_session_observable_result(
                session, 0, &disabled_result) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                !disabled_result.valid,
            "classic disabled observable getter");
      CHECK(fabs(summary.final_energy - summary.coefficient_gevp.energy) <=
                    1.0e-3 &&
                fabs(summary.corrected_variance) <= 1.0e-3,
            "classic p2 final state is inconsistent with the frozen GEVP: "
            "gevp=%.17g final=%.17g variance=%.17g",
            summary.coefficient_gevp.energy, summary.final_energy,
            summary.corrected_variance);
    }
    mvmc_power_lanczos_production_session_destroy(session);
  }
  mvmc_power_lanczos_observable_plan_destroy(&observable_plan);
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

static int scaled_equal(const MVMCScaledComplex *left,
                        const MVMCScaledComplex *right) {
  return left->state == right->state &&
         creal(left->phase) == creal(right->phase) &&
         cimag(left->phase) == cimag(right->phase) &&
         left->log_abs == right->log_abs &&
         left->log_abs_error_bound == right->log_abs_error_bound &&
         left->max_input_log_abs == right->max_input_log_abs &&
         left->cancellation_log_abs == right->cancellation_log_abs &&
         left->cancellation_ratio == right->cancellation_ratio;
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

static void check_corrected_dispatch(
    const MVMCPowerLanczosClassicView *view,
    MVMCKrylovBoundedCommunicator chain_communicator, size_t split_size,
    int power_step, MVMCPowerLanczosCorrectedBootstrapMode bootstrap_mode,
    const char *label) {
  const uint64_t root =
      view != NULL && view->pure_spin ? UINT64_C(9) : UINT64_C(5);
  static const MVMCPowerLanczosStabilizationOutputIdentity identity = {
      "classic-dispatch-anchor",
      "0123456789012345678901234567890123456789",
      "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
      "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
      "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
      "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
      "focused-anchor",
      "STABILIZATION_SEED_1"};
  MVMCPowerLanczosCorrectedDispatchInput input;
  MVMCPowerLanczosCorrectedDispatchResult result;
  MVMCKrylovStatus status;
  char output_directory[] = "./power_lanczos_dispatch_XXXXXX";
  char output_path[256];
  int path_length = -1;
  int directory_ready = 1;
  if (world_rank == 0 && mkdtemp(output_directory) == NULL) {
    directory_ready = 0;
  }
#ifdef _mpi_use
  MPI_Bcast(&directory_ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(output_directory, (int)sizeof(output_directory), MPI_CHAR, 0,
            MPI_COMM_WORLD);
#endif
  CHECK(directory_ready, "%s corrected dispatch output directory", label);
  if (!directory_ready) return;
  memset(&input, 0, sizeof(input));
  input.classic_view = view;
  input.world_communicator = world_communicator();
  input.chain_communicator = chain_communicator;
  input.power_step = power_step;
  input.resolved_base_seed = UINT64_C(0x53544142494c3031);
  input.run_index = UINT64_C(1);
  input.mpi_world_rank = (size_t)world_rank;
  input.mpi_world_size = (size_t)world_size;
  input.split_size = split_size;
  input.base_generation = UINT64_C(23);
  input.bootstrap_mode = bootstrap_mode;
  if (bootstrap_mode ==
      MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL) {
    input.root_terminal_words = world_rank == 0 ? &root : NULL;
    input.root_terminal_word_count = world_rank == 0 ? 1u : 0u;
    input.terminal_proposal_counter = UINT64_C(4096);
  }
  input.controls.coefficient_warm_up = 16;
  input.controls.coefficient_sample_count = 64;
  input.controls.coefficient_interval = 1;
  input.controls.final_warm_up = 16;
  input.controls.final_sample_count = 64;
  input.controls.final_interval = 1;
  input.root_identity = world_rank == 0 ? &identity : NULL;
  input.root_output_directory = world_rank == 0 ? output_directory : NULL;
  input.root_output_basename =
      world_rank == 0 ? "stabilization.json" : NULL;
  status = mvmc_power_lanczos_corrected_dispatch_run(&input, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            result.bootstrap_mode == bootstrap_mode &&
            result.bootstrap_proposal_counter ==
                (bootstrap_mode ==
                         MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_UNIFORM_SECTOR
                     ? UINT64_C(1)
                     : UINT64_C(4096)) &&
            result.output_status ==
                MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK &&
            result.output_size > 0 && strlen(result.output_sha256) == 64 &&
            isfinite(result.energy) && isfinite(result.variance),
        "%s corrected dispatch status=%s decision=%s output=%s", label,
        mvmc_krylov_status_string(status),
        mvmc_power_lanczos_stabilization_decision_string(result.decision),
        result.output_sha256);
  if (world_rank == 0) {
    path_length = snprintf(output_path, sizeof(output_path), "%s/%s",
                           output_directory, "stabilization.json");
    CHECK(path_length > 0 && (size_t)path_length < sizeof(output_path) &&
              access(output_path, F_OK) == 0,
          "%s corrected dispatch atomic output exists", label);
    if (path_length > 0 && (size_t)path_length < sizeof(output_path)) {
      (void)unlink(output_path);
    }
    (void)rmdir(output_directory);
  }
}

static void run_production_anchor(void) {
  enum { qp_total = 8, singular_qp = 3, nproj = 20 };
  int transfer_rows[4][4] = {
      {1, 0, 0, 0},
      {0, 0, 1, 0},
      {1, 1, 0, 1},
      {0, 1, 1, 1},
  };
  int transfer_rows_snapshot[4][4];
  int *transfer_indices[4] = {
      transfer_rows[0], transfer_rows[1], transfer_rows[2],
      transfer_rows[3]};
  double complex transfer_parameters[4] = {0.75, 0.75, 0.50, 0.50};
  double complex transfer_parameters_snapshot[4];
  int intra_indices[1] = {0};
  int intra_indices_snapshot[1];
  double intra_parameters[1] = {0.40};
  double intra_parameters_snapshot[1];
  MVMCPowerLanczosClassicView view;
  MVMCPowerLanczosClassicBridge *bridge = NULL;
  MVMCPowerLanczosClassicBridge *complex_bridge = NULL;
  MVMCPowerLanczosClassicBridgeSummary bridge_summary;
  const MVMCKrylovFockModel *model = NULL;
  double slater[qp_total * 16];
  double slater_snapshot[qp_total * 16];
  double complex complex_slater[qp_total * 16];
  double complex weights[qp_total];
  double complex weights_snapshot[qp_total];
  double complex parameters[nproj];
  double complex parameters_snapshot[nproj];
  int gutzwiller[2] = {0, 1};
  int gutzwiller_snapshot[2];
  int jastrow_row_0[2] = {-1, 0};
  int jastrow_row_1[2] = {0, -1};
  int *jastrow[2] = {jastrow_row_0, jastrow_row_1};
  int spin_row_0[2] = {-1, 0};
  int spin_row_1[2] = {0, -1};
  int *spin_jastrow[2] = {spin_row_0, spin_row_1};
  int dh2_row[4] = {1, 1, 0, 0};
  int *dh2[1] = {dh2_row};
  int dh4_row[8] = {1, 1, 1, 1, 0, 0, 0, 0};
  int *dh4[1] = {dh4_row};
  int binding_snapshot[20];
  MVMCKrylovScaledAmplitudeCallback amplitude = NULL;
  void *amplitude_context = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *bounded_workspace = NULL;
  MVMCKrylovBoundedResult result;
  MVMCKrylovBoundedResult immutable_result;
  uint64_t root = UINT64_C(5);
  const uint64_t root_snapshot = root;
  uint64_t generation_hash = 0;
  uint64_t plan_hash = 0;
  MVMCKrylovStatus status;
  int index;
  const size_t split_size = world_size > 1 ? 2u : 1u;
  int chain_rank = 0;
  int chain_size = 1;
#ifdef _mpi_use
  MPI_Comm chain_communicator = MPI_COMM_NULL;
  CHECK(MPI_Comm_split(MPI_COMM_WORLD,
                       world_rank / (int)split_size, world_rank,
                       &chain_communicator) == MPI_SUCCESS,
        "production session chain communicator");
  MPI_Comm_rank(chain_communicator, &chain_rank);
  MPI_Comm_size(chain_communicator, &chain_size);
#else
  MVMCKrylovBoundedCommunicator chain_communicator = 0;
#endif

  fill_slater(slater, qp_total, singular_qp);
  for (index = 0; index < qp_total; ++index) {
    weights[index] = 1.0 / (double)(index + 1);
  }
  for (index = 0; index < qp_total * 16; ++index) {
    complex_slater[index] = slater[index];
  }
  for (index = 0; index < nproj; ++index) {
    parameters[index] = (double)(index + 1) / 100.0;
  }
  memcpy(slater_snapshot, slater, sizeof(slater));
  memcpy(weights_snapshot, weights, sizeof(weights));
  memcpy(parameters_snapshot, parameters, sizeof(parameters));
  memcpy(transfer_rows_snapshot, transfer_rows, sizeof(transfer_rows));
  memcpy(transfer_parameters_snapshot, transfer_parameters,
         sizeof(transfer_parameters));
  memcpy(intra_indices_snapshot, intra_indices, sizeof(intra_indices));
  memcpy(intra_parameters_snapshot, intra_parameters,
         sizeof(intra_parameters));
  memcpy(gutzwiller_snapshot, gutzwiller, sizeof(gutzwiller));
  memcpy(binding_snapshot, jastrow_row_0, sizeof(jastrow_row_0));
  memcpy(binding_snapshot + 2, jastrow_row_1, sizeof(jastrow_row_1));
  memcpy(binding_snapshot + 4, spin_row_0, sizeof(spin_row_0));
  memcpy(binding_snapshot + 6, spin_row_1, sizeof(spin_row_1));
  memcpy(binding_snapshot + 8, dh2_row, sizeof(dh2_row));
  memcpy(binding_snapshot + 12, dh4_row, sizeof(dh4_row));

  memset(&view, 0, sizeof(view));
  view.site_count = 2;
  view.up_electron_count = 1;
  view.down_electron_count = 1;
  view.arithmetic = MVMC_POWER_LANCZOS_CLASSIC_REAL;
  view.transfer_count = 4;
  view.transfer_indices = transfer_indices;
  view.transfer_parameters = transfer_parameters;
  view.coulomb_intra_count = 1;
  view.coulomb_intra_indices = intra_indices;
  view.coulomb_intra_parameters = intra_parameters;
  view.qp_total = qp_total;
  view.qp_start = qp_total * chain_rank / chain_size;
  view.qp_end = qp_total * (chain_rank + 1) / chain_size;
  view.nproj = nproj;
  view.ngutzwiller_idx = 2;
  view.njastrow_idx = 1;
  view.nspin_jastrow_idx = 1;
  view.ndoublon_holon_2site_idx = 1;
  view.ndoublon_holon_4site_idx = 1;
  view.gutzwiller_idx = gutzwiller;
  view.jastrow_idx = jastrow;
  view.spin_jastrow_idx = spin_jastrow;
  view.doublon_holon_2site_idx = dh2;
  view.doublon_holon_4site_idx = dh4;
  view.projection_parameters = parameters;
  view.qp_weights = weights;
  view.slater_real = slater;

  status = mvmc_power_lanczos_classic_bridge_create(
      &view, world_communicator(), chain_communicator, &bridge);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && bridge != NULL,
        "production classic bridge failed with status %d", (int)status);
  if (status != MVMC_KRYLOV_STATUS_OK || bridge == NULL) goto cleanup;
  model = mvmc_power_lanczos_classic_bridge_model(bridge);
  CHECK(model != NULL && model->term_count == 5 &&
            model->operator_count == 12,
        "production model expansion is unexpected");
  if (model == NULL) goto cleanup;
  amplitude = mvmc_power_lanczos_classic_bridge_amplitude(bridge);
  amplitude_context =
      mvmc_power_lanczos_classic_bridge_amplitude_context(bridge);
  generation_hash =
      mvmc_power_lanczos_classic_bridge_generation_hash(bridge);
  CHECK(generation_hash != 0,
        "production amplitude generation hash is zero");
  memset(&bridge_summary, 0, sizeof(bridge_summary));
  CHECK(mvmc_power_lanczos_classic_bridge_summary(
            bridge, &bridge_summary) == MVMC_KRYLOV_STATUS_OK &&
            bridge_summary.valid &&
            bridge_summary.version ==
                MVMC_POWER_LANCZOS_CLASSIC_BRIDGE_VERSION &&
            bridge_summary.arithmetic == MVMC_POWER_LANCZOS_CLASSIC_REAL &&
            bridge_summary.model_bytes > 0 &&
            bridge_summary.amplitude_bytes > 0 &&
            bridge_summary.allocated_bytes >=
                bridge_summary.model_bytes + bridge_summary.amplitude_bytes &&
            bridge_summary.term_count == 5 &&
            bridge_summary.operator_count == 12 &&
            bridge_summary.amplitude_generation_hash == generation_hash,
        "production classic bridge summary is invalid");

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
      bounded_workspace, &root, 1, amplitude, amplitude_context, &result);
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
            memcmp(transfer_rows, transfer_rows_snapshot,
                   sizeof(transfer_rows)) == 0 &&
            memcmp(transfer_parameters, transfer_parameters_snapshot,
                   sizeof(transfer_parameters)) == 0 &&
            memcmp(intra_indices, intra_indices_snapshot,
                   sizeof(intra_indices)) == 0 &&
            memcmp(intra_parameters, intra_parameters_snapshot,
                   sizeof(intra_parameters)) == 0 &&
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

  if (status == MVMC_KRYLOV_STATUS_OK && result.valid) {
    int values_equal = 1;
    immutable_result = result;
    transfer_rows[0][0] = 0;
    transfer_parameters[0] = -7.0 + 3.0 * I;
    intra_indices[0] = 1;
    intra_parameters[0] = -9.0;
    slater[2] += 17.0;
    weights[0] = -2.0 + 5.0 * I;
    parameters[0] += 3.0;
    gutzwiller[0] = 1;
    jastrow_row_0[0] = 0;
    spin_row_0[0] = 0;
    dh2_row[0] = 0;
    dh4_row[0] = 0;
    status = mvmc_bounded_krylov_evaluate(
        bounded_workspace, &root, 1, amplitude, amplitude_context, &result);
    for (index = 0; index <= 3; ++index) {
      values_equal =
          values_equal &&
          scaled_equal(&immutable_result.value[index], &result.value[index]);
    }
    CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid && values_equal &&
              result.statistics.trace_hash ==
                  immutable_result.statistics.trace_hash &&
              result.statistics.terminal_amplitude_calls ==
                  immutable_result.statistics.terminal_amplitude_calls &&
              mvmc_power_lanczos_classic_bridge_generation_hash(bridge) ==
                  generation_hash &&
              model->term_count == 5 && model->operator_count == 12,
          "classic bridge retained caller-owned production buffers");
    memcpy(transfer_rows, transfer_rows_snapshot, sizeof(transfer_rows));
    memcpy(transfer_parameters, transfer_parameters_snapshot,
           sizeof(transfer_parameters));
    memcpy(intra_indices, intra_indices_snapshot, sizeof(intra_indices));
    memcpy(intra_parameters, intra_parameters_snapshot,
           sizeof(intra_parameters));
    memcpy(slater, slater_snapshot, sizeof(slater));
    memcpy(weights, weights_snapshot, sizeof(weights));
    memcpy(parameters, parameters_snapshot, sizeof(parameters));
    memcpy(gutzwiller, gutzwiller_snapshot, sizeof(gutzwiller));
    memcpy(jastrow_row_0, binding_snapshot, sizeof(jastrow_row_0));
    memcpy(jastrow_row_1, binding_snapshot + 2, sizeof(jastrow_row_1));
    memcpy(spin_row_0, binding_snapshot + 4, sizeof(spin_row_0));
    memcpy(spin_row_1, binding_snapshot + 6, sizeof(spin_row_1));
    memcpy(dh2_row, binding_snapshot + 8, sizeof(dh2_row));
    memcpy(dh4_row, binding_snapshot + 12, sizeof(dh4_row));
  }

  run_session_anchor(model, amplitude, amplitude_context, generation_hash,
                     chain_communicator, split_size);

  check_corrected_dispatch(
      &view, chain_communicator, split_size, 1,
      MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL,
                           "electronic real order-1");

  {
    MVMCPowerLanczosClassicView invalid_view = view;
    MVMCPowerLanczosClassicBridge *invalid_bridge = NULL;
    invalid_view.transfer_count = -1;
    status = mvmc_power_lanczos_classic_bridge_create(
        &invalid_view, world_communicator(), chain_communicator,
        &invalid_bridge);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              invalid_bridge == NULL,
          "classic bridge negative count boundary status=%s",
          mvmc_krylov_status_string(status));
    mvmc_power_lanczos_classic_bridge_destroy(invalid_bridge);

    invalid_view = view;
    invalid_view.unsupported_amplitude_features =
        MVMC_CLASSIC_KRYLOV_UNSUPPORTED_RBM;
    status = mvmc_power_lanczos_classic_bridge_create(
        &invalid_view, world_communicator(), chain_communicator,
        &invalid_bridge);
    CHECK(status == MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL &&
              invalid_bridge == NULL,
          "classic bridge unsupported feature boundary status=%s",
          mvmc_krylov_status_string(status));
    mvmc_power_lanczos_classic_bridge_destroy(invalid_bridge);

    invalid_view = view;
    invalid_view.qp_start = -1;
    status = mvmc_power_lanczos_classic_bridge_create(
        &invalid_view, world_communicator(), chain_communicator,
        &invalid_bridge);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              invalid_bridge == NULL,
          "classic bridge QP partition boundary status=%s",
          mvmc_krylov_status_string(status));
    mvmc_power_lanczos_classic_bridge_destroy(invalid_bridge);

    invalid_view = view;
    if (world_rank == 0) invalid_view.transfer_indices = NULL;
    status = mvmc_power_lanczos_classic_bridge_create(
        &invalid_view, world_communicator(), chain_communicator,
        &invalid_bridge);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              invalid_bridge == NULL,
          "classic bridge root model pointer boundary status=%s",
          mvmc_krylov_status_string(status));
    mvmc_power_lanczos_classic_bridge_destroy(invalid_bridge);
  }

  {
    MVMCPowerLanczosClassicView complex_view = view;
    MVMCKrylovScaledAmplitudeCallback complex_amplitude;
    void *complex_context;
    uint64_t complex_generation_hash;
    complex_view.arithmetic = MVMC_POWER_LANCZOS_CLASSIC_COMPLEX;
    complex_view.slater_real = NULL;
    complex_view.slater_complex = complex_slater;
    status = mvmc_power_lanczos_classic_bridge_create(
        &complex_view, world_communicator(), chain_communicator,
        &complex_bridge);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && complex_bridge != NULL,
          "complex classic bridge status=%s",
          mvmc_krylov_status_string(status));
    if (status == MVMC_KRYLOV_STATUS_OK && complex_bridge != NULL) {
      complex_amplitude =
          mvmc_power_lanczos_classic_bridge_amplitude(complex_bridge);
      complex_context =
          mvmc_power_lanczos_classic_bridge_amplitude_context(
              complex_bridge);
      complex_generation_hash =
          mvmc_power_lanczos_classic_bridge_generation_hash(complex_bridge);
      memset(&bridge_summary, 0, sizeof(bridge_summary));
      CHECK(mvmc_power_lanczos_classic_bridge_summary(
                complex_bridge, &bridge_summary) ==
                    MVMC_KRYLOV_STATUS_OK &&
                bridge_summary.valid &&
                bridge_summary.arithmetic ==
                    MVMC_POWER_LANCZOS_CLASSIC_COMPLEX &&
                bridge_summary.term_count == 5 &&
                bridge_summary.operator_count == 12 &&
                bridge_summary.amplitude_generation_hash ==
                    complex_generation_hash &&
                complex_generation_hash != 0,
            "complex classic bridge summary");
      status = mvmc_bounded_krylov_evaluate(
          bounded_workspace, &root, 1, complex_amplitude,
          complex_context, &result);
      CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
                result.evaluated_order == 3 &&
                scaled_is_finite_nonzero(&result.value[0]),
            "complex classic bridge bounded evaluation status=%s",
            mvmc_krylov_status_string(status));
      check_corrected_dispatch(
          &complex_view, chain_communicator, split_size, 2,
          MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL,
          "electronic complex order-2 provided terminal");
      check_corrected_dispatch(
          &complex_view, chain_communicator, split_size, 2,
          MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_UNIFORM_SECTOR,
          "electronic complex order-2 corrected bootstrap");
    }
  }

  {
    int inter_row[2] = {0, 1};
    int hund_row[2] = {0, 1};
    int exchange_row[2] = {1, 0};
    int *inter_indices[1] = {inter_row};
    int *hund_indices[1] = {hund_row};
    int *exchange_indices[1] = {exchange_row};
    double inter_parameters[1] = {-0.25};
    double hund_parameters[1] = {-0.5};
    double exchange_parameters[1] = {-0.5};
    MVMCPowerLanczosClassicView pure_spin_view = view;
    MVMCPowerLanczosClassicBridge *pure_spin_bridge = NULL;
    const MVMCKrylovFockModel *pure_spin_model;
    pure_spin_view.pure_spin = 1;
    pure_spin_view.transfer_count = 0;
    pure_spin_view.transfer_indices = NULL;
    pure_spin_view.transfer_parameters = NULL;
    pure_spin_view.coulomb_intra_count = 0;
    pure_spin_view.coulomb_intra_indices = NULL;
    pure_spin_view.coulomb_intra_parameters = NULL;
    pure_spin_view.coulomb_inter_count = 1;
    pure_spin_view.coulomb_inter_indices = inter_indices;
    pure_spin_view.coulomb_inter_parameters = inter_parameters;
    pure_spin_view.hund_count = 1;
    pure_spin_view.hund_indices = hund_indices;
    pure_spin_view.hund_parameters = hund_parameters;
    pure_spin_view.exchange_count = 1;
    pure_spin_view.exchange_indices = exchange_indices;
    pure_spin_view.exchange_parameters = exchange_parameters;
    status = mvmc_power_lanczos_classic_bridge_create(
        &pure_spin_view, world_communicator(), chain_communicator,
        &pure_spin_bridge);
    pure_spin_model =
        mvmc_power_lanczos_classic_bridge_model(pure_spin_bridge);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && pure_spin_bridge != NULL &&
              pure_spin_model != NULL && pure_spin_model->pure_spin == 1 &&
              pure_spin_model->term_count == 8 &&
              pure_spin_model->operator_count == 32 &&
              mvmc_power_lanczos_classic_bridge_generation_hash(
                  pure_spin_bridge) != 0,
          "pure-spin classic bridge status=%s",
          mvmc_krylov_status_string(status));
    if (status == MVMC_KRYLOV_STATUS_OK && pure_spin_bridge != NULL) {
      check_corrected_dispatch(
          &pure_spin_view, chain_communicator, split_size, 2,
          MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL,
          "pure-spin real order-2 provided terminal");
      check_corrected_dispatch(
          &pure_spin_view, chain_communicator, split_size, 2,
          MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_UNIFORM_SECTOR,
          "pure-spin real order-2 corrected bootstrap");
    }
    mvmc_power_lanczos_classic_bridge_destroy(pure_spin_bridge);
  }

cleanup:
  mvmc_power_lanczos_classic_bridge_destroy(complex_bridge);
  mvmc_bounded_krylov_workspace_destroy(bounded_workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_power_lanczos_classic_bridge_destroy(bridge);
#ifdef _mpi_use
  if (chain_communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&chain_communicator);
  }
#endif
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
