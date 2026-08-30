#include "power_lanczos_production_session.h"
#include "power_lanczos_production_session_internal.h"

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
static int TerminalBroadcastCorruptionEnabled = 0;
static size_t TerminalBroadcastCorruptionRank = 0;
static size_t TerminalBroadcastCorruptionWord = 0;
static uint64_t TerminalBroadcastCorruptionMask = 0;

MVMCKrylovStatus
mvmc_power_lanczos_production_session_testing_corrupt_terminal_after_broadcast(
    size_t target_world_rank, size_t word_index, uint64_t xor_mask) {
  if (xor_mask == 0) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  TerminalBroadcastCorruptionRank = target_world_rank;
  TerminalBroadcastCorruptionWord = word_index;
  TerminalBroadcastCorruptionMask = xor_mask;
  TerminalBroadcastCorruptionEnabled = 1;
  return MVMC_KRYLOV_STATUS_OK;
}
#endif

MVMCKrylovStatus mvmc_power_lanczos_production_dimension(
    int power_step, size_t *basis_count, size_t *upper_count) {
  size_t dimension;
  if (basis_count == NULL || upper_count == NULL ||
      (power_step != 1 && power_step != 2)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  dimension = (size_t)power_step + 1;
  *basis_count = dimension;
  *upper_count = dimension * (dimension + 1) / 2;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_fail_internal(
    MVMCPowerLanczosProductionSession *session,
    MVMCKrylovStatus status) {
  if (session != NULL) {
    session->valid = 0;
    session->status = status;
    session->state = MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FAILED;
  }
  return status;
}

static MVMCKrylovStatus create_verified_local(
    const MVMCPowerLanczosProductionSessionConfig *config,
    const uint64_t *canonical_base_words, size_t word_count,
    MVMCPowerLanczosProductionSession **session) {
  MVMCPowerLanczosProductionSession *candidate;
  MVMCPowerLanczosSnapshotSummary snapshot_summary;
  MVMCKrylovStatus status;
  size_t basis_count = 0;
  size_t upper_count = 0;
  if (session == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *session = NULL;
  if (config == NULL || canonical_base_words == NULL || word_count == 0 ||
      config->base_generation == 0 ||
      config->configuration_bit_count == 0 ||
      config->mpi_world_size == 0 || config->split_size == 0 ||
      config->mpi_world_rank >= config->mpi_world_size) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_production_dimension(
      config->power_step, &basis_count, &upper_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  candidate = (MVMCPowerLanczosProductionSession *)calloc(
      1, sizeof(*candidate));
  if (candidate == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
#ifdef _mpi_use
  candidate->world_communicator = MPI_COMM_NULL;
  candidate->chain_communicator = MPI_COMM_NULL;
  candidate->leader_communicator = MPI_COMM_NULL;
#endif
  candidate->maximum_leave_one_projective_distance = NAN;
  candidate->final_energy = NAN;
  candidate->final_energy_imaginary = NAN;
  candidate->corrected_variance = NAN;
  candidate->corrected_variance_imaginary = NAN;
  candidate->final_local_energy_abs_squared_diagnostic = NAN;
  candidate->state = MVMC_POWER_LANCZOS_PRODUCTION_SESSION_CREATED;
  candidate->config = *config;
  candidate->basis_count = basis_count;
  candidate->upper_count = upper_count;
  status = mvmc_power_lanczos_snapshot_create(
      canonical_base_words, word_count, config->configuration_bit_count,
      config->base_generation, &candidate->snapshot);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_rng_derive(
        config->resolved_base_seed, config->run_index,
        MVMC_POWER_LANCZOS_RNG_STAGE_BASE, config->mpi_world_rank,
        config->mpi_world_size, config->split_size, &candidate->base_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_rng_derive(
        config->resolved_base_seed, config->run_index,
        MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT, config->mpi_world_rank,
        config->mpi_world_size, config->split_size,
        &candidate->coefficient_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_rng_derive(
        config->resolved_base_seed, config->run_index,
        MVMC_POWER_LANCZOS_RNG_STAGE_FINAL, config->mpi_world_rank,
        config->mpi_world_size, config->split_size, &candidate->final_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (candidate->base_rng.derived_seed ==
           candidate->coefficient_rng.derived_seed ||
       candidate->base_rng.derived_seed == candidate->final_rng.derived_seed ||
       candidate->coefficient_rng.derived_seed ==
           candidate->final_rng.derived_seed)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_summary(
        candidate->snapshot, &snapshot_summary);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    mvmc_power_lanczos_production_session_destroy(candidate);
    return status;
  }
  if (snapshot_summary.allocated_bytes > SIZE_MAX - sizeof(*candidate)) {
    mvmc_power_lanczos_production_session_destroy(candidate);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  candidate->allocated_bytes =
      sizeof(*candidate) + snapshot_summary.allocated_bytes;
  candidate->valid = 1;
  candidate->status = MVMC_KRYLOV_STATUS_OK;
  candidate->state =
      MVMC_POWER_LANCZOS_PRODUCTION_SESSION_SNAPSHOTS_VERIFIED;
  *session = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

static size_t configuration_particle_count(const uint64_t *words,
                                           size_t word_count) {
  size_t total = 0;
  size_t index;
  for (index = 0; index < word_count; ++index) {
    uint64_t word = words[index];
    while (word != 0) {
      word &= word - UINT64_C(1);
      ++total;
    }
  }
  return total;
}

static MVMCKrylovStatus synchronize_creation_status(
    MVMCKrylovBoundedCommunicator communicator,
    MVMCKrylovStatus local_status) {
#ifdef _mpi_use
  int local = (int)local_status;
  int global = (int)MVMC_KRYLOV_STATUS_OK;
  if (MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global < (int)MVMC_KRYLOV_STATUS_OK ||
      global > (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (local_status != MVMC_KRYLOV_STATUS_OK &&
      global == (int)MVMC_KRYLOV_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  return (MVMCKrylovStatus)global;
#else
  (void)communicator;
  return local_status;
#endif
}

static MVMCKrylovStatus validate_root_terminal(
    const MVMCPowerLanczosProductionSessionConfig *config,
    const MVMCPowerLanczosTerminalSnapshotInput *terminal,
    unsigned char digest[MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES]) {
  MVMCKrylovStatus status;
  if (config == NULL || terminal == NULL || digest == NULL ||
      terminal->canonical_words == NULL || terminal->word_count == 0 ||
      terminal->configuration_bit_count == 0 || terminal->generation == 0 ||
      terminal->terminal_proposal_counter == 0 ||
      terminal->expected_particle_count == 0 ||
      terminal->expected_particle_count > terminal->configuration_bit_count ||
      terminal->generation != config->base_generation ||
      terminal->configuration_bit_count != config->configuration_bit_count ||
      !isfinite(terminal->log_sampling_support)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_snapshot_digest(
      terminal->canonical_words, terminal->word_count,
      terminal->configuration_bit_count, terminal->generation, digest);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      configuration_particle_count(terminal->canonical_words,
                                   terminal->word_count) !=
          terminal->expected_particle_count) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return status;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_create(
    const MVMCPowerLanczosProductionSessionConfig *config,
    MVMCKrylovBoundedCommunicator world_communicator,
    const MVMCPowerLanczosTerminalSnapshotInput *rank_0_terminal,
    MVMCPowerLanczosProductionSession **session) {
  uint64_t metadata[5] = {0, 0, 0, 0, 0};
  unsigned char expected_digest[MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES] = {0};
  unsigned char actual_digest[MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES] = {0};
  uint64_t *selected_words = NULL;
  size_t word_count = 0;
  size_t bit_count = 0;
  size_t expected_particle_count = 0;
  size_t word_bytes = 0;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  MVMCKrylovStatus root_status = MVMC_KRYLOV_STATUS_OK;
  MVMCPowerLanczosProductionSession *candidate = NULL;
#ifdef _mpi_use
  int initialized = 0;
  int actual_rank = -1;
  int actual_size = 0;
#endif
  if (session == NULL) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  } else {
    *session = NULL;
  }
  if (config == NULL || config->mpi_world_size == 0 ||
      config->mpi_world_rank >= config->mpi_world_size) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (config != NULL && config->mpi_world_rank != 0 &&
      rank_0_terminal != NULL) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#ifdef _mpi_use
  if (world_communicator == MPI_COMM_NULL ||
      MPI_Initialized(&initialized) != MPI_SUCCESS || !initialized) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (MPI_Comm_rank(world_communicator, &actual_rank) != MPI_SUCCESS ||
      MPI_Comm_size(world_communicator, &actual_size) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (actual_rank < 0 || actual_size <= 0 ||
      (config != NULL &&
       ((size_t)actual_rank != config->mpi_world_rank ||
        (size_t)actual_size != config->mpi_world_size))) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#else
  (void)world_communicator;
  if (config != NULL &&
      (config->mpi_world_rank != 0 || config->mpi_world_size != 1)) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#endif
  status = synchronize_creation_status(world_communicator, status);
  if (status != MVMC_KRYLOV_STATUS_OK || config == NULL || session == NULL) {
    return status != MVMC_KRYLOV_STATUS_OK
               ? status : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }

  if (config->mpi_world_rank == 0) {
    root_status = validate_root_terminal(
        config, rank_0_terminal, expected_digest);
    if (root_status == MVMC_KRYLOV_STATUS_OK) {
      metadata[0] = rank_0_terminal->terminal_proposal_counter;
      metadata[1] = rank_0_terminal->generation;
      metadata[2] = (uint64_t)rank_0_terminal->word_count;
      metadata[3] = (uint64_t)rank_0_terminal->configuration_bit_count;
      metadata[4] = (uint64_t)rank_0_terminal->expected_particle_count;
      if ((size_t)metadata[2] != rank_0_terminal->word_count ||
          (size_t)metadata[3] != rank_0_terminal->configuration_bit_count ||
          (size_t)metadata[4] != rank_0_terminal->expected_particle_count) {
        root_status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
    }
  }
#ifdef _mpi_use
  {
    int root_status_value = (int)root_status;
    if (MPI_Bcast(&root_status_value, 1, MPI_INT, 0,
                  world_communicator) != MPI_SUCCESS) {
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    if (root_status_value < (int)MVMC_KRYLOV_STATUS_OK ||
        root_status_value >
            (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    root_status = (MVMCKrylovStatus)root_status_value;
  }
#endif
  if (root_status != MVMC_KRYLOV_STATUS_OK) return root_status;

#ifdef _mpi_use
  if (MPI_Bcast(metadata, 5, MPI_UINT64_T, 0, world_communicator) !=
          MPI_SUCCESS ||
      MPI_Bcast(expected_digest, MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES,
                MPI_UNSIGNED_CHAR, 0, world_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  word_count = (size_t)metadata[2];
  bit_count = (size_t)metadata[3];
  expected_particle_count = (size_t)metadata[4];
  if ((uint64_t)word_count != metadata[2] ||
      (uint64_t)bit_count != metadata[3] ||
      (uint64_t)expected_particle_count != metadata[4] ||
      metadata[0] == 0 || metadata[1] == 0 || word_count == 0 ||
      bit_count == 0 || expected_particle_count == 0 ||
      metadata[1] != config->base_generation ||
      bit_count != config->configuration_bit_count ||
      word_count > SIZE_MAX / sizeof(*selected_words)) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = synchronize_creation_status(world_communicator, status);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;

  word_bytes = word_count * sizeof(*selected_words);
  selected_words = (uint64_t *)malloc(word_bytes);
  if (selected_words == NULL) {
    status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  } else if (config->mpi_world_rank == 0) {
    memcpy(selected_words, rank_0_terminal->canonical_words, word_bytes);
  }
  status = synchronize_creation_status(world_communicator, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    free(selected_words);
    return status;
  }

#ifdef _mpi_use
  {
    size_t offset = 0;
    while (offset < word_count) {
      const size_t remaining = word_count - offset;
      const int chunk = remaining > (size_t)INT_MAX
                            ? INT_MAX
                            : (int)remaining;
      if (MPI_Bcast(selected_words + offset, chunk, MPI_UINT64_T, 0,
                    world_communicator) != MPI_SUCCESS) {
        status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
        break;
      }
      offset += (size_t)chunk;
    }
  }
#endif
#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
  if (TerminalBroadcastCorruptionEnabled) {
    if (config->mpi_world_rank == TerminalBroadcastCorruptionRank) {
      if (TerminalBroadcastCorruptionWord < word_count) {
        selected_words[TerminalBroadcastCorruptionWord] ^=
            TerminalBroadcastCorruptionMask;
      } else {
        status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
    }
    TerminalBroadcastCorruptionEnabled = 0;
    TerminalBroadcastCorruptionRank = 0;
    TerminalBroadcastCorruptionWord = 0;
    TerminalBroadcastCorruptionMask = 0;
  }
#endif
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_digest(
        selected_words, word_count, bit_count, metadata[1], actual_digest);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (memcmp(actual_digest, expected_digest, sizeof(actual_digest)) != 0 ||
       configuration_particle_count(selected_words, word_count) !=
           expected_particle_count)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  status = synchronize_creation_status(world_communicator, status);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = create_verified_local(config, selected_words, word_count,
                                   &candidate);
  }
  status = synchronize_creation_status(world_communicator, status);
  free(selected_words);
  if (status != MVMC_KRYLOV_STATUS_OK || candidate == NULL) {
    mvmc_power_lanczos_production_session_destroy(candidate);
    return status != MVMC_KRYLOV_STATUS_OK
               ? status : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }

  candidate->snapshot_selection.valid = 1;
  candidate->snapshot_selection.version =
      MVMC_POWER_LANCZOS_TERMINAL_SNAPSHOT_SELECTION_VERSION;
  candidate->snapshot_selection.selection_kind =
      MVMC_POWER_LANCZOS_TERMINAL_SNAPSHOT_WORLD_RANK_0_BASE_TERMINAL_AFTER_LAST_PROPOSAL;
  candidate->snapshot_selection.owner_world_rank = 0;
  candidate->snapshot_selection.terminal_proposal_counter = metadata[0];
  candidate->snapshot_selection.generation = metadata[1];
  candidate->snapshot_selection.word_count = word_count;
  candidate->snapshot_selection.configuration_bit_count = bit_count;
  candidate->snapshot_selection.expected_particle_count =
      expected_particle_count;
  memcpy(candidate->snapshot_selection.snapshot_sha256, expected_digest,
         sizeof(expected_digest));
  *session = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_production_session_destroy(
    MVMCPowerLanczosProductionSession *session) {
  if (session == NULL) return;
  if (session->coefficient_trace == session->coefficient_local_trace) {
    session->coefficient_local_trace = NULL;
  }
  if (session->final_trace == session->final_local_trace) {
    session->final_local_trace = NULL;
  }
  mvmc_power_lanczos_primitive_trace_destroy(session->final_trace);
  mvmc_power_lanczos_primitive_trace_destroy(session->coefficient_trace);
  mvmc_power_lanczos_primitive_trace_destroy(session->final_local_trace);
  mvmc_power_lanczos_primitive_trace_destroy(
      session->coefficient_local_trace);
  mvmc_power_lanczos_observable_confirmation_destroy(
      session->observable_confirmation);
  mvmc_power_lanczos_observable_block_destroy(
      session->observable_final_blocks);
  mvmc_power_lanczos_observable_block_destroy(
      session->observable_coefficient_blocks);
  mvmc_power_lanczos_observable_evaluator_workspace_destroy(
      session->observable_evaluator);
  mvmc_power_lanczos_block_collective_destroy(
      session->block_collective_workspace);
  mvmc_power_lanczos_chain_rng_destroy(session->chain_rng_workspace);
  mvmc_bounded_krylov_workspace_destroy(session->bounded_workspace);
  mvmc_bounded_krylov_plan_destroy(session->bounded_plan);
  free(session->observable_leave_one_projective_distances);
  free(session->observable_leave_one_coefficients);
  free(session->observable_final_block_means);
  free(session->observable_leave_one_estimates);
  free(session->observable_final_estimates);
  free(session->observable_coefficient_estimates);
  free(session->observable_block_counts);
  free(session->observable_final_block_sums);
  free(session->observable_coefficient_block_sums);
  free(session->observable_final_sample);
  free(session->observable_coefficient_sample);
  free(session->observable_numeric_evidence_sample);
  free(session->primitive_trace_support);
  free(session->primitive_trace_bounds);
  free(session->primitive_trace_values);
  free(session->chain_provenance_digests);
  free(session->chain_digest_buffer);
  free(session->leave_one_projective_distances);
  free(session->leave_one_second_moments);
  free(session->leave_one_energies);
  free(session->leave_one_coefficients);
  free(session->final_evidence_block_sums);
  free(session->final_evidence_block_counts);
  free(session->coefficient_evidence_block_sums);
  free(session->coefficient_evidence_block_counts);
  free(session->global_block_sums);
  free(session->local_block_sums);
  free(session->global_block_counts);
  free(session->local_block_counts);
  free(session->proposal_words);
  free(session->final_local_accumulators);
  free(session->coefficient_hamiltonian_squared);
  free(session->coefficient_hamiltonian_adjoint);
  free(session->coefficient_hamiltonian);
  free(session->coefficient_overlap);
  free(session->coefficient_blocks);
  free(session->model_operators);
  free(session->model_terms);
#ifdef _mpi_use
  if (session->leader_communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&session->leader_communicator);
  }
  if (session->chain_communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&session->chain_communicator);
  }
  if (session->world_communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&session->world_communicator);
  }
#endif
  mvmc_power_lanczos_observable_plan_destroy(&session->observable_plan);
  mvmc_power_lanczos_snapshot_destroy(session->snapshot);
  memset(session, 0, sizeof(*session));
  free(session);
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_verify(
    MVMCPowerLanczosProductionSession *session) {
  MVMCPowerLanczosSnapshotSummary snapshot_summary;
  MVMCKrylovStatus status;
  int matches = 0;
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state !=
          MVMC_POWER_LANCZOS_PRODUCTION_SESSION_SNAPSHOTS_VERIFIED) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_snapshot_verify(session->snapshot);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_summary(session->snapshot,
                                                 &snapshot_summary);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (!session->snapshot_selection.valid ||
       session->snapshot_selection.version !=
           MVMC_POWER_LANCZOS_TERMINAL_SNAPSHOT_SELECTION_VERSION ||
       session->snapshot_selection.selection_kind !=
           MVMC_POWER_LANCZOS_TERMINAL_SNAPSHOT_WORLD_RANK_0_BASE_TERMINAL_AFTER_LAST_PROPOSAL ||
       session->snapshot_selection.owner_world_rank != 0 ||
       session->snapshot_selection.terminal_proposal_counter == 0 ||
       session->snapshot_selection.generation != snapshot_summary.generation ||
       session->snapshot_selection.word_count != snapshot_summary.word_count ||
       session->snapshot_selection.configuration_bit_count !=
           snapshot_summary.configuration_bit_count ||
       session->snapshot_selection.expected_particle_count == 0 ||
       session->snapshot_selection.expected_particle_count >
           snapshot_summary.configuration_bit_count ||
       memcmp(session->snapshot_selection.snapshot_sha256,
              snapshot_summary.base_sha256,
              MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES) != 0)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_rng_domain_matches(
        &session->base_rng, session->config.resolved_base_seed,
        session->config.run_index, MVMC_POWER_LANCZOS_RNG_STAGE_BASE,
        session->config.mpi_world_rank, session->config.mpi_world_size,
        session->config.split_size, &matches);
    if (status == MVMC_KRYLOV_STATUS_OK && !matches) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_rng_domain_matches(
        &session->coefficient_rng, session->config.resolved_base_seed,
        session->config.run_index,
        MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT,
        session->config.mpi_world_rank, session->config.mpi_world_size,
        session->config.split_size, &matches);
    if (status == MVMC_KRYLOV_STATUS_OK && !matches) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_rng_domain_matches(
        &session->final_rng, session->config.resolved_base_seed,
        session->config.run_index, MVMC_POWER_LANCZOS_RNG_STAGE_FINAL,
        session->config.mpi_world_rank, session->config.mpi_world_size,
        session->config.split_size, &matches);
    if (status == MVMC_KRYLOV_STATUS_OK && !matches) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(
        session, status);
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_summary(
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosProductionSessionSummary *summary) {
  MVMCPowerLanczosProductionSessionSummary candidate;
  MVMCKrylovStatus status;
  if (summary == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(summary, 0, sizeof(*summary));
  summary->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (session == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = session->valid;
  candidate.status = session->status;
  candidate.version = MVMC_POWER_LANCZOS_PRODUCTION_SESSION_VERSION;
  candidate.state = session->state;
  candidate.power_step = session->config.power_step;
  candidate.basis_count = session->basis_count;
  candidate.upper_count = session->upper_count;
  candidate.allocated_bytes = session->allocated_bytes;
  candidate.execution_prepared = session->execution_prepared;
  candidate.chain_rank = session->chain_rank;
  candidate.chain_size = session->chain_size;
  candidate.chain_count = session->chain_count;
  candidate.block_count = session->execution.block_count;
  candidate.execution_fingerprint = session->execution_fingerprint;
  candidate.scale_pilot_enabled = session->execution.scale_pilot_enabled;
  candidate.scale_pilot_warm_up =
      (uint64_t)session->execution.scale_pilot_warm_up;
  candidate.scale_pilot_sample_count_per_chain =
      (uint64_t)session->execution.scale_pilot_sample_count;
  candidate.scale_pilot_proposals = session->scale_pilot_proposals;
  candidate.scale_pilot_accepted_steps =
      session->scale_pilot_accepted_steps;
  candidate.scale_pilot_sample_count = session->scale_pilot_sample_count;
  candidate.eta_relative = session->execution.eta_relative;
  candidate.resolved_eta =
      session->execution.coefficient_guide_policy.eta;
  memcpy(candidate.log_basis_scale,
         session->execution.coefficient_guide_policy.log_basis_scale,
         sizeof(candidate.log_basis_scale));
  candidate.coefficient_proposals = session->coefficient_proposals;
  candidate.coefficient_accepted_steps =
      session->coefficient_accepted_steps;
  candidate.final_proposals = session->final_proposals;
  candidate.final_accepted_steps = session->final_accepted_steps;
  candidate.coefficient_sample_count = session->coefficient_sample_count;
  candidate.final_sample_count = session->final_sample_count;
  candidate.observable_enabled = session->execution.observable_enabled;
  candidate.observable_request_count = session->observable_request_count;
  memcpy(candidate.observable_census_sha256,
         session->observable_plan.raw_observable_census_sha256,
         sizeof(candidate.observable_census_sha256));
  memcpy(candidate.observable_semantic_sha256,
         session->observable_plan.semantic_observable_census_sha256,
         sizeof(candidate.observable_semantic_sha256));
  candidate.maximum_leave_one_projective_distance =
      session->maximum_leave_one_projective_distance;
  candidate.final_energy = session->final_energy;
  candidate.final_energy_imaginary = session->final_energy_imaginary;
  candidate.corrected_variance = session->corrected_variance;
  candidate.corrected_variance_imaginary =
      session->corrected_variance_imaginary;
  candidate.final_local_energy_abs_squared_diagnostic =
      session->final_local_energy_abs_squared_diagnostic;
  candidate.coefficient_provenance_hash =
      session->coefficient_provenance_hash;
  candidate.final_policy_hash = session->final_policy_hash;
  memcpy(candidate.coefficient_provenance_sha256,
         session->coefficient_provenance_sha256,
         sizeof(candidate.coefficient_provenance_sha256));
  candidate.coefficient_gevp = session->coefficient_gevp;
  candidate.snapshot_selection = session->snapshot_selection;
  candidate.base_rng = session->base_rng;
  candidate.coefficient_rng = session->coefficient_rng;
  candidate.final_rng = session->final_rng;
  status = mvmc_power_lanczos_snapshot_summary(
      session->snapshot, &candidate.snapshot);
  if (status != MVMC_KRYLOV_STATUS_OK &&
      session->status == MVMC_KRYLOV_STATUS_OK) {
    return status;
  }
  *summary = candidate;
  return candidate.status;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_observable_result(
    const MVMCPowerLanczosProductionSession *session, size_t request_index,
    MVMCPowerLanczosProductionObservableResult *result) {
  MVMCPowerLanczosProductionObservableResult candidate;
  const MVMCPowerLanczosObservableRecord *record;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state != MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED ||
      !session->execution.observable_enabled ||
      request_index >= session->observable_request_count ||
      session->observable_plan.records == NULL ||
      session->observable_coefficient_estimates == NULL ||
      session->observable_final_estimates == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  record = &session->observable_plan.records[request_index];
  candidate.valid = 1;
  candidate.request_index = request_index;
  candidate.raw_ordinal = record->raw_ordinal;
  candidate.family = record->family;
  memcpy(candidate.canonical_operator_id, record->canonical_operator_id,
         sizeof(candidate.canonical_operator_id));
  candidate.coefficient_estimate =
      session->observable_coefficient_estimates[request_index];
  candidate.final_estimate = session->observable_final_estimates[request_index];
  if (!isfinite(creal(candidate.coefficient_estimate)) ||
      !isfinite(cimag(candidate.coefficient_estimate)) ||
      !isfinite(creal(candidate.final_estimate)) ||
      !isfinite(cimag(candidate.final_estimate))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus
mvmc_power_lanczos_production_session_observable_block_result(
    const MVMCPowerLanczosProductionSession *session, size_t block_index,
    size_t request_index,
    MVMCPowerLanczosProductionObservableBlockResult *result) {
  MVMCPowerLanczosProductionObservableBlockResult candidate;
  size_t request_offset;
  size_t coefficient_offset;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state != MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED ||
      !session->execution.observable_enabled ||
      block_index >= session->execution.block_count ||
      request_index >= session->observable_request_count ||
      session->observable_leave_one_estimates == NULL ||
      session->observable_final_block_means == NULL ||
      session->observable_leave_one_coefficients == NULL ||
      session->observable_leave_one_projective_distances == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  request_offset =
      block_index * session->observable_request_count + request_index;
  coefficient_offset =
      block_index * MVMC_POWER_LANCZOS_CONFIRMATION_DIMENSION;
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.block_index = block_index;
  candidate.request_index = request_index;
  candidate.coefficient_leave_one_estimate =
      session->observable_leave_one_estimates[request_offset];
  candidate.final_block_mean =
      session->observable_final_block_means[request_offset];
  memcpy(candidate.leave_one_coefficient,
         session->observable_leave_one_coefficients + coefficient_offset,
         sizeof(candidate.leave_one_coefficient));
  candidate.leave_one_projective_distance =
      session->observable_leave_one_projective_distances[block_index];
  if (!isfinite(creal(candidate.coefficient_leave_one_estimate)) ||
      !isfinite(cimag(candidate.coefficient_leave_one_estimate)) ||
      !isfinite(creal(candidate.final_block_mean)) ||
      !isfinite(cimag(candidate.final_block_mean)) ||
      !isfinite(creal(candidate.leave_one_coefficient[0])) ||
      !isfinite(cimag(candidate.leave_one_coefficient[0])) ||
      !isfinite(creal(candidate.leave_one_coefficient[1])) ||
      !isfinite(cimag(candidate.leave_one_coefficient[1])) ||
      !isfinite(candidate.leave_one_projective_distance)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_matrix_result(
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosProductionMatrixResult *result) {
  MVMCPowerLanczosProductionMatrixResult candidate;
  size_t entry;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state <
          MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_FROZEN ||
      session->state > MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED ||
      !session->coefficient_evidence_valid || session->basis_count > 3 ||
      session->upper_count > 6 || session->coefficient_sample_count == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.basis_count = session->basis_count;
  candidate.upper_count = session->upper_count;
  candidate.sample_count = session->coefficient_sample_count;
  memcpy(candidate.overlap, session->coefficient_overlap_mean,
         session->upper_count * sizeof(*candidate.overlap));
  memcpy(candidate.hamiltonian, session->coefficient_hamiltonian_mean,
         session->upper_count * sizeof(*candidate.hamiltonian));
  memcpy(candidate.hamiltonian_adjoint,
         session->coefficient_hamiltonian_adjoint_mean,
         session->upper_count * sizeof(*candidate.hamiltonian_adjoint));
  memcpy(candidate.hamiltonian_squared,
         session->coefficient_hamiltonian_squared_mean,
         session->upper_count * sizeof(*candidate.hamiltonian_squared));
  for (entry = 0; entry < session->upper_count; ++entry) {
    if (!isfinite(creal(candidate.overlap[entry])) ||
        !isfinite(cimag(candidate.overlap[entry])) ||
        !isfinite(creal(candidate.hamiltonian[entry])) ||
        !isfinite(cimag(candidate.hamiltonian[entry])) ||
        !isfinite(creal(candidate.hamiltonian_adjoint[entry])) ||
        !isfinite(cimag(candidate.hamiltonian_adjoint[entry])) ||
        !isfinite(creal(candidate.hamiltonian_squared[entry])) ||
        !isfinite(cimag(candidate.hamiltonian_squared[entry]))) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_block_result(
    const MVMCPowerLanczosProductionSession *session, size_t block_index,
    MVMCPowerLanczosProductionCoefficientBlockResult *result) {
  MVMCPowerLanczosProductionCoefficientBlockResult candidate;
  size_t entry;
  size_t offset;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state <
          MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_FROZEN ||
      session->state > MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED ||
      !session->coefficient_evidence_valid ||
      block_index >= session->execution.block_count ||
      session->basis_count > 3 || session->upper_count > 6 ||
      session->coefficient_evidence_block_counts == NULL ||
      session->coefficient_evidence_block_sums == NULL ||
      session->leave_one_coefficients == NULL ||
      session->leave_one_energies == NULL ||
      session->leave_one_second_moments == NULL ||
      session->leave_one_projective_distances == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.block_index = block_index;
  candidate.basis_count = session->basis_count;
  candidate.upper_count = session->upper_count;
  candidate.sample_count =
      session->coefficient_evidence_block_counts[block_index];
  if (candidate.sample_count >= session->coefficient_sample_count) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  candidate.leave_one_sample_count =
      session->coefficient_sample_count - candidate.sample_count;
  offset = block_index * 4 * session->upper_count;
  memcpy(candidate.overlap_sum,
         session->coefficient_evidence_block_sums + offset,
         session->upper_count * sizeof(*candidate.overlap_sum));
  memcpy(candidate.hamiltonian_sum,
         session->coefficient_evidence_block_sums + offset +
             session->upper_count,
         session->upper_count * sizeof(*candidate.hamiltonian_sum));
  memcpy(candidate.hamiltonian_adjoint_sum,
         session->coefficient_evidence_block_sums + offset +
             2 * session->upper_count,
         session->upper_count *
             sizeof(*candidate.hamiltonian_adjoint_sum));
  memcpy(candidate.hamiltonian_squared_sum,
         session->coefficient_evidence_block_sums + offset +
             3 * session->upper_count,
         session->upper_count * sizeof(*candidate.hamiltonian_squared_sum));
  memcpy(candidate.leave_one_coefficient,
         session->leave_one_coefficients + block_index * session->basis_count,
         session->basis_count * sizeof(*candidate.leave_one_coefficient));
  candidate.leave_one_energy = session->leave_one_energies[block_index];
  candidate.leave_one_second_moment =
      session->leave_one_second_moments[block_index];
  candidate.leave_one_projective_distance =
      session->leave_one_projective_distances[block_index];
  if (candidate.sample_count == 0 ||
      candidate.leave_one_sample_count == 0 ||
      !isfinite(candidate.leave_one_energy) ||
      !isfinite(creal(candidate.leave_one_second_moment)) ||
      !isfinite(cimag(candidate.leave_one_second_moment)) ||
      !isfinite(candidate.leave_one_projective_distance)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  for (entry = 0; entry < session->upper_count; ++entry) {
    if (!isfinite(creal(candidate.overlap_sum[entry])) ||
        !isfinite(cimag(candidate.overlap_sum[entry])) ||
        !isfinite(creal(candidate.hamiltonian_sum[entry])) ||
        !isfinite(cimag(candidate.hamiltonian_sum[entry])) ||
        !isfinite(creal(candidate.hamiltonian_adjoint_sum[entry])) ||
        !isfinite(cimag(candidate.hamiltonian_adjoint_sum[entry])) ||
        !isfinite(creal(candidate.hamiltonian_squared_sum[entry])) ||
        !isfinite(cimag(candidate.hamiltonian_squared_sum[entry]))) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  for (entry = 0; entry < session->basis_count; ++entry) {
    if (!isfinite(creal(candidate.leave_one_coefficient[entry])) ||
        !isfinite(cimag(candidate.leave_one_coefficient[entry]))) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_final_block_result(
    const MVMCPowerLanczosProductionSession *session, size_t block_index,
    MVMCPowerLanczosProductionFinalBlockResult *result) {
  MVMCPowerLanczosProductionFinalBlockResult candidate;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      (session->state !=
           MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_READY &&
       session->state !=
           MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED) ||
      !session->final_evidence_valid ||
      block_index >= session->execution.block_count ||
      session->final_evidence_block_counts == NULL ||
      session->final_evidence_block_sums == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.block_index = block_index;
  candidate.sample_count = session->final_evidence_block_counts[block_index];
  candidate.energy_sum = session->final_evidence_block_sums[2 * block_index];
  candidate.local_energy_abs_squared_sum =
      session->final_evidence_block_sums[2 * block_index + 1];
  if (candidate.sample_count == 0 ||
      !isfinite(creal(candidate.energy_sum)) ||
      !isfinite(cimag(candidate.energy_sum)) ||
      !isfinite(creal(candidate.local_energy_abs_squared_sum)) ||
      !isfinite(cimag(candidate.local_energy_abs_squared_sum))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_trace_summary(
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosPrimitiveTraceSummary *summary) {
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->coefficient_trace == NULL) {
    if (summary != NULL) memset(summary, 0, sizeof(*summary));
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_primitive_trace_summary(
      session->coefficient_trace, summary);
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_final_trace_summary(
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosPrimitiveTraceSummary *summary) {
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->final_trace == NULL) {
    if (summary != NULL) memset(summary, 0, sizeof(*summary));
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_primitive_trace_summary(
      session->final_trace, summary);
}

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_trace_group(
    const MVMCPowerLanczosProductionSession *session,
    size_t group_ordinal, MVMCPowerLanczosPrimitiveTraceGroup *group) {
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->coefficient_trace == NULL) {
    if (group != NULL) memset(group, 0, sizeof(*group));
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_primitive_trace_group(
      session->coefficient_trace, group_ordinal, group);
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_final_trace_group(
    const MVMCPowerLanczosProductionSession *session,
    size_t group_ordinal, MVMCPowerLanczosPrimitiveTraceGroup *group) {
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->final_trace == NULL) {
    if (group != NULL) memset(group, 0, sizeof(*group));
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_primitive_trace_group(
      session->final_trace, group_ordinal, group);
}

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_trace_sample(
    const MVMCPowerLanczosProductionSession *session,
    size_t group_ordinal, size_t sample_index, double complex *values,
    size_t value_capacity, double *absolute_numeric_bounds,
    size_t bound_capacity, uint8_t *support_flags,
    size_t support_capacity, double *tail_magnitude) {
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->coefficient_trace == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_primitive_trace_sample(
      session->coefficient_trace, group_ordinal, sample_index, values,
      value_capacity, absolute_numeric_bounds, bound_capacity,
      support_flags, support_capacity, tail_magnitude);
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_final_trace_sample(
    const MVMCPowerLanczosProductionSession *session,
    size_t group_ordinal, size_t sample_index, double complex *values,
    size_t value_capacity, double *absolute_numeric_bounds,
    size_t bound_capacity, uint8_t *support_flags,
    size_t support_capacity, double *tail_magnitude) {
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->final_trace == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_primitive_trace_sample(
      session->final_trace, group_ordinal, sample_index, values,
      value_capacity, absolute_numeric_bounds, bound_capacity,
      support_flags, support_capacity, tail_magnitude);
}

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_trace_export_scalars(
    const MVMCPowerLanczosProductionSession *session, double *scalars,
    size_t scalar_capacity) {
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->coefficient_trace == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_primitive_trace_export_scalars(
      session->coefficient_trace, scalars, scalar_capacity);
}

MVMCKrylovStatus
mvmc_power_lanczos_production_session_final_trace_export_scalars(
    const MVMCPowerLanczosProductionSession *session, double *scalars,
    size_t scalar_capacity) {
  if (session == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->final_trace == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_primitive_trace_export_scalars(
      session->final_trace, scalars, scalar_capacity);
}

const char *mvmc_power_lanczos_production_session_state_string(
    MVMCPowerLanczosProductionSessionState state) {
  switch (state) {
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_CREATED:
      return "created";
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_SNAPSHOTS_VERIFIED:
      return "snapshots_verified";
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_RUNNING:
      return "coefficient_running";
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_READY:
      return "coefficient_ready";
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_FROZEN:
      return "coefficient_frozen";
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_RUNNING:
      return "final_running";
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_READY:
      return "final_ready";
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED:
      return "finalized";
    case MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FAILED:
      return "failed";
    default:
      return "invalid";
  }
}

#endif
