#include "power_lanczos_production_session_internal.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE) ||                    \
    !defined(MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
#error "production session stages require bounded Krylov and P5 core"
#endif

#include "krylov_final_state_chain.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { PROPOSAL_DIGEST_BYTES = 32, PROVENANCE_RECORD_CAPACITY = 512 };

static int valid_status(MVMCKrylovStatus status) {
  return status >= MVMC_KRYLOV_STATUS_OK &&
         status <= MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static MVMCKrylovStatus confirmation_status(
    MVMCPowerLanczosObservableConfirmationStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_CONFIRMATION_OK:
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    case MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT:
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    case MVMC_POWER_LANCZOS_CONFIRMATION_ALLOCATION_FAILURE:
      return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    case MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE:
      return MVMC_KRYLOV_STATUS_NONFINITE;
    case MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_STATE:
    case MVMC_POWER_LANCZOS_CONFIRMATION_BLOCK_CONTRACT:
    case MVMC_POWER_LANCZOS_CONFIRMATION_GEVP_FAILURE:
    case MVMC_POWER_LANCZOS_CONFIRMATION_INTERNAL_FAILURE:
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    default:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
}

static MVMCKrylovStatus synchronize_status(
    MVMCPowerLanczosProductionSession *session,
    MVMCKrylovStatus local_status, int use_world) {
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#ifdef _mpi_use
  MPI_Comm communicator =
      use_world ? session->world_communicator : session->chain_communicator;
  struct {
    int value;
    int index;
  } local, global;
  local.value = (int)effective;
  local.index = use_world ? (int)session->config.mpi_world_rank
                          : (int)session->chain_rank;
  if (MPI_Allreduce(&local, &global, 1, MPI_2INT, MPI_MAXLOC,
                    communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  effective = (MVMCKrylovStatus)global.value;
#else
  (void)session;
  (void)use_world;
#endif
  return effective;
}

static MVMCKrylovStatus reduce_chain_leader_counter(
    MVMCPowerLanczosProductionSession *session, uint64_t local_counter,
    uint64_t *global_counter) {
  uint64_t contribution;
  if (session == NULL || global_counter == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  contribution = session->chain_rank == 0 ? local_counter : UINT64_C(0);
#ifdef _mpi_use
  if (MPI_Allreduce(&contribution, global_counter, 1, MPI_UINT64_T, MPI_SUM,
                    session->world_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  *global_counter = contribution;
#endif
  return *global_counter >= contribution
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
}

static uint64_t double_bits(double value) {
  uint64_t bits = 0;
  memcpy(&bits, &value, sizeof(bits));
  return bits;
}

static void hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    *hash ^= value & UINT64_C(0xff);
    *hash *= UINT64_C(1099511628211);
    value >>= 8;
  }
}

static void digest_initialize(uint64_t hash[4], uint64_t domain) {
  int lane;
  for (lane = 0; lane < 4; ++lane) {
    hash[lane] = UINT64_C(1469598103934665603);
    hash_u64(&hash[lane], domain ^
                              (UINT64_C(0x9e3779b97f4a7c15) *
                               (uint64_t)(lane + 1)));
  }
}

static void digest_u64(uint64_t hash[4], uint64_t value) {
  int lane;
  for (lane = 0; lane < 4; ++lane) {
    hash_u64(&hash[lane], value ^
                              (UINT64_C(0xd6e8feb86659fd93) *
                               (uint64_t)(lane + 1)));
  }
}

static void store_u64_be(uint64_t value, unsigned char output[8]) {
  int index;
  for (index = 7; index >= 0; --index) {
    output[index] = (unsigned char)value;
    value >>= 8;
  }
}

static void digest_finish(const uint64_t hash[4],
                          unsigned char digest[PROPOSAL_DIGEST_BYTES]) {
  int lane;
  for (lane = 0; lane < 4; ++lane) {
    store_u64_be(hash[lane], digest + 8 * lane);
  }
}

static MVMCKrylovStatus audit_chain_digest(
    MVMCPowerLanczosProductionSession *session,
    const unsigned char digest[PROPOSAL_DIGEST_BYTES]) {
#ifdef _mpi_use
  size_t rank;
  if (MPI_Allgather(digest, PROPOSAL_DIGEST_BYTES, MPI_BYTE,
                    session->chain_digest_buffer,
                    PROPOSAL_DIGEST_BYTES, MPI_BYTE,
                    session->chain_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  for (rank = 1; rank < session->chain_size; ++rank) {
    if (memcmp(session->chain_digest_buffer,
               session->chain_digest_buffer +
                   rank * PROPOSAL_DIGEST_BYTES,
               PROPOSAL_DIGEST_BYTES) != 0) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
#else
  (void)session;
  (void)digest;
#endif
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus audit_proposal(
    MVMCPowerLanczosProductionSession *session,
    const uint64_t *current_words,
    const MVMCKrylovPositiveSamplerProposalDrawResult *draw,
    uint64_t proposal_ordinal) {
  uint64_t hash[4];
  unsigned char digest[PROPOSAL_DIGEST_BYTES];
  size_t word;
  if (current_words == NULL || draw == NULL || !draw->valid ||
      draw->status != MVMC_KRYLOV_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  digest_initialize(hash, UINT64_C(0x503650524f504f53));
  digest_u64(hash, proposal_ordinal);
  digest_u64(hash, (uint64_t)session->word_count);
  for (word = 0; word < session->word_count; ++word) {
    digest_u64(hash, current_words[word]);
  }
  for (word = 0; word < session->word_count; ++word) {
    digest_u64(hash, session->proposal_words[word]);
  }
  if (session->state ==
      MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_RUNNING) {
    size_t value;
    digest_u64(hash, session->coefficient_sampler.policy_hash);
    digest_u64(hash, session->coefficient_sampler.accepted_generation);
    digest_u64(hash, (uint64_t)(unsigned int)
                         session->coefficient_sampler.krylov.evaluated_order);
    digest_u64(hash,
               double_bits(session->coefficient_sampler.guide.log_guide));
    for (value = 0; value <= session->basis_count; ++value) {
      const MVMCScaledComplex *component =
          &session->coefficient_sampler.krylov.value[value];
      digest_u64(hash, (uint64_t)(unsigned int)component->state);
      digest_u64(hash, double_bits(creal(component->phase)));
      digest_u64(hash, double_bits(cimag(component->phase)));
      digest_u64(hash, double_bits(component->log_abs));
    }
  } else if (session->state ==
             MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_RUNNING) {
    digest_u64(hash, session->final_sampler.policy_hash);
    digest_u64(hash, session->final_sampler.accepted_generation);
    digest_u64(hash, session->final_sampler.configuration_hash);
    digest_u64(hash, session->final_sampler.integrity_hash);
    digest_u64(hash, session->final_sampler.final_state.policy_hash);
  } else {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  digest_u64(hash, (uint64_t)(unsigned int)draw->component);
  digest_u64(hash, (uint64_t)(unsigned int)draw->self_proposal);
  digest_u64(hash, (uint64_t)draw->neighbor_count);
  digest_u64(hash, (uint64_t)draw->selected_neighbor_index);
  digest_u64(hash, (uint64_t)draw->component_draw_count);
  digest_u64(hash, (uint64_t)draw->global_subset_draw_count);
  digest_u64(hash, (uint64_t)draw->shell_draw_count);
  digest_u64(hash, (uint64_t)draw->shell_distance);
  digest_u64(hash, double_bits(draw->uniform));
  digest_u64(hash, double_bits(draw->log_proposal_ratio));
  digest_u64(hash, draw->proposal_policy_hash);
  digest_u64(hash, draw->proposal_model_hash);
  digest_u64(hash, draw->rng_after.state);
  digest_u64(hash, draw->rng_after.stream);
  digest_u64(hash, draw->rng_after.draws);
  digest_finish(hash, digest);
  return audit_chain_digest(session, digest);
}

static double complex streaming_value(
    const MVMCKrylovStreamingComplexSum *sum) {
  return (sum->real.sum + sum->real.compensation) +
         I * (sum->imag.sum + sum->imag.compensation);
}

static MVMCKrylovStatus streaming_add(
    MVMCKrylovStreamingComplexSum *sum, double complex value) {
  double next;
  if (sum == NULL || !isfinite(creal(value)) || !isfinite(cimag(value))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  next = sum->real.sum + creal(value);
  if (!isfinite(next)) return MVMC_KRYLOV_STATUS_NONFINITE;
  if (fabs(sum->real.sum) >= fabs(creal(value))) {
    sum->real.compensation +=
        (sum->real.sum - next) + creal(value);
  } else {
    sum->real.compensation +=
        (creal(value) - next) + sum->real.sum;
  }
  sum->real.sum = next;
  next = sum->imag.sum + cimag(value);
  if (!isfinite(next)) return MVMC_KRYLOV_STATUS_NONFINITE;
  if (fabs(sum->imag.sum) >= fabs(cimag(value))) {
    sum->imag.compensation +=
        (sum->imag.sum - next) + cimag(value);
  } else {
    sum->imag.compensation +=
        (cimag(value) - next) + sum->imag.sum;
  }
  sum->imag.sum = next;
  return isfinite(sum->real.compensation) &&
                 isfinite(sum->imag.compensation)
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

static MVMCKrylovStatus coefficient_proposal_step(
    MVMCPowerLanczosProductionSession *session, uint64_t ordinal) {
  MVMCKrylovPositiveSamplerRng proposal_rng;
  MVMCPowerLanczosChainRngResult rng_result;
  MVMCKrylovPositiveSamplerProposalDrawResult draw;
  MVMCKrylovPositiveSamplerStepResult step;
  MVMCKrylovStatus status;
  memset(&proposal_rng, 0, sizeof(proposal_rng));
  memset(&draw, 0, sizeof(draw));
  status = mvmc_power_lanczos_chain_rng_derive_proposal(
      session->chain_rng_workspace, &session->coefficient_rng,
      &session->coefficient_rank_rng, ordinal, &proposal_rng, &rng_result);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_draw_mixture_rng(
        &session->model, &session->execution.proposal_policy,
        session->coefficient_words, session->word_count, &proposal_rng,
        session->proposal_words, session->word_count, &draw);
  }
  status = synchronize_status(session, status, 0);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = audit_proposal(session, session->coefficient_words, &draw,
                          ordinal);
  status = synchronize_status(session, status, 0);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = mvmc_krylov_positive_sampler_step(
      session->bounded_workspace,
      &session->execution.coefficient_guide_policy,
      session->coefficient_words, session->word_count,
      &session->coefficient_sampler, session->proposal_words,
      session->word_count, draw.log_proposal_ratio, draw.uniform,
      session->execution.amplitude, session->execution.amplitude_context,
      &step);
  return synchronize_status(session, status, 0);
}

static MVMCKrylovStatus run_scale_pilot(
    MVMCPowerLanczosProductionSession *session, uint64_t *ordinal) {
  double *log_squared = NULL;
  double local_max[MVMC_KRYLOV_MAX_ORDER + 1];
  double global_max[MVMC_KRYLOV_MAX_ORDER + 1];
  double local_sum[MVMC_KRYLOV_MAX_ORDER + 1];
  double global_sum[MVMC_KRYLOV_MAX_ORDER + 1];
  const size_t value_count = session->basis_count + 1;
  const size_t sample_count = session->execution.scale_pilot_sample_count;
  size_t storage_count;
  size_t step;
  size_t sample;
  size_t value_index;
  uint64_t policy_hash = UINT64_C(1469598103934665603);
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (!session->execution.scale_pilot_enabled) return status;
  if (ordinal == NULL || value_count > MVMC_KRYLOV_MAX_ORDER + 1 ||
      sample_count == 0 || sample_count > SIZE_MAX / value_count ||
      sample_count * value_count > SIZE_MAX / sizeof(*log_squared)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  storage_count = sample_count * value_count;
  log_squared = (double *)malloc(storage_count * sizeof(*log_squared));
  if (log_squared == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  for (value_index = 0; value_index < value_count; ++value_index) {
    local_max[value_index] = -INFINITY;
    global_max[value_index] = -INFINITY;
    local_sum[value_index] = 0.0;
    global_sum[value_index] = 0.0;
  }
  for (step = 0; step < session->execution.scale_pilot_warm_up; ++step) {
    if (*ordinal == UINT64_MAX) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      break;
    }
    status = coefficient_proposal_step(session, *ordinal);
    if (status != MVMC_KRYLOV_STATUS_OK) break;
    ++*ordinal;
  }
  for (sample = 0;
       status == MVMC_KRYLOV_STATUS_OK && sample < sample_count; ++sample) {
    if (*ordinal == UINT64_MAX) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      break;
    }
    status = coefficient_proposal_step(session, *ordinal);
    if (status != MVMC_KRYLOV_STATUS_OK) break;
    ++*ordinal;
    if (!session->coefficient_sampler.valid ||
        session->coefficient_sampler.status != MVMC_KRYLOV_STATUS_OK ||
        session->coefficient_sampler.krylov.evaluated_order <
            (int)session->basis_count) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      break;
    }
    for (value_index = 0; value_index < value_count; ++value_index) {
      const MVMCScaledComplex *value =
          &session->coefficient_sampler.krylov.value[value_index];
      double squared = -INFINITY;
      if (!mvmc_scaled_complex_is_valid(value)) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
        break;
      }
      if (value->state == MVMC_SCALED_COMPLEX_NONFINITE) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
        break;
      }
      if (value->state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
        squared = 2.0 * value->log_abs;
        if (!isfinite(squared)) {
          status = MVMC_KRYLOV_STATUS_NONFINITE;
          break;
        }
      }
      log_squared[sample * value_count + value_index] = squared;
      if (session->chain_rank == 0 && squared > local_max[value_index]) {
        local_max[value_index] = squared;
      }
    }
  }
  status = synchronize_status(session, status, 1);
#ifdef _mpi_use
  if (status == MVMC_KRYLOV_STATUS_OK &&
      MPI_Allreduce(local_max, global_max, (int)value_count, MPI_DOUBLE,
                    MPI_MAX, session->world_communicator) != MPI_SUCCESS) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memcpy(global_max, local_max, value_count * sizeof(*global_max));
  }
#endif
  if (status == MVMC_KRYLOV_STATUS_OK) {
    for (value_index = 0; value_index < value_count; ++value_index) {
      if (isnan(global_max[value_index]) ||
          global_max[value_index] == INFINITY) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
        break;
      }
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK && session->chain_rank == 0) {
    for (sample = 0; sample < sample_count; ++sample) {
      for (value_index = 0; value_index < value_count; ++value_index) {
        const double squared =
            log_squared[sample * value_count + value_index];
        if (isfinite(squared)) {
          local_sum[value_index] +=
              exp(squared - global_max[value_index]);
        }
      }
    }
  }
#ifdef _mpi_use
  if (status == MVMC_KRYLOV_STATUS_OK &&
      MPI_Allreduce(local_sum, global_sum, (int)value_count, MPI_DOUBLE,
                    MPI_SUM, session->world_communicator) != MPI_SUCCESS) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memcpy(global_sum, local_sum, value_count * sizeof(*global_sum));
  }
#endif
  if (status == MVMC_KRYLOV_STATUS_OK) {
    const double global_count =
        (double)sample_count * (double)session->chain_count;
    const double eta = session->execution.eta_relative * (double)value_count;
    if (!isfinite(global_count) || global_count <= 0.0 || !isfinite(eta) ||
        eta <= 0.0 || sample_count > UINT64_MAX / session->chain_count) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    for (value_index = 0;
         status == MVMC_KRYLOV_STATUS_OK && value_index < value_count;
         ++value_index) {
      double log_mean;
      double scale;
      if (global_max[value_index] == -INFINITY &&
          global_sum[value_index] == 0.0) {
        /* The strictly positive guide floor keeps this component sampleable.
         * A pilot with no nonzero observation supplies no scale estimate, so
         * retain the neutral log scale and let coefficient/GEVP/statistical
         * gates decide from the scored samples. */
        scale = 0.0;
      } else if (!isfinite(global_max[value_index]) ||
                 !isfinite(global_sum[value_index]) ||
                 global_sum[value_index] <= 0.0) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
        break;
      } else {
        log_mean = global_max[value_index] + log(global_sum[value_index]) -
                   log(global_count);
        scale = -0.5 * log_mean;
        if (!isfinite(log_mean) || !isfinite(scale)) {
          status = MVMC_KRYLOV_STATUS_NONFINITE;
          break;
        }
      }
      session->execution.coefficient_guide_policy.lambda[value_index] = 1.0;
      session->execution.coefficient_guide_policy
          .log_basis_scale[value_index] = scale;
      session->execution.matrix_policy.guide_lambda[value_index] = 1.0;
      session->execution.matrix_policy.target_weight[value_index] = 1.0;
      session->execution.matrix_policy.log_basis_scale[value_index] = scale;
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      session->execution.coefficient_guide_policy.eta = eta;
      session->execution.matrix_policy.eta = eta;
      hash_u64(&policy_hash, UINT64_C(0x50365343414c4531));
      hash_u64(&policy_hash, double_bits(session->execution.eta_relative));
      hash_u64(&policy_hash, (uint64_t)value_count);
      for (value_index = 0; value_index < value_count; ++value_index) {
        hash_u64(&policy_hash,
                 double_bits(session->execution.coefficient_guide_policy
                                 .log_basis_scale[value_index]));
      }
      session->execution.coefficient_guide_policy.policy_hash =
          policy_hash == 0 ? UINT64_C(1) : policy_hash;
      session->scale_pilot_proposals = *ordinal;
      session->scale_pilot_sample_count =
          (uint64_t)sample_count * (uint64_t)session->chain_count;
      status = reduce_chain_leader_counter(
          session, session->coefficient_sampler.accepted_generation,
          &session->scale_pilot_accepted_steps);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = mvmc_krylov_positive_sampler_initialize(
          session->bounded_workspace,
          &session->execution.coefficient_guide_policy,
          session->coefficient_words, session->word_count,
          session->execution.amplitude,
          session->execution.amplitude_context,
          &session->coefficient_sampler);
      }
    }
  }
  free(log_squared);
  return synchronize_status(session, status, 1);
}

static MVMCKrylovStatus extract_coefficient_blocks(
    MVMCPowerLanczosProductionSession *session,
    MVMCKrylovStatus local_status) {
  const size_t matrix_entries_per_block = 4 * session->upper_count;
  const size_t entries_per_block = session->coefficient_entries_per_block;
  const size_t observable_entries_per_block =
      session->observable_request_count *
      MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES;
  size_t block;
  size_t entry;
  memset(session->local_block_counts, 0,
         session->execution.block_count *
             sizeof(*session->local_block_counts));
  memset(session->local_block_sums, 0,
         session->execution.block_count * session->block_entry_count *
             sizeof(*session->local_block_sums));
  if (local_status != MVMC_KRYLOV_STATUS_OK) return local_status;
  if (!session->coefficient_accumulator.valid ||
      session->coefficient_accumulator.status != MVMC_KRYLOV_STATUS_OK ||
      session->coefficient_accumulator.sample_count !=
          (uint64_t)session->execution.coefficient_sample_count ||
      entries_per_block > session->block_entry_count) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (session->execution.observable_enabled) {
    MVMCPowerLanczosObservableBlockSummary summary;
    local_status = mvmc_power_lanczos_observable_block_summary(
        session->observable_coefficient_blocks, &summary);
    if (local_status == MVMC_KRYLOV_STATUS_OK &&
        (!summary.valid || summary.status != MVMC_KRYLOV_STATUS_OK ||
         summary.stage != MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT ||
         summary.request_count != session->observable_request_count ||
         summary.entries_per_request !=
             MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES ||
         summary.completed_block_count != session->execution.block_count ||
         summary.completed_sample_count !=
             (uint64_t)session->execution.coefficient_sample_count ||
         summary.current_block_sample_count != 0 ||
         summary.discarded_partial_sample_count != 0)) {
      local_status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (local_status == MVMC_KRYLOV_STATUS_OK) {
      local_status = mvmc_power_lanczos_observable_block_export(
          session->observable_coefficient_blocks,
          session->observable_coefficient_block_sums,
          session->execution.block_count * observable_entries_per_block,
          session->observable_block_counts,
          session->execution.block_count);
    }
    if (local_status != MVMC_KRYLOV_STATUS_OK) return local_status;
  }
  for (block = 0; block < session->execution.block_count; ++block) {
    const size_t storage_offset = block * session->upper_count;
    const size_t output_offset = block * entries_per_block;
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator =
        &session->coefficient_blocks[block];
    session->local_block_counts[block] = accumulator->sample_count;
    if (!accumulator->valid ||
        accumulator->status != MVMC_KRYLOV_STATUS_OK ||
        accumulator->sample_count !=
            (uint64_t)(session->execution.coefficient_sample_count /
                       session->execution.block_count)) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (session->execution.observable_enabled &&
        session->observable_block_counts[block] !=
            accumulator->sample_count) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    for (entry = 0; entry < session->upper_count; ++entry) {
      session->local_block_sums[output_offset + entry] =
          streaming_value(session->coefficient_overlap + storage_offset +
                          entry);
      session->local_block_sums[output_offset + session->upper_count +
                                entry] =
          streaming_value(session->coefficient_hamiltonian + storage_offset +
                          entry);
      session->local_block_sums[output_offset + 2 * session->upper_count +
                                entry] =
          streaming_value(session->coefficient_hamiltonian_adjoint +
                          storage_offset + entry);
      session->local_block_sums[output_offset + 3 * session->upper_count +
                                entry] =
          streaming_value(session->coefficient_hamiltonian_squared +
                          storage_offset + entry);
    }
    if (session->execution.observable_enabled) {
      memcpy(session->local_block_sums + output_offset +
                 matrix_entries_per_block,
             session->observable_coefficient_block_sums +
                 block * observable_entries_per_block,
             observable_entries_per_block *
                 sizeof(*session->local_block_sums));
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus append_numeric_evidence(
    MVMCPowerLanczosProductionSession *session, size_t *primitive,
    const MVMCPowerLanczosNumericEvidence *evidence) {
  if (session == NULL || primitive == NULL || evidence == NULL ||
      !evidence->valid || evidence->status != MVMC_KRYLOV_STATUS_OK ||
      *primitive >= session->maximum_primitive_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  session->primitive_trace_values[*primitive] = evidence->value;
  session->primitive_trace_bounds[*primitive] =
      evidence->absolute_numeric_bound;
  session->primitive_trace_support[*primitive] = evidence->support_flags;
  ++(*primitive);
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus append_coefficient_primitive_trace(
    MVMCPowerLanczosProductionSession *session, size_t sample,
    const MVMCKrylovMatrixMeasurementSampleEvidence *matrix_evidence) {
  size_t primitive = 0;
  size_t entry;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (session->chain_rank != 0) return MVMC_KRYLOV_STATUS_OK;
  if (matrix_evidence == NULL || !matrix_evidence->valid ||
      matrix_evidence->status != MVMC_KRYLOV_STATUS_OK ||
      matrix_evidence->upper_count != session->upper_count ||
      !matrix_evidence->has_hamiltonian_adjoint ||
      session->coefficient_local_trace == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = append_numeric_evidence(
      session, &primitive, &matrix_evidence->denominator);
  for (entry = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  entry < session->upper_count; ++entry) {
    status = append_numeric_evidence(
        session, &primitive, matrix_evidence->overlap + entry);
  }
  for (entry = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  entry < session->upper_count; ++entry) {
    status = append_numeric_evidence(
        session, &primitive, matrix_evidence->hamiltonian + entry);
  }
  for (entry = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  entry < session->upper_count; ++entry) {
    status = append_numeric_evidence(
        session, &primitive,
        matrix_evidence->hamiltonian_adjoint + entry);
  }
  for (entry = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  entry < session->upper_count; ++entry) {
    status = append_numeric_evidence(
        session, &primitive,
        matrix_evidence->hamiltonian_squared + entry);
  }
  for (entry = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  entry < session->observable_request_count *
                              MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES;
       ++entry) {
    status = append_numeric_evidence(
        session, &primitive,
        session->observable_numeric_evidence_sample + entry);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      primitive != session->coefficient_primitive_count) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return mvmc_power_lanczos_primitive_trace_append(
      session->coefficient_local_trace, 0, sample,
      session->primitive_trace_values, session->primitive_trace_bounds,
      session->primitive_trace_support,
      session->coefficient_primitive_count, 0.0);
}

static MVMCKrylovStatus append_final_primitive_trace(
    MVMCPowerLanczosProductionSession *session, size_t sample,
    const MVMCPowerLanczosNumericEvidence *energy_evidence) {
  size_t primitive = 0;
  size_t request;
  double tail_magnitude;
  MVMCKrylovStatus status;
  if (session->chain_rank != 0) return MVMC_KRYLOV_STATUS_OK;
  if (energy_evidence == NULL || !energy_evidence->valid ||
      energy_evidence->status != MVMC_KRYLOV_STATUS_OK ||
      session->final_local_trace == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  session->primitive_trace_values[primitive] = 1.0;
  session->primitive_trace_bounds[primitive] = 0.0;
  session->primitive_trace_support[primitive] =
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO;
  ++primitive;
  status = append_numeric_evidence(session, &primitive, energy_evidence);
  for (request = 0; status == MVMC_KRYLOV_STATUS_OK &&
                    request < session->observable_request_count; ++request) {
    status = append_numeric_evidence(
        session, &primitive,
        session->observable_numeric_evidence_sample + request);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      primitive != session->final_primitive_count) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  tail_magnitude = hypot(creal(energy_evidence->value),
                         cimag(energy_evidence->value));
  if (status != MVMC_KRYLOV_STATUS_OK || !isfinite(tail_magnitude)) {
    return status == MVMC_KRYLOV_STATUS_OK
               ? MVMC_KRYLOV_STATUS_NONFINITE
               : status;
  }
  return mvmc_power_lanczos_primitive_trace_append(
      session->final_local_trace, 0, sample,
      session->primitive_trace_values, session->primitive_trace_bounds,
      session->primitive_trace_support, session->final_primitive_count,
      tail_magnitude);
}

static MVMCKrylovStatus audit_local_blocks(
    MVMCPowerLanczosProductionSession *session, size_t entries_per_block,
    MVMCKrylovStatus local_status) {
  uint64_t hash[4];
  unsigned char digest[PROPOSAL_DIGEST_BYTES];
  size_t index;
  local_status = synchronize_status(session, local_status, 0);
  if (local_status != MVMC_KRYLOV_STATUS_OK) return local_status;
  digest_initialize(hash, UINT64_C(0x5036424c4f434b53));
  digest_u64(hash, (uint64_t)session->execution.block_count);
  digest_u64(hash, (uint64_t)entries_per_block);
  for (index = 0; index < session->execution.block_count; ++index) {
    digest_u64(hash, session->local_block_counts[index]);
  }
  for (index = 0;
       index < session->execution.block_count * entries_per_block; ++index) {
    digest_u64(hash, double_bits(creal(session->local_block_sums[index])));
    digest_u64(hash, double_bits(cimag(session->local_block_sums[index])));
  }
  digest_finish(hash, digest);
  return audit_chain_digest(session, digest);
}

static MVMCKrylovStatus reduce_blocks(
    MVMCPowerLanczosProductionSession *session, size_t entries_per_block,
    MVMCKrylovStatus local_status) {
  MVMCPowerLanczosBlockCollectiveResult result;
  MVMCKrylovStatus status = local_status;
  size_t count_bytes = session->execution.block_count *
                       sizeof(*session->global_block_counts);
  size_t sum_count = session->execution.block_count * entries_per_block;
  size_t sum_bytes = sum_count * sizeof(*session->global_block_sums);
#ifdef _mpi_use
  if (session->chain_rank == 0) {
    status = mvmc_power_lanczos_block_collective_reduce(
        session->block_collective_workspace, local_status,
        session->execution.block_count, entries_per_block,
        session->local_block_counts, session->local_block_sums,
        session->global_block_counts, session->execution.block_count,
        session->global_block_sums, sum_count, &result);
  }
  {
    int encoded = (int)status;
    if (MPI_Bcast(&encoded, 1, MPI_INT, 0,
                  session->chain_communicator) != MPI_SUCCESS) {
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    status = valid_status((MVMCKrylovStatus)encoded)
                 ? (MVMCKrylovStatus)encoded
                 : MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (count_bytes > (size_t)INT_MAX || sum_bytes > (size_t)INT_MAX ||
       MPI_Bcast(session->global_block_counts, (int)count_bytes, MPI_BYTE,
                 0, session->chain_communicator) != MPI_SUCCESS ||
       MPI_Bcast(session->global_block_sums, (int)sum_bytes, MPI_BYTE, 0,
                 session->chain_communicator) != MPI_SUCCESS)) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  status = mvmc_power_lanczos_block_collective_reduce(
      session->block_collective_workspace, local_status,
      session->execution.block_count, entries_per_block,
      session->local_block_counts, session->local_block_sums,
      session->global_block_counts, session->execution.block_count,
      session->global_block_sums, sum_count, &result);
  (void)count_bytes;
  (void)sum_bytes;
#endif
  return synchronize_status(session, status, 0);
}

static MVMCKrylovStatus persist_coefficient_evidence(
    MVMCPowerLanczosProductionSession *session) {
  const size_t evidence_entries_per_block = 4 * session->upper_count;
  size_t block;
  size_t entry;
  if (session->coefficient_evidence_block_counts == NULL ||
      session->coefficient_evidence_block_sums == NULL ||
      session->global_block_counts == NULL ||
      session->global_block_sums == NULL ||
      session->coefficient_entries_per_block < evidence_entries_per_block) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  session->coefficient_evidence_valid = 0;
  for (block = 0; block < session->execution.block_count; ++block) {
    const size_t source = block * session->coefficient_entries_per_block;
    const size_t destination = block * evidence_entries_per_block;
    if (session->global_block_counts[block] == 0) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    session->coefficient_evidence_block_counts[block] =
        session->global_block_counts[block];
    for (entry = 0; entry < evidence_entries_per_block; ++entry) {
      const double complex value = session->global_block_sums[source + entry];
      if (!isfinite(creal(value)) || !isfinite(cimag(value))) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
      session->coefficient_evidence_block_sums[destination + entry] = value;
    }
  }
  session->coefficient_evidence_valid = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus persist_final_evidence(
    MVMCPowerLanczosProductionSession *session) {
  size_t block;
  size_t entry;
  if (session->final_evidence_block_counts == NULL ||
      session->final_evidence_block_sums == NULL ||
      session->global_block_counts == NULL ||
      session->global_block_sums == NULL ||
      session->final_entries_per_block < 2) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  session->final_evidence_valid = 0;
  for (block = 0; block < session->execution.block_count; ++block) {
    const size_t source = block * session->final_entries_per_block;
    if (session->global_block_counts[block] == 0) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    session->final_evidence_block_counts[block] =
        session->global_block_counts[block];
    for (entry = 0; entry < 2; ++entry) {
      const double complex value = session->global_block_sums[source + entry];
      if (!isfinite(creal(value)) || !isfinite(cimag(value))) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
      session->final_evidence_block_sums[2 * block + entry] = value;
    }
  }
  session->final_evidence_valid = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static int trace_packet_bytes(size_t primitive_count, size_t *bytes) {
  size_t value_bytes;
  size_t bound_bytes;
  size_t total;
  if (bytes == NULL || primitive_count == 0 ||
      primitive_count > SIZE_MAX / sizeof(double complex) ||
      primitive_count > SIZE_MAX / sizeof(double)) {
    return 0;
  }
  value_bytes = primitive_count * sizeof(double complex);
  bound_bytes = primitive_count * sizeof(double);
  if (value_bytes > SIZE_MAX - bound_bytes) return 0;
  total = value_bytes + bound_bytes;
  if (total > SIZE_MAX - primitive_count) return 0;
  total += primitive_count;
  if (total > SIZE_MAX - sizeof(double)) return 0;
  *bytes = total + sizeof(double);
  return 1;
}

static void pack_trace_sample(
    unsigned char *packet, size_t primitive_count,
    const double complex *values, const double *bounds,
    const uint8_t *support, double tail) {
  const size_t value_bytes = primitive_count * sizeof(values[0]);
  const size_t bound_bytes = primitive_count * sizeof(bounds[0]);
  memcpy(packet, values, value_bytes);
  memcpy(packet + value_bytes, bounds, bound_bytes);
  memcpy(packet + value_bytes + bound_bytes, support, primitive_count);
  memcpy(packet + value_bytes + bound_bytes + primitive_count,
         &tail, sizeof(tail));
}

static void unpack_trace_sample(
    const unsigned char *packet, size_t primitive_count,
    double complex *values, double *bounds, uint8_t *support,
    double *tail) {
  const size_t value_bytes = primitive_count * sizeof(values[0]);
  const size_t bound_bytes = primitive_count * sizeof(bounds[0]);
  memcpy(values, packet, value_bytes);
  memcpy(bounds, packet + value_bytes, bound_bytes);
  memcpy(support, packet + value_bytes + bound_bytes, primitive_count);
  memcpy(tail, packet + value_bytes + bound_bytes + primitive_count,
         sizeof(*tail));
}

static MVMCKrylovStatus gather_primitive_trace(
    MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosPrimitiveTraceStage stage,
    MVMCPowerLanczosPrimitiveTrace **local_trace_slot,
    MVMCPowerLanczosPrimitiveTrace **global_trace_slot,
    size_t primitive_count, size_t sample_count, int message_tag) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  MVMCPowerLanczosPrimitiveTrace *global_trace = NULL;
  size_t local_bytes = 0;
  if (session == NULL || local_trace_slot == NULL ||
      global_trace_slot == NULL || *global_trace_slot != NULL ||
      primitive_count == 0 || sample_count == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (session->chain_rank == 0) {
    if (*local_trace_slot == NULL) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    } else {
      local_bytes = mvmc_power_lanczos_primitive_trace_allocated_bytes(
          *local_trace_slot);
      status = mvmc_power_lanczos_primitive_trace_freeze(
          *local_trace_slot);
    }
  } else if (*local_trace_slot != NULL) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  status = synchronize_status(session, status, 1);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;

#ifndef _mpi_use
  (void)stage;
  (void)primitive_count;
  (void)sample_count;
  (void)message_tag;
  *global_trace_slot = *local_trace_slot;
  *local_trace_slot = NULL;
  return MVMC_KRYLOV_STATUS_OK;
#else
  {
    unsigned char *packet = NULL;
    size_t packet_bytes = 0;
    int leader_status;
    int effective_leader_status;
    if (!trace_packet_bytes(primitive_count, &packet_bytes) ||
        packet_bytes > (size_t)INT_MAX) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    if (status == MVMC_KRYLOV_STATUS_OK && session->chain_rank == 0) {
      packet = (unsigned char *)malloc(packet_bytes);
      if (packet == NULL) status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        session->config.mpi_world_rank == 0) {
      size_t group;
      status = mvmc_power_lanczos_primitive_trace_create(
          stage, primitive_count, session->chain_count, sample_count,
          session->execution.bounded_limits.max_workspace_bytes,
          &global_trace);
      for (group = 0; status == MVMC_KRYLOV_STATUS_OK &&
                      group < session->chain_count; ++group) {
        const size_t leader_world_rank =
            group * session->config.split_size;
        size_t expected_chain_size =
            session->config.mpi_world_size - leader_world_rank;
        if (expected_chain_size > session->config.split_size) {
          expected_chain_size = session->config.split_size;
        }
        status = mvmc_power_lanczos_primitive_trace_register_group(
            global_trace, group, leader_world_rank,
            expected_chain_size);
      }
    }
    status = synchronize_status(session, status, 1);
    if (status == MVMC_KRYLOV_STATUS_OK && session->chain_rank == 0 &&
        packet == NULL) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    status = synchronize_status(session, status, 1);
    if (status == MVMC_KRYLOV_STATUS_OK && session->chain_rank == 0 &&
        packet != NULL) {
      if (session->leader_rank == 0) {
        size_t group;
        for (group = 0; group < session->chain_count; ++group) {
          size_t sample;
          for (sample = 0; sample < sample_count; ++sample) {
            double tail = 0.0;
            MVMCKrylovStatus sample_status = MVMC_KRYLOV_STATUS_OK;
            if (group == 0) {
              sample_status = mvmc_power_lanczos_primitive_trace_sample(
                  *local_trace_slot, 0, sample,
                  session->primitive_trace_values, primitive_count,
                  session->primitive_trace_bounds, primitive_count,
                  session->primitive_trace_support, primitive_count,
                  &tail);
            } else if (MPI_Recv(
                           packet, (int)packet_bytes, MPI_BYTE,
                           (int)group, message_tag,
                           session->leader_communicator,
                           MPI_STATUS_IGNORE) != MPI_SUCCESS) {
              sample_status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
            } else {
              unpack_trace_sample(
                  packet, primitive_count,
                  session->primitive_trace_values,
                  session->primitive_trace_bounds,
                  session->primitive_trace_support, &tail);
            }
            if (status == MVMC_KRYLOV_STATUS_OK &&
                sample_status == MVMC_KRYLOV_STATUS_OK) {
              status = mvmc_power_lanczos_primitive_trace_append(
                  global_trace, group, sample,
                  session->primitive_trace_values,
                  session->primitive_trace_bounds,
                  session->primitive_trace_support, primitive_count,
                  tail);
            } else if (status == MVMC_KRYLOV_STATUS_OK) {
              status = sample_status;
            }
          }
        }
      } else {
        size_t sample;
        for (sample = 0; sample < sample_count; ++sample) {
          double tail = 0.0;
          MVMCKrylovStatus sample_status =
              mvmc_power_lanczos_primitive_trace_sample(
                  *local_trace_slot, 0, sample,
                  session->primitive_trace_values, primitive_count,
                  session->primitive_trace_bounds, primitive_count,
                  session->primitive_trace_support, primitive_count,
                  &tail);
          if (sample_status == MVMC_KRYLOV_STATUS_OK) {
            pack_trace_sample(
                packet, primitive_count,
                session->primitive_trace_values,
                session->primitive_trace_bounds,
                session->primitive_trace_support, tail);
          } else {
            memset(packet, 0, packet_bytes);
            if (status == MVMC_KRYLOV_STATUS_OK) status = sample_status;
          }
          if (MPI_Send(packet, (int)packet_bytes, MPI_BYTE, 0,
                       message_tag,
                       session->leader_communicator) != MPI_SUCCESS &&
              status == MVMC_KRYLOV_STATUS_OK) {
            status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
          }
        }
      }
      leader_status = (int)status;
      if (MPI_Allreduce(&leader_status, &effective_leader_status, 1,
                        MPI_INT, MPI_MAX,
                        session->leader_communicator) != MPI_SUCCESS) {
        status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
      } else {
        status = (MVMCKrylovStatus)effective_leader_status;
      }
      if (session->leader_rank == 0 &&
          status == MVMC_KRYLOV_STATUS_OK) {
        status = mvmc_power_lanczos_primitive_trace_freeze(global_trace);
      }
    }
    free(packet);
  }
  status = synchronize_status(session, status, 1);
  if (session->chain_rank == 0) {
    mvmc_power_lanczos_primitive_trace_destroy(*local_trace_slot);
    *local_trace_slot = NULL;
    if (session->allocated_bytes >= local_bytes) {
      session->allocated_bytes -= local_bytes;
    } else if (status == MVMC_KRYLOV_STATUS_OK) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  if (session->config.mpi_world_rank == 0) {
    if (status == MVMC_KRYLOV_STATUS_OK) {
      const size_t global_bytes =
          mvmc_power_lanczos_primitive_trace_allocated_bytes(global_trace);
      if (session->allocated_bytes > SIZE_MAX - global_bytes) {
        status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      } else {
        session->allocated_bytes += global_bytes;
        *global_trace_slot = global_trace;
        global_trace = NULL;
      }
    }
    mvmc_power_lanczos_primitive_trace_destroy(global_trace);
  }
  return synchronize_status(session, status, 1);
#endif
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_run_coefficient(
    MVMCPowerLanczosProductionSession *session) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  uint64_t ordinal = 0;
  size_t step;
  size_t sample;
  if (session == NULL || !session->valid || !session->execution_prepared ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state !=
          MVMC_POWER_LANCZOS_PRODUCTION_SESSION_SNAPSHOTS_VERIFIED) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  session->state =
      MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_RUNNING;
  status = mvmc_power_lanczos_snapshot_coefficient_begin(
      session->snapshot, &session->coefficient_words, &session->word_count);
  status = synchronize_status(session, status, 0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        session->bounded_workspace, session->execution.amplitude,
        session->execution.amplitude_context,
        session->execution.amplitude_generation_hash);
  }
  status = synchronize_status(session, status, 0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_rng_seed(
        session->coefficient_rng.derived_seed,
        session->coefficient_rng.stage_tag,
        &session->coefficient_rank_rng);
  }
  status = synchronize_status(session, status, 0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_initialize(
        session->bounded_workspace,
        &session->execution.coefficient_guide_policy,
        session->coefficient_words, session->word_count,
        session->execution.amplitude, session->execution.amplitude_context,
        &session->coefficient_sampler);
  }
  status = synchronize_status(session, status, 0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = run_scale_pilot(session, &ordinal);
  }
  status = synchronize_status(session, status, 1);
  for (step = 0; step < session->execution.coefficient_warm_up; ++step) {
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = coefficient_proposal_step(session, ordinal);
      if (status == MVMC_KRYLOV_STATUS_OK) ++ordinal;
    }
  }
  for (sample = 0; sample < session->execution.coefficient_sample_count;
       ++sample) {
    MVMCKrylovMatrixMeasurementSampleEvidence matrix_evidence;
    memset(&matrix_evidence, 0, sizeof(matrix_evidence));
    for (step = 0; step < session->execution.coefficient_interval; ++step) {
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = coefficient_proposal_step(session, ordinal);
        if (status == MVMC_KRYLOV_STATUS_OK) ++ordinal;
      }
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
      status =
          mvmc_krylov_matrix_measurement_block_accumulator_add_sample_with_evidence(
              &session->coefficient_accumulator,
              &session->execution.matrix_policy,
              session->coefficient_sampler.krylov.value,
              (size_t)session->coefficient_sampler.krylov.evaluated_order + 1,
              &matrix_evidence, &diagnostics);
      status = synchronize_status(session, status, 0);
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        session->execution.observable_enabled) {
      MVMCPowerLanczosObservableSampleDiagnostics diagnostics;
      status =
          mvmc_power_lanczos_observable_coefficient_sample_with_evidence(
              session->observable_evaluator, session->bounded_workspace,
              &session->observable_layout, &session->observable_plan,
              session->coefficient_words, session->word_count,
              session->coefficient_sampler.guide.log_guide,
              session->execution.matrix_policy.log_basis_scale,
              &matrix_evidence.guide,
              session->observable_coefficient_sample,
              session->observable_request_count *
                  MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES,
              session->observable_numeric_evidence_sample,
              session->observable_request_count *
                  MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES,
              &diagnostics);
      if (status == MVMC_KRYLOV_STATUS_OK &&
          (!diagnostics.valid ||
           diagnostics.status != MVMC_KRYLOV_STATUS_OK ||
           diagnostics.request_count != session->observable_request_count)) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = mvmc_power_lanczos_observable_block_add_sample(
            session->observable_coefficient_blocks,
            session->observable_coefficient_sample,
            session->observable_request_count *
                MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES);
      }
      status = synchronize_status(session, status, 0);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = append_coefficient_primitive_trace(
          session, sample, &matrix_evidence);
    }
    status = synchronize_status(session, status, 0);
  }
  session->coefficient_proposals = ordinal;
  if (mvmc_bounded_krylov_session_is_active(session->bounded_workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(session->bounded_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  status = synchronize_status(session, status, 0);
  status = synchronize_status(session, status, 1);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = gather_primitive_trace(
        session, MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_COEFFICIENT,
        &session->coefficient_local_trace, &session->coefficient_trace,
        session->coefficient_primitive_count,
        session->execution.coefficient_sample_count, 601);
  }
  status = extract_coefficient_blocks(session, status);
  status = audit_local_blocks(session, session->coefficient_entries_per_block,
                              status);
  status = reduce_blocks(session, session->coefficient_entries_per_block,
                         status);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = persist_coefficient_evidence(session);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_coefficient_complete(
        session->snapshot);
  }
  status = synchronize_status(session, status, 1);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  session->coefficient_sample_count =
      (uint64_t)session->execution.coefficient_sample_count *
      (uint64_t)session->chain_count;
  status = reduce_chain_leader_counter(
      session, session->coefficient_sampler.accepted_generation,
      &session->coefficient_accepted_steps);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  session->state =
      MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_READY;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus gather_coefficient_chain_digests(
    MVMCPowerLanczosProductionSession *session,
    const unsigned char local_digest[PROPOSAL_DIGEST_BYTES],
    MVMCKrylovStatus local_status) {
  MVMCKrylovStatus status =
      synchronize_status(session, local_status, 0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = audit_chain_digest(session, local_digest);
    status = synchronize_status(session, status, 0);
  }
#ifdef _mpi_use
  if (session->chain_rank == 0) {
    struct {
      int value;
      int index;
    } local, global;
    local.value = (int)status;
    local.index = (int)session->leader_rank;
    if (MPI_Allreduce(&local, &global, 1, MPI_2INT, MPI_MAXLOC,
                      session->leader_communicator) != MPI_SUCCESS) {
      status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    } else {
      status = (MVMCKrylovStatus)global.value;
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        MPI_Allgather(local_digest, PROPOSAL_DIGEST_BYTES, MPI_BYTE,
                      session->chain_provenance_digests,
                      PROPOSAL_DIGEST_BYTES, MPI_BYTE,
                      session->leader_communicator) != MPI_SUCCESS) {
      status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
  }
  {
    int encoded = (int)status;
    const size_t bytes =
        session->chain_count * PROPOSAL_DIGEST_BYTES;
    if (MPI_Bcast(&encoded, 1, MPI_INT, 0,
                  session->chain_communicator) != MPI_SUCCESS) {
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    status = valid_status((MVMCKrylovStatus)encoded)
                 ? (MVMCKrylovStatus)encoded
                 : MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    if (status == MVMC_KRYLOV_STATUS_OK &&
        (bytes > (size_t)INT_MAX ||
         MPI_Bcast(session->chain_provenance_digests, (int)bytes, MPI_BYTE,
                   0, session->chain_communicator) != MPI_SUCCESS)) {
      status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
  }
#else
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memcpy(session->chain_provenance_digests, local_digest,
           PROPOSAL_DIGEST_BYTES);
  }
#endif
  return synchronize_status(session, status, 0);
}

static MVMCKrylovStatus coefficient_matrix_means(
    const MVMCPowerLanczosProductionSession *session, size_t excluded_block,
    double complex overlap[6], double complex forward[6],
    double complex reverse[6], double complex squared[6],
    uint64_t *sample_count) {
  uint64_t count = 0;
  size_t block;
  size_t entry;
  memset(overlap, 0, 6 * sizeof(*overlap));
  memset(forward, 0, 6 * sizeof(*forward));
  memset(reverse, 0, 6 * sizeof(*reverse));
  memset(squared, 0, 6 * sizeof(*squared));
  if (sample_count == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *sample_count = 0;
  if (!session->coefficient_evidence_valid ||
      session->coefficient_evidence_block_counts == NULL ||
      session->coefficient_evidence_block_sums == NULL) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  for (block = 0; block < session->execution.block_count; ++block) {
    const size_t offset = block * 4 * session->upper_count;
    if (block == excluded_block) continue;
    if (session->coefficient_evidence_block_counts[block] == 0 ||
        session->coefficient_evidence_block_counts[block] >
            UINT64_MAX - count) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    count += session->coefficient_evidence_block_counts[block];
    for (entry = 0; entry < session->upper_count; ++entry) {
      overlap[entry] +=
          session->coefficient_evidence_block_sums[offset + entry];
      forward[entry] +=
          session->coefficient_evidence_block_sums[
              offset + session->upper_count + entry];
      reverse[entry] += session->coefficient_evidence_block_sums[
          offset + 2 * session->upper_count + entry];
      squared[entry] += session->coefficient_evidence_block_sums[
          offset + 3 * session->upper_count + entry];
    }
  }
  if (count == 0 || count >
                        MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  for (entry = 0; entry < session->upper_count; ++entry) {
    const double denominator = (double)count;
    overlap[entry] /= denominator;
    forward[entry] /= denominator;
    reverse[entry] /= denominator;
    squared[entry] /= denominator;
    if (!isfinite(creal(overlap[entry])) ||
        !isfinite(cimag(overlap[entry])) ||
        !isfinite(creal(forward[entry])) ||
        !isfinite(cimag(forward[entry])) ||
        !isfinite(creal(reverse[entry])) ||
        !isfinite(cimag(reverse[entry])) ||
        !isfinite(creal(squared[entry])) ||
        !isfinite(cimag(squared[entry]))) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  *sample_count = count;
  return MVMC_KRYLOV_STATUS_OK;
}

static double complex packed_expectation(
    const double complex packed[6], size_t dimension,
    const double complex *left, const double complex *right) {
  double complex value = 0.0;
  size_t row;
  size_t column;
  size_t index = 0;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column) {
      const double complex upper = packed[index++];
      value += conj(left[row]) * upper * right[column];
      if (row != column) {
        value += conj(left[column]) * conj(upper) * right[row];
      }
    }
  }
  return value;
}

static MVMCKrylovStatus projective_distance(
    const double complex overlap[6], size_t dimension,
    const double complex *left, const double complex *right,
    double *distance) {
  const double complex left_norm_value =
      packed_expectation(overlap, dimension, left, left);
  const double complex right_norm_value =
      packed_expectation(overlap, dimension, right, right);
  const double complex cross =
      packed_expectation(overlap, dimension, left, right);
  const double scale = fmax(1.0, fmax(fabs(creal(left_norm_value)),
                                     fabs(creal(right_norm_value))));
  double ratio;
  double radicand;
  if (distance == NULL || !isfinite(creal(left_norm_value)) ||
      !isfinite(cimag(left_norm_value)) ||
      !isfinite(creal(right_norm_value)) ||
      !isfinite(cimag(right_norm_value)) || !isfinite(creal(cross)) ||
      !isfinite(cimag(cross)) || creal(left_norm_value) <= 0.0 ||
      creal(right_norm_value) <= 0.0 ||
      fabs(cimag(left_norm_value)) > 64.0 * DBL_EPSILON * scale ||
      fabs(cimag(right_norm_value)) > 64.0 * DBL_EPSILON * scale) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  ratio = (creal(cross) * creal(cross) + cimag(cross) * cimag(cross)) /
          (creal(left_norm_value) * creal(right_norm_value));
  if (!isfinite(ratio) || ratio < -64.0 * DBL_EPSILON ||
      ratio > 1.0 + 256.0 * DBL_EPSILON) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (ratio < 0.0) ratio = 0.0;
  if (ratio > 1.0) ratio = 1.0;
  radicand = 1.0 - ratio;
  *distance = sqrt(radicand);
  return isfinite(*distance) ? MVMC_KRYLOV_STATUS_OK
                             : MVMC_KRYLOV_STATUS_NONFINITE;
}

static int close_real(double left, double right) {
  const double scale = fmax(1.0, fmax(fabs(left), fabs(right)));
  return isfinite(left) && isfinite(right) &&
         fabs(left - right) <= 1024.0 * DBL_EPSILON * scale;
}

static int close_complex(double complex left, double complex right) {
  return close_real(creal(left), creal(right)) &&
         close_real(cimag(left), cimag(right));
}

static MVMCKrylovStatus freeze_observable_confirmation(
    MVMCPowerLanczosProductionSession *session) {
  MVMCPowerLanczosObservableConfirmationSummary summary;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  const size_t matrix_entries_per_block = 4 * session->upper_count;
  const size_t observable_entry_count =
      session->observable_request_count *
      MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES;
  size_t block;
  if (!session->execution.observable_enabled) {
    return MVMC_KRYLOV_STATUS_OK;
  }
  for (block = 0;
       status == MVMC_KRYLOV_STATUS_OK &&
       block < session->execution.block_count;
       ++block) {
    const double complex *base =
        session->global_block_sums +
        block * session->coefficient_entries_per_block;
    status = confirmation_status(
        mvmc_power_lanczos_observable_confirmation_add_coefficient_block(
            session->observable_confirmation,
            session->global_block_counts[block], base,
            base + session->upper_count,
            base + 2 * session->upper_count,
            base + 3 * session->upper_count,
            base + matrix_entries_per_block, observable_entry_count));
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = confirmation_status(
        mvmc_power_lanczos_observable_confirmation_freeze_coefficient(
            session->observable_confirmation));
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = confirmation_status(
        mvmc_power_lanczos_observable_confirmation_summary(
            session->observable_confirmation, &summary));
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (!summary.valid ||
       summary.status != MVMC_POWER_LANCZOS_CONFIRMATION_OK ||
       summary.state !=
           MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_FROZEN ||
       summary.request_count != session->observable_request_count ||
       summary.block_count != session->execution.block_count ||
       summary.coefficient_blocks_added != session->execution.block_count ||
       summary.coefficient_sample_count != session->coefficient_sample_count ||
       !close_complex(summary.coefficient[0],
                      session->coefficient_gevp.coefficient[0]) ||
       !close_complex(summary.coefficient[1],
                      session->coefficient_gevp.coefficient[1]) ||
       !close_real(summary.energy, session->coefficient_gevp.energy) ||
       summary.maximum_leave_one_projective_distance >
           session->execution.maximum_leave_one_projective_distance)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  return status;
}

static void block_payload_digest(
    const MVMCPowerLanczosProductionSession *session,
    unsigned char digest[PROPOSAL_DIGEST_BYTES]) {
  uint64_t hash[4];
  size_t index;
  const size_t entries_per_block = session->coefficient_entries_per_block;
  digest_initialize(hash, UINT64_C(0x5036434f45464642));
  digest_u64(hash, (uint64_t)session->execution.block_count);
  digest_u64(hash, (uint64_t)entries_per_block);
  for (index = 0; index < session->execution.block_count; ++index) {
    digest_u64(hash, session->global_block_counts[index]);
  }
  for (index = 0;
       index < session->execution.block_count * entries_per_block; ++index) {
    digest_u64(hash, double_bits(creal(session->global_block_sums[index])));
    digest_u64(hash, double_bits(cimag(session->global_block_sums[index])));
  }
  digest_finish(hash, digest);
}

static int append_u64(unsigned char **cursor, const unsigned char *end,
                      uint64_t value) {
  if (cursor == NULL || *cursor == NULL || end == NULL ||
      *cursor > end || (size_t)(end - *cursor) < sizeof(uint64_t)) {
    return 0;
  }
  store_u64_be(value, *cursor);
  *cursor += sizeof(uint64_t);
  return 1;
}

static int append_bytes(unsigned char **cursor, const unsigned char *end,
                        const void *bytes, size_t count) {
  if (cursor == NULL || *cursor == NULL || end == NULL || bytes == NULL ||
      *cursor > end || count > (size_t)(end - *cursor)) {
    return 0;
  }
  memcpy(*cursor, bytes, count);
  *cursor += count;
  return 1;
}

static MVMCKrylovStatus create_coefficient_provenance(
    MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosSnapshotSummary *snapshot_summary) {
  static const unsigned char domain[] =
      "MVMC-P6-PRODUCTION-COEFFICIENT-V1";
  unsigned char record[PROVENANCE_RECORD_CAPACITY];
  unsigned char block_digest[PROPOSAL_DIGEST_BYTES];
  unsigned char chain_digest[PROPOSAL_DIGEST_BYTES];
  unsigned char *cursor = record;
  const unsigned char *end = record + sizeof(record);
  uint64_t hash = 0;
  size_t index;
  memset(record, 0, sizeof(record));
  block_payload_digest(session, block_digest);
  if (!mvmc_power_lanczos_sha256_bytes(
          session->chain_provenance_digests,
          session->chain_count * PROPOSAL_DIGEST_BYTES, chain_digest)) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (!append_bytes(&cursor, end, domain, sizeof(domain)) ||
      !append_u64(&cursor, end,
                  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_VERSION) ||
      !append_u64(&cursor, end,
                  (uint64_t)(unsigned int)session->config.power_step) ||
      !append_u64(&cursor, end, (uint64_t)session->basis_count) ||
      !append_u64(&cursor, end, (uint64_t)session->upper_count) ||
      !append_u64(&cursor, end,
                  (uint64_t)session->execution.block_count) ||
      !append_u64(&cursor, end, session->coefficient_sample_count) ||
      !append_u64(&cursor, end, session->execution_fingerprint) ||
      !append_u64(&cursor, end,
                  mvmc_bounded_krylov_plan_hash(session->bounded_plan)) ||
      !append_u64(&cursor, end,
                  session->execution.amplitude_generation_hash) ||
      !append_u64(&cursor, end,
                  session->execution.coefficient_guide_policy.policy_hash) ||
      !append_u64(&cursor, end,
                  session->execution.proposal_policy.policy_hash) ||
      !append_u64(&cursor, end,
                  (uint64_t)(unsigned int)
                      session->execution.gevp_policy.cutoff_id) ||
      !append_u64(&cursor, end,
                  double_bits(session->execution.gevp_policy
                                  .rank_relative_cutoff)) ||
      !append_bytes(&cursor, end, snapshot_summary->base_sha256,
                    PROPOSAL_DIGEST_BYTES) ||
      !append_bytes(&cursor, end, chain_digest,
                    PROPOSAL_DIGEST_BYTES) ||
      !append_bytes(&cursor, end, block_digest,
                    PROPOSAL_DIGEST_BYTES) ||
      !append_u64(&cursor, end,
                  (uint64_t)(unsigned int)
                      session->coefficient_gevp.retained_rank) ||
      !append_u64(&cursor, end,
                  (uint64_t)(unsigned int)
                      session->coefficient_gevp.root_multiplicity) ||
      !append_u64(&cursor, end,
                  double_bits(session->coefficient_gevp.energy)) ||
      !append_u64(&cursor, end,
                  double_bits(session->coefficient_gevp
                                  .normwise_backward_error))) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  for (index = 0; index < session->basis_count; ++index) {
    if (!append_u64(&cursor, end,
                    double_bits(creal(
                        session->coefficient_gevp.coefficient[index]))) ||
        !append_u64(&cursor, end,
                    double_bits(cimag(
                        session->coefficient_gevp.coefficient[index])))) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  if (!mvmc_power_lanczos_sha256_bytes(
          record, (size_t)(cursor - record),
          session->coefficient_provenance_sha256)) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  for (index = 0; index < 8; ++index) {
    hash = (hash << 8) | session->coefficient_provenance_sha256[index];
  }
  if (hash == 0) {
    hash = UINT64_C(1469598103934665603);
    for (index = 0; index < PROPOSAL_DIGEST_BYTES; ++index) {
      hash ^= session->coefficient_provenance_sha256[index];
      hash *= UINT64_C(1099511628211);
    }
    if (hash == 0) hash = UINT64_C(1);
  }
  session->coefficient_provenance_hash = hash;
  memset(record, 0, sizeof(record));
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus audit_world_digest(
    MVMCPowerLanczosProductionSession *session,
    const unsigned char digest[PROPOSAL_DIGEST_BYTES]) {
#ifdef _mpi_use
  int lane;
  for (lane = 0; lane < 4; ++lane) {
    uint64_t value = 0;
    uint64_t minimum = 0;
    uint64_t maximum = 0;
    int byte;
    for (byte = 0; byte < 8; ++byte) {
      value = (value << 8) | digest[8 * lane + byte];
    }
    if (MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                      session->world_communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                      session->world_communicator) != MPI_SUCCESS) {
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    if (minimum != maximum) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
#else
  (void)session;
  (void)digest;
#endif
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_freeze_coefficient(
    MVMCPowerLanczosProductionSession *session) {
  MVMCPowerLanczosSnapshotSummary snapshot_summary;
  MVMCPowerLanczosGEVPResult leave_one;
  double complex overlap[6];
  double complex forward[6];
  double complex reverse[6];
  double complex squared[6];
  uint64_t sample_count = 0;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  size_t block;
  double maximum_distance = 0.0;
  if (session == NULL || !session->valid || !session->execution_prepared ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state !=
          MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_READY) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_snapshot_summary(session->snapshot,
                                               &snapshot_summary);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      snapshot_summary.state !=
          MVMC_POWER_LANCZOS_SNAPSHOT_COEFFICIENT_COMPLETE) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  {
    unsigned char zero_digest[PROPOSAL_DIGEST_BYTES] = {0};
    const unsigned char *digest =
        status == MVMC_KRYLOV_STATUS_OK
            ? snapshot_summary.coefficient_sha256
            : zero_digest;
    status = gather_coefficient_chain_digests(session, digest, status);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = coefficient_matrix_means(
        session, SIZE_MAX, overlap, forward, reverse, squared,
        &sample_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (sample_count != session->coefficient_sample_count ||
       mvmc_power_lanczos_gevp_solve_complex_packed(
           &session->execution.gevp_policy, (int)session->basis_count,
           overlap, forward, reverse, squared, session->upper_count,
           &session->coefficient_gevp) != MVMC_KRYLOV_GEVP_OK ||
       !session->coefficient_gevp.valid ||
       session->coefficient_gevp.retained_rank < 1 ||
       session->coefficient_gevp.retained_rank >
           (int)session->basis_count ||
       session->coefficient_gevp.root_multiplicity < 1 ||
       session->coefficient_gevp.root_multiplicity >
           session->coefficient_gevp.retained_rank ||
       session->coefficient_gevp.normwise_backward_error >
           session->execution.gevp_policy.maximum_normwise_backward_error)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memcpy(session->coefficient_overlap_mean, overlap,
           session->upper_count * sizeof(*overlap));
    memcpy(session->coefficient_hamiltonian_mean, forward,
           session->upper_count * sizeof(*forward));
    memcpy(session->coefficient_hamiltonian_adjoint_mean, reverse,
           session->upper_count * sizeof(*reverse));
    memcpy(session->coefficient_hamiltonian_squared_mean, squared,
           session->upper_count * sizeof(*squared));
  }
  for (block = 0;
       status == MVMC_KRYLOV_STATUS_OK &&
       block < session->execution.block_count;
       ++block) {
    double distance = NAN;
    double complex leave_one_second_moment = NAN + I * NAN;
    status = coefficient_matrix_means(session, block, overlap, forward,
                                      reverse, squared, &sample_count);
    if (status == MVMC_KRYLOV_STATUS_OK &&
        (mvmc_power_lanczos_gevp_solve_complex_packed(
             &session->execution.gevp_policy, (int)session->basis_count,
             overlap, forward, reverse, squared, session->upper_count,
             &leave_one) != MVMC_KRYLOV_GEVP_OK ||
         !leave_one.valid ||
         leave_one.retained_rank < 1 ||
         leave_one.retained_rank > (int)session->basis_count ||
         leave_one.root_multiplicity < 1 ||
         leave_one.root_multiplicity > leave_one.retained_rank)) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      const double complex normalization = packed_expectation(
          overlap, session->basis_count, leave_one.coefficient,
          leave_one.coefficient);
      const double complex squared_expectation = packed_expectation(
          squared, session->basis_count, leave_one.coefficient,
          leave_one.coefficient);
      if (sample_count == 0 || !isfinite(leave_one.energy) ||
          !isfinite(creal(normalization)) ||
          !isfinite(cimag(normalization)) || creal(normalization) <= 0.0 ||
          !isfinite(creal(squared_expectation)) ||
          !isfinite(cimag(squared_expectation))) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      } else {
        leave_one_second_moment = squared_expectation / normalization;
        if (!isfinite(creal(leave_one_second_moment)) ||
            !isfinite(cimag(leave_one_second_moment))) {
          status = MVMC_KRYLOV_STATUS_NONFINITE;
        }
      }
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = coefficient_matrix_means(
          session, SIZE_MAX, overlap, forward, reverse, squared,
          &sample_count);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = projective_distance(
          overlap, session->basis_count,
          session->coefficient_gevp.coefficient,
          leave_one.coefficient, &distance);
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        distance >
            session->execution.maximum_leave_one_projective_distance) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      size_t index;
      session->leave_one_projective_distances[block] = distance;
      session->leave_one_energies[block] = leave_one.energy;
      session->leave_one_second_moments[block] = leave_one_second_moment;
      if (distance > maximum_distance) maximum_distance = distance;
      for (index = 0; index < session->basis_count; ++index) {
        session->leave_one_coefficients[
            block * session->basis_count + index] =
            leave_one.coefficient[index];
      }
    }
  }
  session->maximum_leave_one_projective_distance = maximum_distance;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = freeze_observable_confirmation(session);
  }
  status = synchronize_status(session, status, 1);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = create_coefficient_provenance(session, &snapshot_summary);
  }
  status = synchronize_status(session, status, 1);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = audit_world_digest(
        session, session->coefficient_provenance_sha256);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_final_state_policy_create_scaled_basis(
        session->config.power_step,
        MVMC_KRYLOV_FINAL_COEFFICIENT_PRODUCTION_NOISY,
        session->coefficient_provenance_hash,
        session->coefficient_gevp.coefficient,
        session->execution.matrix_policy.log_basis_scale,
        session->basis_count,
        &session->final_policy);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (session->final_policy.source !=
           MVMC_KRYLOV_FINAL_COEFFICIENT_PRODUCTION_NOISY ||
       session->final_policy.provenance_hash !=
           session->coefficient_provenance_hash ||
       session->final_policy.policy_hash == 0)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  status = synchronize_status(session, status, 1);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  session->final_policy_hash = session->final_policy.policy_hash;
  session->state =
      MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_FROZEN;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus final_proposal_step(
    MVMCPowerLanczosProductionSession *session, uint64_t ordinal) {
  MVMCKrylovPositiveSamplerRng proposal_rng;
  MVMCPowerLanczosChainRngResult rng_result;
  MVMCKrylovPositiveSamplerProposalDrawResult draw;
  MVMCKrylovFinalStateStepResult step;
  MVMCKrylovStatus status;
  memset(&proposal_rng, 0, sizeof(proposal_rng));
  memset(&draw, 0, sizeof(draw));
  status = mvmc_power_lanczos_chain_rng_derive_proposal(
      session->chain_rng_workspace, &session->final_rng,
      &session->final_rank_rng, ordinal, &proposal_rng, &rng_result);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_draw_mixture_rng(
        &session->model, &session->execution.proposal_policy,
        session->final_words, session->word_count, &proposal_rng,
        session->proposal_words, session->word_count, &draw);
  }
  status = synchronize_status(session, status, 0);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = audit_proposal(session, session->final_words, &draw, ordinal);
  status = synchronize_status(session, status, 0);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = mvmc_krylov_final_state_sampler_step(
      session->bounded_workspace, &session->final_policy,
      session->final_words, session->word_count, &session->final_sampler,
      session->proposal_words, session->word_count,
      draw.log_proposal_ratio, draw.uniform, session->execution.amplitude,
      session->execution.amplitude_context, &step);
  return synchronize_status(session, status, 0);
}

static MVMCKrylovStatus extract_final_blocks(
    MVMCPowerLanczosProductionSession *session,
    MVMCKrylovStatus local_status) {
  const size_t entries_per_block = session->final_entries_per_block;
  size_t block;
  memset(session->local_block_counts, 0,
         session->execution.block_count *
             sizeof(*session->local_block_counts));
  memset(session->local_block_sums, 0,
         session->execution.block_count * session->block_entry_count *
             sizeof(*session->local_block_sums));
  if (local_status != MVMC_KRYLOV_STATUS_OK) return local_status;
  if (session->execution.observable_enabled) {
    MVMCPowerLanczosObservableBlockSummary summary;
    local_status = mvmc_power_lanczos_observable_block_summary(
        session->observable_final_blocks, &summary);
    if (local_status == MVMC_KRYLOV_STATUS_OK &&
        (!summary.valid || summary.status != MVMC_KRYLOV_STATUS_OK ||
         summary.stage != MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL ||
         summary.request_count != session->observable_request_count ||
         summary.entries_per_request != 1 ||
         summary.completed_block_count != session->execution.block_count ||
         summary.completed_sample_count !=
             (uint64_t)session->execution.final_sample_count ||
         summary.current_block_sample_count != 0 ||
         summary.discarded_partial_sample_count != 0)) {
      local_status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (local_status == MVMC_KRYLOV_STATUS_OK) {
      local_status = mvmc_power_lanczos_observable_block_export(
          session->observable_final_blocks,
          session->observable_final_block_sums,
          session->execution.block_count * session->observable_request_count,
          session->observable_block_counts,
          session->execution.block_count);
    }
    if (local_status != MVMC_KRYLOV_STATUS_OK) return local_status;
  }
  for (block = 0; block < session->execution.block_count; ++block) {
    const size_t output_offset = block * entries_per_block;
    session->local_block_counts[block] =
        (uint64_t)(session->execution.final_sample_count /
                   session->execution.block_count);
    if (session->execution.observable_enabled &&
        session->observable_block_counts[block] !=
            session->local_block_counts[block]) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    session->local_block_sums[output_offset] =
        streaming_value(session->final_local_accumulators + 2 * block);
    session->local_block_sums[output_offset + 1] =
        streaming_value(session->final_local_accumulators + 2 * block + 1);
    if (session->execution.observable_enabled) {
      memcpy(session->local_block_sums + output_offset + 2,
             session->observable_final_block_sums +
                 block * session->observable_request_count,
             session->observable_request_count *
                 sizeof(*session->local_block_sums));
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus final_observables(
    MVMCPowerLanczosProductionSession *session) {
  double complex energy_sum = 0.0;
  double complex diagnostic_sum = 0.0;
  uint64_t count = 0;
  double complex normalization;
  double complex squared_expectation;
  double complex quotient;
  double complex variance;
  size_t block;
  if (!session->final_evidence_valid ||
      session->final_evidence_block_counts == NULL ||
      session->final_evidence_block_sums == NULL) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  for (block = 0; block < session->execution.block_count; ++block) {
    const size_t offset = 2 * block;
    if (session->final_evidence_block_counts[block] == 0 ||
        session->final_evidence_block_counts[block] > UINT64_MAX - count) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    count += session->final_evidence_block_counts[block];
    energy_sum += session->final_evidence_block_sums[offset];
    diagnostic_sum += session->final_evidence_block_sums[offset + 1];
  }
  if (count == 0 || count !=
                        (uint64_t)session->execution.final_sample_count *
                            (uint64_t)session->chain_count) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  energy_sum /= (double)count;
  diagnostic_sum /= (double)count;
  session->final_sample_count = count;
  session->final_energy = creal(energy_sum);
  session->final_energy_imaginary = cimag(energy_sum);
  session->final_local_energy_abs_squared_diagnostic =
      creal(diagnostic_sum);
  if (!isfinite(creal(energy_sum)) || !isfinite(cimag(energy_sum)) ||
      !isfinite(creal(diagnostic_sum)) ||
      !isfinite(cimag(diagnostic_sum)) ||
      fabs(cimag(diagnostic_sum)) >
          64.0 * DBL_EPSILON * fmax(1.0, fabs(creal(diagnostic_sum))) ||
      creal(diagnostic_sum) < 0.0) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  normalization = packed_expectation(
      session->coefficient_overlap_mean, session->basis_count,
      session->coefficient_gevp.coefficient,
      session->coefficient_gevp.coefficient);
  squared_expectation = packed_expectation(
      session->coefficient_hamiltonian_squared_mean,
      session->basis_count, session->coefficient_gevp.coefficient,
      session->coefficient_gevp.coefficient);
  if (!isfinite(creal(normalization)) ||
      !isfinite(cimag(normalization)) || creal(normalization) <= 0.0 ||
      !isfinite(creal(squared_expectation)) ||
      !isfinite(cimag(squared_expectation))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  quotient = squared_expectation / normalization;
  variance = quotient - energy_sum * energy_sum;
  if (!isfinite(creal(quotient)) || !isfinite(cimag(quotient)) ||
      !isfinite(creal(variance)) || !isfinite(cimag(variance))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  session->corrected_variance = creal(variance);
  session->corrected_variance_imaginary = cimag(variance);
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus add_final_observable_blocks(
    MVMCPowerLanczosProductionSession *session) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  size_t block;
  if (!session->execution.observable_enabled) {
    return MVMC_KRYLOV_STATUS_OK;
  }
  for (block = 0;
       status == MVMC_KRYLOV_STATUS_OK &&
       block < session->execution.block_count;
       ++block) {
    const double complex *base =
        session->global_block_sums +
        block * session->final_entries_per_block;
    status = confirmation_status(
        mvmc_power_lanczos_observable_confirmation_add_final_block(
            session->observable_confirmation,
            session->global_block_counts[block], base + 2,
            session->observable_request_count));
  }
  return status;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_run_final(
    MVMCPowerLanczosProductionSession *session) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  uint64_t ordinal = 0;
  size_t step;
  size_t sample;
  const size_t block_length =
      session == NULL || !session->execution_prepared
          ? 0
          : session->execution.final_sample_count /
                session->execution.block_count;
  if (session == NULL || !session->valid || !session->execution_prepared ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state !=
          MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_FROZEN ||
      block_length == 0 ||
      session->final_policy.source !=
          MVMC_KRYLOV_FINAL_COEFFICIENT_PRODUCTION_NOISY ||
      session->final_policy.provenance_hash !=
          session->coefficient_provenance_hash ||
      session->final_policy.policy_hash != session->final_policy_hash) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  session->state = MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_RUNNING;
  status = mvmc_power_lanczos_snapshot_final_begin(
      session->snapshot, &session->final_words, &session->word_count);
  status = synchronize_status(session, status, 0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        session->bounded_workspace, session->execution.amplitude,
        session->execution.amplitude_context,
        session->execution.amplitude_generation_hash);
  }
  status = synchronize_status(session, status, 0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_rng_seed(
        session->final_rng.derived_seed, session->final_rng.stage_tag,
        &session->final_rank_rng);
  }
  status = synchronize_status(session, status, 0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_final_state_sampler_initialize(
        session->bounded_workspace, &session->final_policy,
        session->final_words, session->word_count,
        session->execution.amplitude, session->execution.amplitude_context,
        &session->final_sampler);
  }
  status = synchronize_status(session, status, 0);
  for (step = 0; step < session->execution.final_warm_up; ++step) {
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = final_proposal_step(session, ordinal);
      if (status == MVMC_KRYLOV_STATUS_OK) ++ordinal;
    }
  }
  for (sample = 0; sample < session->execution.final_sample_count; ++sample) {
    for (step = 0; step < session->execution.final_interval; ++step) {
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = final_proposal_step(session, ordinal);
        if (status == MVMC_KRYLOV_STATUS_OK) ++ordinal;
      }
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      const size_t block = sample / block_length;
      double complex energy = NAN + I * NAN;
      double complex diagnostic = NAN + I * NAN;
      MVMCPowerLanczosNumericEvidence energy_evidence;
      memset(&energy_evidence, 0, sizeof(energy_evidence));
      status = mvmc_krylov_final_state_observable_measure(
          MVMC_KRYLOV_FINAL_OBSERVABLE_ENERGY_COMPLEX,
          &session->final_sampler.final_state, &energy);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = mvmc_krylov_final_state_local_energy_evidence(
            &session->final_sampler.final_state, &energy_evidence);
      }
      if (status == MVMC_KRYLOV_STATUS_OK &&
          energy_evidence.value != energy) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = mvmc_krylov_final_state_observable_measure(
            MVMC_KRYLOV_FINAL_OBSERVABLE_LOCAL_ENERGY_ABS_SQUARED_DIAGNOSTIC,
            &session->final_sampler.final_state, &diagnostic);
      }
      if (status == MVMC_KRYLOV_STATUS_OK &&
          (block >= session->execution.block_count ||
           streaming_add(session->final_local_accumulators + 2 * block,
                         energy) != MVMC_KRYLOV_STATUS_OK ||
           streaming_add(session->final_local_accumulators + 2 * block + 1,
                         diagnostic) != MVMC_KRYLOV_STATUS_OK)) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      }
      if (status == MVMC_KRYLOV_STATUS_OK &&
          session->execution.observable_enabled) {
        const double raw_basis_log_scale[2] = {0.0, 0.0};
        MVMCPowerLanczosObservableSampleDiagnostics diagnostics;
        status = mvmc_power_lanczos_observable_final_sample_with_evidence(
            session->observable_evaluator, session->bounded_workspace,
            &session->observable_layout, &session->observable_plan,
            session->final_words, session->word_count,
            session->final_sampler.final_state.log_weight,
            raw_basis_log_scale, session->final_policy.coefficient,
            session->observable_final_sample,
            session->observable_request_count,
            session->observable_numeric_evidence_sample,
            session->observable_request_count, &diagnostics);
        if (status == MVMC_KRYLOV_STATUS_OK &&
            (!diagnostics.valid ||
             diagnostics.status != MVMC_KRYLOV_STATUS_OK ||
             diagnostics.request_count != session->observable_request_count)) {
          status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
        }
        if (status == MVMC_KRYLOV_STATUS_OK) {
          status = mvmc_power_lanczos_observable_block_add_sample(
              session->observable_final_blocks,
              session->observable_final_sample,
              session->observable_request_count);
        }
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = append_final_primitive_trace(
            session, sample, &energy_evidence);
      }
      status = synchronize_status(session, status, 0);
    }
  }
  session->final_proposals = ordinal;
  if (mvmc_bounded_krylov_session_is_active(session->bounded_workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(session->bounded_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  status = synchronize_status(session, status, 0);
  status = synchronize_status(session, status, 1);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = gather_primitive_trace(
        session, MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_FINAL,
        &session->final_local_trace, &session->final_trace,
        session->final_primitive_count,
        session->execution.final_sample_count, 602);
  }
  status = extract_final_blocks(session, status);
  status = audit_local_blocks(session, session->final_entries_per_block,
                              status);
  status = reduce_blocks(session, session->final_entries_per_block, status);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = persist_final_evidence(session);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_final_complete(session->snapshot);
  }
  status = synchronize_status(session, status, 1);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = final_observables(session);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = add_final_observable_blocks(session);
  }
  status = synchronize_status(session, status, 1);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  status = reduce_chain_leader_counter(
      session, session->final_sampler.accepted_generation,
      &session->final_accepted_steps);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  session->state = MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_READY;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_finalize(
    MVMCPowerLanczosProductionSession *session) {
  MVMCPowerLanczosSnapshotSummary snapshot_summary;
  MVMCKrylovStatus status;
  uint64_t expected_coefficient_proposals;
  uint64_t expected_final_proposals;
  uint64_t expected_scale_pilot_proposals;
  uint64_t expected_scored_coefficient_proposals;
  uint64_t total_scale_pilot_proposals;
  uint64_t total_scored_coefficient_proposals;
  uint64_t total_final_proposals;
  if (session == NULL || !session->valid || !session->execution_prepared ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state != MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_READY) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  expected_scale_pilot_proposals =
      (uint64_t)session->execution.scale_pilot_warm_up +
      (uint64_t)session->execution.scale_pilot_sample_count;
  expected_scored_coefficient_proposals =
      (uint64_t)session->execution.coefficient_warm_up +
      (uint64_t)session->execution.coefficient_sample_count *
          (uint64_t)session->execution.coefficient_interval;
  expected_coefficient_proposals = expected_scale_pilot_proposals +
                                   expected_scored_coefficient_proposals;
  expected_final_proposals =
      (uint64_t)session->execution.final_warm_up +
      (uint64_t)session->execution.final_sample_count *
          (uint64_t)session->execution.final_interval;
  if ((expected_scale_pilot_proposals != 0 &&
       (uint64_t)session->chain_count >
           UINT64_MAX / expected_scale_pilot_proposals) ||
      (expected_scored_coefficient_proposals != 0 &&
       (uint64_t)session->chain_count >
           UINT64_MAX / expected_scored_coefficient_proposals) ||
      (expected_final_proposals != 0 &&
       (uint64_t)session->chain_count >
           UINT64_MAX / expected_final_proposals)) {
    return mvmc_power_lanczos_production_session_fail_internal(
        session, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  }
  total_scale_pilot_proposals =
      expected_scale_pilot_proposals * (uint64_t)session->chain_count;
  total_scored_coefficient_proposals =
      expected_scored_coefficient_proposals *
      (uint64_t)session->chain_count;
  total_final_proposals =
      expected_final_proposals * (uint64_t)session->chain_count;
  status = mvmc_power_lanczos_snapshot_verify(session->snapshot);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_summary(session->snapshot,
                                                 &snapshot_summary);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (snapshot_summary.state !=
           MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_COMPLETE ||
       session->coefficient_proposals != expected_coefficient_proposals ||
       session->final_proposals != expected_final_proposals ||
       session->coefficient_rank_rng.draws != expected_coefficient_proposals ||
       session->final_rank_rng.draws != expected_final_proposals ||
       session->scale_pilot_accepted_steps > total_scale_pilot_proposals ||
       session->coefficient_accepted_steps >
           total_scored_coefficient_proposals ||
       session->final_accepted_steps > total_final_proposals ||
       expected_coefficient_proposals == 0 || expected_final_proposals == 0 ||
       session->coefficient_rank_rng.stream !=
           session->coefficient_rng.stage_tag ||
       session->final_rank_rng.stream != session->final_rng.stage_tag ||
       session->final_policy.source !=
           MVMC_KRYLOV_FINAL_COEFFICIENT_PRODUCTION_NOISY ||
       session->final_policy.provenance_hash !=
           session->coefficient_provenance_hash ||
       session->final_policy.policy_hash != session->final_policy_hash ||
       session->coefficient_provenance_hash == 0 ||
       !session->coefficient_evidence_valid ||
       !session->final_evidence_valid ||
       !isfinite(session->final_energy) ||
       !isfinite(session->final_energy_imaginary) ||
       !isfinite(session->corrected_variance) ||
       !isfinite(session->corrected_variance_imaginary) ||
       !isfinite(session->final_local_energy_abs_squared_diagnostic))) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      session->execution.observable_enabled) {
    MVMCPowerLanczosObservableConfirmationSummary summary;
    const size_t result_count = session->execution.block_count *
                                session->observable_request_count;
    const size_t coefficient_count =
        session->execution.block_count *
        MVMC_POWER_LANCZOS_CONFIRMATION_DIMENSION;
    status = confirmation_status(
        mvmc_power_lanczos_observable_confirmation_finalize(
            session->observable_confirmation,
            session->observable_coefficient_estimates,
            session->observable_request_count,
            session->observable_final_estimates,
            session->observable_request_count,
            session->observable_leave_one_estimates, result_count,
            session->observable_final_block_means, result_count,
            session->observable_leave_one_coefficients, coefficient_count,
            session->observable_leave_one_projective_distances,
            session->execution.block_count));
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = confirmation_status(
          mvmc_power_lanczos_observable_confirmation_summary(
              session->observable_confirmation, &summary));
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        (!summary.valid ||
         summary.status != MVMC_POWER_LANCZOS_CONFIRMATION_OK ||
         summary.state != MVMC_POWER_LANCZOS_CONFIRMATION_FINALIZED ||
         summary.request_count != session->observable_request_count ||
         summary.block_count != session->execution.block_count ||
         summary.coefficient_blocks_added !=
             session->execution.block_count ||
         summary.final_blocks_added != session->execution.block_count ||
         summary.coefficient_sample_count !=
             session->coefficient_sample_count ||
         summary.final_sample_count != session->final_sample_count)) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  status = synchronize_status(session, status, 1);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  session->state = MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED;
  return MVMC_KRYLOV_STATUS_OK;
}
