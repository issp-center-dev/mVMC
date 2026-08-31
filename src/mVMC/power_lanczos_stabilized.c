#include "power_lanczos_stabilized.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE) ||                     \
    !defined(MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
#error "power_lanczos_stabilized.c requires the production Krylov core"
#endif

#include "krylov_final_state_chain.h"
#include "krylov_gevp_solver.h"
#include "krylov_matrix_measurement.h"
#include "krylov_streaming_statistics.h"

#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum {
  MVMC_POWER_LANCZOS_PACKED_CAPACITY = 6,
  MVMC_POWER_LANCZOS_BOOTSTRAP_ATTEMPTS = 1024
};

typedef struct {
  uint64_t count;
  double complex overlap[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
  double complex hamiltonian[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
  double complex hamiltonian_adjoint[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
  double complex hamiltonian_squared[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
} MVMCPowerLanczosMatrixBlock;

static void invalidate_result(
    MVMCKrylovStatus status, MVMCPowerLanczosStabilizedResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->version = MVMC_POWER_LANCZOS_STABILIZED_VERSION;
  result->status = status;
  result->energy = NAN;
  result->energy_standard_error = NAN;
  result->variance = NAN;
  result->variance_standard_error = NAN;
  result->final_energy_imaginary = NAN;
  result->variance_imaginary = NAN;
  result->condition_estimate = NAN;
  result->gevp_residual = NAN;
  result->energy_tau_int = NAN;
  result->effective_sample_count = NAN;
}

static MVMCKrylovStatus synchronize_status(
    MVMCKrylovStatus local, MVMCClassicPfaffianCommunicator communicator) {
#ifdef _mpi_use
  int encoded = (int)local;
  int reduced = 0;
  if (MPI_Allreduce(&encoded, &reduced, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS ||
      reduced < (int)MVMC_KRYLOV_STATUS_OK ||
      reduced > (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return (MVMCKrylovStatus)reduced;
#else
  (void)communicator;
  return local;
#endif
}

static uint64_t mix_u64(uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static uint64_t domain_seed(uint64_t seed, uint64_t domain) {
  const uint64_t derived = mix_u64(seed ^ domain);
  return derived == 0 ? UINT64_C(1) : derived;
}

static uint64_t hash_u64(uint64_t hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static MVMCKrylovStatus map_gevp_status(MVMCKrylovGEVPStatus status) {
  switch (status) {
    case MVMC_KRYLOV_GEVP_OK:
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_KRYLOV_GEVP_INVALID_ARGUMENT:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    case MVMC_KRYLOV_GEVP_NONFINITE_INPUT:
      return MVMC_KRYLOV_STATUS_NONFINITE;
    case MVMC_KRYLOV_GEVP_LAPACK_FAILURE:
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    default:
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
}

static int checked_add_u64(uint64_t left, uint64_t right, uint64_t *sum) {
  if (sum == NULL || right > UINT64_MAX - left) return 0;
  *sum = left + right;
  return 1;
}

static MVMCKrylovStatus create_limits(
    const MVMCKrylovFockModel *model, int power_step,
    uint64_t generation_hash, MVMCKrylovBoundedLimits *limits) {
  const int maximum_order = power_step + 1;
  uint64_t frontier = 1;
  uint64_t nodes = 1;
  uint64_t transitions = 0;
  int depth;
  if (model == NULL || limits == NULL || model->term_count == 0 ||
      power_step < 1 || power_step > 2 || generation_hash == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (depth = 0; depth < maximum_order; ++depth) {
    uint64_t next;
    if ((uint64_t)model->term_count > UINT64_MAX / frontier) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    next = frontier * (uint64_t)model->term_count;
    if (!checked_add_u64(nodes, next, &nodes) ||
        !checked_add_u64(transitions, next, &transitions)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    frontier = next;
  }
  if (nodes > UINT64_C(10000000) ||
      transitions > UINT64_C(100000000)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  memset(limits, 0, sizeof(*limits));
  limits->amplitude_policy_hash = generation_hash;
  limits->cache_bytes = (size_t)4 * 1024 * 1024;
  limits->max_row_transitions = model->term_count;
  limits->max_workspace_bytes = (size_t)512 * 1024 * 1024;
  limits->max_node_expansions = nodes;
  limits->max_terminal_amplitude_calls = frontier;
  limits->max_total_row_transitions = transitions;
  limits->max_order = maximum_order;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus draw_bounded(
    void *context, size_t bound, size_t *value) {
  return mvmc_krylov_positive_sampler_rng_draw_bounded(
      (MVMCKrylovPositiveSamplerRng *)context, bound, value);
}

static MVMCKrylovStatus draw_uniform_sector(
    const MVMCKrylovFockModel *model,
    MVMCClassicPfaffianCommunicator communicator,
    MVMCKrylovPositiveSamplerRng *rng, uint64_t *words,
    size_t word_count) {
  MVMCKrylovFockUniformProposalResult proposal;
  MVMCKrylovStatus status;
  memset(&proposal, 0, sizeof(proposal));
  status = mvmc_krylov_fock_proposal_draw_uniform_sector(
      model, draw_bounded, rng, words, word_count, &proposal);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (!proposal.valid || proposal.status != MVMC_KRYLOV_STATUS_OK ||
       proposal.word_count != word_count)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  status = synchronize_status(status, communicator);
#ifdef _mpi_use
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (word_count > (size_t)INT_MAX ||
       MPI_Bcast(words, (int)word_count, MPI_UINT64_T, 0,
                 communicator) != MPI_SUCCESS)) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  return synchronize_status(status, communicator);
}

static MVMCKrylovStatus coefficient_step(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *guide,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *words, size_t word_count,
    MVMCKrylovPositiveSamplerSnapshot *snapshot,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words) {
  MVMCKrylovPositiveSamplerProposalStepResult step;
  return mvmc_krylov_positive_sampler_step_mixture_rng(
      workspace, guide, model, proposal_policy, words, word_count, snapshot,
      rng, amplitude, amplitude_context, proposal_words, word_count, &step);
}

static MVMCKrylovStatus final_step(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *words, size_t word_count,
    MVMCKrylovFinalStateSnapshot *snapshot,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words) {
  MVMCKrylovFinalStateChainStepResult step;
  return mvmc_krylov_final_state_chain_step_mixture_rng(
      workspace, policy, model, proposal_policy, words, word_count,
      snapshot, rng, amplitude, amplitude_context, proposal_words,
      word_count, &step);
}

static void log_sum_add(double value, double *maximum, double *scaled_sum) {
  if (value > *maximum) {
    *scaled_sum =
        *maximum == -INFINITY ? 1.0
                             : *scaled_sum * exp(*maximum - value) + 1.0;
    *maximum = value;
  } else {
    *scaled_sum += exp(value - *maximum);
  }
}

static uint64_t guide_policy_hash(
    int order, double eta, const double *log_basis_scale) {
  uint64_t hash = UINT64_C(1469598103934665603);
  uint64_t bits = 0;
  int index;
  hash = hash_u64(hash, UINT64_C(0x504c534c494d4731));
  hash = hash_u64(hash, (uint64_t)(unsigned)order);
  memcpy(&bits, &eta, sizeof(bits));
  hash = hash_u64(hash, bits);
  for (index = 0; index <= order; ++index) {
    memcpy(&bits, log_basis_scale + index, sizeof(bits));
    hash = hash_u64(hash, bits);
  }
  return hash == 0 ? UINT64_C(1) : hash;
}

static MVMCKrylovStatus scale_pilot(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    size_t pilot_warm_up, size_t pilot_sample_count,
    uint64_t *words, size_t word_count,
    MVMCKrylovPositiveSamplerSnapshot *snapshot,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, MVMCKrylovPositiveGuidePolicy *guide,
    MVMCKrylovMatrixMeasurementPolicy *matrix) {
  double maximum[MVMC_KRYLOV_MAX_ORDER + 1];
  double scaled_sum[MVMC_KRYLOV_MAX_ORDER + 1];
  size_t step;
  size_t sample;
  int index;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  for (index = 0; index <= guide->order; ++index) {
    maximum[index] = -INFINITY;
    scaled_sum[index] = 0.0;
  }
  for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                 step < pilot_warm_up; ++step) {
    status = coefficient_step(
        workspace, guide, model, proposal_policy, words, word_count,
        snapshot, rng, amplitude, amplitude_context, proposal_words);
  }
  for (sample = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   sample < pilot_sample_count; ++sample) {
    status = coefficient_step(
        workspace, guide, model, proposal_policy, words, word_count,
        snapshot, rng, amplitude, amplitude_context, proposal_words);
    for (index = 0; status == MVMC_KRYLOV_STATUS_OK &&
                    index <= guide->order; ++index) {
      const MVMCScaledComplex *value = snapshot->krylov.value + index;
      if (!mvmc_scaled_complex_is_valid(value) ||
          value->state == MVMC_SCALED_COMPLEX_NONFINITE) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      } else if (value->state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
        const double log_squared = 2.0 * value->log_abs;
        if (!isfinite(log_squared)) {
          status = MVMC_KRYLOV_STATUS_NONFINITE;
        } else {
          log_sum_add(log_squared, maximum + index, scaled_sum + index);
        }
      }
    }
  }
  for (index = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  index <= guide->order; ++index) {
    double scale = 0.0;
    if (maximum[index] != -INFINITY) {
      const double log_mean =
          maximum[index] + log(scaled_sum[index]) -
          log((double)pilot_sample_count);
      scale = -0.5 * log_mean;
      if (!isfinite(log_mean) || !isfinite(scale)) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      }
    }
    guide->lambda[index] = 1.0;
    guide->log_basis_scale[index] = scale;
    matrix->guide_lambda[index] = 1.0;
    matrix->target_weight[index] = 1.0;
    matrix->log_basis_scale[index] = scale;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    const double eta = 0x1p-40 * (double)(guide->order + 1);
    guide->eta = eta;
    matrix->eta = eta;
    guide->policy_hash =
        guide_policy_hash(guide->order, eta, guide->log_basis_scale);
    status = mvmc_krylov_positive_sampler_initialize(
        workspace, guide, words, word_count, amplitude, amplitude_context,
        snapshot);
  }
  return status;
}

static MVMCKrylovStatus initialize_matrix_blocks(
    size_t block_count, int dimension, size_t upper_count,
    MVMCKrylovMatrixMeasurementAccumulator *accumulators,
    MVMCKrylovStreamingComplexSum overlap[][6],
    MVMCKrylovStreamingComplexSum hamiltonian[][6],
    MVMCKrylovStreamingComplexSum hamiltonian_adjoint[][6],
    MVMCKrylovStreamingComplexSum hamiltonian_squared[][6]) {
  size_t block;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  for (block = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  block < block_count; ++block) {
    status =
        mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
            dimension, overlap[block], hamiltonian[block],
            hamiltonian_adjoint[block], hamiltonian_squared[block],
            upper_count, accumulators + block);
  }
  return status;
}

static MVMCKrylovStatus extract_matrix_block(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    MVMCPowerLanczosMatrixBlock *block) {
  size_t row;
  size_t column;
  size_t entry = 0;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (block == NULL || accumulator == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(block, 0, sizeof(*block));
  block->count = accumulator->sample_count;
  for (row = 0; status == MVMC_KRYLOV_STATUS_OK &&
                row < accumulator->dimension; ++row) {
    for (column = row; status == MVMC_KRYLOV_STATUS_OK &&
                         column < accumulator->dimension;
         ++column, ++entry) {
      MVMCKrylovJackknifeBlock extracted;
      status = mvmc_krylov_matrix_measurement_extract_block(
          accumulator, MVMC_KRYLOV_MATRIX_OVERLAP, row, column,
          &extracted);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        block->overlap[entry] = extracted.numerator;
        status = mvmc_krylov_matrix_measurement_extract_block(
            accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN, row, column,
            &extracted);
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        block->hamiltonian[entry] = extracted.numerator;
        status = mvmc_krylov_matrix_measurement_extract_block(
            accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN_ADJOINT,
            row, column, &extracted);
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        block->hamiltonian_adjoint[entry] = extracted.numerator;
        status = mvmc_krylov_matrix_measurement_extract_block(
            accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN_SQUARED,
            row, column, &extracted);
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        block->hamiltonian_squared[entry] = extracted.numerator;
      }
    }
  }
  return status;
}

static void matrix_block_add(
    MVMCPowerLanczosMatrixBlock *total,
    const MVMCPowerLanczosMatrixBlock *block, size_t upper_count) {
  size_t entry;
  total->count += block->count;
  for (entry = 0; entry < upper_count; ++entry) {
    total->overlap[entry] += block->overlap[entry];
    total->hamiltonian[entry] += block->hamiltonian[entry];
    total->hamiltonian_adjoint[entry] +=
        block->hamiltonian_adjoint[entry];
    total->hamiltonian_squared[entry] +=
        block->hamiltonian_squared[entry];
  }
}

static MVMCPowerLanczosMatrixBlock matrix_block_without(
    const MVMCPowerLanczosMatrixBlock *total,
    const MVMCPowerLanczosMatrixBlock *removed, size_t upper_count) {
  MVMCPowerLanczosMatrixBlock result = *total;
  size_t entry;
  result.count -= removed->count;
  for (entry = 0; entry < upper_count; ++entry) {
    result.overlap[entry] -= removed->overlap[entry];
    result.hamiltonian[entry] -= removed->hamiltonian[entry];
    result.hamiltonian_adjoint[entry] -=
        removed->hamiltonian_adjoint[entry];
    result.hamiltonian_squared[entry] -=
        removed->hamiltonian_squared[entry];
  }
  return result;
}

static double complex packed_quadratic(
    const double complex *packed, int dimension,
    const double complex *coefficient) {
  double complex value = 0.0;
  size_t entry = 0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column, ++entry) {
      value += conj(coefficient[row]) * packed[entry] *
               coefficient[column];
      if (row != column) {
        value += conj(coefficient[column]) * conj(packed[entry]) *
                 coefficient[row];
      }
    }
  }
  return value;
}

static MVMCKrylovStatus solve_gevp(
    const MVMCKrylovGEVPPolicy *policy,
    const MVMCPowerLanczosMatrixBlock *matrix, int dimension,
    size_t upper_count, MVMCKrylovGEVPResult *result) {
  MVMCKrylovGEVPStatus status =
      mvmc_krylov_gevp_solve_complex_packed(
          policy, dimension, matrix->overlap, matrix->hamiltonian,
          matrix->hamiltonian_adjoint, matrix->hamiltonian_squared,
          upper_count, result);
  return map_gevp_status(status);
}

static MVMCKrylovStatus physical_second_moment(
    const MVMCPowerLanczosMatrixBlock *matrix,
    const MVMCKrylovGEVPResult *gevp, int dimension,
    double complex *second_moment) {
  const double complex normalization =
      packed_quadratic(matrix->overlap, dimension, gevp->coefficient);
  const double complex numerator =
      packed_quadratic(
          matrix->hamiltonian_squared, dimension, gevp->coefficient);
  if (second_moment == NULL ||
      !isfinite(creal(normalization)) ||
      !isfinite(cimag(normalization)) ||
      !(creal(normalization) > 0.0) ||
      fabs(cimag(normalization)) >
          1.0e-10 * fmax(1.0, fabs(creal(normalization))) ||
      !isfinite(creal(numerator)) || !isfinite(cimag(numerator))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  *second_moment = numerator / normalization;
  return isfinite(creal(*second_moment)) &&
                 isfinite(cimag(*second_moment))
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

static double jackknife_standard_error(
    const double *leave_one, size_t block_count) {
  double mean = 0.0;
  double sum = 0.0;
  size_t block;
  for (block = 0; block < block_count; ++block) mean += leave_one[block];
  mean /= (double)block_count;
  for (block = 0; block < block_count; ++block) {
    const double delta = leave_one[block] - mean;
    sum += delta * delta;
  }
  return sqrt(((double)block_count - 1.0) /
              (double)block_count * sum);
}

static double block_mean_standard_error(
    const double complex *block_sum, const uint64_t *block_sample_count,
    size_t block_count, double overall_mean) {
  double squared_deviation = 0.0;
  size_t block;
  for (block = 0; block < block_count; ++block) {
    const double block_mean =
        creal(block_sum[block]) / (double)block_sample_count[block];
    const double delta = block_mean - overall_mean;
    squared_deviation += delta * delta;
  }
  return sqrt(squared_deviation /
              ((double)block_count * ((double)block_count - 1.0)));
}

MVMCKrylovStatus mvmc_power_lanczos_stabilized_run(
    const MVMCPowerLanczosStabilizedInput *input,
    MVMCPowerLanczosStabilizedResult *result) {
  MVMCPowerLanczosClassicBridge *bridge = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovScaledAmplitudeCallback amplitude = NULL;
  void *amplitude_context = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovPositiveGuidePolicy guide;
  MVMCKrylovMatrixMeasurementPolicy matrix_policy;
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovPositiveSamplerRng bootstrap_rng;
  MVMCKrylovPositiveSamplerRng coefficient_rng;
  MVMCKrylovPositiveSamplerRng final_rng;
  MVMCKrylovPositiveSamplerSnapshot coefficient_snapshot;
  MVMCKrylovFinalStateSnapshot final_snapshot;
  MVMCKrylovFinalStatePolicy final_policy;
  MVMCKrylovGEVPPolicy gevp_policy;
  MVMCKrylovGEVPResult gevp;
  MVMCKrylovMatrixMeasurementAccumulator accumulators[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  MVMCKrylovStreamingComplexSum overlap_entries[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS][6];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS][6];
  MVMCKrylovStreamingComplexSum hamiltonian_adjoint_entries[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS][6];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS][6];
  MVMCPowerLanczosMatrixBlock matrix_blocks[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  MVMCPowerLanczosMatrixBlock matrix_total;
  double complex final_energy_block_sum[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  uint64_t final_energy_block_count[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  double leave_one_second_moment[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  double *energy_series = NULL;
  uint64_t *words = NULL;
  uint64_t *proposal_words = NULL;
  size_t dimension = 0;
  size_t upper_count = 0;
  size_t word_count = 0;
  size_t block_count = 0;
  size_t pilot_warm_up = 0;
  size_t pilot_sample_count = 0;
  uint64_t generation_hash = 0;
  uint64_t provenance_hash = 0;
  double complex energy_sum = 0.0;
  double complex energy = NAN + I * NAN;
  double complex second_moment = NAN + I * NAN;
  double complex variance = NAN + I * NAN;
  double energy_standard_error = NAN;
  double second_moment_standard_error = NAN;
  MVMCKrylovTauIntResult tau;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  size_t step;
  size_t sample;
  size_t block;
  int index;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (input == NULL || input->classic_view == NULL ||
      (input->power_step != 1 && input->power_step != 2) ||
      input->seed == 0 || input->sample_count < 8 || input->interval == 0 ||
      input->sample_count > SIZE_MAX / sizeof(*energy_series)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  dimension = (size_t)input->power_step + 1;
  status = mvmc_krylov_streaming_upper_count(dimension, &upper_count);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      upper_count > MVMC_POWER_LANCZOS_PACKED_CAPACITY) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  block_count =
      mvmc_power_lanczos_stabilized_block_count(input->sample_count);
  if (block_count == 0) status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;

  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_classic_bridge_create(
        input->classic_view, input->communicator, input->communicator,
        &bridge);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    model = mvmc_power_lanczos_classic_bridge_model(bridge);
    amplitude = mvmc_power_lanczos_classic_bridge_amplitude(bridge);
    amplitude_context =
        mvmc_power_lanczos_classic_bridge_amplitude_context(bridge);
    generation_hash =
        mvmc_power_lanczos_classic_bridge_generation_hash(bridge);
    if (model == NULL || amplitude == NULL || amplitude_context == NULL ||
        generation_hash == 0) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    word_count = mvmc_krylov_fock_word_count(model->site_count);
    if (word_count == 0 ||
        word_count > SIZE_MAX / sizeof(*words)) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    words = (uint64_t *)calloc(word_count, sizeof(*words));
    proposal_words =
        (uint64_t *)calloc(word_count, sizeof(*proposal_words));
    energy_series =
        (double *)calloc(input->sample_count, sizeof(*energy_series));
    if (words == NULL || proposal_words == NULL || energy_series == NULL) {
      status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = create_limits(
        model, input->power_step, generation_hash, &limits);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_plan_create(model, &limits, &plan);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_proposal_policy_create(
        1, 16, &proposal_policy);
  }
  memset(&guide, 0, sizeof(guide));
  memset(&matrix_policy, 0, sizeof(matrix_policy));
  guide.order = (int)dimension;
  guide.eta = 0x1p-40;
  matrix_policy.order = (int)dimension;
  matrix_policy.eta = 0x1p-40;
  for (index = 0; index <= (int)dimension; ++index) {
    const double initial_weight = index == 0 ? 0x1p-40 : 1.0;
    guide.lambda[index] = initial_weight;
    matrix_policy.guide_lambda[index] = initial_weight;
    matrix_policy.target_weight[index] = 1.0;
  }
  guide.policy_hash =
      guide_policy_hash(guide.order, guide.eta, guide.log_basis_scale);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_rng_seed(
        domain_seed(input->seed, UINT64_C(0x504c424f4f545354)),
        UINT64_C(0x504c424f4f545354), &bootstrap_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = draw_uniform_sector(
        model, input->communicator, &bootstrap_rng, words, word_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        workspace, amplitude, amplitude_context, generation_hash);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_rng_seed(
        domain_seed(input->seed, UINT64_C(0x504c434f45464631)),
        UINT64_C(0x504c434f45464631), &coefficient_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_initialize(
        workspace, &guide, words, word_count, amplitude,
        amplitude_context, &coefficient_snapshot);
  }
  pilot_warm_up = input->warm_up < 256 ? input->warm_up : 256;
  pilot_sample_count =
      input->sample_count < 256 ? input->sample_count : 256;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = scale_pilot(
        workspace, model, &proposal_policy, pilot_warm_up,
        pilot_sample_count, words, word_count, &coefficient_snapshot,
        &coefficient_rng, amplitude, amplitude_context, proposal_words,
        &guide, &matrix_policy);
  }
  for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                 step < input->warm_up; ++step) {
    status = coefficient_step(
        workspace, &guide, model, &proposal_policy, words, word_count,
        &coefficient_snapshot, &coefficient_rng, amplitude,
        amplitude_context, proposal_words);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = initialize_matrix_blocks(
        block_count, (int)dimension, upper_count, accumulators,
        overlap_entries, hamiltonian_entries,
        hamiltonian_adjoint_entries, hamiltonian_squared_entries);
  }
  for (sample = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   sample < input->sample_count; ++sample) {
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    block = sample * block_count / input->sample_count;
    for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   step < input->interval; ++step) {
      status = coefficient_step(
          workspace, &guide, model, &proposal_policy, words, word_count,
          &coefficient_snapshot, &coefficient_rng, amplitude,
          amplitude_context, proposal_words);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_krylov_matrix_measurement_accumulator_add_sample(
          accumulators + block, &matrix_policy,
          coefficient_snapshot.krylov.value,
          (size_t)coefficient_snapshot.krylov.evaluated_order + 1,
          &diagnostics);
    }
  }
  if (workspace != NULL &&
      mvmc_bounded_krylov_session_is_active(workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  status = synchronize_status(status, input->communicator);

  memset(&matrix_total, 0, sizeof(matrix_total));
  for (block = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  block < block_count; ++block) {
    status = extract_matrix_block(
        accumulators + block, matrix_blocks + block);
    if (status == MVMC_KRYLOV_STATUS_OK &&
        matrix_blocks[block].count == 0) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      matrix_block_add(
          &matrix_total, matrix_blocks + block, upper_count);
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      matrix_total.count != (uint64_t)input->sample_count) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      mvmc_krylov_gevp_default_policy(0x1p-40, &gevp_policy) !=
          MVMC_KRYLOV_GEVP_OK) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = solve_gevp(
        &gevp_policy, &matrix_total, (int)dimension, upper_count, &gevp);
  }
  provenance_hash = hash_u64(generation_hash, input->seed);
  if (provenance_hash == 0) provenance_hash = UINT64_C(1);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_final_state_policy_create_scaled_basis(
        input->power_step,
        MVMC_KRYLOV_FINAL_COEFFICIENT_PRODUCTION_NOISY,
        provenance_hash, gevp.coefficient, matrix_policy.log_basis_scale,
        dimension, &final_policy);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_rng_seed(
        domain_seed(input->seed, UINT64_C(0x504c46494e414c31)),
        UINT64_C(0x504c46494e414c31), &final_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        workspace, amplitude, amplitude_context, generation_hash);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    size_t attempt;
    status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    for (attempt = 0;
         status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
         attempt < MVMC_POWER_LANCZOS_BOOTSTRAP_ATTEMPTS;
         ++attempt) {
      status = draw_uniform_sector(
          model, input->communicator, &final_rng, words, word_count);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = mvmc_krylov_final_state_sampler_initialize(
            workspace, &final_policy, words, word_count, amplitude,
            amplitude_context, &final_snapshot);
      }
    }
  }
  for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                 step < input->warm_up; ++step) {
    status = final_step(
        workspace, &final_policy, model, &proposal_policy, words,
        word_count, &final_snapshot, &final_rng, amplitude,
        amplitude_context, proposal_words);
  }
  memset(final_energy_block_sum, 0, sizeof(final_energy_block_sum));
  memset(final_energy_block_count, 0, sizeof(final_energy_block_count));
  for (sample = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   sample < input->sample_count; ++sample) {
    double complex local_energy = NAN + I * NAN;
    MVMCScaledComplexExportStatus export_status;
    block = sample * block_count / input->sample_count;
    for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   step < input->interval; ++step) {
      status = final_step(
          workspace, &final_policy, model, &proposal_policy, words,
          word_count, &final_snapshot, &final_rng, amplitude,
          amplitude_context, proposal_words);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      export_status = mvmc_krylov_final_state_local_energy_export(
          &final_snapshot.final_state, &local_energy);
      if ((export_status != MVMC_SCALED_EXPORT_OK &&
           export_status != MVMC_SCALED_EXPORT_EXACT_ZERO) ||
          !isfinite(creal(local_energy)) ||
          !isfinite(cimag(local_energy))) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      }
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      energy_series[sample] = creal(local_energy);
      energy_sum += local_energy;
      final_energy_block_sum[block] += local_energy;
      ++final_energy_block_count[block];
    }
  }
  if (workspace != NULL &&
      mvmc_bounded_krylov_session_is_active(workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  status = synchronize_status(status, input->communicator);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    energy = energy_sum / (double)input->sample_count;
    status = physical_second_moment(
        &matrix_total, &gevp, (int)dimension, &second_moment);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    variance = second_moment - energy * energy;
    if (!isfinite(creal(variance)) || !isfinite(cimag(variance))) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  for (block = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  block < block_count; ++block) {
    const MVMCPowerLanczosMatrixBlock leave_matrix =
        matrix_block_without(&matrix_total, matrix_blocks + block,
                             upper_count);
    MVMCKrylovGEVPResult leave_gevp;
    double complex leave_second = NAN + I * NAN;
    status = solve_gevp(
        &gevp_policy, &leave_matrix, (int)dimension, upper_count,
        &leave_gevp);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = physical_second_moment(
          &leave_matrix, &leave_gevp, (int)dimension, &leave_second);
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        !isfinite(creal(leave_second))) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    }
    leave_one_second_moment[block] = creal(leave_second);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    energy_standard_error = block_mean_standard_error(
        final_energy_block_sum, final_energy_block_count, block_count,
        creal(energy));
    second_moment_standard_error = jackknife_standard_error(
        leave_one_second_moment, block_count);
    if (!isfinite(energy_standard_error) ||
        !isfinite(second_moment_standard_error)) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_tau_int_geyer_initial_positive(
        energy_series, input->sample_count, &tau);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    result->valid = 1;
    result->status = MVMC_KRYLOV_STATUS_OK;
    result->power_step = input->power_step;
    result->block_count = block_count;
    result->coefficient_samples = (uint64_t)input->sample_count;
    result->final_samples = (uint64_t)input->sample_count;
    result->retained_rank = gevp.retained_rank;
    result->condition_estimate = gevp.condition_estimate;
    result->gevp_residual = gevp.gevp_relative_residual;
    result->energy = creal(energy);
    result->energy_standard_error = energy_standard_error;
    result->variance = creal(variance);
    result->variance_standard_error = hypot(
        second_moment_standard_error,
        2.0 * cabs(energy) * energy_standard_error);
    result->final_energy_imaginary = cimag(energy);
    result->variance_imaginary = cimag(variance);
    result->energy_tau_int = tau.tau_int;
    result->effective_sample_count =
        (double)input->sample_count / (2.0 * tau.tau_int);
  } else {
    invalidate_result(status, result);
  }

  free(energy_series);
  free(proposal_words);
  free(words);
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_power_lanczos_classic_bridge_destroy(bridge);
  return status;
}
