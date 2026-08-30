#ifndef MVMC_POWER_LANCZOS_STABILIZATION_STATISTICS_H
#define MVMC_POWER_LANCZOS_STABILIZATION_STATISTICS_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "power_lanczos_production_session.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_STABILIZATION_RESULT_VERSION UINT64_C(1)
#define MVMC_POWER_LANCZOS_STABILIZATION_IMAGINARY_SIGMA_MULTIPLIER 5.0

typedef enum {
  MVMC_POWER_LANCZOS_STABILIZATION_PASS = 0,
  MVMC_POWER_LANCZOS_STABILIZATION_INCONCLUSIVE = 1,
  MVMC_POWER_LANCZOS_STABILIZATION_FAIL = 2
} MVMCPowerLanczosStabilizationDecision;

enum {
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_BLOCK_COUNT = UINT64_C(1) << 0,
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_COEFFICIENT_TAU = UINT64_C(1) << 1,
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_FINAL_TAU = UINT64_C(1) << 2,
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_GEVP_RANK = UINT64_C(1) << 3,
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_GEVP_RESIDUAL = UINT64_C(1) << 4,
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_GEVP_NORMALIZATION = UINT64_C(1) << 5,
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_NEGATIVE_VARIANCE = UINT64_C(1) << 6,
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_ZERO_SUPPORT = UINT64_C(1) << 7,
  MVMC_POWER_LANCZOS_STABILIZATION_GATE_COMPLEX_RESIDUAL = UINT64_C(1) << 8
};

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  MVMCPowerLanczosStabilizationDecision decision;
  uint64_t failed_gates;
  size_t block_count;
  size_t chain_count;
  uint64_t coefficient_sample_count;
  uint64_t final_sample_count;
  size_t coefficient_block_length;
  size_t final_block_length;
  double coefficient_tau_int_max;
  double final_tau_int_max;
  double coefficient_tau_to_block_length;
  double final_tau_to_block_length;
  double coefficient_effective_block_count;
  double final_effective_block_count;
  double energy;
  double energy_imaginary;
  double energy_standard_error;
  double energy_imaginary_standard_error;
  double coefficient_second_moment;
  double coefficient_second_moment_imaginary;
  double coefficient_second_moment_standard_error;
  double coefficient_second_moment_imaginary_standard_error;
  double variance;
  double variance_imaginary;
  double variance_standard_error;
  double maximum_leave_one_projective_distance;
  int retained_rank;
  double gevp_relative_residual;
  double gevp_normalization_absolute_error;
  uint64_t finite_nonzero_primitive_count;
  uint64_t exact_zero_primitive_count;
  uint64_t numeric_zero_primitive_count;
  uint64_t outside_base_support_samples;
  uint64_t nonzero_outside_support_contributions;
  double maximum_absolute_numeric_bound;
  size_t session_allocated_bytes;
} MVMCPowerLanczosStabilizationResult;

/*
 * Evaluate immutable evidence owned by world rank 0.  Sample-level trace data
 * are consumed in memory and never written.  A non-root session intentionally
 * returns INVALID_ARGUMENT because the gathered trace is root-owned.
 */
MVMCKrylovStatus mvmc_power_lanczos_stabilization_statistics_evaluate(
    const MVMCPowerLanczosProductionSession *session,
    int require_zero_support_evidence,
    MVMCPowerLanczosStabilizationResult *result);

const char *mvmc_power_lanczos_stabilization_decision_string(
    MVMCPowerLanczosStabilizationDecision decision);

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_POWER_LANCZOS_STABILIZATION_STATISTICS_H */
