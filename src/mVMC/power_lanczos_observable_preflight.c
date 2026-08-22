#include "power_lanczos_observable_preflight.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

static int CheckedAddSize(size_t *total, size_t addition) {
  if (addition > SIZE_MAX - *total) return 0;
  *total += addition;
  return 1;
}

static int CheckedMultiplySize(size_t left, size_t right, size_t *product) {
  if (left != 0 && right > SIZE_MAX / left) return 0;
  *product = left * right;
  return 1;
}

static int CheckedAddDouble(double *total, double addition) {
  if (!isfinite(*total) || !isfinite(addition) || addition < 0.0 ||
      addition > 1.0e300 - *total) {
    return 0;
  }
  *total += addition;
  return isfinite(*total);
}

static MVMCPowerLanczosObservablePreflightStatus CommitReject(
    MVMCPowerLanczosObservablePreflightStatus status,
    MVMCPowerLanczosObservablePreflightResult *candidate,
    MVMCPowerLanczosObservablePreflightResult *result) {
  candidate->valid = 1;
  candidate->status = status;
  candidate->corrected_dispatch_enabled = 0;
  *result = *candidate;
  return status;
}

MVMCPowerLanczosObservablePreflightStatus
mvmc_power_lanczos_observable_preflight(
    const MVMCPowerLanczosObservablePreflightInput *input,
    MVMCPowerLanczosObservablePreflightResult *result) {
  MVMCPowerLanczosObservablePreflightResult candidate;
  size_t request_count = 0;
  size_t matrix_bytes;
  size_t target_bytes;
  size_t exact_bytes;
  size_t artifact_entries;
  size_t artifact_bytes;
  size_t rss_bytes;
  double per_source_wall;
  double wall;
  unsigned int requested_mask = 0;
  int family;
  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT;
  if (input == NULL || result == NULL || input->nsite <= 0 ||
      input->nsite > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE ||
      input->nsite_uc <= 0 || input->nsite_uc > input->nsite ||
      input->nsite_uc > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE_UC ||
      input->block_count == 0 || input->saved_source_count == 0 ||
      input->target_cache_bytes_per_target == 0 ||
      !isfinite(input->setup_wall_upper_seconds) ||
      input->setup_wall_upper_seconds < 0.0 ||
      !isfinite(input->per_source_intercept_upper_seconds) ||
      input->per_source_intercept_upper_seconds < 0.0 ||
      !isfinite(input->per_request_slope_upper_seconds) ||
      input->per_request_slope_upper_seconds < 0.0 ||
      (input->block_contract_passed != 0 &&
       input->block_contract_passed != 1)) {
    if (result != NULL) *result = candidate;
    return MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT;
  }
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    size_t family_limit;
    if (input->family_count[family] < 0 ||
        !CheckedAddSize(&request_count,
                        (size_t)input->family_count[family])) {
      *result = candidate;
      return MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT;
    }
    if (input->family_count[family] > 0) {
      requested_mask |= 1u << (unsigned int)family;
    }
    family_limit = (size_t)(family == 0 ? 2 : 6) *
                   (size_t)input->nsite_uc * (size_t)input->nsite;
    if ((size_t)input->family_count[family] > family_limit ||
        (family == 0 &&
         input->family_count[family] >
             MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ONE_BODY) ||
        (family != 0 &&
         input->family_count[family] >
             MVMC_POWER_LANCZOS_OBSERVABLE_MAX_QUARTIC)) {
      candidate.requested_family_mask = requested_mask;
      candidate.request_count = request_count;
      return CommitReject(MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INPUT_LIMIT,
                          &candidate, result);
    }
  }
  candidate.requested_family_mask = requested_mask;
  candidate.request_count = request_count;
  candidate.unique_target_upper = input->unique_target_upper;
  candidate.artifact_file_count = input->artifact_file_count;
  if (request_count == 0 ||
      request_count > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS ||
      input->unique_target_upper > request_count) {
    return CommitReject(MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INPUT_LIMIT,
                        &candidate, result);
  }
  if (!input->block_contract_passed ||
      input->block_count < MVMC_POWER_LANCZOS_OBSERVABLE_MIN_BLOCK_COUNT ||
      input->stage_file_count >
          MVMC_POWER_LANCZOS_OBSERVABLE_MAX_STAGE_FILES ||
      input->artifact_ledger_bytes >
          MVMC_POWER_LANCZOS_OBSERVABLE_MAX_LEDGER_BYTES) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_BLOCK_CONTRACT,
        &candidate, result);
  }
  if ((input->correctness_family_mask & requested_mask) != requested_mask) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_CERTIFICATE_UNAVAILABLE,
        &candidate, result);
  }
  if ((input->scaling_family_mask & requested_mask) != requested_mask) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_SCALING_UNAVAILABLE,
        &candidate, result);
  }
  if (!CheckedMultiplySize(request_count, (size_t)64, &matrix_bytes) ||
      !CheckedMultiplySize(input->unique_target_upper,
                           input->target_cache_bytes_per_target,
                           &target_bytes)) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
        &candidate, result);
  }
  exact_bytes = input->fixed_allocated_bytes;
  if (!CheckedAddSize(&exact_bytes, matrix_bytes) ||
      !CheckedAddSize(&exact_bytes, target_bytes) ||
      !CheckedAddSize(&exact_bytes, input->statistical_workspace_bytes)) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
        &candidate, result);
  }
  rss_bytes = exact_bytes;
  if (!CheckedAddSize(&rss_bytes, input->rss_headroom_bytes) ||
      !CheckedMultiplySize(request_count, input->block_count,
                           &artifact_entries) ||
      !CheckedMultiplySize(artifact_entries, (size_t)80,
                           &artifact_bytes)) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
        &candidate, result);
  }
  candidate.artifact_payload_bytes = artifact_bytes;
  if (!CheckedAddSize(&artifact_bytes, input->artifact_ledger_bytes)) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
        &candidate, result);
  }
  per_source_wall = input->per_source_intercept_upper_seconds;
  if (!CheckedAddDouble(
          &per_source_wall,
          input->per_request_slope_upper_seconds * (double)request_count)) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
        &candidate, result);
  }
  wall = input->setup_wall_upper_seconds;
  if ((double)input->saved_source_count > 1.0e300 /
          fmax(per_source_wall, 1.0) ||
      !CheckedAddDouble(
          &wall, (double)input->saved_source_count * per_source_wall)) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
        &candidate, result);
  }
  candidate.matrix_accumulator_bytes = matrix_bytes;
  candidate.target_cache_bytes = target_bytes;
  candidate.exact_allocated_upper_bytes = exact_bytes;
  candidate.peak_rss_upper_bytes = rss_bytes;
  candidate.artifact_upper_bytes = artifact_bytes;
  candidate.wall_upper_seconds = wall;
  if (exact_bytes > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_BYTES ||
      rss_bytes > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_BYTES ||
      artifact_bytes >
          MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ARTIFACT_BYTES ||
      input->artifact_file_count >
          MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ARTIFACT_FILES ||
      wall > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_WALL_SECONDS) {
    return CommitReject(
        MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
        &candidate, result);
  }
  candidate.valid = 1;
  candidate.status = MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_READY;
  candidate.corrected_dispatch_enabled = 1;
  *result = candidate;
  return candidate.status;
}

const char *mvmc_power_lanczos_observable_preflight_error(
    MVMCPowerLanczosObservablePreflightStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_READY:
      return "";
    case MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT:
      return "P6 INPUT REJECTED: OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT";
    case MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INPUT_LIMIT:
      return "P6 INPUT REJECTED: OBSERVABLE_INPUT_LIMIT";
    case MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_BLOCK_CONTRACT:
      return "P6 INPUT REJECTED: OBSERVABLE_BLOCK_CONTRACT_UNAVAILABLE";
    case MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_CERTIFICATE_UNAVAILABLE:
      return "P6 INPUT REJECTED: OBSERVABLE_CERTIFICATE_UNAVAILABLE";
    case MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_SCALING_UNAVAILABLE:
      return "P6 INPUT REJECTED: OBSERVABLE_SCALING_UNAVAILABLE";
    case MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT:
      return "P6 INPUT REJECTED: OBSERVABLE_RESOURCE_LIMIT";
  }
  return "P6 INPUT REJECTED: OBSERVABLE_PREFLIGHT_UNKNOWN";
}
