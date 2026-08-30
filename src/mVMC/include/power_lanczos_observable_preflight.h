#ifndef MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_H
#define MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_H

#include "power_lanczos_observable_census.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_OBSERVABLE_MIN_BLOCK_COUNT 32
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_STAGE_FILES 6
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_LEDGER_BYTES \
  ((size_t)16 * 1024 * 1024)
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_WALL_SECONDS 3600.0
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_BYTES \
  ((size_t)8 * 1024 * 1024 * 1024)
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ARTIFACT_BYTES \
  ((size_t)100 * 1024 * 1024 * 1024)
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ARTIFACT_FILES ((size_t)25000)
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NQP_FULL 8

typedef enum {
  MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_READY = 0,
  MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT,
  MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INPUT_LIMIT,
  MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_BLOCK_CONTRACT,
  MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_CERTIFICATE_UNAVAILABLE,
  MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_SCALING_UNAVAILABLE,
  MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT
} MVMCPowerLanczosObservablePreflightStatus;

typedef struct {
  int nsite;
  int nsite_uc;
  int nqp_full;
  int family_count[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT];
  size_t unique_target_upper;
  size_t block_count;
  uint64_t saved_source_count;
  size_t fixed_allocated_bytes;
  size_t target_cache_bytes_per_target;
  size_t statistical_workspace_bytes;
  size_t artifact_ledger_bytes;
  size_t artifact_file_count;
  size_t stage_file_count;
  size_t rss_headroom_bytes;
  double setup_wall_upper_seconds;
  double per_source_intercept_upper_seconds;
  double per_request_slope_upper_seconds;
  unsigned int correctness_family_mask;
  unsigned int scaling_family_mask;
  int block_contract_passed;
} MVMCPowerLanczosObservablePreflightInput;

typedef struct {
  int valid;
  MVMCPowerLanczosObservablePreflightStatus status;
  int corrected_dispatch_enabled;
  int nsite;
  int nsite_uc;
  int nqp_full;
  int family_count[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT];
  unsigned int requested_family_mask;
  size_t request_count;
  size_t unique_target_upper;
  size_t block_count;
  uint64_t saved_source_count;
  size_t matrix_accumulator_bytes;
  size_t target_cache_bytes;
  size_t exact_allocated_upper_bytes;
  size_t peak_rss_upper_bytes;
  size_t artifact_payload_bytes;
  size_t artifact_upper_bytes;
  size_t artifact_file_count;
  double wall_upper_seconds;
} MVMCPowerLanczosObservablePreflightResult;

/*
 * Pure, allocation-free corrected-dispatch selector.  A non-READY result
 * always keeps corrected dispatch disabled; no legacy fallback is implied.
 */
MVMCPowerLanczosObservablePreflightStatus
mvmc_power_lanczos_observable_preflight(
    const MVMCPowerLanczosObservablePreflightInput *input,
    MVMCPowerLanczosObservablePreflightResult *result);

const char *mvmc_power_lanczos_observable_preflight_error(
    MVMCPowerLanczosObservablePreflightStatus status);

#endif
