#include "power_lanczos_observable_preflight.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

static void Check(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "FAIL: %s\n", message);
    ++failures;
  }
}

static MVMCPowerLanczosObservablePreflightInput ValidInput(void) {
  MVMCPowerLanczosObservablePreflightInput input;
  memset(&input, 0, sizeof(input));
  input.nsite = 16;
  input.nsite_uc = 16;
  input.nqp_full = 8;
  input.family_count[0] = 32;
  input.family_count[1] = 64;
  input.family_count[2] = 64;
  input.unique_target_upper = 120;
  input.block_count = 32;
  input.saved_source_count = 1024;
  input.fixed_allocated_bytes = 1024 * 1024;
  input.target_cache_bytes_per_target = 256;
  input.statistical_workspace_bytes = 2 * 1024 * 1024;
  input.artifact_ledger_bytes = 1024 * 1024;
  input.artifact_file_count = 6;
  input.stage_file_count = 2;
  input.rss_headroom_bytes = 8 * 1024 * 1024;
  input.setup_wall_upper_seconds = 2.0;
  input.per_source_intercept_upper_seconds = 0.001;
  input.per_request_slope_upper_seconds = 0.00001;
  input.correctness_family_mask = 7;
  input.scaling_family_mask = 7;
  input.block_contract_passed = 1;
  return input;
}

static void ExpectStatus(
    const MVMCPowerLanczosObservablePreflightInput *input,
    MVMCPowerLanczosObservablePreflightStatus expected,
    const char *message) {
  MVMCPowerLanczosObservablePreflightResult result;
  MVMCPowerLanczosObservablePreflightStatus status;
  memset(&result, 0xa5, sizeof(result));
  status = mvmc_power_lanczos_observable_preflight(input, &result);
  Check(status == expected, message);
  if (expected == MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT) {
    Check(!result.valid && !result.corrected_dispatch_enabled, message);
  } else {
    Check(result.valid && result.status == expected, message);
    Check(result.corrected_dispatch_enabled ==
              (expected == MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_READY),
          message);
  }
}

int main(void) {
  MVMCPowerLanczosObservablePreflightInput input = ValidInput();
  MVMCPowerLanczosObservablePreflightResult result;
  MVMCPowerLanczosObservablePreflightStatus status;
  status = mvmc_power_lanczos_observable_preflight(&input, &result);
  Check(status == MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_READY &&
            result.valid && result.corrected_dispatch_enabled,
        "valid preflight not enabled");
  Check(result.request_count == 160 &&
            result.nsite == 16 && result.nsite_uc == 16 &&
            result.nqp_full == 8 && result.family_count[0] == 32 &&
            result.family_count[1] == 64 && result.family_count[2] == 64 &&
            result.block_count == 32 && result.saved_source_count == 1024 &&
            result.matrix_accumulator_bytes == 64 * 160 &&
            result.target_cache_bytes == 120 * 256 &&
            result.artifact_payload_bytes == 80 * 32 * 160,
        "preflight exact formulas");

  input = ValidInput();
  input.nqp_full = 0;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT,
               "zero NQPFull");
  input = ValidInput();
  input.nqp_full = -1;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT,
               "negative NQPFull");
  input = ValidInput();
  input.nqp_full = 9;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT,
               "NQPFull over cap");

  input = ValidInput();
  input.correctness_family_mask = 3;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_CERTIFICATE_UNAVAILABLE,
               "missing correctness certificate");
  Check(strcmp(mvmc_power_lanczos_observable_preflight_error(
                   MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_CERTIFICATE_UNAVAILABLE),
               "P6 INPUT REJECTED: OBSERVABLE_CERTIFICATE_UNAVAILABLE") == 0,
        "frozen certificate error string");

  input = ValidInput();
  input.scaling_family_mask = 5;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_SCALING_UNAVAILABLE,
               "missing scaling identity");
  input = ValidInput();
  input.block_count = 31;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_BLOCK_CONTRACT,
               "short block count");
  input = ValidInput();
  input.family_count[0] = 513;
  ExpectStatus(&input, MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INPUT_LIMIT,
               "family count over cap");
  input = ValidInput();
  input.unique_target_upper = 161;
  ExpectStatus(&input, MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INPUT_LIMIT,
               "unique target over request count");
  input = ValidInput();
  input.fixed_allocated_bytes =
      MVMC_POWER_LANCZOS_OBSERVABLE_MAX_BYTES;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
               "allocation cap");
  input = ValidInput();
  input.per_request_slope_upper_seconds = 1.0;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_RESOURCE_LIMIT,
               "wall cap");
  input = ValidInput();
  input.family_count[0] = -1;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT,
               "negative family count");
  input = ValidInput();
  input.setup_wall_upper_seconds = NAN;
  ExpectStatus(&input,
               MVMC_POWER_LANCZOS_OBSERVABLE_PREFLIGHT_INVALID_ARGUMENT,
               "nonfinite timing");

  if (failures != 0) {
    fprintf(stderr, "%d observable preflight assertion(s) failed\n",
            failures);
    return EXIT_FAILURE;
  }
  puts("power_lanczos_observable_preflight_unit: PASS");
  return EXIT_SUCCESS;
}
