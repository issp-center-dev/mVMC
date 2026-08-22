#include "power_lanczos_selector.h"

#include <limits.h>
#include <string.h>

static const char *const ChainControlNames[] = {
    "NLanczosCoeffWarmUp", "NLanczosCoeffSample",
    "NLanczosCoeffInterval", "NLanczosFinalWarmUp",
    "NLanczosFinalSample", "NLanczosFinalInterval"};

int mvmc_power_lanczos_parse_ascii_nonnegative_int(
    const char *token, int *value) {
  unsigned int parsed = 0;
  size_t index;
  if (token == NULL || value == NULL || token[0] == '\0') return 1;
  if (token[0] == '0') {
    if (token[1] != '\0') return 1;
    *value = 0;
    return 0;
  }
  if (token[0] < '1' || token[0] > '9') return 1;
  for (index = 0; token[index] != '\0'; ++index) {
    const unsigned int digit = (unsigned int)(token[index] - '0');
    if (token[index] < '0' || token[index] > '9' ||
        parsed > ((unsigned int)INT_MAX - digit) / 10u) {
      return 1;
    }
    parsed = parsed * 10u + digit;
  }
  *value = (int)parsed;
  return 0;
}

static MVMCPowerLanczosRuntimeStatus Commit(
    MVMCPowerLanczosRuntimeStatus status, MVMCPowerLanczosRoute route,
    int invalid_chain_control, MVMCPowerLanczosRuntimeResult *result) {
  MVMCPowerLanczosRuntimeResult candidate;
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = status;
  candidate.route = route;
  candidate.invalid_chain_control = invalid_chain_control;
  *result = candidate;
  return status;
}

MVMCPowerLanczosRuntimeStatus mvmc_power_lanczos_runtime_validate(
    const MVMCPowerLanczosRuntimeOptions *options,
    MVMCPowerLanczosRuntimeResult *result) {
  int control;
  int has_chain_control = 0;
  if (result != NULL) {
    memset(result, 0, sizeof(*result));
    result->status = MVMC_POWER_LANCZOS_RUNTIME_INVALID_ARGUMENT;
    result->invalid_chain_control = -1;
  }
  if (options == NULL || result == NULL || options->lanczos_mode < 0) {
    return MVMC_POWER_LANCZOS_RUNTIME_INVALID_ARGUMENT;
  }
  if (options->estimator_mode != 0 && options->estimator_mode != 1) {
    return Commit(MVMC_POWER_LANCZOS_RUNTIME_INVALID_ESTIMATOR_MODE,
                  MVMC_POWER_LANCZOS_ROUTE_DISABLED, -1, result);
  }
  for (control = 0; control < MVMC_POWER_LANCZOS_CHAIN_CONTROL_COUNT;
       ++control) {
    if (options->chain_controls[control] < 0) {
      return Commit(MVMC_POWER_LANCZOS_RUNTIME_INVALID_CHAIN_CONTROL,
                    MVMC_POWER_LANCZOS_ROUTE_DISABLED, control, result);
    }
    has_chain_control |= options->chain_controls[control] != 0;
  }
  if (options->guide_mode != 0) {
    return Commit(MVMC_POWER_LANCZOS_RUNTIME_INVALID_GUIDE_MODE,
                  MVMC_POWER_LANCZOS_ROUTE_DISABLED, -1, result);
  }
  if (options->stat_mode != 0) {
    return Commit(MVMC_POWER_LANCZOS_RUNTIME_INVALID_STAT_MODE,
                  MVMC_POWER_LANCZOS_ROUTE_DISABLED, -1, result);
  }
  if (options->lanczos_mode == 0) {
    if (options->estimator_mode != 0 || has_chain_control) {
      return Commit(MVMC_POWER_LANCZOS_RUNTIME_DISABLED_CONTROL_CONFLICT,
                    MVMC_POWER_LANCZOS_ROUTE_DISABLED, -1, result);
    }
    return Commit(MVMC_POWER_LANCZOS_RUNTIME_OK,
                  MVMC_POWER_LANCZOS_ROUTE_DISABLED, -1, result);
  }
  if (options->estimator_mode == 1) {
    if (has_chain_control) {
      return Commit(MVMC_POWER_LANCZOS_RUNTIME_LEGACY_CONTROL_CONFLICT,
                    MVMC_POWER_LANCZOS_ROUTE_DISABLED, -1, result);
    }
    return Commit(MVMC_POWER_LANCZOS_RUNTIME_OK,
                  MVMC_POWER_LANCZOS_ROUTE_LEGACY, -1, result);
  }
  return Commit(MVMC_POWER_LANCZOS_RUNTIME_OK,
                MVMC_POWER_LANCZOS_ROUTE_CORRECTED, -1, result);
}

const char *mvmc_power_lanczos_chain_control_name(int control) {
  if (control < 0 || control >= MVMC_POWER_LANCZOS_CHAIN_CONTROL_COUNT) {
    return "unknown power-Lanczos chain control";
  }
  return ChainControlNames[control];
}

const char *mvmc_power_lanczos_runtime_error(
    const MVMCPowerLanczosRuntimeResult *result) {
  if (result == NULL || !result->valid) {
    return "P6 INPUT REJECTED: INVALID_RUNTIME_OPTIONS";
  }
  switch (result->status) {
    case MVMC_POWER_LANCZOS_RUNTIME_OK:
      return "";
    case MVMC_POWER_LANCZOS_RUNTIME_INVALID_ARGUMENT:
      return "P6 INPUT REJECTED: INVALID_RUNTIME_OPTIONS";
    case MVMC_POWER_LANCZOS_RUNTIME_INVALID_ESTIMATOR_MODE:
      return "P6 INPUT REJECTED: NLanczosEstimatorMode must be the integer 0 (corrected) or 1 (explicit legacy).";
    case MVMC_POWER_LANCZOS_RUNTIME_INVALID_CHAIN_CONTROL:
      return "P6 INPUT REJECTED: chain controls must be non-negative integers.";
    case MVMC_POWER_LANCZOS_RUNTIME_INVALID_GUIDE_MODE:
      return "P6 INPUT REJECTED: NLanczosGuideMode only supports integer 0 in P6.";
    case MVMC_POWER_LANCZOS_RUNTIME_INVALID_STAT_MODE:
      return "P6 INPUT REJECTED: NLanczosStatMode only supports integer 0 in P6.";
    case MVMC_POWER_LANCZOS_RUNTIME_DISABLED_CONTROL_CONFLICT:
      return "P6 INPUT REJECTED: NLanczosMode=0 requires every P6 estimator and chain control at default zero.";
    case MVMC_POWER_LANCZOS_RUNTIME_LEGACY_CONTROL_CONFLICT:
      return "P6 INPUT REJECTED: explicit legacy estimator mode requires all coefficient and final chain controls to be zero.";
  }
  return "P6 INPUT REJECTED: UNKNOWN_RUNTIME_OPTIONS_STATUS";
}
