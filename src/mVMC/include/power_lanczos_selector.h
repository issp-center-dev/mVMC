#ifndef MVMC_POWER_LANCZOS_SELECTOR_H
#define MVMC_POWER_LANCZOS_SELECTOR_H

#include <stddef.h>

#define MVMC_POWER_LANCZOS_INPUT_REJECTED_EXIT_CODE 20
#define MVMC_POWER_LANCZOS_CHAIN_CONTROL_COUNT 6

typedef enum {
  MVMC_POWER_LANCZOS_ROUTE_DISABLED = 0,
  MVMC_POWER_LANCZOS_ROUTE_CORRECTED,
  MVMC_POWER_LANCZOS_ROUTE_LEGACY
} MVMCPowerLanczosRoute;

typedef enum {
  MVMC_POWER_LANCZOS_CONTROL_COEFF_WARMUP = 0,
  MVMC_POWER_LANCZOS_CONTROL_COEFF_SAMPLE,
  MVMC_POWER_LANCZOS_CONTROL_COEFF_INTERVAL,
  MVMC_POWER_LANCZOS_CONTROL_FINAL_WARMUP,
  MVMC_POWER_LANCZOS_CONTROL_FINAL_SAMPLE,
  MVMC_POWER_LANCZOS_CONTROL_FINAL_INTERVAL
} MVMCPowerLanczosChainControl;

typedef struct {
  int lanczos_mode;
  int estimator_mode;
  int chain_controls[MVMC_POWER_LANCZOS_CHAIN_CONTROL_COUNT];
  int guide_mode;
  int stat_mode;
} MVMCPowerLanczosRuntimeOptions;

typedef enum {
  MVMC_POWER_LANCZOS_RUNTIME_OK = 0,
  MVMC_POWER_LANCZOS_RUNTIME_INVALID_ARGUMENT,
  MVMC_POWER_LANCZOS_RUNTIME_INVALID_ESTIMATOR_MODE,
  MVMC_POWER_LANCZOS_RUNTIME_INVALID_CHAIN_CONTROL,
  MVMC_POWER_LANCZOS_RUNTIME_INVALID_GUIDE_MODE,
  MVMC_POWER_LANCZOS_RUNTIME_INVALID_STAT_MODE,
  MVMC_POWER_LANCZOS_RUNTIME_CORRECTED_OBSERVABLE_UNSUPPORTED,
  MVMC_POWER_LANCZOS_RUNTIME_DISABLED_CONTROL_CONFLICT,
  MVMC_POWER_LANCZOS_RUNTIME_LEGACY_CONTROL_CONFLICT
} MVMCPowerLanczosRuntimeStatus;

typedef struct {
  int valid;
  MVMCPowerLanczosRuntimeStatus status;
  MVMCPowerLanczosRoute route;
  int invalid_chain_control;
} MVMCPowerLanczosRuntimeResult;

/*
 * Frozen P6 integer grammar: exactly "0" or [1-9][0-9]*, with a value no
 * larger than INT_MAX.  Whitespace belongs to the caller's line grammar and
 * is not accepted inside this token API.
 */
int mvmc_power_lanczos_parse_ascii_nonnegative_int(
    const char *token, int *value);

/* Pure, allocation-free selector validation. */
MVMCPowerLanczosRuntimeStatus mvmc_power_lanczos_runtime_validate(
    const MVMCPowerLanczosRuntimeOptions *options,
    MVMCPowerLanczosRuntimeResult *result);

const char *mvmc_power_lanczos_chain_control_name(int control);

const char *mvmc_power_lanczos_runtime_error(
    const MVMCPowerLanczosRuntimeResult *result);

#endif
