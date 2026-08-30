#include "power_lanczos_selector.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition, message)                                             \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "power_lanczos_runtime_unit: %s (line %d)\n",          \
              (message), __LINE__);                                           \
      exit(EXIT_FAILURE);                                                     \
    }                                                                         \
  } while (0)

static MVMCPowerLanczosRuntimeOptions DefaultOptions(void) {
  MVMCPowerLanczosRuntimeOptions options;
  memset(&options, 0, sizeof(options));
  return options;
}

static void TestIntegerGrammar(void) {
  static const char *const rejected[] = {
      "",       "-1",  "+1", "00", "01", "1.0", "1e0",
      "0x1",    "NaN", "Inf", "１", "1 #comment", "1 extra",
      "2147483648"};
  int value = -1;
  size_t index;
  CHECK(mvmc_power_lanczos_parse_ascii_nonnegative_int("0", &value) == 0 &&
            value == 0,
        "zero grammar");
  CHECK(mvmc_power_lanczos_parse_ascii_nonnegative_int("1", &value) == 0 &&
            value == 1,
        "one grammar");
  CHECK(mvmc_power_lanczos_parse_ascii_nonnegative_int("2147483647", &value) ==
                0 &&
            value == INT_MAX,
        "INT_MAX grammar");
  for (index = 0; index < sizeof(rejected) / sizeof(rejected[0]); ++index) {
    value = 17;
    CHECK(mvmc_power_lanczos_parse_ascii_nonnegative_int(rejected[index],
                                                         &value) != 0,
          rejected[index]);
    CHECK(value == 17, "failed integer parse must be atomic");
  }
  CHECK(mvmc_power_lanczos_parse_ascii_nonnegative_int(NULL, &value) != 0,
        "NULL token");
  CHECK(mvmc_power_lanczos_parse_ascii_nonnegative_int("1", NULL) != 0,
        "NULL output");
}

static void TestRoutes(void) {
  MVMCPowerLanczosRuntimeOptions options = DefaultOptions();
  MVMCPowerLanczosRuntimeResult result;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_OK &&
            result.route == MVMC_POWER_LANCZOS_ROUTE_DISABLED,
        "disabled route");
  options.lanczos_mode = 1;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_OK &&
            result.route == MVMC_POWER_LANCZOS_ROUTE_CORRECTED,
        "corrected default route");
  options.estimator_mode = 1;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_OK &&
            result.route == MVMC_POWER_LANCZOS_ROUTE_LEGACY,
        "explicit legacy route");
  options.chain_controls[MVMC_POWER_LANCZOS_CONTROL_COEFF_SAMPLE] = 16;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_LEGACY_CONTROL_CONFLICT,
        "legacy chain control conflict");
  options = DefaultOptions();
  options.lanczos_mode = 1;
  options.chain_controls[MVMC_POWER_LANCZOS_CONTROL_FINAL_SAMPLE] = 64;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_OK &&
            result.route == MVMC_POWER_LANCZOS_ROUTE_CORRECTED,
        "corrected positive control route");
  options = DefaultOptions();
  options.lanczos_mode = 2;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
                MVMC_POWER_LANCZOS_RUNTIME_CORRECTED_OBSERVABLE_UNSUPPORTED &&
            result.route == MVMC_POWER_LANCZOS_ROUTE_DISABLED &&
            strstr(mvmc_power_lanczos_runtime_error(&result),
                   "additional observable production enable is out of scope") !=
                NULL,
        "corrected additional observable is explicitly out of scope");
  options.estimator_mode = 1;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
                MVMC_POWER_LANCZOS_RUNTIME_OK &&
            result.route == MVMC_POWER_LANCZOS_ROUTE_LEGACY,
        "legacy augmented route remains available");
}

static void TestInvalidOptions(void) {
  MVMCPowerLanczosRuntimeOptions options = DefaultOptions();
  MVMCPowerLanczosRuntimeResult result;
  options.estimator_mode = 2;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_INVALID_ESTIMATOR_MODE,
        "invalid estimator");
  options = DefaultOptions();
  options.chain_controls[MVMC_POWER_LANCZOS_CONTROL_FINAL_INTERVAL] = -1;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
                MVMC_POWER_LANCZOS_RUNTIME_INVALID_CHAIN_CONTROL &&
            result.invalid_chain_control ==
                MVMC_POWER_LANCZOS_CONTROL_FINAL_INTERVAL,
        "negative chain control");
  CHECK(strcmp(mvmc_power_lanczos_chain_control_name(
                   result.invalid_chain_control),
               "NLanczosFinalInterval") == 0,
        "negative control identity");
  options = DefaultOptions();
  options.lanczos_mode = 1;
  options.chain_controls[MVMC_POWER_LANCZOS_CONTROL_COEFF_SAMPLE] = 33;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
                MVMC_POWER_LANCZOS_RUNTIME_INVALID_CHAIN_CONTROL &&
            result.invalid_chain_control ==
                MVMC_POWER_LANCZOS_CONTROL_COEFF_SAMPLE,
        "corrected coefficient sample block contract");
  options = DefaultOptions();
  options.lanczos_mode = 1;
  options.chain_controls[MVMC_POWER_LANCZOS_CONTROL_FINAL_SAMPLE] = 16;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
                MVMC_POWER_LANCZOS_RUNTIME_INVALID_CHAIN_CONTROL &&
            result.invalid_chain_control ==
                MVMC_POWER_LANCZOS_CONTROL_FINAL_SAMPLE,
        "corrected final sample minimum");
  options = DefaultOptions();
  options.guide_mode = 1;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_INVALID_GUIDE_MODE,
        "guide hook v1 domain");
  options = DefaultOptions();
  options.stat_mode = 1;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_INVALID_STAT_MODE,
        "stat hook v1 domain");
  options = DefaultOptions();
  options.estimator_mode = 1;
  CHECK(mvmc_power_lanczos_runtime_validate(&options, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_DISABLED_CONTROL_CONFLICT,
        "disabled estimator conflict");
  CHECK(mvmc_power_lanczos_runtime_validate(NULL, &result) ==
            MVMC_POWER_LANCZOS_RUNTIME_INVALID_ARGUMENT,
        "NULL options");
  CHECK(mvmc_power_lanczos_runtime_validate(&options, NULL) ==
            MVMC_POWER_LANCZOS_RUNTIME_INVALID_ARGUMENT,
        "NULL result");
}

int main(void) {
  TestIntegerGrammar();
  TestRoutes();
  TestInvalidOptions();
  puts("power_lanczos_runtime_unit: PASS");
  return EXIT_SUCCESS;
}
