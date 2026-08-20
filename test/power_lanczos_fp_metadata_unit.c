#include <stdio.h>
#include <string.h>

#include "power_lanczos_contract.h"
#include "power_lanczos_fp_metadata.h"

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "PowerLanczosBuildMetadata_Unit FAIL: ");                \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                    \
      failures++;                                                               \
    }                                                                           \
  } while (0)

int main(void) {
  const MVMCPowerLanczosBuildMetadata *metadata =
      mvmc_power_lanczos_build_metadata();

  CHECK(metadata != NULL, "metadata pointer is NULL");
  if (metadata != NULL) {
    CHECK(strcmp(metadata->contract_id,
                 MVMC_POWER_LANCZOS_PRODUCTION_CONTRACT_ID) == 0,
          "contract id mismatch");
    CHECK(strcmp(metadata->strict_fp_option_policy,
                 "strict-fp-options:fno-fast-math,ffp-contract-off") == 0,
          "strict FP policy mismatch");
    CHECK(strcmp(metadata->complex_fma_contract_policy,
                 "complex-fma-contract:c99-complex,record-__FP_FAST_FMA,ffp-contract-off-requested") == 0,
          "complex/FMA policy mismatch");
    CHECK(metadata->fast_math_enabled == 0,
          "fast-math must be disabled for power-Lanczos production core");
    CHECK(metadata->fp_contract_off_requested != 0,
          "ffp-contract=off request was not recorded");
  }
  CHECK(mvmc_power_lanczos_strict_fp_enabled() != 0,
        "strict FP metadata gate did not pass");

  if (failures != 0) {
    fprintf(stderr, "PowerLanczosBuildMetadata_Unit: %d failure(s)\n",
            failures);
    return 1;
  }
  printf("PowerLanczosBuildMetadata_Unit: PASS\n");
  return 0;
}
