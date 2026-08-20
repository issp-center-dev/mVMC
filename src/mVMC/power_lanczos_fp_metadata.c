#include "power_lanczos_contract.h"
#include "power_lanczos_fp_metadata.h"

#define MVMC_POWER_LANCZOS_STRICT_FP_OPTION_POLICY \
  "strict-fp-options:fno-fast-math,ffp-contract-off"
#define MVMC_POWER_LANCZOS_COMPLEX_FMA_CONTRACT_POLICY \
  "complex-fma-contract:c99-complex,record-__FP_FAST_FMA,ffp-contract-off-requested"

const MVMCPowerLanczosBuildMetadata *
mvmc_power_lanczos_build_metadata(void) {
  static const MVMCPowerLanczosBuildMetadata metadata = {
      MVMC_POWER_LANCZOS_PRODUCTION_CONTRACT_ID,
      MVMC_POWER_LANCZOS_STRICT_FP_OPTION_POLICY,
      MVMC_POWER_LANCZOS_COMPLEX_FMA_CONTRACT_POLICY,
#ifdef __FAST_MATH__
      1,
#else
      0,
#endif
#ifdef MVMC_POWER_LANCZOS_FP_CONTRACT_OFF_REQUESTED
      1,
#else
      0,
#endif
#ifdef __FP_FAST_FMA
      1,
#else
      0,
#endif
#ifdef __FP_FAST_FMAF
      1,
#else
      0,
#endif
#ifdef __FP_FAST_FMAL
      1
#else
      0
#endif
  };
  return &metadata;
}

int mvmc_power_lanczos_strict_fp_enabled(void) {
  const MVMCPowerLanczosBuildMetadata *metadata =
      mvmc_power_lanczos_build_metadata();
  return metadata != 0 &&
         metadata->fast_math_enabled == 0 &&
         metadata->fp_contract_off_requested != 0;
}
