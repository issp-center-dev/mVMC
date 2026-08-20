#ifndef MVMC_POWER_LANCZOS_FP_METADATA_H
#define MVMC_POWER_LANCZOS_FP_METADATA_H

typedef struct {
  const char *contract_id;
  const char *strict_fp_option_policy;
  const char *complex_fma_contract_policy;
  int fast_math_enabled;
  int fp_contract_off_requested;
  int fp_fast_fma_defined;
  int fp_fast_fmaf_defined;
  int fp_fast_fmal_defined;
} MVMCPowerLanczosBuildMetadata;

const MVMCPowerLanczosBuildMetadata *
mvmc_power_lanczos_build_metadata(void);

int mvmc_power_lanczos_strict_fp_enabled(void);

#endif
