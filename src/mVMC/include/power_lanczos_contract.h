#ifndef MVMC_POWER_LANCZOS_CONTRACT_H
#define MVMC_POWER_LANCZOS_CONTRACT_H

#include <stddef.h>

#define MVMC_POWER_LANCZOS_PRODUCTION_CONTRACT_ID \
  "power-lanczos-p6a-production-contract-v2"

typedef struct {
  int order;
  int lanczosMode;
  int vmcCalMode;
  int orbitalGeneral;
  int twoSz;
  int nProjBF;
  int flagRBM;
  int reweight;
  int exUpdatePath;
  int nPairHopping;
  int nExchangeCoupling;
  int nInterAll;
  int nNBodyInterAll;
  int nNBodyG;
  int nSpinFlipTransfer;
  int nLocSpn;
  int nsite;
  int ne;
  int nTransfer;
  int nQPFull;
} MVMCPowerLanczosContract;

typedef enum {
  MVMC_POWER_LANCZOS_MODEL_INVALID = 0,
  MVMC_POWER_LANCZOS_MODEL_ELECTRONIC_VK,
  MVMC_POWER_LANCZOS_MODEL_LOCAL_SPIN_EXCHANGE
} MVMCPowerLanczosModelClass;

typedef enum {
  MVMC_POWER_LANCZOS_CONTRACT_OK = 0,
  MVMC_POWER_LANCZOS_CONTRACT_INVALID_ARGUMENT,
  MVMC_POWER_LANCZOS_CONTRACT_INVALID_ORDER,
  MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_LANCZOS_MODE_1,
  MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_VMC_CAL_MODE_1,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TWO_SZ,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_BACKFLOW,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_RBM,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_REWEIGHT,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_UPDATE_PATH,
  MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_PAIR_HOPPING,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_EXCHANGE,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_INTER_ALL,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_G,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER,
  MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN,
  MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TRANSFER,
  MVMC_POWER_LANCZOS_CONTRACT_INVALID_COUNT,
  MVMC_POWER_LANCZOS_CONTRACT_SIZE_OVERFLOW
} MVMCPowerLanczosContractStatus;

typedef struct {
  size_t nsite2;
  size_t nsize;
  size_t one_body_matrix_elements;
  size_t projected_one_body_matrix_elements;
  size_t coefficient_basis_dimension;
  size_t matrix_elements;
  size_t anchor_count;
} MVMCPowerLanczosDerivedSizes;

MVMCPowerLanczosModelClass mvmc_power_lanczos_classify_model(
    const MVMCPowerLanczosContract *contract);

MVMCPowerLanczosContractStatus mvmc_power_lanczos_compute_derived_sizes(
    const MVMCPowerLanczosContract *contract,
    MVMCPowerLanczosDerivedSizes *sizes);

MVMCPowerLanczosContractStatus
mvmc_power_lanczos_validate_common_contract(
    const MVMCPowerLanczosContract *contract);

MVMCPowerLanczosContractStatus
mvmc_power_lanczos_validate_production_contract(
    const MVMCPowerLanczosContract *contract);

const char *mvmc_power_lanczos_contract_error(
    MVMCPowerLanczosContractStatus status);

const char *mvmc_power_lanczos_contract_id(void);

#endif
