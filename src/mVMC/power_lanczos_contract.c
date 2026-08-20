#include <limits.h>
#include <stddef.h>

#include "power_lanczos_contract.h"

static int CountFieldsAreNonnegative(
    const MVMCPowerLanczosContract *contract) {
  return contract->lanczosMode >= 0 &&
         contract->vmcCalMode >= 0 &&
         contract->orbitalGeneral >= 0 &&
         contract->nProjBF >= 0 &&
         contract->flagRBM >= 0 &&
         contract->reweight >= 0 &&
         contract->nPairHopping >= 0 &&
         contract->nExchangeCoupling >= 0 &&
         contract->nInterAll >= 0 &&
         contract->nNBodyInterAll >= 0 &&
         contract->nNBodyG >= 0 &&
         contract->nSpinFlipTransfer >= 0 &&
         contract->nLocSpn >= 0 &&
         contract->nTransfer >= 0;
}

static int MulSize(size_t lhs, size_t rhs, size_t *out) {
  if (out == NULL) return 0;
  if (lhs != 0 && rhs > ((size_t)-1) / lhs) return 1;
  *out = lhs * rhs;
  return 0;
}

MVMCPowerLanczosModelClass mvmc_power_lanczos_classify_model(
    const MVMCPowerLanczosContract *contract) {
  if (contract == NULL) return MVMC_POWER_LANCZOS_MODEL_INVALID;
  if (contract->exUpdatePath == 0) {
    return MVMC_POWER_LANCZOS_MODEL_ELECTRONIC_VK;
  }
  if (contract->exUpdatePath == 2 &&
      contract->nsite > 0 &&
      contract->ne > 0 &&
      contract->ne <= INT_MAX / 2 &&
      contract->nLocSpn == contract->nsite &&
      contract->nLocSpn == 2 * contract->ne) {
    return MVMC_POWER_LANCZOS_MODEL_LOCAL_SPIN_EXCHANGE;
  }
  return MVMC_POWER_LANCZOS_MODEL_INVALID;
}

MVMCPowerLanczosContractStatus mvmc_power_lanczos_compute_derived_sizes(
    const MVMCPowerLanczosContract *contract,
    MVMCPowerLanczosDerivedSizes *sizes) {
  MVMCPowerLanczosDerivedSizes localSizes;
  size_t nsize;
  size_t oneBodyElements;
  size_t basisDimension;

  if (contract == NULL) {
    return MVMC_POWER_LANCZOS_CONTRACT_INVALID_ARGUMENT;
  }
  if (contract->order != 1 && contract->order != 2) {
    return MVMC_POWER_LANCZOS_CONTRACT_INVALID_ORDER;
  }
  if (!CountFieldsAreNonnegative(contract) ||
      contract->nsite <= 0 ||
      contract->ne <= 0 ||
      contract->nQPFull <= 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_INVALID_COUNT;
  }
  if (contract->nsite > INT_MAX / 2 ||
      contract->ne > INT_MAX / 2) {
    return MVMC_POWER_LANCZOS_CONTRACT_SIZE_OVERFLOW;
  }

  localSizes.nsite2 = (size_t)contract->nsite * 2u;
  localSizes.nsize = (size_t)contract->ne * 2u;
  nsize = localSizes.nsize;
  if (MulSize(nsize, nsize, &oneBodyElements) != 0 ||
      MulSize((size_t)contract->nQPFull, oneBodyElements,
              &localSizes.projected_one_body_matrix_elements) != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_SIZE_OVERFLOW;
  }
  localSizes.one_body_matrix_elements = oneBodyElements;

  basisDimension = (size_t)contract->order + 1u;
  localSizes.coefficient_basis_dimension = basisDimension;
  if (MulSize(basisDimension, basisDimension,
              &localSizes.matrix_elements) != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_SIZE_OVERFLOW;
  }
  localSizes.anchor_count = (size_t)contract->order + 2u;

  if (sizes != NULL) *sizes = localSizes;
  return MVMC_POWER_LANCZOS_CONTRACT_OK;
}

MVMCPowerLanczosContractStatus
mvmc_power_lanczos_validate_common_contract(
    const MVMCPowerLanczosContract *contract) {
  return mvmc_power_lanczos_compute_derived_sizes(contract, NULL);
}

MVMCPowerLanczosContractStatus
mvmc_power_lanczos_validate_production_contract(
    const MVMCPowerLanczosContract *contract) {
  MVMCPowerLanczosContractStatus status;
  MVMCPowerLanczosModelClass model;

  status = mvmc_power_lanczos_validate_common_contract(contract);
  if (status != MVMC_POWER_LANCZOS_CONTRACT_OK) return status;

  if (contract->lanczosMode != 1) {
    return MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_LANCZOS_MODE_1;
  }
  if (contract->vmcCalMode != 1) {
    return MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_VMC_CAL_MODE_1;
  }
  if (contract->orbitalGeneral != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL;
  }
  if (contract->twoSz != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TWO_SZ;
  }
  if (contract->nProjBF != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_BACKFLOW;
  }
  if (contract->flagRBM != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_RBM;
  }
  if (contract->reweight != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_REWEIGHT;
  }
  if (contract->exUpdatePath != 0 && contract->exUpdatePath != 2) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_UPDATE_PATH;
  }

  model = mvmc_power_lanczos_classify_model(contract);
  if (contract->nPairHopping != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_PAIR_HOPPING;
  }
  if (contract->nInterAll != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_INTER_ALL;
  }
  if (contract->nNBodyInterAll != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL;
  }
  if (contract->nNBodyG != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_G;
  }
  if (contract->nSpinFlipTransfer != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER;
  }
  if (model == MVMC_POWER_LANCZOS_MODEL_ELECTRONIC_VK &&
      contract->nExchangeCoupling != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_EXCHANGE;
  }
  if (model == MVMC_POWER_LANCZOS_MODEL_ELECTRONIC_VK &&
      contract->nTransfer <= 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER;
  }
  if (contract->exUpdatePath == 2 &&
      model != MVMC_POWER_LANCZOS_MODEL_LOCAL_SPIN_EXCHANGE) {
    return MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN;
  }
  if (model == MVMC_POWER_LANCZOS_MODEL_LOCAL_SPIN_EXCHANGE &&
      contract->nTransfer != 0) {
    return MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TRANSFER;
  }
  return MVMC_POWER_LANCZOS_CONTRACT_OK;
}

const char *mvmc_power_lanczos_contract_error(
    MVMCPowerLanczosContractStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_CONTRACT_OK:
      return "";
    case MVMC_POWER_LANCZOS_CONTRACT_INVALID_ARGUMENT:
      return "power-Lanczos contract is missing";
    case MVMC_POWER_LANCZOS_CONTRACT_INVALID_ORDER:
      return "power-Lanczos order must be 1 or 2";
    case MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_LANCZOS_MODE_1:
      return "power-Lanczos corrected production contract requires NLanczosMode=1";
    case MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_VMC_CAL_MODE_1:
      return "power-Lanczos corrected production contract requires NVMCCalMode=1";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL:
      return "power-Lanczos corrected production contract does not support orbital-general (FSZ) mode";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TWO_SZ:
      return "power-Lanczos corrected production contract requires 2Sz=0";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_BACKFLOW:
      return "power-Lanczos corrected production contract does not support BackFlow";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_RBM:
      return "power-Lanczos corrected production contract does not support RBM";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_REWEIGHT:
      return "power-Lanczos corrected production contract does not support reweight=1";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_UPDATE_PATH:
      return "power-Lanczos corrected production contract supports NExUpdatePath=0, or NExUpdatePath=2 for pure-spin mode";
    case MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER:
      return "power-Lanczos electronic contract requires at least one spin-conserving Transfer term";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_PAIR_HOPPING:
      return "power-Lanczos corrected production contract does not support PairHopping";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_EXCHANGE:
      return "power-Lanczos corrected production contract does not support electronic Exchange";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_INTER_ALL:
      return "power-Lanczos corrected production contract does not support InterAll";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL:
      return "power-Lanczos corrected production contract does not support NBodyInterAll";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_G:
      return "power-Lanczos corrected production contract does not support NBodyG";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER:
      return "power-Lanczos corrected production contract supports only spin-conserving Transfer terms";
    case MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN:
      return "power-Lanczos NExUpdatePath=2 contract requires NLocalSpin=Nsite=2*Ne";
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TRANSFER:
      return "power-Lanczos pure-spin contract requires NTransfer=0";
    case MVMC_POWER_LANCZOS_CONTRACT_INVALID_COUNT:
      return "power-Lanczos corrected production contract requires non-negative counts and positive Nsite, Ne, and NQPFull";
    case MVMC_POWER_LANCZOS_CONTRACT_SIZE_OVERFLOW:
      return "power-Lanczos corrected production contract size arithmetic overflows the supported range";
  }
  return "unknown power-Lanczos contract error";
}

const char *mvmc_power_lanczos_contract_id(void) {
  return MVMC_POWER_LANCZOS_PRODUCTION_CONTRACT_ID;
}
