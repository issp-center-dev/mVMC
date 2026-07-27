#include <stddef.h>

#include "lanczos2_contract.h"

Lanczos2ContractStatus ValidateLanczos2Contract(
    const Lanczos2Contract *contract) {
  if (contract == NULL || (contract->step != 1 && contract->step != 2)) {
    return LANCZOS2_CONTRACT_INVALID_STEP;
  }
  if (contract->step == 1) return LANCZOS2_CONTRACT_OK;

  if (contract->lanczosMode != 1) {
    return LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1;
  }
  if (contract->vmcCalMode != 1) {
    return LANCZOS2_CONTRACT_REQUIRES_VMC_CAL_MODE_1;
  }
  if (contract->orbitalGeneral != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL;
  }
  if (contract->nProjBF != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_BACKFLOW;
  }
  if (contract->flagRBM != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_RBM;
  }
  if (contract->exUpdatePath != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH;
  }
  if (contract->nPairHopping != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_PAIR_HOPPING;
  }
  if (contract->nExchangeCoupling != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_EXCHANGE;
  }
  if (contract->nInterAll != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_INTER_ALL;
  }
  if (contract->nNBodyInterAll != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL;
  }
  if (contract->nNBodyG != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_G;
  }
  if (contract->nSpinFlipTransfer != 0) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER;
  }
  return LANCZOS2_CONTRACT_OK;
}

const char *Lanczos2ContractError(Lanczos2ContractStatus status) {
  switch (status) {
    case LANCZOS2_CONTRACT_OK:
      return "";
    case LANCZOS2_CONTRACT_INVALID_STEP:
      return "NLanczosStep must be 1 or 2";
    case LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1:
      return "NLanczosStep=2 requires NLanczosMode=1";
    case LANCZOS2_CONTRACT_REQUIRES_VMC_CAL_MODE_1:
      return "NLanczosStep=2 requires NVMCCalMode=1";
    case LANCZOS2_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL:
      return "NLanczosStep=2 does not support orbital-general (FSZ) mode";
    case LANCZOS2_CONTRACT_UNSUPPORTED_BACKFLOW:
      return "NLanczosStep=2 does not support BackFlow";
    case LANCZOS2_CONTRACT_UNSUPPORTED_RBM:
      return "NLanczosStep=2 does not support RBM";
    case LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH:
      return "NLanczosStep=2 requires NExUpdatePath=0";
    case LANCZOS2_CONTRACT_UNSUPPORTED_PAIR_HOPPING:
      return "NLanczosStep=2 does not support PairHopping";
    case LANCZOS2_CONTRACT_UNSUPPORTED_EXCHANGE:
      return "NLanczosStep=2 does not support Exchange";
    case LANCZOS2_CONTRACT_UNSUPPORTED_INTER_ALL:
      return "NLanczosStep=2 does not support InterAll";
    case LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL:
      return "NLanczosStep=2 does not support NBodyInterAll";
    case LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_G:
      return "NLanczosStep=2 does not support NBodyG";
    case LANCZOS2_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER:
      return "NLanczosStep=2 supports only spin-conserving Transfer terms";
  }
  return "unknown NLanczosStep contract error";
}
