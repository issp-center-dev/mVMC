#include <limits.h>
#include <stddef.h>

#include "lanczos2_contract.h"
#include "power_lanczos_contract.h"

static MVMCPowerLanczosContract ToPowerContract(
    const Lanczos2Contract *contract) {
  MVMCPowerLanczosContract powerContract;
  powerContract.order = contract->step;
  powerContract.lanczosMode = contract->lanczosMode;
  powerContract.vmcCalMode = contract->vmcCalMode;
  powerContract.orbitalGeneral = contract->orbitalGeneral;
  powerContract.twoSz = contract->twoSz;
  powerContract.nProjBF = contract->nProjBF;
  powerContract.flagRBM = contract->flagRBM;
  powerContract.reweight = contract->reweight;
  powerContract.exUpdatePath = contract->exUpdatePath;
  powerContract.nPairHopping = contract->nPairHopping;
  powerContract.nExchangeCoupling = contract->nExchangeCoupling;
  powerContract.nInterAll = contract->nInterAll;
  powerContract.nNBodyInterAll = contract->nNBodyInterAll;
  powerContract.nNBodyG = contract->nNBodyG;
  powerContract.nSpinFlipTransfer = contract->nSpinFlipTransfer;
  powerContract.nLocSpn = contract->nLocSpn;
  powerContract.nsite = contract->nsite;
  powerContract.ne = contract->ne;
  powerContract.nTransfer = contract->nTransfer;
  powerContract.nQPFull = contract->nQPFull;
  return powerContract;
}

static Lanczos2ContractStatus FromPowerStatus(
    MVMCPowerLanczosContractStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_CONTRACT_OK:
      return LANCZOS2_CONTRACT_OK;
    case MVMC_POWER_LANCZOS_CONTRACT_INVALID_ARGUMENT:
    case MVMC_POWER_LANCZOS_CONTRACT_INVALID_ORDER:
      return LANCZOS2_CONTRACT_INVALID_STEP;
    case MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_LANCZOS_MODE_1:
      return LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1;
    case MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_VMC_CAL_MODE_1:
      return LANCZOS2_CONTRACT_REQUIRES_VMC_CAL_MODE_1;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL:
      return LANCZOS2_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TWO_SZ:
      return LANCZOS2_CONTRACT_UNSUPPORTED_TWO_SZ;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_BACKFLOW:
      return LANCZOS2_CONTRACT_UNSUPPORTED_BACKFLOW;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_RBM:
      return LANCZOS2_CONTRACT_UNSUPPORTED_RBM;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_REWEIGHT:
      return LANCZOS2_CONTRACT_UNSUPPORTED_REWEIGHT;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_UPDATE_PATH:
      return LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH;
    case MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER:
      return LANCZOS2_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TRANSFER:
      return LANCZOS2_CONTRACT_UNSUPPORTED_TRANSFER;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_PAIR_HOPPING:
      return LANCZOS2_CONTRACT_UNSUPPORTED_PAIR_HOPPING;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_EXCHANGE:
      return LANCZOS2_CONTRACT_UNSUPPORTED_EXCHANGE;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_INTER_ALL:
      return LANCZOS2_CONTRACT_UNSUPPORTED_INTER_ALL;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL:
      return LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_G:
      return LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_G;
    case MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER:
      return LANCZOS2_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER;
    case MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN:
      return LANCZOS2_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN;
    case MVMC_POWER_LANCZOS_CONTRACT_INVALID_COUNT:
      return LANCZOS2_CONTRACT_INVALID_COUNT;
    case MVMC_POWER_LANCZOS_CONTRACT_SIZE_OVERFLOW:
      return LANCZOS2_CONTRACT_SIZE_OVERFLOW;
  }
  return LANCZOS2_CONTRACT_INVALID_STEP;
}

Lanczos2ModelClass ClassifyLanczos2Model(
    const Lanczos2Contract *contract) {
  MVMCPowerLanczosContract powerContract;
  if (contract == NULL) return LANCZOS2_MODEL_INVALID;
  powerContract = ToPowerContract(contract);
  switch (mvmc_power_lanczos_classify_model(&powerContract)) {
    case MVMC_POWER_LANCZOS_MODEL_ELECTRONIC_VK:
      return LANCZOS2_MODEL_ELECTRONIC_VK;
    case MVMC_POWER_LANCZOS_MODEL_LOCAL_SPIN_EXCHANGE:
      return LANCZOS2_MODEL_LOCAL_SPIN_EXCHANGE;
    case MVMC_POWER_LANCZOS_MODEL_INVALID:
      return LANCZOS2_MODEL_INVALID;
  }
  return LANCZOS2_MODEL_INVALID;
}

Lanczos2ContractStatus ValidateLanczos2Contract(
    const Lanczos2Contract *contract) {
  MVMCPowerLanczosContract powerContract;
  MVMCPowerLanczosContractStatus powerStatus;

  if (contract == NULL || (contract->step != 1 && contract->step != 2)) {
    return LANCZOS2_CONTRACT_INVALID_STEP;
  }
  powerContract = ToPowerContract(contract);
  powerStatus = mvmc_power_lanczos_validate_common_contract(&powerContract);
  if (powerStatus != MVMC_POWER_LANCZOS_CONTRACT_OK) {
    return FromPowerStatus(powerStatus);
  }
  if (contract->step == 1) return LANCZOS2_CONTRACT_OK;

  powerStatus =
      mvmc_power_lanczos_validate_production_contract(&powerContract);
  if (powerStatus != MVMC_POWER_LANCZOS_CONTRACT_OK) {
    return FromPowerStatus(powerStatus);
  }

  return LANCZOS2_CONTRACT_OK;
}

Lanczos2ContractStatus ValidateLegacyLanczos2ExecutionContract(
    const Lanczos2Contract *contract) {
  const Lanczos2ContractStatus status = ValidateLanczos2Contract(contract);
  if (status != LANCZOS2_CONTRACT_OK) return status;
  if (contract->step == 2 &&
      ClassifyLanczos2Model(contract) ==
          LANCZOS2_MODEL_LOCAL_SPIN_EXCHANGE &&
      contract->nQPFull != 1) {
    return LANCZOS2_CONTRACT_UNSUPPORTED_QUANTUM_PROJECTION;
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
    case LANCZOS2_CONTRACT_UNSUPPORTED_TWO_SZ:
      return "NLanczosStep=2 requires 2Sz=0";
    case LANCZOS2_CONTRACT_UNSUPPORTED_BACKFLOW:
      return "NLanczosStep=2 does not support BackFlow";
    case LANCZOS2_CONTRACT_UNSUPPORTED_RBM:
      return "NLanczosStep=2 does not support RBM";
    case LANCZOS2_CONTRACT_UNSUPPORTED_REWEIGHT:
      return "NLanczosStep=2 does not support reweight=1";
    case LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH:
      return "NLanczosStep=2 supports NExUpdatePath=0, or "
             "NExUpdatePath=2 for pure-spin mode";
    case LANCZOS2_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER:
      return "NLanczosStep=2 electronic mode requires at least one "
             "spin-conserving Transfer term";
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
    case LANCZOS2_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN:
      return "NLanczosStep=2 with NExUpdatePath=2 requires "
             "NLocalSpin=Nsite=2*Ne";
    case LANCZOS2_CONTRACT_UNSUPPORTED_TRANSFER:
      return "NLanczosStep=2 pure-spin mode requires NTransfer=0";
    case LANCZOS2_CONTRACT_UNSUPPORTED_QUANTUM_PROJECTION:
      return "legacy pure-spin Lanczos2 execution currently requires "
             "NQPFull=1; the corrected production adapter handles multi-QP "
             "through its dedicated dispatch";
    case LANCZOS2_CONTRACT_INVALID_COUNT:
      return "NLanczosStep requires non-negative counts and positive Nsite, "
             "Ne, and NQPFull";
    case LANCZOS2_CONTRACT_SIZE_OVERFLOW:
      return "NLanczosStep size arithmetic overflows the supported range";
  }
  return "unknown NLanczosStep contract error";
}
