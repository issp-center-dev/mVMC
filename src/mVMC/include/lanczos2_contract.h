#ifndef _LANCZOS2_CONTRACT_H
#define _LANCZOS2_CONTRACT_H

typedef struct {
  int step;
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
} Lanczos2Contract;

typedef enum {
  LANCZOS2_MODEL_INVALID = 0,
  LANCZOS2_MODEL_ELECTRONIC_VK,
  LANCZOS2_MODEL_LOCAL_SPIN_EXCHANGE
} Lanczos2ModelClass;

typedef enum {
  LANCZOS2_CONTRACT_OK = 0,
  LANCZOS2_CONTRACT_INVALID_STEP,
  LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1,
  LANCZOS2_CONTRACT_REQUIRES_VMC_CAL_MODE_1,
  LANCZOS2_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL,
  LANCZOS2_CONTRACT_UNSUPPORTED_TWO_SZ,
  LANCZOS2_CONTRACT_UNSUPPORTED_BACKFLOW,
  LANCZOS2_CONTRACT_UNSUPPORTED_RBM,
  LANCZOS2_CONTRACT_UNSUPPORTED_REWEIGHT,
  LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH,
  LANCZOS2_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER,
  LANCZOS2_CONTRACT_UNSUPPORTED_PAIR_HOPPING,
  LANCZOS2_CONTRACT_UNSUPPORTED_EXCHANGE,
  LANCZOS2_CONTRACT_UNSUPPORTED_INTER_ALL,
  LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL,
  LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_G,
  LANCZOS2_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER,
  LANCZOS2_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN,
  LANCZOS2_CONTRACT_UNSUPPORTED_TRANSFER,
  LANCZOS2_CONTRACT_UNSUPPORTED_QUANTUM_PROJECTION,
  LANCZOS2_CONTRACT_INVALID_COUNT,
  LANCZOS2_CONTRACT_SIZE_OVERFLOW
} Lanczos2ContractStatus;

Lanczos2ModelClass ClassifyLanczos2Model(
    const Lanczos2Contract *contract);

Lanczos2ContractStatus ValidateLanczos2Contract(
    const Lanczos2Contract *contract);

/* Apply the full-support corrected execution capability gate for either
 * power step.  Unlike the legacy step-1 compatibility route, corrected
 * execution must be representable by the classic full-support bridge. */
Lanczos2ContractStatus ValidateCorrectedPowerLanczosExecutionContract(
    const Lanczos2Contract *contract);

/*
 * Apply the additional execution limit of the legacy pure-spin main path.
 * The general production capability contract intentionally accepts multi-QP;
 * corrected adapters use their dedicated dispatch instead of this gate.
 */
Lanczos2ContractStatus ValidateLegacyLanczos2ExecutionContract(
    const Lanczos2Contract *contract);

const char *Lanczos2ContractError(Lanczos2ContractStatus status);

#endif
