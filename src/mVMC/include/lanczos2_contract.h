#ifndef _LANCZOS2_CONTRACT_H
#define _LANCZOS2_CONTRACT_H

typedef struct {
  int step;
  int lanczosMode;
  int vmcCalMode;
  int orbitalGeneral;
  int nProjBF;
  int flagRBM;
  int exUpdatePath;
  int nPairHopping;
  int nExchangeCoupling;
  int nInterAll;
  int nNBodyInterAll;
  int nNBodyG;
  int nSpinFlipTransfer;
} Lanczos2Contract;

typedef enum {
  LANCZOS2_CONTRACT_OK = 0,
  LANCZOS2_CONTRACT_INVALID_STEP,
  LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1,
  LANCZOS2_CONTRACT_REQUIRES_VMC_CAL_MODE_1,
  LANCZOS2_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL,
  LANCZOS2_CONTRACT_UNSUPPORTED_BACKFLOW,
  LANCZOS2_CONTRACT_UNSUPPORTED_RBM,
  LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH,
  LANCZOS2_CONTRACT_UNSUPPORTED_PAIR_HOPPING,
  LANCZOS2_CONTRACT_UNSUPPORTED_EXCHANGE,
  LANCZOS2_CONTRACT_UNSUPPORTED_INTER_ALL,
  LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL,
  LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_G,
  LANCZOS2_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER
} Lanczos2ContractStatus;

Lanczos2ContractStatus ValidateLanczos2Contract(
    const Lanczos2Contract *contract);

const char *Lanczos2ContractError(Lanczos2ContractStatus status);

#endif
