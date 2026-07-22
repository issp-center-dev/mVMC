/*-------------------------------------------------------------
 * Local Lanczos vector for the no-BackFlow FSZ path.
 *-------------------------------------------------------------*/
#include "lslocgrn_fsz.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "calham_fsz.h"
#include "calham_fsz_real.h"
#include "global.h"
#include "locgrn_fsz.h"
#include "locgrn_fsz_real.h"
#include "matrix.h"
#include "projection.h"
#include "qp.h"
#include "qp_real.h"

static int lsCheckedMul(size_t left, size_t right, size_t *result) {
  if(result == NULL) return 1;
  if(left != 0 && right > SIZE_MAX/left) return 1;
  *result = left*right;
  return 0;
}

static void *lsCheckedCalloc(size_t count, size_t size) {
  if(count == 0) count = 1;
  if(size != 0 && count > SIZE_MAX/size) return NULL;
  return calloc(count, size);
}

static int lsStateCheckEnabled(void) {
  const char *value = getenv("MVMC_BF_LANCZOS_STATE_CHECK");
  return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0;
}

void LSLanczosFSZScratchFree(LSLanczosFSZScratch *scratch) {
  if(scratch == NULL) return;
  free(scratch->eleIdx);
  free(scratch->eleCfg);
  free(scratch->eleNum);
  free(scratch->eleSpn);
  free(scratch->projCnt);
  free(scratch->gfProjCnt);
  free(scratch->stateEleIdx);
  free(scratch->stateEleCfg);
  free(scratch->stateEleNum);
  free(scratch->stateEleSpn);
  free(scratch->stateProjCnt);
  free(scratch->gfComplex);
  free(scratch->baseInvComplex);
  free(scratch->basePfComplex);
  free(scratch->gfReal);
  free(scratch->baseInvReal);
  free(scratch->basePfReal);
  memset(scratch, 0, sizeof(*scratch));
}

int LSLanczosFSZScratchInit(LSLanczosFSZScratch *scratch, int useReal) {
  size_t nsize;
  size_t nsite2;
  size_t nqp;
  size_t nproj;
  size_t matrixCount;
  size_t gfCount;

  if(scratch == NULL) return 1;
  memset(scratch, 0, sizeof(*scratch));
  if(Nsize <= 0 || Nsite2 <= 0 || NQPFull <= 0 || NProj < 0) return 1;

  nsize = (size_t)Nsize;
  nsite2 = (size_t)Nsite2;
  nqp = (size_t)NQPFull;
  nproj = (size_t)NProj;
  if(lsCheckedMul(nsize, nsize, &matrixCount) != 0 ||
     lsCheckedMul(matrixCount, nqp, &matrixCount) != 0 ||
     nqp > SIZE_MAX - 2*nsize) {
    return 1;
  }
  gfCount = nqp + 2*nsize;

  scratch->eleIdx = lsCheckedCalloc(nsize, sizeof(int));
  scratch->eleCfg = lsCheckedCalloc(nsite2, sizeof(int));
  scratch->eleNum = lsCheckedCalloc(nsite2, sizeof(int));
  scratch->eleSpn = lsCheckedCalloc(nsize, sizeof(int));
  scratch->projCnt = lsCheckedCalloc(nproj, sizeof(int));
  scratch->gfProjCnt = lsCheckedCalloc(nproj, sizeof(int));
  scratch->useReal = useReal != 0;
  scratch->stateCheck = lsStateCheckEnabled();
  scratch->invCount = matrixCount;
  scratch->pfCount = nqp;
  scratch->gfCount = gfCount;

  if(scratch->useReal) {
    scratch->gfReal = lsCheckedCalloc(gfCount, sizeof(double));
    scratch->baseInvReal = lsCheckedCalloc(matrixCount, sizeof(double));
    scratch->basePfReal = lsCheckedCalloc(nqp, sizeof(double));
  } else {
    scratch->gfComplex = lsCheckedCalloc(gfCount, sizeof(double complex));
    scratch->baseInvComplex = lsCheckedCalloc(matrixCount, sizeof(double complex));
    scratch->basePfComplex = lsCheckedCalloc(nqp, sizeof(double complex));
  }

  if(scratch->stateCheck) {
    scratch->stateEleIdx = lsCheckedCalloc(nsize, sizeof(int));
    scratch->stateEleCfg = lsCheckedCalloc(nsite2, sizeof(int));
    scratch->stateEleNum = lsCheckedCalloc(nsite2, sizeof(int));
    scratch->stateEleSpn = lsCheckedCalloc(nsize, sizeof(int));
    scratch->stateProjCnt = lsCheckedCalloc(nproj, sizeof(int));
  }

  if(scratch->eleIdx == NULL || scratch->eleCfg == NULL ||
     scratch->eleNum == NULL || scratch->eleSpn == NULL ||
     scratch->projCnt == NULL || scratch->gfProjCnt == NULL ||
     (scratch->useReal && (scratch->gfReal == NULL ||
                           scratch->baseInvReal == NULL ||
                           scratch->basePfReal == NULL)) ||
     (!scratch->useReal && (scratch->gfComplex == NULL ||
                            scratch->baseInvComplex == NULL ||
                            scratch->basePfComplex == NULL)) ||
     (scratch->stateCheck &&
      (scratch->stateEleIdx == NULL || scratch->stateEleCfg == NULL ||
       scratch->stateEleNum == NULL || scratch->stateEleSpn == NULL ||
       scratch->stateProjCnt == NULL))) {
    LSLanczosFSZScratchFree(scratch);
    return 1;
  }

  scratch->initialized = 1;
  return 0;
}

static void lsCopyElectronState(LSLanczosFSZScratch *scratch,
                                const int *eleIdx, const int *eleCfg,
                                const int *eleNum, const int *eleProjCnt,
                                const int *eleSpn) {
  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleSpn, eleSpn, (size_t)Nsize*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->projCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
}

static void lsCaptureState(LSLanczosFSZScratch *scratch,
                           const int *eleIdx, const int *eleCfg,
                           const int *eleNum, const int *eleProjCnt,
                           const int *eleSpn) {
  if(!scratch->stateCheck) return;
  memcpy(scratch->stateEleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->stateEleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->stateEleNum, eleNum, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->stateEleSpn, eleSpn, (size_t)Nsize*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->stateProjCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
}

static int lsVerifyComplexState(LSLanczosFSZScratch *scratch,
                                const int *eleIdx, const int *eleCfg,
                                const int *eleNum, const int *eleProjCnt,
                                const int *eleSpn) {
  if(!scratch->stateCheck) return 0;
  if(memcmp(scratch->stateEleIdx, eleIdx, (size_t)Nsize*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleCfg, eleCfg, (size_t)Nsite2*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleNum, eleNum, (size_t)Nsite2*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleSpn, eleSpn, (size_t)Nsize*sizeof(int)) != 0 ||
     (NProj > 0 && memcmp(scratch->stateProjCnt, eleProjCnt,
                          (size_t)NProj*sizeof(int)) != 0) ||
     memcmp(scratch->baseInvComplex, InvM,
            scratch->invCount*sizeof(double complex)) != 0 ||
     memcmp(scratch->basePfComplex, PfM,
            scratch->pfCount*sizeof(double complex)) != 0) {
    return 1;
  }
  return 0;
}

static int lsVerifyRealState(LSLanczosFSZScratch *scratch,
                             const int *eleIdx, const int *eleCfg,
                             const int *eleNum, const int *eleProjCnt,
                             const int *eleSpn) {
  if(!scratch->stateCheck) return 0;
  if(memcmp(scratch->stateEleIdx, eleIdx, (size_t)Nsize*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleCfg, eleCfg, (size_t)Nsite2*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleNum, eleNum, (size_t)Nsite2*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleSpn, eleSpn, (size_t)Nsize*sizeof(int)) != 0 ||
     (NProj > 0 && memcmp(scratch->stateProjCnt, eleProjCnt,
                          (size_t)NProj*sizeof(int)) != 0) ||
     memcmp(scratch->baseInvReal, InvM_real,
            scratch->invCount*sizeof(double)) != 0 ||
     memcmp(scratch->basePfReal, PfM_real,
            scratch->pfCount*sizeof(double)) != 0) {
    return 1;
  }
  return 0;
}

static int lsApplyOuterHop(LSLanczosFSZScratch *scratch,
                           int ri, int rj, int s, int t,
                           const int *eleProjCnt) {
  const int source = rj + t*Nsite;
  const int target = ri + s*Nsite;
  const int electron = scratch->eleCfg[source];
  if(electron < 0 || electron >= Nsize) return 1;

  scratch->eleIdx[electron] = ri;
  scratch->eleSpn[electron] = s;
  scratch->eleCfg[source] = -1;
  scratch->eleCfg[target] = electron;
  scratch->eleNum[source] = 0;
  scratch->eleNum[target] = 1;
  UpdateProjCnt_fsz(rj, ri, t, s, scratch->projCnt,
                    eleProjCnt, scratch->eleNum);
  return 0;
}

static int lsOuterComplex(int ri, int rj, int s, int t,
                          double complex h1, double complex ip,
                          const int *eleIdx, const int *eleCfg,
                          const int *eleNum, const int *eleProjCnt,
                          const int *eleSpn, LSLanczosFSZScratch *scratch,
                          double complex *value) {
  const int target = ri + s*Nsite;
  const int source = rj + t*Nsite;
  double complex ratio;
  double complex pfRatio;
  double complex result;
  int useRebuild;
  int idx;

  if(target == source) {
    *value = eleNum[target] ? h1 : 0.0;
    return 0;
  }
  if(eleNum[target] != 0 || eleNum[source] == 0) {
    *value = 0.0;
    return 0;
  }

  lsCopyElectronState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);
  ratio = GreenFunc1_fsz2(ri, rj, s, t, ip,
                          scratch->eleIdx, scratch->eleCfg,
                          scratch->eleNum, eleProjCnt, scratch->eleSpn,
                          scratch->gfProjCnt, scratch->gfComplex);
  pfRatio = CalculateIP_fcmp(scratch->gfComplex, 0, NQPFull,
                             MPI_COMM_SELF)/ip;
  useRebuild = cabs(pfRatio) > 1.0e-12;
  if(LSLanczosTestProjectionBranchAuditEnabled() && useRebuild &&
     cabs(pfRatio) <= 1.0e-12 && cabs(ratio) > 1.0e-12) {
    return LSLANCZOS_STATE_FAILURE;
  }
  if(useRebuild) {
    double complex ipNew;
    double complex projRatio;
    double complex energy;
    int info;

    lsCopyElectronState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);
    if(lsApplyOuterHop(scratch, ri, rj, s, t, eleProjCnt) != 0) return 1;
    projRatio = ProjRatio(scratch->projCnt, eleProjCnt);
    info = LSLanczosTestConsumeRebuildFailure()
         ? 1
         : CalculateMAll_fsz(scratch->eleIdx, scratch->eleSpn,
                             0, NQPFull);
    if(info != 0) {
      memcpy(InvM, scratch->baseInvComplex,
             scratch->invCount*sizeof(double complex));
      memcpy(PfM, scratch->basePfComplex,
             scratch->pfCount*sizeof(double complex));
      return LSLANCZOS_NUMERIC_REJECT;
    }
    ipNew = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
    energy = CalculateHamiltonian_fsz(ipNew, scratch->eleIdx,
                                      scratch->eleCfg, scratch->eleNum,
                                      scratch->projCnt, scratch->eleSpn);
    result = energy*conj(projRatio*ipNew/ip);
    memcpy(InvM, scratch->baseInvComplex,
           scratch->invCount*sizeof(double complex));
    memcpy(PfM, scratch->basePfComplex,
           scratch->pfCount*sizeof(double complex));
    *value = result;
    return 0;
  }

  lsCopyElectronState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);
  if(lsApplyOuterHop(scratch, ri, rj, s, t, eleProjCnt) != 0) return 1;
  result = ratio*CalculateHamiltonian0_fsz(scratch->eleNum);
  for(idx=0; idx<NTransfer; idx++) {
    const int rk = Transfer[idx][0];
    const int u = Transfer[idx][1];
    const int rl = Transfer[idx][2];
    const int v = Transfer[idx][3];
    double complex green;
    lsCopyElectronState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);
    green = GreenFunc2_fsz2(rk, rl, ri, rj, u, v, s, t, ip,
                            scratch->eleIdx, scratch->eleCfg,
                            scratch->eleNum, eleProjCnt, scratch->eleSpn,
                            scratch->gfProjCnt, scratch->gfComplex);
    result -= ParaTransfer[idx]*green;
  }
  *value = result;
  return 0;
}

static int lsOuterReal(int ri, int rj, int s, int t,
                       double h1, double ip, const int *eleIdx,
                       const int *eleCfg, const int *eleNum,
                       const int *eleProjCnt, const int *eleSpn,
                       LSLanczosFSZScratch *scratch, double *value) {
  const int target = ri + s*Nsite;
  const int source = rj + t*Nsite;
  double ratio;
  double pfRatio;
  double result;
  int useRebuild;
  int idx;

  if(target == source) {
    *value = eleNum[target] ? h1 : 0.0;
    return 0;
  }
  if(eleNum[target] != 0 || eleNum[source] == 0) {
    *value = 0.0;
    return 0;
  }

  lsCopyElectronState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);
  ratio = GreenFunc1_fsz2_real(ri, rj, s, t, ip,
                               scratch->eleIdx, scratch->eleCfg,
                               scratch->eleNum, eleProjCnt,
                               scratch->eleSpn, scratch->gfProjCnt,
                               scratch->gfReal);
  pfRatio = CalculateIP_real(scratch->gfReal, 0, NQPFull,
                             MPI_COMM_SELF)/ip;
  useRebuild = fabs(pfRatio) > 1.0e-12;
  if(LSLanczosTestProjectionBranchAuditEnabled() && useRebuild &&
     fabs(pfRatio) <= 1.0e-12 && fabs(ratio) > 1.0e-12) {
    return LSLANCZOS_STATE_FAILURE;
  }
  if(useRebuild) {
    double ipNew;
    double projRatio;
    double energy;
    int info;


    lsCopyElectronState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);
    if(lsApplyOuterHop(scratch, ri, rj, s, t, eleProjCnt) != 0) return 1;
    projRatio = ProjRatio(scratch->projCnt, eleProjCnt);
    info = LSLanczosTestConsumeRebuildFailure()
         ? 1
         : CalculateMAll_fsz_real(scratch->eleIdx, scratch->eleSpn,
                                  0, NQPFull);
    if(info != 0) {
      memcpy(InvM_real, scratch->baseInvReal,
             scratch->invCount*sizeof(double));
      memcpy(PfM_real, scratch->basePfReal,
             scratch->pfCount*sizeof(double));
      return LSLANCZOS_NUMERIC_REJECT;
    }
    ipNew = CalculateIP_real(PfM_real, 0, NQPFull, MPI_COMM_SELF);
    energy = CalculateHamiltonian_fsz_real(ipNew, scratch->eleIdx,
                                           scratch->eleCfg, scratch->eleNum,
                                           scratch->projCnt, scratch->eleSpn);
    result = energy*projRatio*ipNew/ip;
    memcpy(InvM_real, scratch->baseInvReal,
           scratch->invCount*sizeof(double));
    memcpy(PfM_real, scratch->basePfReal,
           scratch->pfCount*sizeof(double));
    *value = result;
    return 0;
  }

  lsCopyElectronState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);
  if(lsApplyOuterHop(scratch, ri, rj, s, t, eleProjCnt) != 0) return 1;
  result = ratio*creal(CalculateHamiltonian0_fsz(scratch->eleNum));
  for(idx=0; idx<NTransfer; idx++) {
    const int rk = Transfer[idx][0];
    const int u = Transfer[idx][1];
    const int rl = Transfer[idx][2];
    const int v = Transfer[idx][3];
    double green;
    lsCopyElectronState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);
    green = GreenFunc2_fsz2_real(rk, rl, ri, rj, u, v, s, t, ip,
                                 scratch->eleIdx, scratch->eleCfg,
                                 scratch->eleNum, eleProjCnt,
                                 scratch->eleSpn, scratch->gfProjCnt,
                                 scratch->gfReal);
    result -= creal(ParaTransfer[idx])*green;
  }
  *value = result;
  return 0;
}

int LSLocalQ_fsz(const double complex h1, const double complex ip,
                 const int *eleIdx, const int *eleCfg, const int *eleNum,
                 const int *eleProjCnt, const int *eleSpn,
                 LSLanczosFSZScratch *scratch, double complex *lslq) {
  double complex h2;
  int status;
  int idx;
  if(scratch == NULL || !scratch->initialized || scratch->useReal ||
     lslq == NULL) return 1;

  memcpy(scratch->baseInvComplex, InvM,
         scratch->invCount*sizeof(double complex));
  memcpy(scratch->basePfComplex, PfM,
         scratch->pfCount*sizeof(double complex));
  lsCaptureState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);

  h2 = h1*CalculateHamiltonian0_fsz(eleNum);
  for(idx=0; idx<NTransfer; idx++) {
    double complex outer;
    status = lsOuterComplex(Transfer[idx][0], Transfer[idx][2],
                            Transfer[idx][1], Transfer[idx][3], h1, ip,
                            eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
                            scratch, &outer);
    if(status != LSLANCZOS_OK) {
      memcpy(InvM, scratch->baseInvComplex,
             scratch->invCount*sizeof(double complex));
      memcpy(PfM, scratch->basePfComplex,
             scratch->pfCount*sizeof(double complex));
      if(lsVerifyComplexState(scratch, eleIdx, eleCfg, eleNum,
                              eleProjCnt, eleSpn) != 0) {
        return LSLANCZOS_STATE_FAILURE;
      }
      return status;
    }
    h2 -= ParaTransfer[idx]*outer;
  }

  lslq[0] = 1.0;
  lslq[1] = h1;
  lslq[2] = h1;
  lslq[3] = h2;
  if(lsVerifyComplexState(scratch, eleIdx, eleCfg, eleNum,
                          eleProjCnt, eleSpn) != 0) {
    return LSLANCZOS_STATE_FAILURE;
  }
  return LSLANCZOS_OK;
}

int LSLocalQ_fsz_real(const double h1, const double ip,
                      const int *eleIdx, const int *eleCfg,
                      const int *eleNum, const int *eleProjCnt,
                      const int *eleSpn, LSLanczosFSZScratch *scratch,
                      double *lslq) {
  double h2;
  int status;
  int idx;
  if(scratch == NULL || !scratch->initialized || !scratch->useReal ||
     lslq == NULL) return 1;

  memcpy(scratch->baseInvReal, InvM_real,
         scratch->invCount*sizeof(double));
  memcpy(scratch->basePfReal, PfM_real,
         scratch->pfCount*sizeof(double));
  lsCaptureState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn);

  h2 = h1*creal(CalculateHamiltonian0_fsz(eleNum));
  for(idx=0; idx<NTransfer; idx++) {
    double outer;
    status = lsOuterReal(Transfer[idx][0], Transfer[idx][2],
                         Transfer[idx][1], Transfer[idx][3], h1, ip,
                         eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
                         scratch, &outer);
    if(status != LSLANCZOS_OK) {
      memcpy(InvM_real, scratch->baseInvReal,
             scratch->invCount*sizeof(double));
      memcpy(PfM_real, scratch->basePfReal,
             scratch->pfCount*sizeof(double));
      if(lsVerifyRealState(scratch, eleIdx, eleCfg, eleNum,
                           eleProjCnt, eleSpn) != 0) {
        return LSLANCZOS_STATE_FAILURE;
      }
      return status;
    }
    h2 -= creal(ParaTransfer[idx])*outer;
  }

  lslq[0] = 1.0;
  lslq[1] = h1;
  lslq[2] = h1;
  lslq[3] = h2;
  if(lsVerifyRealState(scratch, eleIdx, eleCfg, eleNum,
                       eleProjCnt, eleSpn) != 0) {
    return LSLANCZOS_STATE_FAILURE;
  }
  return LSLANCZOS_OK;
}
