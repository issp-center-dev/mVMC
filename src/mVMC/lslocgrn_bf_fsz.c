/*-------------------------------------------------------------
 * Local Lanczos vector for the BF-FSZ path.
 *-------------------------------------------------------------*/
#include "lslocgrn_bf_fsz.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "calham_fsz.h"
#include "global.h"
#include "locgrn_fsz.h"
#include "matrix.h"
#include "projection.h"
#include "qp.h"
#include "slater.h"

static int lsbffszCheckedMul(size_t left, size_t right, size_t *result) {
  if(result == NULL) return 1;
  if(left != 0 && right > SIZE_MAX/left) return 1;
  *result = left*right;
  return 0;
}

static void *lsbffszCheckedCalloc(size_t count, size_t size) {
  if(count == 0) count = 1;
  if(size != 0 && count > SIZE_MAX/size) return NULL;
  return calloc(count, size);
}

static int lsbffszStateCheckEnabled(void) {
  const char *value = getenv("MVMC_BF_LANCZOS_STATE_CHECK");
  return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0;
}

void LSLanczosBFFSZScratchFree(LSLanczosBFFSZScratch *scratch) {
  if(scratch == NULL) return;
  free(scratch->eleIdx);
  free(scratch->eleCfg);
  free(scratch->eleNum);
  free(scratch->eleSpn);
  free(scratch->projCnt);
  free(scratch->bfCnt);
  free(scratch->gfProjCnt);
  free(scratch->gfBFCnt);
  free(scratch->affected);
  free(scratch->hopIntWork);
  free(scratch->pfIWork);
  free(scratch->stateEleIdx);
  free(scratch->stateEleCfg);
  free(scratch->stateEleNum);
  free(scratch->stateEleSpn);
  free(scratch->stateProjCnt);
  free(scratch->stateBFCnt);
  free(scratch->pfRWork);
  free(scratch->greenBuffer);
  free(scratch->pfBufM);
  free(scratch->pfWork);
  free(scratch->candidateSlater);
  free(scratch->candidateInv);
  free(scratch->candidatePf);
  free(scratch->baseSlater);
  free(scratch->baseInv);
  free(scratch->basePf);
  memset(scratch, 0, sizeof(*scratch));
}

int LSLanczosBFFSZScratchInit(LSLanczosBFFSZScratch *scratch) {
  size_t nsize;
  size_t nsite;
  size_t nsite2;
  size_t nqp;
  size_t nproj;
  size_t nrange;
  size_t invCount;
  size_t slaterCount;
  size_t bfCount;
  size_t matrixCount;
  size_t greenBufferCount;
  size_t pfUpdateIntCount;
  size_t pfUpdateDoubleCount;
  size_t pfIntCount;
  size_t pfDoubleCount;
  int hopIntWorkSize;

  if(scratch == NULL) return 1;
  memset(scratch, 0, sizeof(*scratch));
  if(Nsize <= 0 || Nsite <= 0 || Nsite2 <= 0 || NQPFull <= 0 ||
     NProj < 0 || Nrange <= 0 || NProjBF <= 0 || LapackLWork < Nsize) {
    return 1;
  }
  if(GetGreenFuncBF_fsz_buffer_work_size(
       &greenBufferCount, &pfUpdateIntCount, &pfUpdateDoubleCount) != 0 ||
     GetSlaterElmBF_fsz_hop_int_work_size(&hopIntWorkSize) != 0 ||
     hopIntWorkSize <= 0) return 1;

  nsize = (size_t)Nsize;
  nsite = (size_t)Nsite;
  nsite2 = (size_t)Nsite2;
  nqp = (size_t)NQPFull;
  nproj = (size_t)NProj;
  nrange = (size_t)Nrange;
  if(lsbffszCheckedMul(nsize, nsize, &matrixCount) != 0 ||
     lsbffszCheckedMul(matrixCount, nqp, &invCount) != 0 ||
     lsbffszCheckedMul(nsite2, nsite2, &slaterCount) != 0 ||
     lsbffszCheckedMul(slaterCount, nqp, &slaterCount) != 0 ||
     lsbffszCheckedMul(nsite, nrange, &bfCount) != 0 ||
     lsbffszCheckedMul(bfCount, 16, &bfCount) != 0) return 1;
  pfIntCount = pfUpdateIntCount > nsize ? pfUpdateIntCount : nsize;
  pfDoubleCount = pfUpdateDoubleCount > (size_t)LapackLWork
                ? pfUpdateDoubleCount : (size_t)LapackLWork;

  scratch->stateCheck = lsbffszStateCheckEnabled();
  scratch->hopIntWorkSize = hopIntWorkSize;
  scratch->invCount = invCount;
  scratch->pfCount = nqp;
  scratch->slaterCount = slaterCount;
  scratch->bfCount = bfCount;
  scratch->greenBufferCount = greenBufferCount;
  scratch->pfIntCount = pfIntCount;
  scratch->pfDoubleCount = pfDoubleCount;
  scratch->eleIdx = lsbffszCheckedCalloc(nsize, sizeof(int));
  scratch->eleCfg = lsbffszCheckedCalloc(nsite2, sizeof(int));
  scratch->eleNum = lsbffszCheckedCalloc(nsite2, sizeof(int));
  scratch->eleSpn = lsbffszCheckedCalloc(nsize, sizeof(int));
  scratch->projCnt = lsbffszCheckedCalloc(nproj, sizeof(int));
  scratch->bfCnt = lsbffszCheckedCalloc(bfCount, sizeof(int));
  scratch->gfProjCnt = lsbffszCheckedCalloc(nproj, sizeof(int));
  scratch->gfBFCnt = lsbffszCheckedCalloc(bfCount, sizeof(int));
  scratch->affected = lsbffszCheckedCalloc(nsize, sizeof(int));
  scratch->hopIntWork = lsbffszCheckedCalloc((size_t)hopIntWorkSize,
                                             sizeof(int));
  scratch->pfIWork = lsbffszCheckedCalloc(pfIntCount, sizeof(int));
  scratch->pfRWork = lsbffszCheckedCalloc(pfDoubleCount, sizeof(double));
  scratch->greenBuffer = lsbffszCheckedCalloc(
      greenBufferCount, sizeof(double complex));
  scratch->pfBufM = lsbffszCheckedCalloc(matrixCount,
                                         sizeof(double complex));
  scratch->pfWork = lsbffszCheckedCalloc((size_t)LapackLWork,
                                         sizeof(double complex));
  scratch->candidateSlater = lsbffszCheckedCalloc(
      slaterCount, sizeof(double complex));
  scratch->candidateInv = lsbffszCheckedCalloc(
      invCount, sizeof(double complex));
  scratch->candidatePf = lsbffszCheckedCalloc(nqp, sizeof(double complex));
  scratch->baseSlater = lsbffszCheckedCalloc(
      slaterCount, sizeof(double complex));
  scratch->baseInv = lsbffszCheckedCalloc(invCount, sizeof(double complex));
  scratch->basePf = lsbffszCheckedCalloc(nqp, sizeof(double complex));

  if(scratch->stateCheck) {
    scratch->stateEleIdx = lsbffszCheckedCalloc(nsize, sizeof(int));
    scratch->stateEleCfg = lsbffszCheckedCalloc(nsite2, sizeof(int));
    scratch->stateEleNum = lsbffszCheckedCalloc(nsite2, sizeof(int));
    scratch->stateEleSpn = lsbffszCheckedCalloc(nsize, sizeof(int));
    scratch->stateProjCnt = lsbffszCheckedCalloc(nproj, sizeof(int));
    scratch->stateBFCnt = lsbffszCheckedCalloc(bfCount, sizeof(int));
  }

  if(scratch->eleIdx == NULL || scratch->eleCfg == NULL ||
     scratch->eleNum == NULL || scratch->eleSpn == NULL ||
     scratch->projCnt == NULL || scratch->bfCnt == NULL ||
     scratch->gfProjCnt == NULL || scratch->gfBFCnt == NULL ||
     scratch->affected == NULL || scratch->hopIntWork == NULL ||
     scratch->pfIWork == NULL || scratch->pfRWork == NULL ||
     scratch->greenBuffer == NULL || scratch->pfBufM == NULL ||
     scratch->pfWork == NULL || scratch->candidateSlater == NULL ||
     scratch->candidateInv == NULL || scratch->candidatePf == NULL ||
     scratch->baseSlater == NULL || scratch->baseInv == NULL ||
     scratch->basePf == NULL ||
     (scratch->stateCheck &&
      (scratch->stateEleIdx == NULL || scratch->stateEleCfg == NULL ||
       scratch->stateEleNum == NULL || scratch->stateEleSpn == NULL ||
       scratch->stateProjCnt == NULL || scratch->stateBFCnt == NULL))) {
    LSLanczosBFFSZScratchFree(scratch);
    return 1;
  }
  scratch->initialized = 1;
  return 0;
}

static void lsbffszCopyState(LSLanczosBFFSZScratch *scratch,
                             const int *eleIdx, const int *eleCfg,
                             const int *eleNum, const int *eleProjCnt,
                             const int *eleSpn,
                             const int *eleProjBFCnt) {
  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleSpn, eleSpn, (size_t)Nsize*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->projCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
  memcpy(scratch->bfCnt, eleProjBFCnt, scratch->bfCount*sizeof(int));
}

static void lsbffszCaptureState(LSLanczosBFFSZScratch *scratch,
                                const int *eleIdx, const int *eleCfg,
                                const int *eleNum, const int *eleProjCnt,
                                const int *eleSpn,
                                const int *eleProjBFCnt) {
  if(!scratch->stateCheck) return;
  memcpy(scratch->stateEleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->stateEleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->stateEleNum, eleNum, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->stateEleSpn, eleSpn, (size_t)Nsize*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->stateProjCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
  memcpy(scratch->stateBFCnt, eleProjBFCnt,
         scratch->bfCount*sizeof(int));
}

static int lsbffszVerifyState(const LSLanczosBFFSZScratch *scratch,
                              const int *eleIdx, const int *eleCfg,
                              const int *eleNum, const int *eleProjCnt,
                              const int *eleSpn,
                              const int *eleProjBFCnt) {
  if(!scratch->stateCheck) return 0;
  if(memcmp(scratch->stateEleIdx, eleIdx,
            (size_t)Nsize*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleCfg, eleCfg,
            (size_t)Nsite2*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleNum, eleNum,
            (size_t)Nsite2*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleSpn, eleSpn,
            (size_t)Nsize*sizeof(int)) != 0 ||
     (NProj > 0 && memcmp(scratch->stateProjCnt, eleProjCnt,
                          (size_t)NProj*sizeof(int)) != 0) ||
     memcmp(scratch->stateBFCnt, eleProjBFCnt,
            scratch->bfCount*sizeof(int)) != 0 ||
     memcmp(scratch->baseSlater, SlaterElmBF,
            scratch->slaterCount*sizeof(double complex)) != 0 ||
     memcmp(scratch->baseInv, InvM,
            scratch->invCount*sizeof(double complex)) != 0 ||
     memcmp(scratch->basePf, PfM,
            scratch->pfCount*sizeof(double complex)) != 0) return 1;
  return 0;
}

static void lsbffszRestoreMatrixState(
    const LSLanczosBFFSZScratch *scratch) {
  memcpy(SlaterElmBF, scratch->baseSlater,
         scratch->slaterCount*sizeof(double complex));
  memcpy(InvM, scratch->baseInv,
         scratch->invCount*sizeof(double complex));
  memcpy(PfM, scratch->basePf,
         scratch->pfCount*sizeof(double complex));
}

static int lsbffszApplyHop(LSLanczosBFFSZScratch *scratch,
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
  if(s == t) {
    UpdateProjCnt(rj, ri, s, scratch->projCnt,
                  eleProjCnt, scratch->eleNum);
  } else {
    UpdateProjCnt_fsz(rj, ri, t, s, scratch->projCnt,
                      eleProjCnt, scratch->eleNum);
  }
  MakeProjBFCnt(scratch->bfCnt, scratch->eleNum);
  return 0;
}

static double complex lsbffszGreen1(
    int ri, int rj, int s, int t, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt, LSLanczosBFFSZScratch *scratch) {
  lsbffszCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt,
                   eleSpn, eleProjBFCnt);
  return GreenFunc1BF_fsz2_workspace(
      ri, rj, s, t, ip, scratch->eleIdx, scratch->eleCfg,
      scratch->eleNum, eleProjCnt, scratch->eleSpn,
      scratch->gfProjCnt, eleProjBFCnt, scratch->gfBFCnt,
      scratch->greenBuffer, scratch->pfBufM, scratch->affected,
      scratch->hopIntWork, scratch->hopIntWorkSize, scratch->pfIWork,
      scratch->pfWork, scratch->pfRWork);
}

static double complex lsbffszGreen2(
    int ri, int rj, int rk, int rl, int s, int t, int u, int v,
    double complex ip, const int *eleIdx, const int *eleCfg,
    const int *eleNum, const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt, LSLanczosBFFSZScratch *scratch) {
  lsbffszCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt,
                   eleSpn, eleProjBFCnt);
  return GreenFunc2BF_fsz2(
      ri, rj, rk, rl, s, t, u, v, ip,
      scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
      eleProjCnt, scratch->eleSpn, scratch->gfProjCnt,
      eleProjBFCnt, scratch->gfBFCnt, scratch->greenBuffer,
      scratch->affected, scratch->hopIntWork, scratch->hopIntWorkSize,
      scratch->pfIWork, scratch->pfRWork, scratch->pfBufM,
      scratch->pfWork);
}

static int lsbffszOuter(int ri, int rj, int s, int t,
                        double complex h1, double complex ip,
                        const int *eleIdx, const int *eleCfg,
                        const int *eleNum, const int *eleProjCnt,
                        const int *eleSpn, const int *eleProjBFCnt,
                        LSLanczosBFFSZScratch *scratch,
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

  ratio = lsbffszGreen1(ri, rj, s, t, ip, eleIdx, eleCfg, eleNum,
                        eleProjCnt, eleSpn, eleProjBFCnt, scratch);
  pfRatio = CalculateIP_fcmp(scratch->greenBuffer, 0, NQPFull,
                             MPI_COMM_SELF)/ip;
  useRebuild = cabs(pfRatio) > 1.0e-12;
  if(LSLanczosTestProjectionBranchAuditEnabled() && useRebuild &&
     cabs(pfRatio) <= 1.0e-12 && cabs(ratio) > 1.0e-12) {
    return LSLANCZOS_STATE_FAILURE;
  }
  if(useRebuild) {
    BF_FSZ_MAllResult mAllResult;
    double complex ipNew;
    double complex projRatio;
    double complex energy;


    lsbffszCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt,
                     eleSpn, eleProjBFCnt);
    if(lsbffszApplyHop(scratch, ri, rj, s, t, eleProjCnt) != 0) return 1;
    projRatio = ProjRatio(scratch->projCnt, eleProjCnt);
    MakeSlaterElmBF_fsz_to_serial(scratch->candidateSlater,
                                  scratch->eleNum, scratch->bfCnt);
    if(LSLanczosTestConsumeRebuildFailure()) {
      return LSLANCZOS_NUMERIC_REJECT;
    }
    mAllResult = CalculateMAll_BF_fsz_from_workspace(
        scratch->candidateSlater, scratch->eleIdx, scratch->eleSpn,
        0, NQPFull, scratch->candidatePf, scratch->candidateInv,
        (size_t)Nsize*(size_t)Nsize, scratch->pfBufM,
        scratch->pfIWork, scratch->pfWork, LapackLWork,
        scratch->pfRWork);
    if(mAllResult.status != BF_FSZ_MALL_OK) {
      return mAllResult.status == BF_FSZ_MALL_INVALID_ARGUMENT
           ? LSLANCZOS_STATE_FAILURE
           : LSLANCZOS_NUMERIC_REJECT;
    }
    ipNew = CalculateIP_fcmp(scratch->candidatePf, 0, NQPFull,
                             MPI_COMM_SELF);
    memcpy(SlaterElmBF, scratch->candidateSlater,
           scratch->slaterCount*sizeof(double complex));
    memcpy(InvM, scratch->candidateInv,
           scratch->invCount*sizeof(double complex));
    memcpy(PfM, scratch->candidatePf,
           scratch->pfCount*sizeof(double complex));
    energy = CalculateHamiltonianBF_fsz(
        ipNew, scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        scratch->projCnt, scratch->eleSpn, scratch->bfCnt);
    result = energy*conj(projRatio*ipNew/ip);
    lsbffszRestoreMatrixState(scratch);
    *value = result;
    return 0;
  }

  lsbffszCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt,
                   eleSpn, eleProjBFCnt);
  if(lsbffszApplyHop(scratch, ri, rj, s, t, eleProjCnt) != 0) return 1;
  result = ratio*CalculateHamiltonian0_fsz(scratch->eleNum);
  for(idx=0; idx<NTransfer; idx++) {
    const double complex green = lsbffszGreen2(
        Transfer[idx][0], Transfer[idx][2], ri, rj,
        Transfer[idx][1], Transfer[idx][3], s, t, ip,
        eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn, eleProjBFCnt,
        scratch);
    result -= ParaTransfer[idx]*green;
  }
  *value = result;
  return 0;
}

int LSLocalQBF_fsz(const double complex h1, const double complex ip,
                   const int *eleIdx, const int *eleCfg, const int *eleNum,
                   const int *eleProjCnt, const int *eleSpn,
                   const int *eleProjBFCnt,
                   LSLanczosBFFSZScratch *scratch, double complex *lslq) {
  double complex h2;
  const int savedBFProfileEnabled = BFProfileEnabled;
  const int savedInvGemmCheckEnabled = BFFSZInvGemmCheckEnabled;
  const int savedInvDetailProfileEnabled = BFFSZInvDetailProfileEnabled;
  const int savedC2DetailProfileEnabled = BFFSZC2DetailProfileEnabled;
  int status = 0;
  int idx;

  if(scratch == NULL || !scratch->initialized || lslq == NULL ||
     !isfinite(creal(ip)) || !isfinite(cimag(ip)) || cabs(ip) == 0.0) {
    return 1;
  }

  memcpy(scratch->baseSlater, SlaterElmBF,
         scratch->slaterCount*sizeof(double complex));
  memcpy(scratch->baseInv, InvM,
         scratch->invCount*sizeof(double complex));
  memcpy(scratch->basePf, PfM,
         scratch->pfCount*sizeof(double complex));
  lsbffszCaptureState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt,
                      eleSpn, eleProjBFCnt);

  BFProfileEnabled = 0;
  BFFSZInvGemmCheckEnabled = 0;
  BFFSZInvDetailProfileEnabled = 0;
  BFFSZC2DetailProfileEnabled = 0;

  h2 = h1*CalculateHamiltonian0_fsz(eleNum);
  for(idx=0; idx<NTransfer; idx++) {
    double complex outer;
    const int outerStatus = lsbffszOuter(
        Transfer[idx][0], Transfer[idx][2],
        Transfer[idx][1], Transfer[idx][3], h1, ip,
        eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
        eleProjBFCnt, scratch, &outer);
    if(outerStatus != LSLANCZOS_OK) {
      status = outerStatus;
      break;
    }
    h2 -= ParaTransfer[idx]*outer;
  }

  lsbffszRestoreMatrixState(scratch);
  BFProfileEnabled = savedBFProfileEnabled;
  BFFSZInvGemmCheckEnabled = savedInvGemmCheckEnabled;
  BFFSZInvDetailProfileEnabled = savedInvDetailProfileEnabled;
  BFFSZC2DetailProfileEnabled = savedC2DetailProfileEnabled;

  if(lsbffszVerifyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt,
                        eleSpn, eleProjBFCnt) != 0) {
    return LSLANCZOS_STATE_FAILURE;
  }
  if(status != LSLANCZOS_OK) return status;

  lslq[0] = 1.0;
  lslq[1] = h1;
  lslq[2] = h1;
  lslq[3] = h2;
  return LSLANCZOS_OK;
}

int LSLocalQBF_fsz_real(const double complex h1, const double complex ip,
                        const int *eleIdx, const int *eleCfg,
                        const int *eleNum, const int *eleProjCnt,
                        const int *eleSpn, const int *eleProjBFCnt,
                        LSLanczosBFFSZScratch *scratch, double *lslq) {
  double complex complexLSLQ[4];
  int status;
  int i;
  if(scratch == NULL || lslq == NULL || NLSHam*NLSHam != 4) return 1;
  scratch->realGateIndex = -1;
  scratch->realGateTolerance = 0.0;
  scratch->realGateValue = 0.0;
  status = LSLocalQBF_fsz(h1, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
                          eleSpn, eleProjBFCnt, scratch, complexLSLQ);
  if(status != 0) return status;
  for(i=0; i<4; i++) {
    const double realPart = creal(complexLSLQ[i]);
    const double imagPart = cimag(complexLSLQ[i]);
    const double tolerance = 1.0e-10 + 1.0e-9*fabs(realPart);
    if(!isfinite(realPart) || !isfinite(imagPart)) {
      lslq[i] = NAN;
      continue;
    }
    if(fabs(imagPart) > tolerance) {
      scratch->realGateIndex = i;
      scratch->realGateTolerance = tolerance;
      scratch->realGateValue = complexLSLQ[i];
      return LSLANCZOS_REAL_GATE_FAILURE;
    }
    lslq[i] = realPart;
  }
  return LSLANCZOS_OK;
}
