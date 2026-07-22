/*-------------------------------------------------------------
 * Local Lanczos vector for the non-FSZ BackFlow path.
 *-------------------------------------------------------------*/
#include "lslocgrn_bf.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "calham.h"
#include "calham_real.h"
#include "global.h"
#include "locgrn.h"
#include "locgrn_real.h"
#include "matrix.h"
#include "projection.h"
#include "qp.h"
#include "qp_real.h"
#include "slater.h"

static int lsbfCheckedMul(size_t left, size_t right, size_t *result) {
  if(result == NULL) return 1;
  if(left != 0 && right > SIZE_MAX/left) return 1;
  *result = left*right;
  return 0;
}

static void *lsbfCheckedCalloc(size_t count, size_t size) {
  if(count == 0) count = 1;
  if(size != 0 && count > SIZE_MAX/size) return NULL;
  return calloc(count, size);
}

static int lsbfStateCheckEnabled(void) {
  const char *value = getenv("MVMC_BF_LANCZOS_STATE_CHECK");
  return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0;
}

void LSLanczosBFScratchFree(LSLanczosBFScratch *scratch) {
  if(scratch == NULL) return;
  free(scratch->eleIdx);
  free(scratch->eleCfg);
  free(scratch->eleNum);
  free(scratch->projCnt);
  free(scratch->bfCnt);
  free(scratch->gfProjCnt);
  free(scratch->gfBFCnt);
  free(scratch->stateEleIdx);
  free(scratch->stateEleCfg);
  free(scratch->stateEleNum);
  free(scratch->stateProjCnt);
  free(scratch->stateBFCnt);
  free(scratch->baseEtaFlag);
  free(scratch->baseEta);
  free(scratch->baseSlaterComplex);
  free(scratch->gfSlaterComplex);
  free(scratch->gfComplex);
  free(scratch->baseInvComplex);
  free(scratch->basePfComplex);
  free(scratch->baseSlaterReal);
  free(scratch->gfSlaterReal);
  free(scratch->gfReal);
  free(scratch->baseInvReal);
  free(scratch->basePfReal);
  memset(scratch, 0, sizeof(*scratch));
}

int LSLanczosBFScratchInit(LSLanczosBFScratch *scratch, int useReal) {
  size_t nsize;
  size_t nsite;
  size_t nsite2;
  size_t nqp;
  size_t nproj;
  size_t nrange;
  size_t invCount;
  size_t slaterCount;
  size_t bfCount;
  size_t etaCount;
  size_t gfCount;

  if(scratch == NULL) return 1;
  memset(scratch, 0, sizeof(*scratch));
  if(Nsize <= 0 || Nsite <= 0 || Nsite2 <= 0 || NQPFull <= 0 ||
     NProj < 0 || Nrange <= 0 || NProjBF <= 0) return 1;

  nsize = (size_t)Nsize;
  nsite = (size_t)Nsite;
  nsite2 = (size_t)Nsite2;
  nqp = (size_t)NQPFull;
  nproj = (size_t)NProj;
  nrange = (size_t)Nrange;
  if(lsbfCheckedMul(nsize, nsize, &invCount) != 0 ||
     lsbfCheckedMul(invCount, nqp, &invCount) != 0 ||
     lsbfCheckedMul(nsite2, nsite2, &slaterCount) != 0 ||
     lsbfCheckedMul(slaterCount, nqp, &slaterCount) != 0 ||
     lsbfCheckedMul(nsite, nrange, &bfCount) != 0 ||
     lsbfCheckedMul(bfCount, 16, &bfCount) != 0 ||
     lsbfCheckedMul(nsite, nsite, &etaCount) != 0 ||
     nqp > SIZE_MAX - 2*nsize) return 1;
  gfCount = nqp + 2*nsize;

  scratch->useReal = useReal != 0;
  scratch->stateCheck = lsbfStateCheckEnabled();
  scratch->invCount = invCount;
  scratch->pfCount = nqp;
  scratch->slaterCount = slaterCount;
  scratch->bfCount = bfCount;
  scratch->etaCount = etaCount;
  scratch->gfCount = gfCount;
  scratch->eleIdx = lsbfCheckedCalloc(nsize, sizeof(int));
  scratch->eleCfg = lsbfCheckedCalloc(nsite2, sizeof(int));
  scratch->eleNum = lsbfCheckedCalloc(nsite2, sizeof(int));
  scratch->projCnt = lsbfCheckedCalloc(nproj, sizeof(int));
  scratch->bfCnt = lsbfCheckedCalloc(bfCount, sizeof(int));
  scratch->gfProjCnt = lsbfCheckedCalloc(nproj, sizeof(int));
  scratch->gfBFCnt = lsbfCheckedCalloc(bfCount, sizeof(int));
  scratch->baseEtaFlag = lsbfCheckedCalloc(etaCount, sizeof(int));
  scratch->baseEta = lsbfCheckedCalloc(etaCount, sizeof(double complex));
  scratch->baseSlaterComplex = lsbfCheckedCalloc(
      slaterCount, sizeof(double complex));

  if(scratch->useReal) {
    scratch->baseSlaterReal = lsbfCheckedCalloc(slaterCount, sizeof(double));
    scratch->gfSlaterReal = lsbfCheckedCalloc(slaterCount, sizeof(double));
    scratch->gfReal = lsbfCheckedCalloc(gfCount, sizeof(double));
    scratch->baseInvReal = lsbfCheckedCalloc(invCount, sizeof(double));
    scratch->basePfReal = lsbfCheckedCalloc(nqp, sizeof(double));
  } else {
    scratch->gfSlaterComplex = lsbfCheckedCalloc(
        slaterCount, sizeof(double complex));
    scratch->gfComplex = lsbfCheckedCalloc(gfCount, sizeof(double complex));
    scratch->baseInvComplex = lsbfCheckedCalloc(
        invCount, sizeof(double complex));
    scratch->basePfComplex = lsbfCheckedCalloc(nqp, sizeof(double complex));
  }

  if(scratch->stateCheck) {
    scratch->stateEleIdx = lsbfCheckedCalloc(nsize, sizeof(int));
    scratch->stateEleCfg = lsbfCheckedCalloc(nsite2, sizeof(int));
    scratch->stateEleNum = lsbfCheckedCalloc(nsite2, sizeof(int));
    scratch->stateProjCnt = lsbfCheckedCalloc(nproj, sizeof(int));
    scratch->stateBFCnt = lsbfCheckedCalloc(bfCount, sizeof(int));
  }

  if(scratch->eleIdx == NULL || scratch->eleCfg == NULL ||
     scratch->eleNum == NULL || scratch->projCnt == NULL ||
     scratch->bfCnt == NULL || scratch->gfProjCnt == NULL ||
     scratch->gfBFCnt == NULL || scratch->baseEtaFlag == NULL ||
     scratch->baseEta == NULL || scratch->baseSlaterComplex == NULL ||
     (scratch->useReal &&
      (scratch->baseSlaterReal == NULL || scratch->gfSlaterReal == NULL ||
       scratch->gfReal == NULL || scratch->baseInvReal == NULL ||
       scratch->basePfReal == NULL)) ||
     (!scratch->useReal &&
      (scratch->gfSlaterComplex == NULL || scratch->gfComplex == NULL ||
       scratch->baseInvComplex == NULL || scratch->basePfComplex == NULL)) ||
     (scratch->stateCheck &&
      (scratch->stateEleIdx == NULL || scratch->stateEleCfg == NULL ||
       scratch->stateEleNum == NULL || scratch->stateProjCnt == NULL ||
       scratch->stateBFCnt == NULL))) {
    LSLanczosBFScratchFree(scratch);
    return 1;
  }
  scratch->initialized = 1;
  return 0;
}

static void lsbfCopyState(LSLanczosBFScratch *scratch,
                          const int *eleIdx, const int *eleCfg,
                          const int *eleNum, const int *eleProjCnt) {
  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->projCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
}

static void lsbfCopyEtaTo(double complex *etaOut, int *flagOut) {
  int i, j;
  for(i=0; i<Nsite; i++) {
    for(j=0; j<Nsite; j++) {
      etaOut[(size_t)i*(size_t)Nsite + (size_t)j] = eta[i][j];
      flagOut[(size_t)i*(size_t)Nsite + (size_t)j] = etaFlag[i][j];
    }
  }
}

static void lsbfRestoreEta(const LSLanczosBFScratch *scratch) {
  int i, j;
  for(i=0; i<Nsite; i++) {
    for(j=0; j<Nsite; j++) {
      eta[i][j] = scratch->baseEta[(size_t)i*(size_t)Nsite + (size_t)j];
      etaFlag[i][j] = scratch->baseEtaFlag[
          (size_t)i*(size_t)Nsite + (size_t)j];
    }
  }
}

static void lsbfCaptureState(LSLanczosBFScratch *scratch,
                             const int *eleIdx, const int *eleCfg,
                             const int *eleNum, const int *eleProjCnt,
                             const int *eleProjBFCnt) {
  if(!scratch->stateCheck) return;
  memcpy(scratch->stateEleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->stateEleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->stateEleNum, eleNum, (size_t)Nsite2*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->stateProjCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
  memcpy(scratch->stateBFCnt, eleProjBFCnt,
         scratch->bfCount*sizeof(int));
}

static int lsbfVerifyCommonState(const LSLanczosBFScratch *scratch,
                                 const int *eleIdx, const int *eleCfg,
                                 const int *eleNum, const int *eleProjCnt,
                                 const int *eleProjBFCnt) {
  size_t idx;
  if(!scratch->stateCheck) return 0;
  if(memcmp(scratch->stateEleIdx, eleIdx, (size_t)Nsize*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleCfg, eleCfg, (size_t)Nsite2*sizeof(int)) != 0 ||
     memcmp(scratch->stateEleNum, eleNum, (size_t)Nsite2*sizeof(int)) != 0 ||
     (NProj > 0 && memcmp(scratch->stateProjCnt, eleProjCnt,
                          (size_t)NProj*sizeof(int)) != 0) ||
     memcmp(scratch->stateBFCnt, eleProjBFCnt,
            scratch->bfCount*sizeof(int)) != 0 ||
     memcmp(scratch->baseSlaterComplex, SlaterElmBF,
            scratch->slaterCount*sizeof(double complex)) != 0) return 1;
  for(idx=0; idx<scratch->etaCount; idx++) {
    const size_t i = idx/(size_t)Nsite;
    const size_t j = idx%(size_t)Nsite;
    if(scratch->baseEta[idx] != eta[i][j] ||
       scratch->baseEtaFlag[idx] != etaFlag[i][j]) return 1;
  }
  return 0;
}

static void lsbfRestoreComplex(LSLanczosBFScratch *scratch) {
  memcpy(InvM, scratch->baseInvComplex,
         scratch->invCount*sizeof(double complex));
  memcpy(PfM, scratch->basePfComplex,
         scratch->pfCount*sizeof(double complex));
  memcpy(SlaterElmBF, scratch->baseSlaterComplex,
         scratch->slaterCount*sizeof(double complex));
  lsbfRestoreEta(scratch);
}

static void lsbfRestoreReal(LSLanczosBFScratch *scratch) {
  memcpy(InvM_real, scratch->baseInvReal,
         scratch->invCount*sizeof(double));
  memcpy(PfM_real, scratch->basePfReal,
         scratch->pfCount*sizeof(double));
  memcpy(SlaterElmBF, scratch->baseSlaterComplex,
         scratch->slaterCount*sizeof(double complex));
  memcpy(SlaterElmBF_real, scratch->baseSlaterReal,
         scratch->slaterCount*sizeof(double));
  lsbfRestoreEta(scratch);
}

static int lsbfApplyHop(LSLanczosBFScratch *scratch,
                        int ri, int rj, int s,
                        const int *eleProjCnt) {
  const int source = rj + s*Nsite;
  const int target = ri + s*Nsite;
  const int electron = scratch->eleCfg[source];
  if(electron < 0 || electron >= Ne) return 1;
  scratch->eleIdx[electron + s*Ne] = ri;
  scratch->eleCfg[source] = -1;
  scratch->eleCfg[target] = electron;
  scratch->eleNum[source] = 0;
  scratch->eleNum[target] = 1;
  UpdateProjCnt(rj, ri, s, scratch->projCnt,
                eleProjCnt, scratch->eleNum);
  MakeProjBFCnt(scratch->bfCnt, scratch->eleNum);
  return 0;
}

static void lsbfCopyComplexSlaterToReal(void) {
  size_t idx;
  const size_t count = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  for(idx=0; idx<count; idx++) SlaterElmBF_real[idx] = creal(SlaterElmBF[idx]);
}

static int lsbfOuterComplex(int ri, int rj, int s,
                            double complex h1, double complex ip,
                            const int *eleIdx, const int *eleCfg,
                            const int *eleNum, const int *eleProjCnt,
                            const int *eleProjBFCnt,
                            LSLanczosBFScratch *scratch,
                            double complex *value) {
  const int target = ri + s*Nsite;
  const int source = rj + s*Nsite;
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

  lsbfCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt);
  StoreSlaterElmBF_fcmp(scratch->gfSlaterComplex);
  ratio = GreenFunc1BF(ri, rj, s, ip, scratch->gfSlaterComplex,
                       scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
                       eleProjCnt, scratch->gfProjCnt, eleProjBFCnt,
                       scratch->gfBFCnt, scratch->gfComplex);
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
    lsbfCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt);
    if(lsbfApplyHop(scratch, ri, rj, s, eleProjCnt) != 0) return 1;
    projRatio = ProjRatio(scratch->projCnt, eleProjCnt);
    MakeSlaterElmBF_fcmp(scratch->eleNum, scratch->bfCnt);
    info = LSLanczosTestConsumeRebuildFailure()
         ? 1
         : CalculateMAll_BF_fcmp(scratch->eleIdx, 0, NQPFull);
    if(info != 0) {
      lsbfRestoreComplex(scratch);
      return LSLANCZOS_NUMERIC_REJECT;
    }
    ipNew = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
    energy = CalculateHamiltonianBF_fcmp(
        ipNew, scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        scratch->projCnt, scratch->bfCnt);
    result = energy*conj(projRatio*ipNew/ip);
    lsbfRestoreComplex(scratch);
    *value = result;
    return 0;
  }

  lsbfCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt);
  if(lsbfApplyHop(scratch, ri, rj, s, eleProjCnt) != 0) return 1;
  result = ratio*CalculateHamiltonian0(scratch->eleNum);
  for(idx=0; idx<NTransfer; idx++) {
    double complex green;
    lsbfCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt);
    StoreSlaterElmBF_fcmp(scratch->gfSlaterComplex);
    green = GreenFunc2BF(
        Transfer[idx][0], Transfer[idx][2], ri, rj,
        Transfer[idx][3], s, ip, scratch->gfSlaterComplex,
        scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        eleProjCnt, scratch->gfProjCnt, eleProjBFCnt,
        scratch->gfBFCnt, scratch->gfComplex);
    result -= ParaTransfer[idx]*green;
  }
  *value = result;
  return 0;
}

static int lsbfOuterReal(int ri, int rj, int s,
                         double h1, double ip,
                         const int *eleIdx, const int *eleCfg,
                         const int *eleNum, const int *eleProjCnt,
                         const int *eleProjBFCnt,
                         LSLanczosBFScratch *scratch, double *value) {
  const int target = ri + s*Nsite;
  const int source = rj + s*Nsite;
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

  lsbfCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt);
  StoreSlaterElmBF_real(scratch->gfSlaterReal);
  ratio = GreenFunc1BF_real(
      ri, rj, s, ip, scratch->gfSlaterReal,
      scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
      eleProjCnt, scratch->gfProjCnt, eleProjBFCnt,
      scratch->gfBFCnt, scratch->gfReal);
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
    lsbfCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt);
    if(lsbfApplyHop(scratch, ri, rj, s, eleProjCnt) != 0) return 1;
    projRatio = ProjRatio(scratch->projCnt, eleProjCnt);
    MakeSlaterElmBF_fcmp(scratch->eleNum, scratch->bfCnt);
    lsbfCopyComplexSlaterToReal();
    info = LSLanczosTestConsumeRebuildFailure()
         ? 1
         : CalculateMAll_BF_real(scratch->eleIdx, 0, NQPFull);
    if(info != 0) {
      lsbfRestoreReal(scratch);
      return LSLANCZOS_NUMERIC_REJECT;
    }
    ipNew = CalculateIP_real(PfM_real, 0, NQPFull, MPI_COMM_SELF);
    energy = CalculateHamiltonianBF_real(
        ipNew, scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        scratch->projCnt, scratch->bfCnt);
    result = energy*projRatio*ipNew/ip;
    lsbfRestoreReal(scratch);
    *value = result;
    return 0;
  }

  lsbfCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt);
  if(lsbfApplyHop(scratch, ri, rj, s, eleProjCnt) != 0) return 1;
  result = ratio*CalculateHamiltonian0_real(scratch->eleNum);
  for(idx=0; idx<NTransfer; idx++) {
    double green;
    lsbfCopyState(scratch, eleIdx, eleCfg, eleNum, eleProjCnt);
    StoreSlaterElmBF_real(scratch->gfSlaterReal);
    green = GreenFunc2BF_real(
        Transfer[idx][0], Transfer[idx][2], ri, rj,
        Transfer[idx][3], s, ip, scratch->gfSlaterReal,
        scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        eleProjCnt, scratch->gfProjCnt, eleProjBFCnt,
        scratch->gfBFCnt, scratch->gfReal);
    result -= creal(ParaTransfer[idx])*green;
  }
  *value = result;
  return 0;
}

int LSLocalQBF(const double complex h1, const double complex ip,
               const int *eleIdx, const int *eleCfg, const int *eleNum,
               const int *eleProjCnt, const int *eleProjBFCnt,
               LSLanczosBFScratch *scratch, double complex *lslq) {
  double complex h2;
  int status;
  int idx;
  if(scratch == NULL || !scratch->initialized || scratch->useReal ||
     lslq == NULL) return 1;

  memcpy(scratch->baseInvComplex, InvM,
         scratch->invCount*sizeof(double complex));
  memcpy(scratch->basePfComplex, PfM,
         scratch->pfCount*sizeof(double complex));
  memcpy(scratch->baseSlaterComplex, SlaterElmBF,
         scratch->slaterCount*sizeof(double complex));
  lsbfCopyEtaTo(scratch->baseEta, scratch->baseEtaFlag);
  lsbfCaptureState(scratch, eleIdx, eleCfg, eleNum,
                   eleProjCnt, eleProjBFCnt);

  h2 = h1*CalculateHamiltonian0(eleNum);
  for(idx=0; idx<NTransfer; idx++) {
    double complex outer;
    status = lsbfOuterComplex(Transfer[idx][0], Transfer[idx][2],
                              Transfer[idx][3], h1, ip,
                              eleIdx, eleCfg, eleNum, eleProjCnt,
                              eleProjBFCnt, scratch, &outer);
    if(status != LSLANCZOS_OK) {
      lsbfRestoreComplex(scratch);
      if(lsbfVerifyCommonState(scratch, eleIdx, eleCfg, eleNum,
                               eleProjCnt, eleProjBFCnt) != 0 ||
         (scratch->stateCheck &&
          (memcmp(scratch->baseInvComplex, InvM,
                  scratch->invCount*sizeof(double complex)) != 0 ||
           memcmp(scratch->basePfComplex, PfM,
                  scratch->pfCount*sizeof(double complex)) != 0))) {
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
  if(lsbfVerifyCommonState(scratch, eleIdx, eleCfg, eleNum,
                           eleProjCnt, eleProjBFCnt) != 0 ||
     (scratch->stateCheck &&
      (memcmp(scratch->baseInvComplex, InvM,
              scratch->invCount*sizeof(double complex)) != 0 ||
       memcmp(scratch->basePfComplex, PfM,
              scratch->pfCount*sizeof(double complex)) != 0))) {
    return LSLANCZOS_STATE_FAILURE;
  }
  return LSLANCZOS_OK;
}

int LSLocalQBF_real(const double h1, const double ip,
                    const int *eleIdx, const int *eleCfg, const int *eleNum,
                    const int *eleProjCnt, const int *eleProjBFCnt,
                    LSLanczosBFScratch *scratch, double *lslq) {
  double h2;
  int status;
  int idx;
  if(scratch == NULL || !scratch->initialized || !scratch->useReal ||
     lslq == NULL) return 1;

  memcpy(scratch->baseInvReal, InvM_real,
         scratch->invCount*sizeof(double));
  memcpy(scratch->basePfReal, PfM_real,
         scratch->pfCount*sizeof(double));
  memcpy(scratch->baseSlaterComplex, SlaterElmBF,
         scratch->slaterCount*sizeof(double complex));
  memcpy(scratch->baseSlaterReal, SlaterElmBF_real,
         scratch->slaterCount*sizeof(double));
  lsbfCopyEtaTo(scratch->baseEta, scratch->baseEtaFlag);
  lsbfCaptureState(scratch, eleIdx, eleCfg, eleNum,
                   eleProjCnt, eleProjBFCnt);

  h2 = h1*CalculateHamiltonian0_real(eleNum);
  for(idx=0; idx<NTransfer; idx++) {
    double outer;
    status = lsbfOuterReal(Transfer[idx][0], Transfer[idx][2],
                           Transfer[idx][3], h1, ip,
                           eleIdx, eleCfg, eleNum, eleProjCnt,
                           eleProjBFCnt, scratch, &outer);
    if(status != LSLANCZOS_OK) {
      lsbfRestoreReal(scratch);
      if(lsbfVerifyCommonState(scratch, eleIdx, eleCfg, eleNum,
                               eleProjCnt, eleProjBFCnt) != 0 ||
         (scratch->stateCheck &&
          (memcmp(scratch->baseSlaterReal, SlaterElmBF_real,
                  scratch->slaterCount*sizeof(double)) != 0 ||
           memcmp(scratch->baseInvReal, InvM_real,
                  scratch->invCount*sizeof(double)) != 0 ||
           memcmp(scratch->basePfReal, PfM_real,
                  scratch->pfCount*sizeof(double)) != 0))) {
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
  if(lsbfVerifyCommonState(scratch, eleIdx, eleCfg, eleNum,
                           eleProjCnt, eleProjBFCnt) != 0 ||
     (scratch->stateCheck &&
      (memcmp(scratch->baseSlaterReal, SlaterElmBF_real,
              scratch->slaterCount*sizeof(double)) != 0 ||
       memcmp(scratch->baseInvReal, InvM_real,
              scratch->invCount*sizeof(double)) != 0 ||
       memcmp(scratch->basePfReal, PfM_real,
              scratch->pfCount*sizeof(double)) != 0))) {
    return LSLANCZOS_STATE_FAILURE;
  }
  return LSLANCZOS_OK;
}
