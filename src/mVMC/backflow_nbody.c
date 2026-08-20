#include "./include/backflow_nbody.h"

#include "./include/global.h"
#include "./include/locgrn.h"
#include "./include/locgrn_fsz.h"
#include "./include/matrix.h"
#include "./include/nbody_operator.h"
#include "./include/projection.h"
#include "./include/qp.h"
#include "./include/sector_projection.h"
#include "./include/slater.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int *eleIdx;
  int *eleCfg;
  int *eleNum;
  int *eleProjCnt;
  int *eleSpn;
  int *eleProjBFCnt;
  double complex *slater;
  double *slaterReal;
  double complex *invM;
  double complex *pfM;
  double complex *etaValues;
  int *etaFlags;
  double complex *baseInv;
  size_t nsizeCount;
  size_t nsite2Count;
  size_t nprojCount;
  size_t projBFCount;
  size_t slaterCount;
  size_t invCount;
  size_t nqpCount;
  size_t etaCount;
} BFNBodyStateSnapshot;

static BFNBodyResult BFNBodyResultValue(BFNBodyStatus status,
                                        BFNBodyStage stage, int detail,
                                        int reducedOrder,
                                        double complex value);

const char *BFNBodyStageName(BFNBodyStage stage) {
  switch(stage) {
  case BF_NBODY_STAGE_NONE: return "none";
  case BF_NBODY_STAGE_REDUCE: return "reduce";
  case BF_NBODY_STAGE_DISPATCH: return "dispatch";
  case BF_NBODY_STAGE_CANDIDATE: return "candidate";
  case BF_NBODY_STAGE_PROJECTION: return "projection";
  case BF_NBODY_STAGE_SLATER: return "slater";
  case BF_NBODY_STAGE_PFAFFIAN: return "pfaffian";
  case BF_NBODY_STAGE_RATIO: return "ratio";
  case BF_NBODY_STAGE_WORKSPACE: return "workspace";
  }
  return "unknown";
}

const char *BFNBodyDetailName(int detail) {
  switch(detail) {
  case BF_NBODY_DETAIL_NONE: return "none";
  case BF_NBODY_DETAIL_REDUCER: return "reducer";
  case BF_NBODY_DETAIL_SPIN_CHANGE: return "spin-change";
  case BF_NBODY_DETAIL_BAD_ELECTRON_LABEL: return "bad-electron-label";
  case BF_NBODY_DETAIL_LEGACY_STATE_RESTORE: return "legacy-state-restore";
  case BF_NBODY_DETAIL_INJECTED: return "injected";
  case BF_NBODY_DETAIL_STATE_CHANGED: return "state-changed";
  case BF_NBODY_DETAIL_BASE_STATE_STALE: return "base-state-stale";
  case BF_NBODY_DETAIL_STATE_SNAPSHOT: return "state-snapshot";
  }
  return "unknown";
}

static int BFStateByteCount(size_t count, size_t width, size_t *bytes) {
  if(bytes == NULL || (width != 0 && count > SIZE_MAX/width)) return -1;
  *bytes = count*width;
  return 0;
}

static int BFStateSizeMul(size_t left, size_t right, size_t *value) {
  if(value == NULL || (right != 0 && left > SIZE_MAX/right)) return -1;
  *value = left*right;
  return 0;
}

static void *BFStateAlloc(size_t count, size_t width) {
  size_t bytes;
  if(count == 0 || BFStateByteCount(count, width, &bytes) != 0) {
    return NULL;
  }
  return malloc(bytes);
}

static void *BFStateCopy(const void *source, size_t count, size_t width) {
  size_t bytes;
  void *copy;
  if(count == 0) return NULL;
  if(source == NULL || BFStateByteCount(count, width, &bytes) != 0) {
    return NULL;
  }
  copy = malloc(bytes);
  if(copy != NULL) memcpy(copy, source, bytes);
  return copy;
}

static void BFFreeNBodyStateSnapshot(BFNBodyStateSnapshot *snapshot) {
  if(snapshot == NULL) return;
  free(snapshot->eleIdx);
  free(snapshot->eleCfg);
  free(snapshot->eleNum);
  free(snapshot->eleProjCnt);
  free(snapshot->eleSpn);
  free(snapshot->eleProjBFCnt);
  free(snapshot->slater);
  free(snapshot->slaterReal);
  free(snapshot->invM);
  free(snapshot->pfM);
  free(snapshot->etaValues);
  free(snapshot->etaFlags);
  free(snapshot->baseInv);
  memset(snapshot, 0, sizeof(*snapshot));
}

static int BFInitNBodyStateSnapshot(
    BFNBodyStateSnapshot *snapshot, int useFsz,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt) {
  size_t siteRangeCount;
  size_t nsite2Square;
  size_t nsizeSquare;
  int site;

  if(snapshot == NULL || eleIdx == NULL || eleCfg == NULL || eleNum == NULL
     || Nsize <= 0 || Nsite <= 0 || Nsite2 <= 0 || Nrange <= 0
     || NQPFull <= 0 || NProj < 0
     || (NProj > 0 && eleProjCnt == NULL)
     || (useFsz && eleSpn == NULL) || eleProjBFCnt == NULL
     || SlaterElmBF == NULL || SlaterElmBF_real == NULL
     || InvM == NULL || PfM == NULL || eta == NULL || etaFlag == NULL) {
    return -1;
  }
  memset(snapshot, 0, sizeof(*snapshot));
  snapshot->nsizeCount = (size_t)Nsize;
  snapshot->nsite2Count = (size_t)Nsite2;
  snapshot->nprojCount = (size_t)NProj;
  snapshot->nqpCount = (size_t)NQPFull;
  if(BFStateSizeMul((size_t)Nsite, (size_t)Nrange,
                    &siteRangeCount) != 0
     || BFStateSizeMul((size_t)16, siteRangeCount,
                       &snapshot->projBFCount) != 0
     || BFStateSizeMul((size_t)Nsite2, (size_t)Nsite2,
                       &nsite2Square) != 0
     || BFStateSizeMul((size_t)NQPFull, nsite2Square,
                       &snapshot->slaterCount) != 0
     || BFStateSizeMul((size_t)Nsize, (size_t)Nsize,
                       &nsizeSquare) != 0
     || BFStateSizeMul((size_t)NQPFull, nsizeSquare,
                       &snapshot->invCount) != 0
     || BFStateSizeMul((size_t)Nsite, (size_t)Nsite,
                       &snapshot->etaCount) != 0) {
    return -1;
  }
  snapshot->eleIdx =
      BFStateCopy(eleIdx, snapshot->nsizeCount, sizeof(int));
  snapshot->eleCfg =
      BFStateCopy(eleCfg, snapshot->nsite2Count, sizeof(int));
  snapshot->eleNum =
      BFStateCopy(eleNum, snapshot->nsite2Count, sizeof(int));
  if(NProj > 0) {
    snapshot->eleProjCnt =
        BFStateCopy(eleProjCnt, snapshot->nprojCount, sizeof(int));
  }
  if(useFsz) {
    snapshot->eleSpn =
        BFStateCopy(eleSpn, snapshot->nsizeCount, sizeof(int));
  }
  snapshot->eleProjBFCnt =
      BFStateCopy(eleProjBFCnt, snapshot->projBFCount, sizeof(int));
  snapshot->slater =
      BFStateCopy(
          SlaterElmBF, snapshot->slaterCount, sizeof(double complex));
  snapshot->slaterReal =
      BFStateCopy(SlaterElmBF_real, snapshot->slaterCount, sizeof(double));
  snapshot->invM =
      BFStateCopy(InvM, snapshot->invCount, sizeof(double complex));
  snapshot->pfM =
      BFStateCopy(PfM, snapshot->nqpCount, sizeof(double complex));
  snapshot->etaValues = (double complex *)BFStateAlloc(
      snapshot->etaCount, sizeof(double complex));
  snapshot->etaFlags =
      (int *)BFStateAlloc(snapshot->etaCount, sizeof(int));
  snapshot->baseInv = (double complex *)BFStateAlloc(
      snapshot->invCount, sizeof(double complex));
  if(snapshot->eleIdx == NULL || snapshot->eleCfg == NULL
     || snapshot->eleNum == NULL
     || (NProj > 0 && snapshot->eleProjCnt == NULL)
     || (useFsz && snapshot->eleSpn == NULL)
     || snapshot->eleProjBFCnt == NULL || snapshot->slater == NULL
     || snapshot->slaterReal == NULL || snapshot->invM == NULL
     || snapshot->pfM == NULL || snapshot->etaValues == NULL
     || snapshot->etaFlags == NULL || snapshot->baseInv == NULL) {
    BFFreeNBodyStateSnapshot(snapshot);
    return -1;
  }
  for(site=0;site<Nsite;site++) {
    if(eta[site] == NULL || etaFlag[site] == NULL) {
      BFFreeNBodyStateSnapshot(snapshot);
      return -1;
    }
    memcpy(snapshot->etaValues+(size_t)site*(size_t)Nsite, eta[site],
           (size_t)Nsite*sizeof(double complex));
    memcpy(snapshot->etaFlags+(size_t)site*(size_t)Nsite, etaFlag[site],
           (size_t)Nsite*sizeof(int));
  }
  return 0;
}

static int BFStateComplexArraysNear(const double complex *left,
                                    const double complex *right,
                                    size_t count) {
  const double tolerance = 1.0e-11;
  size_t idx;
  if(left == NULL || right == NULL) return 0;
  for(idx=0;idx<count;idx++) {
    const double scale = 1.0+cabs(left[idx])+cabs(right[idx]);
    if(!isfinite(creal(left[idx])) || !isfinite(cimag(left[idx]))
       || !isfinite(creal(right[idx])) || !isfinite(cimag(right[idx]))
       || cabs(left[idx]-right[idx]) > tolerance*scale) {
      return 0;
    }
  }
  return 1;
}

static int BFNBodyBaseStateCurrent(
    const BFNBodyStateSnapshot *snapshot, int useFsz,
    const int *eleIdx, const int *eleNum, const int *eleSpn,
    const int *eleProjBFCnt, BFNBodyScratch *scratch) {
  BF_FSZ_MAllResult matrixResult;
  int particle;
  if(snapshot == NULL || eleIdx == NULL || eleNum == NULL
     || eleProjBFCnt == NULL || scratch == NULL
     || !ValidateBFNBodyScratch(scratch, 1, useFsz)) {
    return 0;
  }
  if(useFsz) {
    if(eleSpn == NULL) return 0;
    memcpy(scratch->eleSpn, eleSpn, snapshot->nsizeCount*sizeof(int));
    MakeSlaterElmBF_fsz_to_serial(
        scratch->slater, eleNum, eleProjBFCnt);
  } else {
    for(particle=0;particle<Nsize;particle++) {
      scratch->eleSpn[particle] = particle < Ne ? 0 : 1;
    }
    MakeSlaterElmBF_fcmp_to_serial(
        scratch->slater, eleNum, eleProjBFCnt);
  }
  if(memcmp(scratch->slater, snapshot->slater,
            snapshot->slaterCount*sizeof(double complex)) != 0) {
    return 0;
  }
  matrixResult = CalculateMAll_BF_fsz_from_workspace(
      scratch->slater, eleIdx, scratch->eleSpn, 0, NQPFull,
      scratch->candidatePf, snapshot->baseInv,
      snapshot->invCount/snapshot->nqpCount, scratch->pfBufM,
      scratch->pfIWork, scratch->pfWork, LapackLWork,
      scratch->pfRWork);
  return matrixResult.status == BF_FSZ_MALL_OK
      && BFStateComplexArraysNear(
             scratch->candidatePf, snapshot->pfM, snapshot->nqpCount)
      && BFStateComplexArraysNear(
             snapshot->baseInv, snapshot->invM, snapshot->invCount);
}

static int BFNBodyStateChanged(
    const BFNBodyStateSnapshot *snapshot, int useFsz,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt) {
  int site;

  if(snapshot == NULL
     || memcmp(eleIdx, snapshot->eleIdx,
               snapshot->nsizeCount*sizeof(int)) != 0
     || memcmp(eleCfg, snapshot->eleCfg,
               snapshot->nsite2Count*sizeof(int)) != 0
     || memcmp(eleNum, snapshot->eleNum,
               snapshot->nsite2Count*sizeof(int)) != 0
     || (NProj > 0
         && memcmp(eleProjCnt, snapshot->eleProjCnt,
                   snapshot->nprojCount*sizeof(int)) != 0)
     || (useFsz
         && memcmp(eleSpn, snapshot->eleSpn,
                   snapshot->nsizeCount*sizeof(int)) != 0)
     || memcmp(eleProjBFCnt, snapshot->eleProjBFCnt,
               snapshot->projBFCount*sizeof(int)) != 0
     || memcmp(SlaterElmBF, snapshot->slater,
               snapshot->slaterCount*sizeof(double complex)) != 0
     || memcmp(SlaterElmBF_real, snapshot->slaterReal,
               snapshot->slaterCount*sizeof(double)) != 0
     || memcmp(InvM, snapshot->invM,
               snapshot->invCount*sizeof(double complex)) != 0
     || memcmp(PfM, snapshot->pfM,
               snapshot->nqpCount*sizeof(double complex)) != 0) {
    return 1;
  }
  for(site=0;site<Nsite;site++) {
    if(memcmp(eta[site],
              snapshot->etaValues+(size_t)site*(size_t)Nsite,
              (size_t)Nsite*sizeof(double complex)) != 0
       || memcmp(etaFlag[site],
                 snapshot->etaFlags+(size_t)site*(size_t)Nsite,
                 (size_t)Nsite*sizeof(int)) != 0) {
      return 1;
    }
  }
  return 0;
}

static BFNBodyResult BFInjectedNBodyResult(int injectStage,
                                           int reducedOrder) {
  BFNBodyStatus status = BF_NBODY_INVALID_ARGUMENT;
  BFNBodyStage stage = BF_NBODY_STAGE_CANDIDATE;
  if(injectStage == BF_NBODY_INJECT_WORKSPACE) {
    status = BF_NBODY_WORKSPACE_ERROR;
    stage = BF_NBODY_STAGE_WORKSPACE;
  } else if(injectStage == BF_NBODY_INJECT_PFAFFIAN) {
    status = BF_NBODY_PFAFFIAN_ERROR;
    stage = BF_NBODY_STAGE_PFAFFIAN;
  }
  return BFNBodyResultValue(
      status, stage, BF_NBODY_DETAIL_INJECTED,
      reducedOrder, 0.0+0.0*I);
}

static int BFShouldInjectNBodyFailure(const BFNBodyScratch *scratch,
                                      int injectStage) {
  return BFNBodyInjectStage == injectStage
      && scratch != NULL && scratch->termIndex == BFNBodyInjectTerm;
}

int GetBFNBodyScratchSizes(int maxOrder, int useFsz,
                           BFNBodyScratchSizes *sizes) {
  BFNBodyScratchDimensions dimensions;

  if(sizes == NULL || maxOrder < 1 || (useFsz != 0 && useFsz != 1)) {
    return BF_NBODY_WORKSPACE_ERROR;
  }

  memset(&dimensions, 0, sizeof(dimensions));
  dimensions.maxOrder = maxOrder;
  dimensions.useFsz = useFsz;
  dimensions.nsize = Nsize;
  dimensions.nsite2 = Nsite2;
  dimensions.nproj = NProj;
  dimensions.nsite = Nsite;
  dimensions.nrange = Nrange;
  dimensions.nqpFull = NQPFull;
  dimensions.lapackLWork = LapackLWork;

  if(useFsz
     && (GetGreenFuncBF_fsz_buffer_work_size(
             &dimensions.bfFszGreenBufferCount,
             &dimensions.pfUpdateIntCount,
             &dimensions.pfUpdateDoubleCount) != 0
         || GetSlaterElmBF_fsz_hop_int_work_size(
             &dimensions.bfHopIntSize) != 0)) {
    return BF_NBODY_WORKSPACE_ERROR;
  }

  return ComputeBFNBodyScratchSizes(&dimensions, sizes);
}

static BFNBodyResult BFNBodyResultValue(BFNBodyStatus status,
                                        BFNBodyStage stage, int detail,
                                        int reducedOrder,
                                        double complex value) {
  BFNBodyResult result = {
      status, stage, detail, reducedOrder, value
  };
  return result;
}

int ValidateBFNBodyScratch(const BFNBodyScratch *scratch, int n,
                           int useFsz) {
  BFNBodyScratchSizes needed;

  if(scratch == NULL || n < 1 || scratch->maxOrder < n
     || scratch->useFsz != useFsz
     || scratch->sizes.maxOrder != scratch->maxOrder
     || scratch->sizes.useFsz != scratch->useFsz
     || scratch->inputRsi == NULL || scratch->inputRsj == NULL
     || scratch->rsi == NULL || scratch->rsj == NULL
     || scratch->moved == NULL || scratch->eleIdx == NULL
     || scratch->eleCfg == NULL || scratch->eleNum == NULL
     || scratch->eleSpn == NULL || scratch->projBFCnt == NULL
     || scratch->pfIWork == NULL || scratch->greenBuffer == NULL
     || scratch->slater == NULL || scratch->candidatePf == NULL
     || scratch->pfBufM == NULL || scratch->pfWork == NULL
     || scratch->pfRWork == NULL
     || (useFsz
         && (scratch->affected == NULL || scratch->hopIntWork == NULL))
     || (NProj > 0 && scratch->projCnt == NULL)
     || GetBFNBodyScratchSizes(scratch->maxOrder, useFsz, &needed)
          != BF_NBODY_OK
     || (!useFsz && (Ne <= 0 || Ne > INT_MAX/2 || Nsize != 2*Ne))
     || scratch->sizes.inputRsiCount < needed.inputRsiCount
     || scratch->sizes.inputRsjCount < needed.inputRsjCount
     || scratch->sizes.rsiCount < needed.rsiCount
     || scratch->sizes.rsjCount < needed.rsjCount
     || scratch->sizes.movedCount < needed.movedCount
     || scratch->sizes.eleIdxCount < needed.eleIdxCount
     || scratch->sizes.eleCfgCount < needed.eleCfgCount
     || scratch->sizes.eleNumCount < needed.eleNumCount
     || scratch->sizes.eleSpnCount < needed.eleSpnCount
     || scratch->sizes.projCntCount < needed.projCntCount
     || scratch->sizes.projBFCntCount < needed.projBFCntCount
     || scratch->sizes.pfIWorkCount < needed.pfIWorkCount
     || scratch->sizes.affectedCount < needed.affectedCount
     || scratch->sizes.hopIntWorkCount < needed.hopIntWorkCount
     || scratch->sizes.greenBufferCount < needed.greenBufferCount
     || scratch->sizes.slaterCount < needed.slaterCount
     || scratch->sizes.candidatePfCount < needed.candidatePfCount
     || scratch->sizes.pfBufMCount < needed.pfBufMCount
     || scratch->sizes.pfWorkCount < needed.pfWorkCount
     || scratch->sizes.pfRWorkCount < needed.pfRWorkCount) {
    return 0;
  }
  return 1;
}

static int BFCallerStateAliasesNBodyScratch(
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn, const int *eleProjBFCnt,
    const BFNBodyScratch *scratch) {
  const int *callerState[6] = {
      eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn, eleProjBFCnt
  };
  const int *scratchState[6] = {
      scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
      scratch->projCnt, scratch->eleSpn, scratch->projBFCnt
  };
  int callerIdx;
  int scratchIdx;

  for(callerIdx=0;callerIdx<6;callerIdx++) {
    if(callerState[callerIdx] == NULL) continue;
    for(scratchIdx=0;scratchIdx<6;scratchIdx++) {
      if(callerState[callerIdx] == scratchState[scratchIdx]) return 1;
    }
  }
  return 0;
}

static BFNBodyResult BFDispatchReducedNBody(
    const NBodyReduction *reduction, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleProjBFCnt,
    BFNBodyScratch *scratch) {
  const size_t slaterCount =
      (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const int target = scratch->rsi[0];
  const int source = scratch->rsj[0];
  const int spin = source/Nsite;
  double complex value;
  int greenStatus = BF_PF_OK;

  if(reduction->order != 1) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_DISPATCH,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->projCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
  memcpy(scratch->projBFCnt, eleProjBFCnt,
         (size_t)16*(size_t)Nsite*(size_t)Nrange*sizeof(int));
  memcpy(scratch->slater, SlaterElmBF,
         slaterCount*sizeof(double complex));

  value = GreenFunc1BF(
      target%Nsite, source%Nsite, spin, ip, scratch->slater,
      scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
      eleProjCnt, scratch->projCnt, eleProjBFCnt,
      scratch->projBFCnt, scratch->greenBuffer, &greenStatus);

  if(greenStatus != BF_PF_OK) {
    return BFNBodyResultValue(
        BF_NBODY_PFAFFIAN_ERROR, BF_NBODY_STAGE_PFAFFIAN,
        greenStatus, reduction->order, 0.0+0.0*I);
  }

  if(memcmp(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int)) != 0
     || memcmp(scratch->eleCfg, eleCfg,
               (size_t)Nsite2*sizeof(int)) != 0
     || memcmp(scratch->eleNum, eleNum,
               (size_t)Nsite2*sizeof(int)) != 0
     || memcmp(scratch->slater, SlaterElmBF,
               slaterCount*sizeof(double complex)) != 0) {
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_DISPATCH,
        BF_NBODY_DETAIL_LEGACY_STATE_RESTORE, reduction->order, 0.0+0.0*I);
  }

  value *= (double)reduction->sign;
  if(!isfinite(creal(value)) || !isfinite(cimag(value))) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_DISPATCH, BF_NBODY_DETAIL_NONE,
        reduction->order, 0.0+0.0*I);
  }
  if(cabs(value) == 0.0) {
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_RATIO, BF_NBODY_DETAIL_NONE,
        reduction->order, 0.0+0.0*I);
  }
  return BFNBodyResultValue(
      BF_NBODY_OK, BF_NBODY_STAGE_NONE, BF_NBODY_DETAIL_NONE,
      reduction->order, value);
}

static BFNBodyResult BFRebuildReducedNBody(
    const NBodyReduction *reduction, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleProjBFCnt,
    BFNBodyScratch *scratch) {
  double complex ipNew;
  double complex ratio;
  double projRatio;
  int pfDetail = 0;
  int pfStatus;
  int k;

  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));

  for(k=0;k<reduction->order;k++) {
    const int target = scratch->rsi[k];
    const int source = scratch->rsj[k];
    const int sourceSpin = source/Nsite;
    const int targetSpin = target/Nsite;
    const int localLabel = eleCfg[source];
    const int particle = localLabel+sourceSpin*Ne;
    int l;

    if(sourceSpin != targetSpin) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
          BF_NBODY_DETAIL_SPIN_CHANGE, reduction->order, 0.0+0.0*I);
    }
    if(localLabel < 0 || localLabel >= Ne || particle < 0
       || particle >= Nsize || eleNum[source] != 1
       || eleNum[target] != 0 || eleIdx[particle] != source%Nsite) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
          BF_NBODY_DETAIL_BAD_ELECTRON_LABEL, reduction->order,
          0.0+0.0*I);
    }
    for(l=0;l<k;l++) {
      if(scratch->moved[l] == particle
         || scratch->rsi[l] == target || scratch->rsj[l] == source) {
        return BFNBodyResultValue(
            BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
            BF_NBODY_DETAIL_BAD_ELECTRON_LABEL, reduction->order,
            0.0+0.0*I);
      }
    }
    scratch->moved[k] = particle;
  }

  for(k=0;k<reduction->order;k++) {
    const int source = scratch->rsj[k];
    scratch->eleCfg[source] = -1;
    scratch->eleNum[source] = 0;
  }
  for(k=0;k<reduction->order;k++) {
    const int target = scratch->rsi[k];
    const int particle = scratch->moved[k];
    scratch->eleCfg[target] = particle%Ne;
    scratch->eleNum[target] = 1;
    scratch->eleIdx[particle] = target%Nsite;
  }
  for(k=0;k<Nsize;k++) {
    scratch->eleSpn[k] = k < Ne ? 0 : 1;
  }

  MakeProjCnt(scratch->projCnt, scratch->eleNum);
  MakeProjBFCnt(scratch->projBFCnt, scratch->eleNum);
  if(!IsSectorStateAllowed(scratch->eleNum)) {
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_PROJECTION,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }

  MakeSlaterElmBF_fcmp_to_serial(
      scratch->slater, scratch->eleNum, scratch->projBFCnt);
  if(BFShouldInjectNBodyFailure(
         scratch, BF_NBODY_INJECT_PFAFFIAN)) {
    return BFInjectedNBodyResult(
        BF_NBODY_INJECT_PFAFFIAN, reduction->order);
  }
  pfStatus = CalculatePfM_BF_from_workspace(
      scratch->slater, scratch->eleIdx, scratch->eleSpn, 0, NQPFull,
      scratch->candidatePf, &pfDetail, scratch->pfBufM,
      scratch->pfIWork, scratch->pfWork, LapackLWork, scratch->pfRWork);
  if(pfStatus != BF_PF_OK) {
    const BFNBodyStatus status =
        pfStatus == BF_PF_NONFINITE
        ? BF_NBODY_NONFINITE : BF_NBODY_PFAFFIAN_ERROR;
    return BFNBodyResultValue(
        status, BF_NBODY_STAGE_PFAFFIAN,
        pfDetail != 0 ? pfDetail : pfStatus,
        reduction->order, 0.0+0.0*I);
  }

  ipNew = CalculateIP_fcmp(
      scratch->candidatePf, 0, NQPFull, MPI_COMM_SELF);
  if(!isfinite(creal(ipNew)) || !isfinite(cimag(ipNew))) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_RATIO, BF_NBODY_DETAIL_NONE,
        reduction->order, 0.0+0.0*I);
  }
  if(cabs(ipNew) == 0.0) {
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_RATIO,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }

  projRatio = ProjRatio(scratch->projCnt, eleProjCnt);
  if(!isfinite(projRatio)) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_PROJECTION,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  ratio = (double)reduction->sign*conj(projRatio*ipNew/ip);
  if(!isfinite(creal(ratio)) || !isfinite(cimag(ratio))) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_RATIO, BF_NBODY_DETAIL_NONE,
        reduction->order, 0.0+0.0*I);
  }
  return BFNBodyResultValue(
      BF_NBODY_OK, BF_NBODY_STAGE_NONE, BF_NBODY_DETAIL_NONE,
      reduction->order, ratio);
}

static BFNBodyResult GreenFuncNBFImpl(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleProjBFCnt,
    BFNBodyScratch *scratch) {
  NBodyReduction reduction;
  int reduceStatus;
  int k;

  if(BFShouldInjectNBodyFailure(
         scratch, BF_NBODY_INJECT_WORKSPACE)) {
    return BFInjectedNBodyResult(BF_NBODY_INJECT_WORKSPACE, 0);
  }
  if(!ValidateBFNBodyScratch(scratch, n, 0)) {
    if(n < 1) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_REDUCE,
          BF_NBODY_DETAIL_REDUCER, 0, 0.0+0.0*I);
    }
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_NONE,
        BF_NBODY_DETAIL_NONE, 0, 0.0+0.0*I);
  }
  if(BFCallerStateAliasesNBodyScratch(
         eleIdx, eleCfg, eleNum, eleProjCnt, NULL, eleProjBFCnt,
         scratch)) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
        BF_NBODY_DETAIL_NONE, 0, 0.0+0.0*I);
  }
  reduceStatus = ReduceNBodyTerm(
      n, rsi, rsj, eleNum, Nsite2, scratch->rsi, scratch->rsj,
      scratch->maxOrder, &reduction);
  if(reduceStatus != 0 || reduction.kind == NBODY_REDUCED_INVALID) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_REDUCE,
        reduceStatus != 0 ? reduceStatus : BF_NBODY_DETAIL_REDUCER,
        0, 0.0+0.0*I);
  }
  if(reduction.kind == NBODY_REDUCED_ZERO) {
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_REDUCE,
        BF_NBODY_DETAIL_NONE, 0, 0.0+0.0*I);
  }
  if(reduction.kind == NBODY_REDUCED_SCALAR) {
    return BFNBodyResultValue(
        BF_NBODY_OK, BF_NBODY_STAGE_REDUCE, BF_NBODY_DETAIL_NONE, 0,
        (double)reduction.sign+0.0*I);
  }
  if(eleIdx == NULL || eleCfg == NULL || eleNum == NULL
     || (NProj > 0 && eleProjCnt == NULL) || eleProjBFCnt == NULL) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
        BF_NBODY_DETAIL_NONE, reduction.order, 0.0+0.0*I);
  }
  if(!isfinite(creal(ip)) || !isfinite(cimag(ip))) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_RATIO, BF_NBODY_DETAIL_NONE,
        reduction.order, 0.0+0.0*I);
  }
  if(cabs(ip) == 0.0) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_RATIO,
        BF_NBODY_DETAIL_NONE, reduction.order, 0.0+0.0*I);
  }
  for(k=0;k<reduction.order;k++) {
    const int source = scratch->rsj[k];
    const int target = scratch->rsi[k];
    const int spin = source/Nsite;
    const int localLabel = eleCfg[source];
    const int particle = localLabel+spin*Ne;
    if(scratch->rsi[k]/Nsite != scratch->rsj[k]/Nsite) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT,
          reduction.order == 1 ? BF_NBODY_STAGE_DISPATCH
                               : BF_NBODY_STAGE_CANDIDATE,
          BF_NBODY_DETAIL_SPIN_CHANGE, reduction.order, 0.0+0.0*I);
    }
    if(localLabel < 0 || localLabel >= Ne || particle < 0
       || particle >= Nsize || eleNum[source] != 1
       || eleNum[target] != 0 || eleIdx[particle] != source%Nsite) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
          BF_NBODY_DETAIL_BAD_ELECTRON_LABEL, reduction.order,
          0.0+0.0*I);
    }
  }
  if(BFShouldInjectNBodyFailure(
         scratch, BF_NBODY_INJECT_CANDIDATE)) {
    return BFInjectedNBodyResult(
        BF_NBODY_INJECT_CANDIDATE, reduction.order);
  }
  if(reduction.order == 1
     && BFShouldInjectNBodyFailure(
            scratch, BF_NBODY_INJECT_PFAFFIAN)) {
    return BFInjectedNBodyResult(
        BF_NBODY_INJECT_PFAFFIAN, reduction.order);
  }
  /*
   * Keep the established one-body update. The legacy two-body BackFlow
   * update is not configuration-exact for every non-identity projection,
   * so effective N>=2 uses the checked out-of-place one-shot rebuild.
   */
  if(reduction.order == 1) {
    return BFDispatchReducedNBody(
        &reduction, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
        eleProjBFCnt, scratch);
  }
  return BFRebuildReducedNBody(
      &reduction, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
      eleProjBFCnt, scratch);
}

BFNBodyResult GreenFuncNBF(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleProjBFCnt,
    BFNBodyScratch *scratch) {
  BFNBodyStateSnapshot snapshot;
  BFNBodyResult result;
  int changed;
  if(!BFNBodyStateCheckEnabled) {
    return GreenFuncNBFImpl(
        n, rsi, rsj, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
        eleProjBFCnt, scratch);
  }
  if(BFInitNBodyStateSnapshot(
         &snapshot, 0, eleIdx, eleCfg, eleNum, eleProjCnt,
         NULL, eleProjBFCnt) != 0) {
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_WORKSPACE,
        BF_NBODY_DETAIL_STATE_SNAPSHOT, 0, 0.0+0.0*I);
  }
  if(!BFNBodyBaseStateCurrent(
         &snapshot, 0, eleIdx, eleNum, NULL, eleProjBFCnt, scratch)) {
    BFFreeNBodyStateSnapshot(&snapshot);
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_DISPATCH,
        BF_NBODY_DETAIL_BASE_STATE_STALE, 0, 0.0+0.0*I);
  }
  result = GreenFuncNBFImpl(
      n, rsi, rsj, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
      eleProjBFCnt, scratch);
  changed = BFNBodyStateChanged(
      &snapshot, 0, eleIdx, eleCfg, eleNum, eleProjCnt,
      NULL, eleProjBFCnt);
  BFFreeNBodyStateSnapshot(&snapshot);
#pragma omp atomic update
  BFNBodyStateCheckCount++;
  if(changed) {
#pragma omp atomic update
    BFNBodyStateCheckFailureCount++;
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR,
        result.stage == BF_NBODY_STAGE_NONE
        ? BF_NBODY_STAGE_CANDIDATE : result.stage,
        BF_NBODY_DETAIL_STATE_CHANGED, result.reducedOrder,
        0.0+0.0*I);
  }
  return result;
}

static BFNBodyResult BFDispatchReducedNBodyFSZ(
    const NBodyReduction *reduction, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt, BFNBodyScratch *scratch) {
  double complex value;
  int k;
  int l;

  if(reduction->order < 1 || reduction->order > 2) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_DISPATCH,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleSpn, eleSpn, (size_t)Nsize*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->projCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
  memcpy(scratch->projBFCnt, eleProjBFCnt,
         (size_t)16*(size_t)Nsite*(size_t)Nrange*sizeof(int));

  /*
   * A reduced hop product must have disjoint creation/annihilation
   * orbitals. Pin that reducer invariant before entering legacy kernels
   * whose identity shortcuts precede their general state checks.
   */
  for(k=0;k<reduction->order;k++) {
    if(scratch->rsi[k] == scratch->rsj[k]) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_DISPATCH,
          BF_NBODY_DETAIL_REDUCER, reduction->order, 0.0+0.0*I);
    }
    for(l=0;l<k;l++) {
      if(scratch->rsi[k] == scratch->rsi[l]
         || scratch->rsi[k] == scratch->rsj[l]
         || scratch->rsj[k] == scratch->rsi[l]
         || scratch->rsj[k] == scratch->rsj[l]) {
        return BFNBodyResultValue(
            BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_DISPATCH,
            BF_NBODY_DETAIL_REDUCER, reduction->order, 0.0+0.0*I);
      }
    }
  }

  /*
   * Classify a candidate rejected by t-J/doublon sector projection before
   * the one-/two-body legacy kernels turn it into an undifferentiated zero.
   */
  for(k=0;k<reduction->order;k++) {
    scratch->eleNum[scratch->rsj[k]] = 0;
  }
  for(k=0;k<reduction->order;k++) {
    scratch->eleNum[scratch->rsi[k]] = 1;
  }
  if(!IsSectorStateAllowed(scratch->eleNum)) {
    memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_PROJECTION,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));

  if(reduction->order == 1) {
    const int target = scratch->rsi[0];
    const int source = scratch->rsj[0];
    value = GreenFunc1BF_fsz2_workspace(
        target%Nsite, source%Nsite, target/Nsite, source/Nsite, ip,
        scratch->eleIdx, scratch->eleCfg, scratch->eleNum, eleProjCnt,
        scratch->eleSpn, scratch->projCnt, eleProjBFCnt,
        scratch->projBFCnt, scratch->greenBuffer, scratch->pfBufM,
        scratch->affected, scratch->hopIntWork,
        (int)scratch->sizes.hopIntWorkCount, scratch->pfIWork,
        scratch->pfWork, scratch->pfRWork);
  } else {
    const int target0 = scratch->rsi[0];
    const int source0 = scratch->rsj[0];
    const int target1 = scratch->rsi[1];
    const int source1 = scratch->rsj[1];
    value = GreenFunc2BF_fsz2(
        target0%Nsite, source0%Nsite, target1%Nsite, source1%Nsite,
        target0/Nsite, source0/Nsite, target1/Nsite, source1/Nsite,
        ip, scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        eleProjCnt, scratch->eleSpn, scratch->projCnt, eleProjBFCnt,
        scratch->projBFCnt, scratch->greenBuffer, scratch->affected,
        scratch->hopIntWork, (int)scratch->sizes.hopIntWorkCount,
        scratch->pfIWork, scratch->pfRWork, scratch->pfBufM,
        scratch->pfWork);
  }

  if(memcmp(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int)) != 0
     || memcmp(scratch->eleCfg, eleCfg,
               (size_t)Nsite2*sizeof(int)) != 0
     || memcmp(scratch->eleNum, eleNum,
               (size_t)Nsite2*sizeof(int)) != 0
     || memcmp(scratch->eleSpn, eleSpn,
               (size_t)Nsize*sizeof(int)) != 0) {
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_DISPATCH,
        BF_NBODY_DETAIL_LEGACY_STATE_RESTORE, reduction->order,
        0.0+0.0*I);
  }

  value *= (double)reduction->sign;
  if(!isfinite(creal(value)) || !isfinite(cimag(value))) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_DISPATCH,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  if(cabs(value) == 0.0) {
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_RATIO,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  return BFNBodyResultValue(
      BF_NBODY_OK, BF_NBODY_STAGE_NONE, BF_NBODY_DETAIL_NONE,
      reduction->order, value);
}

static BFNBodyResult BFRebuildReducedNBodyFSZ(
    const NBodyReduction *reduction, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    BFNBodyScratch *scratch) {
  double complex ipNew;
  double complex ratio;
  double projRatio;
  int pfDetail = 0;
  int pfStatus;
  int k;

  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleSpn, eleSpn, (size_t)Nsize*sizeof(int));

  for(k=0;k<reduction->order;k++) {
    const int target = scratch->rsi[k];
    const int source = scratch->rsj[k];
    const int particle = eleCfg[source];
    int l;

    if(particle < 0 || particle >= Nsize || eleNum[source] != 1
       || eleNum[target] != 0 || eleIdx[particle] != source%Nsite
       || eleSpn[particle] != source/Nsite) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
          BF_NBODY_DETAIL_BAD_ELECTRON_LABEL, reduction->order,
          0.0+0.0*I);
    }
    for(l=0;l<k;l++) {
      if(scratch->moved[l] == particle
         || scratch->rsi[l] == target || scratch->rsj[l] == source) {
        return BFNBodyResultValue(
            BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
            BF_NBODY_DETAIL_BAD_ELECTRON_LABEL, reduction->order,
            0.0+0.0*I);
      }
    }
    scratch->moved[k] = particle;
  }

  for(k=0;k<reduction->order;k++) {
    const int source = scratch->rsj[k];
    scratch->eleCfg[source] = -1;
    scratch->eleNum[source] = 0;
  }
  for(k=0;k<reduction->order;k++) {
    const int target = scratch->rsi[k];
    const int particle = scratch->moved[k];
    scratch->eleCfg[target] = particle;
    scratch->eleNum[target] = 1;
    scratch->eleIdx[particle] = target%Nsite;
    scratch->eleSpn[particle] = target/Nsite;
  }

  MakeProjCnt(scratch->projCnt, scratch->eleNum);
  MakeProjBFCnt(scratch->projBFCnt, scratch->eleNum);
  if(!IsSectorStateAllowed(scratch->eleNum)) {
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_PROJECTION,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }

  MakeSlaterElmBF_fsz_to_serial(
      scratch->slater, scratch->eleNum, scratch->projBFCnt);
  if(BFShouldInjectNBodyFailure(
         scratch, BF_NBODY_INJECT_PFAFFIAN)) {
    return BFInjectedNBodyResult(
        BF_NBODY_INJECT_PFAFFIAN, reduction->order);
  }
  pfStatus = CalculatePfM_BF_from_workspace(
      scratch->slater, scratch->eleIdx, scratch->eleSpn, 0, NQPFull,
      scratch->candidatePf, &pfDetail, scratch->pfBufM,
      scratch->pfIWork, scratch->pfWork, LapackLWork, scratch->pfRWork);
  if(pfStatus != BF_PF_OK) {
    const BFNBodyStatus status =
        pfStatus == BF_PF_NONFINITE
        ? BF_NBODY_NONFINITE : BF_NBODY_PFAFFIAN_ERROR;
    return BFNBodyResultValue(
        status, BF_NBODY_STAGE_PFAFFIAN,
        pfDetail != 0 ? pfDetail : pfStatus,
        reduction->order, 0.0+0.0*I);
  }

  ipNew = CalculateIP_fcmp(
      scratch->candidatePf, 0, NQPFull, MPI_COMM_SELF);
  if(!isfinite(creal(ipNew)) || !isfinite(cimag(ipNew))) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_RATIO,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  if(cabs(ipNew) == 0.0) {
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_RATIO,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }

  projRatio = ProjRatio(scratch->projCnt, eleProjCnt);
  if(!isfinite(projRatio)) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_PROJECTION,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  ratio = (double)reduction->sign*conj(projRatio*ipNew/ip);
  if(!isfinite(creal(ratio)) || !isfinite(cimag(ratio))) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_RATIO,
        BF_NBODY_DETAIL_NONE, reduction->order, 0.0+0.0*I);
  }
  return BFNBodyResultValue(
      BF_NBODY_OK, BF_NBODY_STAGE_NONE, BF_NBODY_DETAIL_NONE,
      reduction->order, ratio);
}

static BFNBodyResult GreenFuncNBF_fszImpl(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt, BFNBodyScratch *scratch) {
  NBodyReduction reduction;
  int reduceStatus;
  int k;

  if(BFShouldInjectNBodyFailure(
         scratch, BF_NBODY_INJECT_WORKSPACE)) {
    return BFInjectedNBodyResult(BF_NBODY_INJECT_WORKSPACE, 0);
  }
  if(!ValidateBFNBodyScratch(scratch, n, 1)) {
    if(n < 1) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_REDUCE,
          BF_NBODY_DETAIL_REDUCER, 0, 0.0+0.0*I);
    }
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_NONE,
        BF_NBODY_DETAIL_NONE, 0, 0.0+0.0*I);
  }
  if(BFCallerStateAliasesNBodyScratch(
         eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn, eleProjBFCnt,
         scratch)) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
        BF_NBODY_DETAIL_NONE, 0, 0.0+0.0*I);
  }
  reduceStatus = ReduceNBodyTerm(
      n, rsi, rsj, eleNum, Nsite2, scratch->rsi, scratch->rsj,
      scratch->maxOrder, &reduction);
  if(reduceStatus != 0 || reduction.kind == NBODY_REDUCED_INVALID) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_REDUCE,
        reduceStatus != 0 ? reduceStatus : BF_NBODY_DETAIL_REDUCER,
        0, 0.0+0.0*I);
  }
  if(reduction.kind == NBODY_REDUCED_ZERO) {
    return BFNBodyResultValue(
        BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_REDUCE,
        BF_NBODY_DETAIL_NONE, 0, 0.0+0.0*I);
  }
  if(reduction.kind == NBODY_REDUCED_SCALAR) {
    return BFNBodyResultValue(
        BF_NBODY_OK, BF_NBODY_STAGE_REDUCE, BF_NBODY_DETAIL_NONE, 0,
        (double)reduction.sign+0.0*I);
  }
  if(eleIdx == NULL || eleCfg == NULL || eleNum == NULL
     || (NProj > 0 && eleProjCnt == NULL) || eleSpn == NULL
     || eleProjBFCnt == NULL) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
        BF_NBODY_DETAIL_NONE, reduction.order, 0.0+0.0*I);
  }
  if(!isfinite(creal(ip)) || !isfinite(cimag(ip))) {
    return BFNBodyResultValue(
        BF_NBODY_NONFINITE, BF_NBODY_STAGE_RATIO,
        BF_NBODY_DETAIL_NONE, reduction.order, 0.0+0.0*I);
  }
  if(cabs(ip) == 0.0) {
    return BFNBodyResultValue(
        BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_RATIO,
        BF_NBODY_DETAIL_NONE, reduction.order, 0.0+0.0*I);
  }
  for(k=0;k<reduction.order;k++) {
    const int source = scratch->rsj[k];
    const int target = scratch->rsi[k];
    const int particle = eleCfg[source];
    if(particle < 0 || particle >= Nsize || eleNum[source] != 1
       || eleNum[target] != 0 || eleIdx[particle] != source%Nsite
       || eleSpn[particle] != source/Nsite) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
          BF_NBODY_DETAIL_BAD_ELECTRON_LABEL, reduction.order,
          0.0+0.0*I);
    }
  }
  if(BFShouldInjectNBodyFailure(
         scratch, BF_NBODY_INJECT_CANDIDATE)) {
    return BFInjectedNBodyResult(
        BF_NBODY_INJECT_CANDIDATE, reduction.order);
  }
  if(reduction.order <= 2
     && BFShouldInjectNBodyFailure(
            scratch, BF_NBODY_INJECT_PFAFFIAN)) {
    return BFInjectedNBodyResult(
        BF_NBODY_INJECT_PFAFFIAN, reduction.order);
  }

  if(reduction.order <= 2) {
    return BFDispatchReducedNBodyFSZ(
        &reduction, ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
        eleProjBFCnt, scratch);
  }
  return BFRebuildReducedNBodyFSZ(
      &reduction, ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
      scratch);
}

BFNBodyResult GreenFuncNBF_fsz(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt, BFNBodyScratch *scratch) {
  BFNBodyStateSnapshot snapshot;
  BFNBodyResult result;
  int changed;
  if(!BFNBodyStateCheckEnabled) {
    return GreenFuncNBF_fszImpl(
        n, rsi, rsj, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
        eleSpn, eleProjBFCnt, scratch);
  }
  if(BFInitNBodyStateSnapshot(
         &snapshot, 1, eleIdx, eleCfg, eleNum, eleProjCnt,
         eleSpn, eleProjBFCnt) != 0) {
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_WORKSPACE,
        BF_NBODY_DETAIL_STATE_SNAPSHOT, 0, 0.0+0.0*I);
  }
  if(!BFNBodyBaseStateCurrent(
         &snapshot, 1, eleIdx, eleNum, eleSpn, eleProjBFCnt, scratch)) {
    BFFreeNBodyStateSnapshot(&snapshot);
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_DISPATCH,
        BF_NBODY_DETAIL_BASE_STATE_STALE, 0, 0.0+0.0*I);
  }
  result = GreenFuncNBF_fszImpl(
      n, rsi, rsj, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
      eleSpn, eleProjBFCnt, scratch);
  changed = BFNBodyStateChanged(
      &snapshot, 1, eleIdx, eleCfg, eleNum, eleProjCnt,
      eleSpn, eleProjBFCnt);
  BFFreeNBodyStateSnapshot(&snapshot);
#pragma omp atomic update
  BFNBodyStateCheckCount++;
  if(changed) {
#pragma omp atomic update
    BFNBodyStateCheckFailureCount++;
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR,
        result.stage == BF_NBODY_STAGE_NONE
        ? BF_NBODY_STAGE_CANDIDATE : result.stage,
        BF_NBODY_DETAIL_STATE_CHANGED, result.reducedOrder,
        0.0+0.0*I);
  }
  return result;
}
