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
#include <string.h>

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

static int BFValidateNBodyScratch(const BFNBodyScratch *scratch, int n) {
  BFNBodyScratchSizes needed;

  if(scratch == NULL || n < 1 || scratch->maxOrder < n
     || scratch->useFsz != 0 || scratch->sizes.maxOrder != scratch->maxOrder
     || scratch->sizes.useFsz != scratch->useFsz
     || scratch->rsi == NULL || scratch->rsj == NULL
     || scratch->moved == NULL || scratch->eleIdx == NULL
     || scratch->eleCfg == NULL || scratch->eleNum == NULL
     || scratch->eleSpn == NULL || scratch->projBFCnt == NULL
     || scratch->pfIWork == NULL || scratch->greenBuffer == NULL
     || scratch->slater == NULL || scratch->candidatePf == NULL
     || scratch->pfBufM == NULL || scratch->pfWork == NULL
     || scratch->pfRWork == NULL
     || (NProj > 0 && scratch->projCnt == NULL)
     || GetBFNBodyScratchSizes(scratch->maxOrder, 0, &needed)
          != BF_NBODY_OK
     || Ne <= 0 || Ne > INT_MAX/2 || Nsize != 2*Ne
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

static BFNBodyResult BFDispatchReducedNBody(
    const NBodyReduction *reduction, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleProjBFCnt,
    BFNBodyScratch *scratch) {
  const size_t slaterCount =
      (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  double complex value;

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

  if(reduction->order == 1) {
    const int target = scratch->rsi[0];
    const int source = scratch->rsj[0];
    const int spin = source/Nsite;
    value = GreenFunc1BF(
        target%Nsite, source%Nsite, spin, ip, scratch->slater,
        scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        eleProjCnt, scratch->projCnt, eleProjBFCnt,
        scratch->projBFCnt, scratch->greenBuffer);
  } else {
    const int target0 = scratch->rsi[0];
    const int source0 = scratch->rsj[0];
    const int target1 = scratch->rsi[1];
    const int source1 = scratch->rsj[1];
    value = GreenFunc2BF(
        target0%Nsite, source0%Nsite,
        target1%Nsite, source1%Nsite,
        source0/Nsite, source1/Nsite, ip, scratch->slater,
        scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        eleProjCnt, scratch->projCnt, eleProjBFCnt,
        scratch->projBFCnt, scratch->greenBuffer);
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

BFNBodyResult GreenFuncNBF(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleProjBFCnt,
    BFNBodyScratch *scratch) {
  NBodyReduction reduction;
  int reduceStatus;
  int k;

  if(!BFValidateNBodyScratch(scratch, n)) {
    if(n < 1) {
      return BFNBodyResultValue(
          BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_REDUCE,
          BF_NBODY_DETAIL_REDUCER, 0, 0.0+0.0*I);
    }
    return BFNBodyResultValue(
        BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_NONE,
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
          reduction.order <= 2 ? BF_NBODY_STAGE_DISPATCH
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
  if(reduction.order <= 2) {
    return BFDispatchReducedNBody(
        &reduction, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
        eleProjBFCnt, scratch);
  }
  return BFRebuildReducedNBody(
      &reduction, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
      eleProjBFCnt, scratch);
}
