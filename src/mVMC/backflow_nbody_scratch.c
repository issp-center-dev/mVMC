#include "./include/backflow_nbody.h"

#include <limits.h>
#include <stdint.h>
#include <string.h>

static int BFCheckedSizeAdd(size_t a, size_t b, size_t *out) {
  if(out == NULL || a > SIZE_MAX-b) return -1;
  *out = a+b;
  return 0;
}

static int BFCheckedSizeMul(size_t a, size_t b, size_t *out) {
  if(out == NULL || (b != 0 && a > SIZE_MAX/b)) return -1;
  *out = a*b;
  return 0;
}

static int BFAddSlice(size_t *total, size_t count) {
  size_t next;
  if(BFCheckedSizeAdd(*total, count, &next) != 0) return -1;
  *total = next;
  return 0;
}

int ComputeBFNBodyScratchSizes(const BFNBodyScratchDimensions *dimensions,
                               BFNBodyScratchSizes *sizes) {
  BFNBodyScratchSizes value;
  size_t nsiteRange;
  size_t nsite2Square;
  size_t nsizeSquare;
  size_t baseGreenCount;

  if(dimensions == NULL || sizes == NULL
     || dimensions->maxOrder < 1
     || (dimensions->useFsz != 0 && dimensions->useFsz != 1)
     || dimensions->nsize <= 0 || dimensions->nsite2 <= 0
     || dimensions->nproj < 0 || dimensions->nsite <= 0
     || dimensions->nrange <= 0 || dimensions->nqpFull <= 0
     || dimensions->lapackLWork <= 0
     || dimensions->nsite > INT_MAX/2
     || dimensions->nsite2 != 2*dimensions->nsite
     || dimensions->nsize > dimensions->nsite2
     || (dimensions->useFsz
         && (dimensions->bfFszGreenBufferCount == 0
             || dimensions->pfUpdateIntCount == 0
             || dimensions->pfUpdateDoubleCount == 0
             || dimensions->bfHopIntSize <= 0))) {
    return BF_NBODY_WORKSPACE_ERROR;
  }

  memset(&value, 0, sizeof(value));
  value.maxOrder = dimensions->maxOrder;
  value.useFsz = dimensions->useFsz;
  value.inputRsiCount = (size_t)dimensions->maxOrder;
  value.inputRsjCount = (size_t)dimensions->maxOrder;
  value.rsiCount = (size_t)dimensions->maxOrder;
  value.rsjCount = (size_t)dimensions->maxOrder;
  value.movedCount = (size_t)dimensions->maxOrder;
  value.eleIdxCount = (size_t)dimensions->nsize;
  value.eleCfgCount = (size_t)dimensions->nsite2;
  value.eleNumCount = (size_t)dimensions->nsite2;
  value.eleSpnCount = (size_t)dimensions->nsize;
  value.projCntCount = (size_t)dimensions->nproj;

  if(BFCheckedSizeMul((size_t)dimensions->nsite,
                      (size_t)dimensions->nrange, &nsiteRange) != 0
     || BFCheckedSizeMul(nsiteRange, 16u, &value.projBFCntCount) != 0
     || BFCheckedSizeMul((size_t)dimensions->nsite2,
                         (size_t)dimensions->nsite2, &nsite2Square) != 0
     || BFCheckedSizeMul((size_t)dimensions->nsize,
                         (size_t)dimensions->nsize, &nsizeSquare) != 0
     || BFCheckedSizeMul((size_t)dimensions->nqpFull, nsite2Square,
                         &value.slaterCount) != 0
     || BFCheckedSizeMul((size_t)dimensions->nsize, 2u,
                         &baseGreenCount) != 0
     || BFCheckedSizeAdd((size_t)dimensions->nqpFull, baseGreenCount,
                         &baseGreenCount) != 0) {
    return BF_NBODY_WORKSPACE_ERROR;
  }

  value.pfIWorkCount = (size_t)dimensions->nsize;
  value.greenBufferCount = baseGreenCount;
  value.pfRWorkCount = (size_t)dimensions->lapackLWork;
  if(dimensions->useFsz) {
    if(dimensions->pfUpdateIntCount > value.pfIWorkCount) {
      value.pfIWorkCount = dimensions->pfUpdateIntCount;
    }
    if(dimensions->bfFszGreenBufferCount > value.greenBufferCount) {
      value.greenBufferCount = dimensions->bfFszGreenBufferCount;
    }
    if(dimensions->pfUpdateDoubleCount > value.pfRWorkCount) {
      value.pfRWorkCount = dimensions->pfUpdateDoubleCount;
    }
    value.affectedCount = (size_t)dimensions->nsize;
    value.hopIntWorkCount = (size_t)dimensions->bfHopIntSize;
  }
  value.candidatePfCount = (size_t)dimensions->nqpFull;
  value.pfBufMCount = nsizeSquare;
  value.pfWorkCount = (size_t)dimensions->lapackLWork;

  if(BFAddSlice(&value.intCount, value.inputRsiCount) != 0
     || BFAddSlice(&value.intCount, value.inputRsjCount) != 0
     || BFAddSlice(&value.intCount, value.rsiCount) != 0
     || BFAddSlice(&value.intCount, value.rsjCount) != 0
     || BFAddSlice(&value.intCount, value.movedCount) != 0
     || BFAddSlice(&value.intCount, value.eleIdxCount) != 0
     || BFAddSlice(&value.intCount, value.eleCfgCount) != 0
     || BFAddSlice(&value.intCount, value.eleNumCount) != 0
     || BFAddSlice(&value.intCount, value.eleSpnCount) != 0
     || BFAddSlice(&value.intCount, value.projCntCount) != 0
     || BFAddSlice(&value.intCount, value.projBFCntCount) != 0
     || BFAddSlice(&value.intCount, value.pfIWorkCount) != 0
     || BFAddSlice(&value.intCount, value.affectedCount) != 0
     || BFAddSlice(&value.intCount, value.hopIntWorkCount) != 0
     || BFAddSlice(&value.complexCount, value.greenBufferCount) != 0
     || BFAddSlice(&value.complexCount, value.slaterCount) != 0
     || BFAddSlice(&value.complexCount, value.candidatePfCount) != 0
     || BFAddSlice(&value.complexCount, value.pfBufMCount) != 0
     || BFAddSlice(&value.complexCount, value.pfWorkCount) != 0) {
    return BF_NBODY_WORKSPACE_ERROR;
  }
  value.doubleCount = value.pfRWorkCount;

  if(value.intCount > INT_MAX || value.complexCount > INT_MAX
     || value.doubleCount > INT_MAX) {
    return BF_NBODY_WORKSPACE_ERROR;
  }

  *sizes = value;
  return BF_NBODY_OK;
}

static int BFValidateScratchSizes(const BFNBodyScratchSizes *sizes) {
  size_t intTotal = 0;
  size_t complexTotal = 0;

  if(sizes == NULL || sizes->maxOrder < 1
     || (sizes->useFsz != 0 && sizes->useFsz != 1)
     || sizes->intCount == 0 || sizes->complexCount == 0
     || sizes->doubleCount == 0
     || sizes->inputRsiCount != (size_t)sizes->maxOrder
     || sizes->inputRsjCount != (size_t)sizes->maxOrder
     || sizes->rsiCount != (size_t)sizes->maxOrder
     || sizes->rsjCount != (size_t)sizes->maxOrder
     || sizes->movedCount != (size_t)sizes->maxOrder
     || (!sizes->useFsz
         && (sizes->affectedCount != 0 || sizes->hopIntWorkCount != 0))) {
    return -1;
  }

  if(BFAddSlice(&intTotal, sizes->inputRsiCount) != 0
     || BFAddSlice(&intTotal, sizes->inputRsjCount) != 0
     || BFAddSlice(&intTotal, sizes->rsiCount) != 0
     || BFAddSlice(&intTotal, sizes->rsjCount) != 0
     || BFAddSlice(&intTotal, sizes->movedCount) != 0
     || BFAddSlice(&intTotal, sizes->eleIdxCount) != 0
     || BFAddSlice(&intTotal, sizes->eleCfgCount) != 0
     || BFAddSlice(&intTotal, sizes->eleNumCount) != 0
     || BFAddSlice(&intTotal, sizes->eleSpnCount) != 0
     || BFAddSlice(&intTotal, sizes->projCntCount) != 0
     || BFAddSlice(&intTotal, sizes->projBFCntCount) != 0
     || BFAddSlice(&intTotal, sizes->pfIWorkCount) != 0
     || BFAddSlice(&intTotal, sizes->affectedCount) != 0
     || BFAddSlice(&intTotal, sizes->hopIntWorkCount) != 0
     || BFAddSlice(&complexTotal, sizes->greenBufferCount) != 0
     || BFAddSlice(&complexTotal, sizes->slaterCount) != 0
     || BFAddSlice(&complexTotal, sizes->candidatePfCount) != 0
     || BFAddSlice(&complexTotal, sizes->pfBufMCount) != 0
     || BFAddSlice(&complexTotal, sizes->pfWorkCount) != 0
     || intTotal != sizes->intCount
     || complexTotal != sizes->complexCount
     || sizes->pfRWorkCount != sizes->doubleCount) {
    return -1;
  }
  return 0;
}

int BindBFNBodyScratch(const BFNBodyScratchSizes *sizes,
                       int *intBase, size_t intCount,
                       double complex *complexBase, size_t complexCount,
                       double *doubleBase, size_t doubleCount,
                       BFNBodyScratch *scratch) {
  BFNBodyScratch value;
  size_t intOffset = 0;
  size_t complexOffset = 0;

  if(BFValidateScratchSizes(sizes) != 0 || intBase == NULL
     || complexBase == NULL || doubleBase == NULL || scratch == NULL
     || intCount < sizes->intCount || complexCount < sizes->complexCount
     || doubleCount < sizes->doubleCount) {
    return BF_NBODY_WORKSPACE_ERROR;
  }

  memset(&value, 0, sizeof(value));
  value.maxOrder = sizes->maxOrder;
  value.useFsz = sizes->useFsz;
#define BF_BIND_INT(field, countField)                                         \
  do {                                                                         \
    value.field = sizes->countField == 0 ? NULL : intBase+intOffset;           \
    intOffset += sizes->countField;                                             \
  } while(0)
  BF_BIND_INT(inputRsi, inputRsiCount);
  BF_BIND_INT(inputRsj, inputRsjCount);
  BF_BIND_INT(rsi, rsiCount);
  BF_BIND_INT(rsj, rsjCount);
  BF_BIND_INT(moved, movedCount);
  BF_BIND_INT(eleIdx, eleIdxCount);
  BF_BIND_INT(eleCfg, eleCfgCount);
  BF_BIND_INT(eleNum, eleNumCount);
  BF_BIND_INT(eleSpn, eleSpnCount);
  BF_BIND_INT(projCnt, projCntCount);
  BF_BIND_INT(projBFCnt, projBFCntCount);
  BF_BIND_INT(pfIWork, pfIWorkCount);
  BF_BIND_INT(affected, affectedCount);
  BF_BIND_INT(hopIntWork, hopIntWorkCount);
#undef BF_BIND_INT

#define BF_BIND_COMPLEX(field, countField)                                     \
  do {                                                                         \
    value.field = sizes->countField == 0 ? NULL : complexBase+complexOffset;   \
    complexOffset += sizes->countField;                                         \
  } while(0)
  BF_BIND_COMPLEX(greenBuffer, greenBufferCount);
  BF_BIND_COMPLEX(slater, slaterCount);
  BF_BIND_COMPLEX(candidatePf, candidatePfCount);
  BF_BIND_COMPLEX(pfBufM, pfBufMCount);
  BF_BIND_COMPLEX(pfWork, pfWorkCount);
#undef BF_BIND_COMPLEX
  value.pfRWork = sizes->pfRWorkCount == 0 ? NULL : doubleBase;
  value.sizes = *sizes;

  *scratch = value;
  return BF_NBODY_OK;
}
