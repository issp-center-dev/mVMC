#ifndef _BACKFLOW_NBODY
#define _BACKFLOW_NBODY

#include <complex.h>
#include <stddef.h>

typedef enum {
  BF_NBODY_OK = 0,
  BF_NBODY_PHYSICAL_ZERO = 1,
  BF_NBODY_INVALID_ARGUMENT = 2,
  BF_NBODY_WORKSPACE_ERROR = 3,
  BF_NBODY_PFAFFIAN_ERROR = 4,
  BF_NBODY_NONFINITE = 5
} BFNBodyStatus;

typedef enum {
  BF_NBODY_STAGE_NONE = 0,
  BF_NBODY_STAGE_REDUCE = 1,
  BF_NBODY_STAGE_DISPATCH = 2,
  BF_NBODY_STAGE_CANDIDATE = 3,
  BF_NBODY_STAGE_PROJECTION = 4,
  BF_NBODY_STAGE_SLATER = 5,
  BF_NBODY_STAGE_PFAFFIAN = 6,
  BF_NBODY_STAGE_RATIO = 7,
  BF_NBODY_STAGE_WORKSPACE = 8
} BFNBodyStage;

typedef struct {
  BFNBodyStatus status;
  BFNBodyStage stage;
  /* Interpret detail in the context of stage; numeric domains may overlap. */
  int detail;
  int reducedOrder;
  double complex value;
} BFNBodyResult;

typedef enum {
  BF_NBODY_DETAIL_NONE = 0,
  BF_NBODY_DETAIL_REDUCER = -1,
  BF_NBODY_DETAIL_SPIN_CHANGE = 1,
  BF_NBODY_DETAIL_BAD_ELECTRON_LABEL = 2,
  BF_NBODY_DETAIL_LEGACY_STATE_RESTORE = 3,
  BF_NBODY_DETAIL_INJECTED = 4,
  BF_NBODY_DETAIL_STATE_CHANGED = 5,
  BF_NBODY_DETAIL_BASE_STATE_STALE = 6,
  BF_NBODY_DETAIL_STATE_SNAPSHOT = 7
} BFNBodyDetail;

typedef struct {
  int maxOrder;
  int useFsz;
  size_t intCount;
  size_t complexCount;
  size_t doubleCount;
  size_t inputRsiCount;
  size_t inputRsjCount;
  size_t rsiCount;
  size_t rsjCount;
  size_t movedCount;
  size_t eleIdxCount;
  size_t eleCfgCount;
  size_t eleNumCount;
  size_t eleSpnCount;
  size_t projCntCount;
  size_t projBFCntCount;
  size_t pfIWorkCount;
  size_t affectedCount;
  size_t hopIntWorkCount;
  size_t greenBufferCount;
  size_t slaterCount;
  size_t candidatePfCount;
  size_t pfBufMCount;
  size_t pfWorkCount;
  size_t pfRWorkCount;
} BFNBodyScratchSizes;

typedef struct {
  int maxOrder;
  int useFsz;
  int termIndex;
  int *inputRsi;
  int *inputRsj;
  int *rsi;
  int *rsj;
  int *moved;
  int *eleIdx;
  int *eleCfg;
  int *eleNum;
  int *eleSpn;
  int *projCnt;
  int *projBFCnt;
  int *pfIWork;
  int *affected;
  int *hopIntWork;
  double complex *greenBuffer;
  double complex *slater;
  double complex *candidatePf;
  double complex *pfBufM;
  double complex *pfWork;
  double *pfRWork;
  BFNBodyScratchSizes sizes;
} BFNBodyScratch;

const char *BFNBodyStageName(BFNBodyStage stage);
const char *BFNBodyDetailName(int detail);

typedef struct {
  int maxOrder;
  int useFsz;
  int nsize;
  int nsite2;
  int nproj;
  int nsite;
  int nrange;
  int nqpFull;
  int lapackLWork;
  size_t bfFszGreenBufferCount;
  size_t pfUpdateIntCount;
  size_t pfUpdateDoubleCount;
  int bfHopIntSize;
} BFNBodyScratchDimensions;

int ComputeBFNBodyScratchSizes(const BFNBodyScratchDimensions *dimensions,
                               BFNBodyScratchSizes *sizes);
int GetBFNBodyScratchSizes(int maxOrder, int useFsz,
                           BFNBodyScratchSizes *sizes);
int ValidateBFNBodyScratch(const BFNBodyScratch *scratch, int n,
                           int useFsz);
int BindBFNBodyScratch(const BFNBodyScratchSizes *sizes,
                       int *intBase, size_t intCount,
                       double complex *complexBase, size_t complexCount,
                       double *doubleBase, size_t doubleCount,
                       BFNBodyScratch *scratch);

/*
 * The read-only caller-state arrays must not alias writable state slices in
 * scratch. Both evaluators reject exact pointer aliases before evaluation.
 */
BFNBodyResult GreenFuncNBF(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleProjBFCnt,
    BFNBodyScratch *scratch);

BFNBodyResult GreenFuncNBF_fsz(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt, BFNBodyScratch *scratch);

#endif
