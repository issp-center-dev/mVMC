#ifndef _LSLOCGRN_BF_FSZ_H
#define _LSLOCGRN_BF_FSZ_H

#include <complex.h>
#include <stddef.h>

#include "lslocgrn_status.h"

typedef struct {
  int initialized;
  int stateCheck;
  int hopIntWorkSize;
  int realGateIndex;
  double realGateTolerance;
  double complex realGateValue;
  size_t invCount;
  size_t pfCount;
  size_t slaterCount;
  size_t bfCount;
  size_t greenBufferCount;
  size_t pfIntCount;
  size_t pfDoubleCount;
  int *eleIdx;
  int *eleCfg;
  int *eleNum;
  int *eleSpn;
  int *projCnt;
  int *bfCnt;
  int *gfProjCnt;
  int *gfBFCnt;
  int *affected;
  int *hopIntWork;
  int *pfIWork;
  int *stateEleIdx;
  int *stateEleCfg;
  int *stateEleNum;
  int *stateEleSpn;
  int *stateProjCnt;
  int *stateBFCnt;
  double *pfRWork;
  double complex *greenBuffer;
  double complex *pfBufM;
  double complex *pfWork;
  double complex *candidateSlater;
  double complex *candidateInv;
  double complex *candidatePf;
  double complex *baseSlater;
  double complex *baseInv;
  double complex *basePf;
} LSLanczosBFFSZScratch;

int LSLanczosBFFSZScratchInit(LSLanczosBFFSZScratch *scratch);
void LSLanczosBFFSZScratchFree(LSLanczosBFFSZScratch *scratch);

int LSLocalQBF_fsz(const double complex h1, const double complex ip,
                   const int *eleIdx, const int *eleCfg, const int *eleNum,
                   const int *eleProjCnt, const int *eleSpn,
                   const int *eleProjBFCnt,
                   LSLanczosBFFSZScratch *scratch, double complex *lslq);

int LSLocalQBF_fsz_real(const double complex h1, const double complex ip,
                        const int *eleIdx, const int *eleCfg,
                        const int *eleNum, const int *eleProjCnt,
                        const int *eleSpn, const int *eleProjBFCnt,
                        LSLanczosBFFSZScratch *scratch, double *lslq);

#endif
