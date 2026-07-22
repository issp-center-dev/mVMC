#ifndef _LSLOCGRN_BF_H
#define _LSLOCGRN_BF_H

#include <complex.h>
#include <stddef.h>

#include "lslocgrn_status.h"

typedef struct {
  int initialized;
  int useReal;
  int stateCheck;
  size_t invCount;
  size_t pfCount;
  size_t slaterCount;
  size_t bfCount;
  size_t etaCount;
  size_t gfCount;
  int *eleIdx;
  int *eleCfg;
  int *eleNum;
  int *projCnt;
  int *bfCnt;
  int *gfProjCnt;
  int *gfBFCnt;
  int *stateEleIdx;
  int *stateEleCfg;
  int *stateEleNum;
  int *stateProjCnt;
  int *stateBFCnt;
  int *baseEtaFlag;
  double complex *baseEta;
  double complex *baseSlaterComplex;
  double complex *gfSlaterComplex;
  double complex *gfComplex;
  double complex *baseInvComplex;
  double complex *basePfComplex;
  double *baseSlaterReal;
  double *gfSlaterReal;
  double *gfReal;
  double *baseInvReal;
  double *basePfReal;
} LSLanczosBFScratch;

int LSLanczosBFScratchInit(LSLanczosBFScratch *scratch, int useReal);
void LSLanczosBFScratchFree(LSLanczosBFScratch *scratch);

int LSLocalQBF(const double complex h1, const double complex ip,
               const int *eleIdx, const int *eleCfg, const int *eleNum,
               const int *eleProjCnt, const int *eleProjBFCnt,
               LSLanczosBFScratch *scratch, double complex *lslq);

int LSLocalQBF_real(const double h1, const double ip,
                    const int *eleIdx, const int *eleCfg, const int *eleNum,
                    const int *eleProjCnt, const int *eleProjBFCnt,
                    LSLanczosBFScratch *scratch, double *lslq);

#endif
