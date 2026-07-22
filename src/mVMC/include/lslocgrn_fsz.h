#ifndef _LSLOCGRN_FSZ_H
#define _LSLOCGRN_FSZ_H

#include <complex.h>
#include <stddef.h>

#include "lslocgrn_status.h"

typedef struct {
  int initialized;
  int useReal;
  int stateCheck;
  size_t invCount;
  size_t pfCount;
  size_t gfCount;
  int *eleIdx;
  int *eleCfg;
  int *eleNum;
  int *eleSpn;
  int *projCnt;
  int *gfProjCnt;
  int *stateEleIdx;
  int *stateEleCfg;
  int *stateEleNum;
  int *stateEleSpn;
  int *stateProjCnt;
  double complex *gfComplex;
  double complex *baseInvComplex;
  double complex *basePfComplex;
  double *gfReal;
  double *baseInvReal;
  double *basePfReal;
} LSLanczosFSZScratch;

int LSLanczosFSZScratchInit(LSLanczosFSZScratch *scratch, int useReal);
void LSLanczosFSZScratchFree(LSLanczosFSZScratch *scratch);

int LSLocalQ_fsz(const double complex h1, const double complex ip,
                 const int *eleIdx, const int *eleCfg, const int *eleNum,
                 const int *eleProjCnt, const int *eleSpn,
                 LSLanczosFSZScratch *scratch, double complex *lslq);

int LSLocalQ_fsz_real(const double h1, const double ip,
                      const int *eleIdx, const int *eleCfg,
                      const int *eleNum, const int *eleProjCnt,
                      const int *eleSpn, LSLanczosFSZScratch *scratch,
                      double *lslq);

#endif
