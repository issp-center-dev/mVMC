#ifndef _LOCGRN_FSZ_H
#define _LOCGRN_FSZ_H

#include <complex.h>
#include <stddef.h>

double complex GreenFunc1_fsz2(const int ri, const int rj,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn,
                  int *projCntNew, double complex *buffer);

double complex GreenFunc2_fsz2(const int ri, const int rj,
                  const int rk, const int rl, const int s, const int t,
                  const int u, const int v, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn,
                  int *projCntNew, double complex *buffer);

double complex GreenFunc1BF_fsz2_workspace(
                  const int ri, const int rj, const int s, const int t,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn, int *projCntNew,
                  const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, double complex *pfBufM,
                  int *affected, int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double complex *pfWork, double *pfRWork);

double complex GreenFunc2BF_fsz2(
                  const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const int u, const int v,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn, int *projCntNew,
                  const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

int GetGreenFuncBF_fsz_buffer_work_size(
    size_t *bufferComplexCount, size_t *pfUpdateIntCount,
    size_t *pfUpdateDoubleCount);

#endif
