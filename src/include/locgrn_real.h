#ifndef _LOCGRN_REAL
#define _LOCGRN_REAL
#include <complex.h>
#include "global.h"
#include "../projection.c"
#include "../qp_real.c"

double GreenFunc1_real(const int ri, const int rj, const int s, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double *buffer);
double GreenFunc2_real(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double *buffer);

double GreenFuncN_real(const int n, int *rsi, int *rsj, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  double *buffer, int *bufferInt);

#endif
