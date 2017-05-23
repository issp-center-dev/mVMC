#ifndef _LOCGRN_CMP
#define _LOCGRN_CMP

#include <complex.h>

double complex GreenFunc1(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double complex *buffer);
double complex GreenFunc2(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double complex *buffer);

double complex GreenFuncN(const int n, int *rsi, int *rsj, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  double complex *buffer, int *bufferInt);

double complex GreenFunc1BF(const int ri, const int rj, const int s, const double complex ip, double complex* bufM,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew, double complex* buffer);
#endif
