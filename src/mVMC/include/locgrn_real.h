#ifndef _LOCGRN_REAL
#define _LOCGRN_REAL

#include <complex.h>

double GreenFunc1_real(const int ri, const int rj, const int s, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, const double complex *rbmCnt,
                  double complex *rbmCntNew, double *buffer);
double GreenFunc2_real(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, const double complex *rbmCnt,
                  double complex *rbmCntNew, double *buffer);

double GreenFuncN_real(const int n, int *rsi, int *rsj, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  double *buffer, int *bufferInt);

double GreenFunc1BF_real(const int ri, const int rj, const int s, const double ip, double *bufM,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew,
                    double *buffer, int *greenStatus);
int GreenFunc1BF_real_prepare(const int ri, const int rj, const int s, double *greenValue,
                    double *projRatio, double *vecM, int *msaTmp, int *icount,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew);
void GreenFunc1BF_real_finish_batch(const int batchSize, const double ip,
                    const double *projRatio, const int *icount, const int *msaTmp,
                    const double *vecM, double *greenValue, double *pfMNew,
                    double *vecStack, double *wStack);
int GreenFunc2BF_real_ws(const int ri, const int rj, const int rk, const int rl,
                    const int s, const int t, const double ip,
                    double *vecTmp0, double *vecTmp1,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew, double *buffer,
                    double *value);
double GreenFunc2BF_real(const int ri, const int rj, const int rk, const int rl,
                    const int s, const int t, const double ip, double *bufM,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew,
                    double *buffer, int *greenStatus);

#endif
