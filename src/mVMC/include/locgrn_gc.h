#ifndef _LOCGRN_GC_H
#define _LOCGRN_GC_H

#include <complex.h>

typedef struct {
  int maxN;
  int pfLWork;
  int *projCntNew;
  int *rsi;
  int *rsj;
  int *msa;
  int *pfIWork;
  double complex *pfMNew;
  double complex *vec;
  double complex *w;
  double complex *smallMat;
  double complex *pfWork;
  double *pfRWork;
} GCGreenScratch;

double complex GreenFunc1GC(int rsi, int rsj, double complex ip, int ncur,
                            int *eleIdx, int *eleCfg, int *eleNum,
                            const int *eleProjCnt,
                            GCGreenScratch *scratch);
double complex GreenFunc2GC(int rsi, int rsj, int rsk, int rsl,
                            double complex ip, int ncur, int *eleIdx,
                            int *eleCfg, int *eleNum,
                            const int *eleProjCnt,
                            GCGreenScratch *scratch);
double complex GreenFuncNGC(int n, int *rsi, int *rsj, double complex ip,
                            int ncur, int *eleIdx, int *eleCfg, int *eleNum,
                            const int *eleProjCnt,
                            GCGreenScratch *scratch);
double complex CalculateHamiltonianGC(double complex ip, int ncur,
                                      int *eleIdx, int *eleCfg, int *eleNum,
                                      const int *eleProjCnt);
double CalculateSzGC(const int *eleNum);
void CalculateGreenFuncGC(double w, double complex ip, int ncur,
                          int *eleIdx, int *eleCfg, int *eleNum,
                          int *eleProjCnt);

#endif
