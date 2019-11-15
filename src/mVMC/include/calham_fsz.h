#ifndef _CALHAM_FSZ
#define _CALHAM_FSZ
#include <complex.h>

double CalculateSz_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum,const int *eleProjCnt,int *eleSpn);
double complex CalculateHamiltonian_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum,const int *eleProjCnt,int *eleSpn);
double complex CalculateHamiltonian0_fsz(const int *eleNum);
double complex CalculateHamiltonian1_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum,const int *eleProjCnt,int *eleSpn);
double complex CalculateHamiltonian2_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn);


#endif
