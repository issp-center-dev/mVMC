#ifndef _CALHAM
#define _CALHAM
#include <complex.h>

double complex CalculateHamiltonian(const double complex ip, int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt);

double complex CalculateHamiltonian0(const int *eleNum);

double CalculateDoubleOccupation(int *eleIdx, const int *eleCfg,
                                 int *eleNum, const int *eleProjCnt);

double complex CalculateHamiltonianBF_fcmp(const double complex ip, int *eleIdx, const int *eleCfg,
                                   int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt);

#endif
