#ifndef _CALHAM_REAL
#define _CALHAM_REAL

double CalculateHamiltonian_real(const double ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt);

double CalculateHamiltonianBF_real(const double ip, int *eleIdx, const int *eleCfg,
                              int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt);


double CalculateHamiltonian0_real(const int *eleNum);

#endif
