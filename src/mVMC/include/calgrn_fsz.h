#ifndef _CALGRN_FSZ
#define _CALGRN_FSZ
#include <complex.h>

void CalculateGreenFunc_fsz(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                         int *eleNum, int *eleSpn, int *eleProjCnt);
#endif
