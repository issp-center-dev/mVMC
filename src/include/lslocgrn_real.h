#ifndef _LSLOCGRN_REAL
#define _LSLOCGRN_REAL
#include <complex.h>
#include "../workspace.c"
#include "../calham_real.c"
#include "../pfupdate_real.c"
#include "../qp_real.c"

void LSLocalQ_real(const double h1, const double ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt);

double calculateHK_real(const double h1, const double ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt);

double calHCA_real(const int ri, const int rj, const int s,
                   const double h1, const double ip, int *eleIdx, int *eleCfg,
                   int *eleNum, int *eleProjCnt);

double checkGF1_real(const int ri, const int rj, const int s, const double ip,
                int *eleIdx, const int *eleCfg, int *eleNum);

double calHCA1_real(const int ri, const int rj, const int s,
               const double ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt);
double calHCA2_real(const int ri, const int rj, const int s,
                       const double ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt);


#endif
