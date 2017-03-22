#ifndef _PROJECTION_
#define _PROJECTION_

#include <complex.h>

double LogProjVal(const int *projCnt);
double LogProjRatio(const int *projCntNew, const int *projCntOld);
double ProjRatio(const int *projCntNew, const int *projCntOld);
void MakeProjBFCnt(int *projCnt, const int *eleNum);

#endif
