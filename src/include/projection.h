#ifndef _PROJECTION_
#define _PROJECTION_

double LogProjVal(const int *projCnt);
double LogProjRatio(const int *projCntNew, const int *projCntOld);
double ProjRatio(const int *projCntNew, const int *projCntOld);
void MakeProjCnt(int *projCnt, const int *eleNum);
void UpdateProjCnt(const int ri, const int rj, const int s,
				   int *projCntNew, const int *projCntOld,
				   const int *eleNum);

void MakeProjBFCnt(int *projCnt, const int *eleNum);

#endif
