#ifndef _SLATER
#define _SLATER
void UpdateSlaterElm_fcmp();
void SlaterElmDiff_fcmp(double complex *srOptO, const double complex ip, int *eleIdx);
void MakeSlaterElmBF_fcmp(const int *eleNum, const int *eleProjBFCnt);
void UpdateSlaterElmBF_fcmp(const int ma, const int ra, const int rb, const int u,
                       const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double complex*sltElmTmp);

#endif


