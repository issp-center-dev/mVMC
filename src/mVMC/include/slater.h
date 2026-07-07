#ifndef _SLATER
#define _SLATER
void UpdateSlaterElm_fcmp();
void SlaterElmDiff_fcmp(double complex *srOptO, const double complex ip, int *eleIdx);

void SlaterElmBFDiff_fcmp(double complex*srOptO, const double complex ip, int *eleIdx, int *eleNum, int *eleCfg, int *eleProjConst,const int * eleProjBFCnt);

void BackFlowDiff_fcmp(complex double *srOptO, const double complex ip, int *eleIdx, const int *eleNum, int *eleProjConst,
                       const int *eleProjBFCnt);

void MakeSlaterElmBF_fcmp(const int *eleNum, const int *eleProjBFCnt);
void MakeSlaterElmBF_fsz(const int *eleNum, const int *eleProjBFCnt);
void UpdateSlaterElmBF_fcmp(const int ma, const int ra, const int rb, const int u,
                       const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double complex*sltElmTmp);
void UpdateSlaterElmBF_real(const int ma, const int ra, const int rb, const int u,
                       const int *eleCfg, const int *eleNum, const int *eleProjBFCntOld,
                       const int *eleProjBFCnt, int *msa, int *hopNum, double *sltElmTmp);
void UpdateSlaterElmBFGrn(const int ma, const int ra, const int rb, const int u,
                          const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double complex*sltElmTmp);

void UpdateSlaterElmBFGrn_real(const int ma, const int ra, const int rb, const int u,
                          const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double *sltElmTmp,
                          int *restoreRows, int *restoreRowCount);
void UpdateSlaterElmBFGrnVec_real(const int ma, const int ra, const int rb, const int u,
                          const int *eleIdx, const int *eleCfg, const int *eleNum,
                          const int *eleProjBFCntOld, const int *eleProjBFCnt,
                          int *msa, int *hopNum, double *vecTmp);

void StoreSlaterElmBF_fcmp(complex double *bufM);

void StoreSlaterElmBF_real(double *bufM_real);

void StoreSlaterElmBF_real_rows(double *bufM_real, const int *restoreRows, const int *restoreRowCount);

#endif
