#ifndef _RBM_
#define _RBM_

extern inline double complex WeightRBM(const double complex *rbmCnt);
extern inline double complex LogWeightRBM(const double complex *rbmCnt);
extern inline double complex RBMRatio(const double complex *rbmCntNew, const double complex *rbmCntOld);
extern inline double complex LogRBMRatio(const double complex *rbmCntNew, const double complex *rbmCntOld);

void copyFromBurnSampleRBM(double complex *rbmCnt);
void copyToBurnSampleRBM(double complex *rbmCnt);
void saveRBMCnt(const int sample, const double complex *rbmCnt);

void RBMDiff(double complex *srOptO, const double complex *rbmCnt, const int *eleNum);


void MakeRBMCnt(double complex *rbmCnt, const int *eleNum);
void UpdateRBMCnt(const int ri, const int rj, const int s,
                   double complex *rbmCntNew, const double complex *rbmCntOld, const int *eleNum);

void UpdateRBMCnt_fsz(const int ri, const int rj, const int s, const int t,
                   double complex *rbmCntNew, const double complex *rbmCntOld) ;
#endif
