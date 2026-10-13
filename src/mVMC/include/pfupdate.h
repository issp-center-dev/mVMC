#ifndef _PFUPDATE
#define _PFUPDATE
#include <complex.h>
#include "backflow_stable.h"
void CalculateNewPfM(const int mi, const int s, double complex *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd);
void CalculateNewPfM2(const int mi, const int s, double complex *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd);
void UpdateMAll(const int mi, const int s, const int *eleIdx,
                const int qpStart, const int qpEnd);
void updateMAll_child(const int ma, const int s, const int *eleIdx,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double complex *vec1, double complex *vec2);

int CalculateNewPfMBFChecked(const int *icount, const int *msaTmp, double complex *pfMNew, const int *eleIdx, int qpStart, int qpEnd, const double complex *bufM);
int CalculateNewPfMBFCheckedWorkspace(const int *icount, const int *msaTmp,
    double complex *pfMNew, const int *eleIdx, int qpStart, int qpEnd,
    const double complex *bufM, BFStableWorkspaceFcmp *scratch);
void CalculateNewPfMBF(const int *icount, const int *msaTmp,double complex*pfMNew, const int *eleIdx,
                       const int qpStart, const int qpEnd, const double complex*bufM) ;
void CalculateNewPfMBFWithStride(const int *icount, const int *msaTmp, const int msaStride,
                       double complex*pfMNew, const int *eleIdx,
                       const int qpStart, const int qpEnd, const double complex*bufM) ;
void CalculateNewPfMBFWithStrideWorkspace(const int *icount,
    const int *msaTmp, int msaStride, double complex *pfMNew,
    const int *eleIdx, int qpStart, int qpEnd, const double complex *bufM,
    BFStableWorkspaceFcmp *scratch);

double complex calculateNewPfMBFN4_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                 const int *eleIdx, const double complex* bufM, int *status);

int UpdateMAll_BF_fcmp(const int *icount, const int *msaTmp, double complex *pfMNew, const int *eleIdx, int qpStart, int qpEnd, double complex *candidateInv);
int UpdateMAll_BF_fcmpWorkspace(const int *icount, const int *msaTmp,
    double complex *pfMNew, const int *eleIdx, int qpStart, int qpEnd,
    double complex *candidateInv, BFStableWorkspaceFcmp *scratch);


#endif
