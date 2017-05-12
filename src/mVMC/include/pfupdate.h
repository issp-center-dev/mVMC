#ifndef _PFUPDATE
#define _PFUPDATE
#include <complex.h>
void CalculateNewPfM(const int mi, const int s, double complex *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd);
void CalculateNewPfM2(const int mi, const int s, double complex *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd);
void UpdateMAll(const int mi, const int s, const int *eleIdx,
                const int qpStart, const int qpEnd);
void updateMAll_child(const int ma, const int s, const int *eleIdx,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double complex *vec1, double complex *vec2);

void CalculateNewPfMBF(const int *icount, const int *msaTmp,double complex*pfMNew, const int *eleIdx,
                       const int qpStart, const int qpEnd, const double complex*bufM) ;

double complex calculateNewPfMBFN4_child(const int qpidx, const int n, const int *msa,
                                 const int *eleIdx, const double complex* bufM);

void UpdateMAll_BF_fcmp(const int *icount, const int *msaTmp,
                        double complex* pfMNew, const int *eleIdx,
                        const int qpStart, const int qpEnd) ;

#endif


