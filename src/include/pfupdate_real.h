#ifndef _PFUPDATE_REAL
#define _PFUPDATE_REAL
void CalculateNewPfMBF_real(const int *icount, const int *msaTmp,
                            double *pfMNew, const int *eleIdx,
                            const int qpStart, const int qpEnd, const double *bufM);

void UpdateMAll_BF_real(const int *icount, const int *msaTmp,
                  double *pfMNew, const int *eleIdx,
                  const int qpStart, const int qpEnd) ;

#endif


