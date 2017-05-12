#ifndef _PFUPDATE_TWO_FCMP
#define _PFUPDATE_TWO_FCMP
#include <complex.h>

void CalculateNewPfMTwo_fcmp(const int mk, const int t, const int mi, const int s, 
                        double complex *pfMNew, const int *eleIdx,
                        const int qpStart, const int qpEnd, double complex *buffer);
void CalculateNewPfMTwo2_fcmp(const int ma, const int s, const int mb, const int t,
                         double complex *pfMNew, const int *eleIdx,
                         const int qpStart, const int qpEnd);
void UpdateMAllTwo_fcmp(const int ma, const int s, const int mb, const int t,
                   const int raOld, const int rbOld,
                   const int *eleIdx, const int qpStart, const int qpEnd);

#endif
