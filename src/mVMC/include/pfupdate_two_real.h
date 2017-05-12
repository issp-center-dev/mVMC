#ifndef _PFUPDATE_TWO_REAL
#define _PFUPDATE_TWO_REAL

void CalculateNewPfMTwo_real(const int mk, const int t, const int mi, const int s, 
                        double *pfMNew_real, const int *eleIdx,
                        const int qpStart, const int qpEnd, double *buffer);
void CalculateNewPfMTwo2_real(const int ma, const int s, const int mb, const int t,
                         double *pfMNew_real, const int *eleIdx,
                         const int qpStart, const int qpEnd);
void UpdateMAllTwo_real(const int ma, const int s, const int mb, const int t,
                   const int raOld, const int rbOld,
                   const int *eleIdx, const int qpStart, const int qpEnd);
#endif
