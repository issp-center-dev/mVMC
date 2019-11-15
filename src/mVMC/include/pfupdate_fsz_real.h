#ifndef _PFUPDATE_FSZ_REAL
#define _PFUPDATE_FSZ_REAL

void CalculateNewPfM_fsz_real(const int mi, const int s, double *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd);
void CalculateNewPfM2_fsz_real(const int ma, const int s, double *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd);
void UpdateMAll_fsz_real(const int mi, const int s, const int *eleIdx,const int *eleSpn,
                const int qpStart, const int qpEnd);
void updateMAll_child_fsz_real(const int ma, const int s, const int *eleIdx,const int *eleSpn,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double *vec1, double *vec2);


#endif


