#ifndef _PFUPDATE_FSZ
#define _PFUPDATE_FSZ

void CalculateNewPfM_fsz(const int mi, const int s, double complex *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd);
void CalculateNewPfM2_fsz(const int ma, const int s, double complex *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd);
void UpdateMAll_fsz(const int mi, const int s, const int *eleIdx,const int *eleSpn,
                const int qpStart, const int qpEnd);
void updateMAll_child_fsz(const int ma, const int s, const int *eleIdx,const int *eleSpn,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double complex *vec1, double complex *vec2);



#endif


