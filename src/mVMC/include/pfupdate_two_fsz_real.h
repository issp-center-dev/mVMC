#ifndef _PFUPDATE_TWO_FSZ_REAL
#define _PFUPDATE_TWO_FSZ_REAL

void CalculateNewPfMTwo_fsz_real(const int mk, const int t, const int mi, const int s, 
                        double *pfMNew, const int *eleIdx,const int *eleSpn,
                        const int qpStart, const int qpEnd, double *buffer);
void CalculateNewPfMTwo2_fsz_real(const int ma, const int s, const int mb, const int t,
                         double *pfMNew, const int *eleIdx,const int *eleSpn,
                         const int qpStart, const int qpEnd);
void calculateNewPfMTwo_child_fsz_real(const int ma, const int s, const int mb, const int t,
                              double *pfMNew, const int *eleIdx,const int *eleSpn,
                              const int qpStart, const int qpEnd, const int qpidx,
                              double *vec_a, double *vec_b);
void UpdateMAllTwo_fsz_real(const int ma, const int s, const int mb, const int t,
                   const int raOld, const int rbOld,
                   const int *eleIdx,const int*eleSpn, const int qpStart, const int qpEnd);
void updateMAllTwo_child_fsz_real(const int ma, const int s, const int mb, const int t,
                         const int raOld, const int rbOld,
                         const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd, const int qpidx,
                         double *vecP, double *vecQ, double *vecS, double *vecT);


#endif
