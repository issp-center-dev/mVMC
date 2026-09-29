#ifndef _PFUPDATE_REAL
#define _PFUPDATE_REAL
void CalculateNewPfM_real(const int mi, const int s, double *pfMNew_real, const int *eleIdx,
                     const int qpStart, const int qpEnd);
void CalculateNewPfM2_real(const int mi, const int s, double *pfMNew_real, const int *eleIdx,
                     const int qpStart, const int qpEnd);
void UpdateMAll_real(const int mi, const int s, const int *eleIdx,
                const int qpStart, const int qpEnd);

int CalculateNewPfMBF_realChecked(const int *icount, const int *msaTmp, double *pfMNew, const int *eleIdx, int qpStart, int qpEnd, const double *bufM);
void CalculateNewPfMBF_real(const int *icount, const int *msaTmp,
                            double *pfMNew, const int *eleIdx,
                            const int qpStart, const int qpEnd, const double *bufM);
void CalculateNewPfMBFWithStride_real(const int *icount, const int *msaTmp, const int msaStride,
                            double *pfMNew, const int *eleIdx,
                            const int qpStart, const int qpEnd, const double *bufM);
void CalculateNewPfMBFVecWithStride_real(const int *icount, const int *msaTmp, const int msaStride,
                            double *pfMNew, const int qpStart, const int qpEnd,
                            const double *vecM, const int vecStride, const int *eleIdx);
void CalculateNewPfMBFVec_real(const int *icount, const int *msaTmp,
                            double *pfMNew, const int qpStart, const int qpEnd,
                            const double *vecM, const int *eleIdx);
void CalculateNewPfMBFVecBatched_real(const int batchSize, const int *icount, const int *msaTmp,
                            double *pfMNew, const int qpStart, const int qpEnd,
                            const double *vecM, const int *eleIdx);

int UpdateMAll_BF_real(const int *icount, const int *msaTmp, double *pfMNew, const int *eleIdx, int qpStart, int qpEnd, double *candidateInv);


#endif
