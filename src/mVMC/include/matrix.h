#ifndef _MATRIX
#define _MATRIX

#define BF_FSZ_PF_OK 0
#define BF_FSZ_PF_LAPACK_FAILURE 1
#define BF_FSZ_PF_NONFINITE 2

int CalculateMAll_fcmp(const int *eleIdx, const int qpStart, const int qpEnd);
int CalculateMAll_real(const int *eleIdx, const int qpStart, const int qpEnd);

int CalculateMAll_BF_real(const int *eleIdx, const int qpStart, const int qpEnd);
int CalculateMAll_BF_fcmp(const int *eleIdx, const int qpStart, const int qpEnd);
int CalculateMAll_BF_fsz(const int *eleIdx, const int *eleSpn, const int qpStart, const int qpEnd);
int CalculatePfM_BF_fsz(const int *eleIdx, const int *eleSpn, const int qpStart, const int qpEnd,
    double complex *pfMOut, int *failureDetail);

int CalculateMAll_fsz(const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd);
int CalculateMAll_fsz_real(const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd);

int calculateMAll_child_real(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double *bufM, int *iwork, double *work, int lwork, double* pfM_real, double *invM_real);
int calculateMAll_child_fcmp(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double complex *bufM, int *iwork, double complex *work, int lwork,double *rwork);

int calculateMAll_BF_real_child(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double *bufM, int *iwork, double *work, int lwork, double* pfM_real, double *invM_real);
int calculateMAll_BF_fcmp_child(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double complex*bufM, int *iwork, double complex*work, int lwork, double *rwork, double complex* pfM, double complex*invM);
int calculateMAll_BF_fsz_child(const int *eleIdx, const int *eleSpn, const int qpStart, const int qpEnd,
    const int qpidx, double complex *bufM, int *iwork, double complex *work, int lwork,
    double *rwork, double complex *pfM, double complex *invM);
int calculatePfM_BF_fsz_child(const int *eleIdx, const int *eleSpn, const int qpStart,
    const int qpidx, double complex *bufM, int *iwork, double complex *work, int lwork,
    double *rwork, double complex *pfMOut, int *failureDetail);

int calculateMAll_child_fsz(const int *eleIdx,const int *elesSpn, const int qpStart, const int qpEnd, const int qpidx,
    double complex *bufM, int *iwork, double complex *work, int lwork,double *rwork);
// note: CalculateMAll_fsz,calculateMAll_child_fsz will be merged with *fcmp
int calculateMAll_child_fsz_real(const int *eleIdx, const int *elesSpn, const int qpStart, const int qpEnd, const int qpidx,
    double *bufM, int *iwork, double *work, int lwork);


#endif
