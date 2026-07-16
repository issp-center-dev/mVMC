#ifndef _PFUPDATE_FSZ
#define _PFUPDATE_FSZ

#define BF_FSZ_PF_UPDATE_OK 0
#define BF_FSZ_PF_UPDATE_EXACT_ZERO 1
#define BF_FSZ_PF_UPDATE_LAPACK_FAILURE 2
#define BF_FSZ_PF_UPDATE_NONFINITE 3
#define BF_FSZ_PF_UPDATE_INVALID_ARGUMENT 4

#define BF_FSZ_INV_UPDATE_OK 0
#define BF_FSZ_INV_UPDATE_LAPACK_FAILURE 1
#define BF_FSZ_INV_UPDATE_NONFINITE 2
#define BF_FSZ_INV_UPDATE_ANTISYMMETRY 3
#define BF_FSZ_INV_UPDATE_RESIDUAL 4
#define BF_FSZ_INV_UPDATE_INVALID_ARGUMENT 5

#define BF_FSZ_INV_STAGE_NONE 0
#define BF_FSZ_INV_STAGE_GETRF 1
#define BF_FSZ_INV_STAGE_GETRI 2
#define BF_FSZ_INV_STAGE_FINITE 3
#define BF_FSZ_INV_STAGE_ANTISYMMETRY 4
#define BF_FSZ_INV_STAGE_RESIDUAL 5

typedef struct {
  int status;
  int stage;
  int qpidx;
  int lapackInfo;
  double antisymmetryResidual;
  double affectedResidual;
  double checkSeconds;
} BF_FSZ_InvUpdateResult;
void CalculateNewPfM_fsz(const int mi, const int s, double complex *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd);
void CalculateNewPfM2_fsz(const int ma, const int s, double complex *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd);
void UpdateMAll_fsz(const int mi, const int s, const int *eleIdx,const int *eleSpn,
                const int qpStart, const int qpEnd);
void updateMAll_child_fsz(const int ma, const int s, const int *eleIdx,const int *eleSpn,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double complex *vec1, double complex *vec2);



/* F3-a BF-FSZ candidate Pfaffian update. oldPfM and oldInvM use rank-local
   qp indexing, while candidateSlater uses the global qp index. The workspace
   API performs no allocation and never modifies global PfM/InvM state. */
int GetCalculateNewPfMBF_fsz_rows_work_size(size_t *complexCount,
    size_t *intCount, size_t *doubleCount);
int BF_FSZ_ShouldUseFullPfaffian(const int nAffected);
int CalculateNewPfMBF_fsz_rows(
    const int nAffected, const int *affected,
    double complex *pfMNew, const double complex *oldPfM,
    const double complex *oldInvM, const size_t invMQpStride,
    const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater, int *failureDetail);
int CalculateNewPfMBF_fsz_rows_workspace(
    const int nAffected, const int *affected,
    double complex *pfMNew, const double complex *oldPfM,
    const double complex *oldInvM, const size_t invMQpStride,
    const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater, int *failureDetail,
    double complex *complexWork, const size_t complexWorkCount,
    int *iwork, const size_t intWorkCount,
    double *rwork, const size_t rworkCount);

/* F3-b transactional accepted inverse prepare. This computes rank-local
   candidate inverse matrices into invMNew without modifying global state. */
int GetPrepareInvMBF_fsz_rows_work_size(size_t *complexCount,
    size_t *intCount);
BF_FSZ_InvUpdateResult PrepareInvMBF_fsz_rows_workspace(
    const int nAffected, const int *affected,
    const double complex *oldInvM, const size_t oldInvMQpStride,
    double complex *invMNew, const size_t newInvMQpStride,
    const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater,
    double complex *complexWork, const size_t complexWorkCount,
    int *iwork, const size_t intWorkCount);

#endif
