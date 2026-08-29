#ifndef _SLATER
#define _SLATER

#define BF_FSZ_ROW_BUILD_OK 0
#define BF_FSZ_ROW_BUILD_NEEDS_FULL 1
#define BF_FSZ_ROW_BUILD_INVALID_ARGUMENT 2
void UpdateSlaterElm_fcmp();
void SlaterElmDiff_fcmp(double complex *srOptO, const double complex ip, int *eleIdx);

void SlaterElmBFDiff_fcmp(double complex*srOptO, const double complex ip, int *eleIdx, int *eleNum, int *eleCfg, int *eleProjConst,const int * eleProjBFCnt);

void BackFlowDiff_fcmp(complex double *srOptO, const double complex ip, int *eleIdx, const int *eleNum, int *eleProjConst,
                       const int *eleProjBFCnt);
void SlaterElmBFDiff_fsz(double complex *srOptO, const double complex ip,
                         const int *eleIdx, const int *eleSpn,
                         const int *eleNum, const int *eleProjBFCnt);
void BackFlowDiff_fsz(double complex *srOptO, const double complex ip,
                      const int *eleIdx, const int *eleSpn,
                      const int *eleNum, const int *eleProjBFCnt);

void MakeSlaterElmBF_fcmp(const int *eleNum, const int *eleProjBFCnt);
void MakeSlaterElmBF_fcmp_to_serial(double complex *sltElmBF,
                                    const int *eleNum,
                                    const int *eleProjBFCnt);
int RebuildSlaterMAllBF_fcmp(
    const int *eleIdx, const int *eleNum, const int *eleProjBFCnt,
    int qpStart, int qpEnd, double complex *sltElmBF,
    double complex *pfMOut, double complex *invMOut);
int RebuildSlaterMAllBF_real(
    const int *eleIdx, const int *eleNum, const int *eleProjBFCnt,
    int qpStart, int qpEnd, double complex *sltElmBF,
    double complex *pfMComplexOut, double complex *invMComplexOut,
    double *sltElmBFReal, double *pfMRealOut, double *invMRealOut);
int CalculateBFCanonicalPf_fcmp(
    const int *eleIdx, const int *eleNum, const int *eleProjBFCnt,
    int qpStart, int qpEnd, double complex *pfMOut);
int CalculateBFCanonicalPf_real(
    const int *eleIdx, const int *eleNum, const int *eleProjBFCnt,
    int qpStart, int qpEnd, double *pfMOut);
void MakeSlaterElmBF_fsz(const int *eleNum, const int *eleProjBFCnt);
void MakeSlaterElmBF_fsz_to(double complex *sltElmBF, const int *eleNum, const int *eleProjBFCnt);
void MakeSlaterElmBF_fsz_to_serial(double complex *sltElmBF, const int *eleNum,
                       const int *eleProjBFCnt);
void MakeSlaterElmBF_fsz_hop_to(double complex *sltElmBF, const double complex *baseSltElmBF,
                                const int *oldEleProjBFCnt, const int *newEleProjBFCnt);
void MakeSlaterElmBF_fsz_hop_to_serial(double complex *sltElmBF,
                                       const double complex *baseSltElmBF,
                                       const int *oldEleProjBFCnt,
                                       const int *newEleProjBFCnt);
int MakeSlaterElmBF_fsz_hop_to_with_rows(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected);
int MakeSlaterElmBF_fsz_hop_to_with_rows_serial(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected);
/* The returned size excludes affected[Nsize]. The affected and intWork ranges
 * must not overlap; sltElmBF and baseSltElmBF must not alias. */
int GetSlaterElmBF_fsz_hop_int_work_size(int *workSize);
int MakeSlaterElmBF_fsz_hop_to_with_rows_workspace(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize);
int MakeSlaterElmBF_fsz_hop_to_with_rows_workspace_serial(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize);
int GetSlaterElmBF_fsz_hop_row_work_size(size_t *complexCount,
    int qpNum, int rowCapacity);
/* Multi-move candidate rows use the same layout and status contract as the
 * one-move hop API. movedParticles must contain nMoved distinct particles.
 * candidateRows must not overlap the base matrix, integer workspace,
 * affected output, moved list, or particle/projection inputs. INVALID_ARGUMENT
 * and NEEDS_FULL leave candidateRows unchanged; NEEDS_FULL still returns the
 * complete affected list and counts. A non-NULL candidateRows buffer and its
 * QP stride must hold rowCapacity rows per local QP. */
int MakeSlaterElmBF_fsz_multi_move_rows_workspace(
    double complex *candidateRows, size_t rowQpStride,
    size_t rowAffectedStride, int rowCapacity,
    const double complex *baseSltElmBF,
    int nMoved, const int *movedParticles,
    const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int qpStart, int qpEnd,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize);
int MakeSlaterElmBF_fsz_multi_move_rows_workspace_serial(
    double complex *candidateRows, size_t rowQpStride,
    size_t rowAffectedStride, int rowCapacity,
    const double complex *baseSltElmBF,
    int nMoved, const int *movedParticles,
    const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int qpStart, int qpEnd,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize);
int MakeSlaterElmBF_fsz_hop_rows_workspace(
    double complex *candidateRows, size_t rowQpStride,
    size_t rowAffectedStride, int rowCapacity,
    const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int qpStart, int qpEnd,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize);
int MakeSlaterElmBF_fsz_hop_rows_workspace_serial(
    double complex *candidateRows, size_t rowQpStride,
    size_t rowAffectedStride, int rowCapacity,
    const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int qpStart, int qpEnd,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize);
int CommitSlaterElmBF_fsz_hop_workspace(
    double complex *sltElmBF,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *intWork, int intWorkSize);
int CommitSlaterElmBF_fsz_hop_workspace_serial(
    double complex *sltElmBF,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *intWork, int intWorkSize);
extern long long BFRowSelectAdds;
extern long long BFRowSelectLowCountAdds;
extern int BFRowSelectCountersEnabled;
void UpdateSlaterElmBF_fcmp(const int ma, const int ra, const int rb, const int u,
                       const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double complex*sltElmTmp);
void UpdateSlaterElmBF_real(const int ma, const int ra, const int rb, const int u,
                       const int *eleCfg, const int *eleNum, const int *eleProjBFCntOld,
                       const int *eleProjBFCnt, int *msa, int *hopNum, double *sltElmTmp);
void UpdateSlaterElmBFGrn(const int ma, const int ra, const int rb, const int u,
                          const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double complex*sltElmTmp);

void UpdateSlaterElmBFGrn_real(const int ma, const int ra, const int rb, const int u,
                          const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double *sltElmTmp,
                          int *restoreRows, int *restoreRowCount);
void UpdateSlaterElmBFGrnVec_real(const int ma, const int ra, const int rb, const int u,
                          const int *eleIdx, const int *eleCfg, const int *eleNum,
                          const int *eleProjBFCntOld, const int *eleProjBFCnt,
                          int *msa, int *hopNum, double *vecTmp);

void StoreSlaterElmBF_fcmp(complex double *bufM);

void StoreSlaterElmBF_real(double *bufM_real);

void StoreSlaterElmBF_real_rows(double *bufM_real, const int *restoreRows, const int *restoreRowCount);

#endif
