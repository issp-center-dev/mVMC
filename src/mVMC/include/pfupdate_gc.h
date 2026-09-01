#ifndef _PFUPDATE_GC_H
#define _PFUPDATE_GC_H

#include <complex.h>

/* Candidate functions are called after the proposal has temporarily changed
 * eleIdx (and the matching eleCfg/eleNum state), but before shared InvM/PfM
 * are updated. They read the old inverse/Pfaffian and the new ordering. */
void CalculateNewPfMHopGC(int ma, double complex *pfMNew,
                          const int *eleIdx, int ncur, int qpStart,
                          int qpEnd);
void UpdateMAllHopGC(int ma, const int *eleIdx, int ncur, int qpStart,
                     int qpEnd);
void CalculateNewPfMTwoHopGC(int ma, int mb, double complex *pfMNew,
                             const int *eleIdx, int ncur, int qpStart,
                             int qpEnd);
void CalculateNewPfMTwoHopGCWorkspace(
    int ma, int mb, double complex *pfMNew, const int *eleIdx, int ncur,
    int qpStart, int qpEnd, double complex *vecA, double complex *vecB);

/* n-row candidate used by GC Green functions. msa[k] is the moved particle
 * position and rsa[k] its destination fused orbital. All scratch belongs to
 * the caller. vec has n*NsizeMax elements, w has NsizeMax, smallMat/pfRWork
 * have (2*n)^2, pfIWork has 2*n, and pfLWork must be at least (2*n)^2. */
double complex CalculateNewPfMNGC(
    int qpidx, int n, const int *msa, const int *rsa, const int *eleIdx,
    int ncur, double complex *vec, double complex *w,
    double complex *smallMat, int *pfIWork, double complex *pfWork,
    int pfLWork, double *pfRWork);

void CalculateNewPfMAddGC(int rsa, int rsb, double complex *pfMNew,
                          const int *eleIdx, int ncurOld, int qpStart,
                          int qpEnd);
void UpdateMAllAddGC(int rsa, int rsb, const int *eleIdx, int ncurOld,
                     int qpStart, int qpEnd);
void CalculateNewPfMRemoveGC(int pos0, int pos1,
                             double complex *pfMNew, const int *eleIdx,
                             int ncurOld, int qpStart, int qpEnd);
void UpdateMAllRemoveGC(int pos0, int pos1, const int *eleIdx,
                        int ncurOld, int qpStart, int qpEnd);

#endif
