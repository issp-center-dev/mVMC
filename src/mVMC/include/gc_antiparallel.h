#ifndef _GC_ANTIPARALLEL_H
#define _GC_ANTIPARALLEL_H

#include <complex.h>
#include <stdint.h>

/* Grand-canonical sampling of a pair wave function built only from
 * anti-parallel (up-down) pairs.  These helpers do not touch MPI or global
 * state so that the proposal can be verified in isolation. */

/* Nonzero when NGrandCanonical is set and no OrbitalGeneral/Parallel input
 * is present (the orbital table is the Nsite x Nsite anti-parallel form). */
int GCAntiEnabled(int grandCanonical, int orbitalGeneral);

/* Reverse/forward proposal ratio for adding one up-down pair to a state with
 * m pairs on nsite sites.  pAdd is the add-class probability at m and
 * pRemoveNext the remove-class probability at m+1.  Invalid input gives 0. */
double GCAntiRatioAdd(int m, int nsite, double pAdd, double pRemoveNext);

/* Reverse/forward ratio for removing one up-down pair from m pairs. */
double GCAntiRatioRemove(int m, int nsite, double pRemove, double pAddPrev);

/* Unbiased draw in [0, bound) by rejection; UINT32_MAX on invalid input. */
uint32_t GCAntiDrawBelow(uint32_t bound, uint32_t (*next32)(void));

/* Fused orbital index (site + spin*nsite) of the k-th occupied site of the
 * given spin in ascending site order, or -1. */
int GCAntiFindOccupied(const int *eleNum, int nsite, int spin, int k);

/* 1 when the first ncur entries of eleIdx, eleCfg and eleNum describe a
 * consistent state with equal up and down counts, 0 otherwise.  Entries of
 * eleIdx beyond ncur are scratch and are not inspected. */
int GCAntiValidateConfig(const int *eleIdx, const int *eleCfg,
                         const int *eleNum, int nsite, int ncur);

/* Metropolis decision in log space for weight
 *   exp(2*(logProjRatio + Re newLog - Re oldLog)) * proposalRatio.
 * Returns 0 on success with *accepted set; a zero candidate (Re newLog=-Inf)
 * is a valid rejection.  Nonfinite or out-of-range input returns 1. */
int GCAntiAcceptLog(double complex oldLog, double complex newLog,
                    double logProjRatio, double proposalRatio, double draw,
                    int *accepted);

#endif
