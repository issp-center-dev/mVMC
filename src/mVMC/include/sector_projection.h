#ifndef MVMC_SECTOR_PROJECTION_H
#define MVMC_SECTOR_PROJECTION_H

/* Sector projection helpers for doublon-only (NExUpdatePath=6) and t-J (4/5).
 * NExUpdatePath and Nsite must be in scope before this header is included. */

/* Pre-hop 1-hop sector check. Called before eleNum is modified. */
static inline int IsSectorPreserved_1hopPre(int ri, int rj, int s,
                                            const int *eleNum) {
  if (NExUpdatePath == 6)
    return eleNum[rj + (1 - s) * Nsite] != 1;
  if (NExUpdatePath == 4 || NExUpdatePath == 5)
    return eleNum[ri + (1 - s) * Nsite] != 1;
  return 1;
}

/* Post-hop sector check. Called after eleNum is modified and is read-only. */
static inline int IsSectorStateAllowed(const int *eleNum) {
  int r;
  if (NExUpdatePath == 6) {
    for (r = 0; r < Nsite; r++) {
      if (eleNum[r] != eleNum[r + Nsite]) return 0;
    }
  } else if (NExUpdatePath == 4 || NExUpdatePath == 5) {
    for (r = 0; r < Nsite; r++) {
      if (eleNum[r] * eleNum[r + Nsite] == 1) return 0;
    }
  }
  return 1;
}

/* Pre-hop 2-hop sector check. Temporarily mutates eleNum and restores it. */
static inline int WouldPreserve_2hopPre(int rsi, int rsj, int rsk, int rsl,
                                        int *eleNum) {
  int ok;
  eleNum[rsl] = 0;
  eleNum[rsk] = 1;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  ok = IsSectorStateAllowed(eleNum);
  eleNum[rsl] = 1;
  eleNum[rsk] = 0;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;
  return ok;
}

#endif
