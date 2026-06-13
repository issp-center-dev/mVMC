#pragma once
/* Sector projection helpers for doublon-only (NExUpdatePath=6) and t-J (4/5).
 * Requires NExUpdatePath and Nsite from global.h to be in scope.
 * Include this file from locgrn.c and locgrn_real.c (after global.h is in scope).
 * lslocgrn.c and lslocgrn_real.c share the same translation unit, so the helpers
 * defined here are visible to them without a separate include. */

/* Pre-hop 1-hop sector check: returns 1 if c†_{ri,s} a_{rj,s} preserves the sector.
 * Called at hop-BEFORE position (eleNum not yet modified).
 * Returns 1 (no-op) when NExUpdatePath is not in {4,5,6}. */
static inline int IsSectorPreserved_1hopPre(int ri, int rj, int s, const int *eleNum) {
  if (NExUpdatePath == 6)                        /* doublon-only: hop away breaks doublon at rj */
    return eleNum[rj + (1-s)*Nsite] != 1;
  if (NExUpdatePath == 4 || NExUpdatePath == 5)  /* t-J: hop into ri must not create doublon */
    return eleNum[ri + (1-s)*Nsite] != 1;
  return 1;
}

/* Post-hop sector check: returns 1 if current eleNum is a valid sector state.
 * Called at hop-AFTER position (eleNum already modified). Read-only.
 * Returns 1 (no-op) when NExUpdatePath is not in {4,5,6}. */
static inline int IsSectorStateAllowed(const int *eleNum) {
  int r;
  if (NExUpdatePath == 6) {                        /* doublon-only: every site empty or doublon */
    for (r = 0; r < Nsite; r++)
      if (eleNum[r] != eleNum[r + Nsite]) return 0;
  } else if (NExUpdatePath == 4 || NExUpdatePath == 5) { /* t-J: no doubly-occupied site */
    for (r = 0; r < Nsite; r++)
      if (eleNum[r] * eleNum[r + Nsite] == 1) return 0;
  }
  return 1;
}

/* Pre-hop 2-hop sector check: returns 1 if applying c†_{rsi} a_{rsj} c†_{rsk} a_{rsl}
 * would keep eleNum inside the sector. Uses spin-site indices rsi=ri+s*Nsite etc.
 * Temporarily mutates eleNum and reverts — must be called from single-threaded context.
 * Called from calHCACA/calHCACA_real at hop-BEFORE position. */
static inline int WouldPreserve_2hopPre(int rsi, int rsj, int rsk, int rsl, int *eleNum) {
  int ok;
  eleNum[rsl] = 0; eleNum[rsk] = 1;
  eleNum[rsj] = 0; eleNum[rsi] = 1;
  ok = IsSectorStateAllowed(eleNum);
  eleNum[rsl] = 1; eleNum[rsk] = 0;
  eleNum[rsj] = 1; eleNum[rsi] = 0;
  return ok;
}
