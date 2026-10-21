#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "include/gc_antiparallel.h"

int GCAntiEnabled(const int grandCanonical, const int orbitalGeneral) {
  return grandCanonical != 0 && orbitalGeneral == 0;
}

int GCAntiSelectBridge(const double draw) {
  if (draw < 1.0 / 64.0) return GC_ANTI_GLOBAL;
  if (draw < 9.0 / 64.0) return GC_ANTI_ADD_TWO;
  if (draw < 17.0 / 64.0) return GC_ANTI_REMOVE_TWO;
  return -1;
}

int GCAntiSelectFiniteSpinBridge(const double draw) {
  if (draw < 1.0 / 64.0) return GC_ANTI_GLOBAL;
  if (draw < 9.0 / 64.0) return GC_ANTI_ADD_TWO;
  if (draw < 17.0 / 64.0) return GC_ANTI_REMOVE_TWO;
  if (draw < 25.0 / 64.0) return GC_ANTI_SWAP;
  return -1;
}

double GCAntiRatioTwoPairs(const int m, const int nsite, const int adding) {
  double ratio;
  if (nsite < 2 || m < 0 || m > nsite) return 0.0;
  if (adding) {
    if (nsite - m < 2) return 0.0;
    ratio = ((double)(nsite - m) / (double)(m + 2)) *
            ((double)(nsite - m - 1) / (double)(m + 1));
  } else {
    if (m < 2) return 0.0;
    ratio = ((double)m / (double)(nsite - m + 2)) *
            ((double)(m - 1) / (double)(nsite - m + 1));
  }
  return ratio * ratio;
}

int GCAntiDrawBridge(const enum GCAntiBridgeMove move, const int *current,
                     const int nsite, const int m, int *candidate,
                     uint32_t (*next32)(void), double *ratio) {
  int spin;
  if (nsite <= 0 || nsite > INT_MAX / 2 || m < 0 || m > nsite ||
      current == NULL || candidate == NULL || current == candidate ||
      next32 == NULL || ratio == NULL) return -1;
  if (move == GC_ANTI_GLOBAL) {
    /* Choose L labels uniformly among U_i and Dbar_i. Every balanced
     * configuration has exactly one such subset, so q(x->y) is constant. */
    int label, remaining = nsite;
    for (label = 0; label < 2 * nsite; ++label) {
      const int taken = (int)GCAntiDrawBelow((uint32_t)(2 * nsite - label),
                                            next32) < remaining;
      remaining -= taken;
      candidate[label] = label < nsite ? taken : !taken;
    }
    *ratio = 1.0;
    return remaining == 0 ? 1 : -1;
  }
  if (move == GC_ANTI_SWAP) {
    int site, ups = 0, downs = 0, upPick, downPick, upSite = -1, downSite = -1;
    for (site = 0; site < nsite; ++site) {
      ups += current[site] && !current[site + nsite];
      downs += !current[site] && current[site + nsite];
    }
    *ratio = 1.0;
    if (ups == 0 || downs == 0) return 0;
    upPick = (int)GCAntiDrawBelow((uint32_t)ups, next32);
    downPick = (int)GCAntiDrawBelow((uint32_t)downs, next32);
    for (site = 0; site < nsite; ++site) {
      if (current[site] && !current[site + nsite] && upPick-- == 0)
        upSite = site;
      if (!current[site] && current[site + nsite] && downPick-- == 0)
        downSite = site;
    }
    if (upSite < 0 || downSite < 0) return -1;
    memcpy(candidate, current, (size_t)(2 * nsite) * sizeof(*candidate));
    candidate[upSite] = 0;
    candidate[upSite + nsite] = 1;
    candidate[downSite] = 1;
    candidate[downSite + nsite] = 0;
    return 1;
  }
  if (move != GC_ANTI_ADD_TWO && move != GC_ANTI_REMOVE_TWO) return -1;
  *ratio = GCAntiRatioTwoPairs(m, nsite, move == GC_ANTI_ADD_TWO);
  if (*ratio == 0.0) return 0;
  memcpy(candidate, current, (size_t)(2 * nsite) * sizeof(*candidate));
  for (spin = 0; spin < 2; ++spin) {
    const int occupied = move == GC_ANTI_REMOVE_TWO;
    const int count = occupied ? m : nsite - m;
    int first = (int)GCAntiDrawBelow((uint32_t)count, next32);
    int second = (int)GCAntiDrawBelow((uint32_t)(count - 1), next32);
    int site, k = 0, changed = 0;
    if (second >= first) ++second;
    for (site = 0; site < nsite; ++site) {
      const int rs = spin * nsite + site;
      if (current[rs] == occupied) {
        if (k == first || k == second) {
          candidate[rs] = !occupied;
          ++changed;
        }
        ++k;
      }
    }
    if (k != count || changed != 2) return -1;
  }
  return 1;
}

/* Mark the k-th site whose up and down orbitals are both empty. */
static int GCAntiMarkKthEmptySite(int *eleNum, const int nsite, int k,
                                  const int up, const int down) {
  int site;
  for (site = 0; site < nsite; ++site) {
    if (eleNum[site] == 0 && eleNum[site + nsite] == 0 && k-- == 0) {
      eleNum[site] = up;
      eleNum[site + nsite] = down;
      return 1;
    }
  }
  return 0;
}

int GCAntiDrawSpinConfiguration(const int nsite, const int m, const int spin,
                                uint32_t (*next32)(void), int *eleNum) {
  int low, high, doubles, singles, k, used = 0;
  if (nsite <= 0 || nsite > INT_MAX / 2 || m < 0 || m > nsite ||
      spin < 0 || next32 == NULL || eleNum == NULL)
    return -1;
  low = 2 * m - nsite > 0 ? 2 * m - nsite : 0;
  high = m - spin;
  if (low > high) return 0;
  doubles = low + (int)GCAntiDrawBelow((uint32_t)(high - low + 1), next32);
  singles = m - doubles;
  memset(eleNum, 0, (size_t)(2 * nsite) * sizeof(*eleNum));
  /* Sequential uniform choices without replacement give uniform disjoint
   * subsets: first the doubles, then up-only, then down-only sites. */
  for (k = 0; k < doubles + 2 * singles; ++k, ++used) {
    const int up = k < doubles + singles;
    const int down = k < doubles || k >= doubles + singles;
    const uint32_t pick = GCAntiDrawBelow((uint32_t)(nsite - used), next32);
    if (pick == UINT32_MAX ||
        !GCAntiMarkKthEmptySite(eleNum, nsite, (int)pick, up, down))
      return -1;
  }
  return 1;
}

double GCAntiRatioAdd(const int m, const int nsite, const double pAdd,
                      const double pRemoveNext) {
  double ratio;
  if (nsite <= 0 || m < 0 || m >= nsite || !isfinite(pAdd) ||
      !isfinite(pRemoveNext) || pAdd <= 0.0 || pRemoveNext < 0.0) {
    return 0.0;
  }
  /* (L-m)^2 add candidates forward, (m+1)^2 remove candidates backward. */
  ratio = (double)(nsite - m) / (double)(m + 1);
  return (pRemoveNext / pAdd) * ratio * ratio;
}

double GCAntiRatioRemove(const int m, const int nsite, const double pRemove,
                         const double pAddPrev) {
  double ratio;
  if (nsite <= 0 || m <= 0 || m > nsite || !isfinite(pRemove) ||
      !isfinite(pAddPrev) || pRemove <= 0.0 || pAddPrev < 0.0) {
    return 0.0;
  }
  ratio = (double)m / (double)(nsite - m + 1);
  return (pAddPrev / pRemove) * ratio * ratio;
}

uint32_t GCAntiDrawBelow(const uint32_t bound, uint32_t (*next32)(void)) {
  uint32_t threshold;
  uint32_t x;
  if (bound == 0 || next32 == NULL) return UINT32_MAX;
  /* 2^32 mod bound words are discarded so every residue is equally likely. */
  threshold = (uint32_t)(0U - bound) % bound;
  do {
    x = next32();
  } while (x < threshold);
  return x % bound;
}

int GCAntiFindOccupied(const int *eleNum, const int nsite, const int spin,
                       const int k) {
  int site;
  int count = 0;
  if (eleNum == NULL || nsite <= 0 || nsite > INT_MAX / 2 || spin < 0 ||
      spin > 1 || k < 0) {
    return -1;
  }
  for (site = 0; site < nsite; site++) {
    const int rs = site + spin * nsite;
    if (eleNum[rs] != 0) {
      if (count == k) return rs;
      count++;
    }
  }
  return -1;
}

int GCAntiValidateConfig(const int *eleIdx, const int *eleCfg,
                         const int *eleNum, const int nsite,
                         const int ncur) {
  long long nsite2;
  int up = 0;
  int down = 0;
  int rs;
  int position;
  if (eleIdx == NULL || eleCfg == NULL || eleNum == NULL || nsite <= 0) {
    return 0;
  }
  nsite2 = 2LL * (long long)nsite;
  if (nsite2 > INT_MAX || ncur < 0 || ncur % 2 != 0 ||
      (long long)ncur > nsite2) {
    return 0;
  }
  for (rs = 0; rs < (int)nsite2; rs++) {
    if (eleNum[rs] != 0 && eleNum[rs] != 1) return 0;
    if (eleNum[rs] == 1) {
      if (rs < nsite) {
        up++;
      } else {
        down++;
      }
      if (eleCfg[rs] < 0 || eleCfg[rs] >= ncur ||
          eleIdx[eleCfg[rs]] != rs) {
        return 0;
      }
    } else if (eleCfg[rs] != -1) {
      return 0;
    }
  }
  if (up != ncur / 2 || down != ncur / 2) return 0;
  for (position = 0; position < ncur; position++) {
    const int orbital = eleIdx[position];
    if (orbital < 0 || orbital >= (int)nsite2 || eleNum[orbital] != 1 ||
        eleCfg[orbital] != position) {
      return 0;
    }
  }
  return 1;
}

int GCAntiAcceptLog(const double complex oldLog, const double complex newLog,
                    const double logProjRatio, const double proposalRatio,
                    const double draw, int *accepted) {
  long double logWeight;
  if (accepted == NULL) return 1;
  *accepted = 0;
  if (!isfinite(creal(oldLog)) || !isfinite(cimag(oldLog)) ||
      !isfinite(cimag(newLog)) || !isfinite(logProjRatio) ||
      !isfinite(proposalRatio) || proposalRatio <= 0.0 || !isfinite(draw) ||
      draw < 0.0 || draw >= 1.0) {
    return 1;
  }
  /* A vanishing candidate amplitude is a legitimate rejection. */
  if (creal(newLog) == -INFINITY) return 0;
  if (!isfinite(creal(newLog))) return 1;
  logWeight = 2.0L * ((long double)logProjRatio +
                      (long double)creal(newLog) -
                      (long double)creal(oldLog)) +
              logl((long double)proposalRatio);
  if (isnan(logWeight)) return 1;
  *accepted = logWeight >= 0.0L || draw == 0.0 ||
              logl((long double)draw) < logWeight;
  return 0;
}
