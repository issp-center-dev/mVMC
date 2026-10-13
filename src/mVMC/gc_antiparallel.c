#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "include/gc_antiparallel.h"

int GCAntiEnabled(const int grandCanonical, const int orbitalGeneral) {
  return grandCanonical != 0 && orbitalGeneral == 0;
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
