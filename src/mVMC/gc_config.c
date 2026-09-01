#include "include/gc_config.h"

static void GCSwapOccupiedPositions(int pos0, int pos1, int *eleIdx,
                                    int *eleCfg) {
  const int rs0 = eleIdx[pos0];
  const int rs1 = eleIdx[pos1];
  eleIdx[pos0] = rs1;
  eleIdx[pos1] = rs0;
  eleCfg[rs1] = pos0;
  eleCfg[rs0] = pos1;
}

void GCAddPair(const int rsa, const int rsb, int *eleIdx, int *eleCfg,
               int *eleNum, int *ncur) {
  const int pos0 = *ncur;
  const int pos1 = pos0 + 1;
  eleIdx[pos0] = rsa;
  eleIdx[pos1] = rsb;
  eleCfg[rsa] = pos0;
  eleCfg[rsb] = pos1;
  eleNum[rsa] = 1;
  eleNum[rsb] = 1;
  *ncur += 2;
}

int GCRemoveParitySign(const int pos0, const int pos1, const int ncurOld) {
  const int nswap = (pos1 != ncurOld - 1) + (pos0 != ncurOld - 2);
  return (nswap % 2 == 0) ? 1 : -1;
}

int GCRemovePair(const int pos0, const int pos1, int *eleIdx, int *eleCfg,
                 int *eleNum, int *ncur) {
  const int ncurOld = *ncur;
  const int parity = GCRemoveParitySign(pos0, pos1, ncurOld);
  int removed0;
  int removed1;

  if (pos1 != ncurOld - 1) {
    GCSwapOccupiedPositions(pos1, ncurOld - 1, eleIdx, eleCfg);
  }
  if (pos0 != ncurOld - 2) {
    GCSwapOccupiedPositions(pos0, ncurOld - 2, eleIdx, eleCfg);
  }
  removed0 = eleIdx[ncurOld - 2];
  removed1 = eleIdx[ncurOld - 1];
  eleCfg[removed0] = -1;
  eleCfg[removed1] = -1;
  eleNum[removed0] = 0;
  eleNum[removed1] = 0;
  *ncur = ncurOld - 2;
  return parity;
}

void GCRestoreRemovedPair(const int pos0, const int pos1, int *eleIdx,
                          int *eleCfg, int *eleNum, int *ncur) {
  const int ncurOld = *ncur + 2;
  int i;
  if (pos0 != ncurOld - 2) {
    const int temporary = eleIdx[pos0];
    eleIdx[pos0] = eleIdx[ncurOld - 2];
    eleIdx[ncurOld - 2] = temporary;
  }
  if (pos1 != ncurOld - 1) {
    const int temporary = eleIdx[pos1];
    eleIdx[pos1] = eleIdx[ncurOld - 1];
    eleIdx[ncurOld - 1] = temporary;
  }
  for (i = 0; i < ncurOld; i++) {
    eleCfg[eleIdx[i]] = i;
    eleNum[eleIdx[i]] = 1;
  }
  *ncur = ncurOld;
}

static double GCChooseTwo(const int count) {
  return 0.5 * (double)count * (double)(count - 1);
}

double GCProposalRatioAdd(const int ncurOld, const int nsite2,
                          const double pAddAvailX,
                          const double pRemoveAvailY) {
  const int emptyOld = nsite2 - ncurOld;
  const int occupiedNew = ncurOld + 2;
  if (emptyOld < 2 || occupiedNew < 2 || pAddAvailX <= 0.0 ||
      pRemoveAvailY < 0.0) {
    return 0.0;
  }
  return (pRemoveAvailY / pAddAvailX) *
         (GCChooseTwo(emptyOld) / GCChooseTwo(occupiedNew));
}

double GCProposalRatioRemove(const int ncurOld, const int nsite2,
                             const double pRemoveAvailX,
                             const double pAddAvailY) {
  const int emptyNew = nsite2 - ncurOld + 2;
  if (ncurOld < 2 || emptyNew < 2 || pRemoveAvailX <= 0.0 ||
      pAddAvailY < 0.0) {
    return 0.0;
  }
  return (pAddAvailY / pRemoveAvailX) *
         (GCChooseTwo(ncurOld) / GCChooseTwo(emptyNew));
}

int GCFindKthEmpty(const int *eleNum, const int nsite2, const int k) {
  int emptyIndex = 0;
  int rs;
  if (eleNum == 0 || nsite2 < 0 || k < 0) return -1;
  for (rs = 0; rs < nsite2; rs++) {
    if (eleNum[rs] == 0) {
      if (emptyIndex == k) return rs;
      emptyIndex++;
    }
  }
  return -1;
}
