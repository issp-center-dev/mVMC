#ifndef _GC_CONFIG_H
#define _GC_CONFIG_H

void GCAddPair(int rsa, int rsb, int *eleIdx, int *eleCfg, int *eleNum,
               int *ncur);
int GCRemovePair(int pos0, int pos1, int *eleIdx, int *eleCfg, int *eleNum,
                 int *ncur);
void GCRestoreRemovedPair(int pos0, int pos1, int *eleIdx, int *eleCfg,
                          int *eleNum, int *ncur);
int GCRemoveParitySign(int pos0, int pos1, int ncurOld);
double GCProposalRatioAdd(int ncurOld, int nsite2, double pAddAvailX,
                          double pRemoveAvailY);
double GCProposalRatioRemove(int ncurOld, int nsite2, double pRemoveAvailX,
                             double pAddAvailY);
int GCFindKthEmpty(const int *eleNum, int nsite2, int k);

#endif
