#ifndef _VMCMAKE_GC_H
#define _VMCMAKE_GC_H

#include <complex.h>
#include <mpi.h>

enum GCMoveClass { GC_MOVE_HOP, GC_MOVE_ADD, GC_MOVE_REMOVE };

void VMCMakeSampleGC(MPI_Comm comm);
int makeInitialSampleGC(int *eleIdx, int *eleCfg, int *eleNum,
                        int *eleProjCnt, int qpStart, int qpEnd,
                        MPI_Comm comm);
void copyFromBurnSampleGC(int *eleIdx, int *eleCfg, int *eleNum,
                          int *eleProjCnt);
void copyToBurnSampleGC(const int *eleIdx, const int *eleCfg,
                        const int *eleNum, const int *eleProjCnt);
void saveEleConfigGC(int sample, double complex logIp, const int *eleIdx,
                     const int *eleCfg, const int *eleNum,
                     const int *eleProjCnt, int ncur);
enum GCMoveClass GCSelectMoveClass(double classDraw, int ncur, int nsite2);
int GCAttemptMove(enum GCMoveClass moveClass, int arg0, int arg1,
                  double acceptDraw, int *eleIdx, int *eleCfg, int *eleNum,
                  int *eleProjCnt, double complex *logIpOld,
                  double complex *pfMNew, int *projCntNew, int qpStart,
                  int qpEnd, MPI_Comm comm);
int GCMakeOneStep(int *eleIdx, int *eleCfg, int *eleNum,
                  int *eleProjCnt, double complex *logIpOld,
                  double complex *pfMNew, int *projCntNew, int qpStart,
                  int qpEnd, MPI_Comm comm);

#endif
