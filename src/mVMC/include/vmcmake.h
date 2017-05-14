#ifndef _VMCMAKE
#define _VMCMAKE
#include <complex.h>
#include <mpi.h>

void VMCMakeSample(MPI_Comm comm);
int makeInitialSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                      const int qpStart, const int qpEnd, MPI_Comm comm);
void copyFromBurnSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt);
void copyToBurnSample(const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt);
void saveEleConfig(const int sample, const double complex logIp,
                   const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt);
void sortEleConfig(int *eleIdx, int *eleCfg, const int *eleNum);
void ReduceCounter(MPI_Comm comm);
void makeCandidate_hopping(int *mi_, int *ri_, int *rj_, int *s_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg);
void makeCandidate_exchange(int *mi_, int *ri_, int *rj_, int *s_, int *rejectFlag_,
                            const int *eleIdx, const int *eleCfg, const int *eleNum);
void updateEleConfig(int mi, int ri, int rj, int s,
                     int *eleIdx, int *eleCfg, int *eleNum);
void revertEleConfig(int mi, int ri, int rj, int s,
                     int *eleIdx, int *eleCfg, int *eleNum);


/*[s] BackFlow */
void VMC_BF_MakeSample(MPI_Comm comm);
int makeInitialSampleBF(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, int *eleProjBFCnt,
                        const int qpStart, const int qpEnd, MPI_Comm comm);
void copyFromBurnSampleBF(int *eleIdx);
void copyToBurnSampleBF(const int *eleIdx);
void saveEleConfigBF(const int sample, const double logIp,
                     const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt);
/*[e] BackFlow */

typedef enum {HOPPING,HOPPING_FSZ,EXCHANGE,LOCALSPINFLIP, NONE} UpdateType;
UpdateType getUpdateType(int path);

#endif


