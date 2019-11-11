#ifndef _VMCMAKE_FSZ
#define _VMCMAKE_FSZ
#include <complex.h>
#include <mpi.h>

void VMCMakeSample_fsz(MPI_Comm comm);
int makeInitialSample_fsz(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,int *eleSpn,
                      const int qpStart, const int qpEnd, MPI_Comm comm);
void copyFromBurnSample_fsz(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,int *eleSpn);
void copyToBurnSample_fsz(const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt,const int *eleSpn);
void saveEleConfig_fsz(const int sample, const double complex logIp,
                   const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt,const int *eleSpn);
//void sortEleConfig(int *eleIdx, int *eleCfg, const int *eleNum);
void makeCandidate_hopping_fsz(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn);
void makeCandidate_hopping_csz(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn);

void makeCandidate_exchange_fsz(int *mi_, int *ri_, int *rj_, int *s_, int *rejectFlag_,
                            const int *eleIdx, const int *eleCfg, const int *eleNum,const int *eleSpn);
void updateEleConfig_fsz(int mi, int org_r, int dst_r, int org_spn,int dst_spn,
                     int *eleIdx, int *eleCfg, int *eleNum, int *eleSpn) ;
void revertEleConfig_fsz(int mi, int org_ri, int dst_r, int org_spn,int dst_spn,
                     int *eleIdx, int *eleCfg, int *eleNum,int *eleSpn);
void CheckEleConfig_fsz(int *eleIdx, int *eleCfg, int *eleNum,int *eleSpn,MPI_Comm comm);
int CheckEleNum_fsz(int *eleIdx, int *eleCfg, int *eleNum,int *eleSpn,MPI_Comm comm);
void makeCandidate_LocalSpinFlip_localspin(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn);
void makeCandidate_LocalSpinFlip_conduction(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn);


#endif


