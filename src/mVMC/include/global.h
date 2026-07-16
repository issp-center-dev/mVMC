/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
*/
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * global variables
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#ifndef _INCLUDE_GLOBAL
#define _INCLUDE_GLOBAL
#include <complex.h>
#include "stdio.h"
#define D_FileNameMax 256

/***** definition *****/
char CDataFileHead[D_FileNameMax]; /* prefix of output files */
char CParaFileHead[D_FileNameMax]; /* prefix for optimized variational parameters */

int NVMCCalMode; /* calculation mode
                    0: optimization of variational paraneters,
                    1: calculation of expectation values */
int NLanczosMode; /* mode of the single Lanczos step
                     0: none, 1: only energy, 2: Green functions */

int NStoreO; /* choice of store O: 0-> normal other-> store  */
int NSRCG; /* choice of solver for Sx=g: 0-> (Sca)LAPACK other-> CG  */
int NSRCGFallback; /* choice of fallback for SR-CG failure: 0-> off, 1-> on */
int NSRCGAbortOnFail; /* abort when SR-CG fails: 0-> warn and continue, 1-> abort */
int reweight; /* 1: reweight in vmccal.c, vmccal_fsz.c, other: no reweight, default 0   */

int NDataIdxStart; /* starting value of the file index */
int NDataQtySmp; /* the number of output files */

int Nsite; /* the number of sites */
int Ne;    /* the number of electrons with up spin */
int Nup;   /* the number of electrons with up spin */
int Nsize; /* the number of electrons = 2*Ne */
int Nsite2; /* 2*Nsite */
int NzBF; /* BF connectivity */
int TwoSz;

int NSPGaussLeg; /* the number of points for the Gauss-Legendre quadrature */
int NSPStot; /* S of Spin projection */
int NMPTrans; /* the number of quantum projection for translation and point group symmetry */
int NQPFull; /* the total number of quantum projection = NSPGaussLeg*NMPTrans*NQPTransOpt */
int NQPFix; /* for QPFixWeight NSPGaussLeg*NMPTranss */

int NSROptItrStep; /* the number of SR method steps */
int NSROptItrSmp; /* the number of SR method steps for calculation of average value */
int NSROptFixSmp; /* the number of SR method steps with fixed samples (1 is recommended) */

double DSROptRedCut; /* SR stabilizing factor for truncation of redundant directions */
double DSROptStaDel; /* SR stabiliaing factor for diagonal element modification */
double DSROptStepDt; /* step width of the SR method */

int NSROptCGMaxIter; /* the number of maximum iterations in SR-CG method */
double DSROptCGTol; /* the tolerance for SR-CG method */

int NVMCWarmUp; /* Monte Carlo steps for warming up */
int NVMCInterval; /* sampling interval [MCS] */
int NVMCSample; /* the number of samples */
int NExUpdatePath; /* update path  0: hopping, 1: hopping+exchange, 2: exchange(spin), 3: KondoGC, 6: pair hopping(doublon-only) */
int NBlockUpdateSize; /* {DEFINED: _pf_block_update} size of block Pfaffian update */

int RndSeed; /* seed for pseudorandom number generator */
int NSplitSize; /* the number of inner MPI processes */

/* total length of def array */
int NTotalDefInt, NTotalDefDouble;

/* zlocspin.def */
int NLocSpn; /* the number of local spin */
int *LocSpn; /* [Nsite] */
/* local spin flag  0: local spin, 1: itinerant electron */

/* for Hamiltonian */
int NTransfer;
int **Transfer; /* [NTransfer][4] */
double complex*ParaTransfer;

int NCoulombIntra;
int *CoulombIntra; /* [NCoulombIntra] */
double *ParaCoulombIntra;

int NCoulombInter;
int **CoulombInter; /* [NCoulombInter][2] */
double *ParaCoulombInter;

int NHundCoupling;
int **HundCoupling; /* [NHundCoupling][2] */
double *ParaHundCoupling;

int NPairHopping;
int **PairHopping; /* [NPairHopping][2] */
double *ParaPairHopping;

int NExchangeCoupling;
int **ExchangeCoupling; /* [NExchangeCoupling][2] */
double *ParaExchangeCoupling;

int NInterAll;
int **InterAll; /* [NInterAll][8] */
double complex*ParaInterAll;
int NNBodyInterAll;
int NBodyInterAllTotalFactors;
int NBodyInterAllMaxN;
int *NBodyInterAllN;        /* [NNBodyInterAll] */
int *NBodyInterAllOffset;   /* [NNBodyInterAll] */
int **NBodyInterAllIdx;     /* [NBodyInterAllTotalFactors][4] */
double complex *ParaNBodyInterAll;

/* for variational parameters */
int NGutzwillerIdx, *GutzwillerIdx; /* [Nsite] */
int NJastrowIdx, **JastrowIdx; /* [Nsite][Nsite] */
int NSpinJastrowIdx, **SpinJastrowIdx; /* [Nsite][Nsite] */
int NDoublonHolon2siteIdx, **DoublonHolon2siteIdx; /* DoublonHolon2siteIdx[idx][2*Nsite] */
int NDoublonHolon4siteIdx, **DoublonHolon4siteIdx; /* DoublonHolon4siteIdx[idx][4*Nsite] */
int NOrbitalIdx, **OrbitalIdx; /* [Nsite][Nsite] */
int **OrbitalSgn; /* OrbitalSgn[Nsite][Nsite] = +1 or -1 */
int iFlgOrbitalGeneral=0;
int iNOrbitalParallel=0;
int iNOrbitalAntiParallel=0;

/* restricted Boltzman Machine for variational parameters */
int Nneuron,NneuronGeneral,NneuronCharge,NneuronSpin;
int NRBM_HiddenLayerIdx,NRBM_PhysLayerIdx,NRBM_PhysHiddenIdx;
int NGeneralRBM_HiddenLayerIdx, *GeneralRBM_HiddenLayerIdx; /* [Nneuron] */
int NGeneralRBM_PhysLayerIdx,   *GeneralRBM_PhysLayerIdx;   /* [Nsite2] */
int NGeneralRBM_PhysHiddenIdx, **GeneralRBM_PhysHiddenIdx;  /* [Nsite2][Nneuron] */
int NChargeRBM_HiddenLayerIdx, *ChargeRBM_HiddenLayerIdx; /* [Nneuron] */
int NChargeRBM_PhysLayerIdx,   *ChargeRBM_PhysLayerIdx;   /* [Nsite] */
int NChargeRBM_PhysHiddenIdx, **ChargeRBM_PhysHiddenIdx;  /* [Nsite][Nneuron] */
int NSpinRBM_HiddenLayerIdx, *SpinRBM_HiddenLayerIdx; /* [Nneuron] */
int NSpinRBM_PhysLayerIdx,   *SpinRBM_PhysLayerIdx;   /* [Nsite] */
int NSpinRBM_PhysHiddenIdx, **SpinRBM_PhysHiddenIdx;  /* [Nsite][Nneuron] */
int NBlockSize_RBMRatio;  /* block size for RBMRatio function. It is Tuning for performance. */

/* zqptransidx.def */
int NQPTrans, **QPTrans, **QPTransInv; /* [NQPTrans][Nsite] */
int **QPTransSgn; /* QPTransSgn[NQPTrans][NSite] = +1 or -1 */
double complex *ParaQPTrans;

/* zqpopttrans.def */
int NQPOptTrans, **QPOptTrans; /* [NQPOptTrans][Nsite] */
int **QPOptTransSgn; /* QPOptTransSgn[NQPOptTrans][NSite] = +1 or -1 */
double *ParaQPOptTrans;

/* for Green functions */
int NCisAjs,         **CisAjsIdx;         /* [NCisAjs][3] */
int NCisAjsCktAlt,   **CisAjsCktAltIdx;   /* [NCisAjsCktAlt][2] */
int NCisAjsCktAltDC, **CisAjsCktAltDCIdx; /* [NCisAjsCktAltDC][6] */
int NNBodyG;
int NBodyGTotalFactors;
int NBodyGMaxN;
int *NBodyGN;        /* [NNBodyG] */
int *NBodyGOffset;   /* [NNBodyG] */
int **NBodyGIdx;     /* [NBodyGTotalFactors][4] */
int **iOneBodyGIdx; /* For GF2 indirect measurement */

/* Optimization flag */
int *OptFlag; /* [NPara]  1: optimized, 0 or 2: fixed */
int AllComplexFlag;/* 0 -> all real variables, !=0-> including complex variables*/

/* flag for RBM */
int FlagRBM=0;

/* flag for anti-periodic boundry condition */
int APFlag; /* 0: periodic, 1: anti-periodic */

/* flag for shift of correlation factors */
/* 0: no shift, 1: shift. Set in ReadDefFileIdxPara(). */
int FlagShiftGJ=0;
int FlagShiftDH2=0;
int FlagShiftDH4=0;

/* flag for OptTrans mode */
int FlagOptTrans=0;
/* flag for Binary mode */
/* output zvo_var.dat (FileVar) as binary data */
int FlagBinary=0;

/* flag for file flush */
int NFileFlushInterval=1;

/***** Variational Parameters *****/
int NPara; /* the total number of variational prameters NPara= NProj + NSlater+ NOptTrans */
int NProj;    /* the number of correlation factor */
int NRBM, NRBM_PhysLayerIdx, NRBM_HiddenLayerIdx;
int NProjBF;    /* the number of correlation factor */
int NSlater;  /* the number of pair orbital (f_ij) = NOrbitalIdx */
int NOptTrans; /* the number of weights for OptTrans. This is used only for variatonal parameters */
               /* NOptTrans = 0 (not OptTrans mode) or NQPOptTrans (OptTrans mode) */
int **etaFlag;   /* Back Flow correlation factor (eta = 1.0 or ProjBF[0])*/
double complex *Para;   /* variatonal parameters */
double complex *Proj;   /* correlation factor (Proj    =Para) */
double complex *RBM;   /*  (Proj    =Para) */
double complex *ProjBF; /* Back flow correlation factor (Proj    =Para) */
double complex *Slater; /* pair orbital       (Slater  =Para+NProj) */
double complex *OptTrans; /* weights          (OptTrans=Para+NProj+NSlater) */
double complex **eta;   /* Back Flow correlation factor (eta = 1.0 or ProjBF[0])*/

/***** Back Flow ******/
int NBackFlowIdx, **BackFlowIdx; /* [Nsite] */
int Nrange, **PosBF, **RangeIdx; /* [Nsite] */
int NBFIdxTotal,NrangeIdx;
int **BFSubIdx; /* [Nsite] */
double *BFRealProj;   /* BFRealProj[NrangeIdx][NrangeIdx] = -creal(ProjBF[BFSubIdx]) */
double *BFRealSlater; /* BFRealSlater[Nsite][Nsite] = creal(Slater[OrbitalIdx]) */
double *BFRealSlaterSign; /* BFRealSlaterSign[Nsite][Nsite] = OrbitalSgn */
double BFRealEta;

/***** Electron Configuration ******/
int *EleIdx;     /* EleIdx[sample][mi+si*Ne] */
int *EleCfg;     /* EleCfg[sample][ri+si*Nsite] */
int *EleNum;     /* EleIdx[sample][ri+si*Nsite] */
int *EleProjCnt; /* EleProjCnt[sample][proj] */
//[s] MERGE BY TM
int *EleSpn;     /* EleIdx[sample][mi+si*Ne] */ //fsz
double complex *RBMCnt;
int *EleProjBFCnt; /* EleProjCnt[sample][proj] */
//[e] MERGE BY TM
double *logSqPfFullSlater; /* logSqPfFullSlater[sample] */
//double complex *SmpSltElmBF; /* logSqPfFullSlater[sample] */
double *SmpSltElmBF_real; /* logSqPfFullSlater[sample] */

int    *SmpEtaFlag; /* logSqPfFullSlater[sample] */
double *SmpEta; /* logSqPfFullSlater[sample] */

int *TmpEleIdx;
int *TmpEleCfg;
int *TmpEleNum;
int *TmpEleProjCnt;
//[s] MERGE BY TM
int *TmpEleSpn;
int *TmpEleProjBFCnt;
//[e] MERGE BY TM
double complex *TmpRBMCnt;

int *BurnEleIdx;
int *BurnEleCfg;
int *BurnEleNum;
int *BurnEleProjCnt;
int *BurnEleSpn;
double complex *BurnRBMCnt;
int BurnFlag=0; /* 0: off, 1: on */

/***** Slater Elements ******/
double complex *SlaterElm; /* SlaterElm[QPidx][ri+si*Nsite][rj+sj*Nsite] */
double complex *InvM; /* InvM[QPidx][mi+si*Ne][mj+sj*Ne] */
double complex *PfM; /* PfM[QPidx] */
// TBC only for real
double *SlaterElm_real; /* SlaterElm[QPidx][ri+si*Nsite][rj+sj*Nsite] */
double *InvM_real; /* InvM[QPidx][mi+si*Ne][mj+sj*Ne] */
double *PfM_real; /* PfM[QPidx] */
double complex *SlaterElmBF; /* SlaterElm[QPidx][ri+si*Nsite][rj+sj*Nsite] */
double complex *InvM; /* InvM[QPidx][mi+si*Ne][mj+sj*Ne] */
double complex *PfM; /* PfM[QPidx] */
// TBC only for real
double *SlaterElm_real; /* SlaterElm[QPidx][ri+si*Nsite][rj+sj*Nsite] */
double *SlaterElmBF_real;/* SlaterElm[QPidx][ri+si*Nsite][rj+sj*Nsite] */
double *InvM_real; /* InvM[QPidx][mi+si*Ne][mj+sj*Ne] */
double *PfM_real; /* PfM[QPidx] */
/***** Quantum Projection *****/
double complex *QPFullWeight; /* QPFullWeight[NQPFull] */
double complex *QPFixWeight; /* QPFixWeight[NQPFix] */
double complex *SPGLCos, *SPGLSin; /* [NSPGaussLeg]  cos(beta/2) and sin(beta/2) */
double complex *SPGLCosSin, *SPGLCosCos, *SPGLSinSin; /* [NSPGaussLeg] */

/***** Stocastic Reconfiguration *****/
int    SROptSize; /* 1+NPara */
double complex *SROptOO; /* [SROptSize*SROptSize] <O^\dagger O> */
double complex *SROptHO; /* [SROptSize]            < HO > */
double complex *SROptO;  /* [SROptSize] calculation buffar */
double complex *SROptO_Store;  /* [SROptSize*NVMCSample] calculation buffer */
//for real
double *SROptOO_real; /* [SROptSize*SROptSize] <O^\dagger O> */ //TBC
double *SROptHO_real; /* [SROptSize]            < HO > */       //TBC
double *SROptO_real;  /* [SROptSize] calculation buffar */      //TBC
double *SROptO_Store_real;  /* [SROptSize*NVMCSample] calculation buffer */

double complex *SROptData; /* [2+NPara] storage for energy and variational parameters */

/***** Physical Quantity *****/
double complex Wc; /* Weight for correlation sampling = <psi|x> */
double complex Etot; /* <H> */
double complex Etot2; /* <H^2> */
double complex Dbtot;
double complex Dbtot2;

double complex *PhysCisAjs; /* [NCisAjs] */
double complex *PhysCisAjsCktAlt; /* [NCisAjsCktAlt] */
double complex *PhysCisAjsCktAltDC; /* [NCisAjsCktAltDC] */
double complex *LocalCisAjs; /* [NCisAjs] */
double complex *LocalCisAjsCktAltDC; /* [NCisAjsCktAltDC] */

double complex Sztot,Sztot2; /* <Sz>,<Sz^2> */


double complex *PhysCisAjs; /* [NCisAjs] */
double complex *PhysCisAjsCktAlt; /* [NCisAjsCktAlt] */
double complex *PhysCisAjsCktAltDC; /* [NCisAjsCktAltDC] */
double complex *LocalCisAjs; /* [NCisAjs] */

/* for Lattice index */
int Nx, Ny, Nz, Norb;
int **LatticeIdx;         /* [Nsite][4] */

/* for Twist operator */
int NTwist, **TwistIdx;         /* TwistIdx -> SiteIdx, SpinIdx */
double **ParaTwist;         /* [NTwist][3*Nsite*2] */
double complex *PhysTwist; /* [NTwist] */
double complex *PhysNBodyG; /* [NNBodyG] */



const int NLSHam = 2; /* 0: I, 1: H */
double complex *QQQQ; /* QQQQ[NLSHam][NLSHam][NLSHam][NLSHam]*/  //TBC
double complex *LSLQ; /* [NLSHam][NLSHam]*/                      //TBC
double *QQQQ_real; /* QQQQ[NLSHam][NLSHam][NLSHam][NLSHam]*/  //TBC
double *LSLQ_real; /* [NLSHam][NLSHam]*/                      //TBC

double complex *QCisAjsQ; /* QCisAjsQ[NLSHam][NLSHam][NCisAjs]*/ //TBC
double complex *QCisAjsCktAltQ; /* QCisAjsCktAltQ[NLSHam][NLSHam][NCisAjsCktAlt]*/ //TBC
double complex *QCisAjsCktAltQDC; /* QCisAjsCktAltQ[NLSHam][NLSHam][NCisAjsCktAlt]
                                     DC Lanczos Calculation */
double complex *LSLCisAjs; /* [NLSHam][NCisAjs]*/                //TBC

double *QCisAjsQ_real; /* QCisAjsQ[NLSHam][NLSHam][NCisAjs]*/ //TBC
double *QCisAjsCktAltQ_real; /* QCisAjsCktAltQ[NLSHam][NLSHam][NCisAjsCktAlt]*/ //TBC
double *QCisAjsCktAltQDC_real; /* QCisAjsCktAltQ[NLSHam][NLSHam][NCisAjsCktAlt]*/
double *LSLCisAjs_real; /* [NLSHam][NCisAjs]*/                //TBC

/***** Output File *****/
/* FILE *FileCfg; */
FILE *FileOut;
FILE *FileVar;
FILE *FileTime;
FILE *FileSRinfo; /* zvo_SRinfo.dat */
FILE *FileCisAjs;
FILE *FileCisAjsCktAlt;
FILE *FileCisAjsCktAltDC;
FILE *FileTwist;
FILE *FileNBodyG;
FILE *FileLS;
FILE *FileLSQQQQ;
FILE *FileLSQCisAjsQ;
FILE *FileLSQCisAjsCktAltQ;
FILE *FileLSQCisAjsCktAltQ;
FILE *FileLSCisAjs;
FILE *FileLSCisAjsCktAlt;
FILE *FileLSCisAjsCktAltDC;


/* FILE *FileTimerList; */
/* FILE *FileOpt;    /\* zqp_opt *\/ */

/***** HitachiTimer *****/
const int NTimer=1000;
double Timer[1000], TimerStart[1000];
double ccc[100];

/* flag for  SROptimization*/
int SRFlag; /* 0: periodic, 1: Diagonalization */

/***** openMP *****/
int NThread;

/***** for DGETRI and DSKPFA in CalculateMAll *****/
int LapackLWork;

/***** counter for vmcMake *****/
int Counter[6] = {0,0,0,0,0,0};
int Counter_max = 6;
/* 0: hopping, 1: hopping accept, 2: exchange try, 3: exchange accept */
/* 4: local spin flip try, 5 local spin flip accept*/

/***** optional BackFlow profiling counters *****/
#define NBFProfileCounter 48
#define BFPROF_SAMPLE_ROW_REQUEST        0
#define BFPROF_SAMPLE_ROW_RECOMPUTE      1
#define BFPROF_SAMPLE_ROW_REUSE          2
#define BFPROF_SAMPLE_PAIR_REQUEST       3
#define BFPROF_SAMPLE_PAIR_RECOMPUTE     4
#define BFPROF_GREEN_ROW_REQUEST         5
#define BFPROF_GREEN_ROW_RECOMPUTE       6
#define BFPROF_GREEN_ROW_REUSE           7
#define BFPROF_GREEN_PAIR_REQUEST        8
#define BFPROF_GREEN_PAIR_RECOMPUTE      9
#define BFPROF_BFCNT_SNAPSHOT           10
#define BFPROF_BFCNT_TOTAL_ENTRY        11
#define BFPROF_BFCNT_GROUP0_NNZ         12
#define BFPROF_BFCNT_GROUP1_NNZ         13
#define BFPROF_BFCNT_GROUP2_NNZ         14
#define BFPROF_BFCNT_GROUP3_NNZ         15
#define BFPROF_BFCNT_STATE0_NNZ         16
#define BFPROF_BFCNT_STATE1_NNZ         17
#define BFPROF_BFCNT_STATE2_NNZ         18
#define BFPROF_BFCNT_STATE3_NNZ         19
#define BFPROF_HOP_TRY                  20
#define BFPROF_HOP_CANDIDATE_REJECT     21
#define BFPROF_HOP_VALID                22
#define BFPROF_HOP_ACCEPT               23
#define BFPROF_HOP_METROPOLIS_REJECT    24
#define BFPROF_EXCHANGE_TRY             25
#define BFPROF_EXCHANGE_CANDIDATE_REJECT 26
#define BFPROF_EXCHANGE_VALID           27
#define BFPROF_EXCHANGE_ACCEPT          28
#define BFPROF_EXCHANGE_METROPOLIS_REJECT 29
#define BFPROF_SAMPLE_TERM_DENSE_CANDIDATE 30
#define BFPROF_SAMPLE_TERM_GEOMETRY_VALID  31
#define BFPROF_SAMPLE_TERM_SPARSE_PAIR     32
#define BFPROF_SAMPLE_TERM_ACTUAL_ADD      33
#define BFPROF_GREEN_TERM_DENSE_CANDIDATE  34
#define BFPROF_GREEN_TERM_GEOMETRY_VALID   35
#define BFPROF_GREEN_TERM_SPARSE_PAIR      36
#define BFPROF_GREEN_TERM_ACTUAL_ADD       37
#define BFPROF_SAMPLE_DELTA_CNT_TOTAL      38
#define BFPROF_SAMPLE_DELTA_CNT_CHANGED    39
#define BFPROF_SAMPLE_DELTA_PAIR_NEW       40
#define BFPROF_SAMPLE_DELTA_PAIR_OLD       41
#define BFPROF_SAMPLE_DELTA_PAIR_TOTAL     42
#define BFPROF_GREEN_DELTA_CNT_TOTAL       43
#define BFPROF_GREEN_DELTA_CNT_CHANGED     44
#define BFPROF_GREEN_DELTA_PAIR_NEW        45
#define BFPROF_GREEN_DELTA_PAIR_OLD        46
#define BFPROF_GREEN_DELTA_PAIR_TOTAL      47
int BFProfileEnabled = 0;
long long BFProfileCounter[NBFProfileCounter];

#define BFFSZ_PROFILE_SAMPLE 0
#define BFFSZ_PROFILE_GREEN  1
#define NBFFSZProfileSource  2
#define NBFFSZProfileHist    22
#define NBFFSZProfileRatioHist 8
#define BFFSZ_PF_PATH_OPTIMIZED 0
#define BFFSZ_PF_PATH_DIRECT_FULL 1
#define BFFSZ_PF_PATH_FALLBACK 2
#define NBFFSZPfPath 3
#define BF_FSZ_PF_UPDATE_KFULL_DEFAULT 32
int BFFSZGreenRebuildCheckEnabled = 0;
int BFFSZAffectedCheckEnabled = 0;
int BFFSZPfUpdateCheckEnabled = 0;
int BFFSZPfUpdateForceFallback = 0;
int BFFSZPfUpdateInjectedStatus = 0;
int BFFSZPfUpdateInjectedRank = -1;
int BFFSZPfUpdateExplicitStateCheckEnabled = 0;
int BFFSZPfUpdateArgumentCheckEnabled = 0;
int BFFSZPfUpdateKFull = BF_FSZ_PF_UPDATE_KFULL_DEFAULT;
int BFFSZPermuteParticleLabelsCheckEnabled = 0;
int BFFSZInvUpdateCheckEnabled = 0;
int BFFSZInvUpdateForceFallback = 0;
int BFFSZInvUpdateInjectedStage = 0;
int BFFSZInvUpdateInjectedRank = -1;
int BFFSZInvUpdateExplicitStateCheckEnabled = 0;
int BFFSZInvUpdateArgumentCheckEnabled = 0;
long long BFFSZProfileCall[NBFFSZProfileSource];
long long BFFSZProfileChangedSum[NBFFSZProfileSource];
long long BFFSZProfileChangedMax[NBFFSZProfileSource];
long long BFFSZProfileAffectedSum[NBFFSZProfileSource];
long long BFFSZProfileAffectedMax[NBFFSZProfileSource];
long long BFFSZProfileChangedHist[NBFFSZProfileSource][NBFFSZProfileHist];
long long BFFSZProfileAffectedHist[NBFFSZProfileSource][NBFFSZProfileHist];
long long BFFSZProfileAffectedRatioHist[NBFFSZProfileSource][NBFFSZProfileRatioHist];
long long BFFSZProfilePfPath[NBFFSZProfileSource][NBFFSZPfPath];
long long BFFSZProfileInvPath[NBFFSZPfPath];
double BFFSZProfileInvCheckSeconds = 0.0;
double BFFSZProfileInvAntisymmetryMax = 0.0;
double BFFSZProfileInvResidualMax = 0.0;

static inline void AddBFProfileCounter(int idx, long long value) {
  if(!BFProfileEnabled || value == 0) return;
#pragma omp atomic
  BFProfileCounter[idx] += value;
}

static inline int BFFSZProfileHistBin(int value) {
  if(value <= 0) return 0;
  if(value <= 16) return value;
  if(value <= 32) return 17;
  if(value <= 64) return 18;
  if(value <= 128) return 19;
  if(value <= 256) return 20;
  return 21;
}

static inline int BFFSZProfileRatioHistBin(int value, int total) {
  if(value <= 0 || total <= 0) return 0;
  if(16LL*value <= total) return 1;
  if(8LL*value <= total) return 2;
  if(4LL*value <= total) return 3;
  if(2LL*value <= total) return 4;
  if(4LL*value <= 3LL*total) return 5;
  if(value < total) return 6;
  return 7;
}

static inline void RecordBFFSZProfile(int source, int nChanged, int nAffected) {
  int changedBin, affectedBin, ratioBin;
  if(!BFProfileEnabled) return;
  if(source < 0 || source >= NBFFSZProfileSource) return;
  changedBin = BFFSZProfileHistBin(nChanged);
  affectedBin = BFFSZProfileHistBin(nAffected);
  ratioBin = BFFSZProfileRatioHistBin(nAffected, Nsize);
#pragma omp atomic
  BFFSZProfileCall[source]++;
#pragma omp atomic
  BFFSZProfileChangedSum[source] += nChanged;
#pragma omp atomic
  BFFSZProfileAffectedSum[source] += nAffected;
#pragma omp atomic
  BFFSZProfileChangedHist[source][changedBin]++;
#pragma omp atomic
  BFFSZProfileAffectedHist[source][affectedBin]++;
#pragma omp atomic
  BFFSZProfileAffectedRatioHist[source][ratioBin]++;
#pragma omp critical(BFFSZProfileMax)
  {
    if(nChanged > BFFSZProfileChangedMax[source]) {
      BFFSZProfileChangedMax[source] = nChanged;
    }
    if(nAffected > BFFSZProfileAffectedMax[source]) {
      BFFSZProfileAffectedMax[source] = nAffected;
    }
  }
}

static inline void RecordBFFSZPfPath(int source, int path) {
  if(!BFProfileEnabled) return;
  if(source < 0 || source >= NBFFSZProfileSource
      || path < 0 || path >= NBFFSZPfPath) return;
#pragma omp atomic
  BFFSZProfilePfPath[source][path]++;
}

static inline void RecordBFFSZInvPath(int path) {
  if(!BFProfileEnabled || path < 0 || path >= NBFFSZPfPath) return;
  BFFSZProfileInvPath[path]++;
}

static inline void RecordBFFSZInvChecks(double seconds,
                                        double antisymmetry,
                                        double residual) {
  if(!BFProfileEnabled) return;
  BFFSZProfileInvCheckSeconds += seconds;
  if(antisymmetry > BFFSZProfileInvAntisymmetryMax) {
    BFFSZProfileInvAntisymmetryMax = antisymmetry;
  }
  if(residual > BFFSZProfileInvResidualMax) {
    BFFSZProfileInvResidualMax = residual;
  }
}

int useDiagScale=0;
int RescaleSmat=0;

#endif /*  _INCLUDE_GLOBAL */
