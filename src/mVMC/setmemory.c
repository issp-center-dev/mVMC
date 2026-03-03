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
 * Allocate and free memory for global array
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/


#include <complex.h>
#include "global.h"
#include "setmemory.h"

#ifndef _SRC_SETMEMORY
#define _SRC_SETMEMORY

void SetMemoryDef() {  
  int i, j;
  int *pInt;
  double *pDouble;

  /* Int */
  LocSpn = (int*)malloc(sizeof(int)*NTotalDefInt);
  pInt = LocSpn + Nsite;

  Transfer = (int**)malloc(sizeof(int*)*NTransfer);
  for(i=0;i<NTransfer;i++) {
    Transfer[i] = pInt;
    pInt += 4;
  }

  CoulombIntra = pInt;
  pInt += NCoulombIntra;

  CoulombInter = (int**)malloc(sizeof(int*)*NCoulombInter);
  for(i=0;i<NCoulombInter;i++) {
    CoulombInter[i] = pInt;
    pInt += 2;
  }

  HundCoupling = (int**)malloc(sizeof(int*)*NHundCoupling);
  for(i=0;i<NHundCoupling;i++) {
    HundCoupling[i] = pInt;
    pInt += 2;
  }

  PairHopping = (int**)malloc(sizeof(int*)*NPairHopping);
  for(i=0;i<NPairHopping;i++) {
    PairHopping[i] = pInt;
    pInt += 2;
  }

  ExchangeCoupling = (int**)malloc(sizeof(int*)*NExchangeCoupling);
  for(i=0;i<NExchangeCoupling;i++) {
    ExchangeCoupling[i] = pInt;
    pInt += 2;
  }

  GutzwillerIdx = pInt;
  pInt += Nsite;

  JastrowIdx = (int**)malloc(sizeof(int*)*Nsite);
  for(i=0;i<Nsite;i++) {
    JastrowIdx[i] = pInt;
    pInt += Nsite;
  }

  DoublonHolon2siteIdx = (int**)malloc(sizeof(int*)*NDoublonHolon2siteIdx);
  for(i=0;i<NDoublonHolon2siteIdx;i++) {
    DoublonHolon2siteIdx[i] = pInt;
    pInt += 2*Nsite;
  }

  DoublonHolon4siteIdx = (int**)malloc(sizeof(int*)*NDoublonHolon4siteIdx);
  for(i=0;i<NDoublonHolon4siteIdx;i++) {
    DoublonHolon4siteIdx[i] = pInt;
    pInt += 4*Nsite;
  }

//RBM
  if (FlagRBM) {
    ChargeRBM_PhysLayerIdx = pInt;
    pInt += Nsite;
    SpinRBM_PhysLayerIdx = pInt;
    pInt += Nsite;
    GeneralRBM_PhysLayerIdx = pInt;
    pInt += Nsite2;
  
    ChargeRBM_HiddenLayerIdx = pInt;
    pInt += NneuronCharge;
    SpinRBM_HiddenLayerIdx = pInt;
    pInt += NneuronSpin;
    GeneralRBM_HiddenLayerIdx = pInt;
    pInt += NneuronGeneral;
  
    ChargeRBM_PhysHiddenIdx = (int**)malloc(sizeof(int*)*Nsite);
    for(i=0;i<Nsite;i++) {
      ChargeRBM_PhysHiddenIdx[i] = pInt;
      pInt += NneuronCharge;
    }
    SpinRBM_PhysHiddenIdx = (int**)malloc(sizeof(int*)*Nsite);
    for(i=0;i<Nsite;i++) {
      SpinRBM_PhysHiddenIdx[i] = pInt;
      pInt += NneuronSpin;
    }
    GeneralRBM_PhysHiddenIdx = (int**)malloc(sizeof(int*)*Nsite2);
    for(i=0;i<Nsite2;i++) {
      GeneralRBM_PhysHiddenIdx[i] = pInt;
      pInt += NneuronGeneral;
    }
  }
//RBM

 /*[s] For BackFlow */
  if(NBackFlowIdx>0) {
    PosBF = (int**)malloc(sizeof(int*)*Nsite);
    for(i=0;i<Nsite;i++) {
      PosBF[i] = pInt;
      pInt += Nrange;
    }
    RangeIdx = (int**)malloc(sizeof(int*)*Nsite);
    for(i=0;i<Nsite;i++) {
      RangeIdx[i] = pInt;
      pInt += Nsite;
    }
    BackFlowIdx = (int**)malloc(sizeof(int*)*Nsite*Nsite);
    for(i=0;i<Nsite*Nsite;i++) {
      BackFlowIdx[i] = pInt;
      pInt += Nsite*Nsite;
    }
  }
  /*[e] For BackFlow */

  int NOrbit;
  iFlgOrbitalGeneral==0 ? (NOrbit=Nsite): (NOrbit=2*Nsite);
  OrbitalIdx = (int**)malloc(sizeof(int*)*NOrbit);
  for(i=0;i<NOrbit;i++) {
    OrbitalIdx[i] = pInt;
    pInt += NOrbit;
    for(j=0;j<NOrbit;j++) {
      OrbitalIdx[i][j]=0;
    }
  }
  OrbitalSgn = (int**)malloc(sizeof(int*)*NOrbit);
  for(i=0;i<NOrbit;i++) {
    OrbitalSgn[i] = pInt;
    pInt += NOrbit;
    for(j=0;j<NOrbit;j++) {
      OrbitalSgn[i][j]=0;
    }
  }

  QPTrans = (int**)malloc(sizeof(int*)*NQPTrans);
  for(i=0;i<NQPTrans;i++) {
    QPTrans[i] = pInt;
    pInt += Nsite;
  }

  QPTransInv = (int**)malloc(sizeof(int*)*NQPTrans);
  for(i=0;i<NQPTrans;i++) {
    QPTransInv[i] = pInt;
    pInt += Nsite;
  }

  QPTransSgn = (int**)malloc(sizeof(int*)*NQPTrans);
  for(i=0;i<NQPTrans;i++) {
    QPTransSgn[i] = pInt;
    pInt += Nsite;
  }

  CisAjsIdx = (int**)malloc(sizeof(int*)*NCisAjs);
  for(i=0;i<NCisAjs;i++) {
    CisAjsIdx[i] = pInt;
    pInt += 4;
  }

  CisAjsCktAltIdx = (int**)malloc(sizeof(int*)*NCisAjsCktAlt);
  for(i=0;i<NCisAjsCktAlt;i++) {
    CisAjsCktAltIdx[i] = pInt;
    pInt += 2;
  }

  CisAjsCktAltDCIdx = (int**)malloc(sizeof(int*)*NCisAjsCktAltDC);
  for(i=0;i<NCisAjsCktAltDC;i++) {
    CisAjsCktAltDCIdx[i] = pInt;
    pInt += 8;
  }

  LatticeIdx = (int**)malloc(sizeof(int*)*Nsite);
  for(i=0;i<Nsite;i++) {
    LatticeIdx[i] = pInt;
    pInt += 4;
  }

  TwistIdx = (int**)malloc(sizeof(int*)*NTwist);
  for(i=0;i<NTwist;i++) {
    TwistIdx[i] = pInt;
    pInt += 2*Nsite*2;
  }
  
  InterAll = (int**)malloc(sizeof(int*)*NInterAll);
  for(i=0;i<NInterAll;i++) {
    InterAll[i] = pInt;
    pInt += 8;
  }

  QPOptTrans = (int**)malloc(sizeof(int*)*NQPOptTrans);
  for(i=0;i<NQPOptTrans;i++) {
    QPOptTrans[i] = pInt;
    pInt += Nsite;
  }

  QPOptTransSgn = (int**)malloc(sizeof(int*)*NQPOptTrans);
  for(i=0;i<NQPOptTrans;i++) {
    QPOptTransSgn[i] = pInt;
    pInt += Nsite;
  }

  OptFlag = pInt;

  ParaTransfer = (double complex*)malloc(sizeof(double complex)*(NTransfer+NInterAll));  
  ParaInterAll = ParaTransfer+NTransfer;

  ParaCoulombIntra = (double*)malloc(sizeof(double)*(NTotalDefDouble));
  pDouble = ParaCoulombIntra +NCoulombIntra; 

  ParaCoulombInter = pDouble;
  pDouble += NCoulombInter;

  ParaHundCoupling = pDouble;
  pDouble += NHundCoupling;

  ParaPairHopping = pDouble;
  pDouble +=  NPairHopping;

  ParaExchangeCoupling = pDouble;
  pDouble +=  NExchangeCoupling;
  
//  ParaQPTrans = pDouble;
//  pDouble +=  NQPTrans;

  ParaTwist = (double**)malloc(sizeof(double*)*NTwist);
  for(i=0;i<NTwist;i++) {
    ParaTwist[i] = pDouble;
    pDouble += 3*Nsite*2;
  }

  ParaQPOptTrans = pDouble;
  ParaQPTrans = (double complex*)malloc(sizeof(double complex)*(NQPTrans));

  return;
}

void FreeMemoryDef() {
  int i, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==0 && NCisAjsCktAlt>0) {
    for(i=0;i<2*Nsite;i++)
      free(iOneBodyGIdx[i]);
    free(iOneBodyGIdx);
  } 

  free(QPOptTransSgn);
  free(QPOptTrans);
  free(InterAll);
  free(CisAjsCktAltDCIdx);
  free(CisAjsCktAltIdx);
  free(CisAjsIdx);
  free(QPTransSgn);
  free(QPTrans);
  free(QPTransInv);
  free(OrbitalIdx);
  free(OrbitalSgn);
  free(DoublonHolon4siteIdx);
  free(DoublonHolon2siteIdx);
  free(JastrowIdx);
  free(ParaTransfer);
  free(ParaCoulombIntra);
  free(ParaQPTrans);
  free(ExchangeCoupling);
  free(PairHopping);
  free(HundCoupling);
  free(CoulombInter);
  free(Transfer);
  free(LocSpn);
  free(PosBF);
  free(RangeIdx);
  free(BackFlowIdx);
  return;
}

void SetMemory() {
  int i;

  /***** Variational Parameters *****/
  //printf("DEBUG:opt=%d %d %d %d %d Ne=%d\n", AllComplexFlag,NPara,NProj,NSlater,NOrbitalIdx,Ne);
  Para     = (double complex*)malloc(sizeof(double complex)*(NPara));

  Proj     = Para;
  RBM      = Para + NProj;
  ProjBF   = Para + NProj + FlagRBM*NRBM;
  Slater   = Para + NProj + FlagRBM*NRBM + NProjBF;
  OptTrans = Para + NProj + FlagRBM*NRBM + NProjBF + NSlater;

  /***** Electron Configuration ******/
  EleIdx            = (int*)malloc(sizeof(int)*( NVMCSample*2*Ne ));
  EleCfg            = (int*)malloc(sizeof(int)*( NVMCSample*2*Nsite ));
  EleNum            = (int*)malloc(sizeof(int)*( NVMCSample*2*Nsite ));
  EleProjCnt        = (int*)malloc(sizeof(int)*( NVMCSample*NProj ));
//[s] MERGE BY TM
  EleSpn            = (int*)malloc(sizeof(int)*( NVMCSample*2*Ne ));//fsz
//[e] MERGE BY TM
    logSqPfFullSlater = (double*)malloc(sizeof(double)*(NVMCSample));
  //EleProjBFCnt = (int*)malloc(sizeof(int)*( NVMCSample*4*4*Nsite*Nrange));
  if (NBackFlowIdx > 0) {
    EleProjBFCnt = (int*)malloc(sizeof(int)*( NVMCSample*4*4*Nsite*Nrange));
    SmpSltElmBF_real = (double *)malloc(sizeof(double)*(NVMCSample*NQPFull*(2*Nsite)*(2*Nsite)));
    SmpEta = (double*)malloc(sizeof(double*)*NVMCSample*NQPFull*Nsite*Nsite);
    SmpEtaFlag = (int*)malloc(sizeof(int*)*NVMCSample*NQPFull*Nsite*Nsite);
        SlaterElmBF_real = (double*)malloc( sizeof(double)*(NQPFull*(2*Nsite)*(2*Nsite)) );
    eta = (double complex**)malloc(sizeof(double complex*)*Nsite);
      for(i=0;i<Nsite;i++) {
          eta[i] = (double complex*)malloc(sizeof(double complex)*Nsite);
      }
      etaFlag = (int**)malloc(sizeof(int*)*Nsite);
      for(i=0;i<Nsite;i++) {
          etaFlag[i] = (int*)malloc(sizeof(int)*Nsite);
      }
      BFSubIdx = (int**)malloc(sizeof(int*)*NrangeIdx);
      for(i=0;i<NrangeIdx;i++) {
          BFSubIdx[i] = (int*)malloc(sizeof(int)*NrangeIdx);
      }
  }

  TmpEleIdx         = (int*)malloc(sizeof(int)*(2*Ne+2*Nsite+2*Nsite+NProj+2*Ne));//fsz
  TmpEleCfg         = TmpEleIdx + 2*Ne;
  TmpEleNum         = TmpEleCfg + 2*Nsite;
  TmpEleProjCnt     = TmpEleNum + 2*Nsite;
//[s] MERGE BY TM
  TmpEleSpn         = TmpEleProjCnt + NProj; //fsz
  TmpEleProjBFCnt = TmpEleProjCnt + NProj;
//[e] MERGE BY TM

  BurnEleIdx        = (int*)malloc(sizeof(int)*(2*Ne+2*Nsite+2*Nsite+NProj+2*Ne)); //fsz
  BurnEleCfg        = BurnEleIdx + 2*Ne;
  BurnEleNum        = BurnEleCfg + 2*Nsite;
  BurnEleProjCnt    = BurnEleNum + 2*Nsite;
  BurnEleSpn        = BurnEleProjCnt + NProj; //fsz

//RBM
  if (FlagRBM) {
    RBMCnt = (double complex*)malloc( sizeof(double complex)*( NVMCSample*(NRBM_PhysLayerIdx+Nneuron) ) );
    TmpRBMCnt = (double complex*)malloc( sizeof(double complex)*(NRBM_PhysLayerIdx+Nneuron) );
    BurnRBMCnt = (double complex*)malloc( sizeof(double complex)*(NRBM_PhysLayerIdx+Nneuron) );
  }
//RBM

  /***** Slater Elements ******/
  SlaterElm = (double complex*)malloc( sizeof(double complex)*(NQPFull*(2*Nsite)*(2*Nsite)) );
  InvM = (double complex*)malloc( sizeof(double complex)*(NQPFull*(Nsize*Nsize+1)) );
  PfM = InvM + NQPFull*Nsize*Nsize;
// for real TBC
  if (AllComplexFlag == 0){
    SlaterElm_real = (double*)malloc(sizeof(double)*(NQPFull*(2*Nsite)*(2*Nsite)) );
    InvM_real      = (double*)malloc(sizeof(double)*(NQPFull*(Nsize*Nsize+1)) );
    PfM_real       = InvM_real + NQPFull*Nsize*Nsize;
  }

  if (FlagRBM) {
    SlaterElmBF_real = (double*)malloc( sizeof(double)*(NQPFull*(2*Nsite)*(2*Nsite)) );
    eta = (double complex**)malloc(sizeof(double complex*)*Nsite);
    for(i=0;i<Nsite;i++) {
      eta[i] = (double complex*)malloc(sizeof(double complex)*Nsite);
    }
    etaFlag = (int**)malloc(sizeof(int*)*Nsite);
    for(i=0;i<Nsite;i++) {
      etaFlag[i] = (int*)malloc(sizeof(int)*Nsite);
    }
    BFSubIdx = (int**)malloc(sizeof(int*)*NrangeIdx);
    for(i=0;i<NrangeIdx;i++) {
      BFSubIdx[i] = (int*)malloc(sizeof(int)*NrangeIdx);
    }
  }

  /***** Quantum Projection *****/
  QPFullWeight = (double complex*)malloc(sizeof(double complex)*(NQPFull+NQPFix+5*NSPGaussLeg));
  QPFixWeight= QPFullWeight + NQPFull;
  SPGLCos    = QPFullWeight + NQPFull + NQPFix;
  SPGLSin    = SPGLCos + NSPGaussLeg;
  SPGLCosSin = SPGLCos + 2*NSPGaussLeg;
  SPGLCosCos = SPGLCos + 3*NSPGaussLeg;
  SPGLSinSin = SPGLCos + 4*NSPGaussLeg;

  /***** Stocastic Reconfiguration *****/
  if(NVMCCalMode==0){
    //SR components are described by real and complex components of O
    if(NSRCG==0){
      SROptOO = (double complex*)malloc( sizeof(double complex)*((2*SROptSize)*(2*SROptSize+2))) ; //TBC
      SROptHO = SROptOO + (2*SROptSize)*(2*SROptSize); //TBC
      SROptO  = SROptHO + (2*SROptSize);  //TBC
    }else{
      // OO contains only <O_i> and <O_i O_i> in SR-CG
      SROptOO = (double complex*)malloc( sizeof(double complex)*(2*SROptSize)*4) ; //TBC
      SROptHO = SROptOO + 2*SROptSize*2; //TBC
      SROptO  = SROptHO + 2*SROptSize;  //TBC
    }
//for real
    if (AllComplexFlag == 0){
      if(NSRCG==0){
        SROptOO_real = (double*)malloc( sizeof(double )*SROptSize*(SROptSize+2)) ; //TBC
        SROptHO_real = SROptOO_real + (SROptSize)*(SROptSize); //TBC
        SROptO_real  = SROptHO_real + (SROptSize);  //TBC
      }else{
        // OO contains only <O_i> and <O_i O_i> in SR-CG
        SROptOO_real = (double*)malloc( sizeof(double )*SROptSize*4) ; //TBC
        SROptHO_real = SROptOO_real + SROptSize*2; //TBC
        SROptO_real  = SROptHO_real + SROptSize;  //TBC
      }
    }

    if(NSRCG==1 || NStoreO!=0){
      //if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){ //real & sz=0
      if(AllComplexFlag==0){ //real & sz=0
        SROptO_Store_real = (double *)malloc(sizeof(double)*(SROptSize*NVMCSample) );
      }else{
        SROptO_Store      = (double complex*)malloc( sizeof(double complex)*(2*SROptSize*NVMCSample) );
      }
    }
    SROptData = (double complex*)malloc( sizeof(double complex)*(NSROptItrSmp*(2+NPara)) );
  }

  /***** Physical Quantity *****/
  if(NVMCCalMode==1){
    PhysCisAjs  = (double complex*)malloc(sizeof(double complex)
                    *(NCisAjs+NCisAjsCktAlt+NCisAjsCktAltDC+NCisAjs+NCisAjsCktAltDC + NTwist));
    PhysCisAjsCktAlt   = PhysCisAjs       + NCisAjs;
    PhysCisAjsCktAltDC = PhysCisAjsCktAlt + NCisAjsCktAlt;
    PhysTwist = PhysCisAjsCktAltDC + NCisAjsCktAltDC;
    LocalCisAjs = PhysTwist + NTwist;
    //LocalCisAjs = PhysCisAjsCktAltDC + NCisAjsCktAltDC;
    LocalCisAjsCktAltDC = LocalCisAjs + NCisAjs;

    if(NLanczosMode>0){
      QQQQ = (double complex*)malloc(sizeof(double complex)
        *(NLSHam*NLSHam*NLSHam*NLSHam + NLSHam*NLSHam) );
      LSLQ = QQQQ + NLSHam*NLSHam*NLSHam*NLSHam;
      //for real
      QQQQ_real = (double*)malloc(sizeof(double)
      *(NLSHam*NLSHam*NLSHam*NLSHam + NLSHam*NLSHam) );
      LSLQ_real = QQQQ_real + NLSHam*NLSHam*NLSHam*NLSHam;

      if(NLanczosMode>1){
        QCisAjsQ = (double complex*)malloc(sizeof(double complex)
          *(NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*(NCisAjsCktAltDC+NCisAjsCktAlt) + NLSHam*NCisAjs) );
        QCisAjsCktAltQ = QCisAjsQ + NLSHam*NLSHam*NCisAjs;
        QCisAjsCktAltQDC = QCisAjsCktAltQ + NLSHam*NLSHam*NCisAjsCktAlt;
        LSLCisAjs = QCisAjsCktAltQDC + NLSHam*NLSHam*NCisAjsCktAltDC;
        //for real
        QCisAjsQ_real = (double *)malloc(sizeof(double )
        *(NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*(NCisAjsCktAltDC+NCisAjsCktAlt) + NLSHam*NCisAjs) );
        QCisAjsCktAltQ_real = QCisAjsQ_real + NLSHam*NLSHam*NCisAjs;
        QCisAjsCktAltQDC_real = QCisAjsCktAltQ_real + NLSHam*NLSHam*NCisAjsCktAlt;
        LSLCisAjs_real = QCisAjsCktAltQDC_real + NLSHam*NLSHam*NCisAjsCktAltDC;

      }
    }
  }

  initializeWorkSpaceAll();
  return;
}

void FreeMemory() {
  FreeWorkSpaceAll();

  if(NVMCCalMode==1){
    free(PhysCisAjs);
    if(NLanczosMode>0){
      free(QQQQ);
      free(QQQQ_real);
      if(NLanczosMode>1){
        free(QCisAjsQ);
        free(QCisAjsQ_real);
      }
    }
  }

  if(NVMCCalMode==0){
    free(SROptData);
    free(SROptOO);
    //for real
    if (AllComplexFlag == 0){
        free(SROptOO_real);
    }
    if(NSRCG==1 || NStoreO!=0){
      if(AllComplexFlag==0){ //real & sz=0
        free(SROptO_Store_real);
      }else{
        free(SROptO_Store);
      }
    }
  }

  free(QPFullWeight);

  free(InvM);
  free(SlaterElm);

  if (AllComplexFlag == 0){
    free(InvM_real);
    free(SlaterElm_real);
  }

  free(BurnEleIdx);
  free(TmpEleIdx);
  free(logSqPfFullSlater);
  free(EleProjCnt);
  free(EleIdx);
  free(EleNum);
  free(EleSpn);
  free(EleCfg);

  free(Para);

  return;
}

#endif
