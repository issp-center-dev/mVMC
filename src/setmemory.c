/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

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
    pInt += 8;
  }

  CisAjsCktAltDCIdx = (int**)malloc(sizeof(int*)*NCisAjsCktAltDC);
  for(i=0;i<NCisAjsCktAltDC;i++) {
    CisAjsCktAltDCIdx[i] = pInt;
    pInt += 8;
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

  
  ParaQPOptTrans = pDouble;
  ParaQPTrans = (double complex*)malloc(sizeof(double complex)*(NQPTrans));

  // LanczosGreen
  if(NLanczosMode>1){
    CisAjsCktAltLzIdx = malloc(sizeof(int*)*NCisAjsCktAltDC);
    for(i=0;i<NCisAjsCktAltDC;i++) {
      CisAjsCktAltLzIdx[i] = malloc(sizeof(int) * 2);
    }
  }

  return;
}

void FreeMemoryDef() {
  free(ParaTransfer);

  free(QPOptTransSgn);
  free(QPOptTrans);
  free(InterAll);
  free(CisAjsCktAltDCIdx);
  free(CisAjsCktAltIdx);
  free(CisAjsIdx);
  free(QPTransSgn);
  free(QPTrans);
  free(OrbitalIdx);
  free(DoublonHolon4siteIdx);
  free(DoublonHolon2siteIdx);
  free(JastrowIdx);
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
  int i,j,n;

  /***** Variational Parameters *****/
  //printf("DEBUG:opt=%d %d %d %d %d Ne=%d\n", AllComplexFlag,NPara,NProj,NSlater,NOrbitalIdx,Ne);
  Para     = (double complex*)malloc(sizeof(double complex)*(NPara));

  Proj     = Para;
  ProjBF   = Para + NProj;
  Slater   = Para + NProj + NProjBF;
  OptTrans = Para + NProj + NProjBF + NSlater;

  /***** Electron Configuration ******/
  EleIdx            = (int*)malloc(sizeof(int)*( NVMCSample*2*Ne ));
  EleCfg            = (int*)malloc(sizeof(int)*( NVMCSample*2*Nsite ));
  EleNum            = (int*)malloc(sizeof(int)*( NVMCSample*2*Nsite ));
  EleProjCnt        = (int*)malloc(sizeof(int)*( NVMCSample*NProj ));
//[s] MERGE BY TM
  EleSpn            = (int*)malloc(sizeof(int)*( NVMCSample*2*Ne ));//fsz
  EleProjBFCnt = (int*)malloc(sizeof(int)*( NVMCSample*4*4*Nsite*Nrange));
//[e] MERGE BY TM
  logSqPfFullSlater = (double*)malloc(sizeof(double)*(NVMCSample));
  SmpSltElmBF_real = (double *)malloc(sizeof(double)*(NVMCSample*NQPFull*(2*Nsite)*(2*Nsite)));
  SmpEta = (double*)malloc(sizeof(double*)*NVMCSample*NQPFull*Nsite*Nsite);
  SmpEtaFlag = (int*)malloc(sizeof(int*)*NVMCSample*NQPFull*Nsite*Nsite);

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

  /***** Slater Elements ******/
  SlaterElm = (double complex*)malloc( sizeof(double complex)*(NQPFull*(2*Nsite)*(2*Nsite)) );
  InvM = (double complex*)malloc( sizeof(double complex)*(NQPFull*(Nsize*Nsize+1)) );
  PfM = InvM + NQPFull*Nsize*Nsize;
// for real TBC
  SlaterElm_real = (double*)malloc(sizeof(double)*(NQPFull*(2*Nsite)*(2*Nsite)) );
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
  InvM_real      = (double*)malloc(sizeof(double)*(NQPFull*(Nsize*Nsize+1)) );
  PfM_real       = InvM_real + NQPFull*Nsize*Nsize;

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
    SROptOO = (double complex*)malloc( sizeof(double complex)*((2*SROptSize)*(2*SROptSize+2))) ; //TBC
    SROptHO = SROptOO + (2*SROptSize)*(2*SROptSize); //TBC
    SROptO  = SROptHO + (2*SROptSize);  //TBC
//for real
    SROptOO_real = (double*)malloc( sizeof(double )*((SROptSize)*(SROptSize+2))) ; //TBC
    SROptHO_real = SROptOO_real + (SROptSize)*(SROptSize); //TBC
    SROptO_real  = SROptHO_real + (SROptSize);  //TBC

    if(NStoreO!=0){
      if(AllComplexFlag==0){
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
                    *(2*NCisAjs+NCisAjsCktAlt+NCisAjsCktAltDC));
    PhysCisAjsCktAlt   = PhysCisAjs       + NCisAjs;
    PhysCisAjsCktAltDC = PhysCisAjsCktAlt + NCisAjsCktAlt;
    LocalCisAjs = PhysCisAjsCktAltDC + NCisAjsCktAltDC;

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
          *(NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*NCisAjsCktAltDC + NLSHam*NCisAjs) );
        QCisAjsCktAltQ = QCisAjsQ + NLSHam*NLSHam*NCisAjs;
        LSLCisAjs = QCisAjsCktAltQ + NLSHam*NLSHam*NCisAjsCktAltDC;
        //for real
        QCisAjsQ_real = (double *)malloc(sizeof(double )
        *(NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*NCisAjsCktAltDC + NLSHam*NCisAjs) );
        QCisAjsCktAltQ_real = QCisAjsQ_real + NLSHam*NLSHam*NCisAjs;
        LSLCisAjs_real = QCisAjsCktAltQ_real + NLSHam*NLSHam*NCisAjsCktAltDC;

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
  }

  free(QPFullWeight);

  free(InvM);
  free(SlaterElm);

  free(BurnEleIdx);
  free(TmpEleIdx);
  free(logSqPfFullSlater);
  free(EleProjCnt);
  free(EleIdx);
  free(EleCfg);

  free(Para);

  return;
}

#endif
