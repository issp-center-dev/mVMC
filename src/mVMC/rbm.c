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
 * RBM factors
 *-------------------------------------------------------------*/
#include "./include/global.h"
#include "math.h"
#include "./include/rbm.h"

inline double complex WeightRBM(const double complex *rbmCnt) {
  int idx,hi;
  double complex z=0.0;

  #pragma omp parallel for default(shared) private(hi)  firstprivate(NRBM_PhysLayerIdx,Nneuron) reduction(+:z)
  for(idx=0;idx<NRBM_PhysLayerIdx;idx++) {
    z += RBM[idx]*rbmCnt[idx];
  }

  #pragma omp parallel for default(shared) private(hi)  firstprivate(NRBM_PhysLayerIdx,Nneuron) reduction(+:z)
  for(hi=0;hi<Nneuron;hi++) {
    z += clog( ccosh( rbmCnt[hi+NRBM_PhysLayerIdx] ) );
  }
  
  return cexp(z);
}

inline double complex LogWeightRBM(const double complex *rbmCnt) {
  int idx,hi;
  double complex z=0.0;

  #pragma omp parallel for default(shared) private(hi)  firstprivate(NRBM_PhysLayerIdx,Nneuron) reduction(+:z)
  for(idx=0;idx<NRBM_PhysLayerIdx;idx++) {
    z += RBM[idx]*rbmCnt[idx];
  }

  #pragma omp parallel for default(shared) private(hi)  firstprivate(NRBM_PhysLayerIdx,Nneuron) reduction(+:z)
  for(hi=0;hi<Nneuron;hi++) {
    z += clog( ccosh( rbmCnt[hi+NRBM_PhysLayerIdx] ) );
  }

  return z;
}

inline double complex RBMRatio(const double complex *rbmCntNew, const double complex *rbmCntOld) {
  int icnt;
  int idx,hi;
  int Nblk,iblk;
  int hist,hiend;
  double complex z=0.0, zz=1.0;
  double complex RbmCntNew[NBlockSize_RBMRatio], RbmCntOld[NBlockSize_RBMRatio], zzTmp[NBlockSize_RBMRatio];

  #pragma omp parallel for default(shared) private(hi) firstprivate(NRBM_PhysLayerIdx,Nneuron) reduction(+:z)
  for(idx=0;idx<NRBM_PhysLayerIdx;idx++) {
    z += RBM[idx]*(rbmCntNew[idx]-rbmCntOld[idx]);
  }

  Nblk = (Nneuron-1)/NBlockSize_RBMRatio + 1;
  for(iblk=0;iblk<Nblk;iblk++) { 
    hist = iblk*NBlockSize_RBMRatio;
    hiend= hist+NBlockSize_RBMRatio;
    if(hiend > Nneuron) hiend = Nneuron; 
    zz=1.0;
#pragma omp parallel default(shared) private(hi,RbmCntNew,RbmCntOld,zzTmp,icnt) firstprivate(NRBM_PhysLayerIdx,Nneuron) reduction(+:z) reduction(*:zz)
    {
      icnt=0;
#pragma loop swp
#pragma omp for
      for(hi=hist;hi<hiend;hi++) {
        RbmCntNew[icnt] = rbmCntNew[hi+NRBM_PhysLayerIdx];
        RbmCntOld[icnt] = rbmCntOld[hi+NRBM_PhysLayerIdx];
        if(creal(RbmCntNew[icnt]) <= 0.0 ) RbmCntNew[icnt] = -RbmCntNew[icnt];
        if(creal(RbmCntOld[icnt]) <= 0.0 ) RbmCntOld[icnt] = -RbmCntOld[icnt];
        z+=( RbmCntNew[icnt] - RbmCntOld[icnt] );
        zzTmp[icnt] = (1.0+cexp(-2.0*RbmCntNew[icnt])) / (1.0+cexp(-2.0*RbmCntOld[icnt]));
        icnt++;
      }

// Changing for SIMD tuning.
// If the loop length is multiple of 8, It is enabled SIMD.
// If the loop length is not multiple of 8, It is not SIMD.
      if((icnt%8)==0) {
        icnt=0;
#pragma omp for
        for(hi=hist;hi<hiend;hi=hi+8) {
          zz*=zzTmp[icnt]*zzTmp[icnt+1]*zzTmp[icnt+2]*zzTmp[icnt+3]*
              zzTmp[icnt+4]*zzTmp[icnt+5]*zzTmp[icnt+6]*zzTmp[icnt+7];
          icnt+=8;
        }
      } else {
        icnt=0;
#pragma omp for
        for(hi=hist;hi<hiend;hi++) {
          zz*=zzTmp[icnt];
          icnt++;
        }
      }
    }
    zz = clog(zz) ;
    z += zz ;
  }

  return cexp(z);
}

inline double complex LogRBMRatio(const double complex *rbmCntNew, const double complex *rbmCntOld) {
  int icnt;
  int idx,hi;
  int Nblk,iblk;
  int hist,hiend;
  double complex z=0.0, zz=1.0;
  double complex RbmCntNew[NBlockSize_RBMRatio], RbmCntOld[NBlockSize_RBMRatio], zzTmp[NBlockSize_RBMRatio];

  #pragma omp parallel for default(shared) private(hi) firstprivate(NRBM_PhysLayerIdx,Nneuron) reduction(+:z)
  for(idx=0;idx<NRBM_PhysLayerIdx;idx++) {
    z += RBM[idx]*(rbmCntNew[idx]-rbmCntOld[idx]);
  }

  Nblk = (Nneuron-1)/NBlockSize_RBMRatio + 1;
  for(iblk=0;iblk<Nblk;iblk++) { 
    hist = iblk*NBlockSize_RBMRatio;
    hiend= hist+NBlockSize_RBMRatio;
    if(hiend > Nneuron) hiend = Nneuron; 
    zz=1.0;
#pragma omp parallel default(shared) private(hi,RbmCntNew,RbmCntOld,zzTmp,icnt) firstprivate(NRBM_PhysLayerIdx,Nneuron) reduction(+:z) reduction(*:zz)
    {
      icnt=0;
#pragma loop swp
#pragma omp for
      for(hi=hist;hi<hiend;hi++) {
        RbmCntNew[icnt] = rbmCntNew[hi+NRBM_PhysLayerIdx];
        RbmCntOld[icnt] = rbmCntOld[hi+NRBM_PhysLayerIdx];
        if(creal(RbmCntNew[icnt]) <= 0.0 ) RbmCntNew[icnt] = -RbmCntNew[icnt];
        if(creal(RbmCntOld[icnt]) <= 0.0 ) RbmCntOld[icnt] = -RbmCntOld[icnt];
        z+=( RbmCntNew[icnt] - RbmCntOld[icnt] );
        zzTmp[icnt] = (1.0+cexp(-2.0*RbmCntNew[icnt])) / (1.0+cexp(-2.0*RbmCntOld[icnt]));
        icnt++;
      }

// Changing for SIMD tuning.
// If the loop length is multiple of 8, It is enabled SIMD.
// If the loop length is not multiple of 8, It is not SIMD.
      if((icnt%8)==0) {
        icnt=0;
#pragma omp for
        for(hi=hist;hi<hiend;hi=hi+8) {
          zz*=zzTmp[icnt]*zzTmp[icnt+1]*zzTmp[icnt+2]*zzTmp[icnt+3]*
              zzTmp[icnt+4]*zzTmp[icnt+5]*zzTmp[icnt+6]*zzTmp[icnt+7];
          icnt+=8;
        }
      } else {
        icnt=0;
#pragma omp for
        for(hi=hist;hi<hiend;hi++) {
          zz*=zzTmp[icnt];
          icnt++;
        }
      }
    }
    zz = clog(zz) ;
    z += zz ;
  }

  return z;
}


void MakeRBMCnt(double complex *rbmCnt, const int *eleNum) {
  const int *n=eleNum;
  const int *n0=eleNum; //up-spin
  const int *n1=eleNum+Nsite; //down-spin
  int idx;
  int ri,hi,hidx,xi;
  double complex ctmp;
  /* optimization for Kei */
  const int nRBM=NRBM_PhysLayerIdx + Nneuron;
  //const int nRBM=NRBM;
  const int nSite=Nsite;
  const int nSite2=Nsite2;
  const int nNeuronGeneral=NneuronGeneral;
  const int nNeuronCharge=NneuronCharge;
  const int nNeuronSpin=NneuronSpin;
  const double complex *RBM_Hidden     = RBM+NRBM_PhysLayerIdx;
  const double complex *RBM_PhysHidden = RBM+NRBM_PhysLayerIdx+NRBM_HiddenLayerIdx;
  int offset,offset2;

  /* initialization */
  for(idx=0;idx<nRBM;idx++) rbmCnt[idx] = 0.0;

  /* Potential on Physical Layer */
  if(NChargeRBM_PhysLayerIdx>0) {
    for(ri=0;ri<nSite;ri++) {
      rbmCnt[ ChargeRBM_PhysLayerIdx[ri] ] += n0[ri]+n1[ri]-1;
    }
  }
  if(NSpinRBM_PhysLayerIdx>0) {
    offset = NChargeRBM_PhysLayerIdx;
    for(ri=0;ri<nSite;ri++) {
      rbmCnt[ SpinRBM_PhysLayerIdx[ri] + offset] += n0[ri]-n1[ri];
    }
  }
  if(NGeneralRBM_PhysLayerIdx>0) {
    offset = NChargeRBM_PhysLayerIdx + NSpinRBM_PhysLayerIdx;
    for(ri=0;ri<nSite2;ri++) {
      rbmCnt[ GeneralRBM_PhysLayerIdx[ri] + offset] += 2*n0[ri]-1;
    }
  }

  /* Potential on Hidden Layer */
  if(NChargeRBM_HiddenLayerIdx>0) {
    for(hi=0;hi<nNeuronCharge;hi++) {
      hidx = ChargeRBM_HiddenLayerIdx[hi];
      rbmCnt[hi+NRBM_PhysLayerIdx] += RBM_Hidden[hidx];
    }
  }
  if(NSpinRBM_HiddenLayerIdx>0) {
    for(hi=0;hi<nNeuronSpin;hi++) {
      hidx = SpinRBM_HiddenLayerIdx[hi];
      rbmCnt[hi+NRBM_PhysLayerIdx + NneuronCharge] += RBM_Hidden[hidx + NChargeRBM_HiddenLayerIdx];
    }
  }
  if(NGeneralRBM_HiddenLayerIdx>0) {
    offset  = NRBM_PhysLayerIdx + NneuronCharge + NneuronSpin;
    offset2 = NChargeRBM_HiddenLayerIdx + NSpinRBM_HiddenLayerIdx;
    for(hi=0;hi<nNeuronGeneral;hi++) {
      hidx = GeneralRBM_HiddenLayerIdx[hi];
      rbmCnt[hi+offset] += RBM_Hidden[hidx + offset2];
    }
  }

  /* Coupling between Phys-Hidden Layers */
  if(NChargeRBM_PhysHiddenIdx>0) {
    for(hi=0;hi<nNeuronCharge;hi++) {
      ctmp = 0.0;
      for(ri=0;ri<nSite;ri++) {
        xi = n0[ri]+n1[ri]-1;
        hidx = ChargeRBM_PhysHiddenIdx[ri][hi];
        ctmp += RBM_PhysHidden[hidx]*xi;
      }
      rbmCnt[hi+NRBM_PhysLayerIdx] += ctmp;
    }
  }
  if(NSpinRBM_PhysHiddenIdx>0) {
    for(hi=0;hi<nNeuronSpin;hi++) {
      ctmp = 0.0;
      for(ri=0;ri<nSite;ri++) {
        xi = n0[ri]-n1[ri];
        hidx = SpinRBM_PhysHiddenIdx[ri][hi];
        ctmp += RBM_PhysHidden[hidx + NChargeRBM_PhysHiddenIdx]*xi;
      }
      rbmCnt[hi+NRBM_PhysLayerIdx + NneuronCharge] += ctmp;
    }
  }
  if(NGeneralRBM_PhysHiddenIdx>0) {
    offset  = NRBM_PhysLayerIdx + NneuronCharge + NneuronSpin;
    offset2 = NChargeRBM_PhysHiddenIdx+NSpinRBM_PhysHiddenIdx;
    for(hi=0;hi<nNeuronGeneral;hi++) {
      ctmp = 0.0;
      for(ri=0;ri<nSite2;ri++) {
        xi = 2*n[ri]-1;
        hidx = GeneralRBM_PhysHiddenIdx[ri][hi];
        ctmp += RBM_PhysHidden[hidx + offset2]*xi;
      }
      rbmCnt[hi+offset] += ctmp;
    }
  }
 
  return;
}

void copyFromBurnSampleRBM(double complex *rbmCnt) {
  int i,n;
  const double complex *burnRBMCnt = BurnRBMCnt;
  n = NRBM_PhysLayerIdx + Nneuron;
  #pragma loop noalias
  for(i=0;i<n;i++) rbmCnt[i] = burnRBMCnt[i];
  return;
}

void copyToBurnSampleRBM(double complex *rbmCnt) {
  int i,n;
  double complex *burnRBMCnt = BurnRBMCnt;
  n = NRBM_PhysLayerIdx + Nneuron;
  #pragma loop noalias
  for(i=0;i<n;i++) burnRBMCnt[i] = rbmCnt[i];
  return;
}

void saveRBMCnt(const int sample, const double complex *rbmCnt) {
  int i,offset;
  double complex x;
  const int n=NRBM_PhysLayerIdx + Nneuron;

  offset = sample*n;
  #pragma loop noalias
  for(i=0;i<n;i++) RBMCnt[offset+i] = rbmCnt[i];
  
  x = LogWeightRBM(rbmCnt);
  logSqPfFullSlater[sample] += 2.0*cabs(x);//TBC
  
  return;
}

void RBMDiff(double complex *srOptO, const double complex *rbmCnt, const int *eleNum) {
  int idx,hi,xi,ri;
  const int *n=eleNum; //up-spin
  const int *n0=eleNum; //up-spin
  const int *n1=eleNum+Nsite; //down-spin
  int offset  = 2*NRBM_PhysLayerIdx;
  int offset2 = 2*NRBM_PhysLayerIdx + 2*NRBM_HiddenLayerIdx;
  double complex ctmp=0;
  
  for(idx=0;idx<2*NRBM;idx++) {
    srOptO[idx] = 0.0;
  }

  for(idx=0;idx<NRBM_PhysLayerIdx;idx++) {
    ctmp = rbmCnt[idx];
    srOptO[2*idx]   =   ctmp;
    srOptO[2*idx+1] = I*ctmp;
  }

  /* Coupling between Phys-Hidden Layers */
  for(hi=0;hi<NneuronCharge;hi++) {
    idx = ChargeRBM_HiddenLayerIdx[hi];
    ctmp = ctanh( rbmCnt[hi+NRBM_PhysLayerIdx]);
    srOptO[2*idx + offset]     +=   ctmp;
    srOptO[2*idx + offset + 1] += I*ctmp;
    for(ri=0;ri<Nsite;ri++) {
      xi = n0[ri]+n1[ri]-1;
      idx = ChargeRBM_PhysHiddenIdx[ri][hi];
      srOptO[2*idx + offset2]     +=   xi*ctmp;
      srOptO[2*idx + offset2 + 1] += I*xi*ctmp;
    }
  }
  for(hi=0;hi<NneuronSpin;hi++) {
    idx = SpinRBM_HiddenLayerIdx[hi];
    ctmp = ctanh( rbmCnt[hi+NRBM_PhysLayerIdx+NneuronCharge] );
    srOptO[2*idx + offset + 2*NChargeRBM_HiddenLayerIdx]     +=   ctmp;
    srOptO[2*idx + offset + 2*NChargeRBM_HiddenLayerIdx + 1] += I*ctmp;
    for(ri=0;ri<Nsite;ri++) {
      xi = n0[ri]-n1[ri];
      idx = SpinRBM_PhysHiddenIdx[ri][hi];
      srOptO[2*idx + offset2 + 2*NChargeRBM_PhysHiddenIdx]     +=   xi*ctmp;
      srOptO[2*idx + offset2 + 2*NChargeRBM_PhysHiddenIdx + 1] += I*xi*ctmp;
    }
  }
  for(hi=0;hi<NneuronGeneral;hi++) {
    idx = GeneralRBM_HiddenLayerIdx[hi];
    ctmp = ctanh( rbmCnt[hi+NRBM_PhysLayerIdx+NneuronCharge+NneuronSpin] );
    srOptO[2*idx + offset + 2*NChargeRBM_HiddenLayerIdx + 2*NSpinRBM_HiddenLayerIdx]     +=   ctmp;
    srOptO[2*idx + offset + 2*NChargeRBM_HiddenLayerIdx + 2*NSpinRBM_HiddenLayerIdx + 1] += I*ctmp;
    for(ri=0;ri<Nsite2;ri++) {
      xi = 2*n[ri]-1;
      idx = GeneralRBM_PhysHiddenIdx[ri][hi];
      srOptO[2*idx + offset2 + 2*NChargeRBM_PhysHiddenIdx + 2*NSpinRBM_PhysHiddenIdx]     +=   xi*ctmp;
      srOptO[2*idx + offset2 + 2*NChargeRBM_PhysHiddenIdx + 2*NSpinRBM_PhysHiddenIdx + 1] += I*xi*ctmp;
    }
  }
 
    
  return;
}


/* An electron with spin s hops from ri to rj. */
void UpdateRBMCnt(const int ri, const int rj, const int s,
                   double complex *rbmCntNew, const double complex *rbmCntOld, const int *eleNum) {
//  const int *n0=eleNum; //up-spin
//  const int *n1=eleNum+Nsite; //down-spin
  int idx,offset;
  int hi,hidx,xi;
  double complex ctmp;
  /* optimization for Kei */
  const int nRBM=NRBM_PhysLayerIdx + Nneuron;
  const int nSite=Nsite;
  const int nSite2=Nsite2;
  const int nNeuronGeneral=NneuronGeneral;
  const int nNeuronCharge=NneuronCharge;
  const int nNeuronSpin=NneuronSpin;
  const double complex *RBM_PhysHidden = RBM+NRBM_PhysLayerIdx+NRBM_HiddenLayerIdx;
  const int rsi=ri+s*Nsite;
  const int rsj=rj+s*Nsite;
  const int t = 1-s;
  const int nit = eleNum[ri + t*Nsite];
  const int njt = eleNum[rj + t*Nsite];
  
  if(rbmCntNew != rbmCntOld) {
    for(idx=0;idx<nRBM;idx++) rbmCntNew[idx] = rbmCntOld[idx];
  }
  if(ri==rj) return;

  /* Potential on Physical Layer */
  if(NChargeRBM_PhysLayerIdx>0) {
    rbmCntNew[ ChargeRBM_PhysLayerIdx[ri] ] += -1;
    rbmCntNew[ ChargeRBM_PhysLayerIdx[rj] ] +=  1;
  }
  if(NSpinRBM_PhysLayerIdx>0) {
    rbmCntNew[ SpinRBM_PhysLayerIdx[ri] + NChargeRBM_PhysLayerIdx ] += 2*s-1;
    rbmCntNew[ SpinRBM_PhysLayerIdx[rj] + NChargeRBM_PhysLayerIdx ] += 1-2*s;
  }
  if(NGeneralRBM_PhysLayerIdx>0) {
    rbmCntNew[ GeneralRBM_PhysLayerIdx[rsi] + NChargeRBM_PhysLayerIdx + NSpinRBM_PhysLayerIdx] += -2;
    rbmCntNew[ GeneralRBM_PhysLayerIdx[rsj] + NChargeRBM_PhysLayerIdx + NSpinRBM_PhysLayerIdx] +=  2;
  }

  /* Potential on Hidden Layer */
  // No Change

  /* Coupling between Phys-Hidden Layers */
  if(NChargeRBM_PhysHiddenIdx>0) {
    for(hi=0;hi<nNeuronCharge;hi++) {
      hidx = ChargeRBM_PhysHiddenIdx[ri][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx] -= RBM_PhysHidden[hidx];
    }
    
    for(hi=0;hi<nNeuronCharge;hi++) {
      hidx = ChargeRBM_PhysHiddenIdx[rj][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx] += RBM_PhysHidden[hidx];
    }
  }
  if(NSpinRBM_PhysHiddenIdx>0) {
    for(hi=0;hi<nNeuronSpin;hi++) {
      hidx = SpinRBM_PhysHiddenIdx[ri][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge] -= (1-2*s)*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx];
    }
    
    for(hi=0;hi<nNeuronSpin;hi++) {
      hidx = SpinRBM_PhysHiddenIdx[rj][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge] += (1-2*s)*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx];
    }
  }
  if(NGeneralRBM_PhysHiddenIdx>0) {
    for(hi=0;hi<nNeuronGeneral;hi++) {
      hidx = GeneralRBM_PhysHiddenIdx[rsi][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge+nNeuronSpin] -= 2*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx+NSpinRBM_PhysHiddenIdx];
    }
    
    for(hi=0;hi<nNeuronGeneral;hi++) {
      hidx = GeneralRBM_PhysHiddenIdx[rsj][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge+nNeuronSpin] += 2*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx+NSpinRBM_PhysHiddenIdx];
    }
  }



  return;
}

/* An electron hops from ri,s to rj,t. */
void UpdateRBMCnt_fsz(const int ri, const int rj, const int s, const int t,
                   double complex *rbmCntNew, const double complex *rbmCntOld) {
//  const int *n0=eleNum; //up-spin
//  const int *n1=eleNum+Nsite; //down-spin
  int idx,offset;
  int hi,hidx,xi;
  double complex ctmp;
  /* optimization for Kei */
  const int nRBM=NRBM_PhysLayerIdx + Nneuron;
  const int nSite=Nsite;
  const int nNeuronCharge=NneuronCharge;
  const int nNeuronSpin=NneuronSpin;
  const int nNeuronGeneral=NneuronGeneral;
  const double complex *RBM_PhysHidden = RBM+NRBM_PhysLayerIdx+NRBM_HiddenLayerIdx;
  int rsi = ri +s*Nsite;
  int rtj = rj +t*Nsite;
  
  if(rbmCntNew != rbmCntOld) {
    for(idx=0;idx<nRBM;idx++) rbmCntNew[idx] = rbmCntOld[idx];
  }
  
  if(ri==rj){
    if(s == t){
      return;
    }else{
      if(NSpinRBM_PhysLayerIdx>0) {
        rbmCntNew[ SpinRBM_PhysLayerIdx[ri] + NChargeRBM_PhysLayerIdx ] += 2*s-1;
        rbmCntNew[ SpinRBM_PhysLayerIdx[rj] + NChargeRBM_PhysLayerIdx ] += 1-2*t;
      }
      if(NSpinRBM_PhysHiddenIdx>0) {
        for(hi=0;hi<nNeuronSpin;hi++) {
          hidx = SpinRBM_PhysHiddenIdx[ri][hi];
          rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge] -= (1-2*s)*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx];
        }
    
        for(hi=0;hi<nNeuronSpin;hi++) {
          hidx = SpinRBM_PhysHiddenIdx[rj][hi];
          rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge] += (1-2*t)*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx];
        }
      }
      if(NGeneralRBM_PhysHiddenIdx>0) {
        for(hi=0;hi<nNeuronGeneral;hi++) {
          hidx = GeneralRBM_PhysHiddenIdx[rsi][hi];
          rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge+nNeuronSpin] -= 2*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx+NSpinRBM_PhysHiddenIdx];
        }
    
        for(hi=0;hi<nNeuronGeneral;hi++) {
          hidx = GeneralRBM_PhysHiddenIdx[rtj][hi];
          rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge+nNeuronSpin] += 2*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx+NSpinRBM_PhysHiddenIdx];
        }
      }

      return;
    }
  }

  /* Potential on Physical Layer */
  if(NChargeRBM_PhysLayerIdx>0) {
    rbmCntNew[ ChargeRBM_PhysLayerIdx[ri] ] += -1;
    rbmCntNew[ ChargeRBM_PhysLayerIdx[rj] ] +=  1;
  }
  if(NSpinRBM_PhysLayerIdx>0) {
    rbmCntNew[ SpinRBM_PhysLayerIdx[ri] + NChargeRBM_PhysLayerIdx ] += 2*s-1;
    rbmCntNew[ SpinRBM_PhysLayerIdx[rj] + NChargeRBM_PhysLayerIdx ] += 1-2*t;
  }
  if(NGeneralRBM_PhysLayerIdx>0) {
    rbmCntNew[ GeneralRBM_PhysLayerIdx[rsi] + NChargeRBM_PhysLayerIdx + NSpinRBM_PhysLayerIdx] += -2;
    rbmCntNew[ GeneralRBM_PhysLayerIdx[rtj] + NChargeRBM_PhysLayerIdx + NSpinRBM_PhysLayerIdx] +=  2;
  }
  /* Potential on Hidden Layer */
  // No Change

  /* Coupling between Phys-Hidden Layers */
  if(NChargeRBM_PhysHiddenIdx>0) {
    for(hi=0;hi<nNeuronCharge;hi++) {
      hidx = ChargeRBM_PhysHiddenIdx[ri][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx] -= RBM_PhysHidden[hidx];
    }
    
    for(hi=0;hi<nNeuronCharge;hi++) {
      hidx = ChargeRBM_PhysHiddenIdx[rj][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx] += RBM_PhysHidden[hidx];
    }
  }
  if(NSpinRBM_PhysHiddenIdx>0) {
    for(hi=0;hi<nNeuronSpin;hi++) {
      hidx = SpinRBM_PhysHiddenIdx[ri][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge] -= (1-2*s)*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx];
    }
    
    for(hi=0;hi<nNeuronSpin;hi++) {
      hidx = SpinRBM_PhysHiddenIdx[rj][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge] += (1-2*t)*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx];
    }
  }
  if(NGeneralRBM_PhysHiddenIdx>0) {
    for(hi=0;hi<nNeuronGeneral;hi++) {
      hidx = GeneralRBM_PhysHiddenIdx[rsi][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge+nNeuronSpin] -= 2*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx+NSpinRBM_PhysHiddenIdx];
    }
    
    for(hi=0;hi<nNeuronGeneral;hi++) {
      hidx = GeneralRBM_PhysHiddenIdx[rtj][hi];
      rbmCntNew[hi+NRBM_PhysLayerIdx+nNeuronCharge+nNeuronSpin] += 2*RBM_PhysHidden[hidx+NChargeRBM_PhysHiddenIdx+NSpinRBM_PhysHiddenIdx];
    }
  }



  return;
}
