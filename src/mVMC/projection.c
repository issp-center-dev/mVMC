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
 * Gutzwiller-Jastrow Projection and DH correlation factors
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "global.h"
#include "math.h"
#include "projection.h"

inline double LogProjVal(const int *projCnt) {
  int idx;
  double z=0;
  for(idx=0;idx<NProj;idx++) {
    z += creal(Proj[idx]) * (double)(projCnt[idx]);
  }
  return z;
}

inline double LogProjRatio(const int *projCntNew, const int *projCntOld) {
  int idx;
  double z=0;
  for(idx=0;idx<NProj;idx++) {
    z += creal(Proj[idx]) * (double)(projCntNew[idx]-projCntOld[idx]); //TBC we assume gutzwiller and jastrow is real
  }
  return z;
}

inline double ProjRatio(const int *projCntNew, const int *projCntOld) {
  int idx;
  double z=0;
  for(idx=0;idx<NProj;idx++) {
    z += creal(Proj[idx]) * (double)(projCntNew[idx]-projCntOld[idx]); // TVC
  }
  return exp(z);
}
void MakeProjCnt(int *projCnt, const int *eleNum) {
  const int *n0=eleNum; //up-spin
  const int *n1=eleNum+Nsite; //down-spin
  int idx,offset;
  int ri,rj;
  int xi,xj,xn,xm;
  int r0,r1,r2,r3;
  const int *dh;
  /* optimization for Kei */
  const int nProj=NProj;
  const int nSite=Nsite;

  /* initialization */
  for(idx=0;idx<nProj;idx++) projCnt[idx] = 0;

  /* Gutzwiller factor */
  if(NGutzwillerIdx>0) {
    for(ri=0;ri<nSite;ri++) {
      projCnt[ GutzwillerIdx[ri] ] += n0[ri]*n1[ri];
    }
  }

  /* Jastrow factor exp(sum {v_ij * (ni-1) * (nj-1)}) */
  if(NJastrowIdx>0) {
    offset = NGutzwillerIdx;
    for(ri=0;ri<nSite;ri++) {
      xi = n0[ri]+n1[ri]-1;
      if(xi==0) continue;
      for(rj=ri+1;rj<nSite;rj++) {
        xj = n0[rj]+n1[rj]-1;
        idx = offset + JastrowIdx[ri][rj];
        projCnt[idx] += xi*xj;
      }
    }
  }

  /* 2-site doublon-holon correlation factor */
  if(NDoublonHolon2siteIdx>0) {
    offset = NGutzwillerIdx + NJastrowIdx;
    #pragma omp parallel for default(shared) private(xn,dh,ri,xi,r0,r1,xm,idx)
    for(xn=0;xn<NDoublonHolon2siteIdx;xn++) {
      dh=DoublonHolon2siteIdx[xn];
      for(ri=0;ri<nSite;ri++) {
        xi = n0[ri]+n1[ri];
        if(xi==1) continue;
        xi = xi/2; /* xi=0: holon, xi=1: doublon */

        r0 = dh[2*ri];
        r1 = dh[2*ri+1];
        if(xi==0) {
          xm = n0[r0]*n1[r0] + n0[r1]*n1[r1]; /* count doublons on r0 and r1 */
        } else { /* doublon */
          xm = (1-n0[r0])*(1-n1[r0]) 
            + (1-n0[r1])*(1-n1[r1]); /* count holons on r0 and r1 */
        }

        idx = offset + xn + (xi+2*xm)*NDoublonHolon2siteIdx;
        projCnt[idx] += 1;
      }
    }
  }

  /* 4-site doublon-holon correlation factor */
  if(NDoublonHolon4siteIdx>0) {
    offset = NGutzwillerIdx + NJastrowIdx + 6*NDoublonHolon2siteIdx;
    #pragma omp parallel for default(shared) private(xn,dh,ri,xi,r0,r1,r2,r3,xm,idx)
    for(xn=0;xn<NDoublonHolon4siteIdx;xn++) {
      dh=DoublonHolon4siteIdx[xn];
      for(ri=0;ri<nSite;ri++) {
        xi = n0[ri]+n1[ri];
        if(xi==1) continue;
        xi = xi/2; /* xi=0: holon, xi=1: doublon */

        r0 = dh[4*ri];
        r1 = dh[4*ri+1];
        r2 = dh[4*ri+2];
        r3 = dh[4*ri+3];
        if(xi==0) {
          /* count doublons on r0, r1, r2, r3 */
          xm= n0[r0]*n1[r0] + n0[r1]*n1[r1]
            + n0[r2]*n1[r2] + n0[r3]*n1[r3];
        } else { /* doublon */
          /* count holons on r0, r1, r2, r3 */
          xm= (1-n0[r0])*(1-n1[r0]) + (1-n0[r1])*(1-n1[r1])
            + (1-n0[r2])*(1-n1[r2]) + (1-n0[r3])*(1-n1[r3]);
        }
        
        idx = offset + xn + (xi+2*xm)*NDoublonHolon4siteIdx;
        projCnt[idx] += 1;
      }
    }
  }

  return;
}

/* An electron with spin s hops from ri to rj. */
void UpdateProjCnt(const int ri, const int rj, const int s,
                   int *projCntNew, const int *projCntOld,
                   const int *eleNum) {
  const int *n0=eleNum;
  const int *n1=eleNum+Nsite;
  int idx,offset;
  int rk;
  int xi,xn,xm;
  int r0,r1,r2,r3;
  const int *dh;
  /* optimization for Kei */
  const int nProj=NProj;
  const int nSite=Nsite;

  if(projCntNew!=projCntOld) {
    for(idx=0;idx<nProj;idx++) projCntNew[idx] = projCntOld[idx];
  }
  if(ri==rj) return;

  if(NGutzwillerIdx>0){
    idx = GutzwillerIdx[ri];
    projCntNew[idx] -= n0[ri]+n1[ri];
    idx = GutzwillerIdx[rj];
    projCntNew[idx] += n0[rj]*n1[rj];
  }

  if(NJastrowIdx>0){
    offset = NGutzwillerIdx; 
    /* update [ri][rj] */
    if(ri<rj) idx = offset + JastrowIdx[ri][rj];
    else  idx = offset + JastrowIdx[rj][ri];
    projCntNew[idx] += n0[ri]+n1[ri]-n0[rj]-n1[rj]+1;
    /* update [ri][rk] (rk != ri, rj) */
    for(rk=0;rk<nSite;rk++) {
      if(rk==rj) continue;
      if(rk==ri) continue;
      if(rk>ri) idx = offset + JastrowIdx[ri][rk];
      else idx = offset + JastrowIdx[rk][ri]; 
      projCntNew[idx] -= n0[rk]+n1[rk]-1;
    }
    /* update [rj][rk] (rk != ri, rj) */
    for(rk=0;rk<nSite;rk++) {
      if(rk==ri) continue;
      if(rk==rj) continue;
      if(rk>rj) idx = offset + JastrowIdx[rj][rk];
      else idx = offset + JastrowIdx[rk][rj]; 
      projCntNew[idx] += n0[rk]+n1[rk]-1;
    }
  }

  if(NDoublonHolon2siteIdx==0 && NDoublonHolon4siteIdx==0) return;

  offset = NGutzwillerIdx + NJastrowIdx;
  for(idx=offset;idx<nProj;idx++) projCntNew[idx] = 0;

  /* 2-site doublon-holon correlation factor */
  offset = NGutzwillerIdx + NJastrowIdx;
  #pragma omp parallel for default(shared) private(xn,dh,rk,xi,r0,r1,xm,idx)
  for(xn=0;xn<NDoublonHolon2siteIdx;xn++) {
    dh=DoublonHolon2siteIdx[xn];
    for(rk=0;rk<Nsite;rk++) {
      xi = n0[rk]+n1[rk];
      if(xi==1) continue;
      xi = xi/2; /* xi=0: holon, xi=1: doublon */
      
      r0 = dh[2*rk];
      r1 = dh[2*rk+1];
      if(xi==0) {
        xm = n0[r0]*n1[r0] + n0[r1]*n1[r1]; /* count doublons on r0 and r1 */
      } else { /* doublon */
        xm = (1-n0[r0])*(1-n1[r0]) 
          + (1-n0[r1])*(1-n1[r1]); /* count holons on r0 and r1 */
      }
      
      idx = offset + xn + (xi+2*xm)*NDoublonHolon2siteIdx;
      projCntNew[idx] += 1;
      
    }
  }

  /* 4-site doublon-holon correlation factor */
  offset = NGutzwillerIdx + NJastrowIdx + 6*NDoublonHolon2siteIdx;
  #pragma omp parallel for default(shared) private(xn,dh,rk,xi,r0,r1,r2,r3,xm,idx)
  for(xn=0;xn<NDoublonHolon4siteIdx;xn++) {
    dh=DoublonHolon4siteIdx[xn];
    for(rk=0;rk<Nsite;rk++) {
      xi = n0[rk]+n1[rk];
      if(xi==1) continue;
      xi = xi/2; /* xi=0: holon, xi=1: doublon */
      
      r0 = dh[4*rk];
      r1 = dh[4*rk+1];
      r2 = dh[4*rk+2];
      r3 = dh[4*rk+3];
      if(xi==0) {
        /* count doublons on r0, r1, r2, r3 */
        xm= n0[r0]*n1[r0] + n0[r1]*n1[r1]
          + n0[r2]*n1[r2] + n0[r3]*n1[r3];
      } else { /* doublon */
        /* count holons on r0, r1, r2, r3 */
        xm= (1-n0[r0])*(1-n1[r0]) + (1-n0[r1])*(1-n1[r1])
          + (1-n0[r2])*(1-n1[r2]) + (1-n0[r3])*(1-n1[r3]);
      }
      
      idx = offset + xn + (xi+2*xm)*NDoublonHolon4siteIdx;
      projCntNew[idx] += 1;
    }
  }
  
  return;

}

//[s] MERGE BY TM
/* An electron with spin s hops from ri to rj with t. */
// (ri,s) -> (rj,t) assuming s!=t
void UpdateProjCnt_fsz(const int ri, const int rj, const int s,const int t,
                   int *projCntNew, const int *projCntOld,
                   const int *eleNum) {
  const int *n0=eleNum;
  const int *n1=eleNum+Nsite;
  int idx,offset;
  int rk;
  int xi,xn,xm;
  int r0,r1,r2,r3;
  const int *dh;
  /* optimization for Kei */
  const int nProj=NProj;
  const int nSite=Nsite;

  if(projCntNew!=projCntOld) {
    for(idx=0;idx<nProj;idx++) projCntNew[idx] = projCntOld[idx];
  }

  if(ri==rj) return; // this part will be changed if we consider spin jastrow

  if(NGutzwillerIdx>0){
    idx = GutzwillerIdx[ri];
    projCntNew[idx] -= n0[ri]+n1[ri]; // fsz tricky !
    idx = GutzwillerIdx[rj];
    projCntNew[idx] += n0[rj]*n1[rj];
  }

  if(NJastrowIdx>0){
    offset = NGutzwillerIdx; 
    /* update [ri][rj] */
    if(ri<rj) idx = offset + JastrowIdx[ri][rj];
    else  idx = offset + JastrowIdx[rj][ri];
    projCntNew[idx] += n0[ri]+n1[ri]-n0[rj]-n1[rj]+1;
    /* update [ri][rk] (rk != ri, rj) */
    for(rk=0;rk<nSite;rk++) {
      if(rk==rj) continue;
      if(rk==ri) continue;
      if(rk>ri) idx = offset + JastrowIdx[ri][rk];
      else idx = offset + JastrowIdx[rk][ri]; 
      projCntNew[idx] -= n0[rk]+n1[rk]-1;
    }
    /* update [rj][rk] (rk != ri, rj) */
    for(rk=0;rk<nSite;rk++) {
      if(rk==ri) continue;
      if(rk==rj) continue;
      if(rk>rj) idx = offset + JastrowIdx[rj][rk];
      else idx = offset + JastrowIdx[rk][rj]; 
      projCntNew[idx] += n0[rk]+n1[rk]-1;
    }
  }

  if(NDoublonHolon2siteIdx==0 && NDoublonHolon4siteIdx==0) return;

  offset = NGutzwillerIdx + NJastrowIdx;
  for(idx=offset;idx<nProj;idx++) projCntNew[idx] = 0;

  /* 2-site doublon-holon correlation factor */
  offset = NGutzwillerIdx + NJastrowIdx;
  #pragma omp parallel for default(shared) private(xn,dh,rk,xi,r0,r1,xm,idx)
  for(xn=0;xn<NDoublonHolon2siteIdx;xn++) {
    dh=DoublonHolon2siteIdx[xn];
    for(rk=0;rk<Nsite;rk++) {
      xi = n0[rk]+n1[rk];
      if(xi==1) continue;
      xi = xi/2; /* xi=0: holon, xi=1: doublon */
      
      r0 = dh[2*rk];
      r1 = dh[2*rk+1];
      if(xi==0) {
        xm = n0[r0]*n1[r0] + n0[r1]*n1[r1]; /* count doublons on r0 and r1 */
      } else { /* doublon */
        xm = (1-n0[r0])*(1-n1[r0]) 
          + (1-n0[r1])*(1-n1[r1]); /* count holons on r0 and r1 */
      }
      
      idx = offset + xn + (xi+2*xm)*NDoublonHolon2siteIdx;
      projCntNew[idx] += 1;
      
    }
  }

  /* 4-site doublon-holon correlation factor */
  offset = NGutzwillerIdx + NJastrowIdx + 6*NDoublonHolon2siteIdx;
  #pragma omp parallel for default(shared) private(xn,dh,rk,xi,r0,r1,r2,r3,xm,idx)
  for(xn=0;xn<NDoublonHolon4siteIdx;xn++) {
    dh=DoublonHolon4siteIdx[xn];
    for(rk=0;rk<Nsite;rk++) {
      xi = n0[rk]+n1[rk];
      if(xi==1) continue;
      xi = xi/2; /* xi=0: holon, xi=1: doublon */
      
      r0 = dh[4*rk];
      r1 = dh[4*rk+1];
      r2 = dh[4*rk+2];
      r3 = dh[4*rk+3];
      if(xi==0) {
        /* count doublons on r0, r1, r2, r3 */
        xm= n0[r0]*n1[r0] + n0[r1]*n1[r1]
          + n0[r2]*n1[r2] + n0[r3]*n1[r3];
      } else { /* doublon */
        /* count holons on r0, r1, r2, r3 */
        xm= (1-n0[r0])*(1-n1[r0]) + (1-n0[r1])*(1-n1[r1])
          + (1-n0[r2])*(1-n1[r2]) + (1-n0[r3])*(1-n1[r3]);
      }
      
      idx = offset + xn + (xi+2*xm)*NDoublonHolon4siteIdx;
      projCntNew[idx] += 1;
    }
  }
  
  return;
}

void MakeProjBFCnt(int *projCnt, const int *eleNum) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  int idx;
  int k;
  int ri, rk;
  int xid, xih, xkd, xkh;
  int xidh, xihd, xkdh, xkhd;
  int *nBF0, *nBF1, *nBF2, *nBF3;
  /* optimization for Kei */
  const int nSite = Nsite;
  const int nRange = Nrange;

  /* initialization */
  for (idx = 0; idx < 16 * nSite * nRange; idx++) projCnt[idx] = 0;

  nBF0 = projCnt;
  nBF1 = projCnt + 4 * nSite * nRange;
  nBF2 = projCnt + 8 * nSite * nRange;
  nBF3 = projCnt + 12 * nSite * nRange;

  /* BackFlow factor */
  if (NBackFlowIdx > 0) {
    for (ri = 0; ri < nSite; ri++) {
      xid = n0[ri] * n1[ri];
      xih = (1 - n0[ri]) * (1 - n1[ri]);
      xidh = n0[ri] * (1 - n1[ri]);
      xihd = n1[ri] * (1 - n0[ri]);
      //posBF = PosBF + ri*nSite;
      for (k = 0; k < nRange; k++) {
        rk = PosBF[ri][k];
        xkd = n0[rk] * n1[rk];
        xkh = (1 - n0[rk]) * (1 - n1[rk]);
        xkdh = n0[rk] * (1 - n1[rk]);
        xkhd = n1[rk] * (1 - n0[rk]);

        if (ri == rk) {
          nBF0[ri * nRange + k] = 1;
          nBF1[ri * nRange + k] = 1;
          nBF2[ri * nRange + k] = 1;
          nBF3[ri * nRange + k] = 1;
        } else {
          nBF0[ri * nRange + k] = 0;
          nBF1[ri * nRange + k] = 0;
          nBF2[ri * nRange + k] = 0;
          nBF3[ri * nRange + k] = 0;
        }
        nBF0[nSite * nRange + ri * nRange + k] = xid * xkh;
        nBF0[2 * nSite * nRange + ri * nRange + k] = xidh * xkhd;
        nBF0[3 * nSite * nRange + ri * nRange + k] = xid * xkhd + xidh * xkh;

        nBF2[nSite * nRange + ri * nRange + k] = xkd * xih;
        nBF2[2 * nSite * nRange + ri * nRange + k] = xkdh * xihd;
        nBF2[3 * nSite * nRange + ri * nRange + k] = xkd * xihd + xkdh * xih;

        nBF1[nSite * nRange + ri * nRange + k] = xid * xkh;
        nBF1[2 * nSite * nRange + ri * nRange + k] = xihd * xkdh;
        nBF1[3 * nSite * nRange + ri * nRange + k] = xid * xkdh + xihd * xkh;

        nBF3[nSite * nRange + ri * nRange + k] = xkd * xih;
        nBF3[2 * nSite * nRange + ri * nRange + k] = xkhd * xidh;
        nBF3[3 * nSite * nRange + ri * nRange + k] = xkd * xidh + xkhd * xih;
      }
    }
  }

  return;
}
//[e] MERGE BY TM
