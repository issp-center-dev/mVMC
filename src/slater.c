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
 * slater elements
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include <complex.h>
#include "./include/global.h"
#include "./include/slater.h"

#pragma once
void SubSlaterElmBF_fcmp(const int tri, const int trj, double complex*slt_ij, int *ijcount, double complex*slt_ji, int *jicount, const int *eleProjBFCnt);

void SubSlaterElmBF_real(const int tri, const int trj, double *slt_ij, int *ijcount, double *slt_ji, int *jicount, const int *eleProjBFCnt);

void UpdateSlaterElm_fcmp() {
  int ri,ori,tri,sgni,rsi0,rsi1;
  int rj,orj,trj,sgnj,rsj0,rsj1;
  int qpidx,mpidx,spidx,optidx;
  double cs,cc,ss;
  double complex slt_ij,slt_ji;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex *sltE,*sltE_i0,*sltE_i1;

  #pragma omp parallel for default(shared)        \
    private(qpidx,optidx,mpidx,spidx,                      \
            xqpOpt,xqpOptSgn,xqp,xqpSgn,cs,cc,ss,sltE,     \
            ri,ori,tri,sgni,rsi0,rsi1,sltE_i0,sltE_i1,      \
            rj,orj,trj,sgnj,rsj0,rsj1,slt_ij,slt_ji)
  #pragma loop noalias
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    // qpidx  = optidx*NQPFix+NSPGaussLeg*mpidx+spidx 
    // NQPFix = NSPGaussLeg*NMPTrans
    optidx    = qpidx / NQPFix;                // optidx -> optrans projection (will not be used ?)
    mpidx     = (qpidx%NQPFix) / NSPGaussLeg;  // mpidx  -> momentum projection
    spidx     = qpidx % NSPGaussLeg;           // spidx  -> spin     projection

    xqpOpt    = QPOptTrans[optidx];    //
    xqpOptSgn = QPOptTransSgn[optidx]; //
    xqp       = QPTrans[mpidx];        //QPTrans[i][j]: i # of trans. op.,j origin of trans. op. 
    xqpSgn    = QPTransSgn[mpidx];
    cs        = SPGLCosSin[spidx];
    cc        = SPGLCosCos[spidx];
    ss        = SPGLSinSin[spidx];
    
    sltE      = SlaterElm + qpidx*Nsite2*Nsite2;
    
    for(ri=0;ri<Nsite;ri++) {
      ori     = xqpOpt[ri];           // ri (OptTrans) -> ori (Trans)-> tri
      tri     = xqp[ori];
      sgni    = xqpSgn[ori]*xqpOptSgn[ri];
      rsi0    = ri;
      rsi1    = ri+Nsite;
      sltE_i0 = sltE + rsi0*Nsite2;
      sltE_i1 = sltE + rsi1*Nsite2;
      
      for(rj=0;rj<Nsite;rj++) {
        orj = xqpOpt[rj];
        trj = xqp[orj];
        sgnj = xqpSgn[orj]*xqpOptSgn[rj];
        rsj0 = rj;
        rsj1 = rj+Nsite;
        
        slt_ij = Slater[ OrbitalIdx[tri][trj] ] * (double)(OrbitalSgn[tri][trj]*sgni*sgnj);
        slt_ji = Slater[ OrbitalIdx[trj][tri] ] * (double)(OrbitalSgn[trj][tri]*sgni*sgnj);
        
        sltE_i0[rsj0] = -(slt_ij - slt_ji)*cs;   // up   - up
        sltE_i0[rsj1] =   slt_ij*cc + slt_ji*ss; // up   - down
        sltE_i1[rsj0] = -slt_ij*ss - slt_ji*cc;  // down - up
        sltE_i1[rsj1] =  (slt_ij - slt_ji)*cs;   // down - down 
      }
    }
  }

  return;
}

// Calculating Tr[Inv[M]*D_k(X)]
void SlaterElmDiff_fcmp(double complex *srOptO, const double complex ip, int *eleIdx) {
  const int nBuf=NSlater*NQPFull;
  const int nsize = Nsize;
  const int ne = Ne;
  const int nQPFull = NQPFull;
  const int nMPTrans = NMPTrans; // number of translational operators
  const int nSlater = NSlater;
  const int nTrans = NMPTrans * NQPOptTrans; //usually NQPOptTrans=1

  const double complex invIP = 1.0/ip;
  int msi,msj,ri,rj,ori,orj,tri,trj,sgni,sgnj;
  int mpidx,spidx,orbidx,qpidx,optidx,i;
  double complex cs,cc,ss; // including Pf
  int *xqp,*xqpSgn,*xqpOpt,*xqpOptSgn;
  double complex *invM,*invM_i;

  int *orbitalIdx_i;
  int *transOrbIdx; /* transOrbIdx[mpidx][msi][msj] */
  int *transOrbSgn; /* transOrbSgn[mpidx][msi][msj] */
  int *tOrbIdx,*tOrbIdx_i;
  int *tOrbSgn,*tOrbSgn_i;
  double complex *buf, *buffer;
  double complex tmp;

  RequestWorkSpaceInt(2*nTrans*Nsize*Nsize);
  RequestWorkSpaceComplex(NQPFull*NSlater);

  transOrbIdx = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbIdx[mpidx][msi][msj] */
  transOrbSgn = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbSgn[mpidx][msi][msj] */
  buffer = GetWorkSpaceComplex(NQPFull*NSlater);

  for(i=0;i<nBuf;i++) buffer[i]=0.0;

  #pragma omp parallel for default(shared)                        \
    private(qpidx,optidx,mpidx,msi,msj,xqp,xqpSgn,xqpOpt,xqpOptSgn,\
            ri,ori,tri,sgni,rj,orj,trj,sgnj,                       \
            tOrbIdx,tOrbIdx_i,tOrbSgn,tOrbSgn_i,orbitalIdx_i)
  #pragma loop noalias
  for(qpidx=0;qpidx<nTrans;qpidx++) {            // nTran=nMPTrans*NQPOptTrans: usually nMPTrans
    optidx    = qpidx / nMPTrans;                // qpidx=mpidx+nMPTrans*optidx
    mpidx     = qpidx % nMPTrans; 
                                                 // QPOptTrans:     NQPOptTrans* Nsite matrix :  usually NQPOptTrans = 1
                                                 // QPOptTransSgn:  NQPOptTrans* Nsite matrix :  usually NQPOptTrans = 1 
    xqpOpt    = QPOptTrans[optidx];              // xqpOpt[]    = QPOptTrans[optidx][]    : usually optidx=0 and not be used
    xqpOptSgn = QPOptTransSgn[optidx];           // xqpOptSgn[] = QPOptTransSgn[optidx][] : usually optidx=0 and not be used
    xqp       = QPTrans[mpidx];                  // xqp    =  QPTrans[mpidx][]
    xqpSgn    = QPTransSgn[mpidx];               // xqpSgn =  QPTransSgn[mpidx][]
    tOrbIdx   = transOrbIdx + qpidx*nsize*nsize; // tOrbIdx : f_ij
    tOrbSgn   = transOrbSgn + qpidx*nsize*nsize; // tOrbSgn : sign of f_ij
    for(msi=0;msi<nsize;msi++) {                 // nsize=2*Ne
      ri           = eleIdx[msi];                //  ri  : postion where the msi-th electron exists  
      ori          = xqpOpt[ri];                 // ori  : 
      tri          = xqp[ori];                   // tri  : 
      sgni         = xqpSgn[ori]*xqpOptSgn[ri];  // sgni :
      tOrbIdx_i    = tOrbIdx + msi*nsize;        // 
      tOrbSgn_i    = tOrbSgn + msi*nsize;        //
      orbitalIdx_i = OrbitalIdx[tri];            //
      for(msj=0;msj<nsize;msj++) { //nsize=2*Ne //
        rj             = eleIdx[msj];           //
        orj            = xqpOpt[rj];         
        trj            = xqp[orj];
        sgnj           = xqpSgn[orj]*xqpOptSgn[rj];
        tOrbIdx_i[msj] = orbitalIdx_i[trj];
        tOrbSgn_i[msj] = sgni*sgnj*OrbitalSgn[tri][trj];
      }
    }
  }
// calculating Tr(X^{-1}*dX/df_{msi,msj})=-2*alpha(sigma(msi),sigma(msj))(X^{-1})_{msi,msj}
  #pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,spidx,cs,cc,ss,                   \
            tOrbIdx,tOrbSgn,invM,buf,msi,msj,             \
            tOrbIdx_i,tOrbSgn_i,invM_i,orbidx)
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) { // nQPFull = NQPFix * NQPOptTrans: usually = NSPGaussLeg * NMPTrans
    mpidx = qpidx / NSPGaussLeg;       // qpidx   = NSPGaussLeg*mpidx + spidx
    spidx = qpidx % NSPGaussLeg;

    cs = PfM[qpidx] * SPGLCosSin[spidx]; // spin rotation + PfM
    cc = PfM[qpidx] * SPGLCosCos[spidx]; // spin rotation + PfM
    ss = PfM[qpidx] * SPGLSinSin[spidx]; // spin rotation + PfM

    tOrbIdx = transOrbIdx + mpidx*nsize*nsize; // tOrbIdx[msi][msj] = transOrbIdx[mpidx][msi][msj]
    tOrbSgn = transOrbSgn + mpidx*nsize*nsize; // tOrbSgn[msi][msj] = transOrbSgn[mpidx][msi][msj]
    invM    = InvM        + qpidx*Nsize*Nsize; // invM  M^-1 and PfM, InvM: NQPFull*(Nsize*Nsize)+PfM
    buf     = buffer      + qpidx*NSlater;     // buf   fij, buffer: NSlater*NQPFull

    #pragma loop norecurrence
    for(msi=0;msi<ne;msi++) {
      tOrbIdx_i = tOrbIdx + msi*nsize;         // tOrbIdx_i[] = tOrbIdx[msi][msj]
      tOrbSgn_i = tOrbSgn + msi*nsize;         // tOrbSgn_i[] = tOrbSgn[msi][msj]
      invM_i    = invM + msi*nsize;            // invM[]      = invM[msi][]
      for(msj=0;msj<ne;msj++) {                // up-up
        /* si=0 sj=0*/
        orbidx       = tOrbIdx_i[msj];         // 
        buf[orbidx] += invM_i[msj]*cs*tOrbSgn_i[msj]; // invM[msi][msj]
      }
      for(msj=ne;msj<nsize;msj++) {            // up-down
        /* si=0 sj=1*/
        orbidx       = tOrbIdx_i[msj];
        buf[orbidx] -= invM_i[msj]*cc*tOrbSgn_i[msj];
      }
    }
    #pragma loop norecurrence
    for(msi=ne;msi<nsize;msi++) { 
      tOrbIdx_i = tOrbIdx + msi*nsize;
      tOrbSgn_i = tOrbSgn + msi*nsize;
      invM_i = invM + msi*nsize;
      for(msj=0;msj<ne;msj++) {    // down-up
        /* si=1 sj=0*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] += invM_i[msj]*ss*tOrbSgn_i[msj];
      }
      for(msj=ne;msj<nsize;msj++) {// down-down
        /* si=1 sj=1*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] -= invM_i[msj]*cs*tOrbSgn_i[msj];
      }
    }
  }

  /* store SROptO[] */
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[2*orbidx]   = 0.0+0.0*I; // 0
    srOptO[2*orbidx+1] = 0.0+0.0*I; // 0
  }
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    tmp = QPFullWeight[qpidx];
    buf = buffer + qpidx*nSlater;
    for(orbidx=0;orbidx<nSlater;orbidx++) {
      srOptO[2*orbidx]   += tmp * buf[orbidx];   //real      TBC
      srOptO[2*orbidx+1] += tmp * buf[orbidx]*I; //imaginary TBC
      //printf("Re DEBUG: tmp=%lf :orbidx=%d srOptO=%lf %lf invIP=%lf %lf \n",tmp,orbidx,creal(srOptO[2*orbidx]),cimag(srOptO[2*orbidx]),creal(invIP),cimag(invIP));
      //printf("Im DEBUG:orbidx=%d srOptO=%lf %lf invIP=%lf %lf \n",orbidx,creal(srOptO[2*orbidx+1]),cimag(srOptO[2*orbidx+1]),creal(invIP),cimag(invIP));
    }
  }
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[2*orbidx]   *= invIP;
    srOptO[2*orbidx+1] *= invIP;
  }

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  return;
}

void SlaterElmBFDiff_fcmp(double complex*srOptO, const double complex ip, int *eleIdx, int *eleNum, int *eleCfg, int *eleProjConst,const int * eleProjBFCnt){
  const int nBuf=NSlater*NQPFull;
  const int nsize = Nsize;
  const int ne = Ne;
  const int nQPFull = NQPFull;
  const int nMPTrans = NMPTrans;
  const int nSlater = NSlater;
  const int nTrans = NMPTrans * NQPOptTrans;
  const double complex invIP = 1.0/ip;
  
  int xi,xj,xk,xl,xn,xm,xid,xih,xjd,xjh,xkd,xkh,xld;
  int xidh,xihd,xjhd,xjdh,xkdh,xkhd,xldh,xlhd;
  int s,idx,idx2,idx3,idx4,rtmp,rsk,rstmp,rsi,rsj,msi0,msi1,msj0,msj1;
  //double eta;
  const int *n0=eleNum;
  const int *n1=eleNum+Nsite;
  int mi,mj,msi,msj,msk,msl,ri,rj,rk,rl;
  int tri,trj,trk,trl;
  int mpidx,spidx,orbidx,orbidx2,orbsgn,qpidx,i;
  double complex cs,cc,ss;
  int *xqp, *xqpInv;
  double complex *invM,*invM_i,*invM_j,*invM_k,*invM_l;
//  const int *bfCnt0=eleProjBFCnt;
//  const int *bfCnt1=eleProjBFCnt + 4*Nsite*Nrange;
  const int *bfCnt2=eleProjBFCnt + 8*Nsite*Nrange;
  const int *bfCnt3=eleProjBFCnt +12*Nsite*Nrange;
  const int *bfCnt2_n,*bfCnt3_m;
  int **posBF = PosBF;

  int *orbitalIdx_i,*orbitalSgn_i;
  int *transOrbIdx; /* transOrbIdx[mpidx][msi][msj] */
  int *transOrbSgn; /* transOrbSgn[mpidx][msi][msj] */
  int *tOrbIdx,*tOrbIdx_i;
  int *tOrbSgn_i;
  double complex *buf,*buffer;
  double complex tmp;
  //double complex bfdhidx[NMultiSlater*NQPFull*Nsize*Nsize];
  //double complex *bfdhidx_i;
  double complex pTrans[Nsite2*Nsite2];
  double complex *pTrans_i;

  double complex **etaTmp;
  int rki,rlj;
  const int nSite=Nsite;
  const int nRange=Nrange;
  const int nSiteRange = nRange*nSite;
  int idx_ik,idx_jl;
  int bfidx;
  int dki,dlj,nidx,midx,xtmp;

  RequestWorkSpaceInt(2*nTrans*Nsize*Nsize);
  RequestWorkSpaceComplex(NQPFull*NSlater);

  transOrbIdx = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbIdx[mpidx][msi][msj] */
  transOrbSgn = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbSgn[mpidx][msi][msj] */
  buffer = GetWorkSpaceComplex(NQPFull*NSlater);
  
  for(i=0;i<nBuf;i++) buffer[i]=0.0;

  #pragma omp parallel for default(shared)        \
    private(mpidx,msi,msj,xqp,ri,tri,rj,trj,      \
            tOrbIdx,tOrbIdx_i,orbitalIdx_i)
  #pragma loop noalias
  for(mpidx=0;mpidx<nMPTrans;mpidx++) {
    xqp = QPTrans[mpidx];
    tOrbIdx = transOrbIdx + mpidx*nsize*nsize;
    for(msi=0;msi<nsize;msi++) {
      ri = eleIdx[msi];
      tri = xqp[ri];
      tOrbIdx_i = tOrbIdx + msi*nsize;
      orbitalIdx_i = OrbitalIdx[tri];
      for(msj=0;msj<nsize;msj++) {
        rj = eleIdx[msj];
        trj = xqp[rj];
        tOrbIdx_i[msj] = orbitalIdx_i[trj];
      }
    }
  }
            
  for(idx=0;idx<Nsite2*Nsite2;idx++){
     pTrans[idx]=0.0;
  }
  for(idx=0;idx<NTransfer;idx++){
    ri = Transfer[idx][0];
    rj = Transfer[idx][2];
    s  = Transfer[idx][1];// C_ris a_rjs
    rsi = ri + Nsite*s;
    rsj = rj + Nsite*s;
    pTrans_i = pTrans + rsi*Nsite2;
    pTrans_i[rsj] = ParaTransfer[idx];
  }
  /*for(idx=0;idx<Nsite2*Nsite2;idx++){
    rsi = idx/Nsite2;
    rsj = idx%Nsite2;
    printf("Transfer[%d][%d]=%.5e,%.5e\n",rsi,rsj,creal(pTrans[idx]),cimag(pTrans[idx]));
  }*/

  #pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,spidx,cs,cc,ss,           \
            tOrbIdx,invM,buf,msi,msj,             \
            ri,rj,rk,tri,trj,idx,xkd,xkh,xid,xjd, \
            tOrbIdx_i,invM_i,orbidx)
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    mpidx = qpidx / NSPGaussLeg;
    spidx = qpidx % NSPGaussLeg;
    xqp = QPTrans[mpidx];
    xqpInv = QPTransInv[mpidx];

    cs = PfM[qpidx] * SPGLCosSin[spidx];
    cc = PfM[qpidx] * SPGLCosCos[spidx];
    ss = PfM[qpidx] * SPGLSinSin[spidx];

    tOrbIdx = transOrbIdx + mpidx*nsize*nsize;
//    tOrbSgn = transOrbSgn + mpidx*nsize*nsize;
    invM = InvM + qpidx*Nsize*Nsize;
    buf = buffer + qpidx*NSlater;
    etaTmp = eta + qpidx*Nsite*Nsite;
    //flagTmp = EtaFlag + qpidx*Nsite*Nsite;

       // if(slt_ij == 0.0){eta = 1.0;}
       // else{eta = ProjBF[0];}
        //slt_ij += eta*Slater[ OrbitalIdx[tri][trj] ];
        
    //BackFlow Correlation eta1
    if(NBackFlowIdx > 0){
      #pragma loop norecurrence
      //for(msi=0;msi<ne;msi++) {
      for(ri=0;ri<Nsite;ri++) {
        mi = eleCfg[ri];
        msi = mi;
        tri = xqp[ri];
        orbitalIdx_i = OrbitalIdx[tri];
        orbitalSgn_i = OrbitalSgn[tri];
        for(rj=0;rj<Nsite;rj++) {
          // si=0 sj=1
          mj = eleCfg[rj+Nsite];
          msj = mj + ne;
          trj = xqp[rj];
          orbidx = orbitalIdx_i[trj];
          orbsgn = orbitalSgn_i[trj];

          if((mi != -1) && (mj != -1)){
            invM_i = invM + msi*nsize;
            //buf[orbidx] -= tEta[tri][trj]*invM_i[msj]*cc;
            //printf("eta[%d][%d]",tri,trj);
            //printf("=%.2e\n",etaTmp[tri][trj]);
            buf[orbidx] -= etaTmp[tri][trj]*invM_i[msj]*cc*orbsgn;
            //buf[orbidx] -= etaTmp[tri][trj]*invM_i[msj]*cc*tOrbSgn_i[msj];
          }

          for(xn=0;xn<4;xn++){
            bfCnt2_n=bfCnt2+xn*nSiteRange;
//            bfCnt3_n=bfCnt3+xn*nSiteRange;
            for(xm=0;xm<4;xm++){
              if(xm==0 && xn == 0) continue;
              bfidx=BFSubIdx[xn][xm];
//              bfCnt2_m=bfCnt2+xm*nSiteRange;
              bfCnt3_m=bfCnt3+xm*nSiteRange;
          
              for(xk=0;xk<nRange;xk++) {
                rki=posBF[tri][xk];
                rk=xqpInv[rki];
                if(eleCfg[rk] == -1) continue;
                msk = eleCfg[rk];
                idx_ik=tri*nRange+xk;
                invM_k = invM + msk*nsize;
              
                dki = RangeIdx[tri][rki];
                xtmp = 4*dki+xn;
                nidx = xtmp-3-dki;
                //if(itmp <= 3 || itmp%4==0){continue;}
                if(xtmp%4==0){nidx=-1;}
                if(xtmp==0){nidx=0;}
                if(nidx<0){continue;}

                for(xl=0;xl<nRange;xl++) {
                  rlj=posBF[trj][xl];
                  rl=xqpInv[rlj];
                  if(eleCfg[rl+Nsite] == -1) continue;
                  msl = eleCfg[rl+Nsite]+ne;
                  idx_jl=trj*nRange+xl;
                
                  dlj = RangeIdx[trj][rlj];
                  xtmp = 4*dlj+xm;
                  midx = xtmp-3-dlj;
                  //if(itmp <= 3 || itmp%4==0){continue;}
                  if(xtmp%4==0){midx=-1;}
                  if(xtmp==0){midx=0;}
                  if(midx<0){continue;}
                  bfidx=BFSubIdx[nidx][midx];
                  orbsgn =  OrbitalSgn[rki][rlj];

                  buf[orbidx] -= -orbsgn*invM_k[msl]*ProjBF[bfidx]*bfCnt2_n[idx_ik]*bfCnt3_m[idx_jl]*PfM[qpidx];
                }
              }
            }
          }
        }
      }
      /*#pragma loop norecurrence
      for(msi=ne;msi<nsize;msi++) {
        tOrbIdx_i = tOrbIdx + msi*nsize;
        invM_i = invM + msi*nsize;
        for(msj=0;msj<ne;msj++) {
          // si=1 sj=0
          orbidx  = tOrbIdx_i[msj];
          buf[orbidx] += invM_i[msj]*cs*tOrbSgn_i[msj];
        }
        for(msj=ne;msj<nsize;msj++) {
          // si=1 sj=1
          orbidx = tOrbIdx_i[msj];
          buf[orbidx] -= invM_i[msj]*cc*tOrbSgn_i[msj];
        }
      }*/
    }else{
      #pragma loop norecurrence
      for(msi=0;msi<ne;msi++) {
        tOrbIdx_i = tOrbIdx + msi*nsize;
        invM_i = invM + msi*nsize;
        for(msj=0;msj<ne;msj++) {
          /* si=0 sj=0*/
          orbidx = tOrbIdx_i[msj];
          buf[orbidx] += invM_i[msj]*ss*tOrbSgn_i[msj];
        }
        for(msj=ne;msj<nsize;msj++) {
          /* si=0 sj=1*/
          orbidx = tOrbIdx_i[msj];
          buf[orbidx] -= invM_i[msj]*cs*tOrbSgn_i[msj];
        }
      }
      #pragma loop norecurrence
      for(msi=ne;msi<nsize;msi++) {
        tOrbIdx_i = tOrbIdx + msi*nsize;
        invM_i = invM + msi*nsize;
        for(msj=0;msj<ne;msj++) {
          /* si=1 sj=0*/
          orbidx = tOrbIdx_i[msj];
          buf[orbidx] += invM_i[msj]*ss*tOrbSgn_i[msj];
        }
        for(msj=ne;msj<nsize;msj++) {
          /* si=1 sj=1*/
          orbidx = tOrbIdx_i[msj];
          buf[orbidx] -= invM_i[msj]*cs*tOrbSgn_i[msj];
        }
      }
    }
  }

  /* store SROptO[] */
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[2*orbidx] = 0.0+0.0*I;
    srOptO[2*orbidx+1] = 0.0+0.0*I;
  }
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    tmp = QPFullWeight[qpidx];
    buf = buffer + qpidx*nSlater;
    for(orbidx=0;orbidx<nSlater;orbidx++) {
      srOptO[2*orbidx] += tmp * buf[orbidx];
      srOptO[2*orbidx+1] += tmp * buf[orbidx]*I;
    }
  }
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[2*orbidx] *= invIP;
    srOptO[2*orbidx+1] *= invIP;
  }

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  return;
}

/* buffer size = NQPFull*NSlater */
/* bufInt size = NMPTrans*Nsize*Nsize */
void BackFlowDiff_fcmp(complex double *srOptO, const double complex ip, int *eleIdx, const int *eleNum, int *eleProjConst,
                   const int *eleProjBFCnt) {

  const double complex invIP = 1.0/ip;
  int s,t,idx,idx2,rk,rtmp,rsk,rstmp,rl,rtmp2,rtl,rttmp,rsi,rsj;
  double complex t_ki, t_lj;
  //double eta;
  const int *n0=eleNum;
  const int *n1=eleNum+Nsite;
  int msi,msj,ri,rj,tri,trj;
  int mpidx,spidx,orbidx,qpidx,i;
  double complex cs,cc,ss;
  int *xqp;
  const double complex *invM,*invM_i;

  int *orbitalIdx_i;
  int *transOrbIdx, *transOrbSgn; /* transOrbIdx[mpidx][msi][msj] */
  int *tOrbIdx,*tOrbIdx_i;
  int *tOrbSgn,*tOrbSgn_i;
  double complex *buf,*buf_i;
  double complex logz,tmp;

  const int nBuf=(Nsize*Nsize)*NQPFull;
  const int nsize = Nsize;
  const int ne = Ne;
  const int nQPFull = NQPFull;
  const int nMPTrans = NMPTrans;
  const int nTrans = NMPTrans * NQPOptTrans;
  const int nSlater = NSlater;
  int islater;

  int m,n,k,lda,ldb,ldc;
  double complex alpha, beta;
  char transa, transb;

  int bfidx;
  const double complex *sltE;
  const double complex *sltE_i,*sltE_j,*sltE_k;
  double complex *buffer;
  //double bufM[NQPFull*16*Nsize*Nsize];
  double complex *bufM;
  double complex *bufM_i, *bufM_i2, *bufM_ni;
  //double complex trM0=0.0,trM1=0.0,trM2=0.0,trM3=0.0,trM4=0.0,trM5=0.0,trM6=0.0;
  double complex trM[NProjBF];
  double complex pfM;
  int *flagTmp;
  int rki,rkj,rli,rlj;
  const int nSite=Nsite;
  const int nRange=Nrange;
  const int nSiteRange = nRange*nSite;
  int idx_ik,idx_jk,idx_il,idx_jl;
  const int *bfCnt0=eleProjBFCnt;
  const int *bfCnt1=eleProjBFCnt+4*Nsite*Nrange;
  const int *bfCnt0_n,*bfCnt0_m,*bfCnt1_n,*bfCnt1_m;
  int **posBF = PosBF;
  int xi,xj,xk,xl,xn,xm,xid,xih,xjd,xjh,xlh,xkh;
  int xidh,xihd,xjhd,xjdh,xkdh,xkhd,xldh,xlhd;
  int dki,dlj,nidx,midx,itmp;

  RequestWorkSpaceInt(2*nTrans*Nsize*Nsize);
  //RequestWorkSpaceDouble(NQPFull*NSlater);
  RequestWorkSpaceComplex(16*NQPFull*Nsize*Nsize);

  transOrbIdx = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbIdx[mpidx][msi][msj] */
  transOrbSgn = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbSgn[mpidx][msi][msj] */
  //buffer = GetWorkSpaceDouble(NQPFull*NSlater);
  bufM = GetWorkSpaceComplex(16*NQPFull*Nsize*Nsize);

    for(bfidx=0;bfidx<NProjBF;bfidx++) {
      srOptO[2*bfidx] = 0.0 + 0.0*I;
      srOptO[2*bfidx+1] = 0.0 + 0.0*I;
    }

  #pragma omp parallel for default(shared)        \
    private(mpidx,msi,msj,xqp,ri,tri,rj,trj,      \
            tOrbIdx,tOrbIdx_i,orbitalIdx_i)
  #pragma loop noalias
  for(mpidx=0;mpidx<nMPTrans;mpidx++) {
    xqp = QPTrans[mpidx];
    tOrbIdx = transOrbIdx + mpidx*nsize*nsize;
    for(msi=0;msi<nsize;msi++) {
      ri = eleIdx[msi];
      tri = xqp[ri];
      tOrbIdx_i = tOrbIdx + msi*nsize;
      orbitalIdx_i = OrbitalIdx[tri];
      for(msj=0;msj<nsize;msj++) {
        rj = eleIdx[msj];
        trj = xqp[rj];
        tOrbIdx_i[msj] = orbitalIdx_i[trj];
      }
    }
  }

  #pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,spidx,cs,cc,ss,           \
            xid,xih,xjd,xjh,xkh,xlh,rtmp,rtmp2,   \
            s,t,rk,rsk,rtl,rl,idx,idx2,           \
            sltE,sltE_i,sltE_j,sltE_k,            \
            rstmp,rttmp,t_ki,t_lj,                \
            tOrbIdx,invM,buf,msi,msj,             \
            tOrbIdx_i,invM_i,orbidx)              
    //reduction(-:trM0,trM1,trM2)
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    mpidx = qpidx / NSPGaussLeg;
    spidx = qpidx % NSPGaussLeg;
    xqp = QPTrans[mpidx];

    tOrbIdx = transOrbIdx + mpidx*nsize*nsize;
    sltE = SlaterElm + qpidx*Nsite2*Nsite2;
    invM = InvM + qpidx*Nsize*Nsize;
    buf = buffer + qpidx*Nsize*Nsize;
    pfM = PfM[qpidx];
    //flagTmp = etaFlag + qpidx*Nsite*Nsite;

    for(bfidx=0;bfidx<NProjBF;bfidx++) {
      trM[bfidx] = 0.0+0.0*I;
    }
    /* store bufM */
    /* Note that bufM is row-major and skew-symmetric. */
    /* bufM[msj][msi] = sltE[rsi][rsj] */
    #pragma loop noalias
    for(msi=0;msi<ne;msi++) {
      ri = eleIdx[msi];
      tri = xqp[ri];
      rsi = ri + (msi/Ne)*Nsite;
      bufM_i = bufM + msi*Nsize;
      sltE_i = sltE + rsi*Nsite2;
      invM_i = invM + msi*Nsize;
      #pragma loop norecurrence
      for(msj=ne;msj<nsize;msj++) {
        rj = eleIdx[msj];
        trj = xqp[rj];
        rsj = rj + (msj/Ne)*Nsite;
        if(etaFlag[tri][trj] == 0){
          bufM_i[msj] = 0.0;
        }else{
          bufM_i[msj] = Slater[ OrbitalIdx[tri][trj]]*OrbitalSgn[tri][trj];
        }
        trM[0] -= 2.0*invM_i[msj]*bufM_i[msj];
      }
    }

     //printf("eta1\n");
    /* Clear bufM */
    for(msi=0;msi<16*nsize*nsize;msi++) {
      bufM[msi]=0.0+0.0*I;
    }
    /* store bufM */
    /* Note that bufM is row-major and skew-symmetric. */
    /* bufM[msi][msj] = sltE[rsi][rsj] */
    #pragma loop noalias
    for(msi=0;msi<ne;msi++) {
      ri = eleIdx[msi];
      tri = xqp[ri];
      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
      bufM_i = bufM + msi*Nsize;
      invM_i = invM + msi*Nsize;
      xid = n0[tri]*n1[tri];
      #pragma loop norecurrence
      for(msj=ne;msj<nsize;msj++) {
        rj = eleIdx[msj];
        trj = xqp[rj];
        rsj = eleIdx[msj] + (msj/Ne)*Nsite;
        sltE_j = sltE + rsj*Nsite2;
        xjd = n0[trj]*n1[trj];

        for(xn=0;xn<4;xn++) {
          for(xm=0;xm<4;xm++){
            if(xm==0 && xn == 0) continue;
            bfCnt0_n=bfCnt0+xn*nSiteRange;
            bfCnt1_n=bfCnt1+xn*nSiteRange;
            bfCnt0_m=bfCnt0+xm*nSiteRange;
            bfCnt1_m=bfCnt1+xm*nSiteRange;
            
            bufM_ni = bufM_i + (4*xn+xm)*Nsize*Nsize;
            for(xk=0;xk<nRange;xk++) {
              rki=posBF[tri][xk];
              rkj=posBF[trj][xk];
              idx_ik=tri*nRange+xk;
              idx_jk=trj*nRange+xk;
            
              dki = RangeIdx[tri][rki];
              itmp = 4*dki+xn;
              nidx = itmp-3-dki;
              //if(itmp <= 3 || itmp%4==0){continue;}
              if((itmp%4)==0){nidx=-1;}
              if(itmp==0){nidx=0;}
              if(nidx<0){continue;}
            
              for(xl=0;xl<nRange;xl++){
                rli=posBF[tri][xl];
                rlj=posBF[trj][xl];
                idx_il=tri*nRange+xl;
                idx_jl=trj*nRange+xl;
            
                dlj = RangeIdx[trj][rlj];
                itmp = 4*dlj+xm;
                midx = itmp-3-dlj;
                if((itmp%4)==0){midx=-1;}
                if(itmp==0){midx=0;}
                if(midx<0){continue;}
                bfidx=BFSubIdx[nidx][midx];
                          
                tmp = -bfCnt0_n[idx_ik]*bfCnt1_m[idx_jl]*Slater[ OrbitalIdx[rki][rlj]]*OrbitalSgn[rki][rlj]
                      -bfCnt1_n[idx_jk]*bfCnt0_m[idx_il]*Slater[ OrbitalIdx[rli][rkj]]*OrbitalSgn[rli][rkj];
                trM[bfidx] -= invM_i[msj]*tmp;
              }
            }
          }
        }

      }
    }
    //TODO: Check
    for(bfidx=0;bfidx<NProjBF;bfidx++) {
      srOptO[2*bfidx]+=0.5*trM[bfidx]*pfM*invIP;
      srOptO[2*bfidx+1]+=0.5*trM[bfidx]*pfM*invIP*I;
    }
  }//QP

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();

  return;
}


void MakeSlaterElmBF_fcmp(const int *eleNum, const int *eleProjBFCnt) {
  const int *n0=eleNum;
  const int *n1=eleNum+Nsite;
  const int *bfCnt0=eleProjBFCnt;
  const int *bfCnt1=eleProjBFCnt+4*Nsite*Nrange;
  const int *bfCnt0_n,*bfCnt0_m,*bfCnt1_n,*bfCnt1_m;
  int **posBF = PosBF;
  int xi,xj,xk,xl,xn,xm,xid,xih,xjd,xjh,xlh,xkh;
  int xidh,xihd,xjhd,xjdh,xkdh,xkhd,xldh,xlhd;
  int s,idx, rk,rtmp;
  int t,idx2,rl,rtmp2;
  int rki,rkj,rli,rlj;
  int icount, jcount;
  //double eta;
  int ri,ori,tri,sgni,rsi0,rsi1;
  int rj,orj,trj,sgnj,rsj0,rsj1;
  int qpidx,mpidx,spidx,optidx;
  double complex cs,cc,ss;
  double complex slt_ij,slt_ji;
  double complex tmp_ij,tmp_ji;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex *sltE,*sltE_i0,*sltE_i1;
  double complex sum_ij,sum_ji;
  int *flagTmp;
  const int nSite=Nsite;
  const int nRange=Nrange;
  const int nSiteRange = nRange*nSite;
  int idx_ik,idx_jk,idx_il,idx_jl;
  int bfidx;
  int dki,dlj,nidx,midx,itmp;


#pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,spidx,xqp,cs,cc,ss,sltE,  \
            icount,jcount,    \
            ri,tri,rsi0,rsi1,sltE_i0,sltE_i1,     \
            rj,trj,rsj0,rsj1,slt_ij,slt_ji)
#pragma loop noalias
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    mpidx = qpidx / NSPGaussLeg;
    spidx = qpidx % NSPGaussLeg;

    xqp = QPTrans[mpidx];
    cs = SPGLCosSin[spidx];
    cc = SPGLCosCos[spidx];
    ss = SPGLSinSin[spidx];

    sltE = SlaterElmBF + qpidx*Nsite2*Nsite2;

    for(ri=0;ri<Nsite;ri++) {
      tri = xqp[ri];
      rsi0 = ri;
      rsi1 = ri+Nsite;
      sltE_i0 = sltE + rsi0*Nsite2;
      sltE_i1 = sltE + rsi1*Nsite2;

      for(rj=0;rj<Nsite;rj++) {
        trj = xqp[rj];
        rsj0 = rj;
        rsj1 = rj+Nsite;

        /* backflow correlation factor */
        SubSlaterElmBF_fcmp(tri,trj,&slt_ij,&icount,&slt_ji,&jcount,eleProjBFCnt);
        //printf("icount=%d, slt_ij=%.2e\n",icount,slt_ij);

        if(icount == 0){eta[tri][trj] = 1.0;  etaFlag[tri][trj]=0;}
        else{eta[tri][trj] = creal(ProjBF[0]); etaFlag[tri][trj]=1;}

        if(jcount == 0){eta[trj][tri] = 1.0; etaFlag[trj][tri]=0;}
        else{eta[trj][tri] = creal(ProjBF[0]); etaFlag[trj][tri]=1;}

        sltE_i0[rsj0] = -(slt_ij - slt_ji)*cs;
        sltE_i0[rsj1] = slt_ij*cc + slt_ji*ss;
        sltE_i1[rsj0] = -slt_ij*ss - slt_ji*cc;
        sltE_i1[rsj1] = (slt_ij - slt_ji)*cs;
      }
    }
  }

  return;
}

void SubSlaterElmBF_fcmp(const int tri, const int trj, double complex *slt_ij, int *ijcount, double complex* slt_ji, int *jicount, const int *eleProjBFCnt){
  int xn,xm,xk,xl;
  int rki,rkj,rli,rlj;
  int idx_ik,idx_jk,idx_il,idx_jl;
  int bfidx;
  int dki,dlj,nidx,midx,xtmp;
  const int nSite=Nsite;
  const int nRange=Nrange;
  const int nSiteRange = nRange*nSite;
  const int *bfCnt0=eleProjBFCnt;
  const int *bfCnt1=eleProjBFCnt+4*Nsite*Nrange;
  const int *bfCnt0_n,*bfCnt0_m,*bfCnt1_n,*bfCnt1_m;
  double eta;
  int **posBF = PosBF;

  *slt_ij = 0.0;
  *slt_ji = 0.0;
  *ijcount = 0;
  *jicount = 0;
  //#pragma omp parallel for reduction(+:slt_ij,slt_ji)
  for(xn=0;xn<4;xn++){
    bfCnt0_n=bfCnt0+xn*nSiteRange;
    bfCnt1_n=bfCnt1+xn*nSiteRange;
    for(xm=0;xm<4;xm++){
      if(xm==0 && xn == 0) continue;
      bfCnt0_m=bfCnt0+xm*nSiteRange;
      bfCnt1_m=bfCnt1+xm*nSiteRange;

      for(xk=0;xk<nRange;xk++) {
        rki=posBF[tri][xk];
        rkj=posBF[trj][xk];
        idx_ik=tri*nRange+xk;
        idx_jk=trj*nRange+xk;

        dki = RangeIdx[tri][rki];
        xtmp = 4*dki+xn;
        nidx = xtmp-3-dki;
        if(xtmp%4==0){nidx=-1;}
        if(xtmp==0){nidx=0;}
        if(nidx<0){continue;}

        *ijcount += bfCnt0[nSiteRange+idx_ik]+bfCnt0[nSiteRange+idx_jk];
        *jicount += bfCnt0[nSiteRange+idx_ik]+bfCnt0[nSiteRange+idx_jk];

        for(xl=0;xl<nRange;xl++){
          rli=posBF[tri][xl];
          rlj=posBF[trj][xl];
          idx_il=tri*nRange+xl;
          idx_jl=trj*nRange+xl;

          dlj = RangeIdx[trj][rlj];
          xtmp = 4*dlj+xm;
          midx = xtmp-3-dlj;
          if(xtmp%4==0){midx=-1;}
          if(xtmp==0){midx=0;}
          if(midx<0){continue;}
          bfidx=BFSubIdx[nidx][midx];

          //printf("ProjBF[%d]=%.2e\n",bfidx,ProjBF[bfidx]);
          //printf("OrbitalSgn[%d][%d]=%d\n",rki,rlj,OrbitalSgn[rki][rlj]);
          //printf("Slater[%d]=%.2e\n",OrbitalIdx[rki][rlj], Slater[ OrbitalIdx[rki][rlj]]);
          //printf("bfCnt0_n[%d]=%d\n",idx_ik,bfCnt0_n[idx_ik]);
          *slt_ij += -ProjBF[bfidx]*bfCnt0_n[idx_ik]*bfCnt1_m[idx_jl]*Slater[ OrbitalIdx[rki][rlj]]*OrbitalSgn[rki][rlj];
          *slt_ji += -ProjBF[bfidx]*bfCnt1_n[idx_ik]*bfCnt0_m[idx_jl]*Slater[ OrbitalIdx[rlj][rki]]*OrbitalSgn[rlj][rki];
        }
      }
    }
  }

  if(*ijcount == 0){eta = 1.0;}
  else{eta = creal(ProjBF[0]);} //TODO: Check
  *slt_ij += eta*Slater[ OrbitalIdx[tri][trj] ]*OrbitalSgn[tri][trj];

  if(*jicount == 0){eta = 1.0;}
  else{eta = creal(ProjBF[0]);}
  *slt_ji += eta*Slater[ OrbitalIdx[trj][tri] ]*OrbitalSgn[trj][tri];

  return;

}


void SubSlaterElmBF_real(const int tri, const int trj, double *slt_ij, int *ijcount, double* slt_ji, int *jicount, const int *eleProjBFCnt){
    int xn,xm,xk,xl;
    int rki,rkj,rli,rlj;
    int idx_ik,idx_jk,idx_il,idx_jl;
    int bfidx;
    int dki,dlj,nidx,midx,xtmp;
    const int nSite=Nsite;
    const int nRange=Nrange;
    const int nSiteRange = nRange*nSite;
    const int *bfCnt0=eleProjBFCnt;
    const int *bfCnt1=eleProjBFCnt+4*Nsite*Nrange;
    const int *bfCnt0_n,*bfCnt0_m,*bfCnt1_n,*bfCnt1_m;
    double eta;
    int **posBF = PosBF;

    *slt_ij = 0.0;
    *slt_ji = 0.0;
    *ijcount = 0;
    *jicount = 0;
    //#pragma omp parallel for reduction(+:slt_ij,slt_ji)
    for(xn=0;xn<4;xn++){
        bfCnt0_n=bfCnt0+xn*nSiteRange;
        bfCnt1_n=bfCnt1+xn*nSiteRange;
        for(xm=0;xm<4;xm++){
            if(xm==0 && xn == 0) continue;
            bfCnt0_m=bfCnt0+xm*nSiteRange;
            bfCnt1_m=bfCnt1+xm*nSiteRange;

            for(xk=0;xk<nRange;xk++) {
                rki=posBF[tri][xk];
                rkj=posBF[trj][xk];
                idx_ik=tri*nRange+xk;
                idx_jk=trj*nRange+xk;

                dki = RangeIdx[tri][rki];
                xtmp = 4*dki+xn;
                nidx = xtmp-3-dki;
                if(xtmp%4==0){nidx=-1;}
                if(xtmp==0){nidx=0;}
                if(nidx<0){continue;}

                *ijcount += bfCnt0[nSiteRange+idx_ik]+bfCnt0[nSiteRange+idx_jk];
                *jicount += bfCnt0[nSiteRange+idx_ik]+bfCnt0[nSiteRange+idx_jk];

                for(xl=0;xl<nRange;xl++){
                    rli=posBF[tri][xl];
                    rlj=posBF[trj][xl];
                    idx_il=tri*nRange+xl;
                    idx_jl=trj*nRange+xl;

                    dlj = RangeIdx[trj][rlj];
                    xtmp = 4*dlj+xm;
                    midx = xtmp-3-dlj;
                    if(xtmp%4==0){midx=-1;}
                    if(xtmp==0){midx=0;}
                    if(midx<0){continue;}
                    bfidx=BFSubIdx[nidx][midx];

                    //printf("ProjBF[%d]=%.2e\n",bfidx,ProjBF[bfidx]);
                    //printf("OrbitalSgn[%d][%d]=%d\n",rki,rlj,OrbitalSgn[rki][rlj]);
                    //printf("Slater[%d]=%.2e\n",OrbitalIdx[rki][rlj], Slater[ OrbitalIdx[rki][rlj]]);
                    //printf("bfCnt0_n[%d]=%d\n",idx_ik,bfCnt0_n[idx_ik]);
                    *slt_ij += -ProjBF[bfidx]*bfCnt0_n[idx_ik]*bfCnt1_m[idx_jl]*creal(Slater[ OrbitalIdx[rki][rlj]])*OrbitalSgn[rki][rlj];
                    *slt_ji += -ProjBF[bfidx]*bfCnt1_n[idx_ik]*bfCnt0_m[idx_jl]*creal(Slater[ OrbitalIdx[rlj][rki]])*OrbitalSgn[rlj][rki];
                }
            }
        }
    }

    if(*ijcount == 0){eta = 1.0;}
    else{eta = creal(ProjBF[0]);} //TODO: Check
    *slt_ij += eta*Slater[ OrbitalIdx[tri][trj] ]*OrbitalSgn[tri][trj];

    if(*jicount == 0){eta = 1.0;}
    else{eta = creal(ProjBF[0]);}
    *slt_ji += eta*Slater[ OrbitalIdx[trj][tri] ]*OrbitalSgn[trj][tri];

    return;

}


void UpdateSlaterElmBF_fcmp(const int ma, const int ra, const int rb, const int u,
                       const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double complex*sltElmTmp){
    const int *n0=eleNum;
    const int *n1=eleNum+Nsite;
    const int rua=ra+Nsite*u, rub=rb+Nsite*u;
    int **posBF = PosBF;
    const int mua = ma + Ne*u;
    int trua, trub, trsi, trsj;
    int s,idx, rk,rtmp;
    int t,idx2,rl,rtmp2;
    double eta;
    int ri,ori,tri,sgni,rsi0,rsi1;
    int rj,orj,trj,sgnj,rsj0,rsj1;
    int qpidx,mpidx,spidx,optidx;
    double complex cs,cc,ss;
    double complex slt_ij,slt_ji;
    double complex tmp_ij,tmp_ji;
    int *xqp, *xqpInv, *xqpSgn, *xqpOpt, *xqpOptSgn;
    double complex*sltE,*sltE_i0,*sltE_i1;
    int rsz[Nsite2];
    int rhop0[Nsite2],rhop1[Nsite2];
    double complex*pTrans_i;
    int zidx,zidx2,rsi,rsj,msi,msj,itmp,hop,icount,mi,mj;
    int jcount[Nsite],jcount1[Nsite];
    int jcountTmp[Nsite][Nsite],jcount1Tmp[Nsite][Nsite];
    int ri0,ri1,mi0,mi1,mj0,mj1,flag,hidx,itmp0=0,itmp1=0;
    int ijcount,jicount;
    int rki,rkj,rli,rlj;
    const int nSite=Nsite;
    const int nRange=Nrange;
    const int nSiteRange = nRange*nSite;
    int idx_ik,idx_jk,idx_il,idx_jl;
    int bfidx;
    int xn,xm,xk,xl;
    int dki,dlj,nidx,midx,xtmp,xtmp2;
    const int *bfCnt0=eleProjBFCnt;
    const int *bfCnt1=eleProjBFCnt+4*Nsite*Nrange;
    const int *bfCnt0_n,*bfCnt0_m,*bfCnt1_n,*bfCnt1_m;
    double complex sltElm[Nsite*Nsite],sltElm2[Nsite*Nsite];

    for(qpidx=0;qpidx<NQPFull;qpidx++) {
        StartTimer(91);
        itmp0=0;
        itmp1=0;
        mpidx = qpidx / NSPGaussLeg;
        spidx = qpidx % NSPGaussLeg;

        xqp = QPTrans[mpidx];
        xqpInv = QPTransInv[mpidx];

        sltE = sltElmTmp + qpidx*Nsite2*Nsite2;

        icount=0;
        trua = xqpInv[ra];
        rsz[icount] = trua;
        icount++;
        for(idx=0;idx<Nrange;idx++) {
            rtmp=posBF[trua][idx];
            flag=0;
            for(zidx=0;zidx<icount;zidx++){
                if(rsz[zidx] == rtmp) flag=1;
            }
            if(flag==1) continue;
            rsz[icount] = rtmp;
            icount++;
        }
        trub = xqpInv[rb];
        rsz[icount] = trub;
        icount++;
        for(idx=0;idx<Nrange;idx++) {
            rtmp=posBF[trub][idx];
            flag=0;
            for(zidx=0;zidx<icount;zidx++){
                if(rsz[zidx] == rtmp) flag=1;
            }
            if(flag==1) continue;
            rsz[icount] = rtmp;
            icount++;
        }
        hop=0;
        msa[qpidx*Nsite+hop]= ma + Ne*u;
        hop++;

        itmp=0;
        for(zidx=0;zidx<icount;zidx++){
            jcount[zidx]=0;
            jcount1[zidx]=0;
            for(rj=0;rj<Nsite;rj++){
                jcountTmp[rj][zidx] = 0;
                jcount1Tmp[rj][zidx] = 0;
            }
        }

#pragma omp parallel for default(shared)
        for(zidx=0;zidx<Nsite*Nsite;zidx++){
            sltElm[zidx] = 0.0;
            sltElm2[zidx] = 0.0;
        }

        StopTimer(91);
        StartTimer(92);
#pragma omp parallel for default(shared) \
      private(zidx,  \
          ri,tri,rsi0,rsi1,rj,trj,rsj0,rsj1, \
           ijcount,jicount,slt_ij,slt_ji)
        for(zidx=0;zidx<icount*Nsite;zidx++){
            ri  = rsz[zidx/Nsite];
            tri = xqp[ri];
            rsi0 = ri;
            rsi1 = ri+Nsite;

            rj = zidx%Nsite;
            trj = xqp[rj];
            rsj0 = rj;
            rsj1 = rj+Nsite;

            SubSlaterElmBF_fcmp(tri,trj,&slt_ij,&ijcount,&slt_ji,&jicount,eleProjBFCnt);

            sltElm[tri*Nsite+trj] = slt_ij;
            sltElm2[tri*Nsite+trj] = slt_ji;
        }

        StopTimer(92);
        StartTimer(93);
        for(zidx=0;zidx<icount*Nsite;zidx++){
            ri  = rsz[zidx/Nsite];
            tri = xqp[ri];
            rsi0 = ri;
            rsi1 = ri+Nsite;

            rj = zidx%Nsite;
            trj = xqp[rj];
            rsj0 = rj;
            rsj1 = rj+Nsite;

            sltE_i0 = sltE + rsi0*Nsite2;
            sltE_i1 = sltE + rsi1*Nsite2;

            slt_ij = sltElm[tri*Nsite+trj];
            slt_ji = sltElm2[tri*Nsite+trj];

            if(sltE_i0[rsj1] != slt_ij){
                sltE_i0[rsj1] = slt_ij;
                jcountTmp[rj][zidx/Nsite]=1;
            }
            if(sltE_i1[rsj0] != -slt_ji){
                sltE_i1[rsj0] = -slt_ji;
                jcount1Tmp[rj][zidx/Nsite]=1;
            }
        }

        for(zidx=0;zidx<icount*Nsite;zidx++){
            jcount[zidx/Nsite] += jcountTmp[zidx%Nsite][zidx/Nsite];
            jcount1[zidx/Nsite]+= jcount1Tmp[zidx%Nsite][zidx/Nsite];
        }

        //#pragma loop noalias
        for(zidx=0;zidx<icount;zidx++){
            ri  = rsz[zidx];
            tri = xqp[ri];
            rsi0 = ri;
            rsi1 = ri+Nsite;
            mi0 = eleCfg[rsi0];
            mi1 = eleCfg[rsi1];

            //if((jcount == Nsite)){
            if((jcount[zidx] >= 1)){
                rhop0[itmp0] = rsi0;
                itmp0++;
            }
            //if((jcount == Nsite) && (mi0 != -1) && (mua != mi0)){
            if((jcount[zidx] >= Ne) && (mi0 != -1) && (mua != mi0)){
                msa[qpidx*Nsite+hop]=mi0;
                hop++;
            }
            //if((jcount1 == Nsite) ){
            if((jcount1[zidx] >= 1) ){
                rhop1[itmp1] = rsi1;
                itmp1++;
            }
            //if((jcount1 == Nsite) && (mi1 != -1) && (mua != mi1+Ne)){
            if((jcount1[zidx] >= Ne) && (mi1 != -1) && (mua != mi1+Ne)){
                msa[qpidx*Nsite+hop]=mi1+Ne;
                hop++;
            }
        }
        for(hidx=0;hidx<itmp0;hidx++){
            rsi0= rhop0[hidx];
            for(rj=0;rj<Nsite;rj++) {
                rsj0 = rj;
                rsj1 = rj+Nsite;
                sltE[rsj1*Nsite2+rsi0]=-sltE[rsi0*Nsite2+rsj1];
            }
        }
        for(hidx=0;hidx<itmp1;hidx++){
            rsi1= rhop1[hidx];
            for(rj=0;rj<Nsite;rj++) {
                rsj0 = rj;
                rsj1 = rj+Nsite;
                sltE[rsj0*Nsite2+rsi1]=-sltE[rsi1*Nsite2+rsj0];
            }
        }

        hopNum[qpidx] = hop;

        StopTimer(93);
    }

    return ;
}

void UpdateSlaterElmBFGrn(const int ma, const int ra, const int rb, const int u,
                          const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double complex* sltElmTmp){
    const int *n0=eleNum;
    const int *n1=eleNum+Nsite;
    const int rua=ra+Nsite*u, rub=rb+Nsite*u;
    int **posBF = PosBF;
    //int rua=ra+Nsite*u, rub=rb+Nsite*u;
    const int mua = ma + Ne*u;
    int trua, trub, trsi, trsj;
    int s,idx, rk,rtmp;
    int t,idx2,rl,rtmp2;
    double eta;
    int ri,ori,tri,sgni,rsi0,rsi1;
    int rj,orj,trj,sgnj,rsj0,rsj1;
    int qpidx,mpidx,spidx,optidx;
    double complex cs,cc,ss;
    double complex slt_ij,slt_ji;
    double complex tmp_ij,tmp_ji;
    int *xqp, *xqpInv, *xqpSgn, *xqpOpt, *xqpOptSgn;
    double complex *sltE,*sltE_i0,*sltE_i1;
    //double complex pTrans[Nsite2*Nsite2];
    //double complex *pTrans_i;
    int rsz[Nsite2];
    int rhop0[Nsite2],rhop1[Nsite2];
    double complex* pTrans_i;
    int zidx,zidx2,rsi,rsj,msi,msj,itmp,hop,icount,mi,mj;
    int jcount[Nsite],jcount1[Nsite];
    int ri0,ri1,mi0,mi1,mj0,mj1,flag,hidx,itmp0=0,itmp1=0;
    int ijcount,jicount;
    int rki,rkj,rli,rlj;
    const int nSite=Nsite;
    const int nRange=Nrange;
    const int nSiteRange = nRange*nSite;

    for(qpidx=0;qpidx<NQPFull;qpidx++) {
        itmp0=0;
        itmp1=0;
        mpidx = qpidx / NSPGaussLeg;
        spidx = qpidx % NSPGaussLeg;

        xqp = QPTrans[mpidx];
        xqpInv = QPTransInv[mpidx];
        cs = SPGLCosSin[spidx];
        cc = SPGLCosCos[spidx];
        ss = SPGLSinSin[spidx];

        sltE = sltElmTmp + qpidx*Nsite2*Nsite2;

        icount=0;
        trua = xqpInv[ra];
        rsz[icount] = trua;
        icount++;
        for(idx=0;idx<Nrange;idx++) {
            rtmp=posBF[trua][idx];
            flag=0;
            for(zidx=0;zidx<icount;zidx++){
                if(rsz[zidx] == rtmp) flag=1;
            }
            if(flag==1) continue;
            rsz[icount] = rtmp;
            icount++;
        }
        trub = xqpInv[rb];
        rsz[icount] = trub;
        icount++;
        for(idx=0;idx<Nrange;idx++) {
            rtmp=posBF[trub][idx];
            flag=0;
            for(zidx=0;zidx<icount;zidx++){
                if(rsz[zidx] == rtmp) flag=1;
            }
            if(flag==1) continue;
            rsz[icount] = rtmp;
            icount++;
        }
        hop=0;
        msa[qpidx*Nsite+hop]= ma + Ne*u;
        hop++;

        itmp=0;
        for(zidx=0;zidx<icount;zidx++){
            jcount[zidx]=0;
            jcount1[zidx]=0;
        }
        //#pragma omp parallel for reduction(+:jcount,jcount1,itmp,itmp0,itmp1,hop)
        //private(zidx,mi0,mi1,jcount,jcount1,itmp,  \
    //        ri,tri,rsi0,rsi1,sltE_i0,sltE_i1,     \
    //        ijcount,jicount,qpidx,hop,itmp1,itmp0,\
    //        rj,trj,rsj0,rsj1,slt_ij,slt_ji) \

        //#pragma omp parallel for default(shared)        \
    //private(zidx,mi0,mi1,  \
    //        ri,tri,rsi0,rsi1,sltE_i0,sltE_i1,     \
    //        ijcount,jicount,qpidx,\
    //        rj,trj,rsj0,rsj1,slt_ij,slt_ji)
        //reduction(+:jcount,jcount1,itmp,itmp0,itmp1,hop)
        //#pragma loop noalias
        for(zidx=0;zidx<icount*Nsite;zidx++) {
            ri  = rsz[zidx/Nsite];
            tri = xqp[ri];
            rsi0 = ri;
            rsi1 = ri+Nsite;
            sltE_i0 = sltE + rsi0*Nsite2;
            sltE_i1 = sltE + rsi1*Nsite2;

            mi0 = eleCfg[rsi0];
            mi1 = eleCfg[rsi1];

            rj = zidx%Nsite;
            trj = xqp[rj];
            rsj0 = rj;
            rsj1 = rj+Nsite;

            SubSlaterElmBF_fcmp(tri,trj,&slt_ij,&ijcount,&slt_ji,&jicount,eleProjBFCnt);

            if(sltE_i0[rsj1] != slt_ij){
                sltE_i0[rsj1] = slt_ij;
                jcount[zidx/Nsite]++;
            }
            if(sltE_i1[rsj0] != -slt_ji){
                sltE_i1[rsj0] = -slt_ji;
                jcount1[zidx/Nsite]++;
            }
            itmp++;
        }
        for(zidx=0;zidx<icount;zidx++){
            ri  = rsz[zidx];
            tri = xqp[ri];
            rsi0 = ri;
            rsi1 = ri+Nsite;
            mi0 = eleCfg[rsi0];
            mi1 = eleCfg[rsi1];

            //if((jcount == Nsite)){
            if((jcount[zidx] >= 1)){
                rhop0[itmp0] = rsi0;
                itmp0++;
            }
            //if((jcount == Nsite) && (mi0 != -1) && (mua != mi0)){
            if((jcount[zidx] >= Ne) && (mi0 != -1) && (mua != mi0)){
                msa[qpidx*Nsite+hop]=mi0;
                hop++;
            }
            //if((jcount1 == Nsite) ){
            if((jcount1[zidx] >= 1) ){
                rhop1[itmp1] = rsi1;
                itmp1++;
            }
            //if((jcount1 == Nsite) && (mi1 != -1) && (mua != mi1+Ne)){
            if((jcount1[zidx] >= Ne) && (mi1 != -1) && (mua != mi1+Ne)){
                msa[qpidx*Nsite+hop]=mi1+Ne;
                hop++;
            }
        }
        for(hidx=0;hidx<itmp0;hidx++){
            rsi0= rhop0[hidx];
            for(rj=0;rj<Nsite;rj++) {
                rsj0 = rj;
                rsj1 = rj+Nsite;
                sltE[rsj1*Nsite2+rsi0]=-sltE[rsi0*Nsite2+rsj1];
            }
        }
        for(hidx=0;hidx<itmp1;hidx++){
            rsi1= rhop1[hidx];
            for(rj=0;rj<Nsite;rj++) {
                rsj0 = rj;
                rsj1 = rj+Nsite;
                sltE[rsj0*Nsite2+rsi1]=-sltE[rsi1*Nsite2+rsj0];
            }
        }

        hopNum[qpidx] = hop;

    }

    return ;
}


void UpdateSlaterElmBFGrn_real(const int ma, const int ra, const int rb, const int u, 
                     const int *eleCfg, const int *eleNum, const int *eleProjBFCnt, int *msa, int *hopNum, double *sltElmTmp){
  const int *n0=eleNum;
  const int *n1=eleNum+Nsite;
  const int rua=ra+Nsite*u, rub=rb+Nsite*u;
  int **posBF = PosBF;
  //int rua=ra+Nsite*u, rub=rb+Nsite*u;
  const int mua = ma + Ne*u;
  int trua, trub, trsi, trsj;
  int s,idx, rk,rtmp;
  int t,idx2,rl,rtmp2;
  double eta;
  int ri,ori,tri,sgni,rsi0,rsi1;
  int rj,orj,trj,sgnj,rsj0,rsj1;
  int qpidx,mpidx,spidx,optidx;
  double cs,cc,ss;
  double slt_ij,slt_ji;
  double tmp_ij,tmp_ji;
  int *xqp, *xqpInv, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double *sltE,*sltE_i0,*sltE_i1;
  //double complex pTrans[Nsite2*Nsite2];
  //double complex *pTrans_i;
  int rsz[Nsite2];
  int rhop0[Nsite2],rhop1[Nsite2];
  double *pTrans_i;
  int zidx,zidx2,rsi,rsj,msi,msj,itmp,hop,icount,mi,mj;
  int jcount[Nsite],jcount1[Nsite];
  int ri0,ri1,mi0,mi1,mj0,mj1,flag,hidx,itmp0=0,itmp1=0;
  int ijcount,jicount;
  int rki,rkj,rli,rlj;
  const int nSite=Nsite;
  const int nRange=Nrange;
  const int nSiteRange = nRange*nSite;

  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    itmp0=0;
    itmp1=0;
    mpidx = qpidx / NSPGaussLeg;
    spidx = qpidx % NSPGaussLeg;

    xqp = QPTrans[mpidx];
    xqpInv = QPTransInv[mpidx];
    cs = SPGLCosSin[spidx];
    cc = SPGLCosCos[spidx];
    ss = SPGLSinSin[spidx];
    
    sltE = sltElmTmp + qpidx*Nsite2*Nsite2;

    icount=0;
    trua = xqpInv[ra];
    rsz[icount] = trua;
    icount++;
    for(idx=0;idx<Nrange;idx++) {
      rtmp=posBF[trua][idx];
      flag=0;
      for(zidx=0;zidx<icount;zidx++){
        if(rsz[zidx] == rtmp) flag=1;
      }
      if(flag==1) continue;
      rsz[icount] = rtmp;
      icount++;
    }
    trub = xqpInv[rb];
    rsz[icount] = trub;
    icount++;
    for(idx=0;idx<Nrange;idx++) {
      rtmp=posBF[trub][idx];
      flag=0;
      for(zidx=0;zidx<icount;zidx++){
        if(rsz[zidx] == rtmp) flag=1;
      }
      if(flag==1) continue;
      rsz[icount] = rtmp;
      icount++;
    }
    hop=0;
    msa[qpidx*Nsite+hop]= ma + Ne*u;
    hop++;

    itmp=0;
    for(zidx=0;zidx<icount;zidx++){
      jcount[zidx]=0;
      jcount1[zidx]=0;
    }
    //#pragma omp parallel for reduction(+:jcount,jcount1,itmp,itmp0,itmp1,hop)
    //private(zidx,mi0,mi1,jcount,jcount1,itmp,  \
    //        ri,tri,rsi0,rsi1,sltE_i0,sltE_i1,     \
    //        ijcount,jicount,qpidx,hop,itmp1,itmp0,\
    //        rj,trj,rsj0,rsj1,slt_ij,slt_ji) \
    
    //#pragma omp parallel for default(shared)        \
    //private(zidx,mi0,mi1,  \
    //        ri,tri,rsi0,rsi1,sltE_i0,sltE_i1,     \
    //        ijcount,jicount,qpidx,\
    //        rj,trj,rsj0,rsj1,slt_ij,slt_ji)
    //reduction(+:jcount,jcount1,itmp,itmp0,itmp1,hop)
    //#pragma loop noalias
    for(zidx=0;zidx<icount*Nsite;zidx++) {
      ri  = rsz[zidx/Nsite];
      tri = xqp[ri];
      rsi0 = ri;
      rsi1 = ri+Nsite;
      sltE_i0 = sltE + rsi0*Nsite2;
      sltE_i1 = sltE + rsi1*Nsite2;
      
      mi0 = eleCfg[rsi0];
      mi1 = eleCfg[rsi1];
      
      rj = zidx%Nsite;
      trj = xqp[rj];
      rsj0 = rj;
      rsj1 = rj+Nsite;
        
      SubSlaterElmBF_real(tri,trj,&slt_ij,&ijcount,&slt_ji,&jicount,eleProjBFCnt);

      if(sltE_i0[rsj1] != slt_ij){
        sltE_i0[rsj1] = slt_ij;
        jcount[zidx/Nsite]++;
      }
      if(sltE_i1[rsj0] != -slt_ji){
        sltE_i1[rsj0] = -slt_ji;
        jcount1[zidx/Nsite]++;
      }
      itmp++;
    } 
    for(zidx=0;zidx<icount;zidx++){
      ri  = rsz[zidx];
      tri = xqp[ri];
      rsi0 = ri;
      rsi1 = ri+Nsite;
      mi0 = eleCfg[rsi0];
      mi1 = eleCfg[rsi1];
      
      //if((jcount == Nsite)){
        if((jcount[zidx] >= 1)){
          rhop0[itmp0] = rsi0;
          itmp0++;
        }
      //if((jcount == Nsite) && (mi0 != -1) && (mua != mi0)){
        if((jcount[zidx] >= Ne) && (mi0 != -1) && (mua != mi0)){
          msa[qpidx*Nsite+hop]=mi0;
          hop++;
        }
      //if((jcount1 == Nsite) ){
        if((jcount1[zidx] >= 1) ){
          rhop1[itmp1] = rsi1;
          itmp1++;
        }
      //if((jcount1 == Nsite) && (mi1 != -1) && (mua != mi1+Ne)){
        if((jcount1[zidx] >= Ne) && (mi1 != -1) && (mua != mi1+Ne)){
          msa[qpidx*Nsite+hop]=mi1+Ne;
          hop++;
        }
    }
    for(hidx=0;hidx<itmp0;hidx++){
      rsi0= rhop0[hidx];
      for(rj=0;rj<Nsite;rj++) {
        rsj0 = rj;
        rsj1 = rj+Nsite;
        sltE[rsj1*Nsite2+rsi0]=-sltE[rsi0*Nsite2+rsj1];
      }
    }
    for(hidx=0;hidx<itmp1;hidx++){
      rsi1= rhop1[hidx];
      for(rj=0;rj<Nsite;rj++) {
        rsj0 = rj;
        rsj1 = rj+Nsite;
        sltE[rsj0*Nsite2+rsi1]=-sltE[rsi1*Nsite2+rsj0];
      }
    }

    hopNum[qpidx] = hop;

  }

  return ;
}


void StoreSlaterElmBF_fcmp(complex double *bufM){
    int qpidx,rsi,rsj;
    const double complex* sltE;
    const double complex* sltE_i;
    double complex*bufM_i, *bufM_i2;

    /* store SlaterElmBF before hopping */
    for(qpidx=0;qpidx<NQPFull;qpidx++) {
        sltE = SlaterElmBF + qpidx*Nsite2*Nsite2;
        for(rsi=0;rsi<Nsite2;rsi++) {
            bufM_i = bufM + rsi*Nsite2 + qpidx*Nsite2*Nsite2;
            sltE_i = sltE + rsi*Nsite2;
            //#pragma loop norecurrence
            for(rsj=0;rsj<Nsite2;rsj++) {
                bufM_i[rsj] = sltE_i[rsj];
                //printf("bufM[%d]=%.5e\n",msi*Nsize+msj,cabs(bufM[msi*Nsize+msj]));
            }
        }
    }

    return;
}

void StoreSlaterElmBF_real(double *bufM){
    int qpidx,rsi,rsj;
    const double * sltE;
    const double * sltE_i;
    double *bufM_i, *bufM_i2;

    /* store SlaterElmBF before hopping */
    for(qpidx=0;qpidx<NQPFull;qpidx++) {
        sltE = SlaterElmBF_real + qpidx*Nsite2*Nsite2;
        for(rsi=0;rsi<Nsite2;rsi++) {
            bufM_i = bufM + rsi*Nsite2 + qpidx*Nsite2*Nsite2;
            sltE_i = sltE + rsi*Nsite2;
            //#pragma loop norecurrence
            for(rsj=0;rsj<Nsite2;rsj++) {
                bufM_i[rsj] = sltE_i[rsj];
                //printf("bufM[%d]=%.5e\n",msi*Nsize+msj,cabs(bufM[msi*Nsize+msj]));
            }
        }
    }
    return;
}
