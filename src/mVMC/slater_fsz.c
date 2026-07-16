/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

his program is developed based on the mVMC-mini program
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

#include <stdint.h>

void UpdateSlaterElm_fsz();
void MakeSlaterElmBF_fsz(const int *eleNum, const int *eleProjBFCnt);
void MakeSlaterElmBF_fsz_to(double complex *sltElmBF, const int *eleNum, const int *eleProjBFCnt);
void MakeSlaterElmBF_fsz_to_serial(double complex *sltElmBF, const int *eleNum, const int *eleProjBFCnt);
void MakeSlaterElmBF_fsz_hop_to(double complex *sltElmBF, const double complex *baseSltElmBF,
                                const int *oldEleProjBFCnt, const int *newEleProjBFCnt);
void MakeSlaterElmBF_fsz_hop_to_serial(double complex *sltElmBF,
                                       const double complex *baseSltElmBF,
                                       const int *oldEleProjBFCnt,
                                       const int *newEleProjBFCnt);
int MakeSlaterElmBF_fsz_hop_to_with_rows(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected);
int MakeSlaterElmBF_fsz_hop_to_with_rows_serial(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected);
int GetSlaterElmBF_fsz_hop_int_work_size(int *workSize);
int MakeSlaterElmBF_fsz_hop_to_with_rows_workspace(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize);
int MakeSlaterElmBF_fsz_hop_to_with_rows_workspace_serial(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize);
void SlaterElmDiff_fsz(double complex *srOptO, const double complex ip, int *eleIdx,int *eleSpn);
void SlaterElmBFDiff_fsz(double complex *srOptO, const double complex ip,
                         const int *eleIdx, const int *eleSpn,
                         const int *eleNum, const int *eleProjBFCnt);
void BackFlowDiff_fsz(double complex *srOptO, const double complex ip,
                      const int *eleIdx, const int *eleSpn,
                      const int *eleNum, const int *eleProjBFCnt);

typedef struct {
  int sparseKeyCount;
  int sparseEntryCapacity;
  int sparseWorkSize;
  int *sparseWork;
  int *etaFlag;
  int *sparseOffset;
  int *sparseCount;
  int *sparseSite;
  int *sparseSubIdx;
  int *sparseThetaCnt;
} BFSparseWorkspace_fsz;

static inline int BFThetaCount_fsz(const int mu, const int centerSite,
                                   const int spin, const int rangeSlot,
                                   const int *eleProjBFCnt) {
  const int nSiteRange = Nsite*Nrange;
  const int channelOffset = (spin == 0) ? 0 : 4*nSiteRange;
  if(mu < 0 || mu >= 4) return 0;
  return eleProjBFCnt[channelOffset + mu*nSiteRange + centerSite*Nrange + rangeSlot];
}

static int BFHasActiveTheta_fsz(const int centerSite, const int spin,
                                const int *eleProjBFCnt) {
  int rangeSlot;
  /* Match the fcmp BackFlow eta gate: only the doublon-holon theta channel activates eta. */
  for(rangeSlot=0;rangeSlot<Nrange;rangeSlot++) {
    if(BFThetaCount_fsz(1, centerSite, spin, rangeSlot, eleProjBFCnt) != 0) {
      return 1;
    }
  }
  return 0;
}

static inline int BFThetaKey_fsz(const int spin, const int mu,
                                 const int centerSite) {
  return (spin*4 + mu)*Nsite + centerSite;
}

static inline double complex SlaterOrbital_fsz(const int ri, const int si,
                                               const int rj, const int sj) {
  const int rsi = ri + si*Nsite;
  const int rsj = rj + sj*Nsite;
  return Slater[OrbitalIdx[rsi][rsj]] * (double)OrbitalSgn[rsi][rsj];
}

static int BFSubIndexFromMuRange_fsz(const int mu, const int centerSite,
                                     const int rangeSite) {
  int d, xtmp, idx;
  d = RangeIdx[centerSite][rangeSite];
  xtmp = 4*d + mu;
  idx = xtmp - 3 - d;
  if(xtmp%4 == 0) idx = -1;
  if(xtmp == 0) idx = 0;
  return idx;
}

static void MakeBFEtaFlag_fsz(int *etaFlag, const int *eleProjBFCnt) {
  int spin, site;

  for(spin=0;spin<2;spin++) {
    for(site=0;site<Nsite;site++) {
      etaFlag[spin*Nsite + site] = BFHasActiveTheta_fsz(site, spin, eleProjBFCnt);
    }
  }
}

static void MakeBFSparseThetaList_fsz(int *sparseOffset, int *sparseCount,
                                      int *sparseSite, int *sparseSubIdx,
                                      int *sparseThetaCnt,
                                      const int *eleProjBFCnt) {
  int spin, mu, site, rangeSlot;
  int key, cursor = 0;
  int rangeSite, subIdx, thetaCnt;

  for(spin=0;spin<2;spin++) {
    for(mu=0;mu<4;mu++) {
      for(site=0;site<Nsite;site++) {
        key = BFThetaKey_fsz(spin, mu, site);
        sparseOffset[key] = cursor;
        sparseCount[key] = 0;

        for(rangeSlot=0;rangeSlot<Nrange;rangeSlot++) {
          thetaCnt = BFThetaCount_fsz(mu, site, spin, rangeSlot, eleProjBFCnt);
          if(thetaCnt == 0) continue;

          rangeSite = PosBF[site][rangeSlot];
          subIdx = BFSubIndexFromMuRange_fsz(mu, site, rangeSite);
          if(subIdx < 0) continue;

          sparseSite[cursor] = rangeSite;
          sparseSubIdx[cursor] = subIdx;
          sparseThetaCnt[cursor] = thetaCnt;
          cursor++;
          sparseCount[key]++;
        }
      }
    }
  }
}

static int GetBFSparseWorkspaceSize_fsz(int *sparseKeyCount,
                                        int *sparseEntryCapacity,
                                        int *sparseWorkSize) {
  int keyCount, entryCapacity, workSize;

  if(sparseKeyCount == NULL || sparseEntryCapacity == NULL || sparseWorkSize == NULL
      || Nsite <= 0 || Nrange <= 0 || Nsite > INT_MAX/8) {
    fprintf(stderr, "error: invalid BF-FSZ sparse workspace dimensions\n");
    return -1;
  }
  keyCount = 8*Nsite;
  if(Nrange > INT_MAX/keyCount) {
    fprintf(stderr, "error: BF-FSZ sparse workspace entry count overflow\n");
    return -1;
  }
  entryCapacity = keyCount*Nrange;
  if(Nsite > INT_MAX/2 || keyCount > (INT_MAX - 2*Nsite)/2) {
    fprintf(stderr, "error: BF-FSZ sparse workspace key size overflow\n");
    return -1;
  }
  workSize = 2*Nsite + 2*keyCount;
  if(entryCapacity > (INT_MAX - workSize)/3) {
    fprintf(stderr, "error: BF-FSZ sparse workspace size overflow\n");
    return -1;
  }
  workSize += 3*entryCapacity;
  *sparseKeyCount = keyCount;
  *sparseEntryCapacity = entryCapacity;
  *sparseWorkSize = workSize;
  return 0;
}

static int InitBFSparseWorkspaceFrom_fsz(BFSparseWorkspace_fsz *work,
                                         const int *eleProjBFCnt,
                                         int *sparseWork,
                                         const int sparseWorkCapacity) {
  if(work == NULL || eleProjBFCnt == NULL || sparseWork == NULL
      || GetBFSparseWorkspaceSize_fsz(&work->sparseKeyCount,
                                      &work->sparseEntryCapacity,
                                      &work->sparseWorkSize) != 0
      || sparseWorkCapacity < work->sparseWorkSize) {
    fprintf(stderr, "error: invalid BF-FSZ sparse workspace\n");
    return -1;
  }

  work->sparseWork = sparseWork;
  work->etaFlag = work->sparseWork;
  work->sparseOffset = work->etaFlag + 2*Nsite;
  work->sparseCount = work->sparseOffset + work->sparseKeyCount;
  work->sparseSite = work->sparseCount + work->sparseKeyCount;
  work->sparseSubIdx = work->sparseSite + work->sparseEntryCapacity;
  work->sparseThetaCnt = work->sparseSubIdx + work->sparseEntryCapacity;

  MakeBFEtaFlag_fsz(work->etaFlag, eleProjBFCnt);
  MakeBFSparseThetaList_fsz(work->sparseOffset, work->sparseCount,
                            work->sparseSite, work->sparseSubIdx,
                            work->sparseThetaCnt, eleProjBFCnt);
  return 0;
}

static void InitBFSparseWorkspace_fsz(BFSparseWorkspace_fsz *work,
                                      const int *eleProjBFCnt) {
  if(GetBFSparseWorkspaceSize_fsz(&work->sparseKeyCount,
                                  &work->sparseEntryCapacity,
                                  &work->sparseWorkSize) != 0) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if((size_t)work->sparseWorkSize > (size_t)-1/sizeof(int)) {
    fprintf(stderr, "error: BF-FSZ sparse workspace byte size overflow\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  work->sparseWork = (int *)malloc(sizeof(int)*(size_t)work->sparseWorkSize);
  if(work->sparseWork == NULL) {
    fprintf(stderr, "error: failed to allocate BF-FSZ sparse workspace\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(InitBFSparseWorkspaceFrom_fsz(work, eleProjBFCnt, work->sparseWork,
                                   work->sparseWorkSize) != 0) {
    free(work->sparseWork);
    work->sparseWork = NULL;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
}

static void FreeBFSparseWorkspace_fsz(BFSparseWorkspace_fsz *work) {
  free(work->sparseWork);
  work->sparseWork = NULL;
}

static double complex SubSlaterElmBF_fsz_sparse(
    const int ri, const int si, const int rj, const int sj,
    const int *xqp, const int *xqpOpt,
    const int *etaFlag, const int *sparseOffset, const int *sparseCount,
    const int *sparseSite, const int *sparseSubIdx,
    const int *sparseThetaCnt) {
  int mu, nu;
  int k, l, kStart, lStart, kCount, lCount;
  int rk, rki, rl, rlj, tri, trj;
  int nidx, midx, bfidx;
  int cntI, cntJ;
  double complex slt = 0.0 + 0.0*I;
  double eta;

  for(mu=0;mu<4;mu++) {
    kCount = sparseCount[BFThetaKey_fsz(si, mu, ri)];
    if(kCount == 0) continue;
    kStart = sparseOffset[BFThetaKey_fsz(si, mu, ri)];

    for(nu=0;nu<4;nu++) {
      if(mu == 0 && nu == 0) continue;
      lCount = sparseCount[BFThetaKey_fsz(sj, nu, rj)];
      if(lCount == 0) continue;
      lStart = sparseOffset[BFThetaKey_fsz(sj, nu, rj)];

      for(k=0;k<kCount;k++) {
        rk = sparseSite[kStart + k];
        rki = xqp[xqpOpt[rk]];
        nidx = sparseSubIdx[kStart + k];
        cntI = sparseThetaCnt[kStart + k];

        for(l=0;l<lCount;l++) {
          rl = sparseSite[lStart + l];
          rlj = xqp[xqpOpt[rl]];
          midx = sparseSubIdx[lStart + l];
          cntJ = sparseThetaCnt[lStart + l];
          bfidx = BFSubIdx[nidx][midx];
          slt += -ProjBF[bfidx] * (double)(cntI*cntJ)
                 * SlaterOrbital_fsz(rki, si, rlj, sj);
        }
      }
    }
  }

  eta = (etaFlag[si*Nsite + ri] || etaFlag[sj*Nsite + rj])
      ? creal(ProjBF[0]) : 1.0;
  tri = xqp[xqpOpt[ri]];
  trj = xqp[xqpOpt[rj]];
  slt += eta * SlaterOrbital_fsz(tri, si, trj, sj);

  return slt;
}

void UpdateSlaterElm_fsz() {
  int ri,ori,tri,sgni,rsi0,rsi1;
  int rj,orj,trj,sgnj,rsj0,rsj1;
  int qpidx,mpidx,optidx;
  int tri0,tri1,trj0,trj1; //fsz
  double complex slt_i0j0,slt_j0i0;
  double complex slt_i0j1,slt_j1i0;
  double complex slt_i1j0,slt_j0i1;
  double complex slt_i1j1,slt_j1i1;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex *sltE,*sltE_i0,*sltE_i1;

  #pragma omp parallel for default(shared)        \
    private(qpidx,optidx,mpidx,                   \
            xqpOpt,xqpOptSgn,xqp,xqpSgn,sltE,     \
            ri,ori,tri,sgni,rsi0,rsi1,sltE_i0,sltE_i1,      \
            rj,orj,trj,sgnj,rsj0,rsj1,slt_i0j0,slt_j0i0,slt_i0j1,slt_j1i0,slt_i1j0,slt_j0i1,slt_i1j1,slt_j1i1,\
            tri0,tri1,trj0,trj1)
  #pragma loop noalias
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    //printf("qpidx=%d \n",qpidx);
    // note: NQPFix = NSPGaussLeg * NMPTrans, NQPFull = NQPFix * NQPOptTrans(=1);
    // qpidx  = optidx*NQPFix+NSPGaussLeg*mpidx+spidx 
    // optidx will not be used
    optidx    = qpidx / NQPFix;              // optidx =0:  
    mpidx     = (qpidx%NQPFix) / NSPGaussLeg;// mpidx !=0: momentum projection
//  spidx     = qpidx % NSPGaussLeg;         // spidx  =0: spin projection

    xqpOpt    = QPOptTrans[optidx];          // xqpOpta   = 1
    xqpOptSgn = QPOptTransSgn[optidx];       // xqpOptSgn = 1
    xqp       = QPTrans[mpidx];              // xqp[*]=QPTrans[mpidx][*] 
    xqpSgn    = QPTransSgn[mpidx];           // xqpSgn[*]=QPTransSgn[mpidx][*] 
//    cs        = SPGLCosSin[spidx]; // spin projection is not used fsz
//    cc        = SPGLCosCos[spidx];
//    ss        = SPGLSinSin[spidx];
    
    sltE      = SlaterElm + qpidx*Nsite2*Nsite2;
    // for fsz: only for translational projection
    // spin projection is not implemented in ver.1.0
    for(ri=0;ri<Nsite;ri++) {
      ori     = xqpOpt[ri]; // ori = ri
      tri     = xqp[ori];   // ori -> tri momentum projection
      tri0    = tri;        // fsz up
      tri1    = tri0+Nsite; // fsz down
      sgni    = xqpSgn[ori]*xqpOptSgn[ri];
      rsi0    = ri;       // up
      rsi1    = ri+Nsite; // down
      sltE_i0 = sltE + rsi0*Nsite2; // sltE_i0[*]=sltE[rsi0*Nsite2][*]
      sltE_i1 = sltE + rsi1*Nsite2; // sltE_i1[*]=sltE[rsi1*Nsite2][*]
      
      for(rj=0;rj<Nsite;rj++) {
        //
        orj  = xqpOpt[rj];
        trj  = xqp[orj];
        //
        trj0 = trj;       // fsz up
        trj1 = trj+Nsite; // fsz down
        sgnj = xqpSgn[orj]*xqpOptSgn[rj];
        //
        rsj0 = rj;
        rsj1 = rj+Nsite;
        //
// for fsz        
        slt_i0j0        = Slater[ OrbitalIdx[tri0][trj0] ] * (double)(OrbitalSgn[tri0][trj0]*sgni*sgnj);
        slt_j0i0        = Slater[ OrbitalIdx[trj0][tri0] ] * (double)(OrbitalSgn[trj0][tri0]*sgni*sgnj);

        slt_i0j1        = Slater[ OrbitalIdx[tri0][trj1] ] * (double)(OrbitalSgn[tri0][trj1]*sgni*sgnj);
        slt_j1i0        = Slater[ OrbitalIdx[trj1][tri0] ] * (double)(OrbitalSgn[trj1][tri0]*sgni*sgnj);

        slt_i1j0        = Slater[ OrbitalIdx[tri1][trj0] ] * (double)(OrbitalSgn[tri1][trj0]*sgni*sgnj);
        slt_j0i1        = Slater[ OrbitalIdx[trj0][tri1] ] * (double)(OrbitalSgn[trj0][tri1]*sgni*sgnj);

        slt_i1j1        = Slater[ OrbitalIdx[tri1][trj1] ] * (double)(OrbitalSgn[tri1][trj1]*sgni*sgnj);
        slt_j1i1        = Slater[ OrbitalIdx[trj1][tri1] ] * (double)(OrbitalSgn[trj1][tri1]*sgni*sgnj);

        // F_{IJ}-F_{JI}
        sltE_i0[rsj0] =  slt_i0j0-slt_j0i0;   // up   - up
        sltE_i0[rsj1] =  slt_i0j1-slt_j1i0;   // up   - down
        sltE_i1[rsj0] =  slt_i1j0-slt_j0i1;   // down - up
        sltE_i1[rsj1] =  slt_i1j1-slt_j1i1;   // down - down 
      }// for rj 
    }// for ri
  }//for qpidx
  return;
}

static void MakeSlaterElmBF_fsz_qp(double complex *sltElmBF, const int qpidx,
                                   const int *etaFlag, const int *sparseOffset,
                                   const int *sparseCount, const int *sparseSite,
                                   const int *sparseSubIdx, const int *sparseThetaCnt) {
  int ri,ori,sgni,rsi;
  int rj,orj,sgnj,rsj;
  int si,sj;
  int mpidx,optidx;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex slt;
  double complex *sltE, *sltE_i;

  optidx    = qpidx / NQPFix;
  mpidx     = (qpidx%NQPFix) / NSPGaussLeg;

  xqpOpt    = QPOptTrans[optidx];
  xqpOptSgn = QPOptTransSgn[optidx];
  xqp       = QPTrans[mpidx];
  xqpSgn    = QPTransSgn[mpidx];

  sltE      = sltElmBF + qpidx*Nsite2*Nsite2;

#pragma loop noalias
  for(rsi=0;rsi<Nsite2;rsi++) {
    ri = rsi % Nsite;
    si = rsi / Nsite;
    ori  = xqpOpt[ri];
    sgni = xqpSgn[ori]*xqpOptSgn[ri];
    sltE_i = sltE + rsi*Nsite2;
    sltE_i[rsi] = 0.0 + 0.0*I;

    for(rsj=rsi+1;rsj<Nsite2;rsj++) {
      rj = rsj % Nsite;
      sj = rsj / Nsite;
      orj  = xqpOpt[rj];
      sgnj = xqpSgn[orj]*xqpOptSgn[rj];

      slt = (SubSlaterElmBF_fsz_sparse(ri, si, rj, sj, xqp, xqpOpt, etaFlag,
                                       sparseOffset, sparseCount,
                                       sparseSite, sparseSubIdx,
                                       sparseThetaCnt)
             - SubSlaterElmBF_fsz_sparse(rj, sj, ri, si, xqp, xqpOpt, etaFlag,
                                         sparseOffset, sparseCount,
                                         sparseSite, sparseSubIdx,
                                         sparseThetaCnt))
            * (double)(sgni*sgnj);
      sltE_i[rsj] = slt;
      sltE[rsj*Nsite2 + rsi] = -slt;
    }
  }
  return;
}

static void MakeSlaterElmBF_fsz_to_core(double complex *sltElmBF, const int *eleNum,
                                        const int *eleProjBFCnt, const int useOMP) {
  int qpidx;
  BFSparseWorkspace_fsz sparseWork;

  (void)eleNum;

  InitBFSparseWorkspace_fsz(&sparseWork, eleProjBFCnt);

  if(useOMP) {
    #pragma omp parallel for default(shared) private(qpidx)
    for(qpidx=0;qpidx<NQPFull;qpidx++) {
      MakeSlaterElmBF_fsz_qp(sltElmBF, qpidx, sparseWork.etaFlag,
                             sparseWork.sparseOffset, sparseWork.sparseCount,
                             sparseWork.sparseSite, sparseWork.sparseSubIdx,
                             sparseWork.sparseThetaCnt);
    }
  } else {
    for(qpidx=0;qpidx<NQPFull;qpidx++) {
      MakeSlaterElmBF_fsz_qp(sltElmBF, qpidx, sparseWork.etaFlag,
                             sparseWork.sparseOffset, sparseWork.sparseCount,
                             sparseWork.sparseSite, sparseWork.sparseSubIdx,
                             sparseWork.sparseThetaCnt);
    }
  }
  FreeBFSparseWorkspace_fsz(&sparseWork);
  return;
}

void MakeSlaterElmBF_fsz_to(double complex *sltElmBF, const int *eleNum, const int *eleProjBFCnt) {
  MakeSlaterElmBF_fsz_to_core(sltElmBF, eleNum, eleProjBFCnt, 1);
  return;
}

void MakeSlaterElmBF_fsz_to_serial(double complex *sltElmBF, const int *eleNum, const int *eleProjBFCnt) {
  MakeSlaterElmBF_fsz_to_core(sltElmBF, eleNum, eleProjBFCnt, 0);
  return;
}

void MakeSlaterElmBF_fsz(const int *eleNum, const int *eleProjBFCnt) {
  MakeSlaterElmBF_fsz_to(SlaterElmBF, eleNum, eleProjBFCnt);
  return;
}

static int MakeBFChangedEndpointList_fsz(int *changed, int *changedList,
                                         const int *oldEleProjBFCnt,
                                         const int *newEleProjBFCnt) {
  int spin, site, mu, rangeSlot;
  int nChanged = 0;

  for(spin=0;spin<2;spin++) {
    for(site=0;site<Nsite;site++) {
      const int rsi = spin*Nsite + site;
      changed[rsi] = 0;
      for(mu=0;mu<4 && changed[rsi] == 0;mu++) {
        for(rangeSlot=0;rangeSlot<Nrange;rangeSlot++) {
          if(BFThetaCount_fsz(mu, site, spin, rangeSlot, oldEleProjBFCnt)
              != BFThetaCount_fsz(mu, site, spin, rangeSlot, newEleProjBFCnt)) {
            changed[rsi] = 1;
            break;
          }
        }
      }
      if(changed[rsi] != 0) changedList[nChanged++] = rsi;
    }
  }
  return nChanged;
}

static int MakeBFAffectedParticleList_fsz(
    const int *changed, int movedParticle,
    const int *eleIdx, const int *eleSpn,
    int *affectedMark, int *affected) {
  int m, rsm, nAffected = 0;

  memset(affectedMark, 0, sizeof(int)*(size_t)Nsize);
  affectedMark[movedParticle] = 1;
  for(m=0;m<Nsize;m++) {
    if(eleIdx[m] < 0 || eleIdx[m] >= Nsite || eleSpn[m] < 0 || eleSpn[m] >= 2) {
      fprintf(stderr,
              "error: invalid BF-FSZ candidate particle %d: site=%d spin=%d\n",
              m, eleIdx[m], eleSpn[m]);
      return -1;
    }
    rsm = eleIdx[m] + eleSpn[m]*Nsite;
    if(changed[rsm] != 0) affectedMark[m] = 1;
  }
  for(m=0;m<Nsize;m++) {
    if(affectedMark[m] != 0) affected[nAffected++] = m;
  }
  if(nAffected < 1 || nAffected > Nsize) {
    fprintf(stderr, "error: invalid BF-FSZ affected particle count %d\n", nAffected);
    return -1;
  }
  return nAffected;
}

static int BF_FSZ_IntWorkspacesOverlap_fsz(const int *workA, const int workACount,
                                           const int *workB, const int workBCount) {
  size_t workABytes, workBBytes;
  uintptr_t workABegin, workBBegin;

  if(workA == NULL || workB == NULL || workACount < 0 || workBCount < 0
      || (size_t)workACount > SIZE_MAX/sizeof(int)
      || (size_t)workBCount > SIZE_MAX/sizeof(int)) {
    return 1;
  }
  workABytes = sizeof(int)*(size_t)workACount;
  workBBytes = sizeof(int)*(size_t)workBCount;
  workABegin = (uintptr_t)(const void *)workA;
  workBBegin = (uintptr_t)(const void *)workB;
  if(workABegin > UINTPTR_MAX - workABytes
      || workBBegin > UINTPTR_MAX - workBBytes) return 1;
  return workABegin < workBBegin + workBBytes
      && workBBegin < workABegin + workABytes;
}

static inline void MakeSlaterElmBF_fsz_hop_pair(
    double complex *sltE, const int rsi, const int rsj,
    const int *xqp, const int *xqpSgn, const int *xqpOpt,
    const int *xqpOptSgn, const BFSparseWorkspace_fsz *sparseWork) {
  const int ri = rsi % Nsite;
  const int si = rsi / Nsite;
  const int rj = rsj % Nsite;
  const int sj = rsj / Nsite;
  const int ori = xqpOpt[ri];
  const int orj = xqpOpt[rj];
  const int sgni = xqpSgn[ori]*xqpOptSgn[ri];
  const int sgnj = xqpSgn[orj]*xqpOptSgn[rj];
  const double complex slt =
      (SubSlaterElmBF_fsz_sparse(ri, si, rj, sj, xqp, xqpOpt,
                                 sparseWork->etaFlag, sparseWork->sparseOffset,
                                 sparseWork->sparseCount, sparseWork->sparseSite,
                                 sparseWork->sparseSubIdx, sparseWork->sparseThetaCnt)
       - SubSlaterElmBF_fsz_sparse(rj, sj, ri, si, xqp, xqpOpt,
                                   sparseWork->etaFlag, sparseWork->sparseOffset,
                                   sparseWork->sparseCount, sparseWork->sparseSite,
                                   sparseWork->sparseSubIdx, sparseWork->sparseThetaCnt))
      * (double)(sgni*sgnj);

  sltE[(size_t)rsi*(size_t)Nsite2 + (size_t)rsj] = slt;
  sltE[(size_t)rsj*(size_t)Nsite2 + (size_t)rsi] = -slt;
}

static void MakeSlaterElmBF_fsz_hop_qp(double complex *sltElmBF,
                                       const double complex *baseSltElmBF,
                                       const int qpidx, const int *changed,
                                       const int *changedList, const int nChanged,
                                       const BFSparseWorkspace_fsz *sparseWork) {
  int changedIdx, rsi, rsj;
  int mpidx, optidx;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  const size_t qpOffset = (size_t)qpidx*(size_t)Nsite2*(size_t)Nsite2;
  double complex *sltE = sltElmBF + qpOffset;
  const double complex *baseSltE = baseSltElmBF + qpOffset;

  memcpy(sltE, baseSltE, sizeof(double complex)*(size_t)Nsite2*(size_t)Nsite2);
  optidx = qpidx / NQPFix;
  mpidx = (qpidx%NQPFix) / NSPGaussLeg;
  xqpOpt = QPOptTrans[optidx];
  xqpOptSgn = QPOptTransSgn[optidx];
  xqp = QPTrans[mpidx];
  xqpSgn = QPTransSgn[mpidx];

  /* Pass 1: pairs whose lower endpoint changed. */
  for(changedIdx=0;changedIdx<nChanged;changedIdx++) {
    rsi = changedList[changedIdx];
    for(rsj=rsi+1;rsj<Nsite2;rsj++) {
      MakeSlaterElmBF_fsz_hop_pair(sltE, rsi, rsj, xqp, xqpSgn, xqpOpt,
                                    xqpOptSgn, sparseWork);
    }
  }

  /* Pass 2: pairs whose upper endpoint changed and lower endpoint did not. */
  for(changedIdx=0;changedIdx<nChanged;changedIdx++) {
    rsj = changedList[changedIdx];
    for(rsi=0;rsi<rsj;rsi++) {
      if(changed[rsi] != 0) continue;
      MakeSlaterElmBF_fsz_hop_pair(sltE, rsi, rsj, xqp, xqpSgn, xqpOpt,
                                    xqpOptSgn, sparseWork);
    }
  }
}

int GetSlaterElmBF_fsz_hop_int_work_size(int *workSize) {
  int sparseKeyCount, sparseEntryCapacity, sparseWorkSize;

  if(workSize == NULL || Nsite2 <= 0 || Nsize <= 0
      || GetBFSparseWorkspaceSize_fsz(&sparseKeyCount, &sparseEntryCapacity,
                                      &sparseWorkSize) != 0
      || Nsite2 != 2*Nsite || Nsize > Nsite2) {
    fprintf(stderr, "error: invalid BF-FSZ hop workspace dimensions\n");
    return -1;
  }
  (void)sparseKeyCount;
  (void)sparseEntryCapacity;
  if(Nsize > INT_MAX - sparseWorkSize
      || Nsite2 > (INT_MAX - sparseWorkSize - Nsize)/2) {
    fprintf(stderr, "error: BF-FSZ hop workspace size overflow\n");
    return -1;
  }
  *workSize = 2*Nsite2 + Nsize + sparseWorkSize;
  return 0;
}

static int MakeSlaterElmBF_fsz_hop_to_core(double complex *sltElmBF,
                                           const double complex *baseSltElmBF,
                                           const int movedParticle,
                                           const int *eleIdx, const int *eleSpn,
                                           const int *oldEleProjBFCnt,
                                           const int *newEleProjBFCnt,
                                           int *affected, int *nChangedOut,
                                           int *nAffectedOut,
                                           int *intWork, const int intWorkSize,
                                           const int useOMP) {
  const size_t matrixSize = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const int collectRows = (affected != NULL || nChangedOut != NULL || nAffectedOut != NULL);
  int qpidx, m, n, rsm, rsn, nChanged, nAffected = 0;
  int requiredWorkSize, sparseWorkOffset;
  size_t idx;
  int *changed, *changedList, *affectedMark, *sparseIntWork;
  double complex *fullSltElmBF;
  BFSparseWorkspace_fsz sparseWork;

  if(sltElmBF == NULL || baseSltElmBF == NULL || oldEleProjBFCnt == NULL
      || newEleProjBFCnt == NULL || sltElmBF == baseSltElmBF) {
    fprintf(stderr, "error: invalid BF-FSZ incremental Slater arguments\n");
    return -1;
  }
  if(collectRows && (affected == NULL || nChangedOut == NULL || nAffectedOut == NULL
      || eleIdx == NULL || eleSpn == NULL || movedParticle < 0 || movedParticle >= Nsize)) {
    fprintf(stderr, "error: invalid BF-FSZ affected-row arguments\n");
    return -1;
  }
  if(BFFSZAffectedCheckEnabled && !collectRows) {
    fprintf(stderr, "error: BF-FSZ affected oracle requires affected-row arguments\n");
    return -1;
  }
  if(GetSlaterElmBF_fsz_hop_int_work_size(&requiredWorkSize) != 0
      || intWork == NULL || intWorkSize < requiredWorkSize) {
    fprintf(stderr, "error: invalid BF-FSZ hop integer workspace\n");
    return -1;
  }
  if(collectRows
      && BF_FSZ_IntWorkspacesOverlap_fsz(affected, Nsize,
                                         intWork, requiredWorkSize)) {
    fprintf(stderr, "error: BF-FSZ affected and hop workspaces overlap\n");
    return -1;
  }

  changed = intWork;
  changedList = changed + Nsite2;
  affectedMark = changedList + Nsite2;
  sparseWorkOffset = 2*Nsite2 + Nsize;
  sparseIntWork = intWork + sparseWorkOffset;

  nChanged = MakeBFChangedEndpointList_fsz(changed, changedList,
                                            oldEleProjBFCnt, newEleProjBFCnt);
  if(collectRows) {
    nAffected = MakeBFAffectedParticleList_fsz(changed, movedParticle, eleIdx, eleSpn,
                                                affectedMark, affected);
    if(nAffected < 0) return -1;
    *nChangedOut = nChanged;
    *nAffectedOut = nAffected;
  }
  if(InitBFSparseWorkspaceFrom_fsz(&sparseWork, newEleProjBFCnt,
                                    sparseIntWork, intWorkSize - sparseWorkOffset) != 0) {
    return -1;
  }

  if(useOMP) {
    #pragma omp parallel for default(shared) private(qpidx)
    for(qpidx=0;qpidx<NQPFull;qpidx++) {
      MakeSlaterElmBF_fsz_hop_qp(sltElmBF, baseSltElmBF, qpidx, changed,
                                  changedList, nChanged, &sparseWork);
    }
  } else {
    for(qpidx=0;qpidx<NQPFull;qpidx++) {
      MakeSlaterElmBF_fsz_hop_qp(sltElmBF, baseSltElmBF, qpidx, changed,
                                  changedList, nChanged, &sparseWork);
    }
  }

  if(BFFSZGreenRebuildCheckEnabled || BFFSZAffectedCheckEnabled) {
    fullSltElmBF = (double complex *)malloc(sizeof(double complex)*matrixSize);
    if(fullSltElmBF == NULL) {
      fprintf(stderr, "error: failed to allocate BF-FSZ full-rebuild oracle\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    MakeSlaterElmBF_fsz_to_core(fullSltElmBF, NULL, newEleProjBFCnt, useOMP);
    if(BFFSZGreenRebuildCheckEnabled) {
      for(idx=0;idx<matrixSize;idx++) {
        if(sltElmBF[idx] != fullSltElmBF[idx]) {
          fprintf(stderr,
                  "error: BF-FSZ incremental Green rebuild mismatch at index %zu: "
                  "incremental=(%.17e,%.17e) full=(%.17e,%.17e)\n",
                  idx, creal(sltElmBF[idx]), cimag(sltElmBF[idx]),
                  creal(fullSltElmBF[idx]), cimag(fullSltElmBF[idx]));
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
      }
    }
    if(BFFSZAffectedCheckEnabled) {
      for(qpidx=0;qpidx<NQPFull;qpidx++) {
        const size_t qpOffset = (size_t)qpidx*(size_t)Nsite2*(size_t)Nsite2;
        const double complex *baseSltE = baseSltElmBF + qpOffset;
        const double complex *fullSltE = fullSltElmBF + qpOffset;
        for(m=0;m<Nsize;m++) {
          rsm = eleIdx[m] + eleSpn[m]*Nsite;
          for(n=m+1;n<Nsize;n++) {
            rsn = eleIdx[n] + eleSpn[n]*Nsite;
            if(baseSltE[rsm*Nsite2 + rsn] != fullSltE[rsm*Nsite2 + rsn]
                && affectedMark[m] == 0 && affectedMark[n] == 0) {
              fprintf(stderr,
                      "error: BF-FSZ affected coverage mismatch at qp=%d particles=(%d,%d) "
                      "endpoints=(%d,%d) nChanged=%d nAffected=%d\n",
                      qpidx, m, n, rsm, rsn, nChanged, nAffected);
              MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
          }
        }
      }
    }
    free(fullSltElmBF);
  }

  return 0;
}

static int MakeSlaterElmBF_fsz_hop_to_allocated(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    const int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected, const int useOMP) {
  int info, intWorkSize;
  int *intWork;

  if(GetSlaterElmBF_fsz_hop_int_work_size(&intWorkSize) != 0
      || (size_t)intWorkSize > (size_t)-1/sizeof(int)) {
    return -1;
  }
  intWork = (int *)malloc(sizeof(int)*(size_t)intWorkSize);
  if(intWork == NULL) {
    fprintf(stderr, "error: failed to allocate BF-FSZ hop integer workspace\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  info = MakeSlaterElmBF_fsz_hop_to_core(
      sltElmBF, baseSltElmBF, movedParticle, eleIdx, eleSpn,
      oldEleProjBFCnt, newEleProjBFCnt, affected, nChanged, nAffected,
      intWork, intWorkSize, useOMP);
  free(intWork);
  return info;
}

void MakeSlaterElmBF_fsz_hop_to(double complex *sltElmBF,
                                const double complex *baseSltElmBF,
                                const int *oldEleProjBFCnt,
                                const int *newEleProjBFCnt) {
  if(MakeSlaterElmBF_fsz_hop_to_allocated(sltElmBF, baseSltElmBF, -1, NULL, NULL,
                                          oldEleProjBFCnt, newEleProjBFCnt,
                                          NULL, NULL, NULL, 1) != 0) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
}

void MakeSlaterElmBF_fsz_hop_to_serial(double complex *sltElmBF,
                                       const double complex *baseSltElmBF,
                                       const int *oldEleProjBFCnt,
                                       const int *newEleProjBFCnt) {
  if(MakeSlaterElmBF_fsz_hop_to_allocated(sltElmBF, baseSltElmBF, -1, NULL, NULL,
                                          oldEleProjBFCnt, newEleProjBFCnt,
                                          NULL, NULL, NULL, 0) != 0) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
}

int MakeSlaterElmBF_fsz_hop_to_with_rows(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected) {
  return MakeSlaterElmBF_fsz_hop_to_allocated(
      sltElmBF, baseSltElmBF, movedParticle, eleIdx, eleSpn,
      oldEleProjBFCnt, newEleProjBFCnt, affected, nChanged, nAffected, 1);
}

int MakeSlaterElmBF_fsz_hop_to_with_rows_serial(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected) {
  return MakeSlaterElmBF_fsz_hop_to_allocated(
      sltElmBF, baseSltElmBF, movedParticle, eleIdx, eleSpn,
      oldEleProjBFCnt, newEleProjBFCnt, affected, nChanged, nAffected, 0);
}

int MakeSlaterElmBF_fsz_hop_to_with_rows_workspace(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize) {
  return MakeSlaterElmBF_fsz_hop_to_core(
      sltElmBF, baseSltElmBF, movedParticle, eleIdx, eleSpn,
      oldEleProjBFCnt, newEleProjBFCnt, affected, nChanged, nAffected,
      intWork, intWorkSize, 1);
}

int MakeSlaterElmBF_fsz_hop_to_with_rows_workspace_serial(
    double complex *sltElmBF, const double complex *baseSltElmBF,
    int movedParticle, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int *affected, int *nChanged, int *nAffected,
    int *intWork, int intWorkSize) {
  return MakeSlaterElmBF_fsz_hop_to_core(
      sltElmBF, baseSltElmBF, movedParticle, eleIdx, eleSpn,
      oldEleProjBFCnt, newEleProjBFCnt, affected, nChanged, nAffected,
      intWork, intWorkSize, 0);
}

static inline void AddSlaterDiff_fsz(double complex *buf,
                                     const double complex pref,
                                     const int ri, const int si,
                                     const int rj, const int sj,
                                     const double complex scale) {
  const int rsi = ri + si*Nsite;
  const int rsj = rj + sj*Nsite;
  const int orbidx = OrbitalIdx[rsi][rsj];
  buf[orbidx] += pref * scale * (double)OrbitalSgn[rsi][rsj];
}

static void AccumulateSlaterBFDiffOriented_fsz(
    double complex *buf, const double complex pref,
    const int ri, const int si, const int rj, const int sj,
    const int *xqp, const int *xqpOpt,
    const BFSparseWorkspace_fsz *sparseWork) {
  int mu, nu;
  int k, l, kStart, lStart, kCount, lCount;
  int rk, rl, rki, rlj, nidx, midx, bfidx;
  int cntI, cntJ;
  int tri, trj;
  double eta;

  for(mu=0;mu<4;mu++) {
    kCount = sparseWork->sparseCount[BFThetaKey_fsz(si, mu, ri)];
    if(kCount == 0) continue;
    kStart = sparseWork->sparseOffset[BFThetaKey_fsz(si, mu, ri)];

    for(nu=0;nu<4;nu++) {
      if(mu == 0 && nu == 0) continue;
      lCount = sparseWork->sparseCount[BFThetaKey_fsz(sj, nu, rj)];
      if(lCount == 0) continue;
      lStart = sparseWork->sparseOffset[BFThetaKey_fsz(sj, nu, rj)];

      for(k=0;k<kCount;k++) {
        rk = sparseWork->sparseSite[kStart + k];
        rki = xqp[xqpOpt[rk]];
        nidx = sparseWork->sparseSubIdx[kStart + k];
        cntI = sparseWork->sparseThetaCnt[kStart + k];

        for(l=0;l<lCount;l++) {
          rl = sparseWork->sparseSite[lStart + l];
          rlj = xqp[xqpOpt[rl]];
          midx = sparseWork->sparseSubIdx[lStart + l];
          cntJ = sparseWork->sparseThetaCnt[lStart + l];
          bfidx = BFSubIdx[nidx][midx];
          AddSlaterDiff_fsz(buf, pref, rki, si, rlj, sj,
                            -ProjBF[bfidx] * (double)(cntI*cntJ));
        }
      }
    }
  }

  eta = (sparseWork->etaFlag[si*Nsite + ri]
      || sparseWork->etaFlag[sj*Nsite + rj]) ? creal(ProjBF[0]) : 1.0;
  tri = xqp[xqpOpt[ri]];
  trj = xqp[xqpOpt[rj]];
  AddSlaterDiff_fsz(buf, pref, tri, si, trj, sj, eta);
}

static void AccumulateBackFlowDiffOriented_fsz(
    double complex *bfReal, double complex *bfImag, const double complex pref,
    const int ri, const int si, const int rj, const int sj,
    const int *xqp, const int *xqpOpt,
    const BFSparseWorkspace_fsz *sparseWork) {
  int mu, nu;
  int k, l, kStart, lStart, kCount, lCount;
  int rk, rl, rki, rlj, nidx, midx, bfidx;
  int cntI, cntJ;
  int tri, trj;
  double complex term;

  for(mu=0;mu<4;mu++) {
    kCount = sparseWork->sparseCount[BFThetaKey_fsz(si, mu, ri)];
    if(kCount == 0) continue;
    kStart = sparseWork->sparseOffset[BFThetaKey_fsz(si, mu, ri)];

    for(nu=0;nu<4;nu++) {
      if(mu == 0 && nu == 0) continue;
      lCount = sparseWork->sparseCount[BFThetaKey_fsz(sj, nu, rj)];
      if(lCount == 0) continue;
      lStart = sparseWork->sparseOffset[BFThetaKey_fsz(sj, nu, rj)];

      for(k=0;k<kCount;k++) {
        rk = sparseWork->sparseSite[kStart + k];
        rki = xqp[xqpOpt[rk]];
        nidx = sparseWork->sparseSubIdx[kStart + k];
        cntI = sparseWork->sparseThetaCnt[kStart + k];

        for(l=0;l<lCount;l++) {
          rl = sparseWork->sparseSite[lStart + l];
          rlj = xqp[xqpOpt[rl]];
          midx = sparseWork->sparseSubIdx[lStart + l];
          cntJ = sparseWork->sparseThetaCnt[lStart + l];
          bfidx = BFSubIdx[nidx][midx];
          term = -(double)(cntI*cntJ) * SlaterOrbital_fsz(rki, si, rlj, sj);
          bfReal[bfidx] += pref * term;
          bfImag[bfidx] += pref * term * I;
        }
      }
    }
  }

  if(sparseWork->etaFlag[si*Nsite + ri]
      || sparseWork->etaFlag[sj*Nsite + rj]) {
    tri = xqp[xqpOpt[ri]];
    trj = xqp[xqpOpt[rj]];
    /* The BF-FSZ builder uses creal(ProjBF[0]) for eta, so eta has no
       derivative with respect to Im(ProjBF[0]); theta terms above still do. */
    bfReal[0] += pref * SlaterOrbital_fsz(tri, si, trj, sj);
  }
}

void SlaterElmBFDiff_fsz(double complex *srOptO, const double complex ip,
                         const int *eleIdx, const int *eleSpn,
                         const int *eleNum, const int *eleProjBFCnt) {
  const int nBuf = NSlater*NQPFull;
  const double complex invIP = 1.0/ip;
  BFSparseWorkspace_fsz sparseWork;
  double complex *buffer;
  int i, qpidx, optidx, mpidx, msi, msj;
  int ri, rj, si, sj, ori, orj, sgni, sgnj;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex *buf, *invM, *invM_i;
  double complex pref, tmp;
  int orbidx;

  (void)eleNum;

  for(orbidx=0;orbidx<NSlater;orbidx++) {
    srOptO[2*orbidx] = 0.0 + 0.0*I;
    srOptO[2*orbidx+1] = 0.0 + 0.0*I;
  }
  if(NSlater <= 0) return;

  InitBFSparseWorkspace_fsz(&sparseWork, eleProjBFCnt);
  RequestWorkSpaceComplex(NQPFull*NSlater);
  buffer = GetWorkSpaceComplex(NQPFull*NSlater);
  for(i=0;i<nBuf;i++) buffer[i] = 0.0 + 0.0*I;

  #pragma omp parallel for default(shared)                         \
    private(qpidx,optidx,mpidx,xqp,xqpSgn,xqpOpt,xqpOptSgn,        \
            invM,buf,msi,msj,invM_i,ri,rj,si,sj,ori,orj,sgni,sgnj, \
            pref)
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    optidx = qpidx / NQPFix;
    mpidx = (qpidx%NQPFix) / NSPGaussLeg;
    xqpOpt = QPOptTrans[optidx];
    xqpOptSgn = QPOptTransSgn[optidx];
    xqp = QPTrans[mpidx];
    xqpSgn = QPTransSgn[mpidx];
    invM = InvM + qpidx*Nsize*Nsize;
    buf = buffer + qpidx*NSlater;

    for(msi=0;msi<Nsize;msi++) {
      ri = eleIdx[msi];
      si = eleSpn[msi];
      ori = xqpOpt[ri];
      sgni = xqpSgn[ori]*xqpOptSgn[ri];
      invM_i = invM + msi*Nsize;

      for(msj=0;msj<Nsize;msj++) {
        rj = eleIdx[msj];
        sj = eleSpn[msj];
        orj = xqpOpt[rj];
        sgnj = xqpSgn[orj]*xqpOptSgn[rj];
        pref = -invM_i[msj] * PfM[qpidx] * (double)(sgni*sgnj);
        AccumulateSlaterBFDiffOriented_fsz(buf, pref, ri, si, rj, sj,
                                           xqp, xqpOpt, &sparseWork);
      }
    }
  }

  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    tmp = QPFullWeight[qpidx];
    buf = buffer + qpidx*NSlater;
    for(orbidx=0;orbidx<NSlater;orbidx++) {
      srOptO[2*orbidx] += tmp * buf[orbidx];
      srOptO[2*orbidx+1] += tmp * buf[orbidx] * I;
    }
  }
  for(orbidx=0;orbidx<NSlater;orbidx++) {
    srOptO[2*orbidx] *= invIP;
    srOptO[2*orbidx+1] *= invIP;
  }

  ReleaseWorkSpaceComplex();
  FreeBFSparseWorkspace_fsz(&sparseWork);
  return;
}

void BackFlowDiff_fsz(double complex *srOptO, const double complex ip,
                      const int *eleIdx, const int *eleSpn,
                      const int *eleNum, const int *eleProjBFCnt) {
  const double complex invIP = 1.0/ip;
  BFSparseWorkspace_fsz sparseWork;
  double complex *buffer;
  double complex *bfReal, *bfImag;
  int i, qpidx, optidx, mpidx, msi, msj;
  int ri, rj, si, sj, ori, orj, sgni, sgnj;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex *invM, *invM_i;
  double complex pref, tmp;
  int bfidx;

  (void)eleNum;

  for(bfidx=0;bfidx<NProjBF;bfidx++) {
    srOptO[2*bfidx] = 0.0 + 0.0*I;
    srOptO[2*bfidx+1] = 0.0 + 0.0*I;
  }
  if(NProjBF <= 0) return;

  InitBFSparseWorkspace_fsz(&sparseWork, eleProjBFCnt);
  RequestWorkSpaceComplex(2*NQPFull*NProjBF);
  buffer = GetWorkSpaceComplex(2*NQPFull*NProjBF);
  for(i=0;i<2*NQPFull*NProjBF;i++) buffer[i] = 0.0 + 0.0*I;

  #pragma omp parallel for default(shared)                         \
    private(qpidx,optidx,mpidx,xqp,xqpSgn,xqpOpt,xqpOptSgn,        \
            invM,bfReal,bfImag,msi,msj,invM_i,ri,rj,si,sj,ori,orj, \
            sgni,sgnj,pref)
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    optidx = qpidx / NQPFix;
    mpidx = (qpidx%NQPFix) / NSPGaussLeg;
    xqpOpt = QPOptTrans[optidx];
    xqpOptSgn = QPOptTransSgn[optidx];
    xqp = QPTrans[mpidx];
    xqpSgn = QPTransSgn[mpidx];
    invM = InvM + qpidx*Nsize*Nsize;
    bfReal = buffer + qpidx*2*NProjBF;
    bfImag = bfReal + NProjBF;

    for(msi=0;msi<Nsize;msi++) {
      ri = eleIdx[msi];
      si = eleSpn[msi];
      ori = xqpOpt[ri];
      sgni = xqpSgn[ori]*xqpOptSgn[ri];
      invM_i = invM + msi*Nsize;

      for(msj=0;msj<Nsize;msj++) {
        rj = eleIdx[msj];
        sj = eleSpn[msj];
        orj = xqpOpt[rj];
        sgnj = xqpSgn[orj]*xqpOptSgn[rj];
        pref = -invM_i[msj] * PfM[qpidx] * (double)(sgni*sgnj);
        AccumulateBackFlowDiffOriented_fsz(bfReal, bfImag, pref, ri, si,
                                           rj, sj, xqp, xqpOpt, &sparseWork);
      }
    }
  }

  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    tmp = QPFullWeight[qpidx] * invIP;
    bfReal = buffer + qpidx*2*NProjBF;
    bfImag = bfReal + NProjBF;
    for(bfidx=0;bfidx<NProjBF;bfidx++) {
      srOptO[2*bfidx] += tmp * bfReal[bfidx];
      srOptO[2*bfidx+1] += tmp * bfImag[bfidx];
    }
  }

  ReleaseWorkSpaceComplex();
  FreeBFSparseWorkspace_fsz(&sparseWork);
  return;
}

// Calculating Tr[Inv[M]*D_k(X)]
void SlaterElmDiff_fsz(double complex *srOptO, const double complex ip, int *eleIdx,int *eleSpn) {
  const int nBuf=NSlater*NQPFull;
  const int nsize = Nsize;       // Nsize=2*Ne
  const int nQPFull = NQPFull;   // note: NQPFix = NSPGaussLeg * NMPTrans, NQPFull = NQPFix * NQPOptTrans(=1);
  const int nMPTrans = NMPTrans; // number of translational operators
  const int nSlater = NSlater;   // number of slater
  const int nTrans = NMPTrans * NQPOptTrans; //usually NQPOptTrans=1

  const double complex invIP = 1.0/ip;
  int msi,msj,ri,rj,ori,orj,tri,trj,sgni,sgnj;
  int mpidx,orbidx,qpidx,optidx,i;
  int si,sj;//fsz
  double complex cs; // including Pf
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
            ri,ori,tri,sgni,rj,orj,trj,sgnj,si,sj,                       \
            tOrbIdx,tOrbIdx_i,tOrbSgn,tOrbSgn_i,orbitalIdx_i)
  #pragma loop noalias
// this loop for making transObrIdx,transOrbSgn
// transOrbIdx (NMPTrans * NQPOptTrans)*(nsize*nsize) nsize=2Ne
// transOrbSgn (NMPTrans * NQPOptTrans)*(nsize*nsize) nsize=2Ne
  for(qpidx=0;qpidx<nTrans;qpidx++) {//only for trans op. // nTran=nMPTrans*(NQPOptTrans=1)=nMPTrans
    optidx    = qpidx / nMPTrans;                // qpidx=mpidx+nMPTrans*optidx
    mpidx     = qpidx % nMPTrans; 
                                                 // QPOptTrans:     NQPOptTrans* Nsite matrix :  usually NQPOptTrans = 1
                                                 // QPOptTransSgn:  NQPOptTrans* Nsite matrix :  usually NQPOptTrans = 1 
    xqpOpt    = QPOptTrans[optidx];              // xqpOpt[]    = QPOptTrans[optidx][]    : usually optidx=0 and not be used
    xqpOptSgn = QPOptTransSgn[optidx];           // xqpOptSgn[] = QPOptTransSgn[optidx][] : usually optidx=0 and not be used
    xqp       = QPTrans[mpidx];                  // xqp    =  QPTrans[mpidx][]
    xqpSgn    = QPTransSgn[mpidx];               // xqpSgn =  QPTransSgn[mpidx][]
    tOrbIdx   = transOrbIdx + qpidx*nsize*nsize; // tOrbIdx[msi][msj]=OrbitalIdx[tri][trj]
    tOrbSgn   = transOrbSgn + qpidx*nsize*nsize; // tOrbSgn : sign of f_ij
    for(msi=0;msi<nsize;msi++) {                 // nsize=2*Ne including spin degrees of freedom
      ri           = eleIdx[msi];                //  ri  : postion where the msi-th electron exists  
      si           = eleSpn[msi];                //  si  : spin of the msi-th electron  :fsz
      ori          = xqpOpt[ri];                 // ori  : ri -> ori -> tri
      tri          = xqp[ori]+si*Nsite;          // tri  : fsz
      sgni         = xqpSgn[ori]*xqpOptSgn[ri];  // sgni :
      tOrbIdx_i    = tOrbIdx + msi*nsize;        // 
      tOrbSgn_i    = tOrbSgn + msi*nsize;        //
      orbitalIdx_i = OrbitalIdx[tri];            //
      for(msj=0;msj<nsize;msj++) { //nsize=2*Ne //
        rj             = eleIdx[msj];           //
        sj             = eleSpn[msj];        // sj spin of the msj-th electron  :fsz
        orj            = xqpOpt[rj];         
        trj            = xqp[orj]+sj*Nsite;       //fsz
        sgnj           = xqpSgn[orj]*xqpOptSgn[rj]; //xqpSgn
        tOrbIdx_i[msj] = orbitalIdx_i[trj];
        tOrbSgn_i[msj] = sgni*sgnj*OrbitalSgn[tri][trj];
      }
    }
  }
// calculating Tr(X^{-1}*dX/df_{msi,msj})=-2*alpha(sigma(msi),sigma(msj))(X^{-1})_{msi,msj}
  #pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,cs,                   \
            tOrbIdx,tOrbSgn,invM,buf,msi,msj,             \
            tOrbIdx_i,tOrbSgn_i,invM_i,orbidx)
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) { // nQPFull = NQPFix * NQPOptTrans: usually = NSPGaussLeg * NMPTrans
    mpidx = qpidx / NSPGaussLeg;       // qpidx   = NSPGaussLeg*mpidx + spidx

    cs = PfM[qpidx];// * SPGLCosSin[spidx]; //  PfM

    tOrbIdx = transOrbIdx + mpidx*nsize*nsize; // tOrbIdx[msi][msj] = transOrbIdx[mpidx][msi][msj]
    tOrbSgn = transOrbSgn + mpidx*nsize*nsize; // tOrbSgn[msi][msj] = transOrbSgn[mpidx][msi][msj]
    invM    = InvM        + qpidx*Nsize*Nsize; // invM  M^-1 and PfM, InvM: NQPFull*(Nsize*Nsize)+PfM
    buf     = buffer      + qpidx*NSlater;     // buf   fij, buffer: NSlater*NQPFull

    #pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) { //fsz
      tOrbIdx_i = tOrbIdx + msi*nsize;         // tOrbIdx_i[*] = tOrbIdx[msi][*]
      tOrbSgn_i = tOrbSgn + msi*nsize;         // tOrbSgn_i[*] = tOrbSgn[msi][*]
      invM_i    = invM + msi*nsize;            // invM[]       = invM[msi][]
      for(msj=0;msj<nsize;msj++) {             // 
        orbidx       = tOrbIdx_i[msj];         // 
        buf[orbidx] += -1.0*invM_i[msj]*cs*tOrbSgn_i[msj]; // invM[msi][msj]// 2*(1/2)=1
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
