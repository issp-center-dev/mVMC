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
 * Stochastic Reconfiguration method by DPOSV
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#ifdef _SYSTEM_A
 #define M_DPOSV  DPOSV
#elif _lapack_small_nounderscore
 #define M_DPOSV  dposv
#else
 #define M_DPOSV  dposv_
#endif

int M_DPOSV(char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);

int StochasticOpt(MPI_Comm comm);
void stcOptInit(double *const s, double *const g, const int nSmat, int *const smatToParaIdx);
int stcOptMain(double *const s, double *const g, const int nSmat);

int StochasticOpt(MPI_Comm comm) {
  double *s; /* the overlap matrix S */
  double *g; /* the energy gradient and the prameter change */
  int nSmat;

  int paraToSmatIdx[NPara], smatToParaIdx[NPara];

  int optNum=0,cutNum=0;
  double sDiag,sDiagMax,sDiagMin;
  double *sDiagElm;

  double diagCutThreshold;

  FILE *fp;
  double rmax;
  int simax;

  int si,sj; /* index for matrix S */
  int pi,pj; /* index for variational parameters */

  int info=0;
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  if(rank!=0) {
    MPI_Bcast(&info, 1, MPI_INT, 0, comm);
    return info;
  }

  RequestWorkSpaceDouble(SROptSize*(SROptSize+1));
  sDiagElm = GetWorkSpaceDouble(NPara);

  StartTimer(50);

  for(pi=0;pi<NPara;pi++) {
    /* sDiagElm is temporarily used for diagonal elements of S */
    /* S[i][i] = OO[pi+1][pi+1] - OO[0][pi+1] * OO[0][pi+1]; */
    sDiagElm[pi] = SROptOO[(pi+1)*SROptSize+pi+1] - SROptOO[pi+1] * SROptOO[pi+1];
  }

  sDiag = sDiagElm[0];
  sDiagMax=sDiag; sDiagMin=sDiag;
  for(pi=0;pi<NPara;pi++) {
    sDiag = sDiagElm[pi];
    if(sDiag>sDiagMax) sDiagMax=sDiag;
    if(sDiag<sDiagMin) sDiagMin=sDiag;
  }
  diagCutThreshold = sDiagMax*DSROptRedCut;

  si = 0;
  for(pi=0;pi<NPara;pi++) {
    if(OptFlag[pi]!=1) { /* fixed by OptFlag */
      paraToSmatIdx[pi] = -1;
      optNum++;
      continue;
    }

    sDiag = sDiagElm[pi];
    if(sDiag < diagCutThreshold) { /* fixed by diagCut */
      paraToSmatIdx[pi] = -1;
      cutNum++;
    } else { /* optimized */
      paraToSmatIdx[pi] = si;
      smatToParaIdx[si] = pi;
      si += 1;
    }
  }
  nSmat = si;
  for(si=nSmat;si<NPara;si++) {
    smatToParaIdx[si] = -1;
  }

  ReleaseWorkSpaceDouble();
  s = GetWorkSpaceDouble(nSmat*nSmat);
  g = GetWorkSpaceDouble(nSmat);

  StopTimer(50);
  StartTimer(51);

  StartTimer(56);
  stcOptInit(s,g,nSmat,smatToParaIdx);
  StopTimer(56);
  StartTimer(57);
  info = stcOptMain(s,g,nSmat);
  /* g is overwritten with the parameter change. */
  StopTimer(57);

  StopTimer(51);
  StartTimer(52);

  if(info!=0) fprintf(stderr, "StcOpt: DPOSV info=%d\n",info);

  /*** print zqp_SRinfo.dat ***/
  rmax = g[0]; simax=0;
  for(si=0;si<nSmat;si++) {
    if(fabs(rmax) < fabs(g[si])) {
      rmax = g[si]; simax=si;
    }
  }
  fprintf(FileSRinfo, "%5d %5d %5d %5d % .5e % .5e % .5e %5d\n",NPara,nSmat,optNum,cutNum,
          sDiagMax,sDiagMin,rmax,smatToParaIdx[simax]);

  /*** check inf and nan ***/
  for(si=0;si<nSmat;si++) {
    if( !isfinite(g[si]) ) {
      fprintf(stderr, "StcOpt: r[%d]=%.10lf\n",si,g[si]);
      info = 1;
      break;
    }
  }

  /* update variatonal parameters */
  if(info==0) {
    for(si=0;si<nSmat;si++) {
      pi = smatToParaIdx[si];
      Para[pi] += g[si];
    }
  }

  MPI_Bcast(&info, 1, MPI_INT, 0, comm);

  StopTimer(52);

  ReleaseWorkSpaceDouble();
  return info;
}

/* calculate the parameter change from s and g by DPOSV */
/* g is overwritten with the parameter change. */
int stcOptMain(double *const s, double *const g, const int nSmat) {
  /* for DPOSV */
  char uplo;
  int n,nrhs,lds,ldg,info;
  int si;
  StartTimer(53);

  uplo='U'; n=nSmat; nrhs=1; lds=n; ldg=n;
  M_DPOSV(&uplo, &n, &nrhs, s, &lds, g, &ldg, &info);

  StopTimer(53);
  return info;
}


void stcOptInit(double *const s, double *const g, const int nSmat, int *const smatToParaIdx) {
  const double ratioDiag = 1.0 + DSROptStaDel;
  int si,sj,pi,pj,idx,offset;
  double tmp;
  
  /* calculate the overlap matrix S */
  /* S[i][j] = OO[i+1][j+1] - OO[0][i+1] * OO[0][j+1]; */
  for(si=0;si<nSmat;++si) {
    pi = smatToParaIdx[si];
    offset = (pi+1)*SROptSize;
    tmp = SROptOO[pi+1];

    for(sj=0;sj<nSmat;++sj) {
      pj = smatToParaIdx[sj];
      idx = si + nSmat*sj; /* column major */
      s[idx] = SROptOO[offset+pj+1] - tmp * SROptOO[pj+1];
    }

    /* modify diagonal elements */
    idx = si + nSmat*si;
    s[idx] *= ratioDiag;
  }

  /* calculate the energy gradient * (-dt) */
  /* energy gradient = 2.0*( HO[i+1] - HO[0] * OO[i+1]) */
  for(si=0;si<nSmat;++si) {
    pi = smatToParaIdx[si];
    g[si] = -DSROptStepDt*2.0*(SROptHO[pi+1] - SROptHO[0] * SROptOO[pi+1]);
  }

  return;
}
