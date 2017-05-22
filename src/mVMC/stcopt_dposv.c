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
#include "stcopt_dposv.h"
#ifndef _SRC_STCOPT_DPOSV
#define _SRC_STCOPT_DPOSV

// #define _DEBUG_STCOPT_DPOSV

int StochasticOpt(MPI_Comm comm) {
  double *s; /* the overlap matrix S */
  double *g; /* the energy gradient and the prameter change */
  int nSmat;

  int smatToParaIdx[2*NPara];

  int optNum=0,cutNum=0;
  double sDiag,sDiagMax,sDiagMin;
  double *sDiagElm;

  double diagCutThreshold;

  double rmax;
  int simax;

  int si; /* index for matrix S */
  int pi; /* index for variational parameters */

  int info=0;
  int rank,size;
// for real
  int int_x,int_y,j,i;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  if(rank!=0) {
    MPI_Bcast(&info, 1, MPI_INT, 0, comm);
    return info;
  }
//[s] for only real variables TBC
  if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){//real & sz=0
    #pragma omp parallel for default(shared) private(i,int_x,int_y,j)
    #pragma loop noalias
    for(i=0;i<2*SROptSize*(2*SROptSize+2);i++){
      int_x  = i%(2*SROptSize);
      int_y  = (i-int_x)/(2*SROptSize);
      if(int_x%2==0 && int_y%2==0){
        j          = int_x/2+(int_y/2)*SROptSize;
        SROptOO[i] = SROptOO_real[j];// only real part TBC
      }else{
        SROptOO[i] = 0.0+0.0*I;
      }
    }
  }
//[e]


  RequestWorkSpaceDouble(2*SROptSize*(2*SROptSize+2));
  sDiagElm = GetWorkSpaceDouble(2*NPara); //NPara -> 2*NPars

  StartTimer(50);

  for(pi=0;pi<2*NPara;pi++) {
    /* sDiagElm is temporarily used for diagonal elements of S */
    /* S[i][i] = OO[pi+1][pi+1] - OO[0][pi+1] * OO[0][pi+1]; */
    sDiagElm[pi]  = creal(SROptOO[(pi+2)*(2*SROptSize)+(pi+2)]) - creal(SROptOO[pi+2]) * creal(SROptOO[pi+2]);

#ifdef _DEBUG_STCOPT_DPOSV
    fprintf(stderr, "DEBUG in %s (%d): sDiagElm[%d] = %lf\n", __FILE__, __LINE__, pi, sDiagElm[pi]);
#endif
  }

  sDiag = sDiagElm[0];
  sDiagMax=sDiag; sDiagMin=sDiag;
  for(pi=0;pi<2*NPara;pi++) {
    sDiag = sDiagElm[pi];
    if(sDiag>sDiagMax) sDiagMax=sDiag;
    if(sDiag<sDiagMin) sDiagMin=sDiag;
  }
  diagCutThreshold = sDiagMax*DSROptRedCut;

  si = 0;
  for(pi=0;pi<2*NPara;pi++) {
    if(OptFlag[pi]!=1) { /* fixed by OptFlag */
//      paraToSmatIdx[pi] = -1;
      optNum++;
      continue;
    }

    sDiag = sDiagElm[pi];
    if(sDiag < diagCutThreshold) { /* fixed by diagCut */
  //    paraToSmatIdx[pi] = -1;
      cutNum++;
    } else { /* optimized */
    //  paraToSmatIdx[pi] = si;
      smatToParaIdx[si] = pi;
      si += 1;
    }
  }
  nSmat = si;
  for(si=nSmat;si<2*NPara;si++) {
    smatToParaIdx[si] = -1;
  }

#ifdef _DEBUG_STCOPT_DPOSV
  printf("DEBUG in %s (%d): diagCutThreshold = %lg\n", __FILE__, __LINE__, diagCutThreshold);
  printf("DEBUG in %s (%d): optNum, cutNum, nSmat, 2*NPara == %d, %d, %d, %d\n", __FILE__, __LINE__, optNum, cutNum, nSmat, 2*NPara);
#endif

  ReleaseWorkSpaceDouble();
  s = GetWorkSpaceDouble(nSmat*nSmat);
  g = GetWorkSpaceDouble(nSmat);

  StopTimer(50);
  StartTimer(51);

  StartTimer(56);
  stcOptInit(s,g,nSmat,smatToParaIdx);
  StopTimer(56);

#ifdef _DEBUG_STCOPT_DPOSV
  for(i=0; i<nSmat; ++i){
    printf("%lg\n", g[i]);
  }
  for(j=0;j<nSmat;j++) {
    for(i=0;i<nSmat;i++) {
      idx = i + j*nSmat; /* local index (row major) */
      printf("%lg ", s[idx]);
    }
    printf("\n");
  }
#endif

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

#ifdef _DEBUG_STCOPT_DPOSV
  for(si=0; si<nSmat; ++si){
    fprintf(stderr, "%lg\n", g[si]);
  }
#endif

  /*** check inf and nan ***/
  for(si=0;si<nSmat;si++) {
    if( !isfinite(g[si]) ) {
      fprintf(stderr, "StcOpt: r[%d]=%.10lf\n",si,g[si]);
      info = 1;
      break;
    }
  }

  /* update variational parameters */
  if(info==0) {
    for(si=0;si<nSmat;si++) {
      pi = smatToParaIdx[si];
      //printf("DEBUG: nSmat=%d: si=%d pi=%d: %lf\n",nSmat,si,pi,g[si]);
      if(pi%2==0){
        Para[pi/2]     += g[si];  // real
        //printf("Real: DEBUG: nSmat=%d: si=%d pi=%d: %lf %lf\n",nSmat,si,pi,creal(Para[pi/2]),cimag(Para[pi/2]));
      }else{
        Para[(pi-1)/2] += g[si]*I; // imag
        //printf("Imag: DEBUG: nSmat=%d: si=%d pi=%d: %lf %lf\n",nSmat,si,pi,creal(Para[(pi-1)/2]),cimag(Para[(pi-1)/2]));
      }
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
    //offset = (pi+1)*SROptSize;
    offset = (pi+2)*(2*SROptSize);
    tmp = creal(SROptOO[pi+2]);

    for(sj=0;sj<nSmat;++sj) {
      pj = smatToParaIdx[sj];
      idx = si + nSmat*sj; /* column major */
      //s[idx] = SROptOO[offset+pj+1] - tmp * SROptOO[pj+1];
      s[idx] = creal(SROptOO[offset+(pj+2)]) - tmp * creal(SROptOO[pj+2]);
    }

    /* modify diagonal elements */
    idx = si + nSmat*si;
    s[idx] *= ratioDiag;
  }

  /* calculate the energy gradient * (-dt) */
  /* energy gradient = 2.0*( HO[i+1] - HO[0] * OO[i+1]) */
  for(si=0;si<nSmat;++si) {
    pi = smatToParaIdx[si];
    g[si] = -DSROptStepDt*2.0*(creal(SROptHO[pi+2]) - creal(SROptHO[0]) * creal(SROptOO[pi+2]));
  }

  return;
}

#endif
