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
 * Stochastic Reconfiguration method
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

// #define _DEBUG_STCOPT

#include "stcopt.h"

int StochasticOpt(MPI_Comm comm) {
  const int nPara=NPara;
  const int srOptSize=SROptSize;
  const double complex *srOptOO=SROptOO;
  //const double         *srOptOO= SROptOO_real;

  double *r; /* the parameter change */
  int nSmat;
  int smatToParaIdx[2*NPara];//TBC

  int cutNum=0,optNum=0;
  double sDiag,sDiagMax,sDiagMin;
  double diagCutThreshold;

  int si; /* index for matrix S */
  int pi; /* index for variational parameters */

  double rmax;
  int simax;
  int info=0;

// for real
  int int_x,int_y,j,i;

  double complex *para=Para;

  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  r = (double*)calloc(2*SROptSize, sizeof(double));

  StartTimer(50);
//[s] for only real variables TBC
  if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){ //real &  sz=0
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
  #pragma omp parallel for default(shared) private(pi)
  #pragma loop noalias
  for(pi=0;pi<2*nPara;pi++) {
    //for(pi=0;pi<nPara;pi++) {
    /* r[i] is temporarily used for diagonal elements of S */
    /* S[i][i] = OO[pi+1][pi+1] - OO[0][pi+1] * OO[0][pi+1]; */
    r[pi]   = creal(srOptOO[(pi+2)*(2*srOptSize)+(pi+2)]) - creal(srOptOO[pi+2]) * creal(srOptOO[pi+2]);
    //printf("DEBUG: pi=%d: %lf %lf \n",pi,creal(srOptOO[pi]),cimag(srOptOO[pi]));
#ifdef _DEBUG_STCOPT
  fprintf(stderr, "DEBUG in %s (%d): r[%d] = %lf\n", __FILE__, __LINE__, pi, r[pi]);
#endif
  }

// search for max and min
  sDiag = r[0];
  sDiagMax=sDiag; sDiagMin=sDiag;
  for(pi=0;pi<2*nPara;pi++) {
    sDiag = r[pi];
    if(sDiag>sDiagMax) sDiagMax=sDiag;
    if(sDiag<sDiagMin) sDiagMin=sDiag;
  }

// threshold
// optNum = number of parameters 
// cutNum: number of paramers that are cut
  diagCutThreshold = sDiagMax*DSROptRedCut;
  si = 0;
  for(pi=0;pi<2*nPara;pi++) {
    //printf("DEBUG: nPara=%d pi=%d OptFlag=%d r=%lf\n",nPara,pi,OptFlag[pi],r[pi]);
    if(OptFlag[pi]!=1) { /* fixed by OptFlag */
      optNum++;
      continue; //skip sDiag
    }
// s:this part will be skipped if OptFlag[pi]!=1
    sDiag = r[pi];
    if(sDiag < diagCutThreshold) { /* fixed by diagCut */
      cutNum++;
    } else { /* optimized */
      smatToParaIdx[si] = pi; // si -> restricted parameters , pi -> full paramer 0 <-> 2*NPara
      si += 1;
    }
// e
  }
  nSmat = si;
  for(si=nSmat;si<2*nPara;si++) {
    smatToParaIdx[si] = -1; // parameters that will not be optimized
  }

#ifdef _DEBUG_STCOPT
  printf("DEBUG in %s (%d): diagCutThreshold = %lg\n", __FILE__, __LINE__, diagCutThreshold);
  printf("DEBUG in %s (%d): optNum, cutNum, nSmat, 2*nPara == %d, %d, %d, %d\n", __FILE__, __LINE__, optNum, cutNum, nSmat, 2*nPara);
#endif

  StopTimer(50);
  StartTimer(51);


  //printf("DEBUG: nSmat=%d \n",nSmat);
  /* calculate r[i]: global vector [nSmat] */
  info = stcOptMain(r, nSmat, smatToParaIdx, comm);

  StopTimer(51);
  StartTimer(52);

  /*** print zqp_SRinfo.dat ***/
  if(rank==0) {
    if(info!=0) fprintf(stderr, "StcOpt: DPOSV info=%d\n",info);
    rmax = r[0]; simax=0;;
    for(si=0;si<nSmat;si++) {
      if(fabs(rmax) < fabs(r[si])) {
        rmax = r[si]; simax=si;
      }
    }

    fprintf(FileSRinfo, "%5d %5d %5d %5d % .5e % .5e % .5e %5d\n",NPara,nSmat,optNum,cutNum,
            sDiagMax,sDiagMin,rmax,smatToParaIdx[simax]);
  }

  /*** check inf and nan ***/
  if(rank==0) {
    for(si=0;si<nSmat;si++) {
      if( !isfinite(r[si]) ) {
        fprintf(stderr, "StcOpt: r[%d]=%.10lf\n",si,r[si]);
        info = 1;
        break;
      }
    }
  }
  MPI_Bcast(&info, 1, MPI_INT, 0, comm);

  /* update variational parameters */
  if(info==0 && rank==0) {
    #pragma omp parallel for default(shared) private(si,pi)
    #pragma loop noalias
    #pragma loop norecurrence para
    for(si=0;si<nSmat;si++) {
      pi = smatToParaIdx[si];
      if(pi%2==0){
        para[pi/2]     += r[si];  // real
      }else{
        para[(pi-1)/2] += r[si]*I; // imag
      }
    }
  }

  free(r);

  StopTimer(52);
  return info;
}
