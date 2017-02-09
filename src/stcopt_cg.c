/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

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

// #define _DEBUG_STCOPT_CG

#ifdef _SYSTEM_A
 #define M_PDGEMV  PDGEMV
 #define M_PZSYEVD PZSYEVD
 #define M_PZGEMV  PZGEMV
 #define M_PZPOSV  PZPOSV
 #define M_ZPOSV  ZPOSV
 #define M_ZGETRS  ZGETRS
 #define M_ZGETRF  ZGETRF
 #define M_PDGEMV  PDGEMV
 #define M_PZHEEVD  PZHEEVD
 #define M_DGEMV  DGEMV
#elif _lapack_small_nounderscore
 #define M_PDGEMV  pdgemv
 #define M_PZSYEVD pzsyevd
 #define M_PZGEMV  pzgemv
 #define M_ZPOSV  zposv
 #define M_ZGETRS  zgetrs
 #define M_ZGETRF  zgetrf
 #define M_PDGEMV  pdgemv
 #define M_PZHEEVD  pzheevd
 #define M_DGEMV  dgemv
#else
 #define M_NUMROC   numroc_
 #define M_DESCINIT descinit_
 #define M_PDPOSV  pdposv_
 #define M_PDGEMV  pdgemv_
 #define M_PZSYEVD  pzsyevd_
 #define M_PZPOSV  pzposv_
 #define M_PZGEMV  pzgemv_
 #define M_ZPOSV  zposv_
 #define M_ZGETRS  zgetrs_
 #define M_ZGETRF  zgetrf_
 #define M_DSYEV  dsyev_
 #define M_PDGEMV  pdgemv_
 #define M_PZHEEVD  pzheevd_
 #define M_DGEMV  dgemv_
#endif

int Csys2blacs_handle(MPI_Comm comm);
int Cblacs_gridinit(int *ictxt, char *order, int nprow, int npcol);
int Cblacs_gridexit(int ictxt);
int Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *myprow, int *mypcol);
int M_NUMROC(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
int M_DESCINIT(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc,
               int *icsrc, int *ictxt, int *lld, int *info);

void M_PDGEMV(char *trans, int *m, int *n, double *alpha,
              double *a, int *ia, int *ja, int *desca,
              double *x, int *ix, int *jx, int *descx, int *incx,
              double *beta,
              double *y, int *iy, int *jy, int *descy, int *incy);
void M_DGEMV(char *trans, int *m, int *n, double *alpha,
              double *a, int *lda, double *x, int *incx, double *beta,
              double *y, int *incy);


int M_PDPOSV(char *uplo, int *n, int *nrhs, double *a, int *ia, int *ja, int *desca,
             double *b, int *ib, int *jb, int *descb, int *info);
int M_PZPOSV(char *uplo, int *n, int *nrhs, double complex *a, int *ia, int *ja, int *desca,
             double complex *b, int *ib, int *jb, int *descb, int *info);
void M_PDGEMV(char *trans, int *m, int *n, double *alpha,
              double *a, int *ia, int *ja, int *desca,
              double *x, int *ix, int *jx, int *descx, int *incx,
              double *beta,
              double *y, int *iy, int *jy, int *descy, int *incy);
void M_PZGEMV(char *trans, int *m, int *n, double complex *alpha,
              double complex *a, int *ia, int *ja, int *desca,
              double complex *x, int *ix, int *jx, int *descx, int *incx,
              double complex *beta,
              double complex *y, int *iy, int *jy, int *descy, int *incy);

int StochasticOptCG(MPI_Comm comm);
void stcOptCG_Init(const int nSmat, int *const smatToParaIdx, double *VecCG); 
int stcOptCG_Main(const int nSmat, double *VecCG, MPI_Comm comm);
inline double xdot(const int n, double * const p, double * const q);

int StochasticOptCG(MPI_Comm comm) {
  const int nPara=NPara;
  const int srOptSize=SROptSize;
  const double complex *srOptO=SROptOO;
  const double complex *srOptOOdiag=SROptOO + 2*srOptSize;

  double sDiagElm[2*nPara]; /* the parameter change */
  int nSmat;
  int smatToParaIdx[2*nPara];

  int cutNum=0,optNum=0;
  double sDiag,sDiagMax,sDiagMin;
  double diagCutThreshold;

  int si; /* index for matrix S */
  int pi; /* index for variational parameters */

  double rmax;
  int simax;
  int info=0;

  complex double *para=Para;
 //double VecCG[SROptSize*(NVMCSample+7) + NVMCSample];
  double *VecCG;
  //double *stcOd, *stcOs, *stcO;
  //double  *r, *g, *y, *d, *x, *process_y;
  double *r;


  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): Start StochasticOptCG\n", __FILE__, __LINE__);
  fprintf(stderr, "DEBUG in %s (%d): nPara = %d, srOptSize = %d\n", __FILE__, __LINE__, nPara, srOptSize);
#endif

  StartTimer(50);

  if(AllComplexFlag==0){
    for(pi=0; pi<srOptSize; pi++){
      SROptOO[2*pi] = SROptOO_real[pi];
      SROptOO[2*pi+1] = 0.0;
      SROptOO[2*pi+2*srOptSize] = SROptOO_real[pi+srOptSize];
      SROptOO[2*pi+1+2*srOptSize] = 0.0;
      SROptHO[2*pi] = SROptHO_real[pi];
      SROptHO[2*pi+1] = 0.0;
    }
  }

  #pragma omp parallel for default(shared) private(pi)
  #pragma loop noalias
  for(pi=0;pi<2*nPara;pi++) {
    /* calculate diagonal elements of S */
    /* S[pi][pi] = OO[pi+2][pi+2] - O[pi+2] * O[pi+2]; */
    sDiagElm[pi] = creal(srOptOOdiag[pi+2] - conj(srOptO[pi+2]) * srOptO[pi+2]);
#ifdef _DEBUG_STCOPT_CG
    fprintf(stderr, "DEBUG in %s (%d): sDiagElm[%d] = %lf (optflag=%d)\n", __FILE__, __LINE__, pi, sDiagElm[pi], OptFlag[pi]);
#endif
  }

  sDiag = sDiagElm[0];
  sDiagMax=sDiag; sDiagMin=sDiag;
  for(pi=0;pi<2*nPara;pi++) {
    sDiag = sDiagElm[pi];
    if(sDiag>sDiagMax) sDiagMax=sDiag;
    if(sDiag<sDiagMin) sDiagMin=sDiag;
  }

  diagCutThreshold = sDiagMax*DSROptRedCut;
  si = 0;
  for(pi=0;pi<2*nPara;pi++) {
    if(OptFlag[pi]!=1) { /* fixed by OptFlag */
      optNum++;
      continue;
    }

    sDiag = sDiagElm[pi];
    if(sDiag < diagCutThreshold) { /* fixed by diagCut */
      cutNum++;
    } else { /* optimized */
      smatToParaIdx[si] = pi;
      si += 1;
    }
  }
  nSmat = si;
  for(si=nSmat;si<2*nPara;si++) {
    smatToParaIdx[si] = -1;
  }

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): diagCutThreshold = %lf\n", __FILE__, __LINE__, diagCutThreshold);
  fprintf(stderr, "DEBUG in %s (%d): optNum, cutNum, nSmat, nPara == %d, %d, %d, %d\n", __FILE__, __LINE__, optNum, cutNum, nSmat, nPara);
#endif

  StopTimer(50);
  StartTimer(51);
  
  RequestWorkSpaceDouble(2*SROptSize*(NVMCSample+7) + NVMCSample);
  VecCG  = GetWorkSpaceDouble(2*SROptSize*(NVMCSample+7) + NVMCSample);
  
  /* calculate r[i]: global vector [nSmat] */
  stcOptCG_Init(nSmat, smatToParaIdx, VecCG);
  info = stcOptCG_Main(nSmat, VecCG, comm);
  
  StopTimer(51);
  StartTimer(52);

  r = VecCG;
  /*** print zqp_SRinfo.dat ***/
  if(rank==0) {
    //if(info!=0) fprintf(stderr, "StcOpt: DPOSV info=%d\n",info);
    rmax = r[0]; simax=0;;
    for(si=0;si<nSmat;si++) {
      if(fabs(rmax) < fabs(r[si])) {
        rmax = r[si]; simax=si;
      }
    }

    fprintf(FileSRinfo, "%5d %5d %5d %5d % .5e % .5e % .5e %5d, %d\n",NPara,nSmat,optNum,cutNum,
            sDiagMax,sDiagMin,rmax,smatToParaIdx[simax], info);
    //fprintf(FileSRinfo, "%5d %5d %5d %5d % .5e %5d, %d\n",NPara,nSmat,optNum,cutNum,
    //        rmax,smatToParaIdx[simax], info);
  }
  info=0;

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

  /* update variatonal parameters */
  if(info==0 && rank==0) {
    #pragma omp parallel for default(shared) private(si,pi)
    #pragma loop noalias
    #pragma loop norecurrence para
    for(si=0;si<nSmat;si++) {
      pi = smatToParaIdx[si];
      para[pi] += r[si];
    }
  }

  StopTimer(52);
  ReleaseWorkSpaceDouble();


#ifdef _DEBUG_STCOPT_CG
  printf("DEBUG in %s (%d): End StochasticOptCG\n", __FILE__, __LINE__);
#endif
  return info;
}

/* calculate the parameter change r[nSmat] from SOpt.
   Solve S*x = g */
int stcOptCG_Main(const int nSmat, double *VecCG, MPI_Comm comm) {
  int si,pi,pj,idx;
  int rank, size, info;
  int iter;
  int max_iter=NSROptCGMaxIter;
  double delta;
  double alpha, beta;
  double cg_thresh = DSROptCGTol*(double)NPara;
  //double cg_thresh = DSROptRedCut;

  const int srOptSize=SROptSize;
  double *stcOd, *stcOs, *stcO;
  double  *r, *g, *q, *d, *x, *process_y;

#ifdef _DEBUG_STCOPT_CG
  printf("DEBUG in %s (%d): Start stcOptCG_Main\n", __FILE__, __LINE__);
#endif
  
  x = VecCG;              //GetWorkSpaceDouble(nSmat)
  g = x  + nSmat;         //GetWorkSpaceDouble(nSmat)
  q = g  + nSmat;         //GetWorkSpaceDouble(nSmat)
  d = q  + nSmat;         //GetWorkSpaceDouble(nSmat)
  r = d  + nSmat;         //GetWorkSpaceDouble(nSmat)
  process_y = r  + nSmat; //GetWorkSpaceDouble(nSmat)
  
  stcO = process_y + nSmat;                     //GetWorkSpaceDouble(nSmat)
  stcOs = stcO  + nSmat;            //GetWorkSpaceDouble(nSmat*NVMCSample)
  stcOd = stcOs + nSmat*NVMCSample; //GetWorkSpaceDouble(NVMCSample)
  
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  operate_by_S(nSmat, x, r, stcO, comm);
  for(si=0;si<nSmat;++si) {
    r[si] = g[si] - r[si];
    d[si] = r[si];
  }

  delta = xdot(nSmat, r, r);

  for(iter=0; iter < max_iter; iter++){
    //check convergence 
    //if(rank==0) printf("%d: %d %.3e\n",rank, iter,delta);
    if (fabs(sqrt(delta)) < cg_thresh) break;

    // compute vector q=S*d
    operate_by_S(nSmat, d, q, stcO, comm);
    alpha = delta/xdot(nSmat,d,q);
  
    // update solution vector x=x+alpha*d
    for(si=0;si<nSmat;++si) {
      x[si] = x[si] + alpha*d[si];
    }
    // update residual vector r=r-alpha*q, q=S*d
    if((iter+1) % 20 == 0){
      operate_by_S(nSmat, x, r, stcO, comm);
      for(si=0;si<nSmat;++si) {
        r[si] = g[si] - r[si];
      }
    }else{
      for(si=0;si<nSmat;++si) {
        r[si] = r[si] - alpha*q[si];
      }
    }
    beta = xdot(nSmat,r,r)/delta;

    //update the norm of residual vector r
    delta = beta*delta;
    // update direction vector d
    for(si=0;si<nSmat;++si) {
      d[si] = r[si] + beta*d[si];
    }
  }

#ifdef _DEBUG_STCOPT_CG
  printf("DEBUG in %s (%d): End stcOptCG_Main\n", __FILE__, __LINE__);
#endif

  return iter;
}

/* calculate  y = S*x */
/* S is the overlap matrix*/
/* S[i][j] = OO[i+1][j+1] - OO[i+1][0] * OO[0][j+1]; */
int operate_by_S(int nSmat, double *x, double *y, double * stcO, MPI_Comm comm) {
  double process_y[nSmat];
  double *stcOd, *stcOs;
  int rank,size,info=0;
  int i,si;
  int incx = 1, incy = 1;
  double a, coef;
  char transT='T', transN='N';
  double alpha = 1.0, beta = 0.0;
  double invW = 1.0/Wc;
  
  stcOs = stcO  + nSmat;            //GetWorkSpaceDouble(nSmat*NVMCSample)
  stcOd = stcOs + nSmat*NVMCSample; //GetWorkSpaceDouble(NVMCSample)
  
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Bcast(x, nSmat, MPI_DOUBLE, 0, comm);
  
  StartTimer(53);
  for(si=0;si<nSmat;++si) {
    y[si] = 0.0;
    process_y[si] = 0.0;
  }
  for(i=0;i<NVMCSample;++i) {
    stcOd[i] = 0.0;
  }

  /* Od[sample] = sum{sidx} x[sidx] * O[sidx][sample] ; */
  M_DGEMV(&transT, &nSmat, &NVMCSample, &alpha, stcOs, &nSmat, x, &incx, &beta, stcOd, &incy);
  /* y[rank][sidx] = sum{sample} O[sidx][sample] * Od[sample] ; */
  M_DGEMV(&transN, &nSmat, &NVMCSample, &alpha, stcOs, &nSmat, stcOd, &incx, &beta, process_y, &incy);
  
  MPI_Barrier(comm);
  /* compute <OO>*x */
  SafeMpiAllReduce(process_y, y, nSmat, comm);
 
  /* compute <O>*x */
  coef = xdot(nSmat, stcO, x);
  
  /* y = S*x */
  for(si=0;si<nSmat;++si) {
    y[si] = invW*y[si] - coef*stcO[si];
    /* modify diagonal elements */
    y[si] += DSROptStaDel * x[si];
  }

  StopTimer(53);
  return info;
}

void stcOptCG_Init(const int nSmat, int *const smatToParaIdx, double *VecCG) {
  const double ratioDiag = 1.0 + DSROptStaDel;
  const double dt = 2.0*DSROptStepDt;
  const double srOptHO_0 = creal(SROptHO[0]);
  int i;
  int si,sj,pi,pj,idx,offset;
  double *stcOd, *stcOs, *stcO;
  double  *r, *g, *y, *d, *x, *process_y;

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): Start stcOptCG_Init\n", __FILE__, __LINE__);
#endif
  
  x = VecCG;              //GetWorkSpaceDouble(nSmat)
  g = x  + nSmat;         //GetWorkSpaceDouble(nSmat)
  y = g  + nSmat;         //GetWorkSpaceDouble(nSmat)
  d = y  + nSmat;         //GetWorkSpaceDouble(nSmat)
  r = d  + nSmat;         //GetWorkSpaceDouble(nSmat)
  process_y = r  + nSmat; //GetWorkSpaceDouble(nSmat)
  
  stcO = process_y + nSmat;         //GetWorkSpaceDouble(nSmat)
  stcOs = stcO  + nSmat;            //GetWorkSpaceDouble(nSmat*NVMCSample)
  stcOd = stcOs + nSmat*NVMCSample; //GetWorkSpaceDouble(NVMCSample)
  
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<2*SROptSize*(NVMCSample+7) + NVMCSample;si++) {
    VecCG[si] = 0.0;
  }

  for(i=0;i<NVMCSample;++i) {
    offset = i*SROptSize;
    for(si=0;si<nSmat;++si) {
      pi = smatToParaIdx[si]/2;
      
      idx = si + i*nSmat; /* column major */
      stcOs[idx] = SROptO_Store_real[offset+pi+1];
    }
  }

  /* calculate the energy gradient * (-dt) */
  /* energy gradient = 2.0*( HO[i+1] - H * O[i+1]) */
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    pi = smatToParaIdx[si];
    
    stcO[si] = SROptO[pi+2];
    g[si] = -dt*(creal(SROptHO[pi+2]) - srOptHO_0 * creal(SROptO[pi+2]));
  }

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): Start stcOptCG_Init\n", __FILE__, __LINE__);
#endif

  return;
}

inline double xdot(const int n, double * const p, double * const q) {
  int i;
  double z=0;
  #pragma loop noalias
  for(i=0;i<n;i++) {
    z += p[i]*q[i];
  }
  return z;
}
