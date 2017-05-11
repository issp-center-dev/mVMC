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

#ifdef MVMC_SRCG_REAL
  #define fn_StochasticOptCG StochasticOptCG_real
  #define fn_StochasticOptCG_Init StochasticOptCG_Init_real
  #define fn_StochasticOptCG_Main StochasticOptCG_Main_real
  #define fn_operate_by_S operate_by_S_real
  #define fn_xdot xdot_real

  #define elemtype double
  #define CONJ(x) (x)
  #define CREAL(x) (x)
  #define CIMAG(x) (0.0)
  #define Gemv M_DGEMV
  #define Gemm M_DGEMM
  #define Ger M_DGER
  #define allreduce SafeMpiAllReduce
  #define request_workspace RequestWorkSpaceDouble
  #define get_workspace GetWorkSpaceDouble
  #define release_workspace ReleaseWorkSpaceDouble

  #define OFFSET 1
#else // MVMC_SRCG_REAL
  #define fn_StochasticOptCG StochasticOptCG_fcmp
  #define fn_StochasticOptCG_Init StochasticOptCG_Init_fcmp
  #define fn_StochasticOptCG_Main StochasticOptCG_Main_fcmp
  #define fn_operate_by_S operate_by_S_fcmp
  #define fn_xdot xdot_fcmp

  #define Gemv M_ZGEMV
  #define Gemm M_ZGEMM
  #define Ger M_ZGERC
  #define allreduce SafeMpiAllReduce_fcmp
  #define request_workspace RequestWorkSpaceComplex
  #define get_workspace GetWorkSpaceComplex
  #define release_workspace ReleaseWorkSpaceComplex

  #define elemtype double complex
  #define CONJ(x) conj(x)
  #define CREAL(x) creal(x)
  #define CIMAG(x) cimag(x)

  #define OFFSET 2
#endif

#define print_val(x) fprintf(stderr, "%lg ", CREAL(x))
#define println_val(x) fprintf(stderr, "%lg\n", CREAL(x))


int fn_StochasticOptCG(MPI_Comm comm);
void fn_StochasticOptCG_Init(const int nSmat, int *const smatToParaIdx, elemtype *VecCG); 
int fn_StochasticOptCG_Main(const int nSmat, elemtype *VecCG, MPI_Comm comm);
int fn_operate_by_S(const int nSmat, elemtype *x, elemtype *y, elemtype * process_y, elemtype * stcO, MPI_Comm comm);
inline elemtype fn_xdot(const int n, elemtype * const p, elemtype * const q);

int fn_StochasticOptCG(MPI_Comm comm) {
  const int nPara=OFFSET*NPara;
  const int srOptSize=SROptSize;
#ifdef MVMC_SRCG_REAL
  const elemtype *srOptO=SROptOO_real;
  const elemtype *srOptOOdiag=SROptOO_real + srOptSize;
#else
  const elemtype *srOptO=SROptOO;
  const elemtype *srOptOOdiag=SROptOO + 2*srOptSize;
#endif

  double sDiagElm[nPara]; /* the parameter change */
  int nSmat;
  int smatToParaIdx[nPara];

  int cutNum=0,optNum=0;
  double sDiag,sDiagMax,sDiagMin;
  double diagCutThreshold;

  int si; /* index for matrix S */
  int pi; /* index for variational parameters */

  elemtype rmax;
  int simax;
  int info=0;

  double complex *para=Para;
 //double VecCG[SROptSize*(NVMCSample+7) + NVMCSample];
  elemtype *VecCG;
  //double *stcOd, *stcOs, *stcO;
  //double  *r, *g, *y, *d, *x, *process_y;
  elemtype *r;


  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): Start StochasticOptCG\n", __FILE__, __LINE__);
  fprintf(stderr, "DEBUG in %s (%d): nPara = %d, srOptSize = %d\n", __FILE__, __LINE__, nPara, srOptSize);
#endif

  StartTimer(50);

  #pragma omp parallel for default(shared) private(pi)
  #pragma loop noalias
  for(pi=0;pi<nPara;pi++) {
    /* calculate diagonal elements of S */
    /* S[pi][pi] = OO[pi+offset][pi+offset] - O[pi+offset] * O[pi+offset]; */
    sDiagElm[pi] = CREAL(srOptOOdiag[pi+OFFSET]) - CREAL(srOptO[pi+OFFSET]) * CREAL(srOptO[pi+OFFSET]);
  }
#ifdef _DEBUG_STCOPT_CG
  for(pi=0;pi<nPara;pi++) {
    fprintf(stderr, "DEBUG in %s (%d): sDiagElm[%d] = %lf\n", __FILE__, __LINE__, pi, sDiagElm[pi]);
  }
#endif

  sDiag = sDiagElm[0];
  sDiagMax=sDiag; sDiagMin=sDiag;
  for(pi=0;pi<nPara;pi++) {
    sDiag = sDiagElm[pi];
    if(sDiag>sDiagMax) sDiagMax=sDiag;
    if(sDiag<sDiagMin) sDiagMin=sDiag;
  }

  diagCutThreshold = sDiagMax*DSROptRedCut;
  si = 0;
  for(pi=0;pi<nPara;pi++) {
#ifdef MVMC_SRCG_REAL
    if(OptFlag[2*pi]!=1) { /* fixed by OptFlag */
      optNum++;
      continue;
    }
#else
    if(OptFlag[pi]!=1) { /* fixed by OptFlag */
      optNum++;
      continue;
    }
#endif

    sDiag = sDiagElm[pi];
    if(sDiag < diagCutThreshold) { /* fixed by diagCut */
      cutNum++;
    } else { /* optimized */
      smatToParaIdx[si] = pi;
      si += 1;
    }
  }
  nSmat = si;
  for(si=nSmat;si<nPara;si++) {
    smatToParaIdx[si] = -1;
  }

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): diagCutThreshold = %lg\n", __FILE__, __LINE__, diagCutThreshold);
  fprintf(stderr, "DEBUG in %s (%d): optNum, cutNum, nSmat, nPara == %d, %d, %d, %d\n", __FILE__, __LINE__, optNum, cutNum, nSmat, nPara);
#endif

  StopTimer(50);
  StartTimer(51);
  
  // RequestWorkSpaceDouble(2*SROptSize*(NVMCSample+7) + NVMCSample);
  // VecCG  = GetWorkSpaceDouble(2*SROptSize*(NVMCSample+7) + NVMCSample);

  request_workspace(nSmat*(NVMCSample+8) + NVMCSample);
  VecCG  = get_workspace(nSmat*(NVMCSample+8) + NVMCSample);

  /*
  x = VecCG;              //GetWorkSpaceDouble(nSmat)
  g = x  + nSmat;         //GetWorkSpaceDouble(nSmat)
  q = g  + nSmat;         //GetWorkSpaceDouble(nSmat)
  d = q  + nSmat;         //GetWorkSpaceDouble(nSmat)
  r = d  + nSmat;         //GetWorkSpaceDouble(nSmat)
  process_y = r  + nSmat; //GetWorkSpaceDouble(nSmat)
  stcO = process_y + nSmat;                     //GetWorkSpaceDouble(nSmat)
  stcOs = stcO  + nSmat;            //GetWorkSpaceDouble(nSmat*NVMCSample)
  stcOd = stcOs + nSmat*NVMCSample; //GetWorkSpaceDouble(NVMCSample)
  sdiag = stcOd + NVMCSample;
  */
  
  fn_StochasticOptCG_Init(nSmat, smatToParaIdx, VecCG);


#ifndef _DEBUG_STCOPT_CG_LAPACK
  info = fn_StochasticOptCG_Main(nSmat, VecCG, comm);
#ifdef _DEBUG_STCOPT_CG
  for(si=0; si<nSmat; ++si){
    println_val(VecCG[si]);
  }
#endif
#else // _DEBUG_STCOPT_CG_LAPACK
  {
    /* for POSV */
    char uplo;
    char transN='N', transT='T';
    int i,j,n,nrhs,lds,ldg,info;
    elemtype alpha=1.0, beta=0.0;
    elemtype *s = (elemtype*)calloc(nSmat*nSmat, sizeof(elemtype));
    elemtype *xs= (elemtype*)calloc(nSmat, sizeof(elemtype));
    elemtype *process_y= (elemtype*)calloc(nSmat, sizeof(elemtype));
    elemtype *g = VecCG + nSmat;
    elemtype *sdiag = VecCG + 7*nSmat + nSmat*NVMCSample + NVMCSample;

    n=nSmat; nrhs=1; lds=n; ldg=n;
    for(i=0; i<nSmat; ++i){
      xs[i] = 1.0;
      fn_operate_by_S(nSmat, xs, s+i*nSmat, process_y, VecCG+6*nSmat, comm);
      xs[i] = 0.0;
    }

#ifdef _DEBUG_STCOPT_CG
    for(i=0; i<nSmat; ++i){
      println_val(VecCG[i+nSmat]);
    }
    for(i=0; i<nSmat; ++i){
      for(j=0; j<nSmat; ++j){
        print_val(s[j+i*nSmat]);
      }
      fprintf(stderr, "\n");
    }
#endif

#ifdef MVMC_SRCG_REAL
    uplo='U';
    M_DPOSV(&uplo, &n, &nrhs, s, &lds, g, &ldg, &info);
#else
    uplo='U';
    M_ZPOSV(&uplo, &n, &nrhs, s, &lds, g, &ldg, &info);
#endif
    for(i=0; i<nSmat; ++i){
      VecCG[i] = g[i];
    }

#ifdef _DEBUG_STCOPT_CG
    for(i=0; i<nSmat; ++i){
      println_val(VecCG[i]);
    }
#endif

    free(process_y);
    free(xs);
    free(s);
  }
#endif // _DEBUG_STCOPT_CG_LAPACK
  
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
            sDiagMax,sDiagMin,CREAL(rmax),smatToParaIdx[simax], info);
    //fprintf(FileSRinfo, "%5d %5d %5d %5d % .5e %5d, %d\n",NPara,nSmat,optNum,cutNum,
    //        rmax,smatToParaIdx[simax], info);
  }
  info=0;

  /*** check inf and nan ***/
  if(rank==0) {
    for(si=0;si<nSmat;si++) {
      if( !isfinite(CREAL(r[si])) ) {
        fprintf(stderr, "StcOpt: r[%d]=%.10lf\n",si,CREAL(r[si]));
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
#ifdef MVMC_SRCG_REAL
      para[pi]     += r[si];
#else
      if(pi%2==0){
        para[pi/2]     += creal(r[si]);  // real
      }else{
        para[(pi-1)/2] += creal(r[si])*I; // imag
      }
#endif
    }
  }

#ifdef _DEBUG_STCOPT_CG
  for(si=0; si<nSmat; ++si){
      pi = smatToParaIdx[si];
      fprintf(stderr, "%d %d %lg\n", si, pi, CREAL(r[si]));
  }
#endif

  StopTimer(52);
  release_workspace();

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): End StochasticOptCG\n", __FILE__, __LINE__);
#endif
  return info;
}

/* calculate the parameter change r[nSmat] from SOpt.
   Solve S*x = g */
int fn_StochasticOptCG_Main(const int nSmat, elemtype *VecCG, MPI_Comm comm) {
  int si,pi,pj,idx;
  int rank, size, info;
  int iter;
  int max_iter=NSROptCGMaxIter;
  double delta, beta;
  elemtype alpha;
  double cg_thresh = DSROptCGTol*DSROptCGTol * (double)nSmat * (double)nSmat;
  //double cg_thresh = DSROptRedCut;

  const int srOptSize=SROptSize;
  elemtype *stcOd, *stcOs, *stcO, *sdiag;
  elemtype  *r, *g, *q, *d, *x, *process_y;

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): Start stcOptCG_Main\n", __FILE__, __LINE__);
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
  sdiag = stcOd + NVMCSample;
  
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    d[si] = r[si] = g[si];
  }

  delta = CREAL(fn_xdot(nSmat, r, r));

  for(iter=0; iter < max_iter; iter++){
    //check convergence 
    //if(rank==0) printf("%d: %d %.3e\n",rank, iter,delta);
#ifdef _DEBUG_STCOPT_CG
    fprintf(stderr, "delta = %lg, cg_thresh = %lg\n", delta, cg_thresh);
#endif
    if (delta < cg_thresh) break;

    // compute vector q=S*d
    fn_operate_by_S(nSmat, d, q, process_y, stcO, comm);
    alpha = delta/fn_xdot(nSmat,d,q);
  
    // update solution vector x=x+alpha*d
    #pragma omp parallel for default(shared) private(si)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      x[si] = x[si] + alpha*d[si];
    }
    // update residual vector r=r-alpha*q, q=S*d
    if((iter+1) % 20 == 0){
      fn_operate_by_S(nSmat, x, r, process_y, stcO, comm);
      #pragma omp parallel for default(shared) private(si)
      #pragma loop noalias
      for(si=0;si<nSmat;++si) {
        r[si] = g[si] - r[si];
      }
    }else{
      #pragma omp parallel for default(shared) private(si)
      #pragma loop noalias
      for(si=0;si<nSmat;++si) {
        r[si] = r[si] - alpha*q[si];
      }
    }
    beta = CREAL(fn_xdot(nSmat,r,r))/delta;

    //update the norm of residual vector r
    delta = beta*delta;
    // update direction vector d
    #pragma omp parallel for default(shared) private(si)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      d[si] = r[si] + beta*d[si];
    }
  }

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): iter = %d\n", __FILE__, __LINE__, iter);
  fprintf(stderr, "DEBUG in %s (%d): End stcOptCG_Main\n", __FILE__, __LINE__);
#endif

  return iter;
}

/* calculate  y = S*x */
/* S is the overlap matrix*/
/* S[i][j] = OO[i+1][j+1] - OO[i+1][0] * OO[0][j+1]; */
int fn_operate_by_S(int nSmat, elemtype *x, elemtype *y, elemtype * process_y, elemtype * stcO, MPI_Comm comm) {
  elemtype *stcOd, *stcOs, *sdiag;
  int rank,size,info=0;
  int i,si;
  int incx = 1, incy = 1;
  double coef;
  elemtype a;
  elemtype alpha = 1.0, beta = 0.0;
  elemtype invW = 1.0/Wc;
#ifdef MVMC_SRCG_REAL
  char trans1='T';
#else
  char trans1='C';
#endif
  char trans2='N';

#ifdef _DEBUG_STCOPT_CG_GEMV
  elemtype malpha = -1.0;
  elemtype *y2;
  elemtype *S, *S2;
#endif
  stcOs = stcO  + nSmat;            //GetWorkSpaceDouble(nSmat*NVMCSample)
  stcOd = stcOs + nSmat*NVMCSample; //GetWorkSpaceDouble(NVMCSample)
  sdiag = stcOd + NVMCSample;            //GetWorkSpaceDouble(nSmat)
  
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Bcast(x, nSmat, MPI_DOUBLE, 0, comm);
  
  StartTimer(53);
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    y[si] = 0.0;
    process_y[si] = 0.0;
  }
  #pragma omp parallel for default(shared) private(i)
  #pragma loop noalias
  for(i=0;i<NVMCSample;++i) {
    stcOd[i] = 0.0;
  }

  /* Od[sample] = sum{sidx} x[sidx] * O[sidx][sample] ; */
  Gemv(&trans1, &nSmat, &NVMCSample, &alpha, stcOs, &nSmat, x, &incx, &beta, stcOd, &incy);
  /* process_y[sidx] = sum{sample} O[sidx][sample] * Od[sample] ; */
  Gemv(&trans2, &nSmat, &NVMCSample, &alpha, stcOs, &nSmat, stcOd, &incx, &beta, process_y, &incy);

#ifndef MVMC_SRCG_REAL
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    process_y[si] = creal(process_y[si]);
  }
#endif
  
  MPI_Barrier(comm);
  /* compute <OO>*x */
  allreduce(process_y, y, nSmat, comm);
 
  /* compute <O>*x */
  coef = CREAL(fn_xdot(nSmat, stcO, x));

  /* y = S*x */
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    y[si] = invW*y[si] - coef*CREAL(stcO[si]);
    /* modify diagonal elements */
    // y[si] += DSROptStaDel * x[si];

    y[si] += sdiag[si]*DSROptStaDel*x[si];
  }

  StopTimer(53);

#ifdef _DEBUG_STCOPT_CG_GEMV
  fprintf(stderr, "DEBUG in %s (%d): Start STCOPT_CG_GEMV\n", __FILE__, __LINE__);
  y2 = (elemtype*)calloc(nSmat, sizeof(elemtype));
  S = (elemtype*)calloc(nSmat*nSmat, sizeof(elemtype));
  S2 = (elemtype*)calloc(nSmat*nSmat, sizeof(elemtype));

  fprintf(stderr, "DEBUG in %s (%d): S = OO * OO'\n", __FILE__, __LINE__);
  Gemm(&trans2, &trans1, &nSmat, &nSmat, &NVMCSample, &alpha, stcOs, &nSmat, stcOs, &nSmat, &beta, S, &nSmat);

  fprintf(stderr, "DEBUG in %s (%d): S2 = O * O'\n", __FILE__, __LINE__);
  Ger(&nSmat, &nSmat, &alpha, stcO, &incx, stcO, &incy, S2, &nSmat);

  fprintf(stderr, "DEBUG in %s (%d): process_y = S*x\n", __FILE__, __LINE__);
  Gemv(&trans2, &nSmat, &nSmat, &alpha, S, &nSmat, x, &incx, &beta, process_y, &incy);
  fprintf(stderr, "DEBUG in %s (%d): allreduce process_y into y2\n", __FILE__, __LINE__);
  allreduce(process_y, y2, nSmat, comm);

  fprintf(stderr, "DEBUG in %s (%d): y2 = -S2*x + invW * y2\n", __FILE__, __LINE__);
  Gemv(&trans2, &nSmat, &nSmat, &malpha, S2, &nSmat, x, &incx, &invW, y2, &incy);
 
  /* modify diagonal elements */
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    y2[si] += DSROptStaDel * x[si];
  }

  fprintf(stderr, "DEBUG in %s (%d): index y1 y2 diff\n", __FILE__, __LINE__);
  for(si=0; si<nSmat; ++si){
    fprintf(stderr, "%d (%lg %lg) (%lg %lg) (%lg %lg)\n", si, CREAL(y[si]), CIMAG(y[si]), CREAL(y2[si]), CIMAG(y2[si]), CREAL(y[si]-y2[si]), CIMAG(y[si]-y2[si]));
  }
  
  free(S2);
  free(S);
  free(y2);
  fprintf(stderr, "DEBUG in %s (%d): End STCOPT_CG_GEMV\n", __FILE__, __LINE__);
#endif // _DEBUG_STCOPT_CG_GEMV

  return info;

}

void fn_StochasticOptCG_Init(const int nSmat, int *const smatToParaIdx, elemtype *VecCG) {
  const double dt = 2.0*DSROptStepDt;
  int i;
  int si,sj,pi,pj,idx,offset;
  elemtype *stcOd, *stcOs, *stcO;
  elemtype  *r, *g, *y, *d, *x, *process_y;
  elemtype *sdiag;
#ifdef MVMC_SRCG_REAL
  const double *srOptO=SROptOO_real;
  const double *srOptOOdiag=SROptOO_real + SROptSize;
  const double *srOptO_Store = SROptO_Store_real;
  const double *srOptHO=SROptHO_real;
#else
  const double complex *srOptO=SROptOO;
  const double complex *srOptOOdiag=SROptOO + 2*SROptSize;
  const double complex *srOptO_Store = SROptO_Store;
  const double complex *srOptHO=SROptHO;
#endif
  const double srOptHO_0 = CREAL(srOptHO[0]);

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
  sdiag = stcOd + NVMCSample;            //GetWorkSpaceDouble(nSmat)
  
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat*(NVMCSample+8) + NVMCSample;si++) {
    VecCG[si] = 0.0;
  }

  for(i=0;i<NVMCSample;++i) {
    offset = i*OFFSET*SROptSize;
    #pragma omp parallel for default(shared) private(si,pi,idx)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      pi = smatToParaIdx[si];
      
      idx = si + i*nSmat; /* column major */
      stcOs[idx] = CONJ(srOptO_Store[offset+pi+OFFSET]);
    }
  }

  /* calculate the energy gradient * (-dt) */
  /* energy gradient = 2.0*( HO[i+1] - H * O[i+1]) */
  #pragma omp parallel for default(shared) private(si,pi)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    pi = smatToParaIdx[si];
    
    stcO[si] = CREAL(srOptO[pi+OFFSET]);
    g[si] = -dt*(CREAL(srOptHO[pi+OFFSET]) - srOptHO_0 * CREAL(srOptO[pi+OFFSET]));

    sdiag[si] = CREAL(srOptOOdiag[pi+OFFSET]) - CREAL(srOptO[pi+OFFSET]) * CREAL(srOptO[pi+OFFSET]);
  }

#ifdef _DEBUG_STCOPT_CG
  for(i=0; i<nSmat; ++i){
    println_val(g[i]);
  }

  fprintf(stderr, "DEBUG in %s (%d): End stcOptCG_Init\n", __FILE__, __LINE__);
#endif

  return;
}

inline elemtype fn_xdot(const int n, elemtype * const p, elemtype * const q) {
  int i;
  elemtype z=0;
  #pragma loop noalias
  for(i=0;i<n;i++) {
    z += CONJ(p[i])*q[i];
  }
  return z;
}

#undef fn_StochasticOptCG
#undef fn_StochasticOptCG_Init
#undef fn_StochasticOptCG_Main
#undef fn_operate_by_S
#undef fn_xdot

#undef elemtype
#undef CONJ
#undef CREAL
#undef CIMAG
#undef Gemv
#undef Gemm
#undef Ger
#undef allreduce
#undef request_workspace
#undef get_workspace
#undef release_workspace
#undef print_val
#undef println_val

#undef OFFSET
