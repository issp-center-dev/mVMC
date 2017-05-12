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
  #define fn_print_Smat_stderr print_Smat_stderr_real

  #define elemtype double
  #define CONJ(x) (x)
  #define CREAL(x) (x)
  #define CIMAG(x) (0.0)

  #define OFFSET (1)
  #define USE_IMAG (0)
  #define SIZE_VecCG (nSmat*8 + NVMCSample*(nSmat+1))
#else // MVMC_SRCG_REAL
  #define fn_StochasticOptCG StochasticOptCG_fcmp
  #define fn_StochasticOptCG_Init StochasticOptCG_Init_fcmp
  #define fn_StochasticOptCG_Main StochasticOptCG_Main_fcmp
  #define fn_operate_by_S operate_by_S_fcmp
  #define fn_print_Smat_stderr print_Smat_stderr_fcmp

  #define elemtype double complex
  #define CONJ(x) conj(x)
  #define CREAL(x) creal(x)
  #define CIMAG(x) cimag(x)

  #define OFFSET (2)
  #define USE_IMAG (1)
  #define SIZE_VecCG (nSmat*8 + 2*NVMCSample*(nSmat+1))
#endif

#define print_val(x) fprintf(stderr, "%lg ", CREAL(x))
#define println_val(x) fprintf(stderr, "%lg\n", CREAL(x))

/*
   x :: nSmat
   g :: nSmat
   sdiag :: nSmat
   stcO :: nSmat
   stcOs_real :: nSmat*NVMCSample
   stcOs_imag :: nSmat*NVMCSample (Complex) or 0 (Real)
   y_real :: NVMCSample
   y_imag :: NVMCSample (Complex) or 0 (Real)
   z_local :: nSmat
   q :: nSmat
   d :: nSmat
   r :: nSmat
*/

int fn_StochasticOptCG(MPI_Comm comm);
void fn_StochasticOptCG_Init(const int nSmat, int *const smatToParaIdx, double *VecCG); 
int fn_StochasticOptCG_Main(const int nSmat, double *VecCG, MPI_Comm comm);
int fn_operate_by_S(const int nSmat, double *x, double *z, double *VecCG, MPI_Comm comm);
void fn_print_Smat_stderr(const int nSmat, double *VecCG, MPI_Comm comm);

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

  double rmax;
  int simax;
  int info=0;

  double complex *para=Para;
  double *VecCG;
  double *r;


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
  
  RequestWorkSpaceDouble(SIZE_VecCG);
  VecCG = GetWorkSpaceDouble(SIZE_VecCG);
  fn_StochasticOptCG_Init(nSmat, smatToParaIdx, VecCG);

#ifdef _DEBUG_STCOPT_CG_PRINT_SMAT
  fn_print_Smat_stderr(nSmat, VecCG, comm);
  abort();
#endif

  info = fn_StochasticOptCG_Main(nSmat, VecCG, comm);
#ifdef _DEBUG_STCOPT_CG
  for(si=0; si<nSmat; ++si){
    println_val(VecCG[si]);
  }
#endif

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
  ReleaseWorkSpaceDouble();

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): End StochasticOptCG\n", __FILE__, __LINE__);
#endif
  return info;
}

/* calculate the parameter change r[nSmat] from SOpt.
   Solve S*x = g */
int fn_StochasticOptCG_Main(const int nSmat, double *VecCG, MPI_Comm comm) {
  int si,pi,pj,idx;
  int rank, size, info;
  int iter;
  int max_iter=NSROptCGMaxIter;
  double delta, beta;
  double alpha;
  double cg_thresh = DSROptCGTol*DSROptCGTol * (double)nSmat * (double)nSmat;
  //double cg_thresh = DSROptRedCut;

  const int srOptSize=SROptSize;
  double *x, *g, *sdiag, *stcO, *stcOs_real;
  double *stcOs_imag, *y_real, *y_imag, *z_local;
  double *q, *d, *r;

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): Start stcOptCG_Main\n", __FILE__, __LINE__);
#endif
  
  x = VecCG;
  g = x + nSmat;
  sdiag = g + nSmat;
  stcO = sdiag + nSmat;
  stcOs_real = stcO + nSmat;
  stcOs_imag = stcOs_real + NVMCSample*nSmat;
  y_real = stcOs_imag + USE_IMAG*NVMCSample*nSmat;
  y_imag = y_real + NVMCSample;
  z_local = y_imag + USE_IMAG*NVMCSample;
  q = z_local + nSmat;
  d = q + nSmat;
  r = d + nSmat;
  
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    d[si] = r[si] = g[si];
  }

  delta = xdot(nSmat, r, r);

  for(iter=0; iter < max_iter; iter++){
    //check convergence 
    //if(rank==0) printf("%d: %d %.3e\n",rank, iter,delta);
#ifdef _DEBUG_STCOPT_CG
    fprintf(stderr, "delta = %lg, cg_thresh = %lg\n", delta, cg_thresh);
#endif
    if (delta < cg_thresh) break;

    // compute vector q=S*d
    fn_operate_by_S(nSmat, d, q, VecCG, comm);
    alpha = delta/xdot(nSmat,d,q);
  
    // update solution vector x=x+alpha*d
    #pragma omp parallel for default(shared) private(si)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      x[si] = x[si] + alpha*d[si];
    }
    // update residual vector r=r-alpha*q, q=S*d
    if((iter+1) % 20 == 0){
      fn_operate_by_S(nSmat, x, r, VecCG, comm);
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
    beta = xdot(nSmat,r,r)/delta;

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

/* calculate  z = S*x */
/* S is the overlap matrix*/
/* S[i][j] = OO[i+1][j+1] - OO[i+1][0] * OO[0][j+1]; */
int fn_operate_by_S(int nSmat, double *x, double *z, double *VecCG, MPI_Comm comm) {
  int rank,size,info=0;
  int i,si;
  int incx = 1, incy = 1;
  double coef;
  double one = 1.0, zero = 0.0;
  double invW = 1.0/Wc;
  char transT='T';
  char transN='N';

  double *sdiag, *stcO, *stcOs_real, *stcOs_imag;
  double *y_real, *y_imag;
  double *z_local;

  sdiag = VecCG + 2*nSmat;
  stcO = sdiag + nSmat;
  stcOs_real = stcO + nSmat;
  stcOs_imag = stcOs_real + NVMCSample*nSmat;
  y_real = stcOs_imag + USE_IMAG*NVMCSample*nSmat;
  y_imag = y_real + NVMCSample;
  z_local = y_imag + USE_IMAG*NVMCSample;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Bcast(x, nSmat, MPI_DOUBLE, 0, comm);
  
  StartTimer(53);
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    z[si] = 0.0;
    z_local[si] = 0.0;
  }

  // y_real[sample] = sum{si} x[si] * O_real[si][sample]
  M_DGEMV(&transT, &nSmat, &NVMCSample, &one, stcOs_real, &nSmat, x, &incx, &zero, y_real, &incy);
#ifndef MVMC_SRCG_REAL
  // y_imag[sample] = sum{si} x[si] * O_imag[si][sample]
  M_DGEMV(&transT, &nSmat, &NVMCSample, &one, stcOs_imag, &nSmat, x, &incx, &zero, y_imag, &incy);
#endif

  // z_local[si] = sum{sample} O_real[si][sample] * y_real[sample] + O_imag[si][sample] * y_imag[sample]
  M_DGEMV(&transN, &nSmat, &NVMCSample, &one, stcOs_real, &nSmat, y_real, &incx, &zero, z_local, &incy);

#ifndef MVMC_SRCG_REAL
  M_DGEMV(&transN, &nSmat, &NVMCSample, &one, stcOs_imag, &nSmat, y_imag, &incx, &one, z_local, &incy);
#endif

  MPI_Barrier(comm);
  /* compute <OO>*x */
  SafeMpiAllReduce(z_local, z, nSmat, comm);
 
  /* compute <O>*x */
  coef = xdot(nSmat, stcO, x);

  /* y = S*x */
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    z[si] = invW*z[si] - coef*stcO[si];
    /* modify diagonal elements */
    // z[si] += DSROptStaDel * x[si];

    z[si] += sdiag[si]*DSROptStaDel*x[si];
  }

  StopTimer(53);

  return info;
}

void fn_StochasticOptCG_Init(const int nSmat, int *const smatToParaIdx, double *VecCG) {
  const double dt = 2.0*DSROptStepDt;
  int i;
  int si,sj,pi,pj,idx,offset;

  double *x, *g;
  double *sdiag, *stcO, *stcOs_real, *stcOs_imag;
  double *y_real, *y_imag, *z_local;
  double *q, *d, *r;

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
  
  x = VecCG;
  g = x + nSmat;
  sdiag = g + nSmat;
  stcO = sdiag + nSmat;
  stcOs_real = stcO + nSmat;
  stcOs_imag = stcOs_real + NVMCSample*nSmat;
  
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<SIZE_VecCG;si++) {
    VecCG[si] = 0.0;
  }

  for(i=0;i<NVMCSample;++i) {
    offset = i*OFFSET*SROptSize;
    #pragma omp parallel for default(shared) private(si,pi,idx)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      pi = smatToParaIdx[si];
      
      idx = si + i*nSmat; /* column major */
      stcOs_real[idx] = CREAL(srOptO_Store[offset+pi+OFFSET]);
#ifndef MVMC_SRCG_REAL
      stcOs_imag[idx] = CIMAG(srOptO_Store[offset+pi+OFFSET]);
#endif
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

void fn_print_Smat_stderr(const int nSmat, double *VecCG, MPI_Comm comm){
  int si, sj;
  double* S = (double*)calloc(nSmat*nSmat, sizeof(double));
  double* xs = (double*)calloc(nSmat, sizeof(double));

  for(si=0; si<nSmat; ++si){
    xs[si] = 1.0;
    fn_operate_by_S(nSmat, xs, S+si*nSmat, VecCG, comm);
    xs[si] = 0.0;
  }

  for(si=0; si<nSmat; ++si){
  for(sj=0; sj<nSmat; ++sj){
    fprintf(stderr, "%lg ", S[si+nSmat*sj]);
  }
  fprintf(stderr, "\n");
  }

  free(xs);
  free(S);
}

#undef fn_StochasticOptCG
#undef fn_StochasticOptCG_Init
#undef fn_StochasticOptCG_Main
#undef fn_operate_by_S
#undef fn_print_Smat

#undef elemtype
#undef CONJ
#undef CREAL
#undef CIMAG
#undef print_val
#undef println_val
#undef OFFSET
#undef USE_IMAG
#undef SIZE_VecCG

