/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Stochastic Reconfiguration method by PDPOSV
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#ifdef _SYSTEM_A
 #define M_PDSYEVD PDSYEVD
 #define M_PDGEMV  PDGEMV
#elif _lapack_small_nounderscore
 #define M_PDSYEVD pdsyevd
 #define M_PDGEMV  pdgemv
#else
 #define M_NUMROC   numroc_
 #define M_DESCINIT descinit_
 #define M_PDSYEVD  pdsyevd_
 #define M_PDPOSV  pdposv_
 #define M_PDGEMV  pdgemv_
#endif

int Csys2blacs_handle(MPI_Comm comm);
int Cblacs_gridinit(int *ictxt, char *order, int nprow, int npcol);
int Cblacs_gridexit(int ictxt);
int Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *myprow, int *mypcol);
int M_NUMROC(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
int M_DESCINIT(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc,
               int *icsrc, int *ictxt, int *lld, int *info);
int M_PDPOSV(char *uplo, int *n, int *nrhs, double *a, int *ia, int *ja, int *desca,
             double *b, int *ib, int *jb, int *descb, int *info);
void M_PDGEMV(char *trans, int *m, int *n, double *alpha,
              double *a, int *ia, int *ja, int *desca,
              double *x, int *ix, int *jx, int *descx, int *incx,
              double *beta,
              double *y, int *iy, int *jy, int *descy, int *incy);

int StochasticOpt(MPI_Comm comm);
int stcOptMain(double *r, const int nSmat, const int *smatToParaIdx, MPI_Comm comm);

int StochasticOpt(MPI_Comm comm) {
  const int nPara=NPara;
  const int srOptSize=SROptSize;
  const double *srOptOO=SROptOO;

  double r[SROptSize]; /* the parameter change */
  int nSmat;
  int smatToParaIdx[NPara];

  int cutNum=0,optNum=0;
  double sDiag,sDiagMax,sDiagMin;
  double diagCutThreshold;

  int si; /* index for matrix S */
  int pi; /* index for variational parameters */

  double rmax;
  int simax;
  int info=0;

  double *para=Para;

  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  StartTimer(50);

  #pragma omp parallel for default(shared) private(pi)
  #pragma loop noalias
  for(pi=0;pi<nPara;pi++) {
    /* r[i] is temporarily used for diagonal elements of S */
    /* S[i][i] = OO[pi+1][pi+1] - OO[0][pi+1] * OO[0][pi+1]; */
    r[pi] = srOptOO[(pi+1)*srOptSize+pi+1] - srOptOO[pi+1] * srOptOO[pi+1];
  }

  sDiag = r[0];
  sDiagMax=sDiag; sDiagMin=sDiag;
  for(pi=0;pi<nPara;pi++) {
    sDiag = r[pi];
    if(sDiag>sDiagMax) sDiagMax=sDiag;
    if(sDiag<sDiagMin) sDiagMin=sDiag;
  }

  diagCutThreshold = sDiagMax*DSROptRedCut;
  si = 0;
  for(pi=0;pi<nPara;pi++) {
    if(OptFlag[pi]!=1) { /* fixed by OptFlag */
      optNum++;
      continue;
    }

    sDiag = r[pi];
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

  StopTimer(50);
  StartTimer(51);

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
  return info;
}

/* calculate the parameter change r[nSmat] from SOpt.
   The result is gathered in rank 0. */
int stcOptMain(double *r, const int nSmat, const int *smatToParaIdx, MPI_Comm comm) {
  /* global vector */
  double *w; /* workspace */
  /* distributed matrix (row x col) */
  double *s; /* overlap matrix (nSmat x nSmat) */
  double *g; /* energy gradient and parameter change (nSmat x 1) */

  int sSize,gSize;
  int wSize=nSmat;

  /* for MPI */
  int rank,size;
  int dims[2]={0,0};
  MPI_Comm comm_col;

  /* index table */
  int *irToSmatIdx; /* table for local indx to Smat index */
  int *irToParaIdx, *icToParaIdx; /* tables for local index to Para index */

  /* for BLACS */
  int ictxt,nprow,npcol,myprow,mypcol;
  char procOrder='R';

  /* for array descriptor of SCALAPACK */
  int m,n;
  int mlld,mlocr,mlocc; /* for matrix */
  int vlld,vlocr,vlocc; /* for vector */
  int mb=64, nb=64; /* blocking factor */
  int irsrc=0, icsrc=0;
  int descs[9],descg[9];

  /* for PDSYEVD */
  char uplo;
  int nrhs, is, js, ig, jg;

  int info;
  const double ratioDiag = 1.0+DSROptStaDel;
  int si,pi,pj,idx;
  int ir,ic;

  const int srOptSize=SROptSize;
  const double dSROptStepDt = DSROptStepDt;
  const double srOptHO_0=SROptHO[0];
  double *srOptOO=SROptOO;
  double *srOptHO=SROptHO;

  StartTimer(55);

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  MPI_Dims_create(size,2,dims);

  /* initialize the BLACS context */
  nprow=dims[0]; npcol=dims[1];
  ictxt = Csys2blacs_handle(comm);
  Cblacs_gridinit(&ictxt, &procOrder, nprow, npcol);
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myprow, &mypcol);

  /* initialize array descriptors for distributed matrix s */
  m=n=nSmat; /* matrix size */
  mlocr = M_NUMROC(&m, &mb, &myprow, &irsrc, &nprow);
  mlocc = M_NUMROC(&n, &nb, &mypcol, &icsrc, &npcol);
  mlld = (mlocr>0) ? mlocr : 1;
  M_DESCINIT(descs, &m, &n, &mb, &nb, &irsrc, &icsrc, &ictxt, &mlld, &info);
  sSize = (mlocr*mlocc>0) ? mlocr*mlocc : 1;

  /* initialize array descriptors for distributed vector g */
  m=nSmat; n=1; /* vector */
  vlocr = mlocr; /* M_NUMROC(&m, &mb, &myprow, &irsrc, &nprow); */
  vlocc = M_NUMROC(&n, &nb, &mypcol, &icsrc, &npcol);
  vlld = (vlocr>0) ? vlocr : 1;
  M_DESCINIT(descg, &m, &n, &mb, &nb, &irsrc, &icsrc, &ictxt, &vlld, &info);
  gSize = (vlocr*vlocc>0) ? vlocr*vlocc : 1;

  /* allocate memory */
  RequestWorkSpaceDouble(sSize+gSize+wSize);
  s = GetWorkSpaceDouble(sSize);
  g = GetWorkSpaceDouble(gSize);
  w = GetWorkSpaceDouble(wSize);

  /* Para indices of the distributed vector and calculate them */
  RequestWorkSpaceInt(2*mlocr+mlocc);
  irToSmatIdx = GetWorkSpaceInt(mlocr);
  irToParaIdx = GetWorkSpaceInt(mlocr);
  icToParaIdx = GetWorkSpaceInt(mlocc);

  #pragma omp parallel for default(shared) private(ir,si)
  for(ir=0;ir<mlocr;ir++) {
    si = (ir/mb)*nprow*mb + myprow*mb + (ir%mb);
    irToSmatIdx[ir] = si;
    irToParaIdx[ir] = smatToParaIdx[si];
  }
  #pragma omp parallel for default(shared) private(ic)
  for(ic=0;ic<mlocc;ic++) {
    icToParaIdx[ic] = smatToParaIdx[(ic/nb)*npcol*nb + mypcol*nb + (ic%nb)];
  }

  StopTimer(55);
  StartTimer(56);
  /* calculate the overlap matrix S */
  #pragma omp parallel for default(shared) private(ic,ir,pi,pj,idx)
  #pragma loop noalias
  for(ic=0;ic<mlocc;ic++) {
    pj = icToParaIdx[ic]; /* Para index (global) */
    for(ir=0;ir<mlocr;ir++) {
      pi = irToParaIdx[ir]; /* Para index (global) */
      idx = ir + ic*mlocr; /* local index (row major) */

      /* S[i][j] = xOO[i+1][j+1] - xOO[0][i+1] * xOO[0][j+1]; */
      s[idx] = srOptOO[(pi+1)*srOptSize+(pj+1)] - srOptOO[pi+1] * srOptOO[pj+1];
      /* modify diagonal elements */
      if(pi==pj) s[idx] *= ratioDiag;
    }
  }

  /* calculate the energy gradient g and multiply (-dt) */
  if(vlocc>0) {
    #pragma omp parallel for default(shared) private(ir,pi)
    #pragma loop noalias
    #pragma loop norecurrence irToParaIdx
    for(ir=0;ir<vlocr;ir++) {
      pi = irToParaIdx[ir]; /* Para index (global) */
      
      /* energy gradient = 2.0*( xHO[i+1] - xHO[0] * xOO[0][i+1]) */
      /* g[i] = -dt * (energy gradient) */
      g[ir] = -dSROptStepDt*2.0*(srOptHO[pi+1] - srOptHO_0 * srOptOO[pi+1]);
    }
  }

  StopTimer(56);
  StartTimer(57);
  
  /***** solve the linear equation S*r=g by PDPOSV *****/
  uplo='U'; n=nSmat; nrhs=1; is=1; js=1; ig=1; jg=1;
  M_PDPOSV(&uplo, &n, &nrhs, s, &is, &js, descs, 
           g, &ig, &jg, descg, &info); /* g is overwritten with the solution. */
  /* error handle */
  if(info!=0) {
    if(rank==0) fprintf(stderr,"error: PDPOSV info=%d\n",info);
    return info;
  }
  /* end of diagonalization */

  StopTimer(57);
  StartTimer(58);

  /* create a communicator whose process has the same value of mypcol */
  MPI_Comm_split(comm,mypcol,myprow,&comm_col); 

  /* clear workspace */
  for(si=0;si<nSmat;si++) w[si] = 0.0;

  if(mypcol==0) {
    /* copy the solution to workspace */
    #pragma omp parallel for default(shared) private(ir,si)
    #pragma loop noalias
    #pragma loop norecurrence irToSmatIdx
    for(ir=0;ir<vlocr;ir++) {
      si = irToSmatIdx[ir]; /* Smat index (global) */
      w[si] = g[ir];
    }

    /* gather the solution to r on rank=0 process */
    SafeMpiReduce(w,r,nSmat,comm_col);
  }

  StopTimer(58);

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceDouble();
  MPI_Comm_free(&comm_col);
  Cblacs_gridexit(ictxt);
  return info;
}
