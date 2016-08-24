/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Variational parameters
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#define D_AmpMax             4.0

void InitParameter();
int ReadInitParameter(char *initFile);
void SyncModifiedParameter(MPI_Comm comm);
void SetFlagShift();

void shiftGJ();
double shiftDH2();
double shiftDH4();

/* initialize variational parameters */
void InitParameter() {
  int i;

  #pragma omp parallel for default(shared) private(i)
  for(i=0;i<NProj;i++) Proj[i] = 0.0;

  for(i=0;i<NSlater;i++){
    if(OptFlag[i+NProj] > 0){
      Slater[i] = genrand_real2(); /* uniform distribution [0,1) */
      printf("DEBUG: i=%d slater=%lf \n",i,Slater[i]);
    } else {
      Slater[i] = 0.0;
    }
  }

  for(i=0;i<NOptTrans;i++){
    OptTrans[i] = ParaQPOptTrans[i];
  }

  return;
}

/* read initial vaules of variational parameters from initFile */
int ReadInitParameter(char *initFile) {
  FILE *fp;
  int i,xi;
  double xtmp;

  fp = fopen(initFile, "r");
  if(fp!=NULL){
    while(fscanf(fp, "%lf ", &xtmp)!=EOF){
      for(i=1;i<4;i++) fscanf(fp, "%lf ", &xtmp);
      for(xi=0;xi<NProj;xi++) {
        fscanf(fp, "%lf %lf ", &(Proj[xi]), &xtmp);
      }
      for(xi=0;xi<NSlater;xi++) {
        fscanf(fp, "%lf %lf ", &(Slater[xi]), &xtmp);
      }
      for(xi=0;xi<NOptTrans;xi++) {
        fscanf(fp, "%lf %lf ", &(OptTrans[xi]), &xtmp);
      }
    }
    fclose(fp);
  } else { fprintf(stderr, "Error: %s does not exist.\n",initFile); }
  
  return 0;
}

/* sync and modify variational parameters */
void SyncModifiedParameter(MPI_Comm comm) {
  double gShift=0.0;
  double xmax,ratio;
  int i;

#ifdef _mpi_use
  int size;
  MPI_Comm_size(comm, &size);
  if(size>1) MPI_Bcast(Para, NPara, MPI_DOUBLE, 0, comm);
#endif /* _mpi_use */

  /***** shift correlation factors *****/
  /* shift the DH correlation factors */
  if(FlagShiftDH2==1) gShift += shiftDH2();
  if(FlagShiftDH4==1) gShift += shiftDH4();
  /* shift the Gutzwiller factors */
  for(i=0;i<NGutzwillerIdx;i++) Proj[i] += gShift;
  /* shift the Gutzwiller-Jastrow factors */
  if(FlagShiftGJ==1) shiftGJ();

  /***** rescale Slater *****/
  xmax = fabs(Slater[0]);
  for(i=1;i<NSlater;i++){
    if(xmax < fabs(Slater[i])) xmax = fabs(Slater[i]);
  }
  ratio = D_AmpMax/xmax;
  #pragma omp parallel for default(shared) private(i)
  for(i=0;i<NSlater;i++) Slater[i] *= ratio;

  /***** normalize OptTrans *****/
  if(FlagOptTrans>0){
    xmax = fabs(OptTrans[0]);
    for(i=1;i<NOptTrans;++i) {
      if(xmax < fabs(OptTrans[i])) xmax = fabs(OptTrans[i]);
    }
    ratio = 1.0/copysign(xmax, OptTrans[0]);
    #pragma omp parallel for default(shared) private(i)
    for(i=0;i<NOptTrans;++i) {
      OptTrans[i] *= ratio;
    }
  }

  return;
}

/* shift Gutzwiller-Jastrow factor */
void shiftGJ() {
  double shift=0.0;
  const int n = NGutzwillerIdx+NJastrowIdx;
  int i;

  if(NGutzwillerIdx==0||NJastrowIdx==0) return;

  for(i=0;i<n;i++) {
    shift += Proj[i];
  }
  shift /= (double)n;

  for(i=0;i<n;i++) {
    Proj[i] -= shift;
  }

  return;
}

/* shift 2-site DH Correlation factor */
/* The return value is gutzwiller shift width */
double shiftDH2() {
  double gShift=0.0;
  double shift;
  int xn,n0,n1,n2,offset;

  if(NDoublonHolon2siteIdx==0) return 0.0;

  /* 2-site doublon-holon correlation factor */
  offset = NGutzwillerIdx + NJastrowIdx;
  for(xn=0;xn<2*NDoublonHolon2siteIdx;xn++) { /* factor 2: d or h */
    /* n = offset + xn + (xdh+2*xm)*xNDoublonHolon2siteIdx; */
    n0 = offset + xn;
    n1 = n0 + 2*NDoublonHolon2siteIdx;
    n2 = n0 + 4*NDoublonHolon2siteIdx;
    shift = (Proj[n0]+Proj[n1]+Proj[n2])/3.0;
    Proj[n0] -= shift;
    Proj[n1] -= shift;
    Proj[n2] -= shift;
    gShift += shift;
  }

  return gShift;
}

/* shift 4-site DH Correlation factor */
/* The return value is gutzwiller shift width */
double shiftDH4() {
  double gShift=0.0;
  double shift;
  int xn,n0,n1,n2,n3,n4,offset;

  if(NDoublonHolon4siteIdx==0) return 0.0;

  offset = NGutzwillerIdx + NJastrowIdx + 6*NDoublonHolon2siteIdx;
  for(xn=0;xn<2*NDoublonHolon4siteIdx;xn++) { /* factor 2: d or h */
    /* n = offset + xn + (xdh+2*xm)*xNDoublonHolon4siteIdx; */
    n0 = offset + xn;
    n1 = n0 + 2*NDoublonHolon4siteIdx;
    n2 = n0 + 4*NDoublonHolon4siteIdx;
    n3 = n0 + 6*NDoublonHolon4siteIdx;
    n4 = n0 + 8*NDoublonHolon4siteIdx;
    shift = (Proj[n0]+Proj[n1]+Proj[n2]+Proj[n3]+Proj[n4])/5.0;
    Proj[n0] -= shift;
    Proj[n1] -= shift;
    Proj[n2] -= shift;
    Proj[n3] -= shift;
    Proj[n4] -= shift;
    gShift += shift;
  }

  return gShift;
}

void SetFlagShift() {
  int i,start,end;

  /* Gutzwiller */
  FlagShiftGJ=0;
  FlagShiftDH2=0;
  FlagShiftDH4=0;
  if(NGutzwillerIdx==0) return;
  start = 0;
  end = start + NGutzwillerIdx;
  for(i=start;i<end;i++) {
    if(OptFlag[i]!=1) return;
  }
  
  /* Jastrow */
  if(NJastrowIdx>0) {
    start = end;
    end   = start + NJastrowIdx;
    FlagShiftGJ=1;
    for(i=start;i<end;i++) {
      if(OptFlag[i]!=1) {
        FlagShiftGJ=0;
        break;
      }
    }
  }
  
  /* 2-site Doublon-Holon */
  if(NDoublonHolon2siteIdx>0) {
    start = end;
    end   = start + 6*NDoublonHolon2siteIdx;
    FlagShiftDH2=1;
    for(i=start;i<end;i++) {
      if(OptFlag[i]!=1) {
        FlagShiftDH2=0;
        break;
      }
    }
  }
  
  /* 4-site Doublon-Holon */
  if(NDoublonHolon4siteIdx>0) {
    start = end;
    end   = start + 10*NDoublonHolon4siteIdx;
    FlagShiftDH4=1;
    for(i=start;i<end;i++) {
      if(OptFlag[i]!=1) {
        FlagShiftDH4=0;
        break;
      }
    }
  }
  
  return;
}
