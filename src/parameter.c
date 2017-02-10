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
    if(OptFlag[2*i+2*NProj] > 0){ //TBC
      Slater[i] =  1*genrand_real2(); /* uniform distribution [0,1) */
      if(AllComplexFlag != 0){
        Slater[i] += 1*I*genrand_real2(); /* uniform distribution [0,1) */
      }
      //printf("DEBUG: i=%d slater=%lf %lf \n",i,creal(Slater[i]),cimag(Slater[i]));
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
  double tmp_real,tmp_comp;


  fp = fopen(initFile, "r");
  if(fp!=NULL){
    while(fscanf(fp, "%lf ", &xtmp)!=EOF){
      //for(i=1;i<4;i++) fscanf(fp, "%lf ", &xtmp);
      for(i=1;i<6;i++) fscanf(fp, "%lf ", &xtmp);
      for(xi=0;xi<NProj;xi++) {
        fscanf(fp, "%lf %lf %lf ", &tmp_real,&tmp_comp, &xtmp);
        Proj[xi] = tmp_real+tmp_comp*I; 
      }
      for(xi=0;xi<NSlater;xi++) {
        fscanf(fp, "%lf %lf %lf ", &tmp_real,&tmp_comp, &xtmp);
        Slater[xi] = tmp_real+tmp_comp*I; 
      }
      for(xi=0;xi<NOptTrans;xi++) {
        fscanf(fp, "%lf %lf %lf ", &tmp_real,&tmp_comp, &xtmp);
        OptTrans[xi] = tmp_real+tmp_comp*I; 
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
  if(size>1) MPI_Bcast(Para, NPara, MPI_DOUBLE_COMPLEX, 0, comm);
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
  xmax = cabs(Slater[0]);
  for(i=1;i<NSlater;i++){
    if(xmax < cabs(Slater[i])) xmax = cabs(Slater[i]);
  }
  ratio = D_AmpMax/xmax;
  #pragma omp parallel for default(shared) private(i)
  for(i=0;i<NSlater;i++) Slater[i] *= ratio;

  /***** normalize OptTrans *****/
  if(FlagOptTrans>0){
    xmax = cabs(OptTrans[0]);
    for(i=1;i<NOptTrans;++i) {
      if(xmax < cabs(OptTrans[i])) xmax = cabs(OptTrans[i]);
    }
    //ratio = 1.0/copysign(xmax, OptTrans[0]);
    ratio = 1.0/xmax; //TBC
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
    shift += creal(Proj[i]); // TBC
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
    shift = (creal(Proj[n0])+creal(Proj[n1])+creal(Proj[n2]))/3.0; //TBC
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
    shift = (creal(Proj[n0])+creal(Proj[n1])+creal(Proj[n2])+creal(Proj[n3])+creal(Proj[n4]))/5.0; //TBC
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
  if(NGutzwillerIdx==0) return; // no Gutz -> do nothing
  start = 0;
  end = start + NGutzwillerIdx;
  for(i=start;i<end;i++) {
    if(OptFlag[2*i]!=1) return;   // fixed Gutx -> do nothing
  }
  
  /* Jastrow */
  if(NJastrowIdx>0) {
    start = end;
    end   = start + NJastrowIdx;
    FlagShiftGJ=1; // unfixed J -> FlagShift on
    for(i=start;i<end;i++) {
      if(OptFlag[2*i]!=1) {       // fixed jast -> do nothing
        FlagShiftGJ=0;
        break;
      }
    }
  }
  
  /* 2-site Doublon-Holon */
  if(NDoublonHolon2siteIdx>0) {
    start = end;
    end   = start + 6*NDoublonHolon2siteIdx;
    FlagShiftDH2=1;         // unfixed D-H -> FlagShift on
    for(i=start;i<end;i++) {
      if(OptFlag[2*i]!=1) { // fixed D-H -> do nothing
        FlagShiftDH2=0;
        break;
      }
    }
  }
  
  /* 4-site Doublon-Holon */
  if(NDoublonHolon4siteIdx>0) {
    start = end;
    end   = start + 10*NDoublonHolon4siteIdx;
    FlagShiftDH4=1;         // unfixed D-H -> FlagShift on
    for(i=start;i<end;i++) { 
      if(OptFlag[2*i]!=1) { // fixed D-H -> do nothing
        FlagShiftDH4=0;
        break;
      }
    }
  }
  
  return;
}
