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
#include "output.h"
#include "mfmemory.c"
#include "SFMT.h"
#include <matrixlapack.h>

int MakeOrbitalFile(struct BindStruct *X);
void cal_cisajs(struct BindStruct *X);

void WriteHeader(char* cNKWidx, int NKWidx, FILE *fp){
  fprintf(fp, "======================\n");
  fprintf(fp, cNKWidx, 0);
  fprintf(fp, "  %d\n", NKWidx);
  fprintf(fp, "======================\n");
  fprintf(fp, "======================\n");
  fprintf(fp, "======================\n");    
}

void Child_OutputOptData(char* cFileName, char* cNKWidx, double complex *Para, int Nidx){
  FILE *fp_out;
  int i;
  fp_out = fopen(cFileName, "w");
  WriteHeader(cNKWidx, Nidx,  fp_out);
  for(i=0; i<Nidx; i++){
    fprintf(fp_out, "%d % .18e % .18e \n",i,creal(Para[i]),cimag(Para[i]));
  }
  fclose(fp_out);
}

void output(struct BindStruct *X){
  
  FILE *fp;
  char sdt[256];
  int i;
  double tmp;


  sprintf(sdt,"%s_result.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  fprintf(fp," energy %.10lf \n",X->Phys.energy);
  fprintf(fp," num    %.10lf \n",X->Phys.num);
  fclose(fp);
  printf("Energy and num are outputted to %s.\n", sdt);

  sprintf(sdt,"%s_eigen.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  for(i=0;i< X->Def.Nsite*2;i++){
    tmp =  X->Large.EigenValues[i];
    fprintf(fp," %d  %.10lf \n",i+1,tmp);
  }
  fclose(fp);
  printf("Eigenvalues are outputted to %s.\n", sdt);

  sprintf(sdt,"%s_gap.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  tmp =  X->Large.EigenValues[X->Def.Ne*2-1];
  tmp =  X->Large.EigenValues[X->Def.Ne*2] - tmp;
  fprintf(fp,"  %.10lf \n",tmp);
  fclose(fp);
  printf("Energy gaps are outputted to %s.\n", sdt);
  cal_cisajs(X);

  MakeOrbitalFile(X);
}

void cal_cisajs(struct BindStruct *X){
    FILE *fp;
    char sdt[256];
    int site_1,site_2,spin_1,spin_2;
    int Ns,t_site_1,t_site_2;
    double complex tmp;
    Ns = X->Def.Nsite; 
    sprintf(sdt, "%s_UHF_cisajs.dat", X->Def.CDataFileHead);
    fp = fopen(sdt,"w");
    for(site_1=0; site_1< Ns; site_1++){
      for(site_2=0; site_2< Ns; site_2++){
        for(spin_1=0;spin_1<2;spin_1++){
          for(spin_2=0;spin_2<2;spin_2++){
            t_site_1    = site_1+Ns*spin_1;
            t_site_2    = site_2+Ns*spin_2;
            tmp         = X->Large.G[t_site_1][t_site_2];
            
			  fprintf(fp, " %4d %4d %4d %4d %.10lf %.10lf\n", site_1, spin_1, site_2, spin_2, creal(tmp), cimag(tmp));
			  
			//  if(t_site_1==t_site_2) {
		  //		  fprintf(stdout, " Debug: %4d %4d %4d %4d %.10lf %.10lf\n", site_1, spin_1, site_2, spin_2, cabs(tmp), carg(tmp));
			//  }
			  
          }
        } 
      }
    }
    fclose(fp);
    printf("Onebody Green's functions are outputted to %s.\n", sdt);
}

int MakeOrbitalFile(struct BindStruct *X){
  int i, j, ispin, jspin, n;
  int isite, jsite;
  double complex **UHF_Fij;
  double complex *ParamOrbital;
  int *CountOrbital;
  int Orbitalidx;
  int int_i,int_j,int_k,int_l,xMsize;
  double complex **tmp_mat,**vec;
  double *r;
  char fileName[256];
  int ini,fin,tmp_i;

/*this part only for anti-parallel
//[s] for anti-pararell, rediag
  xMsize = X->Def.Nsite;  
  c_malloc2(tmp_mat,xMsize,xMsize);
  c_malloc2(vec,xMsize,xMsize);
  d_malloc1(r,xMsize);
  for(int_l = 0; int_l < 2*xMsize; int_l++){
    for(int_k = 0; int_k < 2*X->Def.Ne; int_k++){
      X->Large.R_SLT[int_l][int_k] = 0.0;
    }
  } 
// for up 
  for(int_i = 0;int_i < xMsize; int_i++){
    for(int_j = 0;int_j < xMsize; int_j++){
      tmp_mat[int_i][int_j] = X->Large.Ham[int_i][int_j];
    }
  }
  ZHEEVall(xMsize,tmp_mat,r,vec);
  for(int_k = 0; int_k < X->Def.Ne; int_k++){
    for(int_l = 0; int_l < xMsize; int_l++){
      X->Large.R_SLT[int_l][2*int_k] = (vec[int_k][int_l]);
    }
  }
// for down
  for(int_i = 0;int_i < xMsize; int_i++){
    for(int_j = 0;int_j < xMsize; int_j++){
      tmp_mat[int_i][int_j] = X->Large.Ham[int_i+xMsize][int_j+xMsize];
    }
  }
  ZHEEVall(xMsize,tmp_mat,r,vec);
  for(int_k = 0; int_k < X->Def.Ne; int_k++){
    for(int_l = 0; int_l < xMsize; int_l++){
      X->Large.R_SLT[int_l+xMsize][2*int_k+1] = (vec[int_k][int_l]);
    }
  }
//[e] for anti-pararell
*/
 
  
  if(X->Def.NOrbitalIdx>0){
    c_malloc2(UHF_Fij, X->Def.Nsite*2, X->Def.Nsite*2);
    for(ispin=0; ispin<2; ispin++){
      for(jspin=0; jspin<2; jspin++){
        for(i=0;i< X->Def.Nsite;i++){
          for(j=0;j< X->Def.Nsite;j++){
            isite = i+ispin*X->Def.Nsite;
            jsite = j+jspin*X->Def.Nsite;
            UHF_Fij[isite][jsite]=0;
            for(n=0;n< 2*X->Def.Ne;n+=2){
              UHF_Fij[isite][jsite]   +=  conj(X->Large.R_SLT[isite][n])*conj(X->Large.R_SLT[jsite][n+1])- conj(X->Large.R_SLT[isite][n+1])*conj(X->Large.R_SLT[jsite][n]);
            }
          }
        }
      }
    }
    //printf(" %d %d %d \n",X->Def.NOrbitalAP,X->Def.NOrbitalP,X->Def.NOrbitalIdx);
//[s] For AntiParallel
    c_malloc1(ParamOrbital, X->Def.NOrbitalIdx);
    i_malloc1(CountOrbital, X->Def.NOrbitalIdx);
    //c_malloc1(ParamOrbital, X->Def.NOrbitalAP);
    //i_malloc1(CountOrbital, X->Def.NOrbitalAP);
    for (i = 0; i < X->Def.NOrbitalAP; i++) {
      ParamOrbital[i] = 0;
      CountOrbital[i] = 0;
    }
    for(i = 0; i < X->Def.Nsite; i++) {
      for(j = 0; j < X->Def.Nsite; j++) {
        isite = i + 0 * X->Def.Nsite;
        jsite = j + 1 * X->Def.Nsite;
        Orbitalidx = X->Def.OrbitalIdx[isite][jsite];
        if(Orbitalidx != -1) {
          ParamOrbital[Orbitalidx] += UHF_Fij[isite][jsite];
          CountOrbital[Orbitalidx] += 1;
        }
      }
    }
    for (i = 0; i < X->Def.NOrbitalAP; i++) {
      ParamOrbital[i] /= (double) CountOrbital[i];
      ParamOrbital[i] += genrand_real2() * pow(10.0, -X->Def.eps_int_slater);
    }
    sprintf(fileName, "%s_APOrbital_opt.dat", X->Def.CParaFileHead);
    Child_OutputOptData(fileName, "NOrbitalAP", ParamOrbital, X->Def.NOrbitalAP);

    //c_free1(ParamOrbital, X->Def.NOrbitalAP);
    //i_free1(CountOrbital, X->Def.NOrbitalAP);
//[e] For AntiParallel

//[s] For Parallel
    //c_malloc1(ParamOrbital, X->Def.NOrbitalP);
    //i_malloc1(CountOrbital, X->Def.NOrbitalP);
    ini =  X->Def.NOrbitalAP;
    fin =  X->Def.NOrbitalAP + X->Def.NOrbitalP;
    for (i =  ini; i < fin; i++) {
      ParamOrbital[i] = 0;
      CountOrbital[i] = 0;
    }
    for(i = 0; i < X->Def.Nsite; i++) {
      for(j = i+1; j < X->Def.Nsite; j++) {
        isite = i + 0 * X->Def.Nsite;
        jsite = j + 0 * X->Def.Nsite;
        Orbitalidx = X->Def.OrbitalIdx[isite][jsite];
        if(Orbitalidx != -1) {
          ParamOrbital[Orbitalidx] += UHF_Fij[isite][jsite];
          CountOrbital[Orbitalidx] += 1;
        }
        isite = i + 1 * X->Def.Nsite;
        jsite = j + 1 * X->Def.Nsite;
        Orbitalidx = X->Def.OrbitalIdx[isite][jsite];
        if(Orbitalidx != -1) {
          ParamOrbital[Orbitalidx] += UHF_Fij[isite][jsite];
          CountOrbital[Orbitalidx] += 1;
        }
      }
    }
    for (i =  ini; i < fin; i++) {
      printf("MDEBUG: %d %d\n",i, CountOrbital[i]);
      ParamOrbital[i] /= (double) CountOrbital[i];
      ParamOrbital[i] += genrand_real2() * pow(10.0, -X->Def.eps_int_slater);
    }
    for (i =  ini; i < fin; i++) {
      tmp_i = i-ini;
      ParamOrbital[tmp_i] = ParamOrbital[i];
    }
    sprintf(fileName, "%s_POrbital_opt.dat", X->Def.CParaFileHead);
    Child_OutputOptData(fileName, "NOrbitalP", ParamOrbital, X->Def.NOrbitalP);

    //c_free1(ParamOrbital, X->Def.NOrbitalP);
    //i_free1(CountOrbital, X->Def.NOrbitalP);
//[e] For Parallel

//[s] For general
    //c_malloc1(ParamOrbital, X->Def.NOrbitalIdx);
    //i_malloc1(CountOrbital, X->Def.NOrbitalIdx);
    for (i = 0; i < X->Def.NOrbitalIdx; i++) {
      ParamOrbital[i] = 0;
      CountOrbital[i] = 0;
    }
    for (ispin = 0; ispin < 2; ispin++) {
      for (jspin = 0; jspin < 2; jspin++) {
        for (i = 0; i < X->Def.Nsite; i++) {
          for (j = 0; j < X->Def.Nsite; j++) {
            isite = i + ispin * X->Def.Nsite;
            jsite = j + jspin * X->Def.Nsite;
            Orbitalidx = X->Def.OrbitalIdx[isite][jsite];
            if (Orbitalidx != -1) {
              // ParamOrbital[Orbitalidx]+=UHF_Fij[isite][jsite];
              ParamOrbital[Orbitalidx] += UHF_Fij[isite][jsite];
              CountOrbital[Orbitalidx] += 1;
              //printf("debug: Orbitaidx[%d][%d]=%d, UHF_Fij=%lf, %lf \n", isite, jsite, Orbitalidx, creal(UHF_Fij[isite][jsite]), cimag(UHF_Fij[isite][jsite]));
              //printf("debug: Orbitaidx[%d][%d]=%d, ParamOrbital_Fij=%lf, %lf \n", isite, jsite, Orbitalidx, creal(ParamOrbital[Orbitalidx]), cimag(ParamOrbital[Orbitalidx]));
            }
          }
        }
      }
    }
    for (i = 0; i < X->Def.NOrbitalIdx; i++) {
      ParamOrbital[i] /= (double) CountOrbital[i];
      ParamOrbital[i] += genrand_real2() * pow(10.0, -X->Def.eps_int_slater);
      //printf("debug: Orbital: idx=%d, param=%lf, %lf , count=%d \n", i, creal(ParamOrbital[i]), cimag(ParamOrbital[i]), CountOrbital[i]);
    }
    sprintf(fileName, "%s_GeneralOrbital_opt.dat", X->Def.CParaFileHead);
    Child_OutputOptData(fileName, "NOrbitalIdx", ParamOrbital, X->Def.NOrbitalIdx);

    c_free1(ParamOrbital, X->Def.NOrbitalIdx);
    i_free1(CountOrbital, X->Def.NOrbitalIdx);
//[e] For general
    printf("Fij for mVMC are outputted to %s.\n", fileName);
    c_free2(UHF_Fij, X->Def.Nsite * 2, X->Def.Nsite * 2);
  }
  return 0;
}
