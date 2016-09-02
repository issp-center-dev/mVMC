#include "output.h"
#include "mfmemory.c"

int MakeOrbitalFile(struct BindStruct *X);
void cal_cisajs(struct BindStruct *X);

void WriteHeader(char* cNKWidx, int NKWidx, FILE *fp){
  fprintf(fp, "======================\n");
  fprintf(fp, cNKWidx);
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
    fprintf(fp_out, "%d % .18e % .18e \n",
            i, creal(Para[i]), cimag(Para[i]));
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

  sprintf(sdt,"%s_eigen.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  for(i=0;i< X->Def.Nsite*2;i++){
    tmp =  X->Large.EigenValues[i];
    fprintf(fp," %d  %.10lf \n",i+1,tmp);
  }
  fclose(fp);

  sprintf(sdt,"%s_gap.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  tmp =  X->Large.EigenValues[X->Def.Ne*2-1];
  tmp =  X->Large.EigenValues[X->Def.Ne*2] - tmp;
  fprintf(fp,"  %.10lf \n",tmp);
  fclose(fp);
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
            fprintf(fp," %4d %4d %4d %4d %.10lf %.10lf\n",site_1,site_2,spin_1,spin_2,cabs(tmp),carg(tmp));
          }
        } 
      }
    }
    fclose(fp);

}

int MakeOrbitalFile(struct BindStruct *X){
  int i, j, ispin, jspin, n;
  int isite, jsite;
  double complex **UHF_Fij;
  double complex *ParamOrbital;
  int *CountOrbital;
  int Orbitalidx;
  char fileName[256];
  
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
              UHF_Fij[isite][jsite]   +=  X->Large.R_SLT[isite][n]*X->Large.R_SLT[jsite][n+1];
            }
          }
        }
      }
    }

    c_malloc1(ParamOrbital, X->Def.NOrbitalIdx);
    i_malloc1(CountOrbital, X->Def.NOrbitalIdx);
    for(i=0; i<X->Def.NOrbitalIdx; i++){
      ParamOrbital[i]=0;
      CountOrbital[i]=0;
    }
    
    for(ispin=0; ispin<2; ispin++){
      for(jspin=0; jspin<2; jspin++){
        for(i=0;i< X->Def.Nsite;i++){
          for(j=0;j< X->Def.Nsite;j++){
            isite = i+ispin*X->Def.Nsite;
            jsite = j+jspin*X->Def.Nsite;
            Orbitalidx=X->Def.OrbitalIdx[isite][jsite];
            if(Orbitalidx !=-1){
              ParamOrbital[Orbitalidx]+=UHF_Fij[isite][jsite];
              CountOrbital[Orbitalidx]+=1;
            }
          }
        }
      }
    }

    for(i=0; i<X->Def.NOrbitalIdx; i++){
      ParamOrbital[i] /= (double)CountOrbital[i];
      //printf("debug: Orbital: idx=%d, param=%lf, %lf \n", i, creal(ParamOrbital[i]), cimag(ParamOrbital[i]));
    };
    
    sprintf(fileName, "%s_orbital_opt.dat", X->Def.CParaFileHead);
    Child_OutputOptData(fileName, "NOrbitalIdx", ParamOrbital, X->Def.NOrbitalIdx);
    
    c_free2(UHF_Fij, X->Def.Nsite*2, X->Def.Nsite*2);
    c_free1(ParamOrbital, X->Def.NOrbitalIdx);



  }
  else{
    return 0;
  }

  return 0;
}
