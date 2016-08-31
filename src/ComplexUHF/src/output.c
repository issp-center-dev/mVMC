#include "output.h"

void output_cisajs(struct BindStruct *X);

void output(struct BindStruct *X){
  
  FILE *fp;
  char sdt[256];
  int i,i_max;
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

  sprintf(sdt,"%s_fij.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  /*
  for(i=0;i< X->Def.Nsite;i++){
    for(j=0;j< X->Def.Nsite;j++){
    tmp    =  0.0;
      for(n=0;n< X->Def.Ne;n++){
        tmp   +=  X->Large.R_SLT_0[i][n]*X->Large.R_SLT_1[j][n] ;
      }
    fprintf(fp," %d  %d %.10lf \n",i,j,tmp);
    }
  }
  */
  fclose(fp);
}

void cal_cisajs(struct BindStruct *X){
    time_t start,end;
    FILE *fp;
    char sdt[256];
    int int_i,int_spin,site_1,site_2,spin_1,spin_2;
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
