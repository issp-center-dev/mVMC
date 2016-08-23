#include "output.h"
void output(struct BindStruct *X){
  
  FILE *fp;
  char sdt[256];
  int i,i_max,j,n;
  double tmp;
 
  sprintf(sdt,"%s_result.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  fprintf(fp," energy %.10lf \n",X->Phys.energy);
  fprintf(fp," num    %.10lf \n",X->Phys.num);
  fclose(fp);

  sprintf(sdt,"%s_eigen.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  for(i=0;i< X->Def.Nsite;i++){
    tmp =  X->Large.EigenValues[0][i];
    fprintf(fp," %d  %.10lf \n",i+1,tmp);
  }
  fclose(fp);

  sprintf(sdt,"%s_fij.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  for(i=0;i< X->Def.Nsite;i++){
    for(j=0;j< X->Def.Nsite;j++){
    tmp    =  0.0;
      for(n=0;n< X->Def.Ne;n++){
        tmp   +=  X->Large.R_SLT_0[i][n]*X->Large.R_SLT_1[j][n] ;
      }
    fprintf(fp," %d  %d %.10lf \n",i,j,tmp);
    }
  }
  fclose(fp);

  sprintf(sdt,"%s_gap.dat",X->Def.CDataFileHead);
  fp=fopen(sdt,"w");
  tmp =  X->Large.EigenValues[0][X->Def.Ne-1];
  tmp =  X->Large.EigenValues[0][X->Def.Ne] - tmp;
  fprintf(fp,"  %.10lf \n",tmp);
  fclose(fp);

}
