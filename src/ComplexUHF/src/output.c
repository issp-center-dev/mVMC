#include "output.h"

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

}
