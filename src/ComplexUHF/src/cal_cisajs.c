#include "cal_cisajs.h"

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
