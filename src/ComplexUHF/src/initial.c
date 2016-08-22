#include "initial.h"

void initial(struct BindStruct *X){

    int int_i,int_j;
    int spin_0,spin_1;
    int site_0,site_1;
    int t_site_0,t_site_1;
    int Ns;
    double theta;
    double complex tmp;

    Ns = X->Def.Nsite;

    for(int_i=0; int_i < 2*X->Def.Nsite; int_i++){
      for(int_j=0; int_j < 2*X->Def.Nsite; int_j++){
        X->Large.G[int_i][int_j]   = 0.0;
      }
    }

    for(int_i=0; int_i < X->Def.NInitial; int_i++){
      site_0  = X->Def.Initial[int_i][0]; 
      site_1  = X->Def.Initial[int_i][1]; 
      spin_0  = X->Def.Initial[int_i][2]; 
      spin_1  = X->Def.Initial[int_i][3]; 
      
      theta   = X->Def.ParaInitial_theta[int_i];
      tmp     = X->Def.ParaInitial[int_i]*(cos(theta)+I*sin(theta));
      //tmp     = I;
      t_site_0 = site_0+spin_0*Ns;
      t_site_1 = site_1+spin_1*Ns;
      printf(" Initial: %4d %4d %lf %lf \n",t_site_0,t_site_1,cabs(tmp),carg(tmp));
      X->Large.G[t_site_0][t_site_1]  = tmp ;
    }
}
