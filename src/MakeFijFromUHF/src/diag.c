#include "diag.h"
#include "mfmemory.c"

void diag(struct BindStruct *X){

    time_t start,end;
    FILE *fp;
    char sdt[256];
    int int_i,int_j,int_k,int_l,site_1,site_2,int_spin;
    double tmp,charge;
    double **tmp_mat,*r,**vec;
    int mfint[7],xMsize;

 
    xMsize = X->Def.Nsite;  
    d_malloc2(tmp_mat,xMsize,xMsize);
    d_malloc2(vec,xMsize,xMsize);
    d_malloc1(r,xMsize);

    for(int_spin = 0;int_spin < 2; int_spin++){
      for(int_i = 0;int_i < xMsize; int_i++){
        for(int_j = 0;int_j < xMsize; int_j++){
           tmp_mat[int_i][int_j] = X->Large.Ham[int_spin][int_i][int_j];
        }
      }
      DSEVvector(xMsize,tmp_mat,r,vec);
      for(int_k = 0; int_k < X->Def.Nsite; int_k++){
        X->Large.EigenValues[int_spin][int_k] = r[int_k];
      }
      if(int_spin==0){
        for(int_k = 0; int_k < X->Def.Ne; int_k++){
          for(int_l = 0; int_l < xMsize; int_l++){
            X->Large.R_SLT_0[int_l][int_k] = vec[int_k][int_l];
            X->Large.L_SLT_0[int_k][int_l] = vec[int_k][int_l];
          }
        }
      }else{
        for(int_k = 0; int_k < X->Def.Ne; int_k++){
          for(int_l = 0; int_l < xMsize; int_l++){
            X->Large.R_SLT_1[int_l][int_k] = vec[int_k][int_l];
            X->Large.L_SLT_1[int_k][int_l] = vec[int_k][int_l];
          }
        }
      }
    }
    d_free2(tmp_mat,xMsize,xMsize);
    d_free2(vec,xMsize,xMsize);
    d_free1(r,xMsize);
}
