#include "green.h"
#include "matrixlapack.h"
#include "mfmemory.c"

void green(struct BindStruct *X){

    time_t start,end;
    FILE *fp;
    char sdt[256],jobz,uplo;
    double alpha,beta;
    int    lda;
    double tmp_0,tmp_1;
    double complex **R_Mat;
    int int_i,int_j,int_k,int_l,site_1,site_2;
    int site_i,site_j,int_i_A,int_j_A,int_i_B,int_j_B;
    int mfint[7],xMsize,Ne;

 
    printf("OK \n");

    xMsize = X->Def.Nsite;  
    Ne     = X->Def.Ne;  

    c_malloc2(R_Mat,2*xMsize,2*xMsize);

    for(int_i=0; int_i < 2*xMsize; int_i++){
      for(int_j=0; int_j < 2*xMsize; int_j++){
        X->Large.G_old[int_i][int_j] = X->Large.G[int_i][int_j];
        R_Mat[int_i][int_j] = 0.0;
      }
    }

    cmp_MMProd(2*xMsize,2*Ne,X->Large.R_SLT,X->Large.L_SLT,R_Mat);

    for(int_i = 0; int_i < 2*xMsize; int_i++){
      for(int_j = 0; int_j < 2*xMsize; int_j++){
        X->Large.G[int_i][int_j] = R_Mat[int_i][int_j];
      }
    }

    c_free2(R_Mat,2*xMsize,2*xMsize);
}
