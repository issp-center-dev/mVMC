#include "green.h"
#include "matrixlapack.h"
#include "mfmemory.c"

void green(struct BindStruct *X){

    double complex **R_Mat;
    int int_i,int_j;
    int mfint[7],xMsize,Ne;

 //   printf("OK \n");

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
