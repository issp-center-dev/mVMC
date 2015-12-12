void green(struct BindStruct *X){

    time_t start,end;
    FILE *fp;
    char sdt[256],jobz,uplo;
    double alpha,beta;
    int    lda;
    double tmp_0,tmp_1;
    double **R_Mat_0;
    double **R_Mat_1;
    int int_i,int_j,int_k,int_l,site_1,site_2,int_spin;
    int site_i,site_j,int_i_A,int_j_A,int_i_B,int_j_B;
    int mfint[7],xMsize,Ne;

 
    printf("OK \n");

    xMsize = X->Def.Nsite;  
    Ne     = X->Def.Ne;  

    d_malloc2(R_Mat_0,xMsize,xMsize);
    d_malloc2(R_Mat_1,xMsize,xMsize);

    for(int_spin = 0; int_spin < 2; int_spin++){
      for(int_i=0; int_i < xMsize; int_i++){
        for(int_j=0; int_j < xMsize; int_j++){
          X->Large.G_old[int_spin][int_i][int_j] = X->Large.G[int_spin][int_i][int_j];
          R_Mat_0[int_i][int_j] = 0.0;
          R_Mat_1[int_i][int_j] = 0.0;
        }
      }
    }

    MMProd(xMsize,Ne,X->Large.R_SLT_0,X->Large.L_SLT_0,R_Mat_0);
    MMProd(xMsize,Ne,X->Large.R_SLT_1,X->Large.L_SLT_1,R_Mat_1);

    for(int_i = 0; int_i < xMsize; int_i++){
      for(int_j = 0; int_j < xMsize; int_j++){
        X->Large.G[0][int_i][int_j] = R_Mat_0[int_i][int_j];
        X->Large.G[1][int_i][int_j] = R_Mat_1[int_i][int_j];
      }
    }


    d_free2(R_Mat_0,xMsize,xMsize);
    d_free2(R_Mat_1,xMsize,xMsize);
}
