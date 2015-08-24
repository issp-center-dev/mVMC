void calenergy(struct BindStruct *X){

    time_t start,end;
    FILE *fp;
    char sdt[256];
    int int_i,int_j,int_k,int_l,site_1,site_2,int_spin;
    double tmp,charge_1,charge_2;
    double E_band,E_CoulombIntra,E_CoulombInter,E_HundCoupling;
    double E_ExchangeCoupling,E_PairHopping,E_InterAll;
    double rest,num,mix;
    int site_i,site_j,int_i_A,int_j_A,int_i_B,int_j_B;
    int site_3,site_4;
    int sigma,tau;
    int    xMsize;

 
    xMsize = X->Def.Nsite;  

    E_band = 0.0;
    for(int_spin = 0;int_spin < 2; int_spin++){
      for(int_k = 0; int_k < X->Def.Ne; int_k++){
        E_band += X->Large.EigenValues[int_spin][int_k];
      }
    }
    //printf("E_band=%lf \n",E_band);
    /*Intra U energy*/
    E_CoulombIntra = 0.0;
    for(int_i =0; int_i < X->Def.NCoulombIntra; int_i++){
      site_1  =  X->Def.CoulombIntra[int_i][0];
      tmp     =  X->Def.ParaCoulombIntra[int_i];
      //printf("site_1=%d %lf %lf\n",site_1,tmp,X->Large.G[0][site_1][site_1]);
      E_CoulombIntra += -1.0*tmp*X->Large.G[0][site_1][site_1]*X->Large.G[1][site_1][site_1];
    }
    //printf("E_CoulombIntra=%lf \n",E_CoulombIntra);
    /*Inter U energy*/
    E_CoulombInter = 0.0;
    for(int_i =0; int_i < X->Def.NCoulombInter; int_i++){
      site_1    =  X->Def.CoulombInter[int_i][0];
      site_2    =  X->Def.CoulombInter[int_i][1];
      tmp       =  X->Def.ParaCoulombInter[int_i];
      //printf("site_1=%d %lf %lf\n",site_1,tmp,X->Large.G[0][site_1][site_1]);
      charge_1  =  X->Large.G[0][site_1][site_1]+X->Large.G[1][site_1][site_1];
      charge_2  =  X->Large.G[0][site_2][site_2]+X->Large.G[1][site_2][site_2];
      E_CoulombInter += -1.0*tmp*charge_1*charge_2;

      /*Diagonal Fock term*/
#if Fock==1
      //printf("Diagonal Fock 1\n");
      E_CoulombInter += 1.0*tmp*X->Large.G[0][site_1][site_2]*X->Large.G[0][site_2][site_1];
      E_CoulombInter += 1.0*tmp*X->Large.G[1][site_1][site_2]*X->Large.G[1][site_2][site_1];
#endif
    }
    //printf("E_CoulombInter=%lf \n",E_CoulombInter);
    /*Hund energy*/
    E_HundCoupling = 0.0;
    for(int_i =0; int_i < X->Def.NHundCoupling; int_i++){
      site_1  =  X->Def.HundCoupling[int_i][0];
      site_2  =  X->Def.HundCoupling[int_i][1];
      tmp     =  -X->Def.ParaHundCoupling[int_i];

      E_HundCoupling += -1.0*tmp*X->Large.G[0][site_1][site_1]*X->Large.G[0][site_2][site_2];
      E_HundCoupling += -1.0*tmp*X->Large.G[1][site_1][site_1]*X->Large.G[1][site_2][site_2];
// add 130725
#if Fock==1
      E_HundCoupling += 1.0*tmp*X->Large.G[0][site_1][site_2]*X->Large.G[0][site_2][site_1];
      E_HundCoupling += 1.0*tmp*X->Large.G[1][site_1][site_2]*X->Large.G[1][site_2][site_1];
#endif
    }

    /*Exchange energy*/
    E_ExchangeCoupling = 0.0;
    for(int_i =0; int_i < X->Def.NExchangeCoupling; int_i++){
      site_1  =  X->Def.ExchangeCoupling[int_i][0];
      site_2  =  X->Def.ExchangeCoupling[int_i][1];
      tmp     =  X->Def.ParaExchangeCoupling[int_i];

      /*Diagonal Fock term*/
#if Fock==1
      E_ExchangeCoupling += -1.0*tmp*X->Large.G[0][site_1][site_2]*X->Large.G[1][site_2][site_1];
      E_ExchangeCoupling += -1.0*tmp*X->Large.G[1][site_1][site_2]*X->Large.G[0][site_2][site_1];
#endif
    }

    /*PairHopping energy*/
    E_PairHopping = 0.0;
    for(int_i =0; int_i < X->Def.NPairHopping; int_i++){
      site_1  =  X->Def.PairHopping[int_i][0];
      site_2  =  X->Def.PairHopping[int_i][1];
      tmp     =  X->Def.ParaPairHopping[int_i];

      /*Diagonal Fock term*/
#if Fock==1
      E_PairHopping += -1.0*tmp*X->Large.G[0][site_1][site_2]*X->Large.G[1][site_1][site_2];
#endif
    }
    /*InterAll energy*/
    E_InterAll = 0.0;
    for(int_i =0; int_i < X->Def.NInterAll; int_i++){
      site_1  =  X->Def.InterAll[int_i][0];
      site_2  =  X->Def.InterAll[int_i][1];
      sigma   =  X->Def.InterAll[int_i][2];
      site_3  =  X->Def.InterAll[int_i][3];
      site_4  =  X->Def.InterAll[int_i][4];
      tau     =  X->Def.InterAll[int_i][5];
      tmp     =  X->Def.ParaInterAll[int_i];
      /*Diagonal Fock term*/
#if Fock==1
      E_InterAll += -1.0*tmp*X->Large.G[sigma][site_1][site_2]*X->Large.G[tau][site_3][site_4];
      if(sigma==tau){
        E_InterAll += 1.0*tmp*X->Large.G[sigma][site_1][site_4]*X->Large.G[sigma][site_3][site_2];
      }
#endif
    }
    /*Calculating Total Enegy*/
    //printf ("%lf %lf %lf %lf %lf %lf\n",E_band ,E_CoulombIntra,E_CoulombInter,E_HundCoupling,E_ExchangeCoupling,E_PairHopping);
    X->Phys.energy  = E_band + E_CoulombIntra + E_CoulombInter + E_HundCoupling;
    X->Phys.energy += E_ExchangeCoupling + E_PairHopping+ E_InterAll;

    mix=X->Def.mix;
    X->Phys.rest = 0.0;
    num=0.0;
    for(int_spin = 0; int_spin < 2; int_spin++){
      for(int_i=0; int_i < xMsize; int_i++){
         num          += X->Large.G[int_spin][int_i][int_i];
        for(int_j=0; int_j < xMsize; int_j++){
          tmp           = X->Large.G_old[int_spin][int_i][int_j]-X->Large.G[int_spin][int_i][int_j];
          X->Phys.rest += tmp*tmp;
          X->Large.G[int_spin][int_i][int_j] = X->Large.G_old[int_spin][int_i][int_j]*(1.0-mix)+mix*X->Large.G[int_spin][int_i][int_j];
        }
      }
    }
    X->Phys.num = num;
    X->Phys.rest = sqrt(X->Phys.rest)/(2.0*X->Def.Nsite*X->Def.Nsite);
    //printf("num=%lf \n",num);
    //printf("rest =%.12lf energy = %lf \n",X->Phys.rest,X->Phys.energy);

}
