#include "cal_energy.h"

void cal_energy(struct BindStruct *X){

    //time_t start,end;
    //FILE *fp;
    //char sdt[256];
    int int_i,int_j,int_k,site_1,site_2;
    double tmp,charge_1,charge_2;
    double E_band,E_CoulombIntra,E_CoulombInter,E_HundCoupling;
    double E_ExchangeCoupling,E_PairHopping;
    double num,mix,tmp_num_r,tmp_num_i;
    int    xMsize;
    int    u_site_1,u_site_2;
    int    d_site_1,d_site_2;
    int    Ns;

 
    xMsize = X->Def.Nsite;  
    Ns     = X->Def.Nsite;  

    E_band = 0.0;
    for(int_k = 0; int_k < X->Def.Ne*2; int_k++){
      E_band += X->Large.EigenValues[int_k];
    }
    printf("E_band=%lf \n",E_band);
    /*Intra U energy*/
    E_CoulombIntra = 0.0;
    for(int_i =0; int_i < X->Def.NCoulombIntra; int_i++){
      site_1  =  X->Def.CoulombIntra[int_i][0];
      tmp     =  X->Def.ParaCoulombIntra[int_i];

      u_site_1  = site_1+0*Ns;
      d_site_1  = site_1+1*Ns;
      //printf("site_1=%d %lf %lf\n",site_1,tmp,X->Large.G[0][site_1][site_1]);
      E_CoulombIntra += -1.0*tmp*X->Large.G[u_site_1][u_site_1]*X->Large.G[d_site_1][d_site_1];
#if Fock==1
      /*Off-Diagonal Fock term*/
      E_CoulombIntra += 1.0*tmp*X->Large.G[u_site_1][d_site_1]*X->Large.G[d_site_1][u_site_1];
#endif
    }
    //printf("E_CoulombIntra=%lf \n",E_CoulombIntra);
    /*Inter U energy*/
    E_CoulombInter = 0.0;
    for(int_i =0; int_i < X->Def.NCoulombInter; int_i++){
      site_1    =  X->Def.CoulombInter[int_i][0];
      site_2    =  X->Def.CoulombInter[int_i][1];
      tmp       =  X->Def.ParaCoulombInter[int_i];

      u_site_1  = site_1+0*Ns;
      d_site_1  = site_1+1*Ns;

      u_site_2  = site_2+0*Ns;
      d_site_2  = site_2+1*Ns;
      //printf("site_1=%d %lf %lf\n",site_1,tmp,X->Large.G[0][site_1][site_1]);
      charge_1  =  X->Large.G[u_site_1][u_site_1]+X->Large.G[d_site_1][d_site_1];
      charge_2  =  X->Large.G[u_site_2][u_site_2]+X->Large.G[d_site_2][d_site_2];
      E_CoulombInter += -1.0*tmp*charge_1*charge_2;

      /*Diagonal Fock term*/
#if Fock==1
      //printf("Diagonal Fock 1\n");
      E_CoulombInter += 1.0*tmp*X->Large.G[u_site_1][u_site_2]*X->Large.G[u_site_2][u_site_1];
      E_CoulombInter += 1.0*tmp*X->Large.G[d_site_1][d_site_2]*X->Large.G[d_site_2][d_site_1];
      /*Off-Diagonal Fock term*/
      E_CoulombInter += 1.0*tmp*X->Large.G[u_site_1][d_site_2]*X->Large.G[d_site_2][u_site_1];
      E_CoulombInter += 1.0*tmp*X->Large.G[d_site_1][u_site_2]*X->Large.G[u_site_2][d_site_1];
#endif
    }
    //printf("E_CoulombInter=%lf \n",E_CoulombInter);
    /*Hund energy*/
    E_HundCoupling = 0.0;
    for(int_i =0; int_i < X->Def.NHundCoupling; int_i++){
      site_1  =  X->Def.HundCoupling[int_i][0];
      site_2  =  X->Def.HundCoupling[int_i][1];
      tmp     =  -X->Def.ParaHundCoupling[int_i];

      u_site_1  = site_1+0*Ns;
      d_site_1  = site_1+1*Ns;

      u_site_2  = site_2+0*Ns;
      d_site_2  = site_2+1*Ns;

      E_HundCoupling += -1.0*tmp*X->Large.G[u_site_1][u_site_1]*X->Large.G[u_site_2][u_site_2];
      E_HundCoupling += -1.0*tmp*X->Large.G[d_site_1][d_site_1]*X->Large.G[d_site_2][d_site_2];

#if Fock==1
      /*Diagonal Fock term*/
      E_HundCoupling += 1.0*tmp*X->Large.G[u_site_1][u_site_2]*X->Large.G[u_site_2][u_site_1];
      E_HundCoupling += 1.0*tmp*X->Large.G[d_site_1][d_site_2]*X->Large.G[d_site_2][d_site_1];
#endif
    }

    /*Exchange energy*/
    E_ExchangeCoupling = 0.0;
    for(int_i =0; int_i < X->Def.NExchangeCoupling; int_i++){
      site_1  =  X->Def.ExchangeCoupling[int_i][0];
      site_2  =  X->Def.ExchangeCoupling[int_i][1];
      tmp     =  X->Def.ParaExchangeCoupling[int_i];

      u_site_1  = site_1+0*Ns;
      d_site_1  = site_1+1*Ns;

      u_site_2  = site_2+0*Ns;
      d_site_2  = site_2+1*Ns;

#if Fock==1
      /*Diagonal Fock term*/
      E_ExchangeCoupling += -1.0*tmp*X->Large.G[u_site_1][u_site_2]*X->Large.G[d_site_2][d_site_1];
      E_ExchangeCoupling += -1.0*tmp*X->Large.G[d_site_1][d_site_2]*X->Large.G[u_site_2][u_site_1];

      /*Off-Diagonal Fock term*/
      E_ExchangeCoupling += 1.0*tmp*X->Large.G[u_site_1][d_site_1]*X->Large.G[d_site_2][u_site_2];
      E_ExchangeCoupling += 1.0*tmp*X->Large.G[d_site_1][u_site_1]*X->Large.G[u_site_2][d_site_2];
#endif
    }

    /*PairHopping energy*/
    E_PairHopping = 0.0;
    for(int_i =0; int_i < X->Def.NPairHopping; int_i++){
      site_1  =  X->Def.PairHopping[int_i][0];
      site_2  =  X->Def.PairHopping[int_i][1];
      tmp     =  X->Def.ParaPairHopping[int_i];

      u_site_1  = site_1+0*Ns;
      d_site_1  = site_1+1*Ns;

      u_site_2  = site_2+0*Ns;
      d_site_2  = site_2+1*Ns;
      /*Diagonal Fock term*/
#if Fock==1
      E_PairHopping += -1.0*tmp*X->Large.G[u_site_1][u_site_2]*X->Large.G[d_site_1][d_site_2];
      E_PairHopping +=  1.0*tmp*X->Large.G[u_site_1][d_site_2]*X->Large.G[d_site_1][u_site_2];
#endif
    }
    /*Calculating Total Enegy*/
    X->Phys.energy  = E_band + E_CoulombIntra + E_CoulombInter + E_HundCoupling;
    X->Phys.energy += E_ExchangeCoupling + E_PairHopping;

    mix=X->Def.mix;
    X->Phys.rest = 0.0;
    num=0.0;
    for(int_i=0; int_i < 2*xMsize; int_i++){
      num            += creal(X->Large.G[int_i][int_i]);
      for(int_j=0; int_j < 2*xMsize; int_j++){
        tmp           = cabs(X->Large.G_old[int_i][int_j]-X->Large.G[int_i][int_j]);
        X->Phys.rest += tmp*tmp;
        X->Large.G[int_i][int_j] = X->Large.G_old[int_i][int_j]*(1.0-mix)+mix*X->Large.G[int_i][int_j];
      }
    }
    X->Phys.num = num;
    X->Phys.rest = sqrt(X->Phys.rest)/(2.0*X->Def.Nsite*X->Def.Nsite);

    //for(int_i=0;int_i<xMsize;int_i++){
      int_i           = 1;
      tmp_num_r       = creal(X->Large.G[int_i][int_i+xMsize]);
      tmp_num_i       = cimag(X->Large.G[int_i][int_i+xMsize]);
      printf("int_i=%d tmp_num_r = %lf tmp_num_i = %lf\n",int_i,tmp_num_r,tmp_num_i);
      printf("\n");
      //tmp_num_r       = creal(X->Large.G[int_i+xMsize][int_i]);
      //tmp_num_i       = cimag(X->Large.G[int_i+xMsize][int_i]);
      //printf("int_i=%d tmp_num_r = %lf tmp_num_i = %lf\n",int_i,tmp_num_r,tmp_num_i);
    //} 
    //printf("num=%lf \n",num);
    //printf("rest =%.12lf energy = %lf \n",X->Phys.rest,X->Phys.energy);

}
