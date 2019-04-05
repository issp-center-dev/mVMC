/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/
#include "cal_energy.h"

void cal_energy(struct BindStruct *X) {

  //time_t start,end;
  //FILE *fp;
  //char sdt[256];
  int int_i, int_j, int_k, site_1, site_2, site_3, site_4;
  double tmp, charge_1, charge_2;
  complex double ctmp;
  double E_band, E_CoulombIntra, E_CoulombInter, E_HundCoupling;
  double E_ExchangeCoupling, E_PairHopping, E_InterAll;
  double num, mix;
  int xMsize;
  int u_site_1, u_site_2;
  int d_site_1, d_site_2;
  int s_site_1, s_site_2, s_site_3, s_site_4;
  int spin_1, spin_2, spin_3, spin_4;
  int Ns;


  xMsize = X->Def.Nsite;
  Ns = X->Def.Nsite;

  E_band = 0.0;
  for (int_k = 0; int_k < X->Def.Nsize; int_k++) {
    E_band += X->Large.EigenValues[int_k];
  }
  //printf("E_band=%lf \n",E_band);
  /*Intra U energy*/
  E_CoulombIntra = 0.0;
  for (int_i = 0; int_i < X->Def.NCoulombIntra; int_i++) {
    site_1 = X->Def.CoulombIntra[int_i][0];
    tmp = X->Def.ParaCoulombIntra[int_i];

    u_site_1 = site_1 + 0 * Ns;
    d_site_1 = site_1 + 1 * Ns;
    //printf("site_1=%d %lf %lf\n",site_1,tmp,X->Large.G[0][site_1][site_1]);
    E_CoulombIntra += -1.0 * tmp * X->Large.G[u_site_1][u_site_1] * X->Large.G[d_site_1][d_site_1];
    /*Off-Diagonal Fock term*/
    if (X->Def.iFlg_Fock == 1) {
        E_CoulombIntra += 1.0 * tmp * X->Large.G[u_site_1][d_site_1] * X->Large.G[d_site_1][u_site_1];
    }
  }
  //printf("E_CoulombIntra=%lf \n",E_CoulombIntra);
  /*Inter U energy*/
  E_CoulombInter = 0.0;
  for (int_i = 0; int_i < X->Def.NCoulombInter; int_i++) {
    site_1 = X->Def.CoulombInter[int_i][0];
    site_2 = X->Def.CoulombInter[int_i][1];
    tmp = X->Def.ParaCoulombInter[int_i];

    u_site_1 = site_1 + 0 * Ns;
    d_site_1 = site_1 + 1 * Ns;

    u_site_2 = site_2 + 0 * Ns;
    d_site_2 = site_2 + 1 * Ns;
    //printf("site_1=%d %lf %lf\n",site_1,tmp,X->Large.G[0][site_1][site_1]);
    charge_1 = X->Large.G[u_site_1][u_site_1] + X->Large.G[d_site_1][d_site_1];
    charge_2 = X->Large.G[u_site_2][u_site_2] + X->Large.G[d_site_2][d_site_2];
    E_CoulombInter += -1.0 * tmp * charge_1 * charge_2;

    /*Diagonal Fock term*/
    if (X->Def.iFlg_Fock == 1) {
      //printf("Diagonal Fock 1\n");
      E_CoulombInter += 1.0 * tmp * X->Large.G[u_site_1][u_site_2] * X->Large.G[u_site_2][u_site_1];
      E_CoulombInter += 1.0 * tmp * X->Large.G[d_site_1][d_site_2] * X->Large.G[d_site_2][d_site_1];
      /*Off-Diagonal Fock term*/
      E_CoulombInter += 1.0 * tmp * X->Large.G[u_site_1][d_site_2] * X->Large.G[d_site_2][u_site_1];
      E_CoulombInter += 1.0 * tmp * X->Large.G[d_site_1][u_site_2] * X->Large.G[u_site_2][d_site_1];
    }
  }
  //printf("E_CoulombInter=%lf \n",E_CoulombInter);
  /*Hund energy*/
  E_HundCoupling = 0.0;
  for (int_i = 0; int_i < X->Def.NHundCoupling; int_i++) {
    site_1 = X->Def.HundCoupling[int_i][0];
    site_2 = X->Def.HundCoupling[int_i][1];
    tmp = -X->Def.ParaHundCoupling[int_i];

    u_site_1 = site_1 + 0 * Ns;
    d_site_1 = site_1 + 1 * Ns;

    u_site_2 = site_2 + 0 * Ns;
    d_site_2 = site_2 + 1 * Ns;

    E_HundCoupling += -1.0 * tmp * X->Large.G[u_site_1][u_site_1] * X->Large.G[u_site_2][u_site_2];
    E_HundCoupling += -1.0 * tmp * X->Large.G[d_site_1][d_site_1] * X->Large.G[d_site_2][d_site_2];

    if (X->Def.iFlg_Fock == 1) {
      /*Diagonal Fock term*/
      E_HundCoupling += 1.0 * tmp * X->Large.G[u_site_1][u_site_2] * X->Large.G[u_site_2][u_site_1];
      E_HundCoupling += 1.0 * tmp * X->Large.G[d_site_1][d_site_2] * X->Large.G[d_site_2][d_site_1];
    }
  }

  /*Exchange energy*/
  E_ExchangeCoupling = 0.0;
  if (X->Def.iFlg_Fock == 1) {
    for (int_i = 0; int_i < X->Def.NExchangeCoupling; int_i++) {
      site_1 = X->Def.ExchangeCoupling[int_i][0];
      site_2 = X->Def.ExchangeCoupling[int_i][1];
      tmp = X->Def.ParaExchangeCoupling[int_i];

      u_site_1 = site_1 + 0 * Ns;
      d_site_1 = site_1 + 1 * Ns;

      u_site_2 = site_2 + 0 * Ns;
      d_site_2 = site_2 + 1 * Ns;

      /*Diagonal Fock term*/
      E_ExchangeCoupling += -1.0 * tmp * X->Large.G[u_site_1][u_site_2] * X->Large.G[d_site_2][d_site_1];
      E_ExchangeCoupling += -1.0 * tmp * X->Large.G[d_site_1][d_site_2] * X->Large.G[u_site_2][u_site_1];

      /*Off-Diagonal Fock term*/
      E_ExchangeCoupling += 1.0 * tmp * X->Large.G[u_site_1][d_site_1] * X->Large.G[d_site_2][u_site_2];
      E_ExchangeCoupling += 1.0 * tmp * X->Large.G[d_site_1][u_site_1] * X->Large.G[u_site_2][d_site_2];
    }

    /*PairHopping energy*/
    E_PairHopping = 0.0;
    for (int_i = 0; int_i < X->Def.NPairHopping; int_i++) {
      site_1 = X->Def.PairHopping[int_i][0];
      site_2 = X->Def.PairHopping[int_i][2];

      tmp = X->Def.ParaPairHopping[int_i];

      u_site_1 = site_1 + 0 * Ns;
      d_site_1 = site_1 + 1 * Ns;

      u_site_2 = site_2 + 0 * Ns;
      d_site_2 = site_2 + 1 * Ns;
      /*Diagonal Fock term*/
      E_PairHopping += -1.0 * tmp * X->Large.G[u_site_1][u_site_2] * X->Large.G[d_site_1][d_site_2];
      E_PairHopping += 1.0 * tmp * X->Large.G[u_site_1][d_site_2] * X->Large.G[d_site_1][u_site_2];
    }

    /*InterAll energy */
    E_InterAll = 0.0;
    for (int_i = 0; int_i < X->Def.NInterAll; int_i++) {
      site_1 = X->Def.InterAll[int_i][0];
      spin_1 = X->Def.InterAll[int_i][1];
      site_2 = X->Def.InterAll[int_i][2];
      spin_2 = X->Def.InterAll[int_i][3];
      site_3 = X->Def.InterAll[int_i][4];
      spin_3 = X->Def.InterAll[int_i][5];
      site_4 = X->Def.InterAll[int_i][6];
      spin_4 = X->Def.InterAll[int_i][7];
      ctmp = X->Def.ParaInterAll[int_i];

      s_site_1 = site_1 + spin_1 * Ns;
      s_site_2 = site_2 + spin_2 * Ns;
      s_site_3 = site_3 + spin_3 * Ns;
      s_site_4 = site_4 + spin_4 * Ns;
      E_InterAll += ctmp * X->Large.G[s_site_1][s_site_2] * X->Large.G[s_site_3][s_site_4];
      E_InterAll -= ctmp * X->Large.G[s_site_1][s_site_3] * X->Large.G[s_site_2][s_site_4];
    }
  }
  /*Calculating Total Enegy*/
  X->Phys.energy = E_band + E_CoulombIntra + E_CoulombInter + E_HundCoupling;
  X->Phys.energy += E_ExchangeCoupling + E_PairHopping + E_InterAll;

  mix = X->Def.mix;
  X->Phys.rest = 0.0;
  num = 0.0;
  for (int_i = 0; int_i < 2 * xMsize; int_i++) {
    num += creal(X->Large.G[int_i][int_i]);
    for (int_j = 0; int_j < 2 * xMsize; int_j++) {
      tmp = cabs(X->Large.G_old[int_i][int_j] - X->Large.G[int_i][int_j]);
      X->Phys.rest += tmp * tmp;
      X->Large.G[int_i][int_j] = X->Large.G_old[int_i][int_j] * (1.0 - mix) + mix * X->Large.G[int_i][int_j];
    }
  }
  X->Phys.num = num;
  X->Phys.rest = sqrt(X->Phys.rest) / (2.0 * X->Def.Nsite * X->Def.Nsite);

  //for(int_i=0;int_i<xMsize;int_i++){
  //int_i           = 1;
  //tmp_num_r       = creal(X->Large.G[int_i][int_i+xMsize]);
  //tmp_num_i       = cimag(X->Large.G[int_i][int_i+xMsize]);
  //printf("int_i=%d tmp_num_r = %lf tmp_num_i = %lf\n",int_i,tmp_num_r,tmp_num_i);
  //printf("\n");
  //tmp_num_r       = creal(X->Large.G[int_i+xMsize][int_i]);
  //tmp_num_i       = cimag(X->Large.G[int_i+xMsize][int_i]);
  //printf("int_i=%d tmp_num_r = %lf tmp_num_i = %lf\n",int_i,tmp_num_r,tmp_num_i);
  //}
  //printf("num=%lf \n",num);
}
