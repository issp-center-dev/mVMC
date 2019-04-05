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
#include "makeham.h"

void makeham(struct BindStruct *X) {

  int int_i, int_j, site_1, site_2, site_3, site_4;
  int int_spin1, int_spin2, int_spin3, int_spin4;
  int s_site_1, s_site_2, s_site_3, s_site_4;
  int t_site_1, t_site_2;
  int u_site_1, u_site_2;
  int d_site_1, d_site_2;
  int Ns;
  double charge;
  double complex tmp;

  Ns = X->Def.Nsite;
  /* Initialize */
  for (int_i = 0; int_i < 2 * X->Def.Nsite; int_i++) {
    for (int_j = 0; int_j < 2 * X->Def.Nsite; int_j++) {
      X->Large.Ham[int_i][int_j] = 0.0;
    }
  }

  /*Transfer input*/
  for (int_i = 0; int_i < X->Def.NTransfer; int_i++) {
    site_1 = X->Def.Transfer[int_i][0];
    int_spin1 = X->Def.Transfer[int_i][1];
    site_2 = X->Def.Transfer[int_i][2];
    int_spin2 = X->Def.Transfer[int_i][3];

    tmp = -X->Def.ParaTransfer[int_i];
    t_site_1 = site_1 + int_spin1 * Ns;
    t_site_2 = site_2 + int_spin2 * Ns;

    X->Large.Ham[t_site_1][t_site_2] += tmp;
  }
  /*Intra U input*/
  for (int_i = 0; int_i < X->Def.NCoulombIntra; int_i++) {
    site_1 = X->Def.CoulombIntra[int_i][0];
    tmp = X->Def.ParaCoulombIntra[int_i];

    u_site_1 = site_1 + 0 * Ns;
    d_site_1 = site_1 + 1 * Ns;

    X->Large.Ham[u_site_1][u_site_1] += tmp * X->Large.G[d_site_1][d_site_1];
    X->Large.Ham[d_site_1][d_site_1] += tmp * X->Large.G[u_site_1][u_site_1];
//#if Fock==1
    /*Off-Diagonal Fock term*/
    X->Large.Ham[u_site_1][d_site_1] += -1.0 * tmp * X->Large.G[d_site_1][u_site_1];
    X->Large.Ham[d_site_1][u_site_1] += -1.0 * tmp * X->Large.G[u_site_1][d_site_1];
//#endif
  }
  /*Inter U input*/
  for (int_i = 0; int_i < X->Def.NCoulombInter; int_i++) {
    site_1 = X->Def.CoulombInter[int_i][0];
    site_2 = X->Def.CoulombInter[int_i][1];
    tmp = X->Def.ParaCoulombInter[int_i];

    u_site_1 = site_1 + 0 * Ns;
    d_site_1 = site_1 + 1 * Ns;

    u_site_2 = site_2 + 0 * Ns;
    d_site_2 = site_2 + 1 * Ns;

    charge = X->Large.G[u_site_2][u_site_2] + X->Large.G[d_site_2][d_site_2];
    X->Large.Ham[u_site_1][u_site_1] += tmp * charge;
    X->Large.Ham[d_site_1][d_site_1] += tmp * charge;

    charge = X->Large.G[u_site_1][u_site_1] + X->Large.G[d_site_1][d_site_1];
    X->Large.Ham[u_site_2][u_site_2] += tmp * charge;
    X->Large.Ham[d_site_2][d_site_2] += tmp * charge;
    if (X->Def.iFlg_Fock == 1) {
      /*Diagonal Fock term*/
      X->Large.Ham[u_site_1][u_site_2] += -tmp * X->Large.G[u_site_2][u_site_1];
      X->Large.Ham[u_site_2][u_site_1] += -tmp * X->Large.G[u_site_1][u_site_2];
      X->Large.Ham[d_site_1][d_site_2] += -tmp * X->Large.G[d_site_2][d_site_1];
      X->Large.Ham[d_site_2][d_site_1] += -tmp * X->Large.G[d_site_1][d_site_2];
      /*Off-Diagonal Fock term*/
      X->Large.Ham[u_site_1][d_site_2] += -tmp * X->Large.G[d_site_2][u_site_1];
      X->Large.Ham[d_site_2][u_site_1] += -tmp * X->Large.G[u_site_1][d_site_2];

      X->Large.Ham[u_site_2][d_site_1] += -tmp * X->Large.G[d_site_1][u_site_2];
      X->Large.Ham[d_site_1][u_site_2] += -tmp * X->Large.G[u_site_2][d_site_1];
    }
  }
  /*Hund input*/
  for (int_i = 0; int_i < X->Def.NHundCoupling; int_i++) {
    site_1 = X->Def.HundCoupling[int_i][0];
    site_2 = X->Def.HundCoupling[int_i][1];
    tmp = -X->Def.ParaHundCoupling[int_i];

    u_site_1 = site_1 + 0 * Ns;
    d_site_1 = site_1 + 1 * Ns;

    u_site_2 = site_2 + 0 * Ns;
    d_site_2 = site_2 + 1 * Ns;

    X->Large.Ham[u_site_1][u_site_1] += tmp * X->Large.G[u_site_2][u_site_2];
    X->Large.Ham[d_site_1][d_site_1] += tmp * X->Large.G[d_site_2][d_site_2];

    X->Large.Ham[u_site_2][u_site_2] += tmp * X->Large.G[u_site_1][u_site_1];
    X->Large.Ham[d_site_2][d_site_2] += tmp * X->Large.G[d_site_1][d_site_1];
    if (X->Def.iFlg_Fock == 1) {
      /*Diagonal Fock term*/
      X->Large.Ham[u_site_1][u_site_2] += -tmp * X->Large.G[u_site_2][u_site_1];
      X->Large.Ham[u_site_2][u_site_1] += -tmp * X->Large.G[u_site_1][u_site_2];

      X->Large.Ham[d_site_1][d_site_2] += -tmp * X->Large.G[d_site_2][d_site_1];
      X->Large.Ham[d_site_2][d_site_1] += -tmp * X->Large.G[d_site_1][d_site_2];
    }
  }
  if (X->Def.iFlg_Fock == 1) {

    /*Exchange input*/
    for (int_i = 0; int_i < X->Def.NExchangeCoupling; int_i++) {
      site_1 = X->Def.ExchangeCoupling[int_i][0];
      site_2 = X->Def.ExchangeCoupling[int_i][1];
      tmp = X->Def.ParaExchangeCoupling[int_i];

      u_site_1 = site_1 + 0 * Ns;
      d_site_1 = site_1 + 1 * Ns;

      u_site_2 = site_2 + 0 * Ns;
      d_site_2 = site_2 + 1 * Ns;

      /*Diagonal Fock term*/
      X->Large.Ham[u_site_1][u_site_2] += tmp * X->Large.G[d_site_2][d_site_1];
      X->Large.Ham[d_site_2][d_site_1] += tmp * X->Large.G[u_site_1][u_site_2];

      X->Large.Ham[d_site_1][d_site_2] += tmp * X->Large.G[u_site_2][u_site_1];
      X->Large.Ham[u_site_2][u_site_1] += tmp * X->Large.G[d_site_1][d_site_2];

      /*Off-Diagonal Fock term*/
      X->Large.Ham[u_site_1][d_site_1] += -tmp * X->Large.G[d_site_2][u_site_2];
      X->Large.Ham[u_site_2][d_site_2] += -tmp * X->Large.G[d_site_1][u_site_1];

      X->Large.Ham[d_site_1][u_site_1] += -tmp * X->Large.G[u_site_2][d_site_2];
      X->Large.Ham[d_site_2][u_site_2] += -tmp * X->Large.G[u_site_1][d_site_1];
    }
    /*PariHopping input*/
    for (int_i = 0; int_i < X->Def.NPairHopping; int_i++) {
      site_1 = X->Def.PairHopping[int_i][0];
      site_2 = X->Def.PairHopping[int_i][1];
      tmp = X->Def.ParaPairHopping[int_i];
      u_site_1 = site_1 + 0 * Ns;
      d_site_1 = site_1 + 1 * Ns;

      u_site_2 = site_2 + 0 * Ns;
      d_site_2 = site_2 + 1 * Ns;

      /*Diagonal Fock term*/
      X->Large.Ham[u_site_1][u_site_2] += tmp * X->Large.G[d_site_1][d_site_2];
      X->Large.Ham[d_site_1][d_site_2] += tmp * X->Large.G[u_site_1][u_site_2];

      /*Off-Diagonal Fock term*/
      X->Large.Ham[u_site_1][d_site_2] += -tmp * X->Large.G[d_site_1][u_site_2];
      X->Large.Ham[d_site_1][u_site_2] += -tmp * X->Large.G[u_site_1][d_site_2];
    }

    /*InterAll input*/
    for (int_i = 0; int_i < X->Def.NInterAll; int_i++) {
      site_1 = X->Def.InterAll[int_i][0];
      int_spin1 = X->Def.InterAll[int_i][1];
      site_2 = X->Def.InterAll[int_i][2];
      int_spin2 = X->Def.InterAll[int_i][3];
      site_3 = X->Def.InterAll[int_i][4];
      int_spin3 = X->Def.InterAll[int_i][5];
      site_4 = X->Def.InterAll[int_i][6];
      int_spin4 = X->Def.InterAll[int_i][7];
      tmp = X->Def.ParaInterAll[int_i];

      s_site_1 = site_1 + int_spin1 * Ns;
      s_site_2 = site_2 + int_spin2 * Ns;
      s_site_3 = site_3 + int_spin3 * Ns;
      s_site_4 = site_4 + int_spin4 * Ns;

      /*Diagonal Fock term*/
      X->Large.Ham[s_site_1][s_site_2] += tmp * X->Large.G[s_site_3][s_site_4];
      X->Large.Ham[s_site_3][s_site_4] += tmp * X->Large.G[s_site_1][s_site_2];

      /*Off-Diagonal Fock term*/
      X->Large.Ham[s_site_1][s_site_4] += -tmp * X->Large.G[s_site_2][s_site_3];
      X->Large.Ham[s_site_2][s_site_3] += -tmp * X->Large.G[s_site_1][s_site_4];

      if (s_site_2 == s_site_3){
        X->Large.Ham[s_site_1][s_site_4] += tmp;
      }

    }
  }
}
