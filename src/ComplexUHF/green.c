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
#include "green.h"
#include "matrixlapack.h"
#include "../common/setmemory.h"

void green(struct BindStruct *X) {

  double complex **R_Mat;
  int int_i, int_j;
  int  xMsize;

  xMsize = X->Def.Nsite;

  R_Mat = cd_2d_allocate(2 * xMsize, 2 * xMsize);

  for (int_i = 0; int_i < 2 * xMsize; int_i++) {
    for (int_j = 0; int_j < 2 * xMsize; int_j++) {
      X->Large.G_old[int_i][int_j] = X->Large.G[int_i][int_j];
      R_Mat[int_i][int_j] = 0.0;
    }
  }

  cmp_MMProd(2 * xMsize, X->Def.Nsize, X->Large.R_SLT, X->Large.L_SLT, R_Mat);

  for (int_i = 0; int_i < 2 * xMsize; int_i++) {
    for (int_j = 0; int_j < 2 * xMsize; int_j++) {
      X->Large.G[int_i][int_j] = R_Mat[int_i][int_j];
      //fprintf(stdout, "DEBUG: X->Large.G[%d][%d]=%lf, %lf \n", int_i, int_j,creal(X->Large.G[int_i][int_j]), cimag(X->Large.G[int_i][int_j]));
    }
  }
  free_cd_2d_allocate(R_Mat);
}
