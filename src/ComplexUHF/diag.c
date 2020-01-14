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
#include <matrixlapack.h>
#include "../common/setmemory.h"
#include "diag.h"
//#include "mfmemory.c"

void diag(struct BindStruct *X) {

  int int_i, int_j, int_k, int_l;
  double *r;
  double complex **tmp_mat, **vec;
  //double complex tmp_mlt;
  int xMsize;


  xMsize = X->Def.Nsite;

  tmp_mat = cd_2d_allocate(2 * xMsize, 2 * xMsize);
  vec = cd_2d_allocate(2 * xMsize, 2 * xMsize);
  r = d_1d_allocate(2 * xMsize);

  for (int_i = 0; int_i < 2 * xMsize; int_i++) {
    for (int_j = 0; int_j < 2 * xMsize; int_j++) {
      tmp_mat[int_i][int_j] = X->Large.Ham[int_i][int_j];
      //printf("Ham: %d %d: %lf %lf\n",int_i,int_j,creal(tmp_mat[int_i][int_j]),cimag(tmp_mat[int_i][int_j]));
//[s] MERGE BY TM
/*
        if(fabs(creal(X->Large.Ham[int_i][int_j]))>0.001){
                  fprintf(stdout, "Debug: mat(%d, %d)=(%lf, %lf)\n", int_i, int_j, creal(X->Large.Ham[int_i][int_j]), cimag(X->Large.Ham[int_i][int_j]));
        }
*/
//[e] MERGE BY TM
    }
  }

  ZHEEVall(2 * xMsize, tmp_mat, r, vec);
//[e]check
  for (int_k = 0; int_k < 2 * xMsize; int_k++) {
    X->Large.EigenValues[int_k] = r[int_k];
    //fprintf(stdout, "Debug: Eigen[%d]=%lf\n", int_k, X->Large.EigenValues[int_k]);
  }
  //For zero-temperature to generate pair-orbitals
  for (int_k = 0; int_k < X->Def.Nsize; int_k++) {
    for (int_l = 0; int_l < 2 * xMsize; int_l++) {
      //X->Large.R_SLT[int_l][int_k] = conj(vec[int_k][int_l]);
//[s] MERGE BY TM
      X->Large.R_SLT[int_l][int_k] = conj(vec[int_k][int_l]); //original
      X->Large.L_SLT[int_k][int_l] = (vec[int_k][int_l]);     //original
      // R_SLT = U^{*}
      // L_SLT = U^{T}
//[e] MERGE BY TM
    }
  }
//	exit(1);
  free_cd_2d_allocate(tmp_mat);
  free_cd_2d_allocate(vec);
  free_d_1d_allocate(r);
}
