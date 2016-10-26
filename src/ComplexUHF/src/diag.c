/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

his program is developed based on the mVMC-mini program
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
#include "diag.h"
#include "mfmemory.c"

void diag(struct BindStruct *X){

    int int_i,int_j,int_k,int_l;
    double *r;
    double complex **tmp_mat,**vec;
    int mfint[7],xMsize;

 
    xMsize = X->Def.Nsite;  
    c_malloc2(tmp_mat,2*xMsize,2*xMsize);
    c_malloc2(vec,2*xMsize,2*xMsize);
    d_malloc1(r,2*xMsize);

    for(int_i = 0;int_i < 2*xMsize; int_i++){
      for(int_j = 0;int_j < 2*xMsize; int_j++){
        tmp_mat[int_i][int_j] = X->Large.Ham[int_i][int_j];
        //          fprintf(stdout, "Debug: mat(%d, %d)=(%lf, %lf)\n", int_i, int_j, creal(X->Large.Ham[int_i][int_j]), cimag(X->Large.Ham[int_i][int_j]));
      }
    }
    ZHEEVall(2*xMsize,tmp_mat,r,vec);
    for(int_k = 0; int_k < 2*xMsize; int_k++){
      X->Large.EigenValues[int_k] = r[int_k];
      //fprintf(stdout, "Debug: Eigen[%d]=%lf\n", int_k, X->Large.EigenValues[int_k]);
    }

    //For zero-temperature to generate pair-orbitals
    for(int_k = 0; int_k < 2*X->Def.Ne; int_k++){
      for(int_l = 0; int_l < 2*xMsize; int_l++){
        X->Large.R_SLT[int_l][int_k] = conj(vec[int_k][int_l]);
        X->Large.L_SLT[int_k][int_l] = (vec[int_k][int_l]);
      }
    }

    c_free2(tmp_mat,2*xMsize,2*xMsize);
    c_free2(vec,2*xMsize,2*xMsize);
    d_free1(r,2*xMsize);
}
