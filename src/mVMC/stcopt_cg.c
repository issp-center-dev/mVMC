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

// #define _DEBUG_STCOPT_CG
// #define _DEBUG_STCOPT_CG_PRINT_SMAT

inline double xdot(const int n, const double * const p, const double * const q) {
  int i;
  double z=0;
  #pragma loop noalias
  for(i=0;i<n;i++) {
    z += p[i]*q[i];
  }
  return z;
}
extern inline double xdot(const int n, const double * const p, const double * const q);

#define MVMC_SRCG_REAL
#include "stcopt_cg_impl.c"
#undef MVMC_SRCG_REAL
#include "stcopt_cg_impl.c"

int StochasticOptCG(MPI_Comm comm)
{
  int ret=0;
  if(AllComplexFlag==0){
    ret = StochasticOptCG_real(comm);
  }else{
    ret = StochasticOptCG_fcmp(comm);
  }
  return ret;
}
