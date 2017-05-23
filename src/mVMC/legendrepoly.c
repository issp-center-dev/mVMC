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
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Legendre Polynomial
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "legendrepoly.h"
#ifndef _LEGENDREPOLY_SRC
#define _LEGENDREPOLY_SRC

double LegendrePoly(const double x, const int n){
  double P01, P02, P03;
  int i;
  if(n<=0) return 1.0;
  else if(n==1) return x;
  else{
    P01 = 1.0;
    P02 = x;
    for(i=2;i<=n;i++){ /* recurrence relation */
      P03 = 1.0/(1.0*i) * ( (2.0*i-1.0)*x*P02 - (1.0*i-1.0)*P01 );
      P01 = P02;
      P02 = P03;
    }
    return P03;
  };
}

#endif

