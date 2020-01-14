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
/* 2009.07.08*/
#include "../common/setmemory.h"

X.Bind.Large.Ham = cd_2d_allocate(2*X.Bind.Def.Nsite, 2*X.Bind.Def.Nsite);
X.Bind.Large.G = cd_2d_allocate(2*X.Bind.Def.Nsite, 2*X.Bind.Def.Nsite);
X.Bind.Large.G_old = cd_2d_allocate(2*X.Bind.Def.Nsite, 2*X.Bind.Def.Nsite);
X.Bind.Large.R_SLT = cd_2d_allocate(2*X.Bind.Def.Nsite, X.Bind.Def.Nsize);
X.Bind.Large.L_SLT = cd_2d_allocate(X.Bind.Def.Nsize, 2*X.Bind.Def.Nsite);
X.Bind.Large.EigenValues = d_1d_allocate(2*X.Bind.Def.Nsite);
printf("LARGE ALLOCATE FINISH !\n");
 
