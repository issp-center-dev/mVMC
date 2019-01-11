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
#include "struct.h"
#include "../common/setmemory.h"
void setmem(struct EDMainCalStruct* X) {
    X->Bind.Def.Tpow = i_1d_allocate(2 * X->Bind.Def.Nsite + 1);
    X->Bind.Def.LocSpn = i_1d_allocate(X->Bind.Def.Nsite);
    X->Bind.Def.All_pair = i_2d_allocate(X->Bind.Def.Nsite * X->Bind.Def.Nsite, 2);
    X->Bind.Phys.spin_real_cor = d_1d_allocate(X->Bind.Def.Nsite * X->Bind.Def.Nsite);
    X->Bind.Phys.charge_real_cor = d_1d_allocate(X->Bind.Def.Nsite * X->Bind.Def.Nsite);
    X->Bind.Phys.loc_spin_z = d_1d_allocate(X->Bind.Def.Nsite * X->Bind.Def.Nsite);
    X->Bind.Def.Transfer = i_2d_allocate(X->Bind.Def.NTransfer, 4);
    X->Bind.Def.ParaTransfer = cd_1d_allocate(X->Bind.Def.NTransfer);
    X->Bind.Def.CoulombIntra = i_2d_allocate(X->Bind.Def.NCoulombIntra, 1);
    X->Bind.Def.ParaCoulombIntra = d_1d_allocate(X->Bind.Def.NCoulombIntra);
    X->Bind.Def.CoulombInter = i_2d_allocate(X->Bind.Def.NCoulombInter, 2);
    X->Bind.Def.ParaCoulombInter = d_1d_allocate(X->Bind.Def.NCoulombInter);
    X->Bind.Def.HundCoupling = i_2d_allocate(X->Bind.Def.NHundCoupling, 2);
    X->Bind.Def.ParaHundCoupling = d_1d_allocate(X->Bind.Def.NHundCoupling);
    X->Bind.Def.PairHopping = i_2d_allocate(X->Bind.Def.NPairHopping, 2);
    X->Bind.Def.ParaPairHopping = d_1d_allocate(X->Bind.Def.NPairHopping);
    X->Bind.Def.ExchangeCoupling = i_2d_allocate(X->Bind.Def.NExchangeCoupling, 2);
    X->Bind.Def.ParaExchangeCoupling = d_1d_allocate(X->Bind.Def.NExchangeCoupling);
    X->Bind.Def.Initial = i_2d_allocate(X->Bind.Def.NInitial, 4);
    X->Bind.Def.ParaInitial = cd_1d_allocate(X->Bind.Def.NInitial);
    X->Bind.Def.OrbitalIdx = i_2d_allocate(X->Bind.Def.Nsite * 2, X->Bind.Def.Nsite * 2);
    X->Bind.Def.OrbitalSgn = i_2d_allocate(X->Bind.Def.Nsite * 2, X->Bind.Def.Nsite * 2);
    X->Bind.Def.CisAjs = i_2d_allocate(X->Bind.Def.NCisAjs, 4);
    X->Bind.Def.CisAjsCktAltDC = i_2d_allocate(X->Bind.Def.NCisAjsCktAltDC, 8);
}