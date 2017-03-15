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
/* 2009.07.08*/

	ui_malloc1(X.Bind.Def.Tpow, 2*X.Bind.Def.Nsite+1);
	i_malloc1(X.Bind.Def.LocSpn, X.Bind.Def.Nsite);
	i_malloc2(X.Bind.Def.All_pair, X.Bind.Def.Nsite*X.Bind.Def.Nsite, 2);
	d_malloc1(X.Bind.Phys.spin_real_cor, X.Bind.Def.Nsite*X.Bind.Def.Nsite);
d_malloc1(X.Bind.Phys.charge_real_cor, X.Bind.Def.Nsite*X.Bind.Def.Nsite);
d_malloc1(X.Bind.Phys.loc_spin_z, X.Bind.Def.Nsite*X.Bind.Def.Nsite);
  
	i_malloc2(X.Bind.Def.Transfer, X.Bind.Def.NTransfer, 4);
	c_malloc1(X.Bind.Def.ParaTransfer, X.Bind.Def.NTransfer);
	i_malloc2(X.Bind.Def.CoulombIntra, X.Bind.Def.NCoulombIntra, 1);
	d_malloc1(X.Bind.Def.ParaCoulombIntra, X.Bind.Def.NCoulombIntra);
	i_malloc2(X.Bind.Def.CoulombInter, X.Bind.Def.NCoulombInter, 2);
	d_malloc1(X.Bind.Def.ParaCoulombInter, X.Bind.Def.NCoulombInter);
	i_malloc2(X.Bind.Def.HundCoupling, X.Bind.Def.NHundCoupling, 2);
	d_malloc1(X.Bind.Def.ParaHundCoupling, X.Bind.Def.NHundCoupling);
	i_malloc2(X.Bind.Def.PairHopping, X.Bind.Def.NPairHopping, 2);
	d_malloc1(X.Bind.Def.ParaPairHopping, X.Bind.Def.NPairHopping);
	i_malloc2(X.Bind.Def.ExchangeCoupling, X.Bind.Def.NExchangeCoupling, 2);
	d_malloc1(X.Bind.Def.ParaExchangeCoupling, X.Bind.Def.NExchangeCoupling);
	i_malloc2(X.Bind.Def.Initial, X.Bind.Def.NInitial, 4);
c_malloc1(X.Bind.Def.ParaInitial, X.Bind.Def.NInitial);
//d_malloc1(X.Bind.Def.ParaInitial, X.Bind.Def.NInitial);
//d_malloc1(X.Bind.Def.ParaInitial_theta, X.Bind.Def.NInitial);
i_malloc2(X.Bind.Def.OrbitalIdx, X.Bind.Def.Nsite*2, X.Bind.Def.Nsite*2);

	i_malloc2(X.Bind.Def.CisAjs, X.Bind.Def.NCisAjs, 4);
  	i_malloc2(X.Bind.Def.CisAjsCktAltDC, X.Bind.Def.NCisAjsCktAltDC, 8);
