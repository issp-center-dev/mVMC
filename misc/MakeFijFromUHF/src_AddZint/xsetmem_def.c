/* 2009.07.08*/
	ui_malloc1(X.Bind.Def.Tpow, 2*X.Bind.Def.Nsite+1);
	i_malloc1(X.Bind.Def.LocSpn, X.Bind.Def.Nsite);
  i_malloc2(X.Bind.Def.All_pair, X.Bind.Def.Nsite*X.Bind.Def.Nsite, 2);
  d_malloc1(X.Bind.Phys.spin_real_cor, X.Bind.Def.Nsite*X.Bind.Def.Nsite);
  d_malloc1(X.Bind.Phys.charge_real_cor, X.Bind.Def.Nsite*X.Bind.Def.Nsite);
  d_malloc1(X.Bind.Phys.loc_spin_z, X.Bind.Def.Nsite*X.Bind.Def.Nsite);
  
	i_malloc1(X.Bind.Def.EDChemi, X.Bind.Def.Nsite);
	d_malloc1(X.Bind.Def.EDParaChemi, X.Bind.Def.Nsite);
	i_malloc2(X.Bind.Def.Transfer, X.Bind.Def.NTransfer, 3);
	i_malloc2(X.Bind.Def.EDTransfer, X.Bind.Def.NTransfer, 2);
	d_malloc1(X.Bind.Def.ParaTransfer, X.Bind.Def.NTransfer);
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
	i_malloc2(X.Bind.Def.InterAll, X.Bind.Def.NInterAll, 6);
	d_malloc1(X.Bind.Def.ParaInterAll, X.Bind.Def.NInterAll);
	i_malloc2(X.Bind.Def.Initial, X.Bind.Def.NInitial, 3);
	d_malloc1(X.Bind.Def.ParaInitial, X.Bind.Def.NInitial);

  i_malloc2(X.Bind.Def.CisAjs, X.Bind.Def.NCisAjs, 3);
  i_malloc2(X.Bind.Def.CisAjsCktAltDC, X.Bind.Def.NCisAjsCktAltDC, 6);
