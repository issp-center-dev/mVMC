/* 2009.07.08*/
	//char_malloc2(X.Bind.Large.off_u, X.Bind.Def.EDNTransfer, X.Bind.Check.idim_max/2+1);
	d_malloc3(X.Bind.Large.Ham, 2, X.Bind.Def.Nsite, X.Bind.Def.Nsite);
	d_malloc3(X.Bind.Large.G,   2, X.Bind.Def.Nsite, X.Bind.Def.Nsite);
	d_malloc3(X.Bind.Large.G_old,   2, X.Bind.Def.Nsite, X.Bind.Def.Nsite);
	d_malloc2(X.Bind.Large.R_SLT_0,X.Bind.Def.Nsite, X.Bind.Def.Ne);
	d_malloc2(X.Bind.Large.L_SLT_0,X.Bind.Def.Ne, X.Bind.Def.Nsite);
	d_malloc2(X.Bind.Large.R_SLT_1,X.Bind.Def.Nsite, X.Bind.Def.Ne);
	d_malloc2(X.Bind.Large.L_SLT_1,X.Bind.Def.Ne, X.Bind.Def.Nsite);
	d_malloc2(X.Bind.Large.EigenValues,  2, X.Bind.Def.Nsite);
	d_malloc1(X.Bind.Large.tmp,  X.Bind.Def.Nsite*2);
  printf("LARGE ALLOCATE FINISH !\n");
 
