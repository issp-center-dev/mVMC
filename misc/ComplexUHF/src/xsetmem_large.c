/* 2009.07.08*/
	//char_malloc2(X.Bind.Large.off_u, X.Bind.Def.EDNTransfer, X.Bind.Check.idim_max/2+1);
	c_malloc2(X.Bind.Large.Ham, 2*X.Bind.Def.Nsite, 2*X.Bind.Def.Nsite);
	c_malloc2(X.Bind.Large.G,   2*X.Bind.Def.Nsite, 2*X.Bind.Def.Nsite);
	c_malloc2(X.Bind.Large.G_old, 2*X.Bind.Def.Nsite, 2*X.Bind.Def.Nsite);
	c_malloc2(X.Bind.Large.R_SLT,2*X.Bind.Def.Nsite, 2*X.Bind.Def.Ne);
	c_malloc2(X.Bind.Large.L_SLT,2*X.Bind.Def.Ne, 2*X.Bind.Def.Nsite);
	d_malloc1(X.Bind.Large.EigenValues,  2*X.Bind.Def.Nsite);
  printf("LARGE ALLOCATE FINISH !\n");
 
