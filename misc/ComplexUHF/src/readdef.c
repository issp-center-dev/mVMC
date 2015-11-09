/*-------------------------------------------------------------
 *[ver.2008.11.4]
 *  Read Definition files
 *-------------------------------------------------------------
 * Copyright (C) 2007-2009 Daisuke Tahara. All rights reserved.
 *-------------------------------------------------------------*/


/*=================================================================================================*/
#define D_FileNameMaxReadDef 200
#define D_CharTmpReadDef     200

int ReadDefFileError(
	char *defname
){
	printf("Error: %s (Broken file or Not exist)\n", defname);
	return 0;
}


int ReadDefFileNInt(
	char *xNameListFile, 
	struct DefineList *X
){
	FILE *fp, *fplist;
	char defname[D_FileNameMaxReadDef];
	char ctmp[D_CharTmpReadDef];
	int itmp;

	fplist = fopen(xNameListFile, "r");
	if(fplist==NULL) return ReadDefFileError(xNameListFile);
	/*=======================================================================*/
	/*modpara.def---------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &itmp); //2
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //3
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //4
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //5
	fscanf(fp,"%s %s\n", ctmp, X->CDataFileHead); //6
	fscanf(fp,"%s %s\n", ctmp, X->CParaFileHead); //7
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);   //8
	fscanf(fp,"%s %d\n", ctmp, &(X->Nsite));      //9
	fscanf(fp,"%s %d\n", ctmp, &(X->Ne));         //10
	fscanf(fp,"%s %lf\n", ctmp, &(X->mix));//11
	fscanf(fp,"%s %d\n", ctmp, &(X->eps_int));//12
	fscanf(fp,"%s %d\n", ctmp, &(X->print));//13
	fscanf(fp,"%s %d\n", ctmp, &(X->IterationMax));//14
	fclose(fp);
	/*locspn.def----------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &(X->NLocSpn));
	fclose(fp);
	/*transfer.def--------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &(X->NTransfer));
	fclose(fp);
	/*coulombintra.def----------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &(X->NCoulombIntra));
	fclose(fp);
	/*coulombinter.def----------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &(X->NCoulombInter));
	fclose(fp);
	/*hund.def------------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &(X->NHundCoupling));
	fclose(fp);
	/*pairhop.def---------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &(X->NPairHopping));
	fclose(fp);
	/*exchange.def--------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &(X->NExchangeCoupling));
	fclose(fp);
  /*cisajs.def----------------------------------------*/
  if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
  fp = fopen(defname, "r");
  if(fp==NULL) return ReadDefFileError(defname);
  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
  fscanf(fp,"%s %d\n", ctmp, &(X->NCisAjs));
  fclose(fp);
  /*cisajscktaltdc.def--------------------------------*/
  if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
  fp = fopen(defname, "r");
  if(fp==NULL) return ReadDefFileError(defname);
  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
  fscanf(fp,"%s %d\n", ctmp, &(X->NCisAjsCktAltDC));
  fclose(fp);
	/*initial.def----------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	fscanf(fp,"%s %d\n", ctmp, &(X->NInitial));
	fclose(fp);
	/*=======================================================================*/
	fclose(fplist);

	X->Nsize   = 2*X->Ne;
	X->fidx = 0;

	return 1;
}


int ReadDefFileIdxPara(
	char *xNameListFile, 
	struct DefineList *X
){
	FILE *fp, *fplist;
	char defname[D_FileNameMaxReadDef];
	char ctmp[D_CharTmpReadDef];
	int itmp;

	int i,j,n,idx;
	int xitmp[8];

	fplist = fopen(xNameListFile, "r");
	if(fplist==NULL) return ReadDefFileError(xNameListFile);
	/*=======================================================================*/
	/*modpara.def---------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	fclose(fp);
	/*locspn.def----------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	fp = fopen(defname, "r");
	if(fp==NULL) return ReadDefFileError(defname);
	for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	idx = 0;
	while( fscanf(fp, "%d %d\n", &(xitmp[0]), &(xitmp[1]) )!=EOF){
		X->LocSpn[xitmp[0]] = xitmp[1];
		idx++;
	}
	if(X->NLocSpn>2*X->Ne){
		printf("Error: 2*Ne must be (2*Ne >= NLocalSpin).\n");
		return 0;
	}
	//if(X->NLocSpn>0 && X->NExUpdatePath==0){
	//	printf("Error: NExUpdatePath (in modpara.def) must be 1.\n");
	//	return 0;
	//}
	if(idx!=X->Nsite) return ReadDefFileError(defname);
	fclose(fp);
	/*transfer.def--------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	if(X->NTransfer>0){
		fp = fopen(defname, "r");
		if(fp==NULL) return ReadDefFileError(defname);
		for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
		idx = 0;
		while( fscanf(fp, "%d %d %d %lf\n", 
					&(X->Transfer[idx][0]),
					&(X->Transfer[idx][1]),
					&(X->Transfer[idx][2]),
					&(X->ParaTransfer[idx])
				)!=EOF
		){
			idx++;
		}
		if(idx!=X->NTransfer) return ReadDefFileError(defname);
		fclose(fp);
	}
	/*coulombintra.def----------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	if(X->NCoulombIntra>0){
		fp = fopen(defname, "r");
		if(fp==NULL) return ReadDefFileError(defname);
		for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
		idx = 0;
		while( fscanf(fp, "%d %lf\n", 
					&(X->CoulombIntra[idx][0]),
					&(X->ParaCoulombIntra[idx])
				)!=EOF
		){
			idx++;
		}
		if(idx!=X->NCoulombIntra) return ReadDefFileError(defname);
		fclose(fp);
	}
	/*coulombinter.def----------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	if(X->NCoulombInter>0){
		fp = fopen(defname, "r");
		if(fp==NULL) return ReadDefFileError(defname);
		for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
		idx = 0;
		while( fscanf(fp, "%d %d %lf\n", 
					&(X->CoulombInter[idx][0]),
					&(X->CoulombInter[idx][1]),
					&(X->ParaCoulombInter[idx])
				)!=EOF
		){
			idx++;
		}
		if(idx!=X->NCoulombInter) return ReadDefFileError(defname);
		fclose(fp);
	}
	/*hund.def------------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	if(X->NHundCoupling>0){
		fp = fopen(defname, "r");
		if(fp==NULL) return ReadDefFileError(defname);
		for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
		idx = 0;
		while( fscanf(fp, "%d %d %lf\n", 
					&(X->HundCoupling[idx][0]),
					&(X->HundCoupling[idx][1]),
					&(X->ParaHundCoupling[idx])
				)!=EOF
		){
			idx++;
		}
		if(idx!=X->NHundCoupling) return ReadDefFileError(defname);
		fclose(fp);
	}
	/*pairhop.def---------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	if(X->NPairHopping>0){
		fp = fopen(defname, "r");
		if(fp==NULL) return ReadDefFileError(defname);
		for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
		idx = 0;
		while( fscanf(fp, "%d %d %lf\n", 
					&(X->PairHopping[idx][0]),
					&(X->PairHopping[idx][1]),
					&(X->ParaPairHopping[idx])
				)!=EOF
		){
			idx++;
		}
		if(idx!=X->NPairHopping) return ReadDefFileError(defname);
		fclose(fp);
	}
	/*exchange.def--------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	if(X->NExchangeCoupling>0){
		fp = fopen(defname, "r");
		if(fp==NULL) return ReadDefFileError(defname);
		for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
		idx = 0;
		while( fscanf(fp, "%d %d %lf\n", 
					&(X->ExchangeCoupling[idx][0]),
					&(X->ExchangeCoupling[idx][1]),
					&(X->ParaExchangeCoupling[idx])
				)!=EOF
		){
			idx++;
		}
		if(idx!=X->NExchangeCoupling) return ReadDefFileError(defname);
		fclose(fp);
	}
 /*cisajs.def----------------------------------------*/
 if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
 if(X->NCisAjs>0){
   fp = fopen(defname, "r");
   if(fp==NULL) return ReadDefFileError(defname);
   for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
   idx = 0;
   while( fscanf(fp, "%d %d %d %d\n", &(xitmp[0]), &(xitmp[1]), &(xitmp[2]), &(xitmp[3])) != EOF){
     X->CisAjs[ xitmp[0] ][0] = xitmp[1];
     X->CisAjs[ xitmp[0] ][1] = xitmp[2];
     X->CisAjs[ xitmp[0] ][2] = xitmp[3];
     idx++;
   }
   if(idx!=X->NCisAjs) return ReadDefFileError(defname);
   fclose(fp);
 }
 /*cisajscktaltdc.def--------------------------------*/
 if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
 if(X->NCisAjsCktAltDC>0){
   fp = fopen(defname, "r");
   if(fp==NULL) return ReadDefFileError(defname);
   for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
   idx = 0;
   while( fscanf(fp, "%d %d %d %d %d %d\n",
     &(xitmp[0]), &(xitmp[1]), &(xitmp[2]), &(xitmp[3]), &(xitmp[4]), &(xitmp[5])
   ) != EOF
   ){
     X->CisAjsCktAltDC[idx][0] = xitmp[0];
     X->CisAjsCktAltDC[idx][1] = xitmp[1];
     X->CisAjsCktAltDC[idx][2] = xitmp[2];
     X->CisAjsCktAltDC[idx][3] = xitmp[3];
     X->CisAjsCktAltDC[idx][4] = xitmp[4];
     X->CisAjsCktAltDC[idx][5] = xitmp[5];
     idx++;
   }
   if(idx!=X->NCisAjsCktAltDC) return ReadDefFileError(defname);
   fclose(fp);
 }
	/*intial.def--------------------------------------*/
	if(fscanf(fplist, "%s\n", defname)==EOF) return ReadDefFileError(xNameListFile);
	if(X->NTransfer>0){
		fp = fopen(defname, "r");
		if(fp==NULL) return ReadDefFileError(defname);
		for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
		idx = 0;
		while( fscanf(fp, "%d %d %d %d %lf %lf\n", 
					&(X->Initial[idx][0]),
					&(X->Initial[idx][1]),
					&(X->Initial[idx][2]),
					&(X->Initial[idx][3]),
					&(X->ParaInitial[idx]),
					&(X->ParaInitial_theta[idx])
				)!=EOF
		){
			idx++;
		}
		if(idx!=X->NInitial) return ReadDefFileError(defname);
		fclose(fp);
	}
 /*=======================================================================*/
	fclose(fplist);
	return 1;
}
