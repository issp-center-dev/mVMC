/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Read Definition Files
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#include <ctype.h>
#include "./include/readdef.h"

int ReadDefFileError(const char *defname);
int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm);
int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm);

int ReadDefFileError(const char *defname){
  fprintf(stderr, "error: %s (Broken file or Not exist)\n", defname);
  return 1;
}

int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm){
  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  char ctmp2[D_FileNameMax];

  int i=0;
  int itmp,tmp_info;

  int rank, info=0;
  const int nBufInt= ParamIdxInt_End;
  const int nBufDouble= ParamIdxDouble_End;
  const int nBufChar=D_FileNameMax;
  int bufInt[nBufInt];
  double bufDouble[nBufDouble];
  int iKWidx=0;

  MPI_Comm_rank(comm, &rank);

  if(rank==0) {
    cFileNameListFile = malloc(sizeof(char)*D_CharTmpReadDef*KWIdxInt_end);
    fprintf(stdout, "  Read File %s .\n", xNameListFile); 
    if(GetFileName(xNameListFile, cFileNameListFile)!=0){
      fprintf(stderr, "  error: Definition files(*.def) are incomplete.\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    
    for(iKWidx=0; iKWidx< KWIdxInt_end; iKWidx++){ 
      strcpy(defname, cFileNameListFile[iKWidx]);
      if(strcmp(defname,"")==0){
	switch (iKWidx){
	case KWModPara:
	case KWLocSpin:
	case KWOrbital:
	  fprintf(stderr, "  Error: Need to make a def file for %s.\n", cKWListOfFileNameList[iKWidx]);
	  MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	  break;
	default:
	  break;
	}
      }
    } 

    for(iKWidx=0; iKWidx< KWIdxInt_end; iKWidx++){ 
      strcpy(defname, cFileNameListFile[iKWidx]);
      if(strcmp(defname,"")==0) continue;
      fprintf(stdout,  "  Read File '%s' for %s.\n", defname, cKWListOfFileNameList[iKWidx]);
      fp = fopen(defname, "r");
      if(fp==NULL){
	info=ReadDefFileError(defname);
	fclose(fp);
	break;
      }
      else{

	switch(iKWidx){
	case KWModPara:
	  /* Read modpara.def---------------------------------------*/
	  //TODO: add error procedure here when parameters are not enough.
	  SetDefultValuesModPara(bufInt, bufDouble);
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &itmp); //2
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //3
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //4
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //5
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %s\n", ctmp, CDataFileHead); //6
	  fgets(ctmp2,sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %s\n", ctmp, CParaFileHead); //7
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);   //8

	  double dtmp;
	  while(fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp)!=NULL){
	    if(*ctmp2 == '\n' || ctmp2[0] == '-')  continue;
	    sscanf(ctmp2,"%s %lf\n", ctmp, &dtmp);
	    if(CheckWords(ctmp, "NVMCCalMode")==0){
	      bufInt[IdxVMCCalcMode]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NLanczosMode")==0){
	      bufInt[IdxLanczosMode]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NDataIdxStart")==0){
	      bufInt[IdxDataIdxStart]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NDataQtySmp")==0){
	      bufInt[IdxDataQtySmp]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "Nsite")==0){
	      bufInt[IdxNsite]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "Ne")==0 || CheckWords(ctmp, "Nelectron")==0 ){
	      bufInt[IdxNe]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NSPGaussLeg")==0){
	      bufInt[IdxSPGaussLeg]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NSPStot")==0){
	      bufInt[IdxSPStot]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NMPTrans")==0){
	      bufInt[IdxMPTrans]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NSROptItrStep")==0){
	      bufInt[IdxSROptItrStep]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NSROptItrSmp")==0){
	      bufInt[IdxSROptItrSmp]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NSROptFixSmp")==0){
	      bufInt[IdxSROptFixSmp]=(int)dtmp;
	    }	
	    else if(CheckWords(ctmp, "DSROptRedCut")==0){
	      bufDouble[IdxSROptRedCut]=(double)dtmp;
	    }	
	    else if(CheckWords(ctmp, "DSROptStaDel")==0){
	      bufDouble[IdxSROptStaDel]=(double)dtmp;
	    }	
	    else if(CheckWords(ctmp, "DSROptStepDt")==0){
	      bufDouble[IdxSROptStepDt]=(double)dtmp;
	    }	
	    else if(CheckWords(ctmp, "NVMCWarmUp")==0){
	      bufInt[IdxVMCWarmUp]=(int)dtmp;
	    }	
	    else if(CheckWords(ctmp, "NVMCIniterval")==0){
	      bufInt[IdxVMCIniterval]=(int)dtmp;
	    }	
	    else if(CheckWords(ctmp, "NVMCSample")==0){
	      bufInt[IdxVMCSample]=(int)dtmp;
	    }	
	    else if(CheckWords(ctmp, "NExUpdatePath")==0){
	      bufInt[IdxExUpdatePath]=(int)dtmp;
	    }	
	    else if(CheckWords(ctmp, "RndSeed")==0){
	      bufInt[IdxRndSeed] = (int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NSplitSize")==0){
	      bufInt[IdxSplitSize]=(int)dtmp;
	    }
	    else if(CheckWords(ctmp, "NStore")==0){
	      NStoreO=(int)dtmp;
	    }
	    else{
	      fprintf(stderr, "  Error: keyword \" %s \" is incorrect. \n", ctmp);
	      info = ReadDefFileError(defname);
	      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);			 
	    }
	  }
	  if(bufInt[IdxRndSeed]<0) {
	    bufInt[IdxRndSeed] = (int)time(NULL);
	    fprintf(stdout, "  remark: Seed = %d\n", bufInt[IdxRndSeed]);
	  }
	  fclose(fp);
	  break;//modpara file

	case KWLocSpin:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNLocSpin]);
	  fclose(fp);
	  break;

	case KWTrans:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNTrans]);
	  fclose(fp);
	  break;

	case KWCoulombIntra:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNCoulombIntra]);
	  fclose(fp);
	  break;

	case KWCoulombInter:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNCoulombInter]);
	  fclose(fp);
	  break;

	case KWHund:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNHund]);
	  fclose(fp);
	  break;
	  
	case KWPairHop:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNPairHop]);
	  fclose(fp);
	  break;
		  
	case KWExchange:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNExchange]);
	  fclose(fp);
	  break;

	case KWGutzwiller:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNGutz]);
      fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
      sscanf(ctmp2,"%s %d\n", ctmp, &iComplexFlgGutzwiller);
	  fclose(fp);
	  break;

	case KWJastrow:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNJast]);
      fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
      sscanf(ctmp2,"%s %d\n", ctmp, &iComplexFlgJastrow);
	  fclose(fp);
	  break;

	case KWDH2:
      fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
      sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNDH2]);
      fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
      sscanf(ctmp2,"%s %d\n", ctmp, &iComplexFlgDH2);
	  fclose(fp);
	  break;

	case KWDH4:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNDH4]);
      fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
      sscanf(ctmp2,"%s %d\n", ctmp, &iComplexFlgDH4);
	  fclose(fp);
	  break;

	case KWOrbital:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNOrbit]);
      fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
      sscanf(ctmp2,"%s %d\n", ctmp, &iComplexFlgOrbital);
	  fclose(fp);
	  break;

	case KWTransSym:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNQPTrans]);
	  fclose(fp);
	  break;
      
	case KWOneBodyG:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNOneBodyG]);
	  fclose(fp);
	  break;

	case KWTwoBodyG:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNTwoBodyG]);
	  fclose(fp);
	  break;

	case KWTwoBodyGEx:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNTwoBodyGEx]);
	  fclose(fp);
	  break;

	case KWInterAll:
	  fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	  fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	  sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNInterAll]);
	  fclose(fp);
	  break;		

	case KWOptTrans:
	  bufInt[IdxNQPOptTrans]=1;
	  if(FlagOptTrans>0) {
	    fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
	    fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
	    sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNQPOptTrans]);
	    if(bufInt[IdxNQPOptTrans]<1) {
	      fprintf(stderr, "error: NQPOptTrans should be larger than 0.\n");
	      info = ReadDefFileError(defname);
	    }
	  }
	  fclose(fp);
	  break;
	  
	default:
	  fclose(fp);
	  break;
	}//case KW
      }
    }
  }
  
  if(info!=0) {
    if(rank==0) {
      fprintf(stderr, "error: Definition files(*.def) are incomplete.\n");
    }
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

#ifdef _mpi_use
  MPI_Bcast(bufInt, nBufInt, MPI_INT, 0, comm);
  MPI_Bcast(&NStoreO, 1, MPI_INT, 0, comm); // for NStoreO
  MPI_Bcast(bufDouble, nBufDouble, MPI_DOUBLE, 0, comm);
  MPI_Bcast(CDataFileHead, nBufChar, MPI_CHAR, 0, comm);
  MPI_Bcast(CParaFileHead, nBufChar, MPI_CHAR, 0, comm);
#endif /* _mpi_use */

  NVMCCalMode            =  bufInt[IdxVMCCalcMode];
  NLanczosMode           =  bufInt[IdxLanczosMode];
  NDataIdxStart          =  bufInt[IdxDataIdxStart];
  NDataQtySmp            =  bufInt[IdxDataQtySmp];
  Nsite                  =  bufInt[IdxNsite];
  Ne                     =  bufInt[IdxNe];
  NSPGaussLeg            =  bufInt[IdxSPGaussLeg];
  NSPStot                =  bufInt[IdxSPStot];
  NMPTrans               =  bufInt[IdxMPTrans];
  NSROptItrStep          =  bufInt[IdxSROptItrStep];
  NSROptItrSmp           =  bufInt[IdxSROptItrSmp];
  NSROptFixSmp           =  bufInt[IdxSROptFixSmp];
  NVMCWarmUp             =  bufInt[IdxVMCWarmUp];
  NVMCIniterval          =  bufInt[IdxVMCIniterval];
  NVMCSample             =  bufInt[IdxVMCSample];
  NExUpdatePath          =  bufInt[IdxExUpdatePath];
  RndSeed                =  bufInt[IdxRndSeed];
  NSplitSize             =  bufInt[IdxSplitSize];
  NLocSpn                =  bufInt[IdxNLocSpin];
  NTransfer              =  bufInt[IdxNTrans];
  NCoulombIntra          =  bufInt[IdxNCoulombIntra];
  NCoulombInter          =  bufInt[IdxNCoulombInter];
  NHundCoupling          =  bufInt[IdxNHund];
  NPairHopping           =  bufInt[IdxNPairHop];
  NExchangeCoupling      =  bufInt[IdxNExchange];
  NGutzwillerIdx         =  bufInt[IdxNGutz];
  NJastrowIdx            =  bufInt[IdxNJast];
  NDoublonHolon2siteIdx  =  bufInt[IdxNDH2];
  NDoublonHolon4siteIdx  =  bufInt[IdxNDH4];
  NOrbitalIdx            =  bufInt[IdxNOrbit];
  NQPTrans               =  bufInt[IdxNQPTrans];
  NCisAjs                =  bufInt[IdxNOneBodyG];
  NCisAjsCktAlt          =  bufInt[IdxNTwoBodyG];
  NCisAjsCktAltDC        =  bufInt[IdxNTwoBodyGEx];
  NInterAll              =  bufInt[IdxNInterAll];
  NQPOptTrans            =  bufInt[IdxNQPOptTrans];

  DSROptRedCut = bufDouble[IdxSROptRedCut];
  DSROptStaDel = bufDouble[IdxSROptStaDel];
  DSROptStepDt = bufDouble[IdxSROptStepDt];

  if(NMPTrans < 0) {
    APFlag = 1; /* anti-periodic boundary */
    NMPTrans *= -1;
  } else {
    APFlag = 0;
  }

  if(DSROptStepDt < 0) {
    SRFlag = 1; /* diagonalization */
    if(rank==0) fprintf(stderr, "remark: Diagonalization Mode\n");
    DSROptStepDt *= -1;
  } else {
    SRFlag = 0;
  }
  Nsize   = 2*Ne;
  Nsite2  = 2*Nsite;
  NSlater = NOrbitalIdx;
  NProj   = NGutzwillerIdx + NJastrowIdx
    + 2*3*NDoublonHolon2siteIdx
    + 2*5*NDoublonHolon4siteIdx;
  NOptTrans = (FlagOptTrans>0) ? NQPOptTrans : 0;
  NPara   = NProj + NSlater + NOptTrans ; 
  NQPFix = NSPGaussLeg * NMPTrans;
  NQPFull = NQPFix * NQPOptTrans;
  SROptSize = NPara+1;
  
  NTotalDefInt = Nsite /* LocSpn */
    + 4*NTransfer /* Transfer */
    + NCoulombIntra /* CoulombIntra */
    + 2*NCoulombInter /* CoulombInter */
    + 2*NHundCoupling /* HundCoupling */
    + 2*NPairHopping /* PairHopping */
    + 2*NExchangeCoupling /* ExchangeCoupling */
    + Nsite /* GutzwillerIdx */
    + Nsite*Nsite /* JastrowIdx */
    + 2*Nsite*NDoublonHolon2siteIdx /* DoublonHolon2siteIdx */
    + 4*Nsite*NDoublonHolon4siteIdx /* DoublonHolon4siteIdx */
    + Nsite*Nsite /* OrbitalIdx */
    + Nsite*Nsite /* OrbitalSgn */
    + Nsite*NQPTrans /* QPTrans */
    + Nsite*NQPTrans /* QPTransSgn */
    + 4*NCisAjs /* CisAjs */
    + 8*NCisAjsCktAlt /* CisAjsCktAlt */
    + 8*NCisAjsCktAltDC /* CisAjsCktAltDC */
    + 8*NInterAll /* InterAll */
    + Nsite*NQPOptTrans /* QPOptTrans */
    + Nsite*NQPOptTrans /* QPOptTransSgn */
    + 2*NPara; /* OptFlag */ // TBC
  
  NTotalDefDouble =
    NCoulombIntra /* ParaCoulombIntra */
    + NCoulombInter /* ParaCoulombInter */
    + NHundCoupling /* ParaHondCoupling */
    + NPairHopping  /* ParaPairHopping */
    + NExchangeCoupling /* ParaExchangeCoupling */
    + NQPTrans /* ParaQPTrans */
    //+ NInterAll /* ParaInterAll */
    + NQPOptTrans; /* ParaQPTransOpt */

  return 0;
}

int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm){
  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax], ctmp2[256];
  int itmp;
  int iKWidx=0;
  
  int i,j,n,idx,idx0,idx1,info=0;
  int fidx=0; /* index for OptFlag */
    int count_idx=0;
  int x0,x1,x2,x3,x4,x5,x6,x7;
  double dReValue, dImValue;
  int tmp_ispin;
  int rank;
  double tmp_real,tmp_comp;

  MPI_Comm_rank(comm, &rank);
  
  if(rank==0) {
    for(iKWidx=KWLocSpin; iKWidx< KWIdxInt_end; iKWidx++){     
      strcpy(defname, cFileNameListFile[iKWidx]);

      if(strcmp(defname,"")==0) continue;   

      fp = fopen(defname, "r");
      if(fp==NULL){
        info= ReadDefFileError(defname);
        fclose(fp);
        continue;
      }

      /*=======================================================================*/
      for(i=0;i<IgnoreLinesInDef;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
      idx=0;

      switch(iKWidx){
      case KWLocSpin:
	/* Read locspn.def----------------------------------------*/
	while( fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp) != NULL){
	  sscanf(ctmp2, "%d %d\n", &(x0), &(x1) );
	  LocSpn[x0] = x1;
	  idx++;
	}
	if(NLocSpn>2*Ne){
	  fprintf(stderr, "Error: 2*Ne must be (2*Ne >= NLocalSpin).\n");
	  info=1;
	}
	if(NLocSpn>0 && NExUpdatePath==0){
	  fprintf(stderr, "Error: NExUpdatePath (in modpara.def) must be 1.\n");
	  info=1;
	}
	if(idx!=Nsite) info = ReadDefFileError(defname);
	fclose(fp);	  
	break;//locspn
	      
      case KWTrans:
	/* transfer.def--------------------------------------*/
	if(NTransfer>0){
	  while( fscanf(fp, "%d %d %d %d %lf %lf\n",
			&(Transfer[idx][0]),
			&(Transfer[idx][1]),
			&(Transfer[idx][2]),
			&(Transfer[idx][3]),
			&dReValue,
	  		&dImValue)!=EOF){

		  	ParaTransfer[idx]=dReValue+I*dImValue;
		  	if(Transfer[idx][1] != Transfer[idx][3]){
				fprintf(stderr, "  Error:  Sz non-conserved system is not yet supported in mVMC ver.1.0.\n");
				info = ReadDefFileError(defname);
				break;
			}
	    idx++;
	  }
	  if(idx!=NTransfer) info = ReadDefFileError(defname);
	}
	fclose(fp);
	break;

      case KWCoulombIntra:
	/*coulombintra.def----------------------------------*/
	if(NCoulombIntra>0){
	  while( fscanf(fp, "%d %lf\n", 
			&(CoulombIntra[idx]),
			&(ParaCoulombIntra[idx]) )!=EOF){
	    idx++;
	  }
	  if(idx!=NCoulombIntra) info = ReadDefFileError(defname);
	}
	fclose(fp);
	break;

      case KWCoulombInter:
	/*coulombinter.def----------------------------------*/
	if(NCoulombInter>0){
	  while(fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp) != NULL){
	    sscanf(ctmp2, "%d %d %lf\n",
		   &(CoulombInter[idx][0]),
		   &(CoulombInter[idx][1]),
		   &(ParaCoulombInter[idx])
		   );
	    idx++;
	  }
	  if(idx!=NCoulombInter) info=ReadDefFileError(defname);
	}
	fclose(fp);
	break;

      case KWHund:
	/*hund.def------------------------------------------*/
	if(NHundCoupling>0){
	  while( fscanf(fp, "%d %d %lf\n", 
			&(HundCoupling[idx][0]),
			&(HundCoupling[idx][1]),
			&(ParaHundCoupling[idx]) )!=EOF){
	    idx++;
	  }
	  if(idx!=NHundCoupling) info=ReadDefFileError(defname);
	}
	fclose(fp);
	break;

      case KWPairHop:
	/*pairhop.def---------------------------------------*/
	if(NPairHopping>0){
	  while( fscanf(fp, "%d %d %lf\n", 
			&(PairHopping[idx][0]),
			&(PairHopping[idx][1]),
			&(ParaPairHopping[idx]) )!=EOF){
	    idx++;
	  }
	  if(idx!=NPairHopping) info=ReadDefFileError(defname);
	}
	fclose(fp);
	break;
	
      case KWExchange:
	/*exchange.def--------------------------------------*/
	if(NExchangeCoupling>0){
	  while( fscanf(fp, "%d %d %lf\n", 
			&(ExchangeCoupling[idx][0]),
			&(ExchangeCoupling[idx][1]),
			&(ParaExchangeCoupling[idx]) )!=EOF){
	    idx++;
	  }
	  if(idx!=NExchangeCoupling) info=ReadDefFileError(defname);
	}
	fclose(fp);
	break;

      case KWGutzwiller:
	/*gutzwilleridx.def---------------------------------*/
	if(NGutzwillerIdx>0){
	  idx0=idx1=0;
	  while( fscanf(fp, "%d ", &i) != EOF){
	    fscanf(fp, "%d\n", &(GutzwillerIdx[i]));
	    idx0++;
	    if(idx0==Nsite) break;
	  }
      fidx=0;
	  while( fscanf(fp, "%d ", &i) != EOF){
	    fscanf(fp, "%d\n", &(OptFlag[2*fidx])); // TBC real

	    OptFlag[2*fidx+1] = iComplexFlgGutzwiller; //  TBC imaginary
	    fidx++;
	    idx1++;
        count_idx++;
	  }
	  if(idx0!=Nsite || idx1!=NGutzwillerIdx) {
	    info=ReadDefFileError(defname);
	  }
	}
	fclose(fp);
	break;

      case KWJastrow:
	/*jastrowidx.def------------------------------------*/
	if(NJastrowIdx>0){
	  idx0 = idx1 = 0;
	  while( fscanf(fp, "%d %d ", &i, &j) != EOF){
	    if(i==j){
	      fprintf(stderr, "Error in %s: [Condition] i neq j\n", defname);
	      info=1;
	      break;
	    }
	    fscanf(fp, "%d\n", &(JastrowIdx[i][j]));
	    JastrowIdx[i][i] = -1; // This case is Gutzwiller.
	    idx0++;
	    if(idx0==Nsite*(Nsite-1)) break;
	  }
      fidx=NGutzwillerIdx;
	  while( fscanf(fp, "%d ", &i) != EOF){
	    fscanf(fp, "%d\n", &(OptFlag[2*fidx])); // TBC real
	    OptFlag[2*fidx+1] = iComplexFlgJastrow; //  TBC imaginary
	    fidx++;
	    idx1++;
          count_idx++;

      }
	  if(idx0!=Nsite*(Nsite-1) || idx1!=NJastrowIdx) {
	    info=ReadDefFileError(defname);
	  }
	}
	fclose(fp);
	break;
	  
      case KWDH2:
	/*doublonholon2siteidx.def--------------------------*/
	if(NDoublonHolon2siteIdx>0){
	  idx0 = idx1 = 0;
	  while( fscanf(fp, "%d %d %d %d\n", &i, &(x0), &(x1), &n) != EOF){
	    DoublonHolon2siteIdx[n][2*i]   = x0;
	    DoublonHolon2siteIdx[n][2*i+1] = x1;
	    idx0++;
	    if(idx0==Nsite*NDoublonHolon2siteIdx) break;
	  }
      fidx=NGutzwillerIdx+NJastrowIdx;
	  while( fscanf(fp, "%d ", &i) != EOF){
	    fscanf(fp, "%d\n", &(OptFlag[2*fidx]));//TBC real
	    OptFlag[2*fidx+1] = iComplexFlgDH2; //  TBC imaginary
	    fidx++;
	    idx1++;
        count_idx++;
      }
	  if(idx0!=Nsite*NDoublonHolon2siteIdx
	     || idx1!=2*3*NDoublonHolon2siteIdx) {
	    info=ReadDefFileError(defname);
	  }
	}
	fclose(fp);
	break;

      case KWDH4:
	/*doublonholon4siteidx.def--------------------------*/
	if(NDoublonHolon4siteIdx>0){
	  idx0 = idx1 = 0;
	  while( fscanf(fp, "%d %d %d %d %d %d\n",
			&i, &(x0), &(x1), &(x2), &(x3), &n) != EOF){
	    DoublonHolon4siteIdx[n][4*i]   = x0;
	    DoublonHolon4siteIdx[n][4*i+1] = x1;
	    DoublonHolon4siteIdx[n][4*i+2] = x2;
	    DoublonHolon4siteIdx[n][4*i+3] = x3;
	    idx0++;
	    if(idx0==Nsite*NDoublonHolon4siteIdx) break;
	  }
      fidx=NGutzwillerIdx+NJastrowIdx+2*3*NDoublonHolon2siteIdx;
        while( fscanf(fp, "%d ", &i) != EOF){
	    fscanf(fp, "%d\n", &(OptFlag[2*fidx]));
	    OptFlag[2*fidx+1] = 0; //  TBC imaginary
	    fidx++;
	    idx1++;
            count_idx++;
        }
	  if(idx0!=Nsite*NDoublonHolon4siteIdx
	     || idx1!=2*5*NDoublonHolon4siteIdx) {
	    info=ReadDefFileError(defname);
	  }
	}
	fclose(fp);
	break;

      case KWOrbital:
	/*orbitalidx.def------------------------------------*/
	if(NOrbitalIdx>0){
	  idx0 = idx1 = 0;
	  if(APFlag==0) {
	    while( fscanf(fp, "%d %d ", &i, &j) != EOF){
	      fscanf(fp, "%d\n", &(OrbitalIdx[i][j]));
	      OrbitalSgn[i][j] = 1;
	      idx0++;
	      if(idx0==Nsite*Nsite) break;
	    }
	  } else { /* anti-periodic boundary mode */
	    while( fscanf(fp, "%d %d ", &i, &j) != EOF){
	      fscanf(fp, "%d %d\n", &(OrbitalIdx[i][j]), &(OrbitalSgn[i][j]));
	      idx0++;
	      if(idx0==Nsite*Nsite) break;
	    }
	  }

        fidx=NProj;
        while( fscanf(fp, "%d ", &i) != EOF){
	    fscanf(fp, "%d\n", &(OptFlag[2*fidx]));
	    OptFlag[2*fidx+1] = iComplexFlgOrbital; //  TBC imaginary
	    fidx ++;
	    idx1++;
        count_idx++;
        }
	  if(idx0!=Nsite*Nsite || idx1!=NOrbitalIdx) {
	    info=ReadDefFileError(defname);
	  }
	}
	fclose(fp);
	break;

      case KWTransSym:
	/*qptransidx.def------------------------------------*/
	if(NQPTrans>0){
	  for(i=0;i<NQPTrans;i++){
	    fscanf(fp, "%d ", &itmp);
	    fscanf(fp, "%lf %lf\n", &dReValue, &dImValue);
		ParaQPTrans[itmp]=dReValue+I*dImValue;
	  }
	  idx = 0;
	  if(APFlag==0) {
	    while( fscanf(fp, "%d %d ", &i, &j) != EOF){
	      fscanf(fp, "%d\n", &(QPTrans[i][j]));
	      QPTransSgn[i][j] = 1;
	      idx++;
	    }
	  } else { /* anti-periodic boundary mode */
	    while( fscanf(fp, "%d %d ", &i, &j) != EOF){
	      fscanf(fp, "%d %d\n", &(QPTrans[i][j]), &(QPTransSgn[i][j]));
	      idx++;
	    }
	  }
	  if(idx!=NQPTrans*Nsite) info=ReadDefFileError(defname);
	}
	fclose(fp);
	break;

      case KWOneBodyG:
	/*cisajs.def----------------------------------------*/
	if(NCisAjs>0){
	  idx = 0;
	  while( fscanf(fp, "%d %d %d %d\n",
			&(x0), &(x1), &(x2), &(x3)) != EOF){
	    CisAjsIdx[idx][0] = x0;
        CisAjsIdx[idx][1] = x1;
	    CisAjsIdx[idx][2] = x2;
	    CisAjsIdx[idx][3] = x3;
        if(x1 != x3){
          fprintf(stderr, "  Error:  Sz non-conserved system is not yet supported in mVMC ver.1.0.\n");
          info = ReadDefFileError(defname);
          break;
        }
	    idx++;
	  }
	  if(idx!=NCisAjs) info=ReadDefFileError(defname);
	}
	fclose(fp);
	break;	  

    
      case KWTwoBodyGEx:
        /*cisajscktalt.def----------------------------------*/
	if(NCisAjsCktAlt>0){
	  idx = 0;
	  while( fscanf(fp, "%d %d %d %d %d %d %d %d\n", 
			&(x0), &(x1), &(x2), &(x3), &(x4),
			&(x5), &(x6), &(x7) ) != EOF ){
	    CisAjsCktAltIdx[idx][0] = x0;
	    CisAjsCktAltIdx[idx][1] = x1;
	    CisAjsCktAltIdx[idx][2] = x2;
	    CisAjsCktAltIdx[idx][3] = x3;
	    CisAjsCktAltIdx[idx][4] = x4;
	    CisAjsCktAltIdx[idx][5] = x5;
	    CisAjsCktAltIdx[idx][6] = x6;
	    CisAjsCktAltIdx[idx][7] = x7;
	    idx++;
	  }
	  if(idx!=NCisAjsCktAlt) info=ReadDefFileError(defname);
	}
	fclose(fp);
	break;

      case KWTwoBodyG:	  
	/*cisajscktaltdc.def--------------------------------*/
	if(NCisAjsCktAltDC>0){
	  idx = 0;	  
	  while( fscanf(fp, "%d %d %d %d %d %d %d %d\n", 
			&(x0), &(x1), &(x2), &(x3), &(x4),
			&(x5), &(x6), &(x7) ) != EOF ){
	    CisAjsCktAltDCIdx[idx][0] = x0;
	    CisAjsCktAltDCIdx[idx][1] = x1;
	    CisAjsCktAltDCIdx[idx][2] = x2;
	    CisAjsCktAltDCIdx[idx][3] = x3;
	    CisAjsCktAltDCIdx[idx][4] = x4;
	    CisAjsCktAltDCIdx[idx][5] = x5;
	    CisAjsCktAltDCIdx[idx][6] = x6;
	    CisAjsCktAltDCIdx[idx][7] = x7;
	    idx++;
        if(x1 != x3 || x5 != x7){
          fprintf(stderr, "  Error:  Sz non-conserved system is not yet supported in mVMC ver.1.0.\n");
          info = ReadDefFileError(defname);
          break;
        }
	  }
	  if(idx!=NCisAjsCktAltDC) info=ReadDefFileError(defname);
	}
	fclose(fp);
	break;

      case KWInterAll:
	/*interall.def---------------------------------------*/
	if(NInterAll>0){
	  idx = 0;
	  while( fscanf(fp, "%d %d %d %d %d %d %d %d %lf %lf\n",
                    &(InterAll[idx][0]),
                    &(InterAll[idx][1]),//ispin1
                    &(InterAll[idx][2]),
                    &(InterAll[idx][3]),//ispin2
                    &(InterAll[idx][4]),
                    &(InterAll[idx][5]),//ispin3
                    &(InterAll[idx][6]),
                    &(InterAll[idx][7]),//ispin4
                    &dReValue,
					&dImValue)!=EOF ){

		  ParaInterAll[idx]=dReValue+I*dImValue;

		  if(!((InterAll[idx][1] == InterAll[idx][3]
		      || InterAll[idx][5] == InterAll[idx][7])
		   ||
		       (InterAll[idx][1] == InterAll[idx][5]
			|| InterAll[idx][3]  == InterAll[idx][7])
		       )
		     )
		    {
		      fprintf(stderr, "  Error:  Sz non-conserved system is not yet supported in mVMC ver.1.0.\n");
		      info = ReadDefFileError(defname);
		      break;
		    }
		  idx++;
	  }
	  if(idx!=NInterAll) info=ReadDefFileError(defname);
	} else {
	  /* do not terminate */
	  /* info=ReadDefFileError(xNameListFile); */
	}
	fclose(fp);
	break;

    case KWOptTrans:
	/*qpopttrans.def------------------------------------*/
	if(FlagOptTrans>0) {
        fidx=NProj+NOrbitalIdx;
        for(i=0;i<NQPOptTrans;i++){
	    fscanf(fp, "%d ", &itmp);
	    fscanf(fp, "%lf\n", &(ParaQPOptTrans[itmp]));
	    OptFlag[fidx] = 1;
	    fidx ++;
        count_idx++;
        }
	  idx = 0;
	  if(APFlag==0) {
	    while( fscanf(fp, "%d %d ", &i, &j) != EOF){
	      fscanf(fp, "%d\n", &(QPOptTrans[i][j]));
	      QPOptTransSgn[i][j] = 1;
	      idx++;
	    }
	  } else { // anti-periodic boundary mode 
	    while( fscanf(fp, "%d %d ", &i, &j) != EOF){
	      fscanf(fp, "%d %d\n", &(QPOptTrans[i][j]), &(QPOptTransSgn[i][j]));
	      idx++;
	    }
	  }
	}
    
	fclose(fp);
	break;
	  
      default:
	fclose(fp);
	break;
      }
    }

    if(count_idx!=NPara){
      fprintf(stderr, "error: OptFlag is incomplete.\n");
      info=1;
    }
  } /* if(rank==0) */

  if(FlagOptTrans<=0){
    ParaQPOptTrans[0]=1.0;
    for(i=0;i<Nsite;++i) {
      QPOptTrans[0][i] = i;
      QPOptTransSgn[0][i] = 1;
    }
  }  
  if(info!=0) {
    if(rank==0) {
      fprintf(stderr, "error: Indices and Parameters of Definition files(*.def) are incomplete.\n");
    }
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  
#ifdef _mpi_use
  SafeMpiBcastInt(LocSpn, NTotalDefInt, comm);
  SafeMpiBcast_fcmp(ParaTransfer, NTransfer+NInterAll, comm);
  SafeMpiBcast(ParaCoulombIntra, NTotalDefDouble, comm);
  SafeMpiBcast_fcmp(ParaQPTrans, NQPTrans, comm);

  /* MPI_Bcast(LocSpn, NTotalDefInt, MPI_INT, 0, comm); */
  /* MPI_Bcast(ParaTransfer, NTotalDefDouble, MPI_DOUBLE, 0, comm); */
#endif /* _mpi_use */
  
  /* set FlagShift */
  if(NVMCCalMode==0) {
    SetFlagShift();
    if(rank==0 && FlagShiftGJ+FlagShiftDH2+FlagShiftDH4>0 ) {
      fprintf(stderr, "remark: FlagShift ( ");
      if(FlagShiftGJ==1)  fprintf(stderr, "GJ ");
      if(FlagShiftDH2==1) fprintf(stderr, "DH2 ");
      if(FlagShiftDH4==1) fprintf(stderr, "DH4 ");
      fprintf(stderr, ") is turned on.\n");
    }
  }
  
  return 0;
}

int ReadInputParameters(char *xNameListFile, MPI_Comm comm)
{
  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax], ctmp2[256];
  int iKWidx=0;
  int i,idx;
  int rank;
  int count=0;
  int info=0;
  double tmp_real, tmp_comp;
  MPI_Comm_rank(comm, &rank);
  
  if(rank==0) {
    for(iKWidx=KWInGutzwiller; iKWidx< KWIdxInt_end; iKWidx++){     
      strcpy(defname, cFileNameListFile[iKWidx]);
      if(strcmp(defname,"")==0) continue;   
      fp = fopen(defname, "r");
      if(fp==NULL){
        info= ReadDefFileError(defname);
        fclose(fp);
        continue;
      }
      /*=======================================================================*/
      idx=0;
		fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
		fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
		sscanf(ctmp2,"%s %d\n", ctmp, &idx);

		switch(iKWidx){
        //get idx

      case KWInGutzwiller:
        if(idx != NGutzwillerIdx){
          info=1;
          continue;
        }
        for(i=0; i<NGutzwillerIdx; i++){
          fscanf(fp, "%d %lf %lf ", &idx, &tmp_real,&tmp_comp);
          Proj[i]=tmp_real+I*tmp_comp;
        }
        break;
        
      case KWInJastrow:
        if(idx != NJastrowIdx){
          info=1;
          continue;
        }
        count= NGutzwillerIdx;
        for(i=count; i<count+NJastrowIdx; i++){
          fscanf(fp, "%d %lf %lf ", &idx, &tmp_real,&tmp_comp);
          Proj[i]=tmp_real+I*tmp_comp;
        }
        break;
        
      case KWInDH2:
        if(idx != NDoublonHolon2siteIdx){
          info=1;
          continue;
        }
        count =NGutzwillerIdx+NJastrowIdx;
        for(i=count; i<count+2*3*NDoublonHolon2siteIdx; i++){
          fscanf(fp, "%d %lf %lf ", &idx, &tmp_real,&tmp_comp);
          Proj[i]=tmp_real+I*tmp_comp;
        }       
        break;
        
      case KWInDH4:
        if(idx != NDoublonHolon4siteIdx){
          info=1;
          continue;
        }
        count =NGutzwillerIdx+NJastrowIdx+2*3*NDoublonHolon2siteIdx;
        for(i=count; i<count+2*5*NDoublonHolon4siteIdx; i++){
          fscanf(fp, "%d %lf %lf ", &idx, &tmp_real,&tmp_comp);
          Proj[i]=tmp_real+I*tmp_comp;
        }        
        break;
        
      case KWInOrbital:
        if(idx != NOrbitalIdx){
          info=1;
          continue;
        }
        for(i=0; i<NSlater; i++){
          fscanf(fp, "%d %lf %lf ", &idx, &tmp_real,&tmp_comp);
          Slater[i]=tmp_real+I*tmp_comp;
        }  
        break;
        
      case KWInOptTrans:
        if(idx != NOptTrans){
          info=1;
          continue;
        }
        for(i=0; i<NOptTrans; i++){
          fscanf(fp, "%d %lf %lf ", &idx, &tmp_real,&tmp_comp);
          OptTrans[i]=tmp_real+I*tmp_comp;
        }  
        break;
        
      default:
        break;
      }
      fclose(fp);
    }
  } /* if(rank==0) */

  if(info!=0) {
    if(rank==0) {
      fprintf(stderr, "error: Indices of %s file is incomplete.\n", defname);
    }
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }  
  return 0;
}



/**********************************************************************/
/* Function checking keywords*/
/**********************************************************************/
/** 
 * 
 * @brief function of checking whether ctmp is same as cKeyWord or not
 * 
 * @param[in] ctmp 
 * @param[in] cKeyWord 
 * @return 0 ctmp is same as cKeyWord
 * 
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int CheckWords(
	       const char* ctmp,
	       const char* cKeyWord
	       )
{

  int i=0;

  char ctmp_small[256]={0};
  char cKW_small[256]={0};
  int n;
  n=strlen(cKeyWord);
  strncpy(cKW_small, cKeyWord, n);
  
  for(i=0; i<n; i++){
    cKW_small[i]=tolower(cKW_small[i]);
  }
  
  n=strlen(ctmp);
  strncpy(ctmp_small, ctmp, n);
  for(i=0; i<n; i++){
    ctmp_small[i]=tolower(ctmp_small[i]);
  }
  if(n<strlen(cKW_small)) n=strlen(cKW_small);
  return(strncmp(ctmp_small, cKW_small, n));
}

/**
 * @brief Function of Checking keyword in NameList file.
 * @param[in] _cKW keyword candidate
 * @param[in] _cKWList Reffercnce of keyword List
 * @param[in] iSizeOfKWidx number of keyword
 * @param[out] _iKWidx index of keyword
 * @retval 0 keyword is correct.
 * @retval -1 keyword is incorrect.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckKW(
	    const char* cKW,
	    char  cKWList[][D_CharTmpReadDef],
	    int iSizeOfKWidx,
	    int* iKWidx
	    ){
  *iKWidx=-1;
  int itmpKWidx;
  int iret=-1;
  for(itmpKWidx=0; itmpKWidx<iSizeOfKWidx; itmpKWidx++){
    if(strcmp(cKW,"")==0){
      break;
    }
    else if(CheckWords(cKW, cKWList[itmpKWidx])==0){
      iret=0;
      *iKWidx=itmpKWidx;
    }
  }
  return iret;
}

/**
 * @brief Function of Validating value.
 * @param[in] icheckValue value to validate.
 * @param[in] ilowestValue lowest value which icheckValue can be set.
 * @param[in] iHighestValue heighest value which icheckValue can be set.
 * @retval 0 value is correct.
 * @retval -1 value is incorrect.
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int ValidateValue(
		  const int icheckValue, 
		  const int ilowestValue, 
		  const int iHighestValue
		  ){

  if(icheckValue < ilowestValue || icheckValue > iHighestValue){
    return(-1);
  }
  return 0;
}

/**
 * @brief Function of Getting keyword and it's variable from characters.
 * @param[in] _ctmpLine characters including keyword and it's variable 
 * @param[out] _ctmp keyword
 * @param[out] _itmp variable for a keyword
 * @retval 0 keyword and it's variable are obtained.
 * @retval 1 ctmpLine is a comment line.
 * @retval -1 format of ctmpLine is incorrect.
 * @version 1.0
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int GetKWWithIdx(
		 char *ctmpLine,
		 char *ctmp,
		 int *itmp
		 )
{
  char *ctmpRead;
  char *cerror;
  char csplit[] = " ,.\t\n";
  if(*ctmpLine=='\n') return(-1);
  ctmpRead = strtok(ctmpLine, csplit);
  if(strncmp(ctmpRead, "=", 1)==0 || strncmp(ctmpRead, "#", 1)==0 || ctmpRead==NULL){
    return(-1);
  }
  strcpy(ctmp, ctmpRead);
    
  ctmpRead = strtok( NULL, csplit );
  *itmp = strtol(ctmpRead, &cerror, 0);
  //if ctmpRead is not integer type
  if(*cerror != '\0'){
    fprintf(stderr, "Error: incorrect format= %s. \n", cerror);
    return(-1);
  }

  ctmpRead = strtok( NULL, csplit );
  if(ctmpRead != NULL){
    fprintf(stderr, "Error: incorrect format= %s. \n", ctmpRead);
    return(-1);
  }    
  return 0;
}

/**
 * @brief Function of Fitting FileName
 * @param[in]  _cFileListNameFile file for getting names of input files.
 * @param[out] _cFileNameList arrays for getting names of input files.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 1.0
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int GetFileName(
		const char* cFileListNameFile,
		char cFileNameList[][D_CharTmpReadDef]
		)
{
  int myrank;
  FILE *fplist;
  char ctmp[D_FileNameMax];
  int itmpKWidx=-1;
  char ctmpFileName[D_FileNameMaxReadDef];
  char ctmpKW[D_CharTmpReadDef], ctmp2[256];
  int i;
  for(i=0; i< KWIdxInt_end; i++){
    strcpy(cFileNameList[i],"");
  }
  
  fplist = fopen(cFileListNameFile, "r");
  if(fplist==NULL) return ReadDefFileError(cFileListNameFile);
  
  while(fgets(ctmp2, 256, fplist) != NULL){ 
    sscanf(ctmp2,"%s %s\n", ctmpKW, ctmpFileName);

    if(strncmp(ctmpKW, "#", 1)==0 || *ctmp2=='\n'){
      continue;
    }
    else if(strcmp(ctmpKW, "")*strcmp(ctmpFileName, "")==0){
      fprintf(stderr, "Error: keyword and filename must be set as a pair in %s.\n", cFileListNameFile);
      fclose(fplist);
      return(-1);
    }
    /*!< Check KW */
    if( CheckKW(ctmpKW, cKWListOfFileNameList, KWIdxInt_end, &itmpKWidx)!=0 ){

      fprintf(stderr, "Error: Wrong keywords '%s' in %s.\n", ctmpKW, cFileListNameFile);
      fprintf(stderr, "%s", "Choose Keywords as follows: \n");
      for(i=0; i<KWIdxInt_end;i++){
        fprintf(stderr, "%s \n", cKWListOfFileNameList[i]);
      }      
      fclose(fplist);
      return(-1);
    }
    /*!< Check cFileNameList to prevent from double registering the file name */    
    if(strcmp(cFileNameList[itmpKWidx], "") !=0){
      fprintf(stderr, "Error: Same keywords exist in %s.\n", cFileListNameFile);
      fclose(fplist);
      return(-1);
    }
    /*!< Copy FileName */
    strcpy(cFileNameList[itmpKWidx], ctmpFileName);
  }
  fclose(fplist);  
  return 0;
}

void SetDefultValuesModPara(int *bufInt, double* bufDouble){
  
  bufInt[IdxVMCCalcMode]=0;
  bufInt[IdxLanczosMode]=0;
  bufInt[IdxDataIdxStart]=0;
  bufInt[IdxDataQtySmp]=1;
  bufInt[IdxNsite]=16;
  bufInt[IdxNe]=8;
  bufInt[IdxSPGaussLeg]=1;
  bufInt[IdxSPStot]=0;
  bufInt[IdxMPTrans]=0;
  bufInt[IdxSROptItrStep]=1000;
  bufInt[IdxSROptItrSmp]=bufInt[IdxSROptItrStep]/10;
  bufInt[IdxVMCWarmUp]=10;
  bufInt[IdxVMCIniterval]=1;
  bufInt[IdxVMCSample]=10;
  bufInt[IdxExUpdatePath]=0;
  bufInt[IdxRndSeed]=11272;
  bufInt[IdxSplitSize]=1;
  bufInt[IdxNLocSpin]=0;
  bufInt[IdxNTrans]=0;
  bufInt[IdxNCoulombIntra]=0;
  bufInt[IdxNCoulombInter]=0;
  bufInt[IdxNHund]=0;
  bufInt[IdxNPairHop]=0;
  bufInt[IdxNExchange]=0;
  bufInt[IdxNGutz]=0;
  bufInt[IdxNJast]=0;
  bufInt[IdxNDH2]=0;
  bufInt[IdxNDH4]=0;
  bufInt[IdxNOrbit]=0;
  bufInt[IdxNQPTrans]=0;
  bufInt[IdxNOneBodyG]=0;
  bufInt[IdxNTwoBodyG]=0;
  bufInt[IdxNTwoBodyGEx]=0;
  bufInt[IdxNInterAll]=0;
  bufInt[IdxNQPOptTrans]=1;
  
  bufDouble[IdxSROptRedCut]=0.001;
  bufDouble[IdxSROptStaDel]=0.02;
  bufDouble[IdxSROptStepDt]=0.02;
  NStoreO=1;
}

/**********************************************************************/
