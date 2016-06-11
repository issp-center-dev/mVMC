/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Read Definition Files
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

int ReadDefFileError(char *defname);
int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm);
int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm);

int ReadDefFileError(char *defname){
  fprintf(stderr, "error: %s (Broken file or Not exist)\n", defname);
  return 1;
}

int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm){
  FILE *fp, *fplist;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  int itmp,tmp_info;

  int rank, info=0;
  const int nBufInt=36;
  const int nBufDouble=3;
  const int nBufChar=D_FileNameMax;
  int bufInt[nBufInt];
  double bufDouble[nBufDouble];

  MPI_Comm_rank(comm, &rank);

  if(rank==0) {
    fplist = fopen(xNameListFile, "r");
    if(fplist!=NULL) {
      /* zmodpara.def */
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &itmp);
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %s\n", ctmp, CDataFileHead);
            fscanf(fp,"%s %s\n", ctmp, CParaFileHead);
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 0])); /* NVMCCalMode */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 1])); /* NLanczosMode */
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 2])); /* NDataIdxStart */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 3])); /* NDataQtySmp */
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 4])); /* Nsite */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 5])); /* Ne */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 6])); /* NSPGaussLeg */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 7])); /* NSPStot */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 8])); /* NMPTrans */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 9])); /* NSROptItrStep */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[10])); /* NSROptItrSmp */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[11])); /* NSROptFixSmp */
            fscanf(fp,"%s %lf\n",ctmp, &(bufDouble[0])); /* DSROptRedCut */
            fscanf(fp,"%s %lf\n",ctmp, &(bufDouble[1])); /* DSROptStaDel */
            fscanf(fp,"%s %lf\n",ctmp, &(bufDouble[2])); /* DSROptStepDt */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[12])); /* NVMCWarmUp */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[13])); /* NVMCIniterval */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[14])); /* NVMCSample */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[15])); /* NExUpdatePath */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[16])); /* RndSeed */
            if(bufInt[16]<0) {
              bufInt[16] = (int)time(NULL);
              fprintf(stderr, "remark: Seed = %d\n", bufInt[16]);
            }
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[17])); /* NSplitSize */
            tmp_info=fscanf(fp,"%s %d\n", ctmp, &(NStoreO)); /* choice of store O: 0-> normal other-> store */
            if(tmp_info==-1){
              NStoreO = 0;
              printf("no input for NStoreO: default value (NStore=0) is used. \n");
            }
            fclose(fp);
          } else { info = ReadDefFileError(defname); }
        } else { info = ReadDefFileError(xNameListFile); }
      }
      /*locspn.def----------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[18])); /* NLocSpn */
            fclose(fp);
          } else { info =  ReadDefFileError(defname); }
        } else { info = ReadDefFileError(xNameListFile); }
      }
      /*transfer.def--------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[19])); /* NTransfer */
            fclose(fp);
          } else { info= ReadDefFileError(defname);}
        } else {info = ReadDefFileError(xNameListFile);}
      }
      /*coulombintra.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[20])); /* NCoulombIntra */
            fclose(fp);
          } else { info =ReadDefFileError(defname);}
        } else { info= ReadDefFileError(xNameListFile);}
      }
      /*coulombinter.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[21])); /* NCoulombInter */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*hund.def------------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[22])); /* NHundCoupling */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*pairhop.def---------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[23])); /* NPairHopping */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*exchange.def--------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[24])); /* NExchangeCoupling */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*gutzwilleridx.def---------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[25])); /* NGutzwillerIdx */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*jastrowidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[26])); /* NJastrowIdx */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*doublonholon2siteidx.def--------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[27])); /* NDoublonHolon2siteIdx */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info =ReadDefFileError(xNameListFile);}
      }
      /*doublonholon4siteidx.def--------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[28])); /* NDoublonHolon4siteIdx */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*orbitalidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[29])); /* NOrbitalIdx */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*qptransidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[30])); /* NQPTrans */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*cisajs.def----------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[31])); /* NCisAjs */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*cisajscktalt.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[32])); /* NCisAjsCktAlt */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*cisajscktaltdc.def--------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[33])); /* NCisAjsCktAltDC */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*interall.def--------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[34])); /* NInterAll */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else {
          ReadDefFileError(xNameListFile);
          fprintf(stderr, "warning: Please add zinterall.def. This job assumes NInterAll=0.\n");
          bufInt[34] = 0;
        }
      }
      /*qpopttrans.def--------------------------------*/
      if(info==0) {
        bufInt[35] = 1;
        if(FlagOptTrans>0) {
          if(fscanf(fplist, "%s\n", defname)!=EOF) {
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              fscanf(fp,"%s %d\n", ctmp, &(bufInt[35])); /* NQPOptTrans */
              fclose(fp);
              if(bufInt[35]<1) {
                fprintf(stderr, "error: NQPOptTrans should be larger than 0.\n");
                info = ReadDefFileError(defname);
              }
            } else { info = ReadDefFileError(defname);}
          } else { info = ReadDefFileError(xNameListFile);}
        }
      }

      fclose(fplist);
    } else { info=ReadDefFileError(xNameListFile);}
  } /* if(rank==0) */

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

  NVMCCalMode            =  bufInt[ 0];
  NLanczosMode           =  bufInt[ 1];
  NDataIdxStart          =  bufInt[ 2];
  NDataQtySmp            =  bufInt[ 3];
  Nsite                  =  bufInt[ 4];
  Ne                     =  bufInt[ 5];
  NSPGaussLeg            =  bufInt[ 6];
  NSPStot                =  bufInt[ 7];
  NMPTrans               =  bufInt[ 8];
  NSROptItrStep          =  bufInt[ 9];
  NSROptItrSmp           =  bufInt[10];
  NSROptFixSmp           =  bufInt[11];
  NVMCWarmUp             =  bufInt[12];
  NVMCIniterval          =  bufInt[13];
  NVMCSample             =  bufInt[14];
  NExUpdatePath          =  bufInt[15];
  RndSeed                =  bufInt[16];
  NSplitSize             =  bufInt[17];
  NLocSpn                =  bufInt[18];
  NTransfer              =  bufInt[19];
  NCoulombIntra          =  bufInt[20];
  NCoulombInter          =  bufInt[21];
  NHundCoupling          =  bufInt[22];
  NPairHopping           =  bufInt[23];
  NExchangeCoupling      =  bufInt[24];
  NGutzwillerIdx         =  bufInt[25];
  NJastrowIdx            =  bufInt[26];
  NDoublonHolon2siteIdx  =  bufInt[27];
  NDoublonHolon4siteIdx  =  bufInt[28];
  NOrbitalIdx            =  bufInt[29];
  NQPTrans               =  bufInt[30];
  NCisAjs                =  bufInt[31];
  NCisAjsCktAlt          =  bufInt[32];
  NCisAjsCktAltDC        =  bufInt[33];
  NInterAll              =  bufInt[34];
  NQPOptTrans            =  bufInt[35];

  DSROptRedCut = bufDouble[0];
  DSROptStaDel = bufDouble[1];
  DSROptStepDt = bufDouble[2];

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
    + 3*NTransfer /* Transfer */
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
    + 3*NCisAjs /* CisAjs */
    + 8*NCisAjsCktAlt /* CisAjsCktAlt */
    + 6*NCisAjsCktAltDC /* CisAjsCktAltDC */
    + 6*NInterAll /* InterAll */
    + Nsite*NQPOptTrans /* QPOptTrans */
    + Nsite*NQPOptTrans /* QPOptTransSgn */
    + 2*NPara; /* OptFlag */ // TBC
  NTotalDefDouble = NTransfer /* ParaTransfer */
    + NCoulombIntra /* ParaCoulombIntra */
    + NCoulombInter /* ParaCoulombInter */
    + NHundCoupling /* ParaHondCoupling */
    + NPairHopping  /* ParaPairHopping */
    + NExchangeCoupling /* ParaExchangeCoupling */
    + NQPTrans /* ParaQPTrans */
    + NInterAll /* ParaInterAll */
    + NQPOptTrans; /* ParaQPTransOpt */

  return 0;
}

int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm){
  FILE *fp, *fplist;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  int itmp;

  int i,j,n,idx,idx0,idx1,info=0;
  int fidx=0; /* index for OptFlag */
  int x0,x1,x2,x3,x4,x5,x6,x7;
  int rank;
  double tmp_real,tmp_comp;

  MPI_Comm_rank(comm, &rank);

  if(rank==0) {
    fplist = fopen(xNameListFile, "r");
    if(fplist!=NULL) {
      /*modpara.def---------------------------------------*/
      fscanf(fplist, "%s\n", defname);

      /*locspn.def----------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            idx = 0;
            while( fscanf(fp, "%d %d\n", &(x0), &(x1) )!=EOF){
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
            if(idx!=Nsite) info=ReadDefFileError(defname);
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }

      /*transfer.def--------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NTransfer>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %lf\n", 
                            &(Transfer[idx][0]),
                            &(Transfer[idx][1]),
                            &(Transfer[idx][2]),
                            &(ParaTransfer[idx]))!=EOF){
                idx++;
              }
              if(idx!=NTransfer) info = ReadDefFileError(defname);
              fclose(fp);
            } else { info = ReadDefFileError(defname);}
          }
        } else { info = ReadDefFileError(xNameListFile);}
      }
    
      /*coulombintra.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NCoulombIntra>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %lf\n", 
                            &(CoulombIntra[idx]),
                            &(ParaCoulombIntra[idx]) )!=EOF){
                idx++;
              }
              if(idx!=NCoulombIntra) info = ReadDefFileError(defname);
              fclose(fp);
            } else { info = ReadDefFileError(defname);}
          }
        } else { info = ReadDefFileError(xNameListFile);}
      }
    
      /*coulombinter.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NCoulombInter>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %lf\n", 
                            &(CoulombInter[idx][0]),
                            &(CoulombInter[idx][1]),
                            &(ParaCoulombInter[idx]) )!=EOF){
                idx++;
              }
              if(idx!=NCoulombInter) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*hund.def------------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NHundCoupling>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %lf\n", 
                            &(HundCoupling[idx][0]),
                            &(HundCoupling[idx][1]),
                            &(ParaHundCoupling[idx]) )!=EOF){
                idx++;
              }
              if(idx!=NHundCoupling) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*pairhop.def---------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NPairHopping>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %lf\n", 
                            &(PairHopping[idx][0]),
                            &(PairHopping[idx][1]),
                            &(ParaPairHopping[idx]) )!=EOF){
                idx++;
              }
              if(idx!=NPairHopping) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*exchange.def--------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NExchangeCoupling>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %lf\n", 
                            &(ExchangeCoupling[idx][0]),
                            &(ExchangeCoupling[idx][1]),
                            &(ParaExchangeCoupling[idx]) )!=EOF){
                idx++;
              }
              if(idx!=NExchangeCoupling) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*gutzwilleridx.def---------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NGutzwillerIdx>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx0 = idx1 = 0;
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(GutzwillerIdx[i]));
                idx0++;
                if(idx0==Nsite) break;
              }
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(OptFlag[2*fidx])); // TBC real
                OptFlag[2*fidx+1] = 0; //  TBC imaginary
                fidx++;
                idx1++;
              }
              if(idx0!=Nsite || idx1!=NGutzwillerIdx) {
                info=ReadDefFileError(defname);
              }
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*jastrowidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NJastrowIdx>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
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
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(OptFlag[2*fidx])); // TBC real
                OptFlag[2*fidx+1] = 0; //  TBC imaginary
                fidx++;
                idx1++;
              }
              if(idx0!=Nsite*(Nsite-1) || idx1!=NJastrowIdx) {
                info=ReadDefFileError(defname);
              }
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*doublonholon2siteidx.def--------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NDoublonHolon2siteIdx>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx0 = idx1 = 0;
              while( fscanf(fp, "%d %d %d %d\n", &i, &(x0), &(x1), &n) != EOF){
                DoublonHolon2siteIdx[n][2*i]   = x0;
                DoublonHolon2siteIdx[n][2*i+1] = x1;
                idx0++;
                if(idx0==Nsite*NDoublonHolon2siteIdx) break;
              }
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(OptFlag[2*fidx]));//TBC real
                OptFlag[2*fidx+1] = 0; //  TBC imaginary
                fidx++;
                idx1++;
              }
              if(idx0!=Nsite*NDoublonHolon2siteIdx
                 || idx1!=2*3*NDoublonHolon2siteIdx) {
                info=ReadDefFileError(defname);
              }
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*doublonholon4siteidx.def--------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NDoublonHolon4siteIdx>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
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
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(OptFlag[2*fidx]));
                OptFlag[2*fidx+1] = 0; //  TBC imaginary
                fidx++;
                idx1++;
              }
              if(idx0!=Nsite*NDoublonHolon4siteIdx
                 || idx1!=2*5*NDoublonHolon4siteIdx) {
                info=ReadDefFileError(defname);
              }
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*orbitalidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NOrbitalIdx>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
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
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(OptFlag[2*fidx]));
                OptFlag[2*fidx+1] = 0; //  TBC imaginary
                fidx += 1;
                idx1++;
              }
              if(idx0!=Nsite*Nsite || idx1!=NOrbitalIdx) {
                info=ReadDefFileError(defname);
              }
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*qptransidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NQPTrans>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              for(i=0;i<NQPTrans;i++){
                fscanf(fp, "%d ",   &itmp);
                fscanf(fp, "%lf ^n",   &(tmp_real));
                //fscanf(fp, "%lf\n", &(tmp_comp));
                ParaQPTrans[itmp] = tmp_real+0*I;
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
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*cisajs.def----------------------------------------*/
      if(info==0){
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NCisAjs>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %d\n",
                            &(x0), &(x1), &(x2), &(x3)) != EOF){
                CisAjsIdx[x0][0] = x1;
                CisAjsIdx[x0][1] = x2;
                CisAjsIdx[x0][2] = x3;
                idx++;
              }
              if(idx!=NCisAjs) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*cisajscktalt.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NCisAjsCktAlt>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
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
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }
    
      /*cisajscktaltdc.def--------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NCisAjsCktAltDC>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %d %d %d\n", 
                            &(x0), &(x1), &(x2), &(x3), &(x4), &(x5) ) != EOF ){
                CisAjsCktAltDCIdx[idx][0] = x0;
                CisAjsCktAltDCIdx[idx][1] = x1;
                CisAjsCktAltDCIdx[idx][2] = x2;
                CisAjsCktAltDCIdx[idx][3] = x3;
                CisAjsCktAltDCIdx[idx][4] = x4;
                CisAjsCktAltDCIdx[idx][5] = x5;
                idx++;
              }
              if(idx!=NCisAjsCktAltDC) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*interall.def---------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NInterAll>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %d %d %d %lf\n", 
                            &(InterAll[idx][0]),
                            &(InterAll[idx][1]),
                            &(InterAll[idx][2]),
                            &(InterAll[idx][3]),
                            &(InterAll[idx][4]),
                            &(InterAll[idx][5]),
                            &(ParaInterAll[idx]) )!=EOF ){
                idx++;
              }
              if(idx!=NInterAll) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else {
          /* do not terminate */
          /* info=ReadDefFileError(xNameListFile); */
        }
      }
    
      /*qpopttrans.def------------------------------------*/
      if(info==0) {
        if(FlagOptTrans>0) {
          if(fscanf(fplist, "%s\n", defname)!=EOF) {
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              for(i=0;i<NQPOptTrans;i++){
                fscanf(fp, "%d ", &itmp);
                fscanf(fp, "%lf\n", &(ParaQPOptTrans[itmp]));
                OptFlag[2*fidx] = 1;
                OptFlag[2*fidx+1] = 0; //  TBC imaginary
                fidx += 1;
              }
              idx = 0;
              if(APFlag==0) {
                while( fscanf(fp, "%d %d ", &i, &j) != EOF){
                  fscanf(fp, "%d\n", &(QPOptTrans[i][j]));
                  QPOptTransSgn[i][j] = 1;
                  idx++;
                }
              } else { /* anti-periodic boundary mode */
                while( fscanf(fp, "%d %d ", &i, &j) != EOF){
                  fscanf(fp, "%d %d\n", &(QPOptTrans[i][j]), &(QPOptTransSgn[i][j]));
                  idx++;
                }
              }
              if(idx!=NQPOptTrans*Nsite) { info=ReadDefFileError(defname);}
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else {
          ParaQPOptTrans[0]=1.0;
          for(i=0;i<Nsite;++i) {
            QPOptTrans[0][i] = i;
            QPOptTransSgn[0][i] = 1;
          }
        }
      }

      fclose(fplist);

    } else { info = ReadDefFileError(xNameListFile);}

    if(fidx!=NPara){
      fprintf(stderr, "error: OptFlag is incomplete.\n");
      info=1;
    }
  } /* if(rank==0) */

  if(info!=0) {
    if(rank==0) {
      fprintf(stderr, "error: Indices and Parameters of Definition files(*.def) are incomplete.\n");
    }
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

#ifdef _mpi_use
  SafeMpiBcastInt(LocSpn, NTotalDefInt, comm);
  SafeMpiBcast(ParaTransfer, NTotalDefDouble, comm);
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
