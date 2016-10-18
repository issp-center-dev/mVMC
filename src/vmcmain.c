/*-------------------------------------------------------------
 * Variational Monte Carlo
 * main program
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
/* #include "fjcoll.h" */
#include "vmcmain.h"

int VMCParaOpt(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2);
int VMCPhysCal(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2);
void outputData();
void printUsageError();
void printOption();
void initMultiDefMode(int nMultiDef, char *fileDirList, MPI_Comm comm_parent, MPI_Comm *comm_child1);

/*main program*/
int main(int argc, char* argv[])
{
  /* input file name */
  char *fileDefList;
  char *fileInitPara;

  int flagReadInitPara=0;
  int info=0;

  /* for MultiDef mode (-m option) */
  int flagMultiDef=0;
  int nMultiDef = 1;
  /* for getopt() */
  int option;
  extern char *optarg;
  extern int optind,opterr,optopt;
  /* for strtol() */
  extern int errno;
  char *endptr;
  long num;

  /* for MPI */
  int rank0=0,size0=1;
  int group1=0,group2=0,rank1=0,rank2=0,size1=1,size2=1;
  MPI_Comm comm0,comm1,comm2;

  MPI_Init(&argc, &argv);
  NThread = omp_get_max_threads();

  InitTimer();
  StartTimer(0);
  StartTimer(1);
  StartTimer(10);

  /* read options */
  while((option=getopt(argc,argv,"bhvm:oF:"))!=-1) {
    switch(option) {
    case 'b': /* BinaryMode */
      FlagBinary=1;
      break;

    case 'h': /* Print Help Message*/
      printUsageError();
      printOption();
      exit(EXIT_SUCCESS);
      break;

    case 'v': /* Print Version*/
      printVersion();
      exit(EXIT_SUCCESS);
      break;

    case 'm': /* MultiDefMode */
      errno = 0;
      num = strtol(optarg,&endptr,10);
      if((errno == ERANGE && (num == LONG_MIN || num == LONG_MAX)) ||
          (errno != 0 && num == 0)) {
        perror("error: -m: strtol()");
        exit(EXIT_FAILURE);
      }
      if(endptr == optarg) {
        fprintf(stderr,"error: -m: No digits were found\n");
        exit(EXIT_FAILURE);
      }
      if(*endptr != '\0') {
        fprintf(stderr,"warning: -m: Futher characters after number: %s\n",endptr);
      }
      if(num > INT_MAX || num < INT_MIN) {
        fprintf(stderr,"error: -m: Numerical result out of range\n");
        exit(EXIT_FAILURE);
      }
      /* strtol() successfully parsed a number */
      flagMultiDef = 1;
      nMultiDef = (int)num;
      break;

    case 'o': /* OptTransMode */
      FlagOptTrans=1;
      break;

    case 'F': /* Flush output file*/
      errno = 0;
      num = strtol(optarg,&endptr,10);
      if((errno == ERANGE && (num == LONG_MIN || num == LONG_MAX)) ||
          (errno != 0 && num == 0)) {
        perror("error: -F: strtol()");
        exit(EXIT_FAILURE);
      }
      if(endptr == optarg) {
        fprintf(stderr,"error: -F: No digits were found\n");
        exit(EXIT_FAILURE);
      }
      if(*endptr != '\0') {
        fprintf(stderr,"warning: -F: Futher characters after number: %s\n",endptr);
      }
      if(num > INT_MAX || num < INT_MIN) {
        fprintf(stderr,"error: -F: Numerical result out of range\n");
        exit(EXIT_FAILURE);
      }
      /* strtol() successfully parsed a number */
      if(num < 1) {
        fprintf(stderr,"error: -F: FileFlushInterval should be natural number.\n");
        exit(EXIT_FAILURE);
      }
      NFileFlushInterval = (int)num;
      break;

    default: /* '?' */
      printUsageError();
      exit(EXIT_FAILURE);
    }
  }

  /* check the number of arguments */
  if((flagMultiDef==0 && argc-optind<1) || (flagMultiDef==1 && argc-optind<2)) {
    fprintf(stderr,"error: Argument count mismatch\n");
    printUsageError();
    exit(EXIT_FAILURE);
  }

  /* set input filename */
  if(flagMultiDef==0) { /* Original mode */
    fileDefList=argv[optind];
    if(argc-optind>1) {
      flagReadInitPara = 1;
      fileInitPara=argv[optind+1];
    }
  } else if(flagMultiDef==1) { /* MultiDef mode */
    fileDefList=argv[optind+1];
    if(argc-optind>2) {
      flagReadInitPara = 1;
      fileInitPara=argv[optind+2];
    }
  }

  if(flagMultiDef==0) { /* Original mode */
    MPI_Comm_dup(MPI_COMM_WORLD,&comm0);
  } else if(flagMultiDef==1) { /* MultiDef mode */
    /* set communicator and change directory */
    initMultiDefMode(nMultiDef,argv[optind],MPI_COMM_WORLD,&comm0);
  }

  MPI_Comm_rank(comm0, &rank0);
  MPI_Comm_size(comm0, &size0);
  StopTimer(10);

  StartTimer(11);
  if(rank0==0) fprintf(stdout,"Start: Read *def files.\n");
  ReadDefFileNInt(fileDefList, comm0);
  if(rank0==0) fprintf(stdout,"End  : Read *def files.\n");
  StopTimer(11);
  
  StartTimer(12);
  SetMemoryDef();
  StopTimer(12);
  
  StartTimer(11);
  if(rank0==0) fprintf(stdout,"Start: Read parameters from *def files.\n");
  ReadDefFileIdxPara(fileDefList, comm0);
  if(rank0==0) fprintf(stdout,"End  : Read parameters from *def files.\n");
  StopTimer(11);
  
  StartTimer(12);
  if(rank0==0) fprintf(stdout,"Start: Set memories.\n");
  SetMemory();
  if(rank0==0) fprintf(stdout,"End  : Set memories.\n");
  StopTimer(12);
  
  /* split MPI coummunicator */
#ifdef _mpi_use
  StartTimer(10);
  group1 = rank0/NSplitSize;
  MPI_Comm_split(comm0,group1,rank0,&comm1);
  MPI_Comm_size(comm1,&size1);
  MPI_Comm_rank(comm1,&rank1);
  group2 = rank1;
  MPI_Comm_split(comm0,group2,rank0,&comm2);
  MPI_Comm_size(comm2,&size2);
  MPI_Comm_rank(comm2,&rank2);

  if(size0%NSplitSize!=0 && rank0==0) {
    fprintf(stderr,"warning: load imbalance. MPI_size0=%d NSplitSize=%d\n",size0,NSplitSize);
  }
  /*   printf("rank=%d group1=%d rank1=%d rank2=%d size1=%d size2=%d\n", */
  /*      rank,group1,rank1,rank2,size1,size2); */
  StopTimer(10);
#endif

  /* initialize Mersenne Twister */
  init_gen_rand(RndSeed+group1);
  /* get the size of work space for LAPACK and PFAPACK */
  LapackLWork = getLWork_fcmp(); //TBC

  StartTimer(13);
  /* initialize variational parameters */
  if(rank0==0) fprintf(stdout,"Start: Initialize parameters.\n");
  InitParameter(); /* Run parallelly for synchronization of random generator */
  if(flagReadInitPara>0 && rank0==0) ReadInitParameter(fileInitPara);
  //[s] add read parameters respectively
  if(rank0==0){
    if(ReadInputParameters(fileDefList, comm0)!=0){
      //[ToDo]: Add Error procedure
      info=1;
    }
  }
  if(rank0==0) fprintf(stdout,"End  : Initialize parameters.\n");
  //[e] add read parameters respectively
  
  SyncModifiedParameter(comm0);
  StopTimer(13);

  /* initialize variables for quantum projection */
  if(rank0==0) fprintf(stdout,"Start: Initialize variables for quantum projection.\n");
  InitQPWeight();
  if(rank0==0) fprintf(stdout,"End  : Initialize variables for quantum projection.\n");
  /* initialize output files */
  if(rank0==0) InitFile(fileDefList, rank0);

  StopTimer(1);

  if(NVMCCalMode==0) {
    StartTimer(2);
    /*-- VMC Parameter Optimization --*/
    if(rank0==0) fprintf(stdout,"Start: Optimize VMC parameters.\n");
    VMCParaOpt(comm0, comm1, comm2);
    if(rank0==0) fprintf(stdout,"End  : Optimize VMC parameters.\n");
    StopTimer(2);
  } else if(NVMCCalMode==1) {
    StartTimer(2);
    /*-- VMC Physical Quantity Calculation --*/
    if(rank0==0) fprintf(stdout,"Start: Calculate VMC physical quantities.\n");
    VMCPhysCal(comm0, comm1, comm2);
    if(rank0==0) fprintf(stdout,"End  : Calculate VMC physical quantities.\n");
    StopTimer(2);
  } else {
    info=1;
    if(rank0==0) fprintf(stderr,"error: NVMCCalMode must be 0 or 1.\n");
  }

  StopTimer(0);
  if(rank0==0) {
    if(NVMCCalMode==0) {
      OutputTimerParaOpt();
    } else if(NVMCCalMode==1) {
      OutputTimerPhysCal();
    } 
  }

  /* close output files */
  if(rank0==0) CloseFile(rank0);

  FreeMemory();
  FreeMemoryDef();
  MPI_Finalize();

  return info;
}

/*-- VMC Parameter Optimization --*/
int VMCParaOpt(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2) {
  int step;
  int info;
  int rank;
  int i,tmp_i;//DEBUG
  MPI_Comm_rank(comm_parent, &rank);

  for(step=0;step<NSROptItrStep;step++) {
    if(rank==0) OutputTime(step);
      StartTimer(20);
    UpdateSlaterElm_fcmp();
    UpdateQPWeight();
      StopTimer(20);
      StartTimer(3);
     if(AllComplexFlag==0){
        // only for real TBC
         StartTimer(69);
         #pragma omp parallel for default(shared) private(tmp_i)
         for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElm_real[tmp_i]= creal(SlaterElm[tmp_i]);
         #pragma omp parallel for default(shared) private(tmp_i)
         for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM_real[tmp_i]= creal(InvM[tmp_i]);
         StopTimer(69);
         // SlaterElm_real will be used in CalculateMAll, note that SlaterElm will not change before SR
         VMCMakeSample_real(comm_child1);
         // only for real TBC
         StartTimer(69);
         #pragma omp parallel for default(shared) private(tmp_i)
         for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM[tmp_i]      = InvM_real[tmp_i]+0.0*I;
         StopTimer(69);
         // only for real TBC
      }else{
        VMCMakeSample(comm_child1);
      } 
      StopTimer(3);
      StartTimer(4);
    VMCMainCal(comm_child1);
      StopTimer(4);
      StartTimer(21);
    WeightAverageWE(comm_parent);
      StartTimer(25);//DEBUG
    if(AllComplexFlag==0){
      WeightAverageSROpt_real(comm_parent);
    }else{
      WeightAverageSROpt(comm_parent);
    }
      StopTimer(25);
    ReduceCounter(comm_child2);
      StopTimer(21);
      StartTimer(22);
    /* output zvo_out and zvo_var */
    if(rank==0) outputData();
      StopTimer(22);
      StartTimer(5);
//DEBUG
    if(rank==0){
      for(i=0;i<2*SROptSize*(2*SROptSize+2);i++){
      //  printf("DEBUG: i=%d SROptOO=%lf +I*%lf\n",i,creal(SROptOO[i]),cimag(SROptOO[i]));
      } 
    }
//DBBUG
    info = StochasticOpt(comm_parent);
    //info = StochasticOptDiag(comm_parent);
    info=0;
      StopTimer(5);

    if(info!=0) {
      if(rank==0) fprintf(stderr, "Error: StcOpt info=%d step=%d\n",info,step);
      return info;
    }

      StartTimer(23);
    SyncModifiedParameter(comm_parent);
      StopTimer(23);

    if(step >= NSROptItrStep-NSROptItrSmp) {
      StoreOptData(step-(NSROptItrStep-NSROptItrSmp));
    }

    FlushFile(step,rank);
  }

  if(rank==0) OutputTime(NSROptItrStep);

  /* output zqp_opt */
  if(rank==0) OutputOptData();

  return 0;
}

/*-- VMC Physical Quantity Calculation --*/
int VMCPhysCal(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2) {
  int ismp;
  int rank;
  MPI_Comm_rank(comm_parent, &rank);

  StartTimer(20);
  UpdateSlaterElm_fcmp();
  StopTimer(20);

  for(ismp=0;ismp<NDataQtySmp;ismp++) {
    if(rank==0) OutputTime(ismp);
    FlushFile(0,rank);

    InitFilePhysCal(ismp, rank);
    
    StartTimer(3);

    VMCMakeSample(comm_child1);

    StopTimer(3);
    StartTimer(4);

    VMCMainCal(comm_child1);

    StopTimer(4);
    StartTimer(21);

    WeightAverageWE(comm_parent);
    WeightAverageGreenFunc(comm_parent);
    ReduceCounter(comm_child2);

    StopTimer(21);
    StartTimer(22);
    /* output zvo_out and green functions */
    if(rank==0) outputData();
    CloseFilePhysCal(rank);

    StopTimer(22);
    StopTimer(5);
  }

  if(rank==0) OutputTime(NDataQtySmp);

  return 0;
}

void outputData() {
  int i;

  /* zvo_out.dat */
 // fprintf(FileOut, "% .18e % .18e % .18e \n", Etot, Etot2, (Etot2 - Etot*Etot)/(Etot*Etot));
    fprintf(FileOut, "% .18e % .18e  % .18e % .18e \n", creal(Etot),cimag(Etot), creal(Etot2), creal((Etot2 - Etot*Etot)/(Etot*Etot)));

  /* zvo_var.dat */
  if(FlagBinary==0) { /* formatted output*/
    fprintf(FileVar, "% .18e % .18e 0.0 % .18e % .18e 0.0 ", creal(Etot), cimag(Etot), creal(Etot2), cimag(Etot2));
    for(i=0;i<NPara;i++)   fprintf(FileVar, "% .18e % .18e 0.0 ", creal(Para[i]),cimag(Para[i]));
    fprintf(FileVar, "\n");
    //for(i=0;i<NPara;i++)  printf("DEBUG:i=%d: % .18e % .18e  \n",i, creal(Para[i]),cimag(Para[i]));
  } else { /* binary output */
    fwrite(Para,sizeof(double),NPara,FileVar);
  }

  if(NVMCCalMode==1) {
    /* zvo_cisajs.dat */
      if(NCisAjs>0) {
          for (i = 0; i < NCisAjs; i++)
              fprintf(FileCisAjs, "% .18e  % .18e ", creal(PhysCisAjs[i]), cimag(PhysCisAjs[i]));
          fprintf(FileCisAjs, "\n");
      }
    /* zvo_cisajscktalt.dat */
      if(NCisAjsCktAlt>0) {
          for (i = 0; i < NCisAjsCktAlt; i++)
              fprintf(FileCisAjsCktAlt, "% .18e  % .18e ", creal(PhysCisAjsCktAlt[i]), cimag(PhysCisAjsCktAlt[i]));
          fprintf(FileCisAjsCktAlt, "\n");
      }

    /* zvo_cisajscktaltdc.dat */
      if(NCisAjsCktAltDC>0) {
          for (i = 0; i < NCisAjsCktAltDC; i++)
              fprintf(FileCisAjsCktAltDC, "% .18e % .18e  ", creal(PhysCisAjsCktAltDC[i]),
                      cimag(PhysCisAjsCktAltDC[i]));
          fprintf(FileCisAjsCktAltDC, "\n");
      }

    if(NLanczosMode>0){
// ignorign Lanczos to be added
      /* zvo_ls.dat */
//      fprintf(FileLS, "% .18e  ", creal(QQQQ[2]));  /* H * I = QQQQ[1],[2],[4],[8] */      //TBC
//      fprintf(FileLS, "% .18e  ", creal(QQQQ[3]));  /* H * H = QQQQ[3],[6],[9],[12] */     //TBC
//      fprintf(FileLS, "% .18e  ", creal(QQQQ[2])); /* H^2 * I = QQQQ[5],[10] */           //TBC
//      fprintf(FileLS, "% .18e  ", creal(QQQQ[11])); /* H^2 * H = QQQQ[7],[11],[13],[14] */ //TBC
//      fprintf(FileLS, "% .18e\n", creal(QQQQ[15])); /* H^2 * H^2 = QQQQ[15] */             //TBC
//
//      /* zvo_ls_qqqq.dat */
//      for(i=0;i<NLSHam*NLSHam*NLSHam*NLSHam;i++) {
//        fprintf(FileLSQQQQ, "% .18e  ", creal(QQQQ[i]));
//      }
//      fprintf(FileLSQQQQ, "\n");
//
//      if(NLanczosMode>1){
//        /* zvo_ls_qcisajsq.dat */
//        for(i=0;i<NLSHam*NLSHam*NCisAjs;i++) {
//          fprintf(FileLSQCisAjsQ, "% .18e  ", QCisAjsQ[i]);
//        }
//        fprintf(FileLSQCisAjsQ, "\n");
//
//        /* zvo_ls_qcisajscktaltq.dat */
//        for(i=0;i<NLSHam*NLSHam*NCisAjsCktAlt;i++) {
//          fprintf(FileLSQCisAjsCktAltQ, "% .18e  ", QCisAjsCktAltQ[i]);
//        }
//        fprintf(FileLSQCisAjsCktAltQ, "\n");
//      }
    }
  }

  return;
}

void printUsageError() {
  fprintf(stderr,"Usage: vmc.out [option] NameListFile [OptParaFile]\n");
  fprintf(stderr,"       vmc.out -m N [option] [--] DirListFile NameListFile [OptParaFile]\n");
  return;
}

void printOption() {
  fprintf(stderr,"  -b     binary mode\n");
  fprintf(stderr,"  -m N   multiDef mode\n");
  fprintf(stderr,"  -o     optTrans mode\n");
  fprintf(stderr,"  -F N   set interval of file flush\n");
  fprintf(stderr,"  -h     show this message\n");
  return;
}

/* This function splits MPI communicator, reads DirName from fileDirList,
   and change current working directory */
void initMultiDefMode(int nMultiDef, char *fileDirList, MPI_Comm comm_parent, MPI_Comm *comm_child1) {
  char dirName[D_FileNameMax];
  char *dirNameList;
  FILE *fp;
  int i;
  int info=0;

  int rank, size;
  int group1, group2, rank1;
  int div, mod, threshold;
  MPI_Comm comm_child2;

  MPI_Comm_rank(comm_parent, &rank);
  MPI_Comm_size(comm_parent, &size);

  /* check MPI size */
  if(size<nMultiDef) {
    if(rank==0) fprintf(stderr,"error: -m: N should be smaller than MPI size.\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else if(size%nMultiDef!=0) {
    if(rank==0) fprintf(stderr,"warning: load imbalance. MPI_size=%d nMultiDef=%d\n",size,nMultiDef);
  }

  /* split MPI communicator */
  div = size / nMultiDef;
  mod = size % nMultiDef;
  threshold = (div+1)*mod;
  if(rank < threshold) {
    group1 = rank / (div+1);
  } else {
    group1 = mod + (rank-threshold)/div;
  }
  MPI_Comm_split(comm_parent,group1,rank,comm_child1);
  MPI_Comm_rank((*comm_child1), &rank1);
  group2 = rank1;
  MPI_Comm_split(comm_parent,group2,rank,&comm_child2);

  /* read fileDirList (only at rank=0 process) */
  if(rank==0) {
    dirNameList = (char*)malloc(nMultiDef*(D_FileNameMax)*sizeof(char));
    if( (fp=fopen(fileDirList, "r")) != NULL ) {
      for(i=0;i<nMultiDef;i++) {
        if(fscanf(fp, "%s\n", (dirNameList + i*D_FileNameMax) )!=1) {
          fprintf(stderr,"error: %s is incomplete.\n",fileDirList);
          info=1;
          break;
        }
      }
      fclose(fp);
    } else {
      fprintf(stderr,"error: DirListFile does not exist.\n");
      info=1;
    }
  }

  /* error handle */
  if(info!=0) {
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

  if(group2==0) { /* rank1==0 */
    /* scatter and broadcast dirName */
    MPI_Scatter(dirNameList,D_FileNameMax,MPI_CHAR,dirName,D_FileNameMax,MPI_CHAR,0,comm_child2);

    /* change current working directory */
    if( chdir(dirName) != 0) {
      /* error handle */
      fprintf(stderr,"error: chdir(): %s: ",dirName);
      perror("");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
  }

  MPI_Comm_free(&comm_child2);
  return;
}
