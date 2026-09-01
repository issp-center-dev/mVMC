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
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * main program
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
/* #include "fjcoll.h" */
#include "vmcmain.h"
#include "physcal_lanczos.h"
#include "physcal_lanczos2.h"
#include "power_lanczos_stabilized.h"

// #define _DEBUG
// #define _DEBUG_DUMP_SROPTO_STORE
// #define _DEBUG_DUMP_SROPTOO
// #define _DEBUG_DUMP_PARA

int VMCParaOpt(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2);
int VMCPhysCal(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2);
void outputData();
void printUsageError();
void printOption();
void initMultiDefMode(int nMultiDef, char *fileDirList, MPI_Comm comm_parent, MPI_Comm *comm_child1);
void StdFace_main(char *fname);

static int RunPowerLanczosStabilized(
    int outputIndex, MPI_Comm communicator) {
  MVMCPowerLanczosClassicView view;
  MVMCPowerLanczosStabilizedInput input;
  MVMCPowerLanczosStabilizedResult result;
  MVMCKrylovStatus status;
  size_t sweep;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &size);
  if (NVMCInterval <= 0 || Nsite <= 0 || NVMCWarmUp < 0 ||
      NVMCSample < 0 ||
      (size_t)NVMCInterval > SIZE_MAX / (size_t)Nsite) {
    if (rank == 0) {
      fprintf(stderr,
              "Error: invalid stabilized power-Lanczos chain controls.\n");
    }
    return 1;
  }
  if (mvmc_power_lanczos_stabilized_block_count(
          (size_t)NVMCSample) == 0) {
    if (rank == 0) {
      fprintf(stderr,
              "Error: stabilized power-Lanczos requires NVMCSample >= 8 "
              "and a divisor between 4 and 16 with at least two samples "
              "per block "
              "(a multiple of 16 is recommended; got %d).\n",
              NVMCSample);
    }
    return 1;
  }
  sweep = (size_t)NVMCInterval * (size_t)Nsite;
  if ((size_t)NVMCWarmUp > SIZE_MAX / sweep) {
    if (rank == 0) {
      fprintf(stderr,
              "Error: stabilized power-Lanczos warm-up overflows size_t.\n");
    }
    return 1;
  }
  memset(&view, 0, sizeof(view));
  view.site_count = Nsite;
  view.up_electron_count = Ne;
  view.down_electron_count = Ne;
  view.pure_spin = NExUpdatePath == 2;
  view.arithmetic = AllComplexFlag == 0
                        ? MVMC_POWER_LANCZOS_CLASSIC_REAL
                        : MVMC_POWER_LANCZOS_CLASSIC_COMPLEX;
  view.transfer_count = NTransfer;
  view.transfer_indices = rank == 0 ? Transfer : NULL;
  view.transfer_parameters = rank == 0 ? ParaTransfer : NULL;
  view.coulomb_intra_count = NCoulombIntra;
  view.coulomb_intra_indices = rank == 0 ? CoulombIntra : NULL;
  view.coulomb_intra_parameters =
      rank == 0 ? ParaCoulombIntra : NULL;
  view.coulomb_inter_count = NCoulombInter;
  view.coulomb_inter_indices = rank == 0 ? CoulombInter : NULL;
  view.coulomb_inter_parameters =
      rank == 0 ? ParaCoulombInter : NULL;
  view.hund_count = NHundCoupling;
  view.hund_indices = rank == 0 ? HundCoupling : NULL;
  view.hund_parameters = rank == 0 ? ParaHundCoupling : NULL;
  view.exchange_count = NExchangeCoupling;
  view.exchange_indices = rank == 0 ? ExchangeCoupling : NULL;
  view.exchange_parameters =
      rank == 0 ? ParaExchangeCoupling : NULL;
  view.pair_hopping_count = NPairHopping;
  view.inter_all_count = NInterAll;
  view.nbody_inter_all_count = NNBodyInterAll;
  view.nbody_g_count = NNBodyG;
  view.qp_total = NQPFull;
  view.qp_start = NQPFull * rank / size;
  view.qp_end = NQPFull * (rank + 1) / size;
  view.scaled_pivot_tolerance = 0.0;
  view.nproj = NProj;
  view.ngutzwiller_idx = NGutzwillerIdx;
  view.njastrow_idx = NJastrowIdx;
  view.nspin_jastrow_idx = NSpinJastrowIdx;
  view.ndoublon_holon_2site_idx = NDoublonHolon2siteIdx;
  view.ndoublon_holon_4site_idx = NDoublonHolon4siteIdx;
  view.gutzwiller_idx = GutzwillerIdx;
  view.jastrow_idx = JastrowIdx;
  view.spin_jastrow_idx = SpinJastrowIdx;
  view.doublon_holon_2site_idx = DoublonHolon2siteIdx;
  view.doublon_holon_4site_idx = DoublonHolon4siteIdx;
  view.projection_parameters = Proj;
  view.qp_weights = QPFullWeight;
  view.slater_real = AllComplexFlag == 0 ? SlaterElm_real : NULL;
  view.slater_complex = AllComplexFlag == 0 ? NULL : SlaterElm;
  if (FlagRBM != 0) {
    view.unsupported_amplitude_features |=
        MVMC_CLASSIC_KRYLOV_UNSUPPORTED_RBM;
  }
  if (NProjBF != 0) {
    view.unsupported_amplitude_features |=
        MVMC_CLASSIC_KRYLOV_UNSUPPORTED_BACKFLOW;
  }
  if (iFlgOrbitalGeneral != 0) {
    view.unsupported_amplitude_features |=
        MVMC_CLASSIC_KRYLOV_UNSUPPORTED_GENERAL_ORBITAL;
  }
#ifdef _pf_block_update
  view.unsupported_amplitude_features |=
      MVMC_CLASSIC_KRYLOV_UNSUPPORTED_BLOCK_UPDATE;
#endif
  memset(&input, 0, sizeof(input));
  input.classic_view = &view;
  input.communicator = communicator;
  input.power_step = NLanczosStep;
  input.seed =
      ((uint64_t)(uint32_t)RndSeed << 32) ^
      (uint64_t)(unsigned)(outputIndex + 1) ^
      UINT64_C(0x504c534c494d5631);
  if (input.seed == 0) input.seed = UINT64_C(1);
  input.warm_up = (size_t)NVMCWarmUp * sweep;
  input.sample_count = (size_t)NVMCSample;
  input.interval = sweep;
  status = mvmc_power_lanczos_stabilized_run(&input, &result);
  if (rank == 0 && status == MVMC_KRYLOV_STATUS_OK) {
    fprintf(FileLS,
            "# stabilized_power_lanczos version=%llu order=%d "
            "coefficient_samples=%llu final_samples=%llu blocks=%zu\n",
            (unsigned long long)result.version, result.power_step,
            (unsigned long long)result.coefficient_samples,
            (unsigned long long)result.final_samples,
            result.block_count);
    fprintf(FileLS,
            "# energy energy_se variance variance_se energy_imag "
            "variance_imag tau_int effective_samples retained_rank "
            "condition gevp_residual\n");
    fprintf(FileLS,
            "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e "
            "%d %.18e %.18e\n",
            result.energy, result.energy_standard_error, result.variance,
            result.variance_standard_error,
            result.final_energy_imaginary, result.variance_imaginary,
            result.energy_tau_int, result.effective_sample_count,
            result.retained_rank, result.condition_estimate,
            result.gevp_residual);
    fprintf(stdout,
            "Stabilized power-Lanczos: E=%.17g +/- %.6g, "
            "variance=%.17g +/- %.6g\n",
            result.energy, result.energy_standard_error, result.variance,
            result.variance_standard_error);
  } else if (rank == 0) {
    fprintf(stderr, "Error: stabilized power-Lanczos failed: %s.\n",
            mvmc_krylov_status_string(status));
  }
  return status == MVMC_KRYLOV_STATUS_OK ? 0 : 1;
}

/*main program*/
int main(int argc, char* argv[])
{
  /* input file name */
  char fileDefList[256];
  char fileInitPara[256];

  int flagReadInitPara=0;
  int info=0;

  /* for MultiDef mode (-m option) */
  int flagMultiDef=0;
  int nMultiDef = 1;
  /* for Standard mode (-s option)*/
  int flagStandard = 0;
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
  while((option=getopt(argc,argv,"bhm:oF:esv"))!=-1) {
    switch(option) {
    case 'b': /* BinaryMode */
      FlagBinary=1;
      break;

    case 'h': /* Print Help Message*/
      printUsageError();
      printOption();
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

    case 'e': /* Expert mode (For compatibility)*/
      /*Nothing to do*/
      flagMultiDef = 0;
      break;

    case 's': /* Standard mode */
      flagMultiDef = 0;
      flagStandard = 1;
      break;

    case 'v': /* Print version */
      printVersion();
      MPI_Finalize();
      return 0;

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
    strcpy(fileDefList, argv[optind]);
    if(argc-optind>1) {
      flagReadInitPara = 1;
      strcpy(fileInitPara, argv[optind + 1]);
    }
  } else if(flagMultiDef==1) { /* MultiDef mode */
    strcpy(fileDefList, argv[optind+1]);
    if(argc-optind>2) {
      flagReadInitPara = 1;
      strcpy(fileInitPara, argv[optind + 2]);
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
  /*
   Standard mode: generating input files
  */
  if (flagStandard == 1) {
    if (rank0 == 0) {
      StdFace_main(fileDefList);
    }
    strcpy(fileDefList, "namelist.def");
  }
  MPI_Barrier(comm0);

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

  if(NProjBF > 0 && AllComplexFlag == 0 && NQPFull > 1) {
    /* The real BackFlow MultiQP path is kept serial for OpenMP stability. */
    omp_set_num_threads(1);
    NThread = 1;
  }

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
  LapackLWork = getLWork_fcmp_dim(NsizeMax); //TBC

  StartTimer(13);
  /* initialize variational parameters */
  if(rank0==0) fprintf(stdout,"Start: Initialize parameters.\n");
  InitParameter(); /* Run parallelly for synchronization of random generator */
  if(flagReadInitPara>0 && rank0==0) ReadInitParameter(fileInitPara);
  //[s] add read parameters respectively
  if(rank0==0){
    if(!ReadInputParameters(fileDefList, comm0)==0){
      //[ToDo]: Add Error procedure
      info=1;
    }
  }
  // free fileDefList
  if(rank0==0) free(cFileNameListFile);
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
    info = VMCPhysCal(comm0, comm1, comm2);
    if(rank0==0) fprintf(stdout,"End  : Calculate VMC physical quantities.\n");
    StopTimer(2);
  } else {
    info=1;
    if(rank0==0) fprintf(stderr,"error: NVMCCalMode must be 0 or 1.\n");
  }

  StopTimer(0);
  ReduceBFProfileCounter(comm0);
  if(rank0==0) {
    if(NVMCCalMode==0) {
      OutputTimerParaOpt();
    } else if(NVMCCalMode==1) {
      OutputTimerPhysCal();
    }
  }

  /* close output files */
  if(rank0==0) CloseFile(rank0);

  if(rank0==0) fprintf(stdout,"Start: Free Memory.\n");
  FreeMemory();
  FreeMemoryDef();
  if(rank0==0) fprintf(stdout,"End: Free Memory.\n");

#ifdef _mpi_use
  MPI_Comm_free(&comm2);
  MPI_Comm_free(&comm1);
#endif
  MPI_Comm_free(&comm0);
  MPI_Finalize();
  if(rank0==0) fprintf(stdout,"Finish calculation.\n");

  return info;
}

/*-- VMC Parameter Optimization --*/
int VMCParaOpt(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2) {
  int step;
  int info;
  int rank;
  int tmp_i;//DEBUG
  int iprogress;
  MPI_Comm_rank(comm_parent, &rank);

  for(step=0;step<NSROptItrStep;step++) {
    //printf("0 DUBUG make:step=%d TwoSz=%d\n",step,TwoSz);
    if(rank==0){
      OutputTime(step);
      if(NSROptItrStep<20){
        iprogress = (int) (100.0*step/NSROptItrStep);
        printf("Progress of Optimization: %d %%.\n", iprogress);
      }
      else if(step%(NSROptItrStep/20)==0){
        iprogress = (int) (100.0*step/NSROptItrStep);
        printf("Progress of Optimization: %d %%.\n", iprogress);
      }
    }

    StartTimer(20);
    //printf("1 DUBUG make:step=%d \n",step);
    if(iFlgOrbitalGeneral==0){//sz is conserved
      UpdateSlaterElm_fcmp();
    }else{
      UpdateSlaterElm_fsz();
    }
    //printf("2 DUBUG make:step=%d \n",step);
    UpdateQPWeight();
    StopTimer(20);
    StartTimer(3);
#ifdef _DEBUG_DETAIL
    printf("Debug: step %d, MakeSample.\n", step);
#endif
    if(FlagGrandCanonical != 0) {
      VMCMakeSampleGC(comm_child1);
    } else if(NProjBF ==0) {
      if(AllComplexFlag==0){ // real
        StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElm_real[tmp_i]= creal(SlaterElm[tmp_i]);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM_real[tmp_i]= creal(InvM[tmp_i]);
        StopTimer(69);
        if(iFlgOrbitalGeneral==0){
          VMCMakeSample_real(comm_child1);
        }else{
          VMCMakeSample_fsz_real(comm_child1);
        }
        StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM[tmp_i]      = InvM_real[tmp_i]+0.0*I;
        StopTimer(69);
      }else{// complex
        if(iFlgOrbitalGeneral==0){// sz =0 & complex
          VMCMakeSample(comm_child1);//VMCMakeSample(comm_child1);
        }else{
          VMCMakeSample_fsz(comm_child1);//VMCMakeSample(comm_child1);
        }
      }
    } else if(iFlgOrbitalGeneral==0) {
      if(AllComplexFlag==0){
        StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElm_real[tmp_i]= creal(SlaterElm[tmp_i]);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM_real[tmp_i]= creal(InvM[tmp_i]);
        StopTimer(69);
        VMC_BF_MakeSample_real(comm_child1);
        StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM[tmp_i]      = InvM_real[tmp_i]+0.0*I;
        StopTimer(69);
      }else{
        VMC_BF_MakeSample(comm_child1);
      }
    } else {
      VMC_BF_MakeSample_fsz(comm_child1);
    }
    StopTimer(3);
    StartTimer(4);
#ifdef _DEBUG_DETAIL
    printf("Debug: step %d, MainCal.\n", step);
#endif
    if(FlagGrandCanonical != 0) {
      VMCMainCalGC(comm_parent, comm_child1);
    } else if(NProjBF ==0) {
      if(iFlgOrbitalGeneral==0){//sz is conserved
        VMCMainCal(comm_parent, comm_child1);
      }else{//fsz
        VMCMainCal_fsz(comm_parent, comm_child1);
      }
    }else if(iFlgOrbitalGeneral==0){
      VMC_BF_MainCal(comm_parent, comm_child1);
    }else{
      VMC_BF_MainCal_fsz(comm_parent, comm_child1);
    }
    StopTimer(4);
    StartTimer(21);
#ifdef _DEBUG_DETAIL
    printf("Debug: step %d, AverageWE.\n", step);
#endif
    WeightAverageWE(comm_parent);
    StartTimer(25);//DEBUG
#ifdef _DEBUG_DETAIL
    printf("Debug: step %d, SROpt.\n", step);
#endif
    //if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){ //real & sz =0
    if(AllComplexFlag==0){ //real
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

#ifdef _DEBUG_DUMP_SROPTO_STORE
    if(rank==0){
      if(AllComplexFlag==0){ //real & sz=0
      //if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){ //real & sz=0
        for(i=0;i<SROptSize*NVMCSample;i++){
          fprintf(stderr, "DEBUG: SROptO_Store_real[%d]=%lf +I*%lf\n",i,creal(SROptO_Store_real[i]),cimag(SROptO_Store_real[i]));
        }
      }else{
        for(i=0;i<2*SROptSize*NVMCSample;i++){
          fprintf(stderr, "DEBUG: SROptO_Store[%d]=%lf +I*%lf\n",i,creal(SROptO_Store[i]),cimag(SROptO_Store[i]));
        }
      }
    }
#endif

#ifdef _DEBUG_DUMP_SROPTOO
    if(rank==0){
      if(AllComplexFlag==0){ //real
      //if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){ //real & sz=0
        for(i=0;i<(NSRCG==0 ? SROptSize*SROptSize: SROptSize*2);i++){
          fprintf(stderr, "DEBUG: SROptOO_real[%d]=%lf +I*%lf\n",i,creal(SROptOO_real[i]),cimag(SROptOO_real[i]));
        }
        for(i=0;i<SROptSize;i++){
          fprintf(stderr, "DEBUG: SROptHO_real[%d]=%lf +I*%lf\n",i,creal(SROptHO_real[i]),cimag(SROptHO_real[i]));
        }
        for(i=0;i<SROptSize;i++){
          fprintf(stderr, "DEBUG: SROptO_real[%d]=%lf +I*%lf\n",i,creal(SROptO_real[i]),cimag(SROptO_real[i]));
        }
      }else{
        for(i=0;i<(NSRCG==0 ? 2*SROptSize*(2*SROptSize): 2*SROptSize*2);i++){
          fprintf(stderr, "DEBUG: SROptOO[%d]=%lf +I*%lf\n",i,creal(SROptOO[i]),cimag(SROptOO[i]));
        }
        for(i=0;i<2*SROptSize;i++){
          fprintf(stderr, "DEBUG: SROptHO[%d]=%lf +I*%lf\n",i,creal(SROptHO[i]),cimag(SROptHO[i]));
        }
        for(i=0;i<2*SROptSize;i++){
          fprintf(stderr, "DEBUG: SROptO[%d]=%lf +I*%lf\n",i,creal(SROptO[i]),cimag(SROptO[i]));
        }
      }
    }
#endif

    StartTimer(5);
    if(NSRCG!=0){
      info = StochasticOptCG(comm_parent);
    }else{
      info = StochasticOpt(comm_parent);
    }
    //info = StochasticOptDiag(comm_parent);
    StopTimer(5);

#ifdef _DEBUG_DUMP_PARA
    for(int i=0; i<NPara; ++i){
      fprintf(stderr, "DEBUG: Para[%d] = %lf %lf\n", i, creal(Para[i]), cimag(Para[i]));
    }
#endif

    // DEBUG
    // abort();

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
  if(rank==0) {
    fprintf(stdout, "Start: Output opt params.\n");
    OutputOptData();
    fprintf(stdout, "End: Output opt params.\n");
  }

  return 0;
}

/*-- VMC Physical Quantity Calculation --*/
int VMCPhysCal(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2) {
  int ismp, tmp_i;
  int correctedInfo;
  int rank;
  MPI_Comm_rank(comm_parent, &rank);

  if(rank==0) fprintf(stdout, "Start: UpdateSlaterElm.\n");
  StartTimer(20);
  if(iFlgOrbitalGeneral==0){//sz is conserved
    UpdateSlaterElm_fcmp();
  }else{
    UpdateSlaterElm_fsz();
  }
  StopTimer(20);
  if(rank==0) fprintf(stdout, "End  : UpdateSlaterElm.\n");

  if (NLanczosMode == 1 && NLanczosEstimatorMode == 1 &&
      AllComplexFlag == 0) {
#pragma omp parallel for default(shared) private(tmp_i)
    for (tmp_i = 0;
         tmp_i < NQPFull * (2 * Nsite) * (2 * Nsite); ++tmp_i) {
      SlaterElm_real[tmp_i] = creal(SlaterElm[tmp_i]);
    }
  }

  if(rank==0) fprintf(stdout, "Start: Sampling.\n");
  for(ismp=0;ismp<NDataQtySmp;ismp++) {
    if(rank==0) OutputTime(ismp);
    FlushFile(0,rank);
    InitFilePhysCal(ismp, rank);
    if (NLanczosMode == 1 && NLanczosEstimatorMode == 1) {
      correctedInfo = RunPowerLanczosStabilized(
          ismp + NDataIdxStart, comm_parent);
      CloseFilePhysCal(rank);
      if (correctedInfo != 0) return correctedInfo;
      continue;
    }
    StartTimer(3);
    if(FlagGrandCanonical != 0) {
      VMCMakeSampleGC(comm_child1);
    } else if(NProjBF ==0) {
      //if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){//real & sz=0
      if(AllComplexFlag==0){//real
        // only for real TBC
        StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElm_real[tmp_i]= creal(SlaterElm[tmp_i]);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM_real[tmp_i]= creal(InvM[tmp_i]);
        StopTimer(69);
        // SlaterElm_real will be used in CalculateMAll, note that SlaterElm will not change before SR
        if(iFlgOrbitalGeneral==0){
          VMCMakeSample_real(comm_child1);
        }else{
          VMCMakeSample_fsz_real(comm_child1);
        }
        // only for real TBC
        StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM[tmp_i]      = InvM_real[tmp_i]+0.0*I;
        StopTimer(69);
        // only for real TBC
      }else{
        if(iFlgOrbitalGeneral==0){
          VMCMakeSample(comm_child1);
        }else{
          VMCMakeSample_fsz(comm_child1);
        }
      }
    }else if(iFlgOrbitalGeneral==0){ // NProjBF != 0
      if(AllComplexFlag==0){
        StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElm_real[tmp_i]= creal(SlaterElm[tmp_i]);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM_real[tmp_i]= creal(InvM[tmp_i]);
        StopTimer(69);
        VMC_BF_MakeSample_real(comm_child1);
        StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
        for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)     InvM[tmp_i]      = InvM_real[tmp_i]+0.0*I;
        StopTimer(69);
      }else{
        VMC_BF_MakeSample(comm_child1);
      }
    }else{
      VMC_BF_MakeSample_fsz(comm_child1);
    }

    StopTimer(3);
    StartTimer(4);
    if(rank==0) fprintf(stdout, "End  : Sampling.\n");

    if(rank==0) fprintf(stdout, "Start: Main calculation.\n");
    if(FlagGrandCanonical != 0) {
      VMCMainCalGC(comm_parent, comm_child1);
    } else if(NProjBF ==0) {
      if(iFlgOrbitalGeneral==0){
        VMCMainCal(comm_parent, comm_child1);
      }else{
        VMCMainCal_fsz(comm_parent, comm_child1);
      }
    }else if(iFlgOrbitalGeneral==0){
      VMC_BF_MainCal(comm_parent, comm_child1);
    }else{
      VMC_BF_MainCal_fsz(comm_parent, comm_child1);
    }
    if(rank==0) fprintf(stdout, "End  : Main calculation.\n");
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
  int i, k, offset;
  int lanczosStatus;

  /* zvo_out.dat */
//[s] MERGE BY TM
 // fprintf(FileOut, "% .18e % .18e % .18e \n", Etot, Etot2, (Etot2 - Etot*Etot)/(Etot*Etot));
 //   fprintf(FileOut, "% .18e % .18e  % .18e % .18e \n", creal(Etot),cimag(Etot), creal(Etot2), creal((Etot2 - Etot*Etot)/(Etot*Etot)));
   fprintf(FileOut, "% .18e % .18e  % .18e % .18e %.18e %.18e\n", creal(Etot),cimag(Etot), creal(Etot2), creal((Etot2 - Etot*Etot)/(Etot*Etot)),creal(Sztot),creal(Sztot2));
   if(FlagGrandCanonical != 0 && FileGC != NULL) {
     fprintf(FileGC, "% .18e % .18e % .18e\n",
             creal(Ntot), creal(Ntot2), creal(Ntot2 - Ntot*Ntot));
     fflush(FileGC);
   }
  // fprintf(FileOut, "% .18e % .18e % .18e \n", Etot, Etot2, (Etot2 - Etot*Etot)/(Etot*Etot));
 // fprintf(FileOut, "% .18e % .18e  % .18e % .18e \n", creal(Etot), cimag(Etot), creal(Etot2),
 //         creal((Etot2 - Etot * Etot) / (Etot * Etot)));
//[e] MERGE BY TM

  /* zvo_var.dat */
  if (FlagBinary == 0) { /* formatted output*/
    fprintf(FileVar, "% .18e % .18e 0.0 % .18e % .18e 0.0 ", creal(Etot), cimag(Etot), creal(Etot2), cimag(Etot2));
    for (i = 0; i < NPara; i++) fprintf(FileVar, "% .18e % .18e 0.0 ", creal(Para[i]), cimag(Para[i]));
    fprintf(FileVar, "\n");
    //for(i=0;i<NPara;i++)  printf("DEBUG:i=%d: % .18e % .18e  \n",i, creal(Para[i]),cimag(Para[i]));
  } else { /* binary output */
    fwrite(Para, sizeof(double), NPara, FileVar);
  }

  if (NVMCCalMode == 1) {
    /* zvo_cisajs.dat */
    if (NCisAjs > 0) {
      for (i = 0; i < NCisAjs; i++)
        fprintf(FileCisAjs, "%d %d %d %d % .18e  % .18e \n", CisAjsIdx[i][0], CisAjsIdx[i][1], CisAjsIdx[i][2],
                CisAjsIdx[i][3], creal(PhysCisAjs[i]), cimag(PhysCisAjs[i]));
      fprintf(FileCisAjs, "\n");
    }
    /* zvo_cisajscktaltex.dat */
    if (NCisAjsCktAlt > 0) {
      for (i = 0; i < NCisAjsCktAlt; i++)
        fprintf(FileCisAjsCktAlt, "% .18e  % .18e ", creal(PhysCisAjsCktAlt[i]), cimag(PhysCisAjsCktAlt[i]));
      fprintf(FileCisAjsCktAlt, "\n");
    }

    /* zvo_cisajscktalt.dat */
    if (NCisAjsCktAltDC > 0) {
      for (i = 0; i < NCisAjsCktAltDC; i++) {
        fprintf(FileCisAjsCktAltDC, "%d %d %d %d %d %d %d %d % .18e % .18e\n",
                CisAjsCktAltDCIdx[i][0], CisAjsCktAltDCIdx[i][1], CisAjsCktAltDCIdx[i][2], CisAjsCktAltDCIdx[i][3],
                CisAjsCktAltDCIdx[i][4], CisAjsCktAltDCIdx[i][5], CisAjsCktAltDCIdx[i][6], CisAjsCktAltDCIdx[i][7],
                creal(PhysCisAjsCktAltDC[i]), cimag(PhysCisAjsCktAltDC[i]));
      }
      fprintf(FileCisAjsCktAltDC, "\n");
    }

    /* zvo_twist.dat */
    if (NTwist > 0) {
      for (i = 0; i < NTwist; i++){
        fprintf(FileTwist, "%.18e  %.18e  %.18e  %.18e \n", creal(PhysTwist[i]), cimag(PhysTwist[i]),cabs(PhysTwist[i]),carg(PhysTwist[i]));
        //fprintf(FileTwist, "% .18e  % .18e \n", creal(PhysTwist[i]), cimag(PhysTwist[i]));
      }
    }

    /* zvo_NBodyG.dat */
    if (NNBodyG > 0) {
      for (i = 0; i < NNBodyG; i++) {
        fprintf(FileNBodyG, "%d", NBodyGN[i]);
        offset = NBodyGOffset[i];
        for (k = 0; k < NBodyGN[i]; k++) {
          fprintf(FileNBodyG, " %d %d %d %d",
                  NBodyGIdx[offset+k][0], NBodyGIdx[offset+k][1],
                  NBodyGIdx[offset+k][2], NBodyGIdx[offset+k][3]);
        }
        fprintf(FileNBodyG, " % .18e  % .18e\n",
                creal(PhysNBodyG[i]), cimag(PhysNBodyG[i]));
      }
      fprintf(FileNBodyG, "\n");
    }

    if (NLanczosMode > 0 && NLanczosEstimatorMode == 0) {
      if (NLanczosStep == 1) {
      if (AllComplexFlag == 0) { //real
        lanczosStatus = PhysCalLanczos_real(
          QQQQ_real, QCisAjsQ_real, QCisAjsCktAltQ_real, QCisAjsCktAltQDC_real,
          NLSHam, Nsite, NCisAjs, NCisAjs, iOneBodyGIdx, CisAjsIdx, NCisAjsCktAlt, NCisAjsCktAltDC, CisAjsCktAltDCIdx,
          NLanczosMode, FileLS, FileLSQQQQ, FileLSQCisAjsQ, FileLSQCisAjsCktAltQ,
          FileLSCisAjs, FileLSCisAjsCktAlt, FileLSCisAjsCktAltDC);
      }else { //complex
        lanczosStatus = PhysCalLanczos_fcmp(
          QQQQ, QCisAjsQ, QCisAjsCktAltQ, QCisAjsCktAltQDC,
          NLSHam, Nsite, NCisAjs, NCisAjs, iOneBodyGIdx, CisAjsIdx, NCisAjsCktAlt, NCisAjsCktAltDC, CisAjsCktAltDCIdx,
          NLanczosMode, FileLS, FileLSQQQQ, FileLSQCisAjsQ, FileLSQCisAjsCktAltQ,
          FileLSCisAjs, FileLSCisAjsCktAlt, FileLSCisAjsCktAltDC);
      }
      if (lanczosStatus != 0) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
      } else {
        Lanczos2Result lanczos2Result = {0};
        Lanczos2SolveStatus lanczos2Status;
        if (AllComplexFlag == 0) {
          lanczos2Status = WriteLanczos2OutputReal(
              FileLS2, FileLS2Moment, LS2Moment_real,
              LS2MomentBasisShift, &lanczos2Result);
        } else {
          lanczos2Status = WriteLanczos2OutputComplex(
              FileLS2, FileLS2Moment, LS2Moment,
              LS2MomentBasisShift, &lanczos2Result);
        }
        if (lanczos2Status != LANCZOS2_SOLVE_OK) {
          fprintf(stderr, "Error: Lanczos2 output failed: %s (LAPACK info=%d).\n",
                  Lanczos2SolveError(lanczos2Status),
                  lanczos2Result.lapackInfo);
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
      }
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
  fprintf(stderr,"  -s     Standard mode\n");
  fprintf(stderr,"  -e     Expert mode\n");
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
