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
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * calculate physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita 
 *-------------------------------------------------------------*/

void VMCMainCal_fsz(MPI_Comm comm);

void VMCMainCal_fsz(MPI_Comm comm) {
  int *eleIdx,*eleCfg,*eleNum,*eleProjCnt,*eleSpn; //fsz
  double complex e,ip;
  double w;
  double sqrtw;
  double complex we;
  double Sz;

  const int qpStart=0;
  const int qpEnd=NQPFull;
  int sample,sampleStart,sampleEnd,sampleSize;
  int i,info,tmp_i;

  /* optimazation for Kei */
  const int nProj=NProj;
  double complex *srOptO = SROptO;
  double         *srOptO_real = SROptO_real;

  int rank,size,int_i;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
#ifdef _DEBUG
  printf("  Debug: SplitLoop\n");
#endif
  SplitLoop(&sampleStart,&sampleEnd,NVMCSample,rank,size);

  /* initialization */
  StartTimer(24);
  clearPhysQuantity();
  StopTimer(24);
  for(sample=sampleStart;sample<sampleEnd;sample++) {
    eleIdx = EleIdx + sample*Nsize;
    eleCfg = EleCfg + sample*Nsite2;
    eleNum = EleNum + sample*Nsite2;
    eleProjCnt = EleProjCnt + sample*NProj;
    eleSpn     = EleSpn + sample*Nsize; //fsz

    StartTimer(40);
#ifdef _DEBUG
    printf("  Debug: sample=%d: CalculateMAll \n",sample);
#endif
    info = CalculateMAll_fsz(eleIdx,eleSpn,qpStart,qpEnd);//info = CalculateMAll_fcmp(eleIdx,qpStart,qpEnd); // InvM,PfM will change
    StopTimer(40);

    if(info!=0) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n",rank,sample,info);
      continue;
    }
#ifdef _DEBUG
    printf("  Debug: sample=%d: CalculateIP \n",sample);
#endif
    ip = CalculateIP_fcmp(PfM,qpStart,qpEnd,MPI_COMM_SELF);
#ifdef _DEBUG
    printf("  Debug: sample=%d: LogProjVal \n",sample);
#endif
    //LogProjVal(eleProjCnt);
    /* calculate reweight */
    w =1.0;
#ifdef _DEBUG
    printf("  Debug: sample=%d: isfinite \n",sample);
#endif
    if( !isfinite(w) ) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d w=%e\n",rank,sample,w);
      continue;
    }

    StartTimer(41);
    /* calculate energy */
#ifdef _DEBUG
    printf("  Debug: sample=%d: calculateHam \n",sample);
#endif
#ifdef _DEBUG
    printf("  Debug: sample=%d: calculateHam_cmp \n",sample);
#endif
    e  = CalculateHamiltonian_fsz(ip,eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn);//fsz
    Sz = CalculateSz_fsz(ip,eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn);//fsz
    StopTimer(41);
    if( !isfinite(creal(e) + cimag(e)) ) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d e=%e\n",rank,sample,creal(e)); //TBC
      continue;
    }

    Wc    += w;
    Etot  += w * e;
    Sztot += w * Sz;
    Etot2 += w * conj(e) * e;
#ifdef _DEBUG
    printf("  Debug: sample=%d: calculateOpt \n",sample);
#endif
    if(NVMCCalMode==0) {
      /* Calculate O for correlation fauctors */
      srOptO[0] = 1.0+0.0*I;//   real 
      srOptO[1] = 0.0+0.0*I;//   real 
      #pragma loop noalias
      for(i=0;i<nProj;i++){ 
        srOptO[(i+1)*2]     = (double)(eleProjCnt[i]); // even real
        srOptO[(i+1)*2+1]   = 0.0+0.0*I;               // odd  comp
      }

      StartTimer(42);
      /* SlaterElmDiff */
      SlaterElmDiff_fsz(SROptO+2*NProj+2,ip,eleIdx,eleSpn) ;//SlaterElmDiff_fcmp(SROptO+2*NProj+2,ip,eleIdx); //TBC: using InvM not InvM_real
      StopTimer(42);
      
      if(FlagOptTrans>0) { // this part will be not used
        calculateOptTransDiff(SROptO+2*NProj+2*NSlater+2, ip); //TBC
      }
      StartTimer(43);
      /* Calculate OO and HO */
      if(NStoreO==0){
        calculateOO(SROptOO,SROptHO,SROptO,w,e,SROptSize);
      }else{
        we    = w*e;
        sqrtw = sqrt(w); 
        #pragma omp parallel for default(shared) private(int_i)
        for(int_i=0;int_i<SROptSize*2;int_i++){
        // SROptO_Store for fortran
          SROptO_Store[int_i+sample*(2*SROptSize)]  = sqrtw*SROptO[int_i];
          SROptHO[int_i]                           += we*SROptO[int_i]; 
        }
      } 
      StopTimer(43);
    } else if(NVMCCalMode==1) {
      StartTimer(42);
      /* Calculate Green Function */
      CalculateGreenFunc_fsz(w,ip,eleIdx,eleCfg,eleNum,eleSpn,eleProjCnt);
      StopTimer(42);

      if(NLanczosMode>0){
        // for sz!=0, Lanczso is not supported
      }
    }
  } /* end of for(sample) */

// calculate OO and HO at NVMCCalMode==0
  if(NStoreO!=0 && NVMCCalMode==0){
    sampleSize=sampleEnd-sampleStart;
    StartTimer(45);
    calculateOO_Store(SROptOO,SROptHO,SROptO_Store,w,e,2*SROptSize,sampleSize);
    StopTimer(45);
  }
  return;
}
