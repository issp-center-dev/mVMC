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
 * initialization of files
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "initfile.h"

#ifndef _INITFILE_SRC
#define _INITFILE_SRC

void InitFile(char *xNameListFile, int rank) {
  char fileName[D_FileNameMax];

  if(rank!=0) return;

  //sprintf(fileName, "%s_cfg_%03d.dat", CDataFileHead, NDataIdxStart);
  //writeConfig(xNameListFile, fileName);

  sprintf(fileName, "%s_time_%03d.dat", CDataFileHead, NDataIdxStart);
  FileTime = fopen(fileName, "w");

  if(NVMCCalMode==0) {
    sprintf(fileName, "%s_SRinfo.dat", CDataFileHead);
    FileSRinfo = fopen(fileName, "w");
    if(SRFlag == 0){
      if (NSRCG != 0) {
        fprintf(FileSRinfo,
              "#Npara Msize optCut diagCut sDiagMax  sDiagMin    absRmax       imax info\n");
      } else {
        fprintf(FileSRinfo,
              "#Npara Msize optCut diagCut sDiagMax  sDiagMin    absRmax       imax\n");
      }
    }else{
      fprintf(FileSRinfo,
            "#Npara Msize optCut diagCut sEigenMax  sEigenMin    absRmax       imax\n");
    }

    sprintf(fileName, "%s_out_%03d.dat", CDataFileHead, NDataIdxStart);
    FileOut = fopen(fileName, "w");

    if(FlagBinary==0) {
      sprintf(fileName, "%s_var_%03d.dat", CDataFileHead, NDataIdxStart);
      FileVar = fopen(fileName, "w");
    } else {
      sprintf(fileName, "%s_varbin_%03d.dat", CDataFileHead, NDataIdxStart);
      FileVar = fopen(fileName, "wb");
      fwrite(&NPara,sizeof(int),1,FileVar);
      fwrite(&NSROptItrStep,sizeof(int),1,FileVar);
    }
  }

  return;
}

void InitFilePhysCal(int i, int rank) {
  char fileName[D_FileNameMax];
  int idx = i+NDataIdxStart;
  int one = 1;

  if(rank!=0) return;

  sprintf(fileName, "%s_out_%03d.dat", CDataFileHead, idx);
  FileOut = fopen(fileName, "w");

  if(FlagBinary==0) {
    sprintf(fileName, "%s_var_%03d.dat", CDataFileHead, idx);
    FileVar = fopen(fileName, "w");
  } else {
    sprintf(fileName, "%s_varbin_%03d.dat", CDataFileHead, idx);
    FileVar = fopen(fileName, "wb");
    fwrite(&NPara,sizeof(int),1,FileVar);
    fwrite(&one,sizeof(int),1,FileVar);
  }

  /* Green function */
  if(NCisAjs>0){
    sprintf(fileName, "%s_cisajs_%03d.dat", CDataFileHead, idx);
    FileCisAjs = fopen(fileName, "w");
  }

  if(NCisAjsCktAlt>0){
    sprintf(fileName, "%s_cisajscktaltex_%03d.dat", CDataFileHead, idx);
    FileCisAjsCktAlt = fopen(fileName, "w");
  }

  if(NCisAjsCktAltDC>0){
    sprintf(fileName, "%s_cisajscktalt_%03d.dat", CDataFileHead, idx);
    FileCisAjsCktAltDC = fopen(fileName, "w");
  }

  if(NTwist>0){
    sprintf(fileName, "%s_twist_%03d.dat", CDataFileHead, idx);
    FileTwist = fopen(fileName, "w");
  }

  if(NNBodyG>0){
    sprintf(fileName, "%s_NBodyG_%03d.dat", CDataFileHead, idx);
    FileNBodyG = fopen(fileName, "w");
  }

  if(NLanczosMode>0 && NLanczosEstimatorMode==1 && NLanczosStep==1){
    sprintf(fileName, "%s_ls_out_%03d.dat", CDataFileHead, idx);
    FileLS = fopen(fileName, "w");

    sprintf(fileName, "%s_ls_qqqq_%03d.dat", CDataFileHead, idx);
    FileLSQQQQ = fopen(fileName, "w");

    if(NLanczosMode>1){
#ifdef _DEBUG
      sprintf(fileName, "%s_ls_qcisajsq_%03d.dat",
              CDataFileHead, idx);
      FileLSQCisAjsQ = fopen(fileName, "w");
      sprintf(fileName, "%s_ls_qcisajscktaltq_%03d.dat",
              CDataFileHead, idx);
      FileLSQCisAjsCktAltQ = fopen(fileName, "w");

 #endif
      sprintf(fileName, "%s_ls_cisajs_%03d.dat",
              CDataFileHead, idx);
      FileLSCisAjs = fopen(fileName, "w");

      /* CACA */
      sprintf(fileName, "%s_ls_cisajscktaltex_%03d.dat",
              CDataFileHead, idx);
      FileLSCisAjsCktAlt = fopen(fileName, "w");

      /* CACADC */
      sprintf(fileName, "%s_ls_cisajscktalt_%03d.dat",
              CDataFileHead, idx);
      FileLSCisAjsCktAltDC = fopen(fileName, "w");
    }
  }
  if(NLanczosMode>0 && NLanczosEstimatorMode==1 && NLanczosStep==2){
    sprintf(fileName, "%s_ls2_out_%03d.dat", CDataFileHead, idx);
    FileLS2 = fopen(fileName, "w");
    if(FileLS2 == NULL){
      fprintf(stderr, "Error: failed to open Lanczos2 output file '%s'.\n",
              fileName);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    sprintf(fileName, "%s_ls2_moment_%03d.dat", CDataFileHead, idx);
    FileLS2Moment = fopen(fileName, "w");
    if(FileLS2Moment == NULL){
      fprintf(stderr, "Error: failed to open Lanczos2 moment file '%s'.\n",
              fileName);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
  if(NLanczosMode>0 && NLanczosEstimatorMode==1){
    sprintf(fileName, "%s_ls_support_%03d.dat", CDataFileHead, idx);
    FileLSSupport = fopen(fileName, "w");
    if(FileLSSupport == NULL){
      fprintf(stderr,
              "Error: failed to open power-Lanczos support file '%s'.\n",
              fileName);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  return;
}

void CloseFile(int rank) {
  if(rank!=0) return;

  fclose(FileTime);

  if(NVMCCalMode==0) {
    fclose(FileSRinfo);
    fclose(FileOut);
    fclose(FileVar);
  }

  return;
}

void CloseFilePhysCal(int rank) {
  if(rank!=0) return;

  fclose(FileOut);
  fclose(FileVar);

  if(NCisAjs>0){
    fclose(FileCisAjs);
  }
  if(NCisAjsCktAlt>0){
    fclose(FileCisAjsCktAlt);
  }
  if(NCisAjsCktAltDC>0){
    fclose(FileCisAjsCktAltDC);
  }
  if(NTwist>0){
    fclose(FileTwist);
  }
  if(NNBodyG>0){
    fclose(FileNBodyG);
  }

  if(NLanczosMode>0 && NLanczosEstimatorMode==1 && NLanczosStep==1){
    fclose(FileLS);
    fclose(FileLSQQQQ);

    if(NLanczosMode>1){
#ifdef _DEBUG
      fclose(FileLSQCisAjsQ);
      fclose(FileLSQCisAjsCktAltQ);
#endif
      fclose(FileLSCisAjs);
      fclose(FileLSCisAjsCktAlt);
      fclose(FileLSCisAjsCktAltDC);
    }
  }
  if(NLanczosMode>0 && NLanczosEstimatorMode==1 && NLanczosStep==2){
    fclose(FileLS2);
    fclose(FileLS2Moment);
  }
  if(NLanczosMode>0 && NLanczosEstimatorMode==1){
    fclose(FileLSSupport);
  }

  return;
}

void FlushFile(int step, int rank) {
  if(rank!=0) return;

  if(step%NFileFlushInterval==0) {
    fflush(FileTime);
    if(NVMCCalMode==0) {
      fflush(FileSRinfo);
      fflush(FileOut);
      fflush(FileVar);
    }
  }
  return;
}

void writeConfig(char *xnamefile, char *fileName) {
  FILE *ofp,*fplist;
  char defname[D_FileNameMax];

  /* clear configfile */
  ofp = fopen(fileName, "w");

  fileCopyAdd(xnamefile,ofp);

  fplist = fopen(xnamefile,"r");
  if(fplist==NULL){
    fprintf(stderr, "Error writecfg.c");
  } else {
    while(fscanf(fplist, "%s\n", defname)!=EOF){
      fileCopyAdd(defname,ofp);
    }
    fclose(fplist);
  }

  fclose(ofp);
  return;
}

int fileCopyAdd(char *inputfileName, FILE *ofp){
  int i;
  FILE *ifp;

  ifp = fopen(inputfileName, "r");
  if(ifp == NULL) return 1;

  fprintf(ofp,"################\n#%s\n################\n",inputfileName);
  while((i=getc(ifp))!=EOF) putc(i,ofp);
  fclose(ifp);

  return 0;
}

#endif
