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
 * initialization of files
 *-------------------------------------------------------------
 * by Satoshi Morita 
 *-------------------------------------------------------------*/

void InitFile(char *xNameListFile, int rank);
void InitFilePhysCal(int i, int rank);
void CloseFile(int rank);
void CloseFilePhysCal(int rank);
void FlushFile(int step, int rank);
void writeConfig(char *xNameFile, char *fileName);
int fileCopyAdd(char *inputFileName, FILE *outputFile);

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
      fprintf(FileSRinfo,
            "#Npara Msize optCut diagCut sDiagMax  sDiagMin    absRmax       imax\n");
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
  
  if(NLanczosMode>0){
    sprintf(fileName, "%s_ls_%03d.dat", CDataFileHead, idx);
    FileLS = fopen(fileName, "w");

    sprintf(fileName, "%s_ls_qqqq_%03d.dat", CDataFileHead, idx);
    FileLSQQQQ = fopen(fileName, "w");
    
    if(NLanczosMode>1){
      sprintf(fileName, "%s_ls_qcisajsq_%03d.dat",
              CDataFileHead, idx);
      FileLSQCisAjsQ = fopen(fileName, "w");

      sprintf(fileName, "%s_ls_qcisajscktaltq_%03d.dat", 
              CDataFileHead, idx);
      FileLSQCisAjsCktAltQ = fopen(fileName, "w");
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
  
  if(NLanczosMode>0){
    fclose(FileLS);
    fclose(FileLSQQQQ);
    
    if(NLanczosMode>1){
      fclose(FileLSQCisAjsQ);
      fclose(FileLSQCisAjsCktAltQ);
    }
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
