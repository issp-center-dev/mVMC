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
 * Read Definition Files
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#include <ctype.h>
#include <stdlib.h>
#include "./include/readdef.h"
#include "./include/global.h"
#include "safempi_fcmp.c"


#define _NOTBACKFLOW

int ReadDefFileError(const char *defname);

int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm);

int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm);


int CheckSite(const int iSite, const int iMaxNum);

int CheckPairSite(const int iSite1, const int iSite2, const int iMaxNum);

int CheckQuadSite(const int iSite1, const int iSite2, const int iSite3, const int iSite4, const int iMaxNum);

int GetTransferInfo(FILE *fp, int **ArrayIdx, double complex *ArrayValue, int Nsite, int NArray, char *defname);

int GetLocSpinInfo(FILE *fp, int *ArrayIdx, int Nsite, char *defname);

int GetInfoCoulombIntra(FILE *fp, int *ArrayIdx, double *ArrayValue, int Nsite, int NArray, char *defname);

int ReadPairHopValue(FILE *fp, int **ArrayIdx, double *ArrayValue, int Nsite, int NArray, char *defname);

int ReadPairDValue(FILE *fp, int **ArrayIdx, double *ArrayValue, int Nsite, int NArray, char *defname);

int GetInfoGutzwiller(FILE *fp, int *ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int Nsite, int NArray,
                      char *defname);

int GetInfoJastrow(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int fidx, int Nsite,
                   int NArray, char *defname);

int
GetInfoDH2(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nsite, int NArray,
           char *defname);

int
GetInfoDH4(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nsite, int NArray,
           char *defname);

int GetInfoTransSym(FILE *fp, int **Array, int **ArraySgn, int **ArrayInv, double complex *ArrayPara, int _APFlag,
                    int Nsite, int NArray, char *defname);

int GetInfoOrbitalParallel(FILE *fp, int **Array, int *ArrayOpt, int **ArraySgn, int *iOptCount,
                           int _fidx, int _iComplexFlag, int _iFlagOrbitalGeneral, int _APFlag, int Nsite, int NArray,
                           int NArrayAP, char *defname);

int GetInfoOrbitalAntiParallel(FILE *fp, int **Array, int *ArrayOpt, int **ArraySgn, int *iOptCount,
                               int _fidx, int _iComplexFlag, int _iFlagOrbitalGeneral, int _APFlag, int Nsite,
                               int NArray, char *defname);

int GetInfoInterAll(FILE *fp, int **ArrayIdx, double complex *ArrayValue,
                    int Nsite, int NArray, char *defname);

int GetInfoOptTrans(FILE *fp, int **Array, double *ArrayPara, int *ArrayOpt, int **ArraySgn,
                    int _iFlagOptTrans, int *iOptCount, int _fidx, int _APFlag, int Nsite, int NArray, char *defname);

int GetInfoTwoBodyG(FILE *fp, int **ArrayIdx, int **ArrayIdxTwoBodyGLz, int **ArrayToIdx, int **ArrayIdxOneBodyG,
                    int _NLanczosMode, int Nsite, int NArray, char *defname);

int GetInfoTwoBodyGEx(FILE *fp, int **ArrayIdx, int Nsite, int NArray, char *defname);

int GetInfoOrbitalGeneral(FILE *fp, int **Array, int *ArrayOpt, int **ArraySgn, int *iOptCount,
                          int _fidx, int _iComplexFlag, int _iFlagOrbitalGeneral, int _APFlag, int Nsite, int NArray,
                          char *defname);

int
GetInfoOneBodyG(FILE *fp, int **ArrayIdx, int **ArrayToIdx, int _NLanczosMode, int Nsite, int NArray, char *defname);

char *ReadBuffInt(FILE *fp, int *iNbuf) {
  char *cerr;
  char ctmp[D_FileNameMax];
  char ctmp2[D_FileNameMax];
  cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
  if (cerr != NULL) {
    cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
    sscanf(ctmp2, "%s %d\n", ctmp, iNbuf); //2
  }
  return cerr;
}

char *ReadBuffIntCmpFlg(FILE *fp, int *iNbuf, int *iComplexFlag) {
  char *cerr;
  char ctmp[D_FileNameMax];
  char ctmp2[D_FileNameMax];
  cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
  if (cerr != NULL) {
    cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
    sscanf(ctmp2, "%s %d\n", ctmp, iNbuf); //2
    if (*iNbuf == 0) {
      cerr = NULL;
      fprintf(stderr, "error: Number of components defined in the 2nd line is set as 0 or illegal header.\n");
      fprintf(stderr, "error: 2nd line %s", ctmp2);
    } else {
      cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
      sscanf(ctmp2, "%s %d\n", ctmp, iComplexFlag);
    }
  }
  return cerr;
}

void SetDefaultValuesModPara(int *buf, double *bufDouble);

int GetInfoFromModPara(int *buf, double *bufDouble);


int ReadDefFileError(const char *defname) {
  fprintf(stderr, "error: %s (Broken file or Not exist)\n", defname);
  return 1;
}

///
/// \param _iFlgOrbitalGeneral Flag of Orbital General
/// \param _iFlgOrbitalAP Flag of Orbital General
/// \param _iFlgOrbitalP
/// \retval 0  normal
/// \retval -1 multiple definition
/// \retval -2 too few definition

int JudgeOrbitalMode(int *_iFlgOrbitalGeneral, const int _iFlgOrbitalAP, const int _iFlgOrbitalP) {

  int iret = 0;
  //(General, AP, P)
  if (*_iFlgOrbitalGeneral == 1) {
    if (_iFlgOrbitalAP == 0 && _iFlgOrbitalP == 0) { //(1, 0, 0)
      iret = 0;
    } else {//(1, 1, 0) or (1, 0, 1) or (1, 1, 1)
      iret = -1;
    }
  } else {
    if (_iFlgOrbitalAP == 1) {
      if (_iFlgOrbitalP == 1) { //(0, 1, 1)
        *_iFlgOrbitalGeneral = 1;
        iret = 0;
      } else { //(0, 1, 0)
        iret = 0;
      }
    } else {// (0, 0, 0) or (0, 0, 1)
      iret = -2;
    }
  }
  if (iret == -1) {
    fprintf(stderr, "error: Multiple definition of Orbital files.\n");
  } else if (iret == -2) {
    fprintf(stderr, "error: Not exist any Orbital file or Need OrbitalAP file.\n");
  }
  return iret;
}

int ReadGreen(char *xNameListFile, int Nca, int **caIdx, int Ncacadc, int **cacaDCIdx, int Ns) {
  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  int iKWidx = 0;
  char *cerr;
  int i, info = 0;

  cFileNameListFile = malloc(sizeof(char) * D_CharTmpReadDef * KWIdxInt_end);
  fprintf(stdout, "  Read File %s .\n", xNameListFile);
  if (GetFileName(xNameListFile, cFileNameListFile) != 0) {
    fprintf(stderr, "  error: Definition files(*.def) are incomplete.\n");
    return -1;
  }

  for (iKWidx = KWLocSpin; iKWidx < KWIdxInt_end; iKWidx++) {
    strcpy(defname, cFileNameListFile[iKWidx]);

    if (strcmp(defname, "") == 0) continue;

    fp = fopen(defname, "r");
    if (fp == NULL) {
      info = ReadDefFileError(defname);
      fclose(fp);
      continue;
    }

    /*=======================================================================*/
    for (i = 0; i < IgnoreLinesInDef; i++)
      cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      if(cerr == NULL) return(-1);
    switch (iKWidx) {
      case KWOneBodyG:
        /*cisajs.def----------------------------------------*/
        if (GetInfoOneBodyG(fp, caIdx, iOneBodyGIdx, 0, Ns, Nca, defname) != 0) return (-1);
        break;
      case KWTwoBodyG:
        /*cisajscktaltdc.def--------------------------------*/
        if (GetInfoTwoBodyG(fp, cacaDCIdx, cacaDCIdx, iOneBodyGIdx, caIdx, 0, Ns, Ncacadc, defname) != 0) return (-1);
        break;
      default:
        break;
    }
    fclose(fp);
  }
  return info;
}

///
/// \param xNameListFile FileNameLists
/// \param Nca Number of CisAjs
/// \param Ncacadc Number of CisAjsCktAltDC
/// \param Ns Number of sites
/// \return Number of calculation target
int CountOneBodyGForLanczos(char *xNameListFile, int Nca, int Ncacadc, int Ns, int **caIdx, int **iFlgOneBodyG) {

  int info = 0;
  int i, j, isite1, isite2;
  int icount = 0;
  int **cacaDCIdx;

  cacaDCIdx = malloc(sizeof(int *) * Ncacadc);
  //pInt=cacaDCIdx[0];
  for (i = 0; i < Ncacadc; i++) {
    cacaDCIdx[i] = malloc(sizeof(int) * 8);
  }

  for (i = 0; i < 2 * Ns; i++) {
    for (j = 0; j < 2 * Ns; j++) {
      iFlgOneBodyG[i][j] = -1;
    }
  }
  info = ReadGreen(xNameListFile, Nca, caIdx, Ncacadc, cacaDCIdx, Ns);
  if (info != 0) {
    free(cacaDCIdx);
    return (info);
  }

  for (i = 0; i < Nca; i++) {
    isite1 = caIdx[i][0] + caIdx[i][1] * Ns;
    isite2 = caIdx[i][2] + caIdx[i][3] * Ns;
    if (iFlgOneBodyG[isite1][isite2] == -1) {
      iFlgOneBodyG[isite1][isite2] = icount;
      icount++;
    }
  }
  //cisajscktalt -> cisajs, cltakt (Note: indecies of the latter Green's function are modified)
  for (i = 0; i < Ncacadc; i++) {
    isite1 = cacaDCIdx[i][0] + cacaDCIdx[i][1] * Ns;
    isite2 = cacaDCIdx[i][2] + cacaDCIdx[i][3] * Ns;
    if (iFlgOneBodyG[isite1][isite2] == -1) {
      iFlgOneBodyG[isite1][isite2] = icount;
      icount++;
    }

    /*
    isite1 = cacaDCIdx[i][4] + cacaDCIdx[i][5] * Ns;
    isite2 = cacaDCIdx[i][6] + cacaDCIdx[i][7] * Ns;
    */
    isite1 = cacaDCIdx[i][6] + cacaDCIdx[i][7] * Ns;
    isite2 = cacaDCIdx[i][4] + cacaDCIdx[i][5] * Ns;
    if (iFlgOneBodyG[isite1][isite2] == -1) {
      iFlgOneBodyG[isite1][isite2] = icount;
      icount++;
    }
    
  }
  free(cacaDCIdx);
  return icount;
}

int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm) {
  FILE *fp;
  char defname[D_FileNameMax];
  char *cerr;

  int rank, info = 0;
  const int nBufInt = ParamIdxInt_End;
  const int nBufDouble = ParamIdxDouble_End;
  const int nBufChar = D_FileNameMax;
  int bufInt[nBufInt];
  double bufDouble[nBufDouble];
  int iKWidx = 0;
  int iret = 0;
  int iFlgOrbitalAntiParallel = 0;
  int iFlgOrbitalParallel = 0;

  int iOrbitalComplex = 0;
  iFlgOrbitalGeneral = 0;
  MPI_Comm_rank(comm, &rank);

  if (rank == 0) {
    cFileNameListFile = malloc(sizeof(char) * D_CharTmpReadDef * KWIdxInt_end);
    fprintf(stdout, "  Read File %s .\n", xNameListFile);
    if (GetFileName(xNameListFile, cFileNameListFile) != 0) {
      fprintf(stderr, "  error: Definition files(*.def) are incomplete.\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (iKWidx = 0; iKWidx < KWIdxInt_end; iKWidx++) {
      strcpy(defname, cFileNameListFile[iKWidx]);
      if (strcmp(defname, "") == 0) {
        switch (iKWidx) {
          case KWModPara:
          case KWLocSpin:
            fprintf(stderr, "  Error: Need to make a def file for %s.\n", cKWListOfFileNameList[iKWidx]);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            break;
          default:
            break;
        }
      }
    }

    SetDefaultValuesModPara(bufInt, bufDouble);
    iret = GetInfoFromModPara(bufInt, bufDouble);
    if (iret != 0) {
      if (rank == 0) {
        fprintf(stderr, "  Error: ModPara file is incomplete.\n");
      }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }


    for (iKWidx = 0; iKWidx < KWIdxInt_end; iKWidx++) {
      strcpy(defname, cFileNameListFile[iKWidx]);

      if (strcmp(defname, "") == 0) continue;
      fprintf(stdout, "  Read File '%s' for %s.\n", defname, cKWListOfFileNameList[iKWidx]);
      fp = fopen(defname, "r");
      if (fp == NULL) {
        info = ReadDefFileError(defname);
        fclose(fp);
        break;
      } else {
        switch (iKWidx) {
          case KWModPara:
            cerr = "";
            break;

          case KWLocSpin:
            cerr = ReadBuffInt(fp, &bufInt[IdxNLocSpin]);
            break;

          case KWTrans:
            cerr = ReadBuffInt(fp, &bufInt[IdxNTrans]);
            break;

          case KWCoulombIntra:
            cerr = ReadBuffInt(fp, &bufInt[IdxNCoulombIntra]);
            break;

          case KWCoulombInter:
            cerr = ReadBuffInt(fp, &bufInt[IdxNCoulombInter]);
            break;

          case KWHund:
            cerr = ReadBuffInt(fp, &bufInt[IdxNHund]);
            break;

          case KWPairHop:
            cerr = ReadBuffInt(fp, &bufInt[IdxNPairHop]);
            break;

          case KWExchange:
            cerr = ReadBuffInt(fp, &bufInt[IdxNExchange]);
            break;

          case KWGutzwiller:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNGutz], &iComplexFlgGutzwiller);
            break;

          case KWJastrow:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNJast], &iComplexFlgJastrow);
            break;

          case KWDH2:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNDH2], &iComplexFlgDH2);
            break;

          case KWDH4:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNDH4], &iComplexFlgDH4);
            break;

          case KWOrbital:
          case KWOrbitalAntiParallel:
            cerr = ReadBuffIntCmpFlg(fp, &iNOrbitalAntiParallel, &iOrbitalComplex);
            iFlgOrbitalAntiParallel = 1;
            bufInt[IdxNOrbit] += iNOrbitalAntiParallel;
            iComplexFlgOrbital += iOrbitalComplex;
            break;

          case KWOrbitalParallel:
            cerr = ReadBuffIntCmpFlg(fp, &iNOrbitalParallel, &iOrbitalComplex);
            iFlgOrbitalParallel = 1;
            bufInt[IdxNOrbit] += 2 * iNOrbitalParallel; //up-up and down-down
            iComplexFlgOrbital += iOrbitalComplex;
            break;

          case KWOrbitalGeneral:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNOrbit], &iOrbitalComplex);
            iFlgOrbitalGeneral = 1;
            iComplexFlgOrbital += iOrbitalComplex;
            break;

          case KWTransSym:
            cerr = ReadBuffInt(fp, &bufInt[IdxNQPTrans]);
            break;

          case KWOneBodyG:
            cerr = ReadBuffInt(fp, &bufInt[IdxNOneBodyG]);
            break;

          case KWTwoBodyG:
            cerr = ReadBuffInt(fp, &bufInt[IdxNTwoBodyG]);
            break;

          case KWTwoBodyGEx:
            cerr = ReadBuffInt(fp, &bufInt[IdxNTwoBodyGEx]);
            break;

          case KWInterAll:
            cerr = ReadBuffInt(fp, &bufInt[IdxNInterAll]);
            break;

          case KWOptTrans:
            bufInt[IdxNQPOptTrans] = 1;
            if (FlagOptTrans > 0) {
              cerr = ReadBuffInt(fp, &bufInt[IdxNQPOptTrans]);
              if (bufInt[IdxNQPOptTrans] < 1) {
                fprintf(stderr, "Error: NQPOptTrans should be larger than 0.\n");
                info = ReadDefFileError(defname);
              }
            }
            break;

          case KWBFRange:
#ifdef _NOTBACKFLOW
            fprintf(stderr, "Error: Back Flow is not supported.\n");
            info = ReadDefFileError(defname);
#else
          fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
          fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
          sscanf(ctmp2,"%s %d %d\n", ctmp, &bufInt[IdxNrange], &bufInt[IdxNNz]);
#endif
            break;

          case KWBF:
#ifdef _NOTBACKFLOW
            fprintf(stderr, "Error: Back Flow is not supported.\n");
            info = ReadDefFileError(defname);
#else
            cerr = ReadBuffInt(fp, &bufInt[IdxNBF]);
#endif
            break;

          default:
            cerr = "";
            break;
        }//case KW
        if (cerr == NULL) {
          info = ReadDefFileError(defname);
        }
        fclose(fp);
      }
    }

    iret = JudgeOrbitalMode(&iFlgOrbitalGeneral, iFlgOrbitalAntiParallel, iFlgOrbitalParallel);
    if (iret < 0) info = iret;

    //TODO: LanczosMode is not supported for Sz not conserved mode.

    //For Lanczos mode: Calculation of Green's function
    if (bufInt[IdxLanczosMode] > 1) {
      //Get info of CisAjs and CisAjsCktAltDC
      int i;
      NCisAjsLz = bufInt[IdxNOneBodyG];
      //bufInt[IdxNTwoBodyGEx] = bufInt[IdxNTwoBodyG];
      CisAjsLzIdx = malloc(sizeof(int *) * NCisAjsLz);
      for (i = 0; i < NCisAjsLz; i++) {
        CisAjsLzIdx[i] = malloc(sizeof(int) * 4);
      }
      iOneBodyGIdx = malloc(sizeof(int *) * (2 * bufInt[IdxNsite])); //For spin
      //pInt=iFlgOneBodyG[0];
      for (i = 0; i < 2 * bufInt[IdxNsite]; i++) {
        iOneBodyGIdx[i] = malloc(sizeof(int) * (2 * bufInt[IdxNsite]));
      }
      bufInt[IdxNOneBodyG] = CountOneBodyGForLanczos(xNameListFile, NCisAjsLz, bufInt[IdxNTwoBodyG],
                                                     bufInt[IdxNsite], CisAjsLzIdx, iOneBodyGIdx);
    }

    //CalcNCond
    if (bufInt[IdxNCond] != -1) {
      if (bufInt[IdxNCond] % 2 != 0) {
        fprintf(stderr, "Error: NCond (in modpara.def) must be even number.\n");
        info = 1;
      } else bufInt[IdxNe] = (bufInt[IdxNLocSpin] + bufInt[IdxNCond]) / 2;
    }

    //CheckGeneral Orbital
    //printf("bufInt[Idx2Sz]=%d \n",bufInt[Idx2Sz]);
    if (bufInt[Idx2Sz] != 0) {
      //if(iOrbitalComplex != 1){
      if (iFlgOrbitalGeneral != 1) {
        fprintf(stderr,
                "Error: OrbitalParallel or OrbitalGeneral files must be needed when 2Sz !=0 (in modpara.def).\n");
        info = 1;
      } else if (bufInt[Idx2Sz] % 2 != 0 && bufInt[Idx2Sz] != -1) {
        fprintf(stderr, "Error: 2Sz (in modpara.def) must be even number.\n");
        info = 1;
      }
    }

    if (iFlgOrbitalGeneral == 1) {
      if (bufInt[IdxSPGaussLeg] > 1) {    //Check NSPGaussLeg
        fprintf(stdout, "Warning: SPGaussLeg (in modpara.def) must be 0 or 1 when orbital is general.\n");
        fprintf(stdout, "         SPGaussLeg set as 1.\n");
        bufInt[IdxSPGaussLeg] = 1;
      } else if (bufInt[IdxLanczosMode] != 0) {
        fprintf(stderr, "Error: Lanczos mode is not supported when orbital is general.\n");
        info = 1;
      }
    }

    //Check LocSpn
    if (bufInt[IdxNLocSpin] > 0) {
      if (bufInt[IdxNLocSpin] == 2 * bufInt[IdxNe] && bufInt[IdxExUpdatePath] != 2) {
        fprintf(stderr, "Error: NExUpdatePath (in modpara.def) must be 2 when 2*Ne = NLocalSpin, i.e. spin system.\n");
        info = 1;
      } else if (bufInt[IdxExUpdatePath] == 0) {
        fprintf(stderr, "Error: NExUpdatePath (in modpara.def) must be 1.\n");
        info = 1;
      } else if (bufInt[IdxNLocSpin] > 2 * bufInt[IdxNe]) {
        fprintf(stderr, "Error: 2*Ne must satisfy the condition, 2*Ne >= NLocalSpin.\n");
        info = 1;
      }
      if (info == 1) {
        fprintf(stderr, "  Error: ModPara file is incomplete.\n");
      }
    }

  }

  if (info != 0) {
    if (rank == 0) {
      fprintf(stderr, "Error: Definition files(*.def) are incomplete.\n");
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if (rank == 0) {
    AllComplexFlag = iComplexFlgGutzwiller + iComplexFlgJastrow + iComplexFlgDH2; //TBC
    AllComplexFlag += iComplexFlgDH4 + iComplexFlgOrbital;//TBC
    //AllComplexFlag  = 1;//DEBUG
    // AllComplexFlag= 0 -> All real, !=0 -> complex
  }

#ifdef _mpi_use
  MPI_Bcast(bufInt, nBufInt, MPI_INT, 0, comm);
  MPI_Bcast(&NStoreO, 1, MPI_INT, 0, comm); // for NStoreO
  MPI_Bcast(&NSRCG, 1, MPI_INT, 0, comm); // for NCG
  MPI_Bcast(&AllComplexFlag, 1, MPI_INT, 0, comm); // for Real
  MPI_Bcast(&iFlgOrbitalGeneral, 1, MPI_INT, 0, comm); // for fsz
  MPI_Bcast(bufDouble, nBufDouble, MPI_DOUBLE, 0, comm);
  MPI_Bcast(CDataFileHead, nBufChar, MPI_CHAR, 0, comm);
  MPI_Bcast(CParaFileHead, nBufChar, MPI_CHAR, 0, comm);
#endif /* _mpi_use */

  NVMCCalMode = bufInt[IdxVMCCalcMode];
  NLanczosMode = bufInt[IdxLanczosMode];
  NDataIdxStart = bufInt[IdxDataIdxStart];
  NDataQtySmp = bufInt[IdxDataQtySmp];
  Nsite = bufInt[IdxNsite];
  Ne = bufInt[IdxNe];
  NSPGaussLeg = bufInt[IdxSPGaussLeg];
  NSPStot = bufInt[IdxSPStot];
  NMPTrans = bufInt[IdxMPTrans];
  NSROptItrStep = bufInt[IdxSROptItrStep];
  NSROptItrSmp = bufInt[IdxSROptItrSmp];
  NSROptFixSmp = bufInt[IdxSROptFixSmp];
  NVMCWarmUp = bufInt[IdxVMCWarmUp];
  NVMCInterval = bufInt[IdxVMCInterval];
  NVMCSample = bufInt[IdxVMCSample];
  NExUpdatePath = bufInt[IdxExUpdatePath];
  RndSeed = bufInt[IdxRndSeed];
  NSplitSize = bufInt[IdxSplitSize];
  NLocSpn = bufInt[IdxNLocSpin];
  NTransfer = bufInt[IdxNTrans];
  NCoulombIntra = bufInt[IdxNCoulombIntra];
  NCoulombInter = bufInt[IdxNCoulombInter];
  NHundCoupling = bufInt[IdxNHund];
  NPairHopping = 2 * bufInt[IdxNPairHop];
  NExchangeCoupling = bufInt[IdxNExchange];
  NGutzwillerIdx = bufInt[IdxNGutz];
  NJastrowIdx = bufInt[IdxNJast];
  NDoublonHolon2siteIdx = bufInt[IdxNDH2];
  NDoublonHolon4siteIdx = bufInt[IdxNDH4];
  NOrbitalIdx = bufInt[IdxNOrbit];
  NQPTrans = bufInt[IdxNQPTrans];
  NCisAjs = bufInt[IdxNOneBodyG];
  NCisAjsCktAlt = bufInt[IdxNTwoBodyGEx];
  NCisAjsCktAltDC = bufInt[IdxNTwoBodyG];
  NInterAll = bufInt[IdxNInterAll];
  NQPOptTrans = bufInt[IdxNQPOptTrans];
  Nrange = bufInt[IdxNrange];
  NBackFlowIdx = bufInt[IdxNBF];
  Nz = bufInt[IdxNNz];
  NSROptCGMaxIter = bufInt[IdxSROptCGMaxIter];
  DSROptRedCut = bufDouble[IdxSROptRedCut];
  DSROptStaDel = bufDouble[IdxSROptStaDel];
  DSROptStepDt = bufDouble[IdxSROptStepDt];
  DSROptCGTol = bufDouble[IdxSROptCGTol];
  TwoSz = bufInt[Idx2Sz];

  if (NMPTrans < 0) {
    APFlag = 1; /* anti-periodic boundary */
    NMPTrans *= -1;
  } else {
    APFlag = 0;
  }

  if (DSROptStepDt < 0) {
    SRFlag = 1; /* diagonalization */
    if (rank == 0) fprintf(stderr, "remark: Diagonalization Mode\n");
    DSROptStepDt *= -1;
  } else {
    SRFlag = 0;
  }

  Nsize = 2 * Ne;
  Nsite2 = 2 * Nsite;
  NSlater = NOrbitalIdx;
  NProj = NGutzwillerIdx + NJastrowIdx
          + 2 * 3 * NDoublonHolon2siteIdx
          + 2 * 5 * NDoublonHolon4siteIdx;
  NOptTrans = (FlagOptTrans > 0) ? NQPOptTrans : 0;

  /* [s] For BackFlow */
  if (NBackFlowIdx > 0) {
    NrangeIdx = 3 * (Nrange - 1) / Nz + 1; //For Nz-conectivity
    NBFIdxTotal = (NrangeIdx - 1) * (NrangeIdx) / 2 + (NrangeIdx);
    NProjBF = NBFIdxTotal * NBackFlowIdx;
  } else {
    NrangeIdx = 0;
    NBFIdxTotal = 0;
    NProjBF = 0;
  }
  /* [e] For BackFlow */

  NPara = NProj + NSlater + NOptTrans + NProjBF;
  NQPFix = NSPGaussLeg * NMPTrans;
  NQPFull = NQPFix * NQPOptTrans;
  SROptSize = NPara + 1;

  NTotalDefInt = Nsite /* LocSpn */
                 + 4 * NTransfer /* Transfer */
                 + NCoulombIntra /* CoulombIntra */
                 + 2 * NCoulombInter /* CoulombInter */
                 + 2 * NHundCoupling /* HundCoupling */
                 + 2 * NPairHopping /* PairHopping */
                 + 2 * NExchangeCoupling /* ExchangeCoupling */
                 + Nsite /* GutzwillerIdx */
                 + Nsite * Nsite /* JastrowIdx */
                 + 2 * Nsite * NDoublonHolon2siteIdx /* DoublonHolon2siteIdx */
                 + 4 * Nsite * NDoublonHolon4siteIdx /* DoublonHolon4siteIdx */
                 //+ (2*Nsite)*(2*Nsite) /* OrbitalIdx */ //fsz
                 //+ (2*Nsite)*(2*Nsite) /* OrbitalSgn */ //fsz
                 + Nsite * NQPTrans /* QPTrans */
                 + Nsite * NQPTrans /* QPTransInv */
                 + Nsite * NQPTrans /* QPTransSgn */
                 + 4 * NCisAjs /* CisAjs */
                 + 8 * NCisAjsCktAlt /* CisAjsCktAlt */
                 + 8 * NCisAjsCktAltDC /* CisAjsCktAltDC */
                 + 8 * NInterAll /* InterAll */
                 + Nsite * NQPOptTrans /* QPOptTrans */
                 + Nsite * NQPOptTrans /* QPOptTransSgn */
                 + 2 * NPara; /* OptFlag */ // TBC

  //Orbitalidx
  if (iFlgOrbitalGeneral == 0) {
    NTotalDefInt += Nsite * Nsite; // OrbitalIdx
    NTotalDefInt += Nsite * Nsite; // OrbitalSgn
  } else if (iFlgOrbitalGeneral == 1) {
    NTotalDefInt += (2 * Nsite) * (2 * Nsite); //OrbitalIdx
    NTotalDefInt += (2 * Nsite) * (2 * Nsite); //OrbitalSgn
  }

  if (NBackFlowIdx > 0) {
    NTotalDefInt += Nsite * Nsite * Nsite * Nsite; /* BackflowIdx */
  }

  NTotalDefDouble =
      NCoulombIntra /* ParaCoulombIntra */
      + NCoulombInter /* ParaCoulombInter */
      + NHundCoupling /* ParaHondCoupling */
      + NPairHopping  /* ParaPairHopping */
      + NExchangeCoupling /* ParaExchangeCoupling */
      //    + NQPTrans /* ParaQPTrans */
      //+ NInterAll /* ParaInterAll */
      + NQPOptTrans; /* ParaQPTransOpt */

  return 0;
}

int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm) {
  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  int iKWidx = 0;
  char *cerr;
  int i, j, n, idx0, idx1, info = 0;
  int fidx = 0; /* index for OptFlag */
  int count_idx = 0;
  int x0, x1;
  int rank;

  int iNOneBodyG;

  MPI_Comm_rank(comm, &rank);

  if (rank == 0) {
    for (iKWidx = KWLocSpin; iKWidx < KWIdxInt_end; iKWidx++) {
      strcpy(defname, cFileNameListFile[iKWidx]);
      if (strcmp(defname, "") == 0) continue;
      fprintf(stdout, "     %s\n", defname);
      fp = fopen(defname, "r");
      if (fp == NULL) {
        info = ReadDefFileError(defname);
        fclose(fp);
        continue;
      }

      /*=======================================================================*/
      for (i = 0; i < IgnoreLinesInDef; i++) fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      idx0 = 0;
      switch (iKWidx) {
        case KWLocSpin: /* Read locspn.def----------------------------------------*/
          if (GetLocSpinInfo(fp, LocSpn, Nsite, defname) != 0) info = 1;
          break;//locspn

        case KWTrans: /* transfer.def--------------------------------------*/
          if (GetTransferInfo(fp, Transfer, ParaTransfer, Nsite, NTransfer, defname) != 0) info = 1;
          break;

        case KWCoulombIntra: /*coulombintra.def----------------------------------*/
          if (GetInfoCoulombIntra(fp, CoulombIntra, ParaCoulombIntra, Nsite, NCoulombIntra, defname) != 0) info = 1;
          break;

        case KWCoulombInter: /*coulombinter.def----------------------------------*/
          if (ReadPairDValue(fp, CoulombInter, ParaCoulombInter, Nsite, NCoulombInter, defname) != 0) info = 1;
          break;

        case KWHund: /*hund.def------------------------------------------*/
          if (ReadPairDValue(fp, HundCoupling, ParaHundCoupling, Nsite, NHundCoupling, defname) != 0) info = 1;
          break;

        case KWExchange: /*exchange.def--------------------------------------*/
          if (ReadPairDValue(fp, ExchangeCoupling, ParaExchangeCoupling, Nsite, NExchangeCoupling, defname) != 0)
            info = 1;
          break;

        case KWPairHop: /*pairhop.def---------------------------------------*/
          if (ReadPairHopValue(fp, PairHopping, ParaPairHopping, Nsite, NPairHopping, defname) != 0) info = 1;
          break;

        case KWGutzwiller: /*gutzwilleridx.def---------------------------------*/
          if (GetInfoGutzwiller(fp, GutzwillerIdx, OptFlag, iComplexFlgGutzwiller, &count_idx, Nsite, NGutzwillerIdx,
                                defname) != 0)
            info = 1;
          break;

        case KWJastrow:
          fidx = NGutzwillerIdx;
          if (GetInfoJastrow(fp, JastrowIdx, OptFlag, iComplexFlgJastrow, &count_idx, fidx, Nsite, NJastrowIdx,
                             defname) != 0)
            info = 1;
          break;

        case KWDH2:
          /*doublonholon2siteidx.def--------------------------*/
          fidx = NGutzwillerIdx + NJastrowIdx;
          if (GetInfoDH2(fp, DoublonHolon2siteIdx, OptFlag, iComplexFlgDH2, &count_idx, fidx, Nsite,
                         NDoublonHolon2siteIdx, defname) != 0)
            info = 1;
          break;

        case KWDH4:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          if (GetInfoDH4(fp, DoublonHolon4siteIdx, OptFlag, iComplexFlgDH4, &count_idx, fidx, Nsite,
                         NDoublonHolon4siteIdx, defname) != 0)
            info = 1;
          break;

        case KWOrbital:
        case KWOrbitalAntiParallel:
          /*orbitalidxs.def------------------------------------*/
          fidx = NProj;
          if (GetInfoOrbitalAntiParallel(fp, OrbitalIdx, OptFlag, OrbitalSgn, &count_idx,
                                         fidx, iComplexFlgOrbital, iFlgOrbitalGeneral, APFlag, Nsite, iNOrbitalAntiParallel,
                                         defname) != 0)
            info = 1;
          break;

        case KWOrbitalGeneral:
          fidx = NProj;
          if (GetInfoOrbitalGeneral(fp, OrbitalIdx, OptFlag, OrbitalSgn, &count_idx,
                                    fidx, iComplexFlgOrbital, iFlgOrbitalGeneral, APFlag, Nsite, NOrbitalIdx,
                                    defname) != 0)
            info = 1;
          break;

        case KWOrbitalParallel:
          /*orbitalidxt.def------------------------------------*/
          fidx = NProj + iNOrbitalAntiParallel;
          if (GetInfoOrbitalParallel(fp, OrbitalIdx, OptFlag, OrbitalSgn, &count_idx,
                                     fidx, iComplexFlgOrbital, iFlgOrbitalGeneral, APFlag, Nsite, iNOrbitalParallel,
                                     iNOrbitalAntiParallel, defname) != 0)
            info = 1;
          break;

        case KWTransSym:
          /*qptransidx.def------------------------------------*/
          if (GetInfoTransSym(fp, QPTrans, QPTransSgn, QPTransInv, ParaQPTrans, APFlag, Nsite, NQPTrans, defname) !=
              0)
            info = 1;
          break;

        case KWOneBodyG:
          /*cisajs.def----------------------------------------*/
          iNOneBodyG = (NLanczosMode < 2) ? NCisAjs : NCisAjsLz;
          if (GetInfoOneBodyG(fp, CisAjsIdx, iOneBodyGIdx, NLanczosMode, Nsite, iNOneBodyG, defname) != 0) info = 1;
          break;

        case KWTwoBodyGEx:
          /*cisajscktalt.def----------------------------------*/
          if (GetInfoTwoBodyGEx(fp, CisAjsCktAltIdx, Nsite, NCisAjsCktAlt, defname) != 0) info = 1;
          break;

        case KWTwoBodyG:
          /*cisajscktaltdc.def--------------------------------*/
          if (GetInfoTwoBodyG(fp, CisAjsCktAltDCIdx, CisAjsCktAltLzIdx, iOneBodyGIdx, CisAjsIdx, NLanczosMode, Nsite,
                              NCisAjsCktAltDC, defname) != 0)
            info = 1;
          break;

        case KWInterAll:
          /*interall.def---------------------------------------*/
          if (GetInfoInterAll(fp, InterAll, ParaInterAll, Nsite, NInterAll, defname) != 0) info = 1;
          break;

        case KWOptTrans:
          /*qpopttrans.def------------------------------------*/
          fidx = NProj + NOrbitalIdx;
          if (GetInfoOptTrans(fp, QPOptTrans, ParaQPOptTrans, OptFlag, QPOptTransSgn, FlagOptTrans, &count_idx, fidx,
                              APFlag, Nsite, NQPOptTrans, defname) != 0)
            info = 1;
          break;

        case KWBFRange:
          /*rangebf.def--------------------------*/
          if (Nrange > 0) {
            for (i = 0; i < 5; i++){
              cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
              if(cerr == NULL) info=1;
            }
            idx0 = idx1 = 0;
            while (fscanf(fp, "%d %d %d\n", &i, &j, &n) != EOF) {
              PosBF[i][idx0 % Nrange] = j;
              RangeIdx[i][j] = n;
              idx0++;
              if (idx0 == Nsite * Nrange) break;
            }
            if (idx0 != Nsite * Nrange) {
              info = ReadDefFileError(defname);
            }
          }
          break;

        case KWBF:
          if (NBackFlowIdx > 0) {
            for (i = 0; i < 5; i++){
              cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
              if(cerr == NULL) info=1;
            }
            idx0 = idx1 = 0;
            for (i = 0; i < Nsite * Nsite; i++) {
              for (j = 0; j < Nsite * Nsite; j++) {
                BackFlowIdx[i][j] = -1;
              }
            }
            while (fscanf(fp, "%d %d %d %d %d\n", &i, &j, &(x0), &(x1), &n) != EOF) {
              BackFlowIdx[i * Nsite + j][x0 * Nsite + x1] = n;
              idx0++;
              //printf("idx0=%d, idx=%d\n",idx0,BackFlowIdx[i]);
              if (idx0 == Nsite * Nsite * Nrange * Nrange) break;
            }
            while (fscanf(fp, "%d ", &i) != EOF) {
              fscanf(fp, "%d\n", &(OptFlag[fidx]));
              //printf("idx1=%d, OptFlag=%d\n",idx1,OptFlag[fidx]);
              fidx++;
              idx1++;
            }
            if (idx0 != Nsite * Nsite * Nrange * Nrange
                || idx1 != NBFIdxTotal * NBackFlowIdx) {
              info = ReadDefFileError(defname);
            }
          }
          break;

        default:
          break;
      }
      fclose(fp);
    }

    if (count_idx != NPara) {
      fprintf(stdout, "count_idx=%d, NPara=%d\n", count_idx, NPara);
      fprintf(stderr, "error: OptFlag is incomplete.\n");
      info = 1;
    }
    fprintf(stdout, "finish reading parameters.\n");
  } /* if(rank==0) */

  if (FlagOptTrans <= 0) { // initialization of QPOptTrans
    ParaQPOptTrans[0] = 1.0;
    for (i = 0; i < Nsite; ++i) {
      QPOptTrans[0][i] = i;
      QPOptTransSgn[0][i] = 1;
    }
  }
  if (info != 0) {
    if (rank == 0) {
      fprintf(stderr, "error: Indices and Parameters of Definition files(*.def) are incomplete.\n");
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //Debug
  /*
  for(i =0; i<2*Nsite; i++){
    for(j=0; j<2*Nsite; j++){
      fprintf(stdout, "Debug: OrbitalIdx[%d][%d]=", i, j);
      fprintf(stdout,"Debug: %d\n", OrbitalIdx[i][j]);
    }
  }
  */
#ifdef _mpi_use
  SafeMpiBcastInt(LocSpn, NTotalDefInt, comm);
  if (NLanczosMode > 1) {
    SafeMpiBcastInt(CisAjsCktAltLzIdx[0], NCisAjsCktAltDC * 2, comm);
  }
  SafeMpiBcast_fcmp(ParaTransfer, NTransfer + NInterAll, comm);
  SafeMpiBcast(ParaCoulombIntra, NTotalDefDouble, comm);
  SafeMpiBcast_fcmp(ParaQPTrans, NQPTrans, comm);
#endif /* _mpi_use */

  /* set FlagShift */
  if (NVMCCalMode == 0) {
    SetFlagShift();
    if (rank == 0 && FlagShiftGJ + FlagShiftDH2 + FlagShiftDH4 > 0) {
      fprintf(stderr, "remark: FlagShift ( ");
      if (FlagShiftGJ == 1) fprintf(stderr, "GJ ");
      if (FlagShiftDH2 == 1) fprintf(stderr, "DH2 ");
      if (FlagShiftDH4 == 1) fprintf(stderr, "DH4 ");
      fprintf(stderr, ") is turned on.\n");
    }
  }

  return 0;
}

int ReadInputParameters(char *xNameListFile, MPI_Comm comm) {
  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax], ctmp2[256];
  int iKWidx = 0;
  int i, idx;
  int rank;
  int count = 0;
  int info = 0;
  int iNTotalIdx;
  double tmp_real, tmp_comp;
  MPI_Comm_rank(comm, &rank);
  char *cerr;
  if (rank == 0) {
    for (iKWidx = KWInGutzwiller; iKWidx < KWIdxInt_end; iKWidx++) {
      strcpy(defname, cFileNameListFile[iKWidx]);
      if (strcmp(defname, "") == 0) continue;
      fp = fopen(defname, "r");
      if (fp == NULL) {
        info = ReadDefFileError(defname);
        fclose(fp);
        continue;
      }
      /*=======================================================================*/
      cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);//1
      if(cerr == NULL) return -1;
      cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);//2
      if(cerr == NULL) return -1;
      sscanf(ctmp2, "%s %d\n", ctmp, &idx);
      cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);//3
      if(cerr == NULL) return -1;
      cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);//4
      if(cerr == NULL) return -1;
      cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);//5
      if(cerr == NULL) return -1;

      switch (iKWidx) {
        //get idx

        case KWInGutzwiller:
          fprintf(stdout, "Read InGutzwiller File. \n");
          if (idx != NGutzwillerIdx) {
            info = 1;
            continue;
          }
          count = 0;
          for (i = 0; i < NGutzwillerIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            Proj[idx+count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInJastrow:
          fprintf(stdout, "Read InJastrow File. \n");
          if (idx != NJastrowIdx) {
            info = 1;
            continue;
          }
          count = NGutzwillerIdx;
          for (i = count; i < count + NJastrowIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            Proj[idx+count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInDH2:
          fprintf(stdout, "Read InDH2 File. \n");
          if (idx != NDoublonHolon2siteIdx) {
            info = 1;
            continue;
          }
          count = NGutzwillerIdx + NJastrowIdx;
          for (i = count; i < count + 2 * 3 * NDoublonHolon2siteIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            Proj[idx+count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInDH4:
          fprintf(stdout, "Read InDH4 File. \n");
          if (idx != NDoublonHolon4siteIdx) {
            info = 1;
            continue;
          }
          count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          for (i = count; i < count + 2 * 5 * NDoublonHolon4siteIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            Proj[idx+count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInOrbital:
        case KWInOrbitalAntiParallel:
          iNTotalIdx = (iFlgOrbitalGeneral > 0) ? iNOrbitalAntiParallel : NOrbitalIdx;
          fprintf(stdout, "Read InOrbital/InOrbitalAntiParallel File. \n");
          if (idx != iNTotalIdx) {
            info = 1;
            continue;
          }
          for (i = 0; i < iNTotalIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            Slater[idx] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInOrbitalParallel:
          fprintf(stdout, "Read iNOrbitalParallel File. \n");
          //printf("MDEBUG: %d %d \n",idx,iNOrbitalParallel);
          if ((idx / 2) != iNOrbitalParallel) {
            info = 1;
            continue;
          }
          for (i = 0; i < 2 * iNOrbitalParallel; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            Slater[iNOrbitalAntiParallel + idx] = tmp_real + I * tmp_comp; //up-up
            //printf("MDEBUG: %d %d %d %lf\n",i,idx,iNOrbitalParallel,tmp_real,tmp_comp);
            //ierr = fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            //Slater[iNOrbitalAntiParallel+idx+1] = tmp_real + I * tmp_comp;//down-down
          }
          break;

        case KWInorbitalGeneral:
          fprintf(stdout, "Read InOrbitalGeneral File. \n");
          if (idx != NOrbitalIdx) {
            info = 1;
            continue;
          }
          for (i = 0; i < NOrbitalIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            Slater[idx] = tmp_real + I * tmp_comp; //up-up
          }
          break;

        case KWInOptTrans:
          fprintf(stdout, "Read InOptTrans File. \n");
          if (idx != NOptTrans) {
            info = 1;
            continue;
          }
          for (i = 0; i < NOptTrans; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            OptTrans[idx] = tmp_real + I * tmp_comp;
          }
          break;

        default:
          info = 0;
          break;
      }
      if (info != 0) {
        fprintf(stderr, "error: Indices of %s file is incomplete.\n", defname);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
      fclose(fp);
    }
  } /* if(rank==0) */

  if (info != 0) {
    if (rank == 0) {
      fprintf(stderr, "error: Indices of %s file is incomplete.\n", defname);
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  return info;
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
    const char *ctmp,
    const char *cKeyWord
) {

  int i = 0;

  char ctmp_small[256] = {0};
  char cKW_small[256] = {0};
  int n;
  n = strlen(cKeyWord);
  strncpy(cKW_small, cKeyWord, n);

  for (i = 0; i < n; i++) {
    cKW_small[i] = tolower(cKW_small[i]);
  }

  n = strlen(ctmp);
  strncpy(ctmp_small, ctmp, n);
  for (i = 0; i < n; i++) {
    ctmp_small[i] = tolower(ctmp_small[i]);
  }
  if (n < strlen(cKW_small)) n = strlen(cKW_small);
  return (strncmp(ctmp_small, cKW_small, n));
}

/**
 * @brief Check Site Number.
 * @param[in] *iSite a site number.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckSite(
    const int iSite,
    const int iMaxNum
) {
  if (iSite >= iMaxNum || iSite < 0) return (-1);
  return 0;
}

/**
 * @brief Check Site Number for a pair -> (siteA, siteB).
 * @param[in] iSite1 a site number on a site A.
 * @param[in] iSite2 a site number on a site B.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckPairSite(
    const int iSite1,
    const int iSite2,
    const int iMaxNum
) {
  if (CheckSite(iSite1, iMaxNum) != 0) {
    return (-1);
  }
  if (CheckSite(iSite2, iMaxNum) != 0) {
    return (-1);
  }
  return 0;
}

/**
 * @brief Check Site Number for a quad -> (siteA, siteB, siteC, siteD).
 * @param[in] iSite1 a site number on site A.
 * @param[in] iSite2 a site number on site B.
 * @param[in] iSite3 a site number on site C.
 * @param[in] iSite4 a site number on site D.
 * @param[in] iMaxNum Max site number.
 * @retval 0 normally finished reading file.
 * @retval -1 unnormally finished reading file.
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int CheckQuadSite(
    const int iSite1,
    const int iSite2,
    const int iSite3,
    const int iSite4,
    const int iMaxNum
) {
  if (CheckPairSite(iSite1, iSite2, iMaxNum) != 0) {
    return (-1);
  }
  if (CheckPairSite(iSite3, iSite4, iMaxNum) != 0) {
    return (-1);
  }
  return 0;
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
    const char *cKW,
    char cKWList[][D_CharTmpReadDef],
    int iSizeOfKWidx,
    int *iKWidx
) {
  *iKWidx = -1;
  int itmpKWidx;
  int iret = -1;
  for (itmpKWidx = 0; itmpKWidx < iSizeOfKWidx; itmpKWidx++) {
    if (strcmp(cKW, "") == 0) {
      break;
    } else if (CheckWords(cKW, cKWList[itmpKWidx]) == 0) {
      iret = 0;
      *iKWidx = itmpKWidx;
    }
  }
  return iret;
}

/**
 * @brief Function of Validating value.
 * @param[in] icheckValue value to validate.
 * @param[in] ilowestValue lowest value which icheckValue can be set.
 * @param[in] iHighestValue highest value which icheckValue can be set.
 * @retval 0 value is correct.
 * @retval -1 value is incorrect.
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 **/
int ValidateValue(
    const int icheckValue,
    const int ilowestValue,
    const int iHighestValue
) {

  if (icheckValue < ilowestValue || icheckValue > iHighestValue) {
    return (-1);
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
) {
  char *ctmpRead;
  char *cerror;
  char csplit[] = " ,.\t\n";
  if (*ctmpLine == '\n') return (-1);
  ctmpRead = strtok(ctmpLine, csplit);
  if (strncmp(ctmpRead, "=", 1) == 0 || strncmp(ctmpRead, "#", 1) == 0 || ctmpRead == NULL) {
    return (-1);
  }
  strcpy(ctmp, ctmpRead);

  ctmpRead = strtok(NULL, csplit);
  *itmp = strtol(ctmpRead, &cerror, 0);
  //if ctmpRead is not integer type
  if (*cerror != '\0') {
    fprintf(stderr, "Error: incorrect format= %s. \n", cerror);
    return (-1);
  }

  ctmpRead = strtok(NULL, csplit);
  if (ctmpRead != NULL) {
    fprintf(stderr, "Error: incorrect format= %s. \n", ctmpRead);
    return (-1);
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
    const char *cFileListNameFile,
    char cFileNameList[][D_CharTmpReadDef]
) {
  FILE *fplist;
  int itmpKWidx = -1;
  char ctmpFileName[D_FileNameMaxReadDef];
  char ctmpKW[D_CharTmpReadDef], ctmp2[256];
  int i;
  for (i = 0; i < KWIdxInt_end; i++) {
    strcpy(cFileNameList[i], "");
  }

  fplist = fopen(cFileListNameFile, "r");
  if (fplist == NULL) return ReadDefFileError(cFileListNameFile);

  while (fgets(ctmp2, 256, fplist) != NULL) {
    sscanf(ctmp2, "%s %s\n", ctmpKW, ctmpFileName);
    if (strncmp(ctmpKW, "#", 1) == 0 || *ctmp2 == '\n' || (strcmp(ctmpKW, "") && strcmp(ctmpFileName, "")) == 0) {
      continue;
    } else if (strcmp(ctmpKW, "") * strcmp(ctmpFileName, "") == 0) {
      fprintf(stderr, "Error: keyword and filename must be set as a pair in %s.\n", cFileListNameFile);
      fclose(fplist);
      return (-1);
    }
    /*!< Check KW */
    if (CheckKW(ctmpKW, cKWListOfFileNameList, KWIdxInt_end, &itmpKWidx) != 0) {

      fprintf(stderr, "Error: Wrong keywords '%s' in %s.\n", ctmpKW, cFileListNameFile);
      fprintf(stderr, "%s", "Choose Keywords as follows: \n");
      for (i = 0; i < KWIdxInt_end; i++) {
        fprintf(stderr, "%s \n", cKWListOfFileNameList[i]);
      }
      fclose(fplist);
      return (-1);
    }
    /*!< Check cFileNameList to prevent from double registering the file name */
    if (strcmp(cFileNameList[itmpKWidx], "") != 0) {
      fprintf(stderr, "Error: Same keywords exist in %s.\n", cFileListNameFile);
      fclose(fplist);
      return (-1);
    }
    /*!< Copy FileName */
    strcpy(cFileNameList[itmpKWidx], ctmpFileName);
  }
  fclose(fplist);
  return 0;
}


/**********************************************************************/
/* Function for ModPara file*/
/**********************************************************************/

/// @brief Function of setting default values
/// \param bufInt buffer for int values
/// \param bufDouble buffer for double values
void SetDefaultValuesModPara(int *bufInt, double *bufDouble) {
  bufInt[IdxVMCCalcMode] = 0;
  bufInt[IdxLanczosMode] = 0;
  bufInt[IdxDataIdxStart] = 0;
  bufInt[IdxDataQtySmp] = 1;
  bufInt[IdxNsite] = 16;
  bufInt[IdxNe] = 8;
  bufInt[IdxSPGaussLeg] = 1;
  bufInt[IdxSPStot] = 0;
  bufInt[IdxMPTrans] = 0;
  bufInt[IdxSROptItrStep] = 1000;
  bufInt[IdxSROptItrSmp] = bufInt[IdxSROptItrStep] / 10;
  bufInt[IdxSROptFixSmp] = 1;
  bufInt[IdxVMCWarmUp] = 10;
  bufInt[IdxVMCInterval] = 1;
  bufInt[IdxVMCSample] = 10;
  bufInt[IdxExUpdatePath] = 0;
  bufInt[IdxRndSeed] = 11272;
  bufInt[IdxSplitSize] = 1;
  bufInt[IdxNLocSpin] = 0;
  bufInt[IdxNTrans] = 0;
  bufInt[IdxNCoulombIntra] = 0;
  bufInt[IdxNCoulombInter] = 0;
  bufInt[IdxNHund] = 0;
  bufInt[IdxNPairHop] = 0;
  bufInt[IdxNExchange] = 0;
  bufInt[IdxNGutz] = 0;
  bufInt[IdxNJast] = 0;
  bufInt[IdxNDH2] = 0;
  bufInt[IdxNDH4] = 0;
  bufInt[IdxNOrbit] = 0;
  bufInt[IdxNQPTrans] = 0;
  bufInt[IdxNOneBodyG] = 0;
  bufInt[IdxNTwoBodyG] = 0;
  bufInt[IdxNTwoBodyGEx] = 0;
  bufInt[IdxNInterAll] = 0;
  bufInt[IdxNQPOptTrans] = 1;
  bufInt[IdxSROptCGMaxIter] = 0;
  bufInt[IdxNBF] = 0;
  bufInt[IdxNrange] = 0;
  bufInt[IdxNNz] = 0;
  bufInt[Idx2Sz] = -1;// -1: sz is not fixed :fsz
  bufInt[IdxNCond] = -1;

  bufDouble[IdxSROptRedCut] = 0.001;
  bufDouble[IdxSROptStaDel] = 0.02;
  bufDouble[IdxSROptStepDt] = 0.02;
  bufDouble[IdxSROptCGTol] = 1.0e-10;
  NStoreO = 1;
  NSRCG = 0;
}

int GetInfoFromModPara(int *bufInt, double *bufDouble) {

  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  char ctmp2[D_FileNameMax];

  int itmp;
  char *cerr;

  int iKWidx = 0;
  int iret = 0;
  fprintf(stdout, "Start: Read ModPara File .\n");
  for (iKWidx = 0; iKWidx < KWIdxInt_end; iKWidx++) {
    strcpy(defname, cFileNameListFile[iKWidx]);
    if (strcmp(defname, "") == 0) continue;
    fp = fopen(defname, "r");
    if (fp == NULL) {
      iret=ReadDefFileError(defname);
      fclose(fp);
      break;
    } else {
      switch (iKWidx) {
        case KWModPara:
          /* Read modpara.def---------------------------------------*/
          cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
          if(cerr == NULL) return -1;
          cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
          if(cerr == NULL) return -1;
          sscanf(ctmp2, "%s %d\n", ctmp, &itmp); //2
          cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //3
          if(cerr == NULL) return -1;
          cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //4
          if(cerr == NULL) return -1;
          cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //5
          if(cerr == NULL) return -1;
          cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
          if(cerr == NULL) return -1;
          sscanf(ctmp2, "%s %s\n", ctmp, CDataFileHead); //6
          cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
          if(cerr == NULL) return -1;
          sprintf(ctmp, "output/%s", CDataFileHead);
          strcpy(CDataFileHead, ctmp);
          sscanf(ctmp2, "%s %s\n", ctmp, CParaFileHead); //7
          sprintf(ctmp, "output/%s", CParaFileHead);
          strcpy(CParaFileHead, ctmp);
          cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);   //8
          if(cerr == NULL) return -1;

          iret = system("mkdir -p output");

          double dtmp;
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (*ctmp2 == '\n' || ctmp2[0] == '-') continue;
            sscanf(ctmp2, "%s %lf\n", ctmp, &dtmp);
            if (CheckWords(ctmp, "NVMCCalMode") == 0) {
              bufInt[IdxVMCCalcMode] = (int) dtmp;
            } else if (CheckWords(ctmp, "NLanczosMode") == 0) {
              bufInt[IdxLanczosMode] = (int) dtmp;
            } else if (CheckWords(ctmp, "NDataIdxStart") == 0) {
              bufInt[IdxDataIdxStart] = (int) dtmp;
            } else if (CheckWords(ctmp, "NDataQtySmp") == 0) {
              bufInt[IdxDataQtySmp] = (int) dtmp;
            } else if (CheckWords(ctmp, "Nsite") == 0) {
              bufInt[IdxNsite] = (int) dtmp;
            } else if (CheckWords(ctmp, "Ne") == 0 || CheckWords(ctmp, "Nelectron") == 0) {
              bufInt[IdxNe] = (int) dtmp;
            } else if (CheckWords(ctmp, "Ncond") == 0) {
              bufInt[IdxNCond] = (int) dtmp;
            } else if (CheckWords(ctmp, "2Sz") == 0) {
              bufInt[Idx2Sz] = (int) dtmp;
              if (bufInt[Idx2Sz] == -1) {
                fprintf(stdout, "Error: 2Sz must be even number.\n");
                return (-1);
              }
            } else if (CheckWords(ctmp, "NSPGaussLeg") == 0) {
              bufInt[IdxSPGaussLeg] = (int) dtmp;
            } else if (CheckWords(ctmp, "NSPStot") == 0) {
              bufInt[IdxSPStot] = (int) dtmp;
            } else if (CheckWords(ctmp, "NMPTrans") == 0) {
              bufInt[IdxMPTrans] = (int) dtmp;
            } else if (CheckWords(ctmp, "NSROptItrStep") == 0) {
              bufInt[IdxSROptItrStep] = (int) dtmp;
            } else if (CheckWords(ctmp, "NSROptItrSmp") == 0) {
              bufInt[IdxSROptItrSmp] = (int) dtmp;
            } else if (CheckWords(ctmp, "DSROptRedCut") == 0) {
              bufDouble[IdxSROptRedCut] = (double) dtmp;
            } else if (CheckWords(ctmp, "DSROptStaDel") == 0) {
              bufDouble[IdxSROptStaDel] = (double) dtmp;
            } else if (CheckWords(ctmp, "DSROptStepDt") == 0) {
              bufDouble[IdxSROptStepDt] = (double) dtmp;
            } else if (CheckWords(ctmp, "NSROptCGMaxIter") == 0) {
              bufInt[IdxSROptCGMaxIter] = (int) dtmp;
            } else if (CheckWords(ctmp, "DSROptCGTol") == 0) {
              bufDouble[IdxSROptCGTol] = (double) dtmp;
            } else if (CheckWords(ctmp, "NVMCWarmUp") == 0) {
              bufInt[IdxVMCWarmUp] = (int) dtmp;
            } else if (CheckWords(ctmp, "NVMCInterval") == 0) {
              bufInt[IdxVMCInterval] = (int) dtmp;
            } else if (CheckWords(ctmp, "NVMCSample") == 0) {
              bufInt[IdxVMCSample] = (int) dtmp;
            } else if (CheckWords(ctmp, "NExUpdatePath") == 0) {
              bufInt[IdxExUpdatePath] = (int) dtmp;
            } else if (CheckWords(ctmp, "RndSeed") == 0) {
              bufInt[IdxRndSeed] = (int) dtmp;
            } else if (CheckWords(ctmp, "NSplitSize") == 0) {
              bufInt[IdxSplitSize] = (int) dtmp;
            } else if (CheckWords(ctmp, "NStore") == 0) {
              NStoreO = (int) dtmp;
            } else if (CheckWords(ctmp, "NSRCG") == 0) {
              NSRCG = (int) dtmp;
            } else {
              fprintf(stderr, "  Error: keyword \" %s \" is incorrect. \n", ctmp);
              iret = ReadDefFileError(defname);
              return iret;
            }
          }
          if (bufInt[IdxRndSeed] < 0) {
            bufInt[IdxRndSeed] = (int) time(NULL);
            fprintf(stdout, "  remark: Seed = %d\n", bufInt[IdxRndSeed]);
          }
          break;//modpara file
        default:
          break;
      }
    }
  }
  fclose(fp);
  fprintf(stdout, "End: Read ModPara File .\n");
  return iret;
}
/**********************************************************************/

/**********************************/
/* [s] Read Parameters from file  */
/**********************************/
int GetTransferInfo(FILE *fp, int **ArrayIdx, double complex *ArrayValue, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0, x2 = 0, x3 = 0;
  double dReValue = 0, dImValue = 0;
  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d %lf %lf\n",
           &x0, &x1, &x2, &x3, &dReValue, &dImValue);
    ArrayIdx[idx][0] = x0;
    ArrayIdx[idx][1] = x1;
    ArrayIdx[idx][2] = x2;
    ArrayIdx[idx][3] = x3;

    if (CheckPairSite(x0, x2, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }
    ArrayValue[idx] = dReValue + I * dImValue;
    idx++;
  }
  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

int GetLocSpinInfo(FILE *fp, int *ArrayIdx, int Nsite, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d \n", &x0, &x1);
    ArrayIdx[x0] = x1;
    if (CheckSite(x0, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect.\n");
      info = 1;
      break;
    }
    idx++;
  }
  if (idx != Nsite) info = ReadDefFileError(defname);
  return info;
}

int GetInfoCoulombIntra(FILE *fp, int *ArrayIdx, double *ArrayValue, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0;
  double dReValue = 0;
  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %lf\n", &x0, &dReValue);
    ArrayIdx[idx] = x0;
    if (CheckSite(x0, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }
    ArrayValue[idx] = dReValue;
    idx++;
  }
  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

int ReadPairHopValue(FILE *fp, int **ArrayIdx, double *ArrayValue, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0;
  double dReValue = 0;
  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %lf\n", &x0, &x1, &dReValue);
    ArrayIdx[2 * idx][0] = x0;
    ArrayIdx[2 * idx][1] = x1;
    ArrayValue[2 * idx] = dReValue;
    if (CheckPairSite(x0, x1, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }
    ArrayIdx[2 * idx + 1][0] = x1;
    ArrayIdx[2 * idx + 1][1] = x0;
    ArrayValue[2 * idx + 1] = dReValue;
    idx++;
  }
  if (2 * idx != NArray) info = ReadDefFileError(defname);
  return info;
}

int ReadPairDValue(FILE *fp, int **ArrayIdx, double *ArrayValue, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0;
  double dReValue = 0;
  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %lf\n", &x0, &x1, &dReValue);
    ArrayIdx[idx][0] = x0;
    ArrayIdx[idx][1] = x1;
    if (CheckPairSite(x0, x1, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }
    ArrayValue[idx] = dReValue;
    idx++;
  }
  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

int GetInfoOptOrbitalParalell(FILE *fp, int *ArrayOpt, int iComplxFlag, int *iTotalOptCount, int fidx) {
  int i;
  int iLocalOptCount = 0;
  int idxOptFlag = 0;
  while (fscanf(fp, "%d ", &i) != EOF) {
    idxOptFlag = 2 * (fidx + 2 * iLocalOptCount);
    fscanf(fp, "%d\n", &(ArrayOpt[idxOptFlag])); // up-up real
    ArrayOpt[idxOptFlag + 1] = iComplxFlag; //  up-up imag
    ArrayOpt[idxOptFlag + 2] = ArrayOpt[idxOptFlag]; //  up-up imag
    ArrayOpt[idxOptFlag + 3] = iComplxFlag; //  up-up imag
    (iLocalOptCount)++;
    (*iTotalOptCount) += 2;
  }
  return (iLocalOptCount);
}

int GetInfoOpt(FILE *fp, int *ArrayOpt, int iComplxFlag, int *iTotalOptCount, int fidx) {
  int i;
  int iLocalOptCount = 0;
  while (fscanf(fp, "%d ", &i) != EOF) {
    fscanf(fp, "%d\n", &(ArrayOpt[2 * fidx])); // TBC real
    if(iComplxFlag>0){
      ArrayOpt[2 * fidx + 1] = ArrayOpt[2 * fidx]; //  TBC imaginary
    }
    fidx++;
    (iLocalOptCount)++;
    (*iTotalOptCount)++;
  }
  return (iLocalOptCount);
}


int GetInfoGutzwiller(FILE *fp, int *ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int Nsite, int NArray,
                      char *defname) {
  int idx0 = 0, idx1 = 0, info = 0;
  int i = 0;
  int fidx = 0;

  if (NArray > 0) {
    idx0 = idx1 = 0;
    while (fscanf(fp, "%d ", &i) != EOF) {
      fscanf(fp, "%d\n", &(ArrayIdx[i]));
      if (CheckSite(i, Nsite) != 0) {
        fprintf(stderr, "Error: Site index is incorrect. \n");
        info = 1;
        break;
      }
      idx0++;
      if (idx0 == Nsite) break;
    }
    fidx = 0;
    idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    if (idx0 != Nsite || idx1 != NArray) {
      info = ReadDefFileError(defname);
    }
  }
  return info;
}

int GetInfoJastrow(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nsite,
                   int NArray, char *defname) {
  int idx0 = 0, idx1 = 0, info = 0;
  int i = 0, j = 0;
  int fidx = _fidx;
  if (NArray > 0) {
    while (fscanf(fp, "%d %d ", &i, &j) != EOF) {
      if (i == j) {
        fprintf(stderr, "Error in %s: [Condition] i neq j\n", defname);
        info = 1;
        break;
      }
      if (CheckPairSite(i, j, Nsite) != 0) {
        fprintf(stderr, "Error: Site index is incorrect. \n");
        info = 1;
        break;
      }

      fscanf(fp, "%d\n", &(ArrayIdx[i][j]));
      ArrayIdx[i][i] = -1; // This case is Gutzwiller.
      idx0++;
      if (idx0 == Nsite * (Nsite - 1)) break;
    }
    idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    if (idx0 != Nsite * (Nsite - 1) || idx1 != NArray) {
      info = ReadDefFileError(defname);
    }
  }
  return info;
}

int
GetInfoDH2(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nsite, int NArray,
           char *defname) {
  int idx0 = 0, idx1 = 0, info = 0;
  int i = 0, x0 = 0, x1 = 0, n = 0;
  int fidx = _fidx;
  if (NArray > 0) {
    idx0 = idx1 = 0;
    while (fscanf(fp, "%d %d %d %d\n", &i, &(x0), &(x1), &n) != EOF) {
      ArrayIdx[n][2 * i] = x0;
      ArrayIdx[n][2 * i + 1] = x1;
      if (CheckSite(i, Nsite) != 0 || CheckPairSite(x0, x1, Nsite) != 0) {
        fprintf(stderr, "Error: Site index is incorrect. \n");
        info = 1;
        break;
      }
      idx0++;
      if (idx0 == Nsite * NArray) break;
    }
    idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    if (idx0 != Nsite * NArray || idx1 != 2 * 3 * NArray) {
      info = ReadDefFileError(defname);
    }
  }
  return info;
}

int
GetInfoDH4(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nsite, int NArray,
           char *defname) {
  int idx0 = 0, idx1 = 0, info = 0;
  int i = 0, x0 = 0, x1 = 0, x2 = 0, x3 = 0, n = 0;
  int fidx = _fidx;
  if (NArray > 0) {
    idx0 = idx1 = 0;
    while (fscanf(fp, "%d %d %d %d %d %d\n",
                  &i, &(x0), &(x1), &(x2), &(x3), &n) != EOF) {
      ArrayIdx[n][4 * i] = x0;
      ArrayIdx[n][4 * i + 1] = x1;
      ArrayIdx[n][4 * i + 2] = x2;
      ArrayIdx[n][4 * i + 3] = x3;
      if (CheckSite(i, Nsite) != 0 || CheckQuadSite(x0, x1, x2, x3, Nsite) != 0) {
        fprintf(stderr, "Error: Site index is incorrect. \n");
        info = 1;
        break;
      }
      idx0++;
      if (idx0 == Nsite * NArray) break;
    }
    idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    if (idx0 != Nsite * NArray || idx1 != 2 * 5 * NArray) {
      info = ReadDefFileError(defname);
    }
  }
  return info;
}

int GetInfoTransSym(FILE *fp, int **Array, int **ArraySgn, int **ArrayInv, double complex *ArrayPara,
                    int _APFlag, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int i = 0, j = 0;
  int itmp = 0, itmpsgn = 0;
  double dReValue = 0, dImValue = 0;

  if (NArray > 0) {
    for (i = 0; i < NArray; i++) {
      itmp = 0;
      dReValue = 0;
      dImValue = 0;
      fgets(ctmp2, D_CharTmpReadDef, fp);
      sscanf(ctmp2, "%d %lf %lf\n", &itmp, &dReValue, &dImValue);
      ArrayPara[itmp] = dReValue + I * dImValue;
    }
    idx = 0;
    while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
      sscanf(ctmp2, "%d %d %d %d\n", &i, &j, &itmp, &itmpsgn);
      Array[i][j] = itmp;
      ArraySgn[i][j] = itmpsgn;
      ArrayInv[i][itmp] = j;
      idx++;
    }
    if (_APFlag == 0) {
      for (i = 0; i < NArray; i++) {
        for (j = 0; j < Nsite; j++) ArraySgn[i][j] = 1;
      }
    }
    if (idx != Nsite * NArray) info = ReadDefFileError(defname);
  }
  return info;
}

int
GetInfoOneBodyG(FILE *fp, int **ArrayIdx, int **ArrayToIdx, int _NLanczosMode, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0, x2 = 0, x3 = 0;
  int isite1 = 0, isite2 = 0;
  int idxLanczos = 0;
  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d\n", &x0, &x1, &x2, &x3);
    if (_NLanczosMode < 2) { // Normal
      ArrayIdx[idx][0] = x0;
      ArrayIdx[idx][1] = x1;
      ArrayIdx[idx][2] = x2;
      ArrayIdx[idx][3] = x3;
    } else { //For Calc Green func by Lanczos mode
      isite1 = x0 + x1 * Nsite;
      isite2 = x2 + x3 * Nsite;
      idxLanczos = ArrayToIdx[isite1][isite2];
      ArrayIdx[idxLanczos][0] = x0;
      ArrayIdx[idxLanczos][1] = x1;
      ArrayIdx[idxLanczos][2] = x2;
      ArrayIdx[idxLanczos][3] = x3;
    }
    if (CheckPairSite(x0, x2, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }
    idx++;
  }

  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

int GetInfoTwoBodyGEx(FILE *fp, int **ArrayIdx, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0, x2 = 0, x3 = 0;
  int x4 = 0, x5 = 0, x6 = 0, x7 = 0;
  //Debug
  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d %d %d %d %d\n", &x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7);
    ArrayIdx[idx][0] = x0; // Index to OneBodyG1
    ArrayIdx[idx][1] = x1; // Index to OneBodyG2
    ArrayIdx[idx][2] = x2; // G1:site i
    ArrayIdx[idx][3] = x3; // G1:site j
    ArrayIdx[idx][4] = x4; // G1:sigma1
    ArrayIdx[idx][5] = x5; // G2:site l
    ArrayIdx[idx][6] = x6; // G2:site k
    ArrayIdx[idx][7] = x7; // G2:sigma2
    if (CheckQuadSite(x2, x3, x5, x6, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }
    idx++;
  }
  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

int GetInfoTwoBodyG(FILE *fp, int **ArrayIdx, int **ArrayIdxTwoBodyGLz, int **ArrayToIdx, int **ArrayIdxOneBodyG,
                    int _NLanczosMode, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0, x2 = 0, x3 = 0;
  int x4 = 0, x5 = 0, x6 = 0, x7 = 0;
  int isite1 = 0, isite2 = 0;
  int idxLanczos = 0;

  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d %d %d %d %d\n", &x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7);
    ArrayIdx[idx][0] = x0; //G1: site i
    ArrayIdx[idx][1] = x1; //G1: sigma i
    ArrayIdx[idx][2] = x2; //G1: site j
    ArrayIdx[idx][3] = x3; //G1: sigma j
    ArrayIdx[idx][4] = x4; //G2: site k
    ArrayIdx[idx][5] = x5; //G2: sigma k
    ArrayIdx[idx][6] = x6; //G2: site l
    ArrayIdx[idx][7] = x7; //G2: sigma l
    if (CheckQuadSite(x0, x2, x4, x6, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }

    if (_NLanczosMode > 1) { //Calc TwoBodyG by Lanczos method
      
      isite1 = x0 + x1 * Nsite;
      isite2 = x2 + x3 * Nsite;      
      idxLanczos = ArrayToIdx[isite1][isite2];
      
      ArrayIdxOneBodyG[idxLanczos][0] = x0;
      ArrayIdxOneBodyG[idxLanczos][1] = x1;
      ArrayIdxOneBodyG[idxLanczos][2] = x2;
      ArrayIdxOneBodyG[idxLanczos][3] = x3;
      /*
      ArrayIdxOneBodyG[idxLanczos][0] = x2;
      ArrayIdxOneBodyG[idxLanczos][1] = x3;
      ArrayIdxOneBodyG[idxLanczos][2] = x0;
      ArrayIdxOneBodyG[idxLanczos][3] = x1;
      */
      ArrayIdxTwoBodyGLz[idx][0] = idxLanczos;

      /*
      isite1 = x4 + x5 * Nsite;
      isite2 = x6 + x7 * Nsite;
      */
      isite1 = x6 + x7 * Nsite;
      isite2 = x4 + x5 * Nsite;
      idxLanczos = ArrayToIdx[isite1][isite2];
      ArrayIdxOneBodyG[idxLanczos][0] = x6;
      ArrayIdxOneBodyG[idxLanczos][1] = x7;
      ArrayIdxOneBodyG[idxLanczos][2] = x4;
      ArrayIdxOneBodyG[idxLanczos][3] = x5;

      /*
      ArrayIdxOneBodyG[idxLanczos][0] = x4;
      ArrayIdxOneBodyG[idxLanczos][1] = x5;
      ArrayIdxOneBodyG[idxLanczos][2] = x6;
      ArrayIdxOneBodyG[idxLanczos][3] = x7;
      */
      ArrayIdxTwoBodyGLz[idx][1] = idxLanczos;
    }
    idx++;
  }
  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

int GetInfoOptTrans(FILE *fp, int **Array, double *ArrayPara, int *ArrayOpt, int **ArraySgn,
                    int _iFlagOptTrans, int *iOptCount, int _fidx, int _APFlag, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int i = 0, j = 0;
  int itmp = 0, itmpsgn = 0;
  int fidx = _fidx;
  double dReValue;
  if (_iFlagOptTrans > 0) {
    for (i = 0; i < NArray; i++) {
      itmp = 0;
      dReValue = 0;
      fgets(ctmp2, D_CharTmpReadDef, fp);
      sscanf(ctmp2, "%d %lf \n", &itmp, &dReValue);
      ArrayPara[itmp] = dReValue;
      ArrayOpt[fidx] = 1;
      fidx++;
      (*iOptCount)++;
    }
    idx = 0;
    while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
      sscanf(ctmp2, "%d %d %d %d\n", &i, &j, &itmp, &itmpsgn);
      Array[i][j] = itmp;
      ArraySgn[i][j] = itmpsgn;
      idx++;
    }
    if (_APFlag == 0) {
      for (i = 0; i < NArray; i++) {
        for (j = 0; j < Nsite; j++) ArraySgn[i][j] = 1;
      }
    }
  }
  return info;
}

int GetInfoInterAll(FILE *fp, int **ArrayIdx, double complex *ArrayValue,
                    int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0, x2 = 0, x3 = 0;
  int x4 = 0, x5 = 0, x6 = 0, x7 = 0;
  double dReValue = 0, dImValue = 0;

  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d %d %d %d %d %lf %lf\n",
           &x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7,
           &dReValue, &dImValue);

    ArrayIdx[idx][0] = x0;
    ArrayIdx[idx][1] = x1;
    ArrayIdx[idx][2] = x2;
    ArrayIdx[idx][3] = x3;
    ArrayIdx[idx][4] = x4;
    ArrayIdx[idx][5] = x5;
    ArrayIdx[idx][6] = x6;
    ArrayIdx[idx][7] = x7;

    if (CheckQuadSite(x0, x2, x4, x6, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }

    ArrayValue[idx] = dReValue + I * dImValue;

    if (TwoSz != -1 && !(x1 == x3 && x5 == x7)) {
      fprintf(stderr, "  Error:  Sz non-conserved system is not yet supported for InterAll.\n");
      info = ReadDefFileError(defname);
      break;
    }
    idx++;
  }
  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

int GetInfoOrbitalAntiParallel(FILE *fp, int **Array, int *ArrayOpt, int **ArraySgn, int *iOptCount,
                               int _fidx, int _iComplexFlag, int _iFlagOrbitalGeneral, int _APFlag, int Nsite,
                               int NArray, char *defname) {
  char ctmp2[256];
  int i, j;
  int idx0 = 0, idx1 = 0;
  int itmp = 0, info = 0;
  int spn_i, spn_j;
  int all_i, all_j;
  int fij = 0, fijSign = 1;
  int fidx = _fidx;

  if (NArray == 0) return 0;

  if (_iFlagOrbitalGeneral == 0) {
    while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
      sscanf(ctmp2, "%d %d %d %d \n", &i, &j, &fij, &fijSign);
      if (CheckPairSite(i, j, Nsite) != 0) {
        fprintf(stderr, "Error: Site index is incorrect. \n");
        return -1;
      }
      idx0++;
      Array[i][j] = fij;
      ArraySgn[i][j] = fijSign;
      if (idx0 == Nsite * Nsite) break;
    }

    if (_APFlag == 0) {
      for (i = 0; i < Nsite; i++) {
        for (j = 0; j < Nsite; j++) {
          ArraySgn[i][j] = 1;
        }
      }
    }
  } else { //_iFlagOrbitalGeneral == 0
    while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
      sscanf(ctmp2, "%d %d %d %d \n", &i, &j, &fij, &fijSign);
      spn_i = 0;
      spn_j = 1;
      all_i = i + spn_i * Nsite; //fsz
      all_j = j + spn_j * Nsite; //fsz
      if (CheckPairSite(i, j, Nsite) != 0) {
        fprintf(stderr, "Error: Site index is incorrect. \n");
        return -1;
      }
      if (all_i >= all_j) itmp = 1;
      idx0++;
      Array[all_i][all_j] = fij;
      ArraySgn[all_i][all_j] = fijSign;
      // Note F_{IJ}=-F_{JI}
      Array[all_j][all_i] = fij;
      ArraySgn[all_j][all_i] = -fijSign;
      if (idx0 == (Nsite * Nsite)) break;
    }

    if (_APFlag == 0) {
      for (i = 0; i < Nsite; i++) {
        for (j = Nsite; j < 2 * Nsite; j++) {
          ArraySgn[i][j] = 1;
          ArraySgn[j][i] = -1;
        }
      }
    }
  }

  idx1 = GetInfoOpt(fp, ArrayOpt, _iComplexFlag, iOptCount, fidx);
  if (idx0 != Nsite * Nsite || idx1 != NArray || itmp == 1) {
    info = ReadDefFileError(defname);
  }

  return info;
}

int GetInfoOrbitalGeneral(FILE *fp, int **Array, int *ArrayOpt, int **ArraySgn, int *iOptCount,
                          int _fidx, int _iComplexFlag, int _iFlagOrbitalGeneral, int _APFlag, int Nsite, int NArray,
                          char *defname) {
  char ctmp2[256];
  int i, j;
  int idx0 = 0, idx1 = 0;
  int itmp = 0, info = 0;
  int spn_i, spn_j;
  int all_i, all_j;
  int fij = 0, fijSign = 1;
  int fidx = _fidx;

  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d %d %d\n", &i, &spn_i, &j, &spn_j, &fij, &fijSign);
    all_i = i + spn_i * Nsite; //fsz
    all_j = j + spn_j * Nsite; //fsz
    if (CheckPairSite(i, j, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      return -1;
    }
    if (all_i >= all_j) itmp = 1;
    idx0++;
    Array[all_i][all_j] = fij;
    ArraySgn[all_i][all_j] = fijSign;
    // Note F_{IJ}=-F_{JI}
    Array[all_j][all_i] = fij;
    ArraySgn[all_j][all_i] = -fijSign;
    if (idx0 == (2 * Nsite * Nsite - Nsite)) break; //2N*(2N-1)/2
  }

  if (_APFlag == 0) {
    for (i = 0; i < 2 * Nsite; i++) {
      for (j = i + 1; j < 2 * Nsite; j++) {
        ArraySgn[i][j] = 1;
        ArraySgn[j][i] = -1;
      }
    }
  }

  idx1 = GetInfoOpt(fp, ArrayOpt, _iComplexFlag, iOptCount, fidx);
  if (idx0 != (2 * Nsite * Nsite - Nsite) || idx1 != NArray || itmp == 1) {
    info = ReadDefFileError(defname);
  }

  return info;
}

int GetInfoOrbitalParallel(FILE *fp, int **Array, int *ArrayOpt, int **ArraySgn, int *iOptCount,
                           int _fidx, int _iComplexFlag, int _iFlagOrbitalGeneral, int _APFlag, int Nsite, int NArray,
                           int NArrayAP, char *defname) {
  char ctmp2[256];
  int i, j;
  int idx0 = 0, idx1 = 0;
  int itmp = 0, info = 0;
  int spn_i;
  int all_i, all_j;
  int fij = 0, fijSign = 1, fij_org;
  int fidx = _fidx;

  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d\n", &i, &j, &fij_org, &fijSign);

    if (CheckPairSite(i, j, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      return -1;
    }

    for (spn_i = 0; spn_i < 2; spn_i++) {
      all_i = i + spn_i * Nsite; //fsz
      all_j = j + spn_i * Nsite; //fsz
      if (all_i >= all_j) itmp = 1;
      idx0++;
      fij = NArrayAP + 2 * fij_org + spn_i;

      Array[all_i][all_j] = fij;
      ArraySgn[all_i][all_j] = fijSign;
      // Note F_{IJ}=-F_{JI}
      Array[all_j][all_i] = fij;
      ArraySgn[all_j][all_i] = -fijSign;
    }
    if (idx0 == (Nsite * (Nsite - 1))) break;
  }

  if (_APFlag == 0) {
    for (spn_i = 0; spn_i < 2; spn_i++) {
      for (i = 0; i < Nsite; i++) {
        for (j = i + 1; j < Nsite; j++) {
          all_i = i + spn_i * Nsite; //fsz
          all_j = j + spn_i * Nsite; //fsz
          ArraySgn[all_i][all_j] = 1;
          ArraySgn[all_j][all_i] = -1;
        }
      }
    }
  }


  idx1 = GetInfoOptOrbitalParalell(fp, ArrayOpt, _iComplexFlag, iOptCount, fidx);
  if (idx0 != (Nsite * (Nsite - 1)) || idx1 != NArray || itmp == 1) {
    info = ReadDefFileError(defname);
  }

  return info;
}
/**********************************/
/* [e] Read Parameters from file  */
/**********************************/
