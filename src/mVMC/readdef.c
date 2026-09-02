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
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "./include/backflow.h"
#include "./include/lanczos2_contract.h"
#include "./include/readdef.h"
#include "./include/global.h"
#include "safempi_fcmp.c"

char (*cFileNameListFile)[D_CharTmpReadDef] = NULL;

#define GC_ANOMALOUS_HERMITE_EPS 1.0e-12

int ReadDefFileError(const char *defname);

int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm);

int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm);


int CheckSite(const int iSite, const int iMaxNum);

int CheckPairSite(const int iSite1, const int iSite2, const int iMaxNum);

int CheckQuadSite(const int iSite1, const int iSite2, const int iSite3, const int iSite4, const int iMaxNum);

int GetTransferInfo(FILE *fp, int **ArrayIdx, double complex *ArrayValue, int Nsite, int NArray, char *defname);

int GetLocSpinInfo(FILE *fp, int *ArrayIdx, int Nsite, int NLocalSpin,
                   char *defname);

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

int GetInfoAnomalousTerm(FILE *fp, int **ArrayIdx,
                         double complex *ArrayValue, int Nsite, int NArray,
                         char *defname);

int GetInfoAnomalousG(FILE *fp, int **ArrayIdx, int Nsite, int NArray,
                      char *defname);

int GetInfoOptTrans(FILE *fp, int **Array, double *ArrayPara, int *ArrayOpt, int **ArraySgn,
                    int _iFlagOptTrans, int *iOptCount, int _fidx, int _APFlag, int Nsite, int NArray, char *defname);

int GetInfoTwoBodyG(FILE *fp, int **ArrayIdx, int Nsite, int NArray, char *defname);

int ReadBuffNBodyG(FILE *fp, int *nbody, int *totalFactors,
                   int *maxN, int nsite, char *defname);

int GetInfoNBodyG(FILE *fp, int *termN, int *termOffset, int **termIdx,
                  int nsite, int nbody, int totalFactors,
                  int maxN, int allowSpinChange, char *defname);

int ReadBuffNBodyInterAll(FILE *fp, int *nbody, int *totalFactors,
                          int *maxN, int nsite, char *defname);

int GetInfoNBodyInterAll(FILE *fp, int *termN, int *termOffset,
                         int **termIdx, double complex *termCoef,
                         int nsite, int nbody, int totalFactors,
                         int maxN, int allowSpinChange, char *defname);

int GetInfoTwoBodyGEx(FILE *fp, int **ArrayIdx, int **ArrayToIdx, int **ArrayIdxOneBodyG,
                      int Nsite, int NArray, char *defname);

int GetInfoOrbitalGeneral(FILE *fp, int **Array, int *ArrayOpt, int **ArraySgn, int *iOptCount,
                          int _fidx, int _iComplexFlag, int _iFlagOrbitalGeneral, int _APFlag, int Nsite, int NArray,
                          char *defname);

int
GetInfoOneBodyG(FILE *fp, int **ArrayIdx, int **ArrayToIdx, int _NLanczosMode, int Nsite, int NArray, char *defname);

//RBM
int GetInfoRBM_Layer(FILE *fp, int *ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nlayer, int NArray,
                      char *defname);
int GetInfoGeneralRBM_Layer(FILE *fp, int *ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nlayer, int NArray,
                      char *defname);
int GetInfoRBM_PhysHidden(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nlayer, int Nlayer2, int NArray,
                      char *defname);
int GetInfoGeneralRBM_PhysHidden(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nlayer, int Nlayer2, int NArray,
                      char *defname);
//RBM

int GetInfoLattice(FILE *fp, int **ArrayIdx, int NArray, int nx, int ny, int nz, int norb, char *defname);
//int GetInfoTwist(FILE *fp, int **ArrayIdx, int NArray, char *defname);
int GetInfoTwist(FILE *fp, int **ArrayIdx, double **ArrayValue, int Nsite, int NTwist, char *defname);

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

static int LineTailIsWhitespace(const char *line, int offset) {
  const unsigned char *cursor;
  if (line == NULL || offset < 0) return 0;
  cursor = (const unsigned char *)line + offset;
  while (*cursor != '\0') {
    if (!isspace(*cursor)) return 0;
    cursor++;
  }
  return 1;
}

static int ParseStrictLongLong(const char *text, long long *value) {
  char *end = NULL;
  long long parsed;
  if (text == NULL || value == NULL || *text == '\0') return 1;
  errno = 0;
  parsed = strtoll(text, &end, 10);
  if (errno == ERANGE || end == text || *end != '\0') return 1;
  *value = parsed;
  return 0;
}

static int ReadAnomalousCount(FILE *fp, int *result,
                              const char *expectedName,
                              const char *defname) {
  char header[D_FileNameMax];
  char line[D_FileNameMax];
  char parsedName[D_FileNameMax];
  char parsedValue[D_FileNameMax];
  long long value;
  int offset = 0;
  if (fp == NULL || result == NULL || expectedName == NULL ||
      fgets(header, sizeof(header), fp) == NULL ||
      fgets(line, sizeof(line), fp) == NULL ||
      sscanf(line, "%255s %255s%n", parsedName, parsedValue, &offset) != 2 ||
      !LineTailIsWhitespace(line, offset) ||
      strcmp(parsedName, expectedName) != 0 ||
      ParseStrictLongLong(parsedValue, &value) != 0 ||
      value < 0 || value > INT_MAX) {
    fprintf(stderr,
            "Error in %s: %s must be an integer in [0, INT_MAX].\n",
            defname, expectedName);
    return 1;
  }
  *result = (int)value;
  return 0;
}

static int ReadPositiveProjectionCount(FILE *fp, int *count,
                                       const char *expectedKeyword,
                                       const char *defname) {
  char header[D_FileNameMax];
  char line[D_FileNameMax];
  char keyword[D_FileNameMax];
  char valueText[D_FileNameMax];
  char extra;
  long long parsed;
  int fields;

  if (fp == NULL || count == NULL || expectedKeyword == NULL) return 1;
  if (fgets(header, sizeof(header), fp) == NULL ||
      fgets(line, sizeof(line), fp) == NULL) {
    fprintf(stderr, "Error: incomplete %s header in %s.\n",
            expectedKeyword, defname);
    return 1;
  }
  fields = sscanf(line, "%255s %255s %c", keyword, valueText, &extra);
  if (fields != 2 || strcmp(keyword, expectedKeyword) != 0 ||
      ParseStrictLongLong(valueText, &parsed) != 0 ||
      parsed <= 0 || parsed > INT_MAX) {
    fprintf(stderr,
            "Error: %s must be an integer in [1, %d] in %s.\n",
            expectedKeyword, INT_MAX, defname);
    return 1;
  }
  *count = (int)parsed;
  return 0;
}

static int CheckedCountTerm(long long *total, long long a, long long b,
                            long long c, const char *name) {
  long long product;
  if (total == NULL || a < 0 || b < 0 || c < 0) {
    fprintf(stderr, "Error: negative size while computing %s.\n", name);
    return 1;
  }
  if ((a != 0 && b > LLONG_MAX / a) ||
      (a * b != 0 && c > LLONG_MAX / (a * b))) {
    fprintf(stderr, "Error: integer overflow while computing %s.\n", name);
    return 1;
  }
  product = a * b * c;
  if (product > LLONG_MAX - *total) {
    fprintf(stderr, "Error: integer overflow while accumulating %s.\n", name);
    return 1;
  }
  *total += product;
  return 0;
}

static int CheckedIntProduct3(int a, int b, int c, const char *name,
                              int *result) {
  long long value;
  if (result == NULL || a < 0 || b < 0 || c < 0) {
    fprintf(stderr, "Error: negative size while computing %s.\n", name);
    return 1;
  }
  value = (long long)a * (long long)b * (long long)c;
  if (value > INT_MAX) {
    fprintf(stderr, "Error: %s exceeds the supported int range.\n", name);
    return 1;
  }
  *result = (int)value;
  return 0;
}

static int CheckProjectionAllocationSizes(int nsite, int nQPTrans) {
  size_t count;
  if (nsite <= 0 || nQPTrans <= 0) return 1;
  count = (size_t)nQPTrans;
  if (count > SIZE_MAX / sizeof(int *) ||
      count > SIZE_MAX / sizeof(double complex) ||
      (size_t)nsite > SIZE_MAX / count) {
    fprintf(stderr,
            "Error: allocation size overflow for translation projection tables.\n");
    return 1;
  }
  return 0;
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

static int ReadBackFlowRangeDefinition(const char *defname) {
  FILE *fp;
  char ctmp[D_FileNameMax];
  int i;
  int info = 0;

  if (defname == NULL || strcmp(defname, "") == 0) return ReadDefFileError("BFRange");
  fp = fopen(defname, "r");
  if (fp == NULL) return ReadDefFileError(defname);
  for (i = 0; i < IgnoreLinesInDef; i++) {
    if (fgets(ctmp, sizeof(ctmp) / sizeof(char), fp) == NULL) {
      info = ReadDefFileError(defname);
      fclose(fp);
      return info;
    }
  }
  info = BFReadRange(fp, defname);
  fclose(fp);
  return info;
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

int ReadGreen(char *xNameListFile, int Nca, int **caIdx, int Ncacadc, int **cacaIdx, int Ns) {
  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  int iKWidx = 0;
  char *cerr;
  int i, info = 0;

  if (cFileNameListFile != NULL) {
    free(cFileNameListFile);
    cFileNameListFile = NULL;
  }
  cFileNameListFile = malloc(sizeof(char) * D_CharTmpReadDef * KWIdxInt_end);
  if (cFileNameListFile == NULL) {
    fprintf(stderr, "  error: Memory allocation failed for filename list.\n");
    return -1;
  }
  // freed at vmcmain.c
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
      continue;
    }

    /*=======================================================================*/
    for (i = 0; i < IgnoreLinesInDef; i++) {
      cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      if (cerr == NULL) {
        fclose(fp);
        return (-1);
      }
    }
    switch (iKWidx) {
      case KWOneBodyG:
        /*cisajs.def----------------------------------------*/
        if (GetInfoOneBodyG(fp, caIdx, iOneBodyGIdx, 0, Ns, Nca, defname) != 0) {
          fclose(fp);
          return (-1);
        }
        break;
      case KWTwoBodyGEx:
        /*cisajscktalt.def----------------------------------*/
        /*load as if it's DC for index rearranging----------*/
        if (GetInfoTwoBodyG(fp, cacaIdx, Ns, Ncacadc, defname) != 0) {
          fclose(fp);
          return (-1);
        }
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
int CountOneBodyGForLanczos(char *xNameListFile, int Nca, int Ncacadc, int Ns, int **iFlgOneBodyG) {

  int info = 0;
  int i, j, isite1, isite2;
  int icount = 0;
  int **cacaIdx;
  int **caIdx;

  cacaIdx = malloc(sizeof(int *) * Ncacadc);
  for (i = 0; i < Ncacadc; i++)
    cacaIdx[i] = malloc(sizeof(int) * 8);
  caIdx = malloc(sizeof(int *) * Nca);
  for (i = 0; i < Nca; i++)
    caIdx[i] = malloc(sizeof(int) * 4);

  for (i = 0; i < 2 * Ns; i++) {
    for (j = 0; j < 2 * Ns; j++) {
      iFlgOneBodyG[i][j] = -1;
    }
  }
  info = ReadGreen(xNameListFile, Nca, caIdx, Ncacadc, cacaIdx, Ns);
  if (info != 0) {
    for (i = 0; i < Ncacadc; i++)
      free(cacaIdx[i]);
    free(cacaIdx);
    for (i = 0; i < Nca; i++)
      free(caIdx[i]);
    free(caIdx);
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
    isite1 = cacaIdx[i][0] + cacaIdx[i][1] * Ns;
    isite2 = cacaIdx[i][2] + cacaIdx[i][3] * Ns;
    if (iFlgOneBodyG[isite1][isite2] == -1) {
      iFlgOneBodyG[isite1][isite2] = icount;
      icount++;
    }

    /*
    isite1 = cacaDCIdx[i][4] + cacaDCIdx[i][5] * Ns;
    isite2 = cacaDCIdx[i][6] + cacaDCIdx[i][7] * Ns;
    */
    isite1 = cacaIdx[i][6] + cacaIdx[i][7] * Ns;
    isite2 = cacaIdx[i][4] + cacaIdx[i][5] * Ns;
    if (iFlgOneBodyG[isite1][isite2] == -1) {
      iFlgOneBodyG[isite1][isite2] = icount;
      icount++;
    }

  }
  for (i = 0; i < Ncacadc; i++)
    free(cacaIdx[i]);
  free(cacaIdx);
  for (i = 0; i < Nca; i++)
    free(caIdx[i]);
  free(caIdx);
  return icount;
}

static int ValidateNBodyCapabilities(int nBackFlow, int nNBodyG,
                                     int nNBodyInterAll,
                                     int allComplexFlag, int lanczosMode,
                                     int reportErrors) {
  int info = 0;

  if (nBackFlow > 0 && nNBodyG > 0) {
    if (allComplexFlag == 0) {
      if (reportErrors) {
        fprintf(stderr,
                "Error: BackFlow NBodyG requires complex variational "
                "parameters.\n");
      }
      info = 1;
    } else if (lanczosMode > 0) {
      if (reportErrors) {
        fprintf(stderr,
                "Error: BackFlow NBodyG is not supported with Lanczos mode.\n");
      }
      info = 1;
    }
  }
  if (nBackFlow > 0 && nNBodyInterAll > 0) {
    if (allComplexFlag == 0) {
      if (reportErrors) {
        fprintf(stderr,
                "Error: BackFlow NBodyInterAll requires complex variational "
                "parameters.\n");
      }
      info = 1;
    } else if (lanczosMode > 0) {
      if (reportErrors) {
        fprintf(stderr,
                "Error: BackFlow NBodyInterAll is not supported with "
                "Lanczos mode.\n");
      }
      info = 1;
    }
  }
  if (nNBodyInterAll > 0 && nBackFlow == 0 && lanczosMode > 0) {
    if (reportErrors) {
      fprintf(stderr,
              "Error: NBodyInterAll is not implemented for Lanczos mode.\n");
    }
    info = 1;
  }
  if (nNBodyInterAll > 0 && nBackFlow == 0 && allComplexFlag == 0) {
    if (reportErrors) {
      fprintf(stderr,
              "Error: NBodyInterAll requires complex variational parameters "
              "because the real local-energy kernels are not implemented "
              "(NNBodyInterAll=%d).\n",
              nNBodyInterAll);
    }
    info = 1;
  }
  return info;
}

int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm) {
  FILE *fp;
  char defname[D_FileNameMax];
  char *cerr;
  char ctmp[D_FileNameMax];
  char ctmp2[D_FileNameMax];
  char updateWeightError[512];

  int rank, size, info = 0;
  const int nBufInt = ParamIdxInt_End;
  const int nBufDouble = ParamIdxDouble_End;
  const int nBufChar = D_FileNameMax;
  int bufInt[nBufInt];
  double bufDouble[nBufDouble];
  int iKWidx = 0;
  int iret = 0;
  int hasLattice = 0;
  int hasBF = 0;
  int hasBFRange = 0;
  int hasAnomalousTerm = 0;
  int hasAnomalousG = 0;
  int iFlgOrbitalAntiParallel = 0;
  int iFlgOrbitalParallel = 0;
  int itmp = 0;

  int iOrbitalComplex = 0;
  iFlgOrbitalGeneral = 0;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (rank == 0) {
    cFileNameListFile = malloc(sizeof(char) * D_CharTmpReadDef * KWIdxInt_end);
    fprintf(stdout, "  Read File %s .\n", xNameListFile);
    if (GetFileName(xNameListFile, cFileNameListFile) != 0) {
      fprintf(stderr, "  error: Definition files(*.def) are incomplete.\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    FlagUpdateWeight =
        (strcmp(cFileNameListFile[KWInUpdateWeight], "") != 0);
    UpdateWeights[UpdateWeightExchange] = 0.0;
    UpdateWeights[UpdateWeightLocalSpinFlip] = 0.0;
    UpdateWeights[UpdateWeightPairSpinFlip] = 0.0;

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
    hasLattice = (strcmp(cFileNameListFile[KWLattice], "") != 0);
    hasBF = (strcmp(cFileNameListFile[KWBF], "") != 0);
    hasBFRange = (strcmp(cFileNameListFile[KWBFRange], "") != 0);
    hasAnomalousTerm =
        (strcmp(cFileNameListFile[KWAnomalousTerm], "") != 0);
    hasAnomalousG = (strcmp(cFileNameListFile[KWAnomalousG], "") != 0);

    SetDefaultValuesModPara(bufInt, bufDouble);
    iret = GetInfoFromModPara(bufInt, bufDouble);
    if (iret != 0) {
      if (rank == 0) {
        fprintf(stderr, "  Error: ModPara file is incomplete.\n");
      }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
 /*
    Block_size is using for Tuning of RBMRatio function.
    It must be multiple of 8.
 */
    if((bufInt[IdxNBlockSize_RBMRatio]%8)>0) {
      itmp=bufInt[IdxNBlockSize_RBMRatio]/8;
      if(itmp==0) itmp=1;
      bufInt[IdxNBlockSize_RBMRatio]=itmp*8;
      printf("\n --- Notice --- \n");
      printf("NBlockSize_RBMRatio Adjust to %d\n\n",bufInt[IdxNBlockSize_RBMRatio]);
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

          case KWSpinJastrow:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNSpinJast], &iComplexFlgSpinJastrow);
            break;

//RBM
          case KWGeneralRBM_HiddenLayer:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNGeneralRBM_HiddenLayer], &iComplexFlgGeneralRBM_HiddenLayer);
            FlagRBM=1;
            break;
          case KWChargeRBM_HiddenLayer:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNChargeRBM_HiddenLayer], &iComplexFlgChargeRBM_HiddenLayer);
            FlagRBM=1;
            break;
          case KWSpinRBM_HiddenLayer:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNSpinRBM_HiddenLayer], &iComplexFlgSpinRBM_HiddenLayer);
            FlagRBM=1;
            break;

          case KWGeneralRBM_PhysLayer:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNGeneralRBM_PhysLayer], &iComplexFlgGeneralRBM_PhysLayer);
            FlagRBM=1;
            break;
          case KWChargeRBM_PhysLayer:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNChargeRBM_PhysLayer], &iComplexFlgChargeRBM_PhysLayer);
            FlagRBM=1;
            break;
          case KWSpinRBM_PhysLayer:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNSpinRBM_PhysLayer], &iComplexFlgSpinRBM_PhysLayer);
            FlagRBM=1;
            break;

          case KWGeneralRBM_PhysHidden:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNGeneralRBM_PhysHidden], &iComplexFlgGeneralRBM_PhysHidden);
            FlagRBM=1;
            break;
          case KWChargeRBM_PhysHidden:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNChargeRBM_PhysHidden], &iComplexFlgChargeRBM_PhysHidden);
            FlagRBM=1;
            break;
          case KWSpinRBM_PhysHidden:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNSpinRBM_PhysHidden], &iComplexFlgSpinRBM_PhysHidden);
            FlagRBM=1;
            break;
//RBM

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
            iComplexFlgOrbitalAntiParallel = iOrbitalComplex;
            iComplexFlgOrbital += iOrbitalComplex;
            break;

          case KWOrbitalParallel:
            cerr = ReadBuffIntCmpFlg(fp, &iNOrbitalParallel, &iOrbitalComplex);
            iFlgOrbitalParallel = 1;
            bufInt[IdxNOrbit] += 2 * iNOrbitalParallel; //up-up and down-down
            iComplexFlgOrbitalParallel = iOrbitalComplex;
            iComplexFlgOrbital += iOrbitalComplex;
            break;

          case KWOrbitalGeneral:
            cerr = ReadBuffIntCmpFlg(fp, &bufInt[IdxNOrbit], &iOrbitalComplex);
            iFlgOrbitalGeneral = 1;
            iComplexFlgOrbitalGeneral = iOrbitalComplex;
            iComplexFlgOrbital += iOrbitalComplex;
            break;

          case KWTransSym:
            cerr = "";
            if (ReadPositiveProjectionCount(fp, &bufInt[IdxNQPTrans],
                                            "NQPTrans", defname) != 0) {
              info = ReadDefFileError(defname);
            }
            break;

          case KWInUpdateWeight:
            cerr = "";
            updateWeightError[0] = '\0';
            if (ReadUpdateWeight(fp, defname, UpdateWeights,
                                 updateWeightError,
                                 sizeof(updateWeightError)) != 0) {
              fprintf(stderr, "Error: %s.\n", updateWeightError);
              info = 1;
            }
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

          case KWNBodyG:
            cerr = "";
            if (ReadBuffNBodyG(fp,
                               &bufInt[IdxNNBodyG],
                               &bufInt[IdxNBodyGTotalFactors],
                               &bufInt[IdxNBodyGMaxN],
                               bufInt[IdxNsite],
                               defname) != 0) {
              info = ReadDefFileError(defname);
            }
            break;

          case KWLattice:
            cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
            if (cerr != NULL) {
              cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
              if (cerr != NULL) {
                int scanned_lat = sscanf(ctmp2, "%s %d %d %d %d\n",
                                         ctmp,
                                         &bufInt[IdxNx], &bufInt[IdxNy],
                                         &bufInt[IdxNz], &bufInt[IdxNorb]);
                if (scanned_lat != 5) {
                  fprintf(stderr,
                          "Error: malformed lattice header (expected 'name Nx Ny Nz Norb'): %s",
                          ctmp2);
                  info = ReadDefFileError(defname);
                } else if (bufInt[IdxNx] <= 0 || bufInt[IdxNy] <= 0 ||
                           bufInt[IdxNz] <= 0 || bufInt[IdxNorb] <= 0) {
                  fprintf(stderr,
                          "Error: lattice dimensions must be positive (got Nx=%d Ny=%d Nz=%d Norb=%d).\n",
                          bufInt[IdxNx], bufInt[IdxNy],
                          bufInt[IdxNz], bufInt[IdxNorb]);
                  info = ReadDefFileError(defname);
                }
              }
            }
            break;

          case KWTwist:
            cerr = ReadBuffInt(fp, &bufInt[IdxNTwist]);
            if (cerr != NULL && bufInt[IdxNTwist] < 0) {
              fprintf(stderr, "Error: NTwist must be non-negative (got %d).\n",
                      bufInt[IdxNTwist]);
              info = ReadDefFileError(defname);
            }
            break;


          case KWInterAll:
            cerr = ReadBuffInt(fp, &bufInt[IdxNInterAll]);
            break;

          case KWLsTrans:
            cerr = ReadBuffInt(fp, &bufInt[IdxNLsTrans]);
            break;

          case KWLsInterAll:
            cerr = ReadBuffInt(fp, &bufInt[IdxNLsInterAll]);
            break;

          case KWNBodyInterAll:
            cerr = "";
            if (ReadBuffNBodyInterAll(fp,
                                      &bufInt[IdxNNBodyInterAll],
                                      &bufInt[IdxNBodyInterAllTotalFactors],
                                      &bufInt[IdxNBodyInterAllMaxN],
                                      bufInt[IdxNsite],
                                      defname) != 0) {
              info = ReadDefFileError(defname);
            }
            break;

          case KWOptTrans:
            bufInt[IdxNQPOptTrans] = 1;
            if (FlagOptTrans > 0) {
              cerr = "";
              if (ReadPositiveProjectionCount(fp,
                                              &bufInt[IdxNQPOptTrans],
                                              "NQPOptTrans", defname) != 0) {
                info = ReadDefFileError(defname);
              }
            }
            break;

          case KWBFRange:
            cerr = fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            if (cerr != NULL) {
              cerr = fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
              if (cerr != NULL) {
                sscanf(ctmp2,"%s %d %d\n", ctmp, &bufInt[IdxNrange], &bufInt[IdxNNz]);
              }
            }
            break;

          case KWBF:
            cerr = ReadBuffInt(fp, &bufInt[IdxNBF]);
            break;

          case KWAnomalousTerm:
            cerr = "";
            if (ReadAnomalousCount(fp, &bufInt[IdxNAnomalousTerm],
                                   "NAnomalousTerm", defname) != 0) {
              info = ReadDefFileError(defname);
            }
            break;

          case KWAnomalousG:
            cerr = "";
            if (ReadAnomalousCount(fp, &bufInt[IdxNAnomalousG],
                                   "NAnomalousG", defname) != 0) {
              info = ReadDefFileError(defname);
            }
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

    {
      const long long rawNMPTrans = (long long)bufInt[IdxMPTrans];
      const long long absNMPTrans =
          rawNMPTrans < 0 ? -rawNMPTrans : rawNMPTrans;
      if (bufInt[IdxNsite] <= 0) {
        fprintf(stderr, "Error: Nsite must be positive (got %d).\n",
                bufInt[IdxNsite]);
        info = 1;
      }
      if (rawNMPTrans == 0 || rawNMPTrans == INT_MIN ||
          absNMPTrans > INT_MAX) {
        fprintf(stderr,
                "Error: NMPTrans must be a nonzero integer with "
                "abs(NMPTrans) <= %d (got %lld).\n",
                INT_MAX, rawNMPTrans);
        info = 1;
      }
      if (bufInt[IdxNQPTrans] <= 0 ||
          (rawNMPTrans != INT_MIN && absNMPTrans > bufInt[IdxNQPTrans])) {
        fprintf(stderr,
                "Error: translation projection requires NQPTrans > 0 and "
                "1 <= abs(NMPTrans) <= NQPTrans "
                "(NMPTrans=%lld, NQPTrans=%d).\n",
                rawNMPTrans, bufInt[IdxNQPTrans]);
        info = 1;
      }
      if (bufInt[IdxNsite] > 0 && bufInt[IdxNQPTrans] > 0 &&
          CheckProjectionAllocationSizes(bufInt[IdxNsite],
                                         bufInt[IdxNQPTrans]) != 0) {
        info = 1;
      }
    }

    iret = JudgeOrbitalMode(&iFlgOrbitalGeneral, iFlgOrbitalAntiParallel, iFlgOrbitalParallel);
    if (iret < 0) info = iret;

    //TODO: LanczosMode is not supported for Sz not conserved mode.

    //For indirect calculation of Green's function
    if (bufInt[IdxLanczosEstimatorMode] == 0 &&
        (bufInt[IdxNTwoBodyGEx] > 0 || bufInt[IdxLanczosMode] > 1)) {
      //Get info of CisAjs and CisAjsCktAlt(GreenTwoEx as if it's DC)
      int i;
      iOneBodyGIdx = malloc(sizeof(int *) * (2 * bufInt[IdxNsite])); //For spin
      for (i = 0; i < 2 * bufInt[IdxNsite]; i++) {
        iOneBodyGIdx[i] = malloc(sizeof(int) * (2 * bufInt[IdxNsite]));
      }
      bufInt[IdxNOneBodyG] = CountOneBodyGForLanczos(xNameListFile,
                                                     bufInt[IdxNOneBodyG], bufInt[IdxNTwoBodyGEx],
                                                     bufInt[IdxNsite], iOneBodyGIdx);
    }

    //CalcNCond
    if (bufInt[IdxNGrandCanonical] == 0) {
      if (bufInt[IdxNCond] != -1) {
        if (bufInt[IdxNCond] % 2 != 0) {
          fprintf(stderr, "Error: NCond (in modpara.def) must be even number.\n");
          info = 1;
        } else {
          bufInt[IdxNe] =
              (bufInt[IdxNLocSpin] + bufInt[IdxNCond]) / 2;
        }
      }
    } else if (bufInt[IdxNCond] != -1) {
      fprintf(stdout,
              "  Note: NCond is ignored when NGrandCanonical=1 "
              "(use NGCInitNelec).\n");
    }

    updateWeightError[0] = '\0';
    if (ValidateUpdateWeightContract(
            FlagUpdateWeight, bufInt[IdxExUpdatePath], bufInt[Idx2Sz],
            bufInt[IdxNsite], bufInt[IdxNLocSpin], bufInt[IdxNe],
            iFlgOrbitalGeneral, UpdateWeights, updateWeightError,
            sizeof(updateWeightError)) != 0) {
      fprintf(stderr, "Error: %s.\n", updateWeightError);
      info = 1;
    } else if (FlagUpdateWeight) {
      fprintf(stdout,
              "  Normalized update weights: Exchange=%.17g "
              "LocalSpinFlip=%.17g PairSpinFlip=%.17g.\n",
              UpdateWeights[UpdateWeightExchange],
              UpdateWeights[UpdateWeightLocalSpinFlip],
              UpdateWeights[UpdateWeightPairSpinFlip]);
    }

    if (bufInt[IdxNsite] > INT_MAX / 2) {
      fprintf(stderr,
              "Error: Nsite is too large: 2*Nsite would overflow int.\n");
      info = 1;
    }
    if (bufInt[IdxNe] > INT_MAX / 2) {
      fprintf(stderr,
              "Error: Ne is too large: 2*Ne would overflow int.\n");
      info = 1;
    }

    if (bufInt[IdxNGrandCanonical] == 0 &&
        (hasAnomalousTerm || hasAnomalousG)) {
      fprintf(stderr,
              "Error: AnomalousTerm/AnomalousG requires NGrandCanonical=1.\n");
      info = 1;
    }
    if (hasAnomalousG && bufInt[IdxVMCCalcMode] != 1) {
      fprintf(stderr, "Error: AnomalousG requires NVMCCalMode=1.\n");
      info = 1;
    }
    if (bufInt[IdxNAnomalousTerm] % 2 != 0) {
      fprintf(stderr,
              "Error: AnomalousTerm terms must appear as adjacent Hermite "
              "pairs.\n");
      info = 1;
    }

    if (bufInt[IdxNGrandCanonical] != 0) {
      fprintf(stdout,
              "  Note: Ne/Nelectron and NCond do not constrain particle "
              "number when NGrandCanonical=1.\n");
      if (bufInt[IdxNGrandCanonical] != 1) {
        fprintf(stderr,
                "Error: NGrandCanonical (in modpara.def) must be 0 or 1.\n");
        info = 1;
      }
      if (iFlgOrbitalGeneral != 1) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 requires OrbitalGeneral.\n");
        info = 1;
      }
      if (bufInt[IdxNLocSpin] > 0) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 does not support LocSpin "
                "(pair add/remove breaks local-spin occupancy).\n");
        info = 1;
      }
      if (bufInt[Idx2Sz] != -1) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 requires 2Sz=-1 "
                "(Sz is not conserved by pair add/remove).\n");
        info = 1;
      }
      if (bufInt[IdxExUpdatePath] != 0) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 requires NExUpdatePath=0 "
                "(GC sampler has its own move classes).\n");
        info = 1;
      }
      if (FlagUpdateWeight) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 does not support the "
                "UpdateWeight input.\n");
        info = 1;
      }
      if (NSRCG != 0) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 supports dense SR only "
                "(NSRCG must be 0).\n");
        info = 1;
      }
      if (bufInt[IdxLanczosMode] > 0) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 does not support Lanczos.\n");
        info = 1;
      }
      if (FlagRBM > 0) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 does not support RBM.\n");
        info = 1;
      }
      if (bufInt[IdxNBF] > 0) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 does not support BackFlow.\n");
        info = 1;
      }
      if (bufInt[IdxSPGaussLeg] != 1 || bufInt[IdxMPTrans] != 1 ||
          bufInt[IdxNQPOptTrans] > 1) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 initially requires "
                "NSPGaussLeg=1, NMPTrans=1, NQPOptTrans<=1.\n");
        info = 1;
      }
      if (FlagOptTrans > 0) {
        fprintf(stderr,
                "Error: NGrandCanonical=1 does not support OptTrans.\n");
        info = 1;
      }
      if (bufInt[IdxNGCInitNelec] != -1 &&
          (bufInt[IdxNGCInitNelec] < 0 ||
           (long long)bufInt[IdxNGCInitNelec] >
               2LL * (long long)bufInt[IdxNsite])) {
        fprintf(stderr,
                "Error: NGCInitNelec must satisfy "
                "0 <= NGCInitNelec <= 2*Nsite.\n");
        info = 1;
      }
    }

    if (bufInt[IdxExUpdatePath] == 4 || bufInt[IdxExUpdatePath] == 5) {
      if (bufInt[IdxNBF] > 0) {
        fprintf(stderr, "Error: NExUpdatePath=4 or 5 (t-J update) does not support BackFlow.\n");
        info = 1;
      }
      if (bufInt[IdxNLocSpin] > 0) {
        fprintf(stderr, "Error: NExUpdatePath=4 or 5 (t-J update) does not support LocSpin.\n");
        info = 1;
      }
      if (bufInt[IdxExUpdatePath] == 4 && 2LL * bufInt[IdxNe] >= bufInt[IdxNsite]) {
        fprintf(stderr,
                "Error: NExUpdatePath=4 (t-J spin hopping) requires 2*Ne < Nsite to keep at least one empty site.\n");
        info = 1;
      }
      if (bufInt[IdxExUpdatePath] == 5 && 2LL * bufInt[IdxNe] > bufInt[IdxNsite]) {
        fprintf(stderr,
                "Error: NExUpdatePath=5 (t-J update) requires 2*Ne <= Nsite to avoid double occupancy.\n");
        info = 1;
      }
      if (FlagRBM > 0) {
        fprintf(stderr, "Error: NExUpdatePath=4 or 5 (t-J update) does not support RBM.\n");
        info = 1;
      }
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

    if (bufInt[IdxLanczosMode] < 0) {
      fprintf(stderr, "Error: NLanczosMode must be non-negative (got %d).\n",
              bufInt[IdxLanczosMode]);
      info = 1;
    }
    if (bufInt[IdxLanczosEstimatorMode] != 0 &&
        bufInt[IdxLanczosEstimatorMode] != 1) {
      fprintf(stderr,
              "Error: NLanczosEstimatorMode must be 0 or 1 (got %d).\n",
              bufInt[IdxLanczosEstimatorMode]);
      info = 1;
    }
    if (bufInt[IdxLanczosMode] > 0 &&
        bufInt[IdxLanczosEstimatorMode] == 1 &&
        bufInt[IdxLanczosMode] != 1) {
      fprintf(stderr,
              "Error: stabilized power-Lanczos currently supports "
              "NLanczosMode=1 (energy/variance) only.\n");
      info = 1;
    }
    if (bufInt[IdxLanczosEstimatorMode] == 1 &&
        bufInt[IdxNTwoBodyGEx] > 0) {
      fprintf(stderr,
              "Error: stabilized power-Lanczos does not support "
              "TwoBodyGEx Green functions.\n");
      info = 1;
    }
    if (iFlgOrbitalGeneral == 1) {
      if (bufInt[IdxSPGaussLeg] > 1) {    //Check NSPGaussLeg
        if (bufInt[IdxLanczosMode] > 0) {
          fprintf(stderr,
                  "Error: FSZ Lanczos requires input SPGaussLeg <= 1 (got %d).\n",
                  bufInt[IdxSPGaussLeg]);
          info = 1;
        } else {
          fprintf(stdout, "Warning: SPGaussLeg (in modpara.def) must be 0 or 1 when orbital is general.\n");
          fprintf(stdout, "         SPGaussLeg set as 1.\n");
          bufInt[IdxSPGaussLeg] = 1;
        }
      }
      if (bufInt[IdxLanczosMode] > 1) {
        fprintf(stderr,
                "Error: FSZ supports only NLanczosMode==1 (got %d).\n",
                bufInt[IdxLanczosMode]);
        info = 1;
      } else if (bufInt[IdxLanczosMode] == 1 &&
                 bufInt[IdxVMCCalcMode] != 1) {
        fprintf(stderr,
                "Error: FSZ NLanczosMode==1 requires NVMCCalMode==1 (got %d).\n",
                bufInt[IdxVMCCalcMode]);
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

    //Check NExUpdatePath=6 (doublon-only pair hopping)
    if (bufInt[IdxExUpdatePath] == 6) {
      if (bufInt[IdxNLocSpin] > 0) {
        fprintf(stderr, "Error: NExUpdatePath=6 (doublon-only) is not compatible with NLocalSpin > 0.\n");
        info = 1;
      }
      if (FlagRBM > 0) {
        fprintf(stderr, "Error: NExUpdatePath=6 (doublon-only) is not compatible with RBM (current implementation).\n");
        info = 1;
      }
      if (bufInt[IdxNBF] > 0) {
        fprintf(stderr, "Error: NExUpdatePath=6 (doublon-only) is not compatible with NBackFlowIdx > 0.\n");
        info = 1;
      }
      if (iFlgOrbitalGeneral == 1) {
        fprintf(stderr, "Error: NExUpdatePath=6 (doublon-only) is not compatible with iFlgOrbitalGeneral=1 (FSZ).\n");
        info = 1;
      }
      if (bufInt[IdxNe] <= 0) {
        fprintf(stderr, "Error: NExUpdatePath=6 (doublon-only) requires Ne > 0.\n");
        info = 1;
      }
      if (bufInt[IdxNe] >= bufInt[IdxNsite]) {
        fprintf(stderr, "Error: NExUpdatePath=6 (doublon-only) requires Ne < Nsite (need empty sites for pair hopping).\n");
        info = 1;
      }
    }

  }//rank 0

  if (rank == 0) {
    FlagLsExplicit =
        strcmp(cFileNameListFile[KWLsTrans], "") != 0 ||
        strcmp(cFileNameListFile[KWLsInterAll], "") != 0;
    AllComplexFlag = iComplexFlgGutzwiller + iComplexFlgJastrow + iComplexFlgSpinJastrow + iComplexFlgDH2; //TBC
    AllComplexFlag += iComplexFlgDH4 + iComplexFlgOrbital;//TBC
    AllComplexFlag += iComplexFlgGeneralRBM_PhysLayer
                    + iComplexFlgGeneralRBM_HiddenLayer
                    + iComplexFlgGeneralRBM_PhysHidden
                    + iComplexFlgChargeRBM_PhysLayer
                    + iComplexFlgChargeRBM_HiddenLayer
                    + iComplexFlgChargeRBM_PhysHidden
                    + iComplexFlgSpinRBM_PhysLayer
                    + iComplexFlgSpinRBM_HiddenLayer
                    + iComplexFlgSpinRBM_PhysHidden;
    //AllComplexFlag  = 1;//DEBUG
    // AllComplexFlag= 0 -> All real, !=0 -> complex
    //if(AllComplexFlag == 0 && iFlgOrbitalGeneral == 1){
    //    fprintf(stderr, "Error: Variational parameters should be complex when orbital is general in this version.\n");
    //    info = 1;
    //}
    if(iComplexFlgOrbital > 0){
      iComplexFlgOrbital = 1;
      fprintf(stderr, "Warning: All the pairings are treated as complex variational parameters.\n");
    }
    if (bufInt[IdxNGrandCanonical] != 0 && AllComplexFlag == 0) {
      fprintf(stderr,
              "Error: NGrandCanonical=1 requires complex variational "
              "parameters.\n");
      info = 1;
    }
    if(FlagRBM != 0 && AllComplexFlag == 0){
      if (iFlgOrbitalGeneral != 0) {
        fprintf(stderr,
                "Error: real-valued RBM does not support orbital-general (FSZ) mode.\n");
        info = 1;
      }
      if (bufInt[IdxLanczosMode] > 0) {
        fprintf(stderr,
                "Error: real-valued RBM does not support Lanczos mode.\n");
        info = 1;
      }
      if (bufInt[IdxExUpdatePath] == 2) {
        fprintf(stderr,
                "Error: real-valued RBM does not support NExUpdatePath=2.\n");
        info = 1;
      }
    }
    /*
     * Single-rank runs exercise the rank-0 pre-broadcast capability gate.
     * Multi-rank runs defer the same helper until after MPI_Bcast so every
     * rank validates the broadcast scalar state.
     */
    if (size == 1
        && ValidateNBodyCapabilities(
             bufInt[IdxNBF], bufInt[IdxNNBodyG],
             bufInt[IdxNNBodyInterAll], AllComplexFlag,
             bufInt[IdxLanczosMode], 1) != 0) {
      info = 1;
    }
  }

  if (info != 0) {
    if (rank == 0) {
      fprintf(stderr, "Error: Definition files(*.def) are incomplete.\n");
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

#ifdef _mpi_use
  MPI_Bcast(bufInt, nBufInt, MPI_INT, 0, comm);
  MPI_Bcast(&FlagRBM, 1, MPI_INT, 0, comm);
  MPI_Bcast(&FlagLsExplicit, 1, MPI_INT, 0, comm);
  MPI_Bcast(&NStoreO, 1, MPI_INT, 0, comm); // for NStoreO
  MPI_Bcast(&NSRCG, 1, MPI_INT, 0, comm); // for NCG
  MPI_Bcast(&NSRCGFallback, 1, MPI_INT, 0, comm); // for SR-CG fallback
  MPI_Bcast(&NSRCGAbortOnFail, 1, MPI_INT, 0, comm); // for SR-CG failure
  MPI_Bcast(&RescaleSmat, 1, MPI_INT, 0, comm); // for Rescale S matrix
  MPI_Bcast(&useDiagScale, 1, MPI_INT, 0, comm); // for Jacobi preconditioned CG
  MPI_Bcast(&reweight, 1, MPI_INT, 0, comm); // for reweight
  MPI_Bcast(&NLanczosSupportMode, 1, MPI_INT, 0, comm);
  MPI_Bcast(&AllComplexFlag, 1, MPI_INT, 0, comm); // for Real
  MPI_Bcast(&iFlgOrbitalGeneral, 1, MPI_INT, 0, comm); // for fsz
  MPI_Bcast(&hasLattice, 1, MPI_INT, 0, comm);
  MPI_Bcast(&hasBF, 1, MPI_INT, 0, comm);
  MPI_Bcast(&hasBFRange, 1, MPI_INT, 0, comm);
  MPI_Bcast(&FlagUpdateWeight, 1, MPI_INT, 0, comm);
  MPI_Bcast(UpdateWeights, UpdateWeightCount, MPI_DOUBLE, 0, comm);
  MPI_Bcast(bufDouble, nBufDouble, MPI_DOUBLE, 0, comm);
  MPI_Bcast(CDataFileHead, nBufChar, MPI_CHAR, 0, comm);
  MPI_Bcast(CParaFileHead, nBufChar, MPI_CHAR, 0, comm);
#endif /* _mpi_use */

  Counter_max = FlagUpdateWeight ? 8 : 6;

  if (ValidateNBodyCapabilities(
        bufInt[IdxNBF], bufInt[IdxNNBodyG], bufInt[IdxNNBodyInterAll],
        AllComplexFlag, bufInt[IdxLanczosMode], rank == 0) != 0) {
    MPI_Abort(comm, EXIT_FAILURE);
  }

  NVMCCalMode = bufInt[IdxVMCCalcMode];
  NLanczosMode = bufInt[IdxLanczosMode];
  NLanczosStep = bufInt[IdxLanczosStep];
  NLanczosEstimatorMode = bufInt[IdxLanczosEstimatorMode];
  FlagGrandCanonical = bufInt[IdxNGrandCanonical];
  NGCInitNelec = bufInt[IdxNGCInitNelec];
  NAnomalousTerm = bufInt[IdxNAnomalousTerm];
  NAnomalousG = bufInt[IdxNAnomalousG];
  if (NLanczosSupportMode != 0 && NLanczosSupportMode != 1) {
    if (rank == 0) {
      fprintf(stderr, "Error: NLanczosSupportMode must be 0 or 1.\n");
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
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
  NSpinJastrowIdx = bufInt[IdxNSpinJast];
  if (NVMCCalMode == 0) {
    if (NSROptItrStep < 1) {
      if (rank == 0) {
        fprintf(stderr,
                "error: NSROptItrStep (= %d) must be a positive integer "
                "when NVMCCalMode=0.\n",
                NSROptItrStep);
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    if (NSROptItrSmp < 1 || NSROptItrSmp > NSROptItrStep) {
      if (rank == 0) {
        fprintf(stderr,
                "error: NSROptItrSmp (= %d) must satisfy "
                "1 <= NSROptItrSmp <= NSROptItrStep (= %d) "
                "when NVMCCalMode=0.\n",
                NSROptItrSmp, NSROptItrStep);
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
  }
  if (NSpinJastrowIdx < 0) {
    if (rank == 0) {
      fprintf(stderr,
              "error: NSpinJastrowIdx (= %d) must be non-negative.\n",
              NSpinJastrowIdx);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
  NDoublonHolon2siteIdx = bufInt[IdxNDH2];
  NDoublonHolon4siteIdx = bufInt[IdxNDH4];
  NOrbitalIdx = bufInt[IdxNOrbit];
  NQPTrans = bufInt[IdxNQPTrans];
  NCisAjs = bufInt[IdxNOneBodyG];
  NCisAjsCktAlt = bufInt[IdxNTwoBodyGEx];
  NCisAjsCktAltDC = bufInt[IdxNTwoBodyG];
  NNBodyG = bufInt[IdxNNBodyG];
  NBodyGTotalFactors = bufInt[IdxNBodyGTotalFactors];
  NBodyGMaxN = bufInt[IdxNBodyGMaxN];
  NInterAll = bufInt[IdxNInterAll];
  NLsTransfer = bufInt[IdxNLsTrans];
  NLsInterAll = bufInt[IdxNLsInterAll];
  NNBodyInterAll = bufInt[IdxNNBodyInterAll];
  NBodyInterAllTotalFactors = bufInt[IdxNBodyInterAllTotalFactors];
  NBodyInterAllMaxN = bufInt[IdxNBodyInterAllMaxN];
  NQPOptTrans = bufInt[IdxNQPOptTrans];
  Nrange = bufInt[IdxNrange];
  NBackFlowIdx = bufInt[IdxNBF];
  NzBF = bufInt[IdxNNz];
  NSROptCGMaxIter = bufInt[IdxSROptCGMaxIter];
  DSROptRedCut = bufDouble[IdxSROptRedCut];
  DSROptStaDel = bufDouble[IdxSROptStaDel];
  DSROptStepDt = bufDouble[IdxSROptStepDt];
  DSROptCGTol = bufDouble[IdxSROptCGTol];
  TwoSz = bufInt[Idx2Sz];

  Nx = bufInt[IdxNx];
  Ny = bufInt[IdxNy];
  Nz = bufInt[IdxNz];
  Norb = bufInt[IdxNorb];
  NTwist  = bufInt[IdxNTwist];

  if (NTwist < 0) {
    if (rank == 0) {
      fprintf(stderr, "Error: NTwist must be non-negative (got %d).\n", NTwist);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (NNBodyG < 0 || NBodyGTotalFactors < 0 || NBodyGMaxN < 0) {
    if (rank == 0) {
      fprintf(stderr,
              "Error: NBodyG sizes must be non-negative "
              "(NNBodyG=%d, NBodyGTotalFactors=%d, NBodyGMaxN=%d).\n",
              NNBodyG, NBodyGTotalFactors, NBodyGMaxN);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (NNBodyInterAll < 0 ||
      NBodyInterAllTotalFactors < 0 ||
      NBodyInterAllMaxN < 0) {
    if (rank == 0) {
      fprintf(stderr,
              "Error: NBodyInterAll sizes must be non-negative "
              "(NNBodyInterAll=%d, NBodyInterAllTotalFactors=%d, "
              "NBodyInterAllMaxN=%d).\n",
              NNBodyInterAll, NBodyInterAllTotalFactors,
              NBodyInterAllMaxN);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (NNBodyInterAll > INT_MAX / 2 ||
      NBodyInterAllTotalFactors > INT_MAX / 4) {
    if (rank == 0) {
      fprintf(stderr,
              "Error: NBodyInterAll sizes are too large "
              "(NNBodyInterAll=%d, NBodyInterAllTotalFactors=%d).\n",
              NNBodyInterAll, NBodyInterAllTotalFactors);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (BFNBodyInjectStage != BF_NBODY_INJECT_NONE
      && BFNBodyInjectTerm >= NNBodyG
      && BFNBodyInjectTerm >= NNBodyInterAll) {
    if (rank == 0) {
      fprintf(stderr,
              "Error: invalid MVMC_BF_NBODY_INJECT_TERM environment value "
              "(term=%d, NNBodyG=%d, NNBodyInterAll=%d).\n",
              BFNBodyInjectTerm, NNBodyG, NNBodyInterAll);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (NTwist > 0) {
    if (!hasLattice) {
      if (rank == 0) {
        fprintf(stderr,
                "Error: NTwist=%d requires a Lattice definition file in namelist.\n",
                NTwist);
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    /* Lattice values are bcast via bufInt, so this is rank-safe.
       The usual lattice.def form has one coordinate per site, so
       Nx*Ny*Nz*Norb should normally match Nsite. Non-matching values are
       accepted for compatibility with inputs that use the header as a
       coordinate bounding box, but warn in all cases. */
    long long lattice_prod = (long long)Nx * (long long)Ny *
                             (long long)Nz * (long long)Norb;
    if (lattice_prod != (long long)Nsite) {
      if (rank == 0) {
        fprintf(stderr,
                "warning: lattice dimensions Nx*Ny*Nz*Norb = %lld do not "
                "match Nsite = %d. The usual lattice.def form has one "
                "coordinate per site, so these values should normally match. "
                "This run will continue for compatibility; make sure the "
                "lattice.def coordinates and twist.def phases encode the "
                "intended geometry.\n",
                lattice_prod, Nsite);
      }
    }
  }

  Nneuron = bufInt[IdxNneuron];
  NneuronGeneral = bufInt[IdxNneuronGeneral];
  NneuronCharge = bufInt[IdxNneuronCharge];
  NneuronSpin = bufInt[IdxNneuronSpin];
  Nneuron += NneuronCharge+NneuronSpin+NneuronGeneral;

  NGeneralRBM_HiddenLayerIdx = bufInt[IdxNGeneralRBM_HiddenLayer];
  NGeneralRBM_PhysLayerIdx   = bufInt[IdxNGeneralRBM_PhysLayer];
  NGeneralRBM_PhysHiddenIdx  = bufInt[IdxNGeneralRBM_PhysHidden];
  NChargeRBM_HiddenLayerIdx = bufInt[IdxNChargeRBM_HiddenLayer];
  NChargeRBM_PhysLayerIdx   = bufInt[IdxNChargeRBM_PhysLayer];
  NChargeRBM_PhysHiddenIdx  = bufInt[IdxNChargeRBM_PhysHidden];
  NSpinRBM_HiddenLayerIdx = bufInt[IdxNSpinRBM_HiddenLayer];
  NSpinRBM_PhysLayerIdx   = bufInt[IdxNSpinRBM_PhysLayer];
  NSpinRBM_PhysHiddenIdx  = bufInt[IdxNSpinRBM_PhysHidden];
  NBlockSize_RBMRatio  = bufInt[IdxNBlockSize_RBMRatio];

  NRBM_PhysHiddenIdx  = NChargeRBM_PhysHiddenIdx + NSpinRBM_PhysHiddenIdx + NGeneralRBM_PhysHiddenIdx;
  NRBM_PhysLayerIdx   = NChargeRBM_PhysLayerIdx + NSpinRBM_PhysLayerIdx + NGeneralRBM_PhysLayerIdx;
  NRBM_HiddenLayerIdx = NChargeRBM_HiddenLayerIdx + NSpinRBM_HiddenLayerIdx + NGeneralRBM_HiddenLayerIdx;
  NRBM = NRBM_HiddenLayerIdx + NRBM_PhysLayerIdx + NRBM_PhysHiddenIdx;

  if (NMPTrans < 0) {
    APFlag = 1; /* anti-periodic boundary */
    NMPTrans = (int)(-(long long)NMPTrans);
  } else {
    APFlag = 0;
  }

  {
    int fszLanczosInfo = 0;
    if (rank == 0 && iFlgOrbitalGeneral != 0 && NLanczosMode == 1) {
      if (NPairHopping > 0 || NExchangeCoupling > 0 || NInterAll > 0 ||
          NNBodyInterAll > 0) {
        fprintf(stderr,
                "Error: FSZ NLanczosMode==1 currently supports only diagonal and Transfer Hamiltonian terms.\n");
        fszLanczosInfo = 1;
      }
      if (FlagRBM != 0 || NTwist > 0 || reweight == 1 ||
          FlagOptTrans > 0 || NExUpdatePath != 0) {
        fprintf(stderr,
                "Error: FSZ NLanczosMode==1 does not support RBM, Twist, "
                "reweight, OptTrans, or special update paths.\n");
        fszLanczosInfo = 1;
      }
    }
#ifdef _mpi_use
    MPI_Bcast(&fszLanczosInfo, 1, MPI_INT, 0, comm);
#endif
    if (fszLanczosInfo != 0) MPI_Abort(comm, EXIT_FAILURE);
  }

  {
    const int backflowSupported = 1;
    int bfInfo = 0;
    if (rank == 0) {
      bfInfo = BFValidateSettings(hasBF, hasBFRange, backflowSupported);
    }
#ifdef _mpi_use
    MPI_Bcast(&bfInfo, 1, MPI_INT, 0, comm);
#endif
    if (bfInfo != 0) {
      MPI_Abort(comm, EXIT_FAILURE);
    }
  }

  if (NSRCG == 2){
    useDiagScale = 1;
    NSRCG = 1;
  //  if (rank == 0) printf("remark: use preconditioned CG (Diag Scale)\n");
  } //else {
  //  useDiagScale = 0;
  //}
  NSRCGFallback = (NSRCGFallback != 0);
  NSRCGAbortOnFail = (NSRCGAbortOnFail != 0);

  if (NSplitSize < 1) {
    if (rank == 0) {
      fprintf(stderr,
              "error: NSplitSize (= %d) must be a positive integer.\n",
              NSplitSize);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }

  if (NSplitSize > 1 && NSRCG != 0) {
    if (rank == 0) {
      fprintf(stderr,
              "error: NSplitSize > 1 cannot be combined with NSRCG != 0 "
              "because the SR-CG stored-O matvec path is not supported for "
              "inner MPI splitting (NSplitSize=%d, NStore=%d, NSRCG=%d).\n",
              NSplitSize, NStoreO, NSRCG);
      fprintf(stderr, "       Use NSplitSize=1, or set NSRCG=0.\n");
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }

  if (useDiagScale){
    if(NSRCG == 1){
      if (rank == 0) printf("remark: use preconditioned CG (Diag Scale)\n");
    }else{
      if (rank == 0) printf("remark: not use preconditioned CG (Diag Scale) because NSRCG=%d != 1. Use direct method instead.\n",NSRCG);
      useDiagScale = 0;
    }
  }

  if (RescaleSmat){
    if (NSRCG != 1) {
      if (rank == 0) {
        fprintf(stderr,
                "error: invalid RescaleSmat configuration. "
                "RescaleSmat=1 requires NSRCG==1, but got NSRCG=%d.\n",
                NSRCG);
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    if (rank == 0) printf("remark: rescale S matrix \n");
  }

  if (DSROptStepDt < 0) {
    SRFlag = 1; /* diagonalization */
    if (rank == 0) fprintf(stderr, "remark: Diagonalization Mode\n");
    DSROptStepDt *= -1;
  } else {
    SRFlag = 0;
  }

  if (SRFlag != 0 && NSRCG != 0) {
    if (rank == 0) {
      fprintf(stderr,
              "error: SR-CG (NSRCG=%d) cannot be combined with diagonalization mode "
              "(DSROptStepDt < 0). Use either DSROptStepDt > 0 with SR-CG, "
              "or DSROptStepDt < 0 with NSRCG=0.\n",
              NSRCG);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }

  Nsize = 2 * Ne;
  Nsite2 = 2 * Nsite;
  NsizeMax = (FlagGrandCanonical != 0) ? Nsite2 : Nsize;
  if (FlagGrandCanonical != 0) {
    int ncur0 = (NGCInitNelec >= 0) ? NGCInitNelec : Nsite;
    if (ncur0 % 2 != 0) ncur0 -= 1;
    if (ncur0 < 0) ncur0 = 0;
    if (ncur0 > Nsite2) ncur0 = Nsite2 - (Nsite2 % 2);
    Ncur = ncur0;
    if (rank == 0) {
      fprintf(stdout, "  Note: initial GC particle number is %d.\n", Ncur);
    }
  }
  NSlater = NOrbitalIdx;
  NProj = NGutzwillerIdx + NJastrowIdx + NSpinJastrowIdx
          + 2 * 3 * NDoublonHolon2siteIdx
          + 2 * 5 * NDoublonHolon4siteIdx;
  NOptTrans = (FlagOptTrans > 0) ? NQPOptTrans : 0;

  if (BFComputeSizes(NBackFlowIdx, Nrange, NzBF, &NrangeIdx, &NBFIdxTotal, &NProjBF) != 0) {
    MPI_Abort(comm, EXIT_FAILURE);
  }

  if (FlagLsExplicit) {
    int independentInfo = 0;
    if (rank == 0 &&
        (NLanczosMode != 1 || NLanczosEstimatorMode != 1 ||
         NLanczosStep != 1 || NVMCCalMode != 1 || AllComplexFlag != 0 ||
         iFlgOrbitalGeneral != 0 || NProjBF != 0 || FlagRBM != 0 ||
         reweight != 0 || NExUpdatePath != 0 || NPairHopping != 0 ||
         NExchangeCoupling != 0 || NNBodyInterAll != 0 || NNBodyG != 0 ||
         NCisAjs != 0 || NCisAjsCktAlt != 0 || NCisAjsCktAltDC != 0 ||
         NCoulombIntra != 0 || NCoulombInter != 0 ||
         NHundCoupling != 0 ||
         (NLsTransfer <= 0 && NLsInterAll <= 0) ||
         (NTransfer <= 0 && NInterAll <= 0))) {
      fprintf(stderr,
              "Error: LsTrans/LsInterAll requires real physical mode, "
              "NLanczosMode=1, NLanczosEstimatorMode=1, "
              "NLanczosStep=1, NVMCCalMode=1, NExUpdatePath=0, and "
              "no Green functions, with nonempty H/H' containing only "
              "Trans+InterAll.\n");
      independentInfo = 1;
    }
#ifdef _mpi_use
    MPI_Bcast(&independentInfo, 1, MPI_INT, 0, comm);
#endif
    if (independentInfo != 0) MPI_Abort(comm, EXIT_FAILURE);
  }

  {
    int lanczos2NQPFull = 0;
    const long long lanczos2MPTrans =
        NMPTrans < 0 ? -(long long)NMPTrans : (long long)NMPTrans;
    if (NSPGaussLeg > 0 && lanczos2MPTrans > 0 &&
        lanczos2MPTrans <= INT_MAX && NQPOptTrans > 0 &&
        NSPGaussLeg <= INT_MAX / (int)lanczos2MPTrans &&
        NSPGaussLeg * (int)lanczos2MPTrans <=
            INT_MAX / NQPOptTrans) {
      lanczos2NQPFull =
          NSPGaussLeg * (int)lanczos2MPTrans * NQPOptTrans;
    }
    const Lanczos2Contract lanczos2Contract = {
        .step = NLanczosStep,
        .lanczosMode = NLanczosMode,
        .vmcCalMode = NVMCCalMode,
        .orbitalGeneral = iFlgOrbitalGeneral,
        .nProjBF = NProjBF,
        .flagRBM = FlagRBM,
        .reweight = reweight,
        .exUpdatePath = NExUpdatePath,
        .nPairHopping = NPairHopping,
        .nExchangeCoupling = NExchangeCoupling,
        .nInterAll = NInterAll,
        .nNBodyInterAll = NNBodyInterAll,
        .nNBodyG = NNBodyG,
        .nSpinFlipTransfer = 0,
        .nLocSpn = NLocSpn,
        .nsite = Nsite,
        .ne = Ne,
        .nTransfer = NTransfer,
        .nQPFull = lanczos2NQPFull};
    int lanczos2StatusCode = LANCZOS2_CONTRACT_OK;
    Lanczos2ContractStatus lanczos2Status;
    if (rank == 0) {
      lanczos2StatusCode =
          (int)ValidateLanczos2Contract(&lanczos2Contract);
    }
#ifdef _mpi_use
    MPI_Bcast(&lanczos2StatusCode, 1, MPI_INT, 0, comm);
#endif
    lanczos2Status = (Lanczos2ContractStatus)lanczos2StatusCode;
    if (lanczos2Status != LANCZOS2_CONTRACT_OK) {
      if (rank == 0) {
        fprintf(stderr, "Error: %s.\n",
                Lanczos2ContractError(lanczos2Status));
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
  }

  /* BFValidateSettings has already rejected unsupported BackFlow Twist and
     reweight inputs before the derived sizes are finalized here. */

  NPara = NProj + NSlater + NOptTrans + NProjBF + NRBM * FlagRBM;
  if (CheckedIntProduct3(NSPGaussLeg, NMPTrans, 1, "NQPFix",
                         &NQPFix) != 0 ||
      CheckedIntProduct3(NQPFix, NQPOptTrans, 1, "NQPFull",
                         &NQPFull) != 0) {
    MPI_Abort(comm, EXIT_FAILURE);
  }
  SROptSize = NPara + 1;

  {
    long long totalInt = 0;
    long long totalDouble = 0;
    const int orbitalDim = iFlgOrbitalGeneral == 0 ? Nsite : Nsite2;
    int sizeInfo = 0;

    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, 1, 1, "LocSpn");
    sizeInfo |= CheckedCountTerm(&totalInt, NTransfer, 4, 1, "Transfer");
    sizeInfo |= CheckedCountTerm(&totalInt, NCoulombIntra, 1, 1, "CoulombIntra");
    sizeInfo |= CheckedCountTerm(&totalInt, NCoulombInter, 2, 1, "CoulombInter");
    sizeInfo |= CheckedCountTerm(&totalInt, NHundCoupling, 2, 1, "HundCoupling");
    sizeInfo |= CheckedCountTerm(&totalInt, NPairHopping, 2, 1, "PairHopping");
    sizeInfo |= CheckedCountTerm(&totalInt, NExchangeCoupling, 2, 1, "ExchangeCoupling");
    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, 1, 1, "GutzwillerIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, Nsite, 1, "JastrowIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, Nsite, 1, "SpinJastrowIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, NDoublonHolon2siteIdx, 2,
                                 "DoublonHolon2siteIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, NDoublonHolon4siteIdx, 4,
                                 "DoublonHolon4siteIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, NQPTrans, 3,
                                 "QPTrans/QPTransInv/QPTransSgn");
    sizeInfo |= CheckedCountTerm(&totalInt, NCisAjs, 4, 1, "CisAjs");
    sizeInfo |= CheckedCountTerm(&totalInt, NCisAjsCktAlt, 2, 1, "CisAjsCktAlt");
    sizeInfo |= CheckedCountTerm(&totalInt, NCisAjsCktAltDC, 8, 1, "CisAjsCktAltDC");
    sizeInfo |= CheckedCountTerm(&totalInt, NNBodyG, 2, 1, "NBodyG metadata");
    sizeInfo |= CheckedCountTerm(&totalInt, NBodyGTotalFactors, 4, 1, "NBodyG factors");
    sizeInfo |= CheckedCountTerm(&totalInt, NInterAll, 8, 1, "InterAll");
    sizeInfo |= CheckedCountTerm(&totalInt, NLsTransfer, 4, 1, "LsTransfer");
    sizeInfo |= CheckedCountTerm(&totalInt, NLsInterAll, 8, 1, "LsInterAll");
    sizeInfo |= CheckedCountTerm(&totalInt, NNBodyInterAll, 2, 1,
                                 "NBodyInterAll metadata");
    sizeInfo |= CheckedCountTerm(&totalInt, NBodyInterAllTotalFactors, 4, 1,
                                 "NBodyInterAll factors");
    sizeInfo |= CheckedCountTerm(&totalInt, NAnomalousTerm, 5, 1,
                                 "AnomalousTerm");
    sizeInfo |= CheckedCountTerm(&totalInt, NAnomalousG, 5, 1,
                                 "AnomalousG");
    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, NQPOptTrans, 2,
                                 "QPOptTrans/QPOptTransSgn");
    sizeInfo |= CheckedCountTerm(&totalInt, Nsite, 4, 1, "LatticeIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, NTwist, Nsite, 4, "TwistIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, NneuronCharge, 1,
                                 "ChargeRBM_HiddenLayerIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, Nsite, 1,
                                 "ChargeRBM_PhysLayerIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, Nsite, NneuronCharge,
                                 "ChargeRBM_PhysHiddenIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, NneuronSpin, 1,
                                 "SpinRBM_HiddenLayerIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, Nsite, 1,
                                 "SpinRBM_PhysLayerIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, Nsite, NneuronSpin,
                                 "SpinRBM_PhysHiddenIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, NneuronGeneral, 1,
                                 "GeneralRBM_HiddenLayerIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, Nsite2, 1,
                                 "GeneralRBM_PhysLayerIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, FlagRBM, Nsite2, NneuronGeneral,
                                 "GeneralRBM_PhysHiddenIdx");
    sizeInfo |= CheckedCountTerm(&totalInt, NPara, 2, 1, "OptFlag");
    sizeInfo |= CheckedCountTerm(&totalInt, orbitalDim, orbitalDim, 2,
                                 "OrbitalIdx/OrbitalSgn");
    if (NBackFlowIdx > 0) {
      sizeInfo |= CheckedCountTerm(&totalInt, Nsite, Nrange, 1,
                                   "PosBF");
      sizeInfo |= CheckedCountTerm(&totalInt, Nsite, Nsite, 1,
                                   "RangeIdx");
    }

    sizeInfo |= CheckedCountTerm(&totalDouble, NCoulombIntra, 1, 1,
                                 "ParaCoulombIntra");
    sizeInfo |= CheckedCountTerm(&totalDouble, NCoulombInter, 1, 1,
                                 "ParaCoulombInter");
    sizeInfo |= CheckedCountTerm(&totalDouble, NHundCoupling, 1, 1,
                                 "ParaHundCoupling");
    sizeInfo |= CheckedCountTerm(&totalDouble, NPairHopping, 1, 1,
                                 "ParaPairHopping");
    sizeInfo |= CheckedCountTerm(&totalDouble, NExchangeCoupling, 1, 1,
                                 "ParaExchangeCoupling");
    sizeInfo |= CheckedCountTerm(&totalDouble, NTwist, Nsite, 6,
                                 "ParaTwist");
    sizeInfo |= CheckedCountTerm(&totalDouble, NQPOptTrans, 1, 1,
                                 "ParaQPOptTrans");

    if (sizeInfo != 0 || totalInt > INT_MAX || totalDouble > INT_MAX ||
        (size_t)totalInt > SIZE_MAX / sizeof(int) ||
        (size_t)totalDouble > SIZE_MAX / sizeof(double)) {
      if (sizeInfo == 0) {
        fprintf(stderr,
                "Error: definition table allocation size exceeds the supported range.\n");
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    NTotalDefInt = (int)totalInt;
    NTotalDefDouble = (int)totalDouble;
  }

  return 0;
}

int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm) {
  FILE *fp;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  int iKWidx = 0;
  int i, info = 0;
  int fidx = 0; /* index for OptFlag */
  int count_idx = 0;
  int bfRangeLoaded = 0;
  int rank;
  char updateWeightError[512];

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
      switch (iKWidx) {
        case KWInUpdateWeight:
          /* Parsed and normalized during ReadDefFileNInt(). */
          break;

        case KWLocSpin: /* Read locspn.def----------------------------------------*/
          if (GetLocSpinInfo(fp, LocSpn, Nsite, NLocSpn, defname) != 0) info = 1;
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

        case KWSpinJastrow:
          fidx = NGutzwillerIdx + NJastrowIdx;
          if (GetInfoJastrow(fp, SpinJastrowIdx, OptFlag, iComplexFlgSpinJastrow, &count_idx, fidx, Nsite, NSpinJastrowIdx,
                             defname) != 0)
            info = 1;
          break;

        case KWDH2:
          /*doublonholon2siteidx.def--------------------------*/
          fidx = NGutzwillerIdx + NJastrowIdx + NSpinJastrowIdx;
          if (GetInfoDH2(fp, DoublonHolon2siteIdx, OptFlag, iComplexFlgDH2, &count_idx, fidx, Nsite,
                         NDoublonHolon2siteIdx, defname) != 0)
            info = 1;
          break;

        case KWDH4:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NGutzwillerIdx + NJastrowIdx + NSpinJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          if (GetInfoDH4(fp, DoublonHolon4siteIdx, OptFlag, iComplexFlgDH4, &count_idx, fidx, Nsite,
                         NDoublonHolon4siteIdx, defname) != 0)
            info = 1;
          break;

//RBM
        case KWChargeRBM_PhysLayer:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj;
          if (GetInfoRBM_Layer(fp, ChargeRBM_PhysLayerIdx, OptFlag, iComplexFlgChargeRBM_PhysLayer,
                         &count_idx, fidx, Nsite, NChargeRBM_PhysLayerIdx, defname) != 0)
            info = 1;
          break;

        case KWSpinRBM_PhysLayer:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj+NChargeRBM_PhysLayerIdx;
          if (GetInfoRBM_Layer(fp, SpinRBM_PhysLayerIdx, OptFlag, iComplexFlgSpinRBM_PhysLayer,
                         &count_idx, fidx, Nsite, NSpinRBM_PhysLayerIdx, defname) != 0)
            info = 1;
          break;

        case KWGeneralRBM_PhysLayer:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj+NChargeRBM_PhysLayerIdx + NSpinRBM_PhysLayerIdx;
          if (GetInfoGeneralRBM_Layer(fp, GeneralRBM_PhysLayerIdx, OptFlag, iComplexFlgGeneralRBM_PhysLayer,
                         &count_idx, fidx, Nsite, NGeneralRBM_PhysLayerIdx, defname) != 0)
            info = 1;
          break;


        case KWChargeRBM_HiddenLayer:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj+NRBM_PhysLayerIdx;
          if (GetInfoRBM_Layer(fp, ChargeRBM_HiddenLayerIdx, OptFlag, iComplexFlgChargeRBM_HiddenLayer,
                         &count_idx, fidx, NneuronCharge, NChargeRBM_HiddenLayerIdx, defname) != 0)
            info = 1;
          break;

        case KWSpinRBM_HiddenLayer:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj+NRBM_PhysLayerIdx+NChargeRBM_HiddenLayerIdx;
          if (GetInfoRBM_Layer(fp, SpinRBM_HiddenLayerIdx, OptFlag, iComplexFlgSpinRBM_HiddenLayer,
                         &count_idx, fidx, NneuronSpin, NSpinRBM_HiddenLayerIdx, defname) != 0)
            info = 1;
          break;

        case KWGeneralRBM_HiddenLayer:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj+NRBM_PhysLayerIdx+NChargeRBM_HiddenLayerIdx+NSpinRBM_HiddenLayerIdx;
          if (GetInfoRBM_Layer(fp, GeneralRBM_HiddenLayerIdx, OptFlag, iComplexFlgGeneralRBM_HiddenLayer,
                         &count_idx, fidx, NneuronGeneral, NGeneralRBM_HiddenLayerIdx, defname) != 0)
            info = 1;
          break;


        case KWChargeRBM_PhysHidden:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj+NRBM_PhysLayerIdx+NRBM_HiddenLayerIdx;
          if (GetInfoRBM_PhysHidden(fp, ChargeRBM_PhysHiddenIdx, OptFlag, iComplexFlgChargeRBM_PhysHidden,
                         &count_idx, fidx, Nsite, NneuronCharge, NChargeRBM_PhysHiddenIdx, defname) != 0)
            info = 1;
          break;

        case KWSpinRBM_PhysHidden:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj+NRBM_PhysLayerIdx+NRBM_HiddenLayerIdx+NChargeRBM_PhysHiddenIdx;
          if (GetInfoRBM_PhysHidden(fp, SpinRBM_PhysHiddenIdx, OptFlag, iComplexFlgSpinRBM_PhysHidden,
                         &count_idx, fidx, Nsite, NneuronSpin, NSpinRBM_PhysHiddenIdx, defname) != 0)
            info = 1;
          break;

        case KWGeneralRBM_PhysHidden:
          /*doublonholon4siteidx.def--------------------------*/
          fidx = NProj+NRBM_PhysLayerIdx+NRBM_HiddenLayerIdx+NChargeRBM_PhysHiddenIdx+NSpinRBM_PhysHiddenIdx;
          if (GetInfoGeneralRBM_PhysHidden(fp, GeneralRBM_PhysHiddenIdx, OptFlag, iComplexFlgGeneralRBM_PhysHidden,
                         &count_idx, fidx, Nsite, NneuronGeneral, NGeneralRBM_PhysHiddenIdx, defname) != 0)
            info = 1;
          break;
//RBM
        case KWOrbital:
        case KWOrbitalAntiParallel:
          /*orbitalidxs.def------------------------------------*/
          fidx = NProj + FlagRBM * NRBM + NProjBF;
          if (GetInfoOrbitalAntiParallel(fp, OrbitalIdx, OptFlag, OrbitalSgn, &count_idx,
                                         fidx, iComplexFlgOrbital, iFlgOrbitalGeneral, APFlag, Nsite, iNOrbitalAntiParallel,
                                         defname) != 0)
            info = 1;
          break;

        case KWOrbitalGeneral:
          fidx = NProj + FlagRBM * NRBM + NProjBF;
          if (GetInfoOrbitalGeneral(fp, OrbitalIdx, OptFlag, OrbitalSgn, &count_idx,
                                    fidx, iComplexFlgOrbital, iFlgOrbitalGeneral, APFlag, Nsite, NOrbitalIdx,
                                    defname) != 0)
            info = 1;
          break;

        case KWOrbitalParallel:
          /*orbitalidxt.def------------------------------------*/
          fidx = NProj + FlagRBM * NRBM + NProjBF + iNOrbitalAntiParallel;
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
          if (GetInfoOneBodyG(fp, CisAjsIdx, iOneBodyGIdx, NCisAjsCktAlt>0, Nsite, NCisAjs, defname) != 0) info = 1;
          break;

        case KWTwoBodyGEx:
          /*cisajscktalt.def----------------------------------*/
          if (GetInfoTwoBodyGEx(fp, CisAjsCktAltIdx, iOneBodyGIdx, CisAjsIdx, Nsite, NCisAjsCktAlt, defname) != 0) info = 1;
          break;

        case KWTwoBodyG:
          /*cisajscktaltdc.def--------------------------------*/
          if (GetInfoTwoBodyG(fp, CisAjsCktAltDCIdx, Nsite, NCisAjsCktAltDC, defname) != 0)
            info = 1;
          break;

        case KWNBodyG:
          /*nbodyg.def----------------------------------------*/
          if (GetInfoNBodyG(fp, NBodyGN, NBodyGOffset, NBodyGIdx,
                            Nsite, NNBodyG, NBodyGTotalFactors,
                            NBodyGMaxN, iFlgOrbitalGeneral,
                            defname) != 0)
            info = 1;
          break;

        case KWLattice:
          /*lattice.def---------------------------------------*/
          if (GetInfoLattice(fp, LatticeIdx, Nsite, Nx, Ny, Nz, Norb, defname) != 0) info = 1;
          break;

        case KWTwist:
          /*twist.def---------------------------------------*/
          if (GetInfoTwist(fp, TwistIdx, ParaTwist, Nsite, NTwist, defname) != 0) info = 1;
          break;


        case KWInterAll:
          /*interall.def---------------------------------------*/
          if (GetInfoInterAll(fp, InterAll, ParaInterAll, Nsite, NInterAll, defname) != 0) info = 1;
          break;

        case KWLsTrans:
          if (GetTransferInfo(fp, LsTransfer, ParaLsTransfer, Nsite,
                              NLsTransfer, defname) != 0) info = 1;
          break;

        case KWLsInterAll:
          if (GetInfoInterAll(fp, LsInterAll, ParaLsInterAll, Nsite,
                              NLsInterAll, defname) != 0) info = 1;
          break;

        case KWNBodyInterAll:
          /*nbodyinterall.def----------------------------------*/
          if (GetInfoNBodyInterAll(fp, NBodyInterAllN,
                                   NBodyInterAllOffset,
                                   NBodyInterAllIdx,
                                   ParaNBodyInterAll,
                                   Nsite, NNBodyInterAll,
                                   NBodyInterAllTotalFactors,
                                   NBodyInterAllMaxN,
                                   iFlgOrbitalGeneral,
                                   defname) != 0)
            info = 1;
          break;

        case KWAnomalousTerm:
          if (GetInfoAnomalousTerm(fp, AnomalousTerm, ParaAnomalousTerm,
                                   Nsite, NAnomalousTerm, defname) != 0)
            info = 1;
          break;

        case KWAnomalousG:
          if (GetInfoAnomalousG(fp, AnomalousG, Nsite, NAnomalousG,
                                defname) != 0)
            info = 1;
          break;

        case KWOptTrans:
          /*qpopttrans.def------------------------------------*/
          fidx = NProj + FlagRBM * NRBM + NProjBF + NSlater;
          if (GetInfoOptTrans(fp, QPOptTrans, ParaQPOptTrans, OptFlag, QPOptTransSgn, FlagOptTrans, &count_idx, fidx,
                              APFlag, Nsite, NQPOptTrans, defname) != 0)
            info = 1;
          break;

        case KWBFRange:
          /*rangebf.def--------------------------*/
          if (Nrange > 0 && !bfRangeLoaded) {
            if (BFReadRange(fp, defname) != 0) info = 1;
            else bfRangeLoaded = 1;
          }
          break;

        case KWBF:
          if (NBackFlowIdx > 0) {
            int bfInfo = 0;
            if (!bfRangeLoaded && Nrange > 0) {
              bfInfo = ReadBackFlowRangeDefinition(cFileNameListFile[KWBFRange]);
              if (bfInfo == 0) bfRangeLoaded = 1;
            }
            if (bfInfo == 0 && BFReadDefinition(fp, OptFlag, &count_idx, defname) != 0) bfInfo = 1;
            if (bfInfo != 0) info = 1;
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
    if (info == 0 && NProjBF > 0 && iFlgOrbitalGeneral != 0 && NVMCCalMode == 0) {
      int hasOptimizedParameter = 0;
      for (i = 0; i < 2 * NPara; i++) {
        if (OptFlag[i] == 1) {
          hasOptimizedParameter = 1;
          break;
        }
      }
      if (!hasOptimizedParameter) {
        fprintf(stderr,
                "Error: BackFlow FSZ optimization requires at least one optimized variational parameter.\n");
        info = 1;
      }
    }
    if (info == 0 && BFValidateFszDefinitionDetails() != 0) info = 1;
    fprintf(stdout, "finish reading parameters.\n");
  } /* if(rank==0) */

  /*
   * Definition arrays are read only on rank 0.  Broadcast the verdict before
   * any rank consumes them so malformed input reaches the same abort path on
   * every rank.
   */
#ifdef _mpi_use
  MPI_Bcast(&info, 1, MPI_INT, 0, comm);
#endif

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

  if (FlagLsExplicit) {
    int lsInfo = 0;
    if (rank == 0) {
      for (i = 0; i < NTransfer; ++i) {
        if (Transfer[i][1] != Transfer[i][3] ||
            !isfinite(creal(ParaTransfer[i])) ||
            cimag(ParaTransfer[i]) != 0.0) {
          fprintf(stderr,
                  "Error: physical Trans row %d must be finite, real, and "
                  "S_z-conserving with LsTrans/LsInterAll.\n", i);
          lsInfo = 1;
          break;
        }
      }
      for (i = 0; !lsInfo && i < NInterAll; ++i) {
        if (InterAll[i][1] != InterAll[i][3] ||
            InterAll[i][5] != InterAll[i][7] ||
            !isfinite(creal(ParaInterAll[i])) ||
            cimag(ParaInterAll[i]) != 0.0) {
          fprintf(stderr,
                  "Error: physical InterAll row %d must be finite, real, "
                  "and S_z-conserving with LsTrans/LsInterAll.\n", i);
          lsInfo = 1;
          break;
        }
      }
      for (i = 0; !lsInfo && i < NLsTransfer; ++i) {
        if (LsTransfer[i][1] < 0 || LsTransfer[i][1] > 1 ||
            LsTransfer[i][3] < 0 || LsTransfer[i][3] > 1 ||
            LsTransfer[i][1] != LsTransfer[i][3] ||
            !isfinite(creal(ParaLsTransfer[i])) ||
            cimag(ParaLsTransfer[i]) != 0.0) {
          fprintf(stderr,
                  "Error: LsTrans row %d must be finite, real, and "
                  "S_z-conserving.\n", i);
          lsInfo = 1;
          break;
        }
      }
      for (i = 0; !lsInfo && i < NLsInterAll; ++i) {
        if (LsInterAll[i][1] < 0 || LsInterAll[i][1] > 1 ||
            LsInterAll[i][3] < 0 || LsInterAll[i][3] > 1 ||
            LsInterAll[i][5] < 0 || LsInterAll[i][5] > 1 ||
            LsInterAll[i][7] < 0 || LsInterAll[i][7] > 1 ||
            LsInterAll[i][1] != LsInterAll[i][3] ||
            LsInterAll[i][5] != LsInterAll[i][7] ||
            !isfinite(creal(ParaLsInterAll[i])) ||
            cimag(ParaLsInterAll[i]) != 0.0) {
          fprintf(stderr,
                  "Error: LsInterAll row %d must be finite, real, and "
                  "S_z-conserving.\n", i);
          lsInfo = 1;
          break;
        }
      }
    }
#ifdef _mpi_use
    MPI_Bcast(&lsInfo, 1, MPI_INT, 0, comm);
#endif
    if (lsInfo != 0) MPI_Abort(comm, EXIT_FAILURE);
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
  SafeMpiBcast_fcmp(ParaTransfer, NTransfer + NInterAll, comm);
  if (NLsTransfer > 0) {
    SafeMpiBcastInt(LsTransfer[0], 4 * NLsTransfer, comm);
  }
  if (NLsInterAll > 0) {
    SafeMpiBcastInt(LsInterAll[0], 8 * NLsInterAll, comm);
  }
  SafeMpiBcast_fcmp(ParaLsTransfer, NLsTransfer, comm);
  SafeMpiBcast_fcmp(ParaLsInterAll, NLsInterAll, comm);
  SafeMpiBcast_fcmp(ParaNBodyInterAll, NNBodyInterAll, comm);
  SafeMpiBcast_fcmp(ParaAnomalousTerm, NAnomalousTerm, comm);
  SafeMpiBcast(ParaCoulombIntra, NTotalDefDouble, comm);
  SafeMpiBcast_fcmp(ParaQPTrans, NQPTrans, comm);
#endif /* _mpi_use */

  updateWeightError[0] = '\0';
  if (ValidateUpdateWeightLocSpin(
          FlagUpdateWeight, Nsite, NLocSpn, LocSpn, UpdateWeights,
          updateWeightError, sizeof(updateWeightError)) != 0) {
    if (rank == 0) fprintf(stderr, "Error: %s.\n", updateWeightError);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  {
    int nSpinFlipTransfer = 0;
    Lanczos2ContractStatus lanczos2Status;
    Lanczos2Contract lanczos2Contract;
    for (i = 0; i < NTransfer; i++) {
      if (Transfer[i][1] != Transfer[i][3]) nSpinFlipTransfer++;
    }
    lanczos2Contract.step = NLanczosStep;
    lanczos2Contract.lanczosMode = NLanczosMode;
    lanczos2Contract.vmcCalMode = NVMCCalMode;
    lanczos2Contract.orbitalGeneral = iFlgOrbitalGeneral;
    lanczos2Contract.nProjBF = NProjBF;
    lanczos2Contract.flagRBM = FlagRBM;
    lanczos2Contract.reweight = reweight;
    lanczos2Contract.exUpdatePath = NExUpdatePath;
    lanczos2Contract.nPairHopping = NPairHopping;
    lanczos2Contract.nExchangeCoupling = NExchangeCoupling;
    lanczos2Contract.nInterAll = NInterAll;
    lanczos2Contract.nNBodyInterAll = NNBodyInterAll;
    lanczos2Contract.nNBodyG = NNBodyG;
    lanczos2Contract.nSpinFlipTransfer = nSpinFlipTransfer;
    lanczos2Contract.nLocSpn = NLocSpn;
    lanczos2Contract.nsite = Nsite;
    lanczos2Contract.ne = Ne;
    lanczos2Contract.nTransfer = NTransfer;
    lanczos2Contract.nQPFull = NQPFull;
    lanczos2Status =
        NLanczosMode > 0 && NLanczosEstimatorMode == 1 && !FlagLsExplicit
            ? ValidateCorrectedPowerLanczosExecutionContract(
                  &lanczos2Contract)
            : ValidateLanczos2Contract(&lanczos2Contract);
    if (lanczos2Status != LANCZOS2_CONTRACT_OK) {
      if (rank == 0) {
        fprintf(stderr, "Error: %s.\n",
                Lanczos2ContractError(lanczos2Status));
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
  }

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

        case KWInSpinJastrow:
          fprintf(stdout, "Read InSpinJastrow File. \n");
          if (idx != NSpinJastrowIdx) {
            info = 1;
            break;
          }
          count = NGutzwillerIdx + NJastrowIdx;
          if (NSpinJastrowIdx > 0) {
            int *seen = (int *)calloc(NSpinJastrowIdx, sizeof(int));
            if (seen == NULL) {
              fprintf(stderr, "Error: memory allocation failed in InSpinJastrow parser.\n");
              info = 1;
              break;
            }
            for (i = 0; i < NSpinJastrowIdx; i++) {
              if (fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp) != 3) {
                fprintf(stderr, "Error in %s: invalid parameter row %d in InSpinJastrow.\n", defname, i + 1);
                info = 1;
                break;
              }
              if (idx < 0 || idx >= NSpinJastrowIdx) {
                fprintf(stderr, "Error in %s: index %d out of range [0, %d) in InSpinJastrow.\n",
                        defname, idx, NSpinJastrowIdx);
                info = 1;
                break;
              }
              if (seen[idx] != 0) {
                fprintf(stderr, "Error in %s: duplicated index %d in InSpinJastrow.\n", defname, idx);
                info = 1;
                break;
              }
              seen[idx] = 1;
              Proj[idx+count] = tmp_real + I * tmp_comp;
            }
            if (info == 0) {
              for (i = 0; i < NSpinJastrowIdx; i++) {
                if (seen[i] == 0) {
                  fprintf(stderr, "Error in %s: missing index %d in InSpinJastrow.\n", defname, i);
                  info = 1;
                  break;
                }
              }
            }
            free(seen);
            if (info != 0) break;
          }
          break;

        case KWInDH2:
          fprintf(stdout, "Read InDH2 File. \n");
          if (idx != NDoublonHolon2siteIdx) {
            info = 1;
            continue;
          }
          count = NGutzwillerIdx + NJastrowIdx + NSpinJastrowIdx;
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
          count = NGutzwillerIdx + NJastrowIdx + NSpinJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          for (i = count; i < count + 2 * 5 * NDoublonHolon4siteIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            Proj[idx+count] = tmp_real + I * tmp_comp;
          }
          break;

//RBM
        case KWInChargeRBM_PhysLayer:
          fprintf(stdout, "Read InChargeRBM_PhysLayer File. \n");
          if (idx != NChargeRBM_PhysLayerIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = 0;
          for (i = count; i < count + NChargeRBM_PhysLayerIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;


        case KWInSpinRBM_PhysLayer:
          fprintf(stdout, "Read InSpinRBM_PhysLayer File. \n");
          if (idx != NSpinRBM_PhysLayerIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = NChargeRBM_PhysLayerIdx;
          for (i = count; i < count + NSpinRBM_PhysLayerIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInGeneralRBM_PhysLayer:
          fprintf(stdout, "Read InGeneralRBM_PhysLayer File. \n");
          if (idx != NGeneralRBM_PhysLayerIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = NChargeRBM_PhysLayerIdx + NSpinRBM_PhysLayerIdx;
          for (i = count; i < count + NGeneralRBM_PhysLayerIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;



        case KWInChargeRBM_HiddenLayer:
          fprintf(stdout, "Read InChargeRBM_HiddenLayer File. \n");
          if (idx != NChargeRBM_HiddenLayerIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = NRBM_PhysLayerIdx;
          for (i = count; i < count + NChargeRBM_HiddenLayerIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInSpinRBM_HiddenLayer:
          fprintf(stdout, "Read InSpinRBM_HiddenLayer File. \n");
          if (idx != NSpinRBM_HiddenLayerIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = NRBM_PhysLayerIdx+NChargeRBM_HiddenLayerIdx;
          for (i = count; i < count + NSpinRBM_HiddenLayerIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInGeneralRBM_HiddenLayer:
          fprintf(stdout, "Read InGeneralRBM_HiddenLayer File. \n");
          if (idx != NGeneralRBM_HiddenLayerIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = NRBM_PhysLayerIdx+NChargeRBM_HiddenLayerIdx+NSpinRBM_HiddenLayerIdx;
          for (i = count; i < count + NGeneralRBM_HiddenLayerIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;


        case KWInChargeRBM_PhysHidden:
          fprintf(stdout, "Read InChargeRBM_PhysHidden File. \n");
          if (idx != NChargeRBM_PhysHiddenIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = NRBM_PhysLayerIdx + NRBM_HiddenLayerIdx;
          for (i = count; i < count + NChargeRBM_PhysHiddenIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInSpinRBM_PhysHidden:
          fprintf(stdout, "Read InSpinRBM_PhysHidden File. \n");
          if (idx != NSpinRBM_PhysHiddenIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = NRBM_PhysLayerIdx + NRBM_HiddenLayerIdx + NChargeRBM_PhysHiddenIdx;
          for (i = count; i < count + NSpinRBM_PhysHiddenIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;

        case KWInGeneralRBM_PhysHidden:
          fprintf(stdout, "Read InGeneralRBM_PhysHidden File. \n");
          if (idx != NGeneralRBM_PhysHiddenIdx) {
            info = 1;
            continue;
          }
          //count = NGutzwillerIdx + NJastrowIdx + 2 * 3 * NDoublonHolon2siteIdx;
          count = NRBM_PhysLayerIdx + NRBM_HiddenLayerIdx + NChargeRBM_PhysHiddenIdx + NSpinRBM_PhysHiddenIdx;
          for (i = count; i < count + NGeneralRBM_PhysHiddenIdx; i++) {
            fscanf(fp, "%d %lf %lf ", &idx, &tmp_real, &tmp_comp);
            RBM[idx + count] = tmp_real + I * tmp_comp;
          }
          break;
//RBM

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
  bufInt[IdxLanczosStep] = 1;
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
  bufInt[IdxNSpinJast] = 0;
  bufInt[IdxNDH2] = 0;
  bufInt[IdxNDH4] = 0;
  bufInt[IdxNOrbit] = 0;
  bufInt[IdxNQPTrans] = 0;
  bufInt[IdxNOneBodyG] = 0;
  bufInt[IdxNTwoBodyG] = 0;
  bufInt[IdxNTwoBodyGEx] = 0;
  bufInt[IdxNNBodyG] = 0;
  bufInt[IdxNBodyGTotalFactors] = 0;
  bufInt[IdxNBodyGMaxN] = 0;
  bufInt[IdxNInterAll] = 0;
  bufInt[IdxNLsTrans] = 0;
  bufInt[IdxNLsInterAll] = 0;
  bufInt[IdxNNBodyInterAll] = 0;
  bufInt[IdxNBodyInterAllTotalFactors] = 0;
  bufInt[IdxNBodyInterAllMaxN] = 0;
  bufInt[IdxNQPOptTrans] = 1;
  bufInt[IdxSROptCGMaxIter] = 0;
  bufInt[IdxNBF] = 0;
  bufInt[IdxNrange] = 0;
  bufInt[IdxNNz] = 0;
  bufInt[Idx2Sz] = -1;// -1: sz is not fixed :fsz
  bufInt[IdxNCond] = -1;
  bufInt[IdxLanczosEstimatorMode] = 0;
  bufInt[IdxNGrandCanonical] = 0;
  bufInt[IdxNGCInitNelec] = -1;
  bufInt[IdxNAnomalousTerm] = 0;
  bufInt[IdxNAnomalousG] = 0;

//RBM
  bufInt[IdxNneuron] = 0;
  bufInt[IdxNneuronGeneral] = 0;
  bufInt[IdxNneuronCharge] = 0;
  bufInt[IdxNneuronSpin] = 0;
  bufInt[IdxNChargeRBM_HiddenLayer] = 0;
  bufInt[IdxNChargeRBM_PhysLayer]   = 0;
  bufInt[IdxNChargeRBM_PhysHidden]  = 0;
  bufInt[IdxNSpinRBM_HiddenLayer] = 0;
  bufInt[IdxNSpinRBM_PhysLayer]   = 0;
  bufInt[IdxNSpinRBM_PhysHidden]  = 0;
  bufInt[IdxNGeneralRBM_HiddenLayer] = 0;
  bufInt[IdxNGeneralRBM_PhysLayer]   = 0;
  bufInt[IdxNGeneralRBM_PhysHidden]  = 0;
  bufInt[IdxNBlockSize_RBMRatio]  = 200;
//RBM

  bufDouble[IdxSROptRedCut] = 0.001;
  bufDouble[IdxSROptStaDel] = 0.02;
  bufDouble[IdxSROptStepDt] = 0.02;
  bufDouble[IdxSROptCGTol] = 1.0e-10;
  NStoreO = 1;
  NSRCG = 0;
  NSRCGFallback = 0;
  NSRCGAbortOnFail = 1;
  RescaleSmat  = 0;
  useDiagScale = 0;
  reweight = 0;
  NLanczosSupportMode = 0;

  bufInt[IdxNx] = 1;
  bufInt[IdxNy] = 1;
  bufInt[IdxNz] = 1;
  bufInt[IdxNorb] = 1;
  bufInt[IdxNTwist] = 0;
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
          char valueText[D_FileNameMax];
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (*ctmp2 == '\n' || ctmp2[0] == '-') continue;
            if (sscanf(ctmp2, "%255s %255s", ctmp, valueText) != 2) {
              fprintf(stderr, "Error: malformed ModPara line: %s", ctmp2);
              return ReadDefFileError(defname);
            }
            dtmp = strtod(valueText, NULL);
            if (CheckWords(ctmp, "NVMCCalMode") == 0) {
              bufInt[IdxVMCCalcMode] = (int) dtmp;
            } else if (CheckWords(ctmp, "NLanczosMode") == 0) {
              bufInt[IdxLanczosMode] = (int) dtmp;
            } else if (CheckWords(ctmp, "NLanczosStep") == 0) {
              /*
               * Do not silently truncate a fractional or non-finite step.
               * Store an invalid sentinel and let the rank-safe capability
               * gate issue the single canonical diagnostic.
               */
              bufInt[IdxLanczosStep] =
                  (dtmp == 1.0 || dtmp == 2.0) ? (int)dtmp : 0;
            } else if (CheckWords(ctmp, "NLanczosEstimatorMode") == 0) {
              bufInt[IdxLanczosEstimatorMode] =
                  (dtmp == 0.0 || dtmp == 1.0) ? (int)dtmp : -1;
            } else if (CheckWords(ctmp, "NGrandCanonical") == 0) {
              bufInt[IdxNGrandCanonical] = (int)dtmp;
            } else if (CheckWords(ctmp, "NGCInitNelec") == 0) {
              bufInt[IdxNGCInitNelec] = (int)dtmp;
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
            } else if (CheckWords(ctmp, "NSPGaussLeg") == 0) {
              bufInt[IdxSPGaussLeg] = (int) dtmp;
            } else if (CheckWords(ctmp, "NSPStot") == 0) {
              bufInt[IdxSPStot] = (int) dtmp;
            } else if (CheckWords(ctmp, "NMPTrans") == 0) {
              long long rawNMPTrans;
              int offset = 0;
              if (sscanf(ctmp2, "%*s %*s %n", &offset) < 0 ||
                  !LineTailIsWhitespace(ctmp2, offset) ||
                  ParseStrictLongLong(valueText, &rawNMPTrans) != 0 ||
                  rawNMPTrans < INT_MIN || rawNMPTrans > INT_MAX) {
                fprintf(stderr,
                        "Error: NMPTrans must be an exact integer in int range: %s",
                        ctmp2);
                return ReadDefFileError(defname);
              }
              bufInt[IdxMPTrans] = (int)rawNMPTrans;
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
            } else if (CheckWords(ctmp, "NSRCGFallback") == 0) {
              NSRCGFallback = (int) dtmp;
            } else if (CheckWords(ctmp, "NSRCGAbortOnFail") == 0) {
              NSRCGAbortOnFail = (int) dtmp;
            } else if (CheckWords(ctmp, "RescaleSmat") == 0) {
              RescaleSmat = (int) dtmp;
            } else if (CheckWords(ctmp, "useDiagScale") == 0) {
              useDiagScale = (int) dtmp;
            } else if (CheckWords(ctmp, "reweight") == 0) {
              reweight = (int) dtmp;
            } else if (CheckWords(ctmp, "NLanczosSupportMode") == 0) {
              NLanczosSupportMode =
                  (dtmp == 0.0 || dtmp == 1.0) ? (int)dtmp : -1;
//RBM
            } else if (CheckWords(ctmp, "Nneuron") == 0) {
              bufInt[IdxNneuron] = (int) dtmp;
            } else if (CheckWords(ctmp, "NneuronCharge") == 0) {
              bufInt[IdxNneuronCharge] = (int) dtmp;
            } else if (CheckWords(ctmp, "NneuronSpin") == 0) {
              bufInt[IdxNneuronSpin] = (int) dtmp;
            } else if (CheckWords(ctmp, "NneuronGeneral") == 0) {
              bufInt[IdxNneuronGeneral] = (int) dtmp;
            } else if (CheckWords(ctmp, "NBlockSize_RBMRatio") == 0) {
              bufInt[IdxNBlockSize_RBMRatio] = (int) dtmp;
//RBM
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
      fclose(fp);
    }
  }
  //fclose(fp);
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

int GetLocSpinInfo(FILE *fp, int *ArrayIdx, int Nsite, int NLocalSpin,
                   char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0, localSpinCount = 0;
  int site;
  int x0 = 0, x1 = 0;
  if (fp == NULL || ArrayIdx == NULL || Nsite <= 0 ||
      NLocalSpin < 0 || NLocalSpin > Nsite) {
    fprintf(stderr, "Error: Invalid LocSpin dimensions.\n");
    return 1;
  }
  for (site = 0; site < Nsite; site++) ArrayIdx[site] = -1;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    if (sscanf(ctmp2, "%d %d", &x0, &x1) != 2) {
      fprintf(stderr, "Error: Malformed LocSpin definition.\n");
      info = 1;
      break;
    }
    if (CheckSite(x0, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect.\n");
      info = 1;
      break;
    }
    if (x1 != 0 && x1 != 1) {
      fprintf(stderr,
              "Error: LocSpin flag must be 0 or 1 (site=%d, value=%d).\n",
              x0, x1);
      info = 1;
      break;
    }
    if (ArrayIdx[x0] != -1) {
      fprintf(stderr, "Error: Duplicate LocSpin site index: %d.\n", x0);
      info = 1;
      break;
    }
    ArrayIdx[x0] = x1;
    localSpinCount += x1;
    idx++;
  }
  if (info == 0) {
    for (site = 0; site < Nsite; site++) {
      if (ArrayIdx[site] == -1) {
        fprintf(stderr, "Error: Missing LocSpin site index: %d.\n", site);
        info = 1;
        break;
      }
    }
  }
  if (info == 0 && idx != Nsite) info = ReadDefFileError(defname);
  if (info == 0 && localSpinCount != NLocalSpin) {
    fprintf(stderr,
            "Error: NLocalSpin header (%d) does not match LocSpin flag "
            "count (%d).\n",
            NLocalSpin, localSpinCount);
    info = 1;
  }
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
    }else{
      ArrayOpt[2 * fidx + 1] = 0;
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
  int i = 0, j = 0, n = 0;
  int a = 0, b = 0;
  int line = 0;
  int idxExpected = Nsite * (Nsite - 1);
  int *pairClass = NULL;
  int fidx = _fidx;
  if (NArray > 0) {
    pairClass = (int *)malloc(sizeof(int) * Nsite * Nsite);
    if (pairClass == NULL) {
      fprintf(stderr, "Error: memory allocation failed in GetInfoJastrow.\n");
      return 1;
    }
    for (i = 0; i < Nsite; i++) {
      for (j = 0; j < Nsite; j++) {
        if (i == j) {
          ArrayIdx[i][j] = -1;
          pairClass[i * Nsite + j] = -1;
        } else {
          ArrayIdx[i][j] = -2;
          pairClass[i * Nsite + j] = -2;
        }
      }
    }

    while (idx0 < idxExpected) {
      if (fscanf(fp, "%d %d %d", &i, &j, &n) != 3) {
        fprintf(stderr, "Error in %s: failed to read Jastrow row %d.\n", defname, line + 1);
        info = 1;
        break;
      }
      line++;
      if (i == j) {
        fprintf(stderr, "Error in %s: [Condition] i neq j at row %d.\n", defname, line);
        info = 1;
        break;
      }
      if (CheckPairSite(i, j, Nsite) != 0) {
        fprintf(stderr, "Error in %s: Site index is incorrect at row %d.\n", defname, line);
        info = 1;
        break;
      }
      if (n < 0 || n >= NArray) {
        fprintf(stderr, "Error in %s: class index %d out of range [0, %d) at row %d.\n",
                defname, n, NArray, line);
        info = 1;
        break;
      }
      if (ArrayIdx[i][j] != -2) {
        fprintf(stderr, "Error in %s: duplicated pair (%d,%d) at row %d.\n", defname, i, j, line);
        info = 1;
        break;
      }
      ArrayIdx[i][j] = n;

      if (i < j) {
        a = i;
        b = j;
      } else {
        a = j;
        b = i;
      }
      if (pairClass[a * Nsite + b] == -2) {
        pairClass[a * Nsite + b] = n;
      } else if (pairClass[a * Nsite + b] != n) {
        fprintf(stderr, "Error in %s: inconsistent class for pair (%d,%d) at row %d.\n",
                defname, a, b, line);
        info = 1;
        break;
      }
      idx0++;
    }
    if (info == 0) {
      for (i = 0; i < Nsite; i++) {
        for (j = 0; j < Nsite; j++) {
          if (i == j) continue;
          if (ArrayIdx[i][j] == -2) {
            fprintf(stderr, "Error in %s: missing pair (%d,%d).\n", defname, i, j);
            info = 1;
            break;
          }
        }
        if (info != 0) break;
      }
    }
    if (info == 0) {
      for (i = 0; i < Nsite; i++) {
        for (j = i + 1; j < Nsite; j++) {
          if (pairClass[i * Nsite + j] == -2) {
            fprintf(stderr, "Error in %s: missing undirected pair (%d,%d).\n", defname, i, j);
            info = 1;
            break;
          }
          ArrayIdx[i][j] = pairClass[i * Nsite + j];
          ArrayIdx[j][i] = pairClass[i * Nsite + j];
        }
        if (info != 0) break;
      }
    }
    if (info == 0) {
      idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    }
    if (info != 0 || idx0 != idxExpected || idx1 != NArray) {
      info = ReadDefFileError(defname);
    }
    free(pairClass);
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

static int ReadProjectionMappings(FILE *fp, int **array, int **arraySgn,
                                  int **arrayInv, int apFlag, int nsite,
                                  int nArray, const char *defname) {
  char line[D_CharTmpReadDef];
  unsigned char *sourceSeen = NULL;
  unsigned char *destinationSeen = NULL;
  size_t expectedRows;
  size_t row = 0;
  int info = 0;

  if (fp == NULL || array == NULL || arraySgn == NULL || nsite <= 0 ||
      nArray <= 0 || (size_t)nsite > SIZE_MAX / (size_t)nArray) {
    fprintf(stderr, "Error: invalid projection mapping dimensions in %s.\n",
            defname);
    return 1;
  }
  expectedRows = (size_t)nsite * (size_t)nArray;
  sourceSeen = (unsigned char *)calloc(expectedRows, sizeof(*sourceSeen));
  destinationSeen =
      (unsigned char *)calloc(expectedRows, sizeof(*destinationSeen));
  if (sourceSeen == NULL || destinationSeen == NULL) {
    fprintf(stderr,
            "Error: failed to allocate projection parser state for %s.\n",
            defname);
    free(sourceSeen);
    free(destinationSeen);
    return 1;
  }

  while (fgets(line, sizeof(line), fp) != NULL) {
    int transform = 0;
    int source = 0;
    int destination = 0;
    int sign = 1;
    int offset = -1;
    int fields;
    size_t sourceKey;
    size_t destinationKey;

    if (row >= expectedRows) {
      fprintf(stderr,
              "Error: extra projection mapping row in %s "
              "(expected %zu rows).\n",
              defname, expectedRows);
      info = 1;
      break;
    }
    fields = sscanf(line, "%d %d %d %d %n", &transform, &source,
                    &destination, &sign, &offset);
    if (fields == 3 && apFlag == 0) {
      offset = -1;
      sign = 1;
      fields = sscanf(line, "%d %d %d %n", &transform, &source,
                      &destination, &offset);
    }
    if ((fields != 4 && !(fields == 3 && apFlag == 0)) ||
        !LineTailIsWhitespace(line, offset)) {
      fprintf(stderr, "Error: malformed projection mapping row in %s: %s",
              defname, line);
      info = 1;
      break;
    }
    if (transform < 0 || transform >= nArray ||
        source < 0 || source >= nsite ||
        destination < 0 || destination >= nsite) {
      fprintf(stderr,
              "Error: projection mapping index out of range in %s "
              "(transform=%d, source=%d, destination=%d).\n",
              defname, transform, source, destination);
      info = 1;
      break;
    }
    if (sign != -1 && sign != 1) {
      fprintf(stderr,
              "Error: projection sign must be +1 or -1 in %s "
              "(transform=%d, source=%d, sign=%d).\n",
              defname, transform, source, sign);
      info = 1;
      break;
    }
    sourceKey = (size_t)transform * (size_t)nsite + (size_t)source;
    destinationKey =
        (size_t)transform * (size_t)nsite + (size_t)destination;
    if (sourceSeen[sourceKey] != 0) {
      fprintf(stderr,
              "Error: duplicate projection source in %s "
              "(transform=%d, source=%d).\n",
              defname, transform, source);
      info = 1;
      break;
    }
    if (destinationSeen[destinationKey] != 0) {
      fprintf(stderr,
              "Error: projection mapping is not a permutation in %s "
              "(transform=%d, duplicate destination=%d).\n",
              defname, transform, destination);
      info = 1;
      break;
    }

    sourceSeen[sourceKey] = 1;
    destinationSeen[destinationKey] = 1;
    array[transform][source] = destination;
    arraySgn[transform][source] = apFlag != 0 ? sign : 1;
    if (arrayInv != NULL) arrayInv[transform][destination] = source;
    row++;
  }

  if (info == 0 && row != expectedRows) {
    fprintf(stderr,
            "Error: projection mapping row count mismatch in %s "
            "(got %zu, expected %zu).\n",
            defname, row, expectedRows);
    info = 1;
  }
  if (info == 0) {
    int transform;
    int site;
    for (transform = 0; transform < nArray && info == 0; transform++) {
      for (site = 0; site < nsite; site++) {
        size_t key = (size_t)transform * (size_t)nsite + (size_t)site;
        if (sourceSeen[key] == 0 || destinationSeen[key] == 0) {
          fprintf(stderr,
                  "Error: incomplete projection permutation in %s "
                  "(transform=%d, site=%d).\n",
                  defname, transform, site);
          info = 1;
          break;
        }
        if (arrayInv != NULL &&
            arrayInv[transform][array[transform][site]] != site) {
          fprintf(stderr,
                  "Error: inconsistent forward/inverse projection mapping "
                  "in %s (transform=%d, site=%d).\n",
                  defname, transform, site);
          info = 1;
          break;
        }
      }
    }
  }

  free(sourceSeen);
  free(destinationSeen);
  return info;
}

static int ReadTransSymParameters(FILE *fp, double complex *arrayPara,
                                  int nArray, const char *defname) {
  char line[D_CharTmpReadDef];
  unsigned char *seen;
  int row;

  if (fp == NULL || arrayPara == NULL || nArray <= 0) return 1;
  seen = (unsigned char *)calloc((size_t)nArray, sizeof(*seen));
  if (seen == NULL) {
    fprintf(stderr,
            "Error: failed to allocate projection parameter parser state "
            "for %s.\n",
            defname);
    return 1;
  }
  for (row = 0; row < nArray; row++) {
    int transform = 0;
    int offset = -1;
    int fields;
    double realPart = 0.0;
    double imagPart = 0.0;
    if (fgets(line, sizeof(line), fp) == NULL) {
      fprintf(stderr,
              "Error: missing projection parameter row in %s "
              "(got %d, expected %d).\n",
              defname, row, nArray);
      free(seen);
      return 1;
    }
    fields = sscanf(line, "%d %lf %lf %n", &transform, &realPart,
                    &imagPart, &offset);
    if (fields == 2) {
      offset = -1;
      imagPart = 0.0;
      fields = sscanf(line, "%d %lf %n", &transform, &realPart, &offset);
    }
    if ((fields != 2 && fields != 3) ||
        !LineTailIsWhitespace(line, offset) || !isfinite(realPart) ||
        !isfinite(imagPart)) {
      fprintf(stderr,
              "Error: malformed projection parameter row in %s: %s",
              defname, line);
      free(seen);
      return 1;
    }
    if (transform < 0 || transform >= nArray) {
      fprintf(stderr,
              "Error: projection parameter index out of range in %s "
              "(index=%d).\n",
              defname, transform);
      free(seen);
      return 1;
    }
    if (seen[transform] != 0) {
      fprintf(stderr,
              "Error: duplicate projection parameter index in %s "
              "(index=%d).\n",
              defname, transform);
      free(seen);
      return 1;
    }
    seen[transform] = 1;
    arrayPara[transform] = realPart + I * imagPart;
  }
  for (row = 0; row < nArray; row++) {
    if (seen[row] == 0) {
      fprintf(stderr,
              "Error: missing projection parameter index in %s "
              "(index=%d).\n",
              defname, row);
      free(seen);
      return 1;
    }
  }
  free(seen);
  return 0;
}

int GetInfoTransSym(FILE *fp, int **Array, int **ArraySgn, int **ArrayInv,
                    double complex *ArrayPara, int _APFlag, int Nsite,
                    int NArray, char *defname) {
  if (NArray <= 0 || Nsite <= 0) {
    fprintf(stderr, "Error: invalid TransSym dimensions in %s.\n", defname);
    return 1;
  }
  if (ReadTransSymParameters(fp, ArrayPara, NArray, defname) != 0) return 1;
  return ReadProjectionMappings(fp, Array, ArraySgn, ArrayInv, _APFlag,
                                Nsite, NArray, defname);
}

int
GetInfoOneBodyG(FILE *fp, int **ArrayIdx, int **ArrayToIdx, int IndirectGFOn, int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0, x2 = 0, x3 = 0;
  int isite1 = 0, isite2 = 0;
  int idxLanczos = 0;
  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d\n", &x0, &x1, &x2, &x3);
    if (!IndirectGFOn) { // Normal
      ArrayIdx[idx][0] = x0;
      ArrayIdx[idx][1] = x1;
      ArrayIdx[idx][2] = x2;
      ArrayIdx[idx][3] = x3;
    } else { //For Calc Green func indirectly
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

  if (!IndirectGFOn && idx != NArray)
    info = ReadDefFileError(defname);
  return info;
}

// Formerly CisAjsCktAlt
int GetInfoTwoBodyGEx(FILE *fp, int **ArrayIdx, int **ArrayToIdx, int **ArrayIdxOneBodyG,
                      int Nsite, int NArray, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0;
  int x0 = 0, x1 = 0, x2 = 0, x3 = 0;
  int x4 = 0, x5 = 0, x6 = 0, x7 = 0;
  int isite1 = 0, isite2 = 0;
  int idxLanczos = 0;

  if (NArray == 0) return 0;
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    sscanf(ctmp2, "%d %d %d %d %d %d %d %d\n", &x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7);

    isite1 = x0 + x1 * Nsite;
    isite2 = x2 + x3 * Nsite;
    idxLanczos = ArrayToIdx[isite1][isite2];
    ArrayIdxOneBodyG[idxLanczos][0] = x0;
    ArrayIdxOneBodyG[idxLanczos][1] = x1;
    ArrayIdxOneBodyG[idxLanczos][2] = x2;
    ArrayIdxOneBodyG[idxLanczos][3] = x3;
    ArrayIdx[idx][0] = idxLanczos;

    isite1 = x6 + x7 * Nsite;
    isite2 = x4 + x5 * Nsite;
    idxLanczos = ArrayToIdx[isite1][isite2];
    /*
    ArrayIdxOneBodyG[idxLanczos][0] = x4;
    ArrayIdxOneBodyG[idxLanczos][1] = x5;
    ArrayIdxOneBodyG[idxLanczos][2] = x6;
    ArrayIdxOneBodyG[idxLanczos][3] = x7;
     */
    ArrayIdxOneBodyG[idxLanczos][0] = x6;
    ArrayIdxOneBodyG[idxLanczos][1] = x7;
    ArrayIdxOneBodyG[idxLanczos][2] = x4;
    ArrayIdxOneBodyG[idxLanczos][3] = x5;
    ArrayIdx[idx][1] = idxLanczos;

    if (CheckQuadSite(x0, x2, x4, x6, Nsite) != 0) {
      fprintf(stderr, "Error: Site index is incorrect. \n");
      info = 1;
      break;
    }
    idx++;
  }
  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

// Formerly CisAjsCktAltDC
int GetInfoTwoBodyG(FILE *fp, int **ArrayIdx, int Nsite, int NArray, char *defname) {
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

    idx++;
  }
  if (idx != NArray) info = ReadDefFileError(defname);
  return info;
}

static int ReadNBodyGLine(FILE *fp, char **line, size_t *cap) {
  int c;
  size_t len = 0;

  if (*line == NULL || *cap == 0) {
    *cap = 256;
    *line = (char *)malloc(*cap);
    if (*line == NULL) return -1;
  }

  while ((c = fgetc(fp)) != EOF) {
    if (len + 1 >= *cap) {
      char *tmp;
      size_t newCap;
      if (*cap > ((size_t)-1) / 2) return -1;
      newCap = (*cap) * 2;
      tmp = (char *)realloc(*line, newCap);
      if (tmp == NULL) return -1;
      *line = tmp;
      *cap = newCap;
    }
    (*line)[len++] = (char)c;
    if (c == '\n') break;
  }

  if (len == 0 && c == EOF) return 0;

  (*line)[len] = '\0';
  return 1;
}

static int IsIgnorableNBodyGLine(const char *line) {
  const unsigned char *p = (const unsigned char *)line;
  while (isspace(*p)) p++;
  return (*p == '\0' || *p == '#');
}

static int ParseNBodyGHeaderCount(const char *line, int *nbody, char *defname) {
  const unsigned char *p = (const unsigned char *)line;
  char *endptr;
  long nLong;

  while (isspace(*p)) p++;
  if (*p == '\0') {
    fprintf(stderr, "Error: missing NBodyG header count in %s.\n", defname);
    return 1;
  }
  while (*p != '\0' && !isspace(*p)) p++;
  while (isspace(*p)) p++;

  errno = 0;
  nLong = strtol((const char *)p, &endptr, 10);
  if ((const char *)p == endptr || errno != 0 ||
      nLong < 0 || nLong > INT_MAX) {
    fprintf(stderr, "Error: NNBodyG must be non-negative in %s: %s", defname, line);
    return 1;
  }
  while (isspace((unsigned char)*endptr)) endptr++;
  if (*endptr != '\0') {
    fprintf(stderr, "Error: extra token in NBodyG header in %s: %s", defname, line);
    return 1;
  }

  *nbody = (int)nLong;
  return 0;
}

static int ReadNBodyGHeader(FILE *fp, int *nbody, char *defname) {
  char *line = NULL;
  size_t cap = 0;
  int i;

  *nbody = 0;
  for (i = 0; i < IgnoreLinesInDef; i++) {
    int status = ReadNBodyGLine(fp, &line, &cap);
    if (status < 0) {
      fprintf(stderr, "Error: failed to read NBodyG header in %s.\n", defname);
      free(line);
      return 1;
    }
    if (status == 0) {
      fprintf(stderr, "Error: incomplete NBodyG header in %s.\n", defname);
      free(line);
      return 1;
    }
    if (i == 1 && ParseNBodyGHeaderCount(line, nbody, defname) != 0) {
      free(line);
      return 1;
    }
  }

  free(line);
  return 0;
}

static int ValidateNBodyGFactors(int n, int *values, int nsite,
                                 int checkSpinChange, int allowSpinChange,
                                 char *defname) {
  int k;
  for (k = 0; k < n; k++) {
    const int siteOut = values[4*k + 0];
    const int spinOut = values[4*k + 1];
    const int siteIn  = values[4*k + 2];
    const int spinIn  = values[4*k + 3];

    if (CheckPairSite(siteOut, siteIn, nsite) != 0) {
      fprintf(stderr,
              "Error: site index of NBodyG is incorrect in %s "
              "(site_out=%d, site_in=%d, Nsite=%d).\n",
              defname, siteOut, siteIn, nsite);
      return 1;
    }
    if (spinOut < 0 || spinOut > 1 || spinIn < 0 || spinIn > 1) {
      fprintf(stderr,
              "Error: spin index of NBodyG is incorrect in %s "
              "(spin_out=%d, spin_in=%d).\n",
              defname, spinOut, spinIn);
      return 1;
    }
    if (checkSpinChange && !allowSpinChange && spinOut != spinIn) {
      fprintf(stderr,
              "Error: spin-changing NBodyG factor requires orbital-general "
              "mode in %s (spin_out=%d, spin_in=%d).\n",
              defname, spinOut, spinIn);
      return 1;
    }
  }
  return 0;
}

static int ParseNBodyGTermLine(const char *line, int nsite,
                               int checkSpinChange, int allowSpinChange,
                               int *nOut, int **valuesOut,
                               char *defname) {
  char *endptr;
  const char *p = line;
  long nLong;
  int n;
  int k;
  int *values = NULL;

  while (isspace((unsigned char)*p)) p++;
  errno = 0;
  nLong = strtol(p, &endptr, 10);
  if (p == endptr || errno != 0 || nLong < 1 ||
      nLong > INT_MAX / 4) {
    fprintf(stderr, "Error: invalid NBodyG term order in %s: %s", defname, line);
    return 1;
  }
  n = (int)nLong;
  p = endptr;

  values = (int *)malloc(sizeof(int) * 4 * (size_t)n);
  if (values == NULL) {
    fprintf(stderr, "Error: memory allocation failed in NBodyG parser.\n");
    return 1;
  }

  for (k = 0; k < 4*n; k++) {
    long v;
    while (isspace((unsigned char)*p)) p++;
    errno = 0;
    v = strtol(p, &endptr, 10);
    if (p == endptr || errno != 0 || v < INT_MIN || v > INT_MAX) {
      fprintf(stderr, "Error: malformed NBodyG term in %s: %s", defname, line);
      free(values);
      return 1;
    }
    values[k] = (int)v;
    p = endptr;
  }

  while (isspace((unsigned char)*p)) p++;
  if (*p != '\0') {
    fprintf(stderr, "Error: extra token in NBodyG term in %s: %s", defname, line);
    free(values);
    return 1;
  }

  if (ValidateNBodyGFactors(n, values, nsite,
                            checkSpinChange, allowSpinChange,
                            defname) != 0) {
    free(values);
    return 1;
  }

  *nOut = n;
  *valuesOut = values;
  return 0;
}

int ReadBuffNBodyG(FILE *fp, int *nbody, int *totalFactors,
                   int *maxN, int nsite, char *defname) {
  char *line = NULL;
  size_t cap = 0;
  int headerN = 0;
  int count = 0;
  int total = 0;
  int maxSeen = 0;
  int status;
  int info = 0;

  *nbody = 0;
  *totalFactors = 0;
  *maxN = 0;

  if (ReadNBodyGHeader(fp, &headerN, defname) != 0) return 1;

  while ((status = ReadNBodyGLine(fp, &line, &cap)) > 0) {
    int n = 0;
    int *values = NULL;

    if (IsIgnorableNBodyGLine(line)) continue;
    if (count >= headerN) {
      fprintf(stderr, "Error: too many NBodyG terms in %s.\n", defname);
      info = 1;
      break;
    }
    if (ParseNBodyGTermLine(line, nsite, 0, 1, &n, &values, defname) != 0) {
      info = 1;
      break;
    }
    if (total > INT_MAX - n) {
      fprintf(stderr, "Error: NBodyG factor count is too large in %s.\n", defname);
      free(values);
      info = 1;
      break;
    }
    total += n;
    if (n > maxSeen) maxSeen = n;
    count++;
    free(values);
  }

  if (status < 0) {
    fprintf(stderr, "Error: failed to read NBodyG terms in %s.\n", defname);
    info = 1;
  }
  free(line);

  if (info == 0 && count != headerN) {
    fprintf(stderr,
            "Error: NBodyG term count mismatch in %s "
            "(header=%d, data=%d).\n",
            defname, headerN, count);
    info = 1;
  }
  if (info != 0) return 1;

  *nbody = headerN;
  *totalFactors = total;
  *maxN = maxSeen;
  return 0;
}

int GetInfoNBodyG(FILE *fp, int *termN, int *termOffset, int **termIdx,
                  int nsite, int nbody, int totalFactors,
                  int maxN, int allowSpinChange, char *defname) {
  char *line = NULL;
  size_t cap = 0;
  int count = 0;
  int offset = 0;
  int maxSeen = 0;
  int status;
  int info = 0;

  while ((status = ReadNBodyGLine(fp, &line, &cap)) > 0) {
    int n = 0;
    int *values = NULL;
    int k;

    if (IsIgnorableNBodyGLine(line)) continue;
    if (count >= nbody) {
      fprintf(stderr, "Error: too many NBodyG terms in %s.\n", defname);
      info = 1;
      break;
    }
    if (ParseNBodyGTermLine(line, nsite, 1, allowSpinChange,
                            &n, &values, defname) != 0) {
      info = 1;
      break;
    }
    if (offset > totalFactors - n) {
      fprintf(stderr, "Error: NBodyG factor count changed in %s.\n", defname);
      free(values);
      info = 1;
      break;
    }

    termN[count] = n;
    termOffset[count] = offset;
    for (k = 0; k < n; k++) {
      termIdx[offset + k][0] = values[4*k + 0];
      termIdx[offset + k][1] = values[4*k + 1];
      termIdx[offset + k][2] = values[4*k + 2];
      termIdx[offset + k][3] = values[4*k + 3];
    }
    offset += n;
    if (n > maxSeen) maxSeen = n;
    count++;
    free(values);
  }

  if (status < 0) {
    fprintf(stderr, "Error: failed to read NBodyG terms in %s.\n", defname);
    info = 1;
  }
  free(line);

  if (info == 0 &&
      (count != nbody || offset != totalFactors || maxSeen != maxN)) {
    fprintf(stderr,
            "Error: NBodyG size mismatch in %s "
            "(terms=%d/%d, factors=%d/%d, maxN=%d/%d).\n",
            defname, count, nbody, offset, totalFactors, maxSeen, maxN);
    info = 1;
  }

  return info;
}

static int ParseNBodyInterAllHeaderCount(const char *line, int *nbody,
                                         char *defname) {
  const unsigned char *p = (const unsigned char *)line;
  char *endptr;
  long nLong;

  while (isspace(*p)) p++;
  if (*p == '\0') {
    fprintf(stderr, "Error: missing NBodyInterAll header count in %s.\n", defname);
    return 1;
  }
  while (*p != '\0' && !isspace(*p)) p++;
  while (isspace(*p)) p++;

  errno = 0;
  nLong = strtol((const char *)p, &endptr, 10);
  if ((const char *)p == endptr || errno != 0 ||
      nLong < 0 || nLong > INT_MAX) {
    fprintf(stderr, "Error: NNBodyInterAll must be non-negative in %s: %s",
            defname, line);
    return 1;
  }
  while (isspace((unsigned char)*endptr)) endptr++;
  if (*endptr != '\0') {
    fprintf(stderr, "Error: extra token in NBodyInterAll header in %s: %s",
            defname, line);
    return 1;
  }

  *nbody = (int)nLong;
  return 0;
}

static int ReadNBodyInterAllHeader(FILE *fp, int *nbody, char *defname) {
  char *line = NULL;
  size_t cap = 0;
  int i;

  *nbody = 0;
  for (i = 0; i < IgnoreLinesInDef; i++) {
    int status = ReadNBodyGLine(fp, &line, &cap);
    if (status < 0) {
      fprintf(stderr, "Error: failed to read NBodyInterAll header in %s.\n",
              defname);
      free(line);
      return 1;
    }
    if (status == 0) {
      fprintf(stderr, "Error: incomplete NBodyInterAll header in %s.\n",
              defname);
      free(line);
      return 1;
    }
    if (i == 1 &&
        ParseNBodyInterAllHeaderCount(line, nbody, defname) != 0) {
      free(line);
      return 1;
    }
  }

  free(line);
  return 0;
}

static int ValidateNBodyInterAllFactors(int n, int *values, int nsite,
                                        int checkSpinChange,
                                        int allowSpinChange,
                                        char *defname) {
  int k;
  for (k = 0; k < n; k++) {
    const int siteOut = values[4*k + 0];
    const int spinOut = values[4*k + 1];
    const int siteIn  = values[4*k + 2];
    const int spinIn  = values[4*k + 3];

    if (CheckPairSite(siteOut, siteIn, nsite) != 0) {
      fprintf(stderr,
              "Error: site index of NBodyInterAll is incorrect in %s "
              "(site_out=%d, site_in=%d, Nsite=%d).\n",
              defname, siteOut, siteIn, nsite);
      return 1;
    }
    if (spinOut < 0 || spinOut > 1 || spinIn < 0 || spinIn > 1) {
      fprintf(stderr,
              "Error: spin index of NBodyInterAll is incorrect in %s "
              "(spin_out=%d, spin_in=%d).\n",
              defname, spinOut, spinIn);
      return 1;
    }
    if (checkSpinChange && !allowSpinChange && spinOut != spinIn) {
      fprintf(stderr,
              "Error: spin-changing NBodyInterAll factor requires "
              "orbital-general mode in %s (spin_out=%d, spin_in=%d).\n",
              defname, spinOut, spinIn);
      return 1;
    }
  }
  return 0;
}

static int ParseNBodyInterAllTermLine(const char *line, int nsite,
                                      int checkSpinChange,
                                      int allowSpinChange, int *nOut,
                                      int **valuesOut,
                                      double complex *coefOut,
                                      char *defname) {
  char *endptr;
  const char *p = line;
  long nLong;
  int n;
  int k;
  int *values = NULL;
  double dReValue;
  double dImValue;

  while (isspace((unsigned char)*p)) p++;
  errno = 0;
  nLong = strtol(p, &endptr, 10);
  if (p == endptr || errno != 0 || nLong < 1 ||
      nLong > INT_MAX / 4) {
    fprintf(stderr, "Error: invalid NBodyInterAll term order in %s: %s",
            defname, line);
    return 1;
  }
  n = (int)nLong;
  p = endptr;

  values = (int *)malloc(sizeof(int) * 4 * (size_t)n);
  if (values == NULL) {
    fprintf(stderr, "Error: memory allocation failed in NBodyInterAll parser.\n");
    return 1;
  }

  for (k = 0; k < 4*n; k++) {
    long v;
    while (isspace((unsigned char)*p)) p++;
    errno = 0;
    v = strtol(p, &endptr, 10);
    if (p == endptr || errno != 0 || v < INT_MIN || v > INT_MAX) {
      fprintf(stderr, "Error: malformed NBodyInterAll term in %s: %s",
              defname, line);
      free(values);
      return 1;
    }
    values[k] = (int)v;
    p = endptr;
  }

  while (isspace((unsigned char)*p)) p++;
  errno = 0;
  dReValue = strtod(p, &endptr);
  if (p == endptr || errno != 0 || !isfinite(dReValue)) {
    fprintf(stderr, "Error: invalid NBodyInterAll coefficient in %s: %s",
            defname, line);
    free(values);
    return 1;
  }
  p = endptr;

  while (isspace((unsigned char)*p)) p++;
  errno = 0;
  dImValue = strtod(p, &endptr);
  if (p == endptr || errno != 0 || !isfinite(dImValue)) {
    fprintf(stderr, "Error: invalid NBodyInterAll coefficient in %s: %s",
            defname, line);
    free(values);
    return 1;
  }
  p = endptr;

  while (isspace((unsigned char)*p)) p++;
  if (*p != '\0') {
    fprintf(stderr, "Error: extra token in NBodyInterAll term in %s: %s",
            defname, line);
    free(values);
    return 1;
  }

  if (ValidateNBodyInterAllFactors(n, values, nsite,
                                   checkSpinChange, allowSpinChange,
                                   defname) != 0) {
    free(values);
    return 1;
  }

  *nOut = n;
  *valuesOut = values;
  *coefOut = dReValue + I * dImValue;
  return 0;
}

int ReadBuffNBodyInterAll(FILE *fp, int *nbody, int *totalFactors,
                          int *maxN, int nsite, char *defname) {
  char *line = NULL;
  size_t cap = 0;
  int headerN = 0;
  int count = 0;
  int total = 0;
  int maxSeen = 0;
  int status;
  int info = 0;

  *nbody = 0;
  *totalFactors = 0;
  *maxN = 0;

  if (ReadNBodyInterAllHeader(fp, &headerN, defname) != 0) return 1;
  if (headerN > INT_MAX / 2) {
    fprintf(stderr, "Error: NNBodyInterAll is too large in %s.\n", defname);
    return 1;
  }

  while ((status = ReadNBodyGLine(fp, &line, &cap)) > 0) {
    int n = 0;
    int *values = NULL;
    double complex coef;

    if (IsIgnorableNBodyGLine(line)) continue;
    if (count >= headerN) {
      fprintf(stderr, "Error: too many NBodyInterAll terms in %s.\n", defname);
      info = 1;
      break;
    }
    if (ParseNBodyInterAllTermLine(line, nsite, 0, 1,
                                   &n, &values, &coef, defname) != 0) {
      info = 1;
      break;
    }
    if (total > INT_MAX - n || total + n > INT_MAX / 4) {
      fprintf(stderr,
              "Error: NBodyInterAll factor count is too large in %s.\n",
              defname);
      free(values);
      info = 1;
      break;
    }
    total += n;
    if (n > maxSeen) maxSeen = n;
    count++;
    free(values);
  }

  if (status < 0) {
    fprintf(stderr, "Error: failed to read NBodyInterAll terms in %s.\n",
            defname);
    info = 1;
  }
  free(line);

  if (info == 0 && count != headerN) {
    fprintf(stderr,
            "Error: NBodyInterAll term count mismatch in %s "
            "(header=%d, data=%d).\n",
            defname, headerN, count);
    info = 1;
  }
  if (info != 0) return 1;

  *nbody = headerN;
  *totalFactors = total;
  *maxN = maxSeen;
  return 0;
}

int GetInfoNBodyInterAll(FILE *fp, int *termN, int *termOffset,
                         int **termIdx, double complex *termCoef,
                         int nsite, int nbody, int totalFactors,
                         int maxN, int allowSpinChange, char *defname) {
  char *line = NULL;
  size_t cap = 0;
  int count = 0;
  int offset = 0;
  int maxSeen = 0;
  int status;
  int info = 0;

  while ((status = ReadNBodyGLine(fp, &line, &cap)) > 0) {
    int n = 0;
    int *values = NULL;
    double complex coef;
    int k;

    if (IsIgnorableNBodyGLine(line)) continue;
    if (count >= nbody) {
      fprintf(stderr, "Error: too many NBodyInterAll terms in %s.\n", defname);
      info = 1;
      break;
    }
    if (ParseNBodyInterAllTermLine(line, nsite, 1, allowSpinChange,
                                   &n, &values, &coef, defname) != 0) {
      info = 1;
      break;
    }
    if (offset > totalFactors - n) {
      fprintf(stderr,
              "Error: NBodyInterAll factor count changed in %s.\n",
              defname);
      free(values);
      info = 1;
      break;
    }

    termN[count] = n;
    termOffset[count] = offset;
    termCoef[count] = coef;
    for (k = 0; k < n; k++) {
      termIdx[offset + k][0] = values[4*k + 0];
      termIdx[offset + k][1] = values[4*k + 1];
      termIdx[offset + k][2] = values[4*k + 2];
      termIdx[offset + k][3] = values[4*k + 3];
    }
    offset += n;
    if (n > maxSeen) maxSeen = n;
    count++;
    free(values);
  }

  if (status < 0) {
    fprintf(stderr, "Error: failed to read NBodyInterAll terms in %s.\n",
            defname);
    info = 1;
  }
  free(line);

  if (info == 0 &&
      (count != nbody || offset != totalFactors || maxSeen != maxN)) {
    fprintf(stderr,
            "Error: NBodyInterAll size mismatch in %s "
            "(terms=%d/%d, factors=%d/%d, maxN=%d/%d).\n",
            defname, count, nbody, offset, totalFactors, maxSeen, maxN);
    info = 1;
  }

  return info;
}

static int ReadOptTransParameters(FILE *fp, double *arrayPara, int nArray,
                                  const char *defname) {
  char line[D_CharTmpReadDef];
  unsigned char *seen;
  int row;

  if (fp == NULL || arrayPara == NULL || nArray <= 0) return 1;
  seen = (unsigned char *)calloc((size_t)nArray, sizeof(*seen));
  if (seen == NULL) {
    fprintf(stderr,
            "Error: failed to allocate OptTrans parameter parser state "
            "for %s.\n",
            defname);
    return 1;
  }
  for (row = 0; row < nArray; row++) {
    int transform = 0;
    int offset = -1;
    double value = 0.0;
    if (fgets(line, sizeof(line), fp) == NULL) {
      fprintf(stderr,
              "Error: missing OptTrans parameter row in %s "
              "(got %d, expected %d).\n",
              defname, row, nArray);
      free(seen);
      return 1;
    }
    if (sscanf(line, "%d %lf %n", &transform, &value, &offset) != 2 ||
        !LineTailIsWhitespace(line, offset) || !isfinite(value)) {
      fprintf(stderr, "Error: malformed OptTrans parameter row in %s: %s",
              defname, line);
      free(seen);
      return 1;
    }
    if (transform < 0 || transform >= nArray) {
      fprintf(stderr,
              "Error: OptTrans parameter index out of range in %s "
              "(index=%d).\n",
              defname, transform);
      free(seen);
      return 1;
    }
    if (seen[transform] != 0) {
      fprintf(stderr,
              "Error: duplicate OptTrans parameter index in %s "
              "(index=%d).\n",
              defname, transform);
      free(seen);
      return 1;
    }
    seen[transform] = 1;
    arrayPara[transform] = value;
  }
  for (row = 0; row < nArray; row++) {
    if (seen[row] == 0) {
      fprintf(stderr,
              "Error: missing OptTrans parameter index in %s "
              "(index=%d).\n",
              defname, row);
      free(seen);
      return 1;
    }
  }
  free(seen);
  return 0;
}

int GetInfoOptTrans(FILE *fp, int **Array, double *ArrayPara,
                    int *ArrayOpt, int **ArraySgn, int _iFlagOptTrans,
                    int *iOptCount, int _fidx, int _APFlag, int Nsite,
                    int NArray, char *defname) {
  int i;
  if (_iFlagOptTrans <= 0) return 0;
  if (NArray <= 0 || Nsite <= 0 || ArrayOpt == NULL || iOptCount == NULL) {
    fprintf(stderr, "Error: invalid OptTrans dimensions in %s.\n", defname);
    return 1;
  }
  if (ReadOptTransParameters(fp, ArrayPara, NArray, defname) != 0) return 1;
  if (ReadProjectionMappings(fp, Array, ArraySgn, NULL, _APFlag, Nsite,
                             NArray, defname) != 0) {
    return 1;
  }
  for (i = 0; i < NArray; i++) ArrayOpt[_fidx + i] = 1;
  *iOptCount += NArray;
  return 0;
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

int GetInfoAnomalousTerm(FILE *fp, int **ArrayIdx,
                         double complex *ArrayValue, int Nsite, int NArray,
                         char *defname) {
  char line[D_FileNameMax];
  int idx = 0;
  int type, s1, sp1, s2, sp2;
  double re, im;
  int nread;

  if (NArray == 0) return 0;
  while (fgets(line, sizeof(line), fp) != NULL) {
    if (idx >= NArray) {
      fprintf(stderr, "Error in %s: too many AnomalousTerm rows.\n",
              defname);
      return 1;
    }
    nread = 0;
    if (sscanf(line, "%d %d %d %d %d %lf %lf%n",
               &type, &s1, &sp1, &s2, &sp2, &re, &im, &nread) != 7) {
      fprintf(stderr, "Error in %s: malformed AnomalousTerm row.\n",
              defname);
      return 1;
    }
    if (!LineTailIsWhitespace(line, nread)) {
      fprintf(stderr,
              "Error in %s: AnomalousTerm line has extra fields.\n",
              defname);
      return 1;
    }
    if (!isfinite(re) || !isfinite(im)) {
      fprintf(stderr,
              "Error in %s: AnomalousTerm coefficient is not finite.\n",
              defname);
      return 1;
    }
    if (type != 0 && type != 1) {
      fprintf(stderr, "Error in %s: AnomalousTerm type must be 0 or 1.\n",
              defname);
      return 1;
    }
    if (CheckSite(s1, Nsite) != 0 || CheckSite(s2, Nsite) != 0) {
      fprintf(stderr,
              "Error in %s: Site index of AnomalousTerm is incorrect.\n",
              defname);
      return 1;
    }
    if (sp1 < 0 || sp1 > 1 || sp2 < 0 || sp2 > 1) {
      fprintf(stderr,
              "Error in %s: Spin index of AnomalousTerm is incorrect.\n",
              defname);
      return 1;
    }
    if (s1 == s2 && sp1 == sp2) {
      fprintf(stderr,
              "Error in %s: AnomalousTerm cannot use the same fermion "
              "operator twice in one pair.\n",
              defname);
      return 1;
    }
    ArrayIdx[idx][0] = type;
    ArrayIdx[idx][1] = s1;
    ArrayIdx[idx][2] = sp1;
    ArrayIdx[idx][3] = s2;
    ArrayIdx[idx][4] = sp2;
    ArrayValue[idx] = re + im * I;
    idx++;
  }
  if (idx != NArray) {
    fprintf(stderr,
            "Error in %s: AnomalousTerm row count does not match header.\n",
            defname);
    return 1;
  }
  for (idx = 0; idx + 1 < NArray; idx += 2) {
    const int *a = ArrayIdx[idx];
    const int *b = ArrayIdx[idx + 1];
    if (b[0] != 1 - a[0] || b[1] != a[3] || b[2] != a[4] ||
        b[3] != a[1] || b[4] != a[2] ||
        cabs(ArrayValue[idx + 1] - conj(ArrayValue[idx])) >
            GC_ANOMALOUS_HERMITE_EPS) {
      fprintf(stderr,
              "Error in %s: AnomalousTerm Hermite pair is inconsistent "
              "at rows %d-%d.\n",
              defname, idx + 1, idx + 2);
      return 1;
    }
  }
  return 0;
}

int GetInfoAnomalousG(FILE *fp, int **ArrayIdx, int Nsite, int NArray,
                      char *defname) {
  char line[D_FileNameMax];
  int idx = 0;
  int type, s1, sp1, s2, sp2;
  int nread;

  if (NArray == 0) return 0;
  while (fgets(line, sizeof(line), fp) != NULL) {
    if (idx >= NArray) {
      fprintf(stderr, "Error in %s: too many AnomalousG rows.\n", defname);
      return 1;
    }
    nread = 0;
    if (sscanf(line, "%d %d %d %d %d%n",
               &type, &s1, &sp1, &s2, &sp2, &nread) != 5) {
      fprintf(stderr, "Error in %s: malformed AnomalousG row.\n", defname);
      return 1;
    }
    if (!LineTailIsWhitespace(line, nread)) {
      fprintf(stderr, "Error in %s: AnomalousG line has extra fields.\n",
              defname);
      return 1;
    }
    if (type != 0 && type != 1) {
      fprintf(stderr, "Error in %s: AnomalousG type must be 0 or 1.\n",
              defname);
      return 1;
    }
    if (CheckSite(s1, Nsite) != 0 || CheckSite(s2, Nsite) != 0) {
      fprintf(stderr, "Error in %s: Site index of AnomalousG is incorrect.\n",
              defname);
      return 1;
    }
    if (sp1 < 0 || sp1 > 1 || sp2 < 0 || sp2 > 1) {
      fprintf(stderr, "Error in %s: Spin index of AnomalousG is incorrect.\n",
              defname);
      return 1;
    }
    if (s1 == s2 && sp1 == sp2) {
      fprintf(stderr,
              "Error in %s: AnomalousG cannot use the same fermion "
              "operator twice in one pair.\n",
              defname);
      return 1;
    }
    ArrayIdx[idx][0] = type;
    ArrayIdx[idx][1] = s1;
    ArrayIdx[idx][2] = sp1;
    ArrayIdx[idx][3] = s2;
    ArrayIdx[idx][4] = sp2;
    idx++;
  }
  if (idx != NArray) {
    fprintf(stderr,
            "Error in %s: AnomalousG row count does not match header.\n",
            defname);
    return 1;
  }
  return 0;
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

int GetInfoRBM_Layer(FILE *fp, int *ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nlayer, int NArray,
                      char *defname) {
  int idx0 = 0, idx1 = 0, info = 0;
  int i = 0;
  int fidx = _fidx;

  if (NArray > 0) {
    idx0 = idx1 = 0;
    while (fscanf(fp, "%d ", &i) != EOF) {
      fscanf(fp, "%d\n", &(ArrayIdx[i]));
      if (CheckSite(i, Nlayer) != 0) {
        fprintf(stderr, "Error: Layer index is incorrect. (i=%d, Nlayer=%d) \n",i,Nlayer);
        info = 1;
        break;
      }
      idx0++;
      if (idx0 == Nlayer) break;
    }

    idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    if (idx0 != Nlayer || idx1 != NArray) {
      info = ReadDefFileError(defname);
    }
  }
  return info;
}

int GetInfoGeneralRBM_Layer(FILE *fp, int *ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nlayer, int NArray,
                      char *defname) {
  int idx0 = 0, idx1 = 0, info = 0;
  int i = 0, s = 0;
  int fidx = _fidx;

  if (NArray > 0) {
    idx0 = idx1 = 0;
    while (fscanf(fp, "%d %d ", &i, &s) != EOF) {
      fscanf(fp, "%d\n", &(ArrayIdx[i + s*Nsite]));
      if (CheckSite(i, Nlayer) != 0 || CheckSite(s, 2) != 0) {
        fprintf(stderr, "Error: Layer index is incorrect. (i=%d, Nlayer=%d) \n",i,Nlayer);
        info = 1;
        break;
      }
      idx0++;
      if (idx0 == 2*Nlayer) break;
    }

    idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    if (idx0 != 2*Nlayer || idx1 != NArray) {
      info = ReadDefFileError(defname);
    }
  }
  return info;
}


int GetInfoRBM_PhysHidden(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nlayer, int Nlayer2, int NArray,
                      char *defname) {
  int idx0 = 0, idx1 = 0, info = 0;
  int i = 0, j = 0;
  int fidx = _fidx;
  if (NArray > 0) {
    while (fscanf(fp, "%d %d ", &i, &j) != EOF) {
      if (CheckSite(i, Nlayer) != 0 || CheckSite(j, Nlayer2)) {
        fprintf(stderr, "Error: Site index is incorrect. \n");
        info = 1;
        break;
      }

      fscanf(fp, "%d\n", &(ArrayIdx[i][j]));
      idx0++;
      if (idx0 == Nlayer * Nlayer2) break;
    }
    idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    if (idx0 != Nlayer * Nlayer2 || idx1 != NArray) {
      info = ReadDefFileError(defname);
    }
  }
  return info;

}

int GetInfoGeneralRBM_PhysHidden(FILE *fp, int **ArrayIdx, int *ArrayOpt, int iComplxFlag, int *iOptCount, int _fidx, int Nlayer, int Nlayer2, int NArray,
                      char *defname) {
  int idx0 = 0, idx1 = 0, info = 0;
  int i = 0, j = 0, s = 0;
  int fidx = _fidx;
  if (NArray > 0) {
    while (fscanf(fp, "%d %d %d ", &i, &s, &j) != EOF) {
      if (CheckSite(i, Nlayer) != 0 || CheckSite(s, 2) != 0 || CheckSite(j, Nlayer2) != 0) {
        fprintf(stderr, "Error: Site index is incorrect. \n");
        info = 1;
        break;
      }

      fscanf(fp, "%d\n", &(ArrayIdx[i+s*Nsite][j]) );
      idx0++;
      if (idx0 == 2*Nlayer * Nlayer2) break;
    }
    idx1 = GetInfoOpt(fp, ArrayOpt, iComplxFlag, iOptCount, fidx);
    if (idx0 != 2*Nlayer * Nlayer2 || idx1 != NArray) {
      info = ReadDefFileError(defname);
    }
  }
  return info;

}

int GetInfoLattice(FILE *fp, int **ArrayIdx, int NArray, int nx, int ny, int nz, int norb, char *defname) {
  char ctmp2[256];
  int idx = 0, info = 0, scanned;
  int x1 = 0, x2 = 0, x3 = 0, x4 = 0;
  int *seen = NULL;
  int i;
  if (NArray == 0) return 0;
  seen = (int *)calloc((size_t)NArray, sizeof(int));
  if (seen == NULL) {
    fprintf(stderr, "Error: GetInfoLattice failed to allocate seen[].\n");
    return 1;
  }
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    scanned = sscanf(ctmp2, "%d %d %d %d %d\n",
                     &idx, &x1, &x2, &x3, &x4);
    if (scanned != 5) {
      fprintf(stderr, "Error: malformed lattice.def line (expected 5 ints): %s",
              ctmp2);
      info = 1;
      break;
    }
    if (idx < 0 || idx >= NArray) {
      fprintf(stderr, "Error: lattice idx=%d out of range [0,%d).\n", idx, NArray);
      info = 1;
      break;
    }
    if (CheckSite(x1, nx) != 0 || CheckSite(x2, ny) != 0 ||
        CheckSite(x3, nz) != 0 || CheckSite(x4, norb) != 0) {
      fprintf(stderr, "Error: Site index for Lattice is incorrect (idx=%d, %d %d %d %d).\n",
              idx, x1, x2, x3, x4);
      info = 1;
      break;
    }
    if (seen[idx] != 0) {
      fprintf(stderr,
              "Error: lattice idx=%d appears more than once. lattice.def "
              "uses one row per site; the last column is an orbital index, "
              "not a spin index.\n",
              idx);
      info = 1;
      break;
    }
    seen[idx] = 1;
    ArrayIdx[idx][0] = x1;
    ArrayIdx[idx][1] = x2;
    ArrayIdx[idx][2] = x3;
    ArrayIdx[idx][3] = x4;
  }
  if (info == 0) {
    for (i = 0; i < NArray; i++) {
      if (seen[i] == 0) {
        fprintf(stderr, "Error: lattice idx=%d is missing.\n", i);
        info = ReadDefFileError(defname);
        break;
      }
    }
  }
  free(seen);
  return info;
}

int GetInfoTwist(FILE *fp, int **ArrayIdx, double **ArrayValue, int Nsite, int NTwist, char *defname) {
  char ctmp2[256];
  int info = 0, scanned;
  int twist_idx = 0, i = 0, s = 0;
  int row, key, t, k;
  double dReValueX = 0;
  double dReValueY = 0;
  double dReValueZ = 0;
  int *count = NULL;       /* per-twist row count */
  char *seen = NULL;       /* flat 2D: seen[t * (2*Nsite) + key] */
  size_t seen_size;
  if (NTwist == 0 || Nsite == 0) return 0;
  count = (int *)calloc((size_t)NTwist, sizeof(int));
  seen_size = (size_t)NTwist * (size_t)(2 * Nsite);
  seen = (char *)calloc(seen_size, sizeof(char));
  if (count == NULL || seen == NULL) {
    fprintf(stderr, "Error: GetInfoTwist failed to allocate validation buffers.\n");
    free(count);
    free(seen);
    return 1;
  }
  while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
    scanned = sscanf(ctmp2, "%d %d %d %lf %lf %lf \n",
                     &twist_idx, &i, &s, &dReValueX, &dReValueY, &dReValueZ);
    if (scanned != 6) {
      fprintf(stderr, "Error: malformed twist.def line (expected 6 fields): %s",
              ctmp2);
      info = 1;
      break;
    }
    if (twist_idx < 0 || twist_idx >= NTwist) {
      fprintf(stderr, "Error: twist_idx=%d out of range [0,%d).\n",
              twist_idx, NTwist);
      info = 1;
      break;
    }
    if (s != 0 && s != 1) {
      fprintf(stderr, "Error: twist spin index s=%d must be 0 or 1.\n", s);
      info = 1;
      break;
    }
    if (CheckSite(i, Nsite) != 0) {
      fprintf(stderr, "Error: twist site index i=%d is out of [0,%d).\n",
              i, Nsite);
      info = 1;
      break;
    }
    key = i + s * Nsite;     /* unique site-spin slot in [0, 2*Nsite) */
    if (seen[(size_t)twist_idx * (size_t)(2 * Nsite) + (size_t)key] != 0) {
      fprintf(stderr, "Error: twist=%d (site=%d, spin=%d) appears more than once.\n",
              twist_idx, i, s);
      info = 1;
      break;
    }
    seen[(size_t)twist_idx * (size_t)(2 * Nsite) + (size_t)key] = 1;
    row = count[twist_idx];  /* row index within this twist (0..2*Nsite-1) */
    ArrayIdx[twist_idx][2*row]   = i;
    ArrayIdx[twist_idx][2*row+1] = s;
    ArrayValue[twist_idx][3*row]   = dReValueX;
    ArrayValue[twist_idx][3*row+1] = dReValueY;
    ArrayValue[twist_idx][3*row+2] = dReValueZ;
    count[twist_idx]++;
  }
  if (info == 0) {
    for (t = 0; t < NTwist; t++) {
      if (count[t] != 2 * Nsite) {
        fprintf(stderr, "Error: twist=%d has %d rows, expected %d.\n",
                t, count[t], 2 * Nsite);
        info = ReadDefFileError(defname);
        break;
      }
      for (k = 0; k < 2 * Nsite; k++) {
        if (seen[(size_t)t * (size_t)(2 * Nsite) + (size_t)k] != 1) {
          fprintf(stderr, "Error: twist=%d missing entry for slot %d (site=%d, spin=%d).\n",
                  t, k, k % Nsite, k / Nsite);
          info = ReadDefFileError(defname);
          break;
        }
      }
      if (info != 0) break;
    }
  }
  free(count);
  free(seen);
  return info;
}



/**********************************/
/* [e] Read Parameters from file  */
/**********************************/
