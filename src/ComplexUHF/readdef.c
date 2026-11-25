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
 *[ver.2008.11.4]
 *  Read Definition files
 *-------------------------------------------------------------
 * Copyright (C) 2007-2009 Daisuke Tahara. All rights reserved.
 *-------------------------------------------------------------*/


/*=================================================================================================*/
#include "readdef.h"
#include "complex.h"
#include <ctype.h>
/**
 * Keyword List in NameListFile.
 **/
static char cKWListOfFileNameList[][D_CharTmpReadDef]={
  "ModPara", "LocSpin",
  "Trans", "CoulombIntra", "CoulombInter",
  "Hund", "PairHop", "Exchange",
  "Gutzwiller", "Jastrow",
  "DH2", "DH4", "Orbital", "OrbitalAntiParallel",
  "OrbitalParallel", "OrbitalGeneral",
  "TransSym", "InGutzwiller", "InJastrow",
  "InDH2", "InDH4", "InOrbital",
  "OneBodyG", "TwoBodyG", "TwoBodyGEx",
  "InterAll", "OptTrans", "InOptTrans",
  "Initial"
};

/**
 * Number of Keyword List in NameListFile for this program.
 **/
enum KWIdxInt {
    KWModPara, KWLocSpin,
    KWTrans, KWCoulombIntra, KWCoulombInter,
    KWHund, KWPairHop, KWExchange,
    KWGutzwiller, KWJastrow,
    KWDH2, KWDH4, KWOrbital, KWOrbitalAntiParallel,
    KWOrbitalParallel, KWOrbitalGeneral,
    KWTransSym, KWInGutzwiller, KWInJastrow,
    KWInDH2, KWInDH4, KWInOrbital,
    KWOneBodyG, KWTwoBodyG, KWTwoBodyGEx,
    KWInterAll, KWOptTrans, KWInOptTrans,
    KWInitial, KWIdxInt_end
};

static char (*cFileNameListFile)[D_CharTmpReadDef];

int ReadDefFileError(
	const char *defname
){
	printf("Error: %s (Broken file or Not exist)\n", defname);
	return 1;
}

int CheckSite(
        const int iSite,
        const int iMaxNum
);

int CheckPairSite(
        const int iSite1,
        const int iSite2,
        const int iMaxNum
);

int CheckQuadSite(
        const int iSite1,
        const int iSite2,
        const int iSite3,
        const int iSite4,
        const int iMaxNum
);

void SetInitialValue(struct DefineList *X);

/**
 * @brief Check that a definition file does not hold more records than its header declares.
 * @retval 0 the record at idx is inside the declared range.
 * @retval -1 the declared number of records has already been read.
 */
static int CheckRecordCount(const int idx, const int iNum, const char *defname) {
  if (idx >= iNum) {
    fprintf(stderr, "Error in %s: more records than declared (%d).\n", defname, iNum);
    return -1;
  }
  return 0;
}

/**
 * @brief Check that a spin index is 0 or 1.
 */
static int CheckSpin(const int iSpin) {
  if (iSpin != 0 && iSpin != 1) return -1;
  return 0;
}

static int ReportFormat(const char *defname, const int row) {
  fprintf(stderr, "Error in %s: incorrect format at record %d.\n", defname, row);
  return 1;
}

static int ReportSite(const char *defname, const int row) {
  fprintf(stderr, "Error in %s: Site index is incorrect at record %d.\n", defname, row);
  return 1;
}

static int ReportSpin(const char *defname, const int row) {
  fprintf(stderr, "Error in %s: Spin index must be 0 or 1 at record %d.\n", defname, row);
  return 1;
}

static int ReportOrbitalIdx(const char *defname, const int row, const int idx, const int iNum) {
  fprintf(stderr, "Error in %s: orbital index %d is outside [0, %d) at record %d.\n", defname, idx, iNum, row);
  return 1;
}

char* ReadBuffInt(FILE *fp, int *iNbuf){
	char *cerr;
	char ctmp[D_FileNameMax];
	char ctmp2[D_FileNameMax];
	cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
	if(cerr != NULL) {
		cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
		if(cerr != NULL) {
			if (sscanf(ctmp2, "%s %d\n", ctmp, iNbuf) != 2 || *iNbuf < 0) {
				fprintf(stderr, "Error: the number of records in a header must be a non-negative integer.\n");
				cerr = NULL;
			}
		}
	}
	return cerr;
}

int JudgeOrbitalMode(int *_iFlgOrbitalGeneral, const int _iFlgOrbitalAP, const int _iFlgOrbitalP,int *OutputMode){

  int iret=0;
  //(General, AP, P)
  if(*_iFlgOrbitalGeneral==1){
    if(_iFlgOrbitalAP==0 && _iFlgOrbitalP==0){ //(1, 0, 0)
      iret=0; 
      *OutputMode = 0;   // only general
    }
    else{//(1, 1, 0) or (1, 0, 1) or (1, 1, 1)
      iret= -1;
    }
  }
  else{
    if(_iFlgOrbitalAP==1){
      if(_iFlgOrbitalP==1){ //(0, 1, 1)
        *_iFlgOrbitalGeneral=1;
        iret= 0;
        *OutputMode = 2; // AP + P
      }
      else{ //(0, 1, 0)
        iret= 0;
        *OutputMode = 1; // AP
      }
    }
    else{// (0, 0, 0) or (0, 0, 1)
      iret= -2;
    }
  }
  if(iret==-1) {
    fprintf(stderr, "error: Multiple definition of Orbital files.\n");
  }
  else if(iret==-2){
    fprintf(stderr, "error: Not exist any Orbital file or Need OrbitalAP file.\n");
  }
  return iret;
}


int ReadDefFileNInt(
	char *xNameListFile, 
	struct DefineList *X
) {
	FILE *fp;
	char defname[D_FileNameMaxReadDef];
	char ctmp[D_CharTmpReadDef];
	char ctmp2[D_FileNameMax];
	int itmp, info;
	int iKWidx = 0;
	info = 0;
	char *cerr;
	int iNOrbitalAP=0;
    int iNOrbitalP=0;
    int iFlgOrbitalAP=0, iFlgOrbitalP=0;
    X->iFlgOrbital = 0;
#if Fock == 1
    X->iFlg_Fock = 1;
#else
    X->iFlg_Fock = 0;
#endif

	cFileNameListFile = malloc(sizeof(char) * D_CharTmpReadDef * KWIdxInt_end);
	fprintf(stdout, "  Read File %s .\n", xNameListFile);
	if (GetFileName(xNameListFile, cFileNameListFile) != 0) {
		fprintf(stderr, "  error: Definition files(*.def) are incomplete.\n");
		return (-1);
	}

	for (iKWidx = 0; iKWidx < KWIdxInt_end; iKWidx++) {
		strcpy(defname, cFileNameListFile[iKWidx]);
		if (strcmp(defname, "") == 0) {
			switch (iKWidx) {
				case KWModPara:
				case KWLocSpin:
					fprintf(stderr, "  Error: Need to make a def file for %s.\n", cKWListOfFileNameList[iKWidx]);
					return (-1);
				default:
					break;
			}
		}
	}

	SetInitialValue(X);

	for (iKWidx = 0; iKWidx < KWIdxInt_end; iKWidx++) {
        strcpy(defname, cFileNameListFile[iKWidx]);
        if (strcmp(defname, "") == 0) continue;
        fprintf(stdout, "  Read File '%s' for %s.\n", defname, cKWListOfFileNameList[iKWidx]);
        cerr = ctmp;
        fp = fopen(defname, "r");
        if (fp == NULL) {
            ReadDefFileError(defname);
            return (-1);
        } else {
            switch (iKWidx) {
                case KWModPara:
                    /* Read modpara.def---------------------------------------*/
                    //TODO: add error procedure here when parameters are not enough.
                    //SetDefultValuesModPara(X);
                    cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
                    cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
                    sscanf(ctmp2, "%s %d\n", ctmp, &itmp); //2
                    cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //3
                    cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //4
                    cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp); //5
                    cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
                    sscanf(ctmp2, "%s %s\n", ctmp, X->CDataFileHead); //6
                    cerr = fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
                    sscanf(ctmp2, "%s %s\n", ctmp, X->CParaFileHead); //7
                    cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);   //8

                    double dtmp;
                    while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
                        if (*ctmp2 == '\n' || ctmp2[0] == '-') continue;
                        sscanf(ctmp2, "%s %lf\n", ctmp, &dtmp);
                        if (CheckWords(ctmp, "NVMCCalMode") == 0 ||
                            CheckWords(ctmp, "NLanczosMode") == 0 ||
                            CheckWords(ctmp, "NDataIdxStart") == 0 ||
                            CheckWords(ctmp, "NDataQtySmp") == 0 ||
                            CheckWords(ctmp, "NDataQtySmp") == 0 ||
                            CheckWords(ctmp, "NDataQtySmp") == 0 ||
                            CheckWords(ctmp, "NSPGaussLeg") == 0 ||
                            CheckWords(ctmp, "NSPStot") == 0 ||
                            CheckWords(ctmp, "NSROptItrStep") == 0 ||
                            CheckWords(ctmp, "NSROptItrSmp") == 0 ||
                            CheckWords(ctmp, "DSROptRedCut") == 0 ||
                            CheckWords(ctmp, "DSROptStaDel") == 0 ||
                            CheckWords(ctmp, "DSROptStepDt") == 0 ||
                            CheckWords(ctmp, "NVMCWarmUp") == 0 ||
                            CheckWords(ctmp, "NVMCInterval") == 0 ||
                            CheckWords(ctmp, "NVMCSample") == 0 ||
                            CheckWords(ctmp, "NExUpdatePath") == 0 ||
                            CheckWords(ctmp, "NSplitSize") == 0 ||
                            CheckWords(ctmp, "NStore") == 0
                                ) {
                            fprintf(stdout, "!! Warning: %s is not used for Hatree Fock Calculation. !!\n", ctmp);

                            continue;
                        } else if (CheckWords(ctmp, "Nsite") == 0) {
                            X->Nsite = (int) dtmp;
                        } else if (CheckWords(ctmp, "Ne") == 0 || CheckWords(ctmp, "Nelectron") == 0) {
                            X->Ne = (int) dtmp;
                        } else if (CheckWords(ctmp, "Ncond") == 0) {
                            X->Ncond = (int) dtmp;
                        } else if (CheckWords(ctmp, "2Sz") == 0) {
                            X->TwoSz = (int) dtmp;
                            if (X->TwoSz == -1) {
                                fprintf(stdout, "Error: 2Sz must be even number.\n");
                                return (-1);
                            }
                        } else if (CheckWords(ctmp, "Mix") == 0) {
                            X->mix = dtmp;
                        } else if (CheckWords(ctmp, "EPS") == 0) {
                            X->eps_int = (int) dtmp;
                        } else if (CheckWords(ctmp, "Print") == 0) {
                            X->print = (int) dtmp;
                        } else if (CheckWords(ctmp, "IterationMax") == 0) {
                            X->IterationMax = (int) dtmp;
                        } else if (CheckWords(ctmp, "RndSeed") == 0) {
                            X->RndSeed = (int) dtmp;
                        } else if (CheckWords(ctmp, "EpsSlater") == 0) {
                            X->eps_int_slater = (int) dtmp;
                        } else if (CheckWords(ctmp, "NMPTrans") == 0) {
                            X->NMPTrans = (int) dtmp;
                        } else {
                            fprintf(stdout, "  Warning: keyword \" %s \" is incorrect. \n", ctmp);
                        }
                    }
                    break;//modpara file

                case KWLocSpin:
                    cerr = ReadBuffInt(fp, &X->NLocSpn);
                    break;

                case KWTrans:
                    cerr = ReadBuffInt(fp, &X->NTransfer);
                    break;

                case KWCoulombIntra:
                    cerr = ReadBuffInt(fp, &X->NCoulombIntra);
                    break;

                case KWCoulombInter:
                    cerr = ReadBuffInt(fp, &X->NCoulombInter);
                    break;

                case KWHund:
                    cerr = ReadBuffInt(fp, &X->NHundCoupling);
                    break;

                case KWPairHop:
                    cerr = ReadBuffInt(fp, &X->NPairHopping);
                    X->NPairHopping *= 2;
                    break;

                case KWExchange:
                    cerr = ReadBuffInt(fp, &X->NExchangeCoupling);
                    break;

                case KWInterAll:
                    cerr = ReadBuffInt(fp, &X->NInterAll);
                    break;

                case KWOrbital:
                case KWOrbitalAntiParallel:
                    X->iFlgOrbital = 1;
                    cerr = ReadBuffInt(fp, &iNOrbitalAP);
                    X->NOrbitalAP = iNOrbitalAP;
                    iFlgOrbitalAP = 1;
                    X->NOrbitalIdx += iNOrbitalAP;
                    break;

                case KWOrbitalParallel:
                    X->iFlgOrbital = 1;
                    cerr = ReadBuffInt(fp, &iNOrbitalP);
                    X->NOrbitalP = 2 * iNOrbitalP;//up-up and down-down
                    iFlgOrbitalP = 1;
                    X->NOrbitalIdx += 2 * iNOrbitalP; //up-up and down-down
                    break;

                case KWOrbitalGeneral:
                    X->iFlgOrbital = 1;
                    cerr = ReadBuffInt(fp, &X->NOrbitalIdx);
                    X->iFlgOrbitalGeneral = 1;
                    break;

                case KWOneBodyG:
                    cerr = ReadBuffInt(fp, &X->NCisAjs);
                    break;

                case KWInitial:
                    cerr = ReadBuffInt(fp, &X->NInitial);
                    break;

                default:
                    break;
            }//case KW
            fclose(fp);
        }
        if (info != 0 || cerr == NULL) {
            fprintf(stderr, "error: Definition files(*.def) are incomplete.\n");
            fprintf(stdout, " Error:  Read File %s .\n", defname);
            return -1;
        }
    }

	if (X->iFlgOrbital == 1) {
        int iret = 0;
        iret = JudgeOrbitalMode(&X->iFlgOrbitalGeneral, iFlgOrbitalAP, iFlgOrbitalP, &X->OrbitalOutputMode);
        if (iret < 0) info = iret;
    }

  if (info != 0) {
		fprintf(stderr, "error: Definition files(*.def) are incomplete.\n");
		fprintf(stdout, " Error:  Read File %s .\n", defname);
		return -1;
	}

  if (X->Nsite <= 0) {
    fprintf(stderr, "Error: Nsite (in modpara.def) must be a positive integer.\n");
    return -1;
  }

  //CalcNCond
  if (X->Ncond == -1) {
    /* Ne is the number of electrons per spin, as in mVMC. */
    if (X->Ne == -1) {
      fprintf(stderr, "Error: Ne or NCond must be given in modpara.def.\n");
      return -1;
    }
    if (X->Ne < 0) {
      fprintf(stderr, "Error: Ne (in modpara.def) must be a non-negative integer.\n");
      return -1;
    }
    X->Ncond = 2 * X->Ne - X->NLocSpn;
    if (X->Ncond < 0) {
      fprintf(stderr, "Error: 2*Ne must be (2*Ne >= NLocalSpin).\n");
      return -1;
    }
  } else if (X->Ncond < 0) {
    fprintf(stderr, "Error: NCond (in modpara.def) must be a non-negative integer.\n");
    return -1;
  }
  if (X->iFlgOrbital == 1) {
    if (X->Ncond % 2 != 0) {
      fprintf(stderr, "Error: NCond (in modpara.def) must be even number.\n");
      return -1;
    }
    X->Ne = (X->NLocSpn + X->Ncond) / 2;
    X->Nsize = 2 * X->Ne;
  } else {
    X->Nsize = X->NLocSpn + X->Ncond;
    if (X->Ne == -1) X->Ne = X->Nsize / 2;
  }
  if (X->Nsize <= 0 || X->Nsize >= 2 * X->Nsite) {
    fprintf(stderr, "Error: the number of electrons NLocalSpin+NCond=%d must satisfy 0 < N < 2*Nsite (Nsite=%d).\n",
            X->Nsize, X->Nsite);
    return -1;
  }

  //CheckGeneral Orbital
  if(X->TwoSz !=0) {
      if (X->iFlgOrbital == 1) {
          if (X->iFlgOrbitalGeneral != 1) {
              fprintf(stderr, "Error:2Sz=%d\n", X->TwoSz);
              fprintf(stderr,
                      "Error: OrbitalParallel or OrbitalGeneral files must be needed when 2Sz !=0 (in modpara.def).\n");
              return -1;
          } else if (X->TwoSz % 2 != 0 && X->TwoSz != -1) {
              fprintf(stderr, "Error: 2Sz (in modpara.def) must be even number.\n");
              return -1;
          }
      }
  }

  if (X->NMPTrans < 0) {
		X->APFlag = 1; /* anti-periodic boundary */
		X->NMPTrans *= -1;
	} else {
		X->APFlag = 0;
	}
	X->fidx = 0;
	return 0;
}


int ReadDefFileIdxPara(
	char *xNameListFile, 
	struct DefineList *X
) {
  FILE *fp;
  char defname[D_FileNameMaxReadDef];
  char ctmp[D_FileNameMax], ctmp2[256];
  int iKWidx = 0;
  int i, j;
  int spn_i, spn_j, all_i, all_j;
  int fij, fijsgn, fij_org;
  int ndata;
  int idx, Orbitalidx, Orbitalsgn;
  int x0, x1, x2, x3;
  int ia[8];
  int info, i_spin, j_spin;
  double dReValue;
  double dImValue;
  char *cerr;
  info = 0;

  //Initialize OrbitalIdx
  // -1: not output
  for (i = 0; i < X->Nsite*2; i++) {
    for (j = 0; j < X->Nsite*2; j++) {
      X->OrbitalIdx[i][j] = -1;
      X->OrbitalSgn[i][j] = 0;
    }
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
    for (i = 0; i < IgnoreLinesInDef; i++){
      cerr = fgets(ctmp, sizeof(ctmp) / sizeof(char), fp);
      if(cerr==NULL) {
        fclose(fp);
        ReadDefFileError(defname);
        return -1;
      }
    }
    idx = 0;

    switch (iKWidx) {
      case KWLocSpin:
        /* Read locspn.def----------------------------------------*/
        while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
          if (CheckRecordCount(idx, X->Nsite, defname) != 0) {
            info = 1;
            break;
          }
          if (sscanf(ctmp2, "%d %d\n", &(x0), &(x1)) != 2) {
            info = ReportFormat(defname, idx);
            break;
          }
          if (CheckSite(x0, X->Nsite) != 0) {
            info = ReportSite(defname, idx);
            break;
          }
          if (x1 != 0 && x1 != 1) {
            fprintf(stderr, "Error in %s: the local-spin flag must be 0 or 1 at record %d.\n", defname, idx);
            info = 1;
            break;
          }
          X->LocSpn[x0] = x1;
          idx++;
        }
        if (X->Nsize < X->NLocSpn) {
          fprintf(stderr, "Error: 2*Ne must be (2*Ne >= NLocalSpin).\n");
          info = 1;
        }
        if (idx != X->Nsite) {
          info = ReadDefFileError(defname);
        }
        break;//locspn

      case KWTrans:
        /* transfer.def--------------------------------------*/
        if (X->NTransfer > 0) {
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (CheckRecordCount(idx, X->NTransfer, defname) != 0) {
              info = 1;
              break;
            }
            dImValue = 0;
            ndata = sscanf(ctmp2, "%d %d %d %d %lf %lf\n", &x0, &x1, &x2, &x3, &dReValue, &dImValue);
            if (ndata < 5) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckPairSite(x0, x2, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            if (CheckSpin(x1) != 0 || CheckSpin(x3) != 0) {
              info = ReportSpin(defname, idx);
              break;
            }
            X->Transfer[idx][0] = x0;
            X->Transfer[idx][1] = x1;
            X->Transfer[idx][2] = x2;
            X->Transfer[idx][3] = x3;
            X->ParaTransfer[idx] = dReValue + dImValue * I;
            idx++;
          }
          if (idx != X->NTransfer) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWCoulombIntra:
        /*coulombintra.def----------------------------------*/
        if (X->NCoulombIntra > 0) {
          idx = 0;
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (CheckRecordCount(idx, X->NCoulombIntra, defname) != 0) {
              info = 1;
              break;
            }
            if (sscanf(ctmp2, "%d %lf\n", &x0, &dReValue) != 2) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckSite(x0, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            X->CoulombIntra[idx][0] = x0;
            X->ParaCoulombIntra[idx] = dReValue;
            idx++;
          }
          if (idx != X->NCoulombIntra) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWCoulombInter:
        /*coulombinter.def----------------------------------*/
        if (X->NCoulombInter > 0) {
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (CheckRecordCount(idx, X->NCoulombInter, defname) != 0) {
              info = 1;
              break;
            }
            if (sscanf(ctmp2, "%d %d %lf\n", &x0, &x1, &dReValue) != 3) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckPairSite(x0, x1, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            X->CoulombInter[idx][0] = x0;
            X->CoulombInter[idx][1] = x1;
            X->ParaCoulombInter[idx] = dReValue;
            idx++;
          }
          if (idx != X->NCoulombInter) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWHund:
        /*hund.def------------------------------------------*/
        if (X->NHundCoupling > 0) {
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (CheckRecordCount(idx, X->NHundCoupling, defname) != 0) {
              info = 1;
              break;
            }
            if (sscanf(ctmp2, "%d %d %lf\n", &x0, &x1, &dReValue) != 3) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckPairSite(x0, x1, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            X->HundCoupling[idx][0] = x0;
            X->HundCoupling[idx][1] = x1;
            X->ParaHundCoupling[idx] = dReValue;
            idx++;
          }
          if (idx != X->NHundCoupling) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWPairHop:
        /*pairhop.def---------------------------------------*/
        if (X->NPairHopping > 0) {
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (CheckRecordCount(idx, X->NPairHopping / 2, defname) != 0) {
              info = 1;
              break;
            }
            if (sscanf(ctmp2, "%d %d %lf\n", &x0, &x1, &dReValue) != 3) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckPairSite(x0, x1, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            X->PairHopping[2*idx][0] = x0;
            X->PairHopping[2*idx][1] = x1;
            X->ParaPairHopping[2*idx] = dReValue;
            X->PairHopping[2*idx+1][0] = x1;
            X->PairHopping[2*idx+1][1] = x0;
            X->ParaPairHopping[2*idx+1] = dReValue;
            idx++;
          }
          if (idx != X->NPairHopping/2) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWExchange:
        /*exchange.def--------------------------------------*/
        if (X->NExchangeCoupling > 0) {
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (CheckRecordCount(idx, X->NExchangeCoupling, defname) != 0) {
              info = 1;
              break;
            }
            if (sscanf(ctmp2, "%d %d %lf\n", &x0, &x1, &dReValue) != 3) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckPairSite(x0, x1, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            X->ExchangeCoupling[idx][0] = x0;
            X->ExchangeCoupling[idx][1] = x1;
            X->ParaExchangeCoupling[idx] = dReValue;
            idx++;
          }
          if (idx != X->NExchangeCoupling) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWInterAll:
        X->iFlg_Fock = 1;
        /*inteall.def--------------------------------------*/
        if (X->NInterAll > 0) {
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (CheckRecordCount(idx, X->NInterAll, defname) != 0) {
              info = 1;
              break;
            }
            dImValue = 0;
            ndata = sscanf(ctmp2, "%d %d %d %d %d %d %d %d %lf %lf\n",
                           &ia[0], &ia[1], &ia[2], &ia[3], &ia[4], &ia[5], &ia[6], &ia[7],
                           &dReValue, &dImValue);
            if (ndata < 9) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckQuadSite(ia[0], ia[2], ia[4], ia[6], X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            if (CheckSpin(ia[1]) != 0 || CheckSpin(ia[3]) != 0 ||
                CheckSpin(ia[5]) != 0 || CheckSpin(ia[7]) != 0) {
              info = ReportSpin(defname, idx);
              break;
            }
            for (i = 0; i < 8; i++) {
              X->InterAll[idx][i] = ia[i];
            }
            X->ParaInterAll[idx] = dReValue + dImValue * I;
            idx++;
          }
          if (idx != X->NInterAll) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWOrbital:
      case KWOrbitalAntiParallel:
        idx = 0;
        if(X->iFlgOrbitalGeneral==0) {
          while( fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp) != NULL){
            ndata = sscanf(ctmp2, "%d %d %d %d\n",
                   &i,
                   &j,
                   &Orbitalidx,
                   &Orbitalsgn);
            if (ndata < 3) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (ndata == 3){ //sgn does not exist in input file
              Orbitalsgn = 1;
            }
            if(CheckPairSite(i, j, X->Nsite) != 0){
              info = ReportSite(defname, idx);
              break;
            }
            if (Orbitalidx < 0 || Orbitalidx >= X->NOrbitalAP) {
              info = ReportOrbitalIdx(defname, idx, Orbitalidx, X->NOrbitalAP);
              break;
            }
            X->OrbitalIdx[i+X->Nsite*0][j+X->Nsite*1]=Orbitalidx;
            X->OrbitalSgn[i+X->Nsite*0][j+X->Nsite*1]=Orbitalsgn;
            idx++;
            if(idx==X->Nsite*X->Nsite) break;
          }
          if(idx!=X->Nsite*X->Nsite) {
            info=ReadDefFileError(defname);
          }
        }else{//general orbital
          while( fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp) != NULL){
            ndata = sscanf(ctmp2, "%d %d %d %d\n", &i, &j, &fij, &fijsgn);
            if (ndata < 3) {
              info = ReportFormat(defname, idx);
              break;
            }
            spn_i = 0;
            spn_j = 1;
            if (ndata == 3){ //sgn does not exist in input file
              fijsgn = 1;
            }
            if(CheckPairSite(i, j, X->Nsite) != 0){
              info = ReportSite(defname, idx);
              break;
            }
            if (fij < 0 || fij >= X->NOrbitalAP) {
              info = ReportOrbitalIdx(defname, idx, fij, X->NOrbitalAP);
              break;
            }
            all_i = i+spn_i*X->Nsite; //fsz
            all_j = j+spn_j*X->Nsite; //fsz
            idx++;
            X->OrbitalIdx[all_i][all_j]=fij;
            X->OrbitalSgn[all_i][all_j]=fijsgn;
            // Note F_{IJ}=-F_{JI}
            X->OrbitalIdx[all_j][all_i]=fij;
            X->OrbitalSgn[all_j][all_i]=-fijsgn;
            if(idx==(X->Nsite*X->Nsite)) break;
          }
        }
        break;

      case KWOrbitalParallel:
        /*orbitalidxt.def------------------------------------*/
          idx = 0;
            while( fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp) != NULL){
              ndata = sscanf(ctmp2, "%d %d %d %d\n", &i, &j, &fij_org, &fijsgn);
              if (ndata < 3) {
                info = ReportFormat(defname, idx);
                break;
              }
              if (ndata == 3){ //sgn does not exist in input file
                fijsgn = 1;
              }
              if (CheckPairSite(i, j, X->Nsite) != 0) {
                info = ReportSite(defname, idx);
                break;
              }
              if (fij_org < 0 || 2 * fij_org + 1 >= X->NOrbitalP) {
                info = ReportOrbitalIdx(defname, idx, fij_org, X->NOrbitalP / 2);
                break;
              }
              for (spn_i = 0; spn_i < 2; spn_i++) {
                all_i = i + spn_i * X->Nsite; //fsz
                all_j = j + spn_i * X->Nsite; //fsz
                idx++;
                fij = X->NOrbitalAP + 2 * fij_org + spn_i;
                X->OrbitalIdx[all_i][all_j] = fij;
                X->OrbitalSgn[all_i][all_j] = fijsgn;
                X->OrbitalIdx[all_j][all_i] = fij;
                X->OrbitalSgn[all_j][all_i] = -fijsgn;
              }
              if (idx == (X->Nsite * (X->Nsite - 1))) {
                break;
              }
            }

        break;

      case KWOrbitalGeneral:
        /*orbitalidx.def------------------------------------*/
        if (X->NOrbitalIdx > 0) {
          idx = 0;
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            //TODO: Replace for spin dependent
            ndata = sscanf(ctmp2, "%d %d %d %d %d %d\n",
                   &i,
                   &i_spin,
                   &j,
                   &j_spin,
                   &Orbitalidx,
                   &Orbitalsgn);
            if (ndata < 5) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (ndata == 5){ //sgn does not exist in input file
              Orbitalsgn = 1;
            }
            if (CheckPairSite(i, j, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            if (CheckSpin(i_spin) != 0 || CheckSpin(j_spin) != 0) {
              info = ReportSpin(defname, idx);
              break;
            }
            if (Orbitalidx < 0 || Orbitalidx >= X->NOrbitalIdx) {
              info = ReportOrbitalIdx(defname, idx, Orbitalidx, X->NOrbitalIdx);
              break;
            }
            (X->OrbitalIdx[i + X->Nsite * i_spin][j + X->Nsite * j_spin]) = Orbitalidx;
            (X->OrbitalSgn[i + X->Nsite * i_spin][j + X->Nsite * j_spin]) = Orbitalsgn;
            idx++;
            if (idx == X->Nsite * (2 * X->Nsite - 1)) break;
          }
          if (idx != X->Nsite * (2 * X->Nsite - 1)) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWOneBodyG:
        /*cisajs.def----------------------------------------*/
        if (X->NCisAjs > 0) {
          idx = 0;
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            if (CheckRecordCount(idx, X->NCisAjs, defname) != 0) {
              info = 1;
              break;
            }
            if (sscanf(ctmp2, "%d %d %d %d\n", &(x0), &(x1), &(x2), &(x3)) != 4) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckPairSite(x0, x2, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            if (CheckSpin(x1) != 0 || CheckSpin(x3) != 0) {
              info = ReportSpin(defname, idx);
              break;
            }
            X->CisAjs[idx][0] = x0;
            X->CisAjs[idx][1] = x1;
            X->CisAjs[idx][2] = x2;
            X->CisAjs[idx][3] = x3;
            if (x1 != x3) {
              fprintf(stderr, "  Error:  Sz non-conserved system is not yet supported in mVMC ver.1.0.\n");
              info = ReadDefFileError(defname);
              break;
            }
            idx++;
          }
          if (idx != X->NCisAjs) {
            info = ReadDefFileError(defname);
          }
        }
        break;

      case KWInitial:
        /*initial.def------------------------------------*/
        if (X->NInitial > 0) {
          idx = 0;
          while (fgets(ctmp2, sizeof(ctmp2) / sizeof(char), fp) != NULL) {
            /* As before, records beyond NInitial are ignored. */
            if (idx >= X->NInitial) break;
            dImValue = 0;
            ndata = sscanf(ctmp2, "%d %d %d %d %lf %lf\n", &x0, &x1, &x2, &x3, &dReValue, &dImValue);
            if (ndata < 5) {
              info = ReportFormat(defname, idx);
              break;
            }
            if (CheckPairSite(x0, x2, X->Nsite) != 0) {
              info = ReportSite(defname, idx);
              break;
            }
            if (CheckSpin(x1) != 0 || CheckSpin(x3) != 0) {
              info = ReportSpin(defname, idx);
              break;
            }
            X->Initial[idx][0] = x0;
            X->Initial[idx][1] = x1;
            X->Initial[idx][2] = x2;
            X->Initial[idx][3] = x3;
            X->ParaInitial[idx] = dReValue + dImValue * I;
            idx++;
          }
        }

        if (idx != X->NInitial) {
          info = ReadDefFileError(defname);
        } else {
          //	 info=ReadDefFileError(xNameListFile);
        }
        break;

      default:
        fprintf(stdout, "!! Warning: %s is not used for Hatree Fock Calculation. !!\n", defname);
        break;
    }
    fclose(fp);
  }

  if (info != 0) {
    fprintf(stderr, "error: Indices and Parameters of Definition files(*.def) are incomplete.\n");
    return -1;
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
	FILE *fplist;
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
        /* strlen on the uninitialized buffers was undefined behavior;
           clearing the first byte is enough for the empty-string checks. */
        ctmpKW[0] = '\0';
        ctmpFileName[0] = '\0';
        sscanf(ctmp2,"%s %s\n", ctmpKW, ctmpFileName);
		if(strncmp(ctmpKW, "#", 1)==0 || *ctmp2=='\n' || (strcmp(ctmpKW, "")&&strcmp(ctmpFileName,""))==0)
        {
			continue;
		}
		else if(strcmp(ctmpKW, "")*strcmp(ctmpFileName, "")==0){
			fprintf(stderr, "Error: keyword and filename must be set as a pair in %s.\n", cFileListNameFile);
			fclose(fplist);
			return(-1);
		}
		/*!< Check KW */
		if( CheckKW(ctmpKW, cKWListOfFileNameList, KWIdxInt_end, &itmpKWidx)!=0 ){
			fprintf(stderr, "Warning: Wrong keywords '%s' in %s.\n", ctmpKW, cFileListNameFile);
			/*
			fprintf(stderr, "%s", "Choose Keywords as follows: \n");
			for(i=0; i<KWIdxInt_end;i++){
				fprintf(stderr, "%s \n", cKWListOfFileNameList[i]);
			}
			fclose(fplist);
			return(-1);
			*/
			/* itmpKWidx is -1 here; skip to avoid writing cFileNameList[-1] */
			continue;
        }
        else {
            /*!< Check cFileNameList to prevent from double registering the file name */
            if (strcmp(cFileNameList[itmpKWidx], "") != 0) {
                fprintf(stderr, "Error: Same keywords exist in %s.\n", cFileListNameFile);
                fclose(fplist);
                return (-1);
            }
        }
		/*!< Copy FileName */
		strcpy(cFileNameList[itmpKWidx], ctmpFileName);
	}
	fclose(fplist);
	return 0;
}

void SetInitialValue(struct DefineList *X){
  X->NTransfer=0;
  X->NCoulombIntra=0;
  X->NCoulombInter=0;
  X->NHundCoupling=0;
  X->NPairHopping=0;
  X->NExchangeCoupling=0;
  X->NOrbitalIdx=0;
  X->NCisAjs=0;
  X->NInitial=0;
  X->mix=0.5;
  X->eps_int=10;
  X->print=0;
  X->IterationMax=2000;
  X->eps_int_slater=6;
  X->Ncond=-1;
  X->Ne=-1;
  X->Nsize=0;
  X->NOrbitalAP=0;
  X->NOrbitalP=0;
  X->iFlgOrbitalGeneral=0;
  X->OrbitalOutputMode=0;
  X->TwoSz=-1;
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
)
{
  if(iSite>=iMaxNum || iSite < 0)  return(-1);
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
)
{
  if(CheckSite(iSite1, iMaxNum)!=0){
    return(-1);
  }
  if(CheckSite(iSite2, iMaxNum)!=0){
    return(-1);
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
)
{
  if(CheckPairSite(iSite1, iSite2, iMaxNum)!=0){
    return(-1);
  }
  if(CheckPairSite(iSite3, iSite4, iMaxNum)!=0){
    return(-1);
  }
  return 0;
}
