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

int ReadDefFileError(
	const char *defname
){
	printf("Error: %s (Broken file or Not exist)\n", defname);
	return 0;
}


int ReadDefFileNInt(
	char *xNameListFile, 
	struct DefineList *X
){
	FILE *fp, *fplist;
	char defname[D_FileNameMaxReadDef];
	char ctmp[D_CharTmpReadDef];
	char ctmp2[D_FileNameMax];
	int itmp, info;
	int iKWidx=0;

	cFileNameListFile = malloc(sizeof(char)*D_CharTmpReadDef*KWIdxInt_end);
	fprintf(stdout, "  Read File %s .\n", xNameListFile);
	if(GetFileName(xNameListFile, cFileNameListFile)!=0){
		fprintf(stderr, "  error: Definition files(*.def) are incomplete.\n");
		return(-1);
	}

	for(iKWidx=0; iKWidx< KWIdxInt_end; iKWidx++){
		strcpy(defname, cFileNameListFile[iKWidx]);
		if(strcmp(defname,"")==0){
			switch (iKWidx){
				case KWModPara:
				case KWLocSpin:
					fprintf(stderr, "  Error: Need to make a def file for %s.\n", cKWListOfFileNameList[iKWidx]);
					return (-1);
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
					//SetDefultValuesModPara(X);
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &itmp); //2
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //3
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //4
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp); //5
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %s\n", ctmp, X->CDataFileHead); //6
					fgets(ctmp2,sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %s\n", ctmp, X->CParaFileHead); //7
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);   //8

					double dtmp;
					while(fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp)!=NULL){
						if(*ctmp2 == '\n' || ctmp2[0] == '-')  continue;
						sscanf(ctmp2,"%s %lf\n", ctmp, &dtmp);
						if(CheckWords(ctmp, "Nsite")==0){
							X->Nsite=(int)dtmp;
						}
						else if(CheckWords(ctmp, "Ne")==0 || CheckWords(ctmp, "Nelectron")==0 ){
							X->Ne=(int)dtmp;
						}
						else if(CheckWords(ctmp, "Mix")==0){
							X->mix=(int)dtmp;
						}
						else if(CheckWords(ctmp, "EPS")==0){
							X->eps_int=(int)dtmp;
						}
						else if(CheckWords(ctmp, "Print")==0){
							X->print=(int)dtmp;
						}
						else if(CheckWords(ctmp, "IterationMax")==0){
							X->IterationMax=(int)dtmp;
						}
						else{
							fprintf(stderr, "  Error: keyword \" %s \" is incorrect. \n", ctmp);
							info = ReadDefFileError(defname);
						}
					}
					fclose(fp);
					break;//modpara file

				case KWLocSpin:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &(X->NLocSpn));
					fclose(fp);
					break;

				case KWTrans:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp,  &(X->NTransfer));
					fclose(fp);
					break;

				case KWCoulombIntra:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &(X->NCoulombIntra));
					fclose(fp);
					break;

				case KWCoulombInter:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &(X->NCoulombInter));
					fclose(fp);
					break;

				case KWHund:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &(X->NHundCoupling));
					fclose(fp);
					break;

				case KWPairHop:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &(X->NPairHopping));
					fclose(fp);
					break;

				case KWExchange:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &(X->NExchangeCoupling));
					fclose(fp);
					break;

				case KWOneBodyG:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &(X->NCisAjs));
					fclose(fp);
					break;

/*
				case KWInterAll:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &bufInt[IdxNInterAll]);
					fclose(fp);
					break;
*/
				case KWInitial:
					fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
					fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp);
					sscanf(ctmp2,"%s %d\n", ctmp, &(X->NInitial));
					fclose(fp);
					break;

				default:
					fclose(fp);
					break;
			}//case KW
		}
	}

	if(info!=0) {
		fprintf(stderr, "error: Definition files(*.def) are incomplete.\n");
		return -1;
	}

	X->Nsize   = 2*X->Ne;
	X->fidx = 0;
	return 1;
}


int ReadDefFileIdxPara(
	char *xNameListFile, 
	struct DefineList *X
){
	FILE *fp, *fplist;
	char defname[D_FileNameMaxReadDef];
	char ctmp[D_FileNameMax], ctmp2[256];
	int iKWidx=0;
	int i, idx;
	int x0,x1,x2,x3,x4,x5,x6,x7;
	int itmp, info;
	double dReValue;
	double dImValue;

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
					X->LocSpn[x0] = x1;
					idx++;
				}
				if(2 * X->Ne < X->NLocSpn){
					fprintf(stderr, "Error: 2*Ne must be (2*Ne >= NLocalSpin).\n");
					info=1;
				}
				if(idx!=X->Nsite) info = ReadDefFileError(defname);
				fclose(fp);
				break;//locspn

			case KWTrans:
				/* transfer.def--------------------------------------*/
				if(X->NTransfer>0){
					while( fscanf(fp, "%d %d %d %d %lf %lf\n",
								  &(X->Transfer[idx][0]),
								  &(X->Transfer[idx][1]),
								  &(X->Transfer[idx][2]),
								  &(X->Transfer[idx][3]),
								  &dReValue,
								  &dImValue)!=EOF){
//						X->ParaTransfer[idx]=dReValue+dImValue*I;
						X->ParaTransfer[idx]=dReValue;
						if(X->Transfer[idx][1] != X->Transfer[idx][3]){
							fprintf(stderr, "  Error:  Sz non-conserved system is not yet supported in mVMC ver.1.0.\n");
							info = ReadDefFileError(defname);
							break;
						}
						idx++;
					}
					if(idx!=X->NTransfer) info = ReadDefFileError(defname);
				}
				fclose(fp);
				break;

			case KWCoulombIntra:
				/*coulombintra.def----------------------------------*/
				if(X->NCoulombIntra>0){
					while( fscanf(fp, "%d %lf\n",
								  &(X->CoulombIntra[idx][0]),
								  &(X->ParaCoulombIntra[idx]) )!=EOF){
						idx++;
					}
					if(idx!=X->NCoulombIntra) info = ReadDefFileError(defname);
				}
				fclose(fp);
				break;

			case KWCoulombInter:
				/*coulombinter.def----------------------------------*/
				if(X->NCoulombInter>0){
					while(fgets(ctmp2, sizeof(ctmp2)/sizeof(char), fp) != NULL){
						sscanf(ctmp2, "%d %d %lf\n",
							   &(X->CoulombInter[idx][0]),
							   &(X->CoulombInter[idx][1]),
							   &(X->ParaCoulombInter[idx])
						);
						idx++;
					}
					if(idx!=X->NCoulombInter) info=ReadDefFileError(defname);
				}
				fclose(fp);
				break;

			case KWHund:
				/*hund.def------------------------------------------*/
				if(X->NHundCoupling>0){
					while( fscanf(fp, "%d %d %lf\n",
								  &(X->HundCoupling[idx][0]),
								  &(X->HundCoupling[idx][1]),
								  &(X->ParaHundCoupling[idx]) )!=EOF){
						idx++;
					}
					if(idx!=X->NHundCoupling) info=ReadDefFileError(defname);
				}
				fclose(fp);
				break;

			case KWPairHop:
				/*pairhop.def---------------------------------------*/
				if(X->NPairHopping>0){
					while( fscanf(fp, "%d %d %lf\n",
								  &(X->PairHopping[idx][0]),
								  &(X->PairHopping[idx][1]),
								  &(X->ParaPairHopping[idx]) )!=EOF){
						idx++;
					}
					if(idx!=X->NPairHopping) info=ReadDefFileError(defname);
				}
				fclose(fp);
				break;

			case KWExchange:
				/*exchange.def--------------------------------------*/
				if(X->NExchangeCoupling>0){
					while( fscanf(fp, "%d %d %lf\n",
								  &(X->ExchangeCoupling[idx][0]),
								  &(X->ExchangeCoupling[idx][1]),
								  &(X->ParaExchangeCoupling[idx]) )!=EOF){
						idx++;
					}
					if(idx!=X->NExchangeCoupling) info=ReadDefFileError(defname);
				}
				fclose(fp);
				break;

			case KWOneBodyG:
				/*cisajs.def----------------------------------------*/
				if(X->NCisAjs>0){
					idx = 0;
					while( fscanf(fp, "%d %d %d %d\n",
								  &(x0), &(x1), &(x2), &(x3)) != EOF){
						X->CisAjs[idx][0] = x0;
						X->CisAjs[idx][1] = x1;
						X->CisAjs[idx][2] = x2;
						X->CisAjs[idx][3] = x3;
						if(x1 != x3){
							fprintf(stderr, "  Error:  Sz non-conserved system is not yet supported in mVMC ver.1.0.\n");
							info = ReadDefFileError(defname);
							break;
						}
						idx++;
					}
					if(idx!=X->NCisAjs) info=ReadDefFileError(defname);
				}
				fclose(fp);
				break;

				/*
			case KWInterAll:
				//
				if(X->NInterAll>0){
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
					// info=ReadDefFileError(xNameListFile);
				}
				fclose(fp);
				break;
		*/
			case KWInitial:
				/*initial.def------------------------------------*/
				if(X->NInitial>0){
					idx = 0;
					while( fscanf(fp, "%d %d %d %d %lf %lf\n",
								  &(X->Initial[idx][0]),
								  &(X->Initial[idx][1]),
								  &(X->Initial[idx][2]),
								  &(X->Initial[idx][3]),
								  &(X->ParaInitial[idx]),
								  &(X->ParaInitial_theta[idx])
					)!=EOF
							){
						idx++;
					}

					if(idx!=X->NInitial) info=ReadDefFileError(defname);
				} else {
					// info=ReadDefFileError(xNameListFile);
				}
				fclose(fp);
				break;

			default:
				fclose(fp);
				break;
		}
	}


	if(info!=0) {
		fprintf(stderr, "error: Indices and Parameters of Definition files(*.def) are incomplete.\n");
		return -1;
	}

	return 1;
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
