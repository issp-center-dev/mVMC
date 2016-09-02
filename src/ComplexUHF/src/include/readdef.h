#pragma once
#include "Def.h"

#define D_FileNameMax 256
#define D_FileNameMaxReadDef 200
#define D_CharTmpReadDef     200

#define IgnoreLinesInDef 5

/**
 * File Name List in NameListFile.
 **/

int CheckWords( const char* ctmp, const char* cKeyWord);
int CheckKW(const char* cKW, char  cKWList[][D_CharTmpReadDef], int iSizeOfKWidx, int* iKWidx);
int GetKWWithIdx(char *ctmpLine, char *ctmp, int *itmp);
int ValidateValue(const int icheckValue, const int ilowestValue, const int iHighestValue);
int GetFileName(const char* cFileListNameFile, char cFileNameList[][D_CharTmpReadDef]);

void SetDefultValuesModPara(struct DefineList *X);

int ReadDefFileNInt(
	char *xNameListFile, 
	struct DefineList *X
                    );

int ReadDefFileIdxPara(
	char *xNameListFile, 
	struct DefineList *X
                       );

int trim(char *s);