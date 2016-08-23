#pragma once
#include "Def.h"

#define D_FileNameMax 256
#define D_FileNameMaxReadDef 200
#define D_CharTmpReadDef     200

#define IgnoreLinesInDef 5

/**
 * Keyword List in NameListFile.
 **/
static char cKWListOfFileNameList[][D_CharTmpReadDef]={
		"ModPara", "LocSpin",
		"Trans", "CoulombIntra", "CoulombInter",
		"Hund", "PairHop", "Exchange",
		"OneBodyG", "InterAll",
		"Initial"
};

/**
 * Number of Keyword List in NameListFile for this program.
 **/
enum KWIdxInt{
	KWModPara, KWLocSpin,
	KWTrans, KWCoulombIntra,KWCoulombInter,
	KWHund, KWPairHop, KWExchange,
	KWOneBodyG, KWInterAll,
	KWInitial,KWIdxInt_end
};

/**
 * File Name List in NameListFile.
 **/
static char (*cFileNameListFile)[D_CharTmpReadDef];

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
