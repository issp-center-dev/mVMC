#pragma once
#include "Def.h"

int ReadDefFileError(
	char *defname
                     );


int ReadDefFileNInt(
	char *xNameListFile, 
	struct DefineList *X
                    );

int ReadDefFileIdxPara(
	char *xNameListFile, 
	struct DefineList *X
                       );
