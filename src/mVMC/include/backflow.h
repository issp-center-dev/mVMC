#ifndef _BACKFLOW
#define _BACKFLOW

#include <stdio.h>

int BFComputeSizes(int nBackFlowIdx, int nrange, int nzBF,
                   int *nrangeIdx, int *nBFIdxTotal, int *nProjBF);
int BFValidateSettings(int hasBF, int hasBFRange, int backflowSupported);
int BFDefIntCount(void);
int BFWorkIntCount(void);
int BFReadRange(FILE *fp, const char *defname);
int BFReadDefinition(FILE *fp, int *optFlag, int *countIdx, const char *defname);
void BFBindDefTables(int **pInt);
void BFFreeDefTables(void);
void BFAllocRuntime(void);
void BFFreeRuntime(void);
void BFInitParameters(void);
void BFSetupIndex(void);
void BFRefreshRealLookupTables(void);

#endif
