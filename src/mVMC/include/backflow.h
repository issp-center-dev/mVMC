#ifndef _BACKFLOW
#define _BACKFLOW

int BFComputeSizes(int nBackFlowIdx, int nrange, int nzBF,
                   int *nrangeIdx, int *nBFIdxTotal, int *nProjBF);
int BFValidateSettings(int hasBF, int hasBFRange, int backflowSupported);
int BFDefIntCount(void);
int BFWorkIntCount(void);
void BFBindDefTables(int **pInt);
void BFFreeDefTables(void);
void BFAllocRuntime(void);
void BFFreeRuntime(void);
void BFSetupIndex(void);

#endif
