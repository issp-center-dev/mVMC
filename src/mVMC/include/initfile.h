#ifndef _INITFILE
#define _INITFILE

void InitFile(char *xNameListFile, int rank);
void InitFilePhysCal(int i, int rank);
void CloseFile(int rank);
void CloseFilePhysCal(int rank);
void FlushFile(int step, int rank);
void writeConfig(char *xNameFile, char *fileName);
int fileCopyAdd(char *inputFileName, FILE *outputFile);

#endif
