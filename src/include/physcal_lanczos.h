#ifndef _PHYSCAL_LANCZOS
#define _PHYSCAL_LANCZOS
#include <complex.h>
#include <stdio.h>

int PhysCalLanczos_real(
 double *_QQQQ_real,
 double *_QCisAjsQ_real,
 double *_QCisAjsCktAltQ_real,
 const int _nLSHam,
 const int _Ns,
 const int _nCisAjs,
 const int _nCisAjsLz,
 int **_iOneBodyGIdx,
 int **_CisAjsLzIdx,
 const int _nCisAjsCktAlt,
 const int _NLanczosmode,
 FILE *_FileLS,
 FILE *_FileLSQQQQ,
 FILE *_FileLSQCisAjsQ,
 FILE *_FileLSQCisAjsCktAltQ,
 FILE *_FileLSCisAjs,
 FILE *_FileLSCisAjsCktAlt
);

int PhysCalLanczos_fcmp(
 double complex* _QQQQ,
 double complex* _QCisAjsQ,
 double complex* _QCisAjsCktAltQ,
 const int _nLSHam,
 const int _Ns,
 const int _nCisAjs,
 const int _nCisAjsLz,
 int **_iOneBodyGIdx,
 int **_CisAjsLzIdx,
 const int _nCisAjsCktAlt,
 const int _NLanczosmode,
 FILE *_FileLS,
 FILE *_FileLSQQQQ,
 FILE *_FileLSQCisAjsQ,
 FILE *_FileLSQCisAjsCktAltQ,
 FILE *_FileLSCisAjs,
 FILE *_FileLSCisAjsCktAlt
);
#endif


