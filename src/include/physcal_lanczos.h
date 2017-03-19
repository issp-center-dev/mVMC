#ifndef _PHYSCAL_LANCZOS
#define _PHYSCAL_LANCZOS
#include <complex.h>
#include <stdio.h>

void PhysCalLanczos_real(
 double *_QQQQ_real,
 double *_QCisAjsQ_real,
 double *_QCisAjsCktAltQ_real,
 const int _nLSHam,
 const int _nCisAjs,
 const int _nCisAjsCktAlt,
 const int _NLanczosmode,
 FILE *_FileLS,
 FILE *_FileLSQQQQ,
 FILE *_FileLSQCisAjsQ,
 FILE *_FileLSQCisAjsCktAltQ
);

void PhysCalLanczos_fcmp(
 double complex* _QQQQ,
 double complex* _QCisAjsQ,
 double complex* _QCisAjsCktAltQ,
 const int _nLSHam,
 const int _nCisAjs,
 const int _nCisAjsCktAlt,
 const int _NLanczosmode,
 FILE *_FileLS,
 FILE *_FileLSQQQQ,
 FILE *_FileLSQCisAjsQ,
 FILE *_FileLSQCisAjsCktAltQ
);
#endif


