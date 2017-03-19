#include "physcal_lanczos.h"

void PhysCalLanczos_real
(
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
 )
{
  int i=0;
  /* zvo_ls.dat */
  fprintf(_FileLS, "% .18e  ", _QQQQ_real[2]);  /* H * I = QQQQ[1],[2],[4],[8] */      //TBC
  fprintf(_FileLS, "% .18e  ", _QQQQ_real[3]);  /* H * H = QQQQ[3],[6],[9],[12] */     //TBC
  fprintf(_FileLS, "% .18e  ", _QQQQ_real[10]); /* H^2 * I = QQQQ[5],[10] */           //TBC
  fprintf(_FileLS, "% .18e  ", _QQQQ_real[11]); /* H^2 * H = QQQQ[7],[11],[13],[14] */ //TBC
  fprintf(_FileLS, "% .18e\n", _QQQQ_real[15]); /* H^2 * H^2 = QQQQ[15] */             //TBC

  /* zvo_ls_qqqq.dat */
  for (i = 0; i < _nLSHam * _nLSHam * _nLSHam * _nLSHam; i++) {
    fprintf(_FileLSQQQQ, "% .18e  ", _QQQQ_real[i]);
  }
  fprintf(_FileLSQQQQ, "\n");

  if (_NLanczosmode > 1) {
    /* zvo_ls_qcisajsq.dat */
    for (i = 0; i < _nLSHam * _nLSHam * _nCisAjs; i++) {
      fprintf(_FileLSQCisAjsQ, "% .18e  ", _QCisAjsQ_real[i]);
    }
    fprintf(_FileLSQCisAjsQ, "\n");

    /* zvo_ls_qcisajscktaltq.dat */
    for (i = 0; i < _nLSHam * _nLSHam * _nCisAjsCktAlt; i++) {
      fprintf(_FileLSQCisAjsCktAltQ, "% .18e  ", _QCisAjsCktAltQ_real[i]);
    }
    fprintf(_FileLSQCisAjsCktAltQ, "\n");
  }
}

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
)
{
  int i=0;
  /* zvo_ls.dat */
  fprintf(_FileLS, "% .18e  ", creal(_QQQQ[2]));  /* H * I = _QQQQ[1],[2],[4],[8] */      //TBC
  fprintf(_FileLS, "% .18e  ", creal(_QQQQ[3]));  /* H * H = _QQQQ[3],[6],[9],[12] */     //TBC
  fprintf(_FileLS, "% .18e  ", creal(_QQQQ[10])); /* H^2 * I = _QQQQ[5],[10] */           //TBC
  fprintf(_FileLS, "% .18e  ", creal(_QQQQ[11])); /* H^2 * H = _QQQQ[7],[11],[13],[14] */ //TBC
  fprintf(_FileLS, "% .18e\n", creal(_QQQQ[15])); /* H^2 * H^2 = _QQQQ[15] */             //TBC

  /* zvo_ls_qqqq.dat */
  for (i = 0; i < _nLSHam * _nLSHam * _nLSHam * _nLSHam; i++) {
    fprintf(_FileLSQQQQ, "% .18e  ", creal(_QQQQ[i]));
  }
  fprintf(_FileLSQQQQ, "\n");

  if (_NLanczosmode > 1) {
    /* zvo_ls_qcisajsq.dat */
    for (i = 0; i < _nLSHam * _nLSHam * _nCisAjs; i++) {
      fprintf(_FileLSQCisAjsQ, "% .18e  ", creal(_QCisAjsQ[i]));
    }
    fprintf(_FileLSQCisAjsQ, "\n");

    /* zvo_ls_qcisajscktaltq.dat */
    for (i = 0; i < _nLSHam * _nLSHam * _nCisAjsCktAlt; i++) {
      fprintf(_FileLSQCisAjsCktAltQ, "% .18e  ",  creal(_QCisAjsCktAltQ[i]));
    }
    fprintf(_FileLSQCisAjsCktAltQ, "\n");
  }
}
