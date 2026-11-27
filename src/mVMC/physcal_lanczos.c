#include "physcal_lanczos.h"
#include <float.h>
#include "math.h"
#include "stdlib.h"

int CalculateEneByAlpha(double H1, double H2_1, double H2_2,
            double H3, double H4, double alpha,
            double *ene, double*ene_V);

int CalculatePhysVal_real(
    double H1, double H2_1, double alpha,
    double *_QPhysQ_real, int NPhys, int NLSHam,
    double *_Phys_LS_real
);

int CalculatePhysVal_fcmp(
    double complex H1, double complex H2_1, double complex alpha,
    double complex *_QPhysQ_fcmp, int NPhys,int NLSHam,
    double complex *_Phys_LS_fcmp
);


int PhysCalLanczos_real
(
 double *_QQQQ_real,
 double *_QCisAjsQ_real,
 double *_QCisAjsCktAltQ_real,
 double *_QCisAjsCktAltQDC_real,
 const int _nLSHam,
 const int _Ns,
 const int _nCisAjs,
 const int _nCisAjsLz,
 int **_iOneBodyGIdx,
 int **_CisAjsLzIdx,
 const int _nCisAjsCktAlt,
 const int _nCisAjsCktAltDC,
 int **_CisAjsCktAltDC,
 const int _NLanczosmode,
 size_t sampleCount,
 FILE *_FileLS,
 FILE *_FileLSQQQQ,
 FILE *_FileLSQCisAjsQ,
 FILE *_FileLSQCisAjsCktAltQ,
 FILE *_FileLSCisAjs,
 FILE *_FileLSCisAjsCktAlt,
 FILE *_FileLSCisAjsCktAltDC
 )
{
  int i=0;
  int idx=0;
  double alpha, ene, ene_v;
  double alpha_p, ene_p, ene_vp;
  double alpha_m, ene_m, ene_vm;
  double *LS_CisAjs_real;
  double *LS_CisAjsCktAlt_real;
  double *LS_CisAjsCktAltDC_real;
  LS_CisAjs_real = (double*)malloc(sizeof(double)*_nCisAjs);
  LS_CisAjsCktAlt_real = (double*)malloc(sizeof(double)*_nCisAjsCktAlt);
  LS_CisAjsCktAltDC_real = (double*)malloc(sizeof(double)*_nCisAjsCktAltDC);

  if (CalculateEne(_QQQQ_real[2], _QQQQ_real[3], _QQQQ_real[10],
                   _QQQQ_real[11], _QQQQ_real[15], sampleCount, &alpha_p, &ene_p,
                   &ene_vp, &alpha_m, &ene_m, &ene_vm) != 0) {
    fprintf(stderr,
            "Error: singular or non-finite first-step Lanczos solution.\n");
    free(LS_CisAjs_real);
    free(LS_CisAjsCktAlt_real);
    free(LS_CisAjsCktAltDC_real);
    return -1;
  }

  //determine alpha
  if (ene_p > ene_m) {
    alpha = alpha_m;
    ene = ene_m;
    ene_v = ene_vm;
  }else {
    alpha = alpha_p;
    ene = ene_p;
    ene_v = ene_vp;
  }

  fprintf(_FileLS, "% .18e  ", ene);
  fprintf(_FileLS, "% .18e  ", ene_v);
  fprintf(_FileLS, "% .18e  ", alpha);

  /* zvo_ls_qqqq.dat */
  for (i = 0; i < _nLSHam * _nLSHam * _nLSHam * _nLSHam; i++) {
    fprintf(_FileLSQQQQ, "% .18e  ", _QQQQ_real[i]);
  }
  fprintf(_FileLSQQQQ, "\n");

  if (_NLanczosmode > 1) {
    //determine alpha
    if (ene_p > ene_m) {
      alpha = alpha_m;
    }else {
      alpha = alpha_p;
    }
    /* zvo_ls_qcisajsq.dat */
#ifdef _DEBUG
    for (i = 0; i < _nLSHam * _nLSHam * _nCisAjs; i++) {
      fprintf(_FileLSQCisAjsQ, "% .18e  ", _QCisAjsQ_real[i]);
      //fprintf(_FileLSQCisAjsQ, "%d % .18e  \n", i, _QCisAjsQ_real[i]);
    }
    fprintf(_FileLSQCisAjsQ, "\n");
#endif

    CalculatePhysVal_real(_QQQQ_real[2], _QQQQ_real[3],
              alpha, _QCisAjsQ_real, _nCisAjs,
              _nLSHam, LS_CisAjs_real);
    /* zvo_ls_cisajs.dat */
    for (i = 0; i < _nCisAjsLz; i++) {
      idx=_iOneBodyGIdx[_CisAjsLzIdx[i][0]+_CisAjsLzIdx[i][1]*_Ns][_CisAjsLzIdx[i][2]+_CisAjsLzIdx[i][3]*_Ns];
      fprintf(_FileLSCisAjs, "%d %d %d %d % .18e 0.0 \n", _CisAjsLzIdx[idx][0], _CisAjsLzIdx[idx][1], _CisAjsLzIdx[idx][2], _CisAjsLzIdx[idx][3], LS_CisAjs_real[idx]);
    }
    fprintf(_FileLSCisAjs, "\n");

    /* zvo_ls_qcisajscktaltq.dat */
#ifdef _DEBUG
    for (i = 0; i < _nLSHam * _nLSHam * _nCisAjsCktAlt; i++) {
      fprintf(_FileLSQCisAjsCktAltQ, "% .18e  ", _QCisAjsCktAltQ_real[i]);
    }
    fprintf(_FileLSQCisAjsCktAltQ, "\n");
#endif

    CalculatePhysVal_real(_QQQQ_real[2], _QQQQ_real[3],
              alpha, _QCisAjsCktAltQDC_real, _nCisAjsCktAltDC,
              _nLSHam, LS_CisAjsCktAltDC_real);
    /* zvo_ls_cisajscktalt.dat */
    for (i = 0; i < _nCisAjsCktAltDC; i++) {
      fprintf(_FileLSCisAjsCktAltDC, "%d %d %d %d %d %d %d %d % .18e 0.0\n",
              _CisAjsCktAltDC[i][0], _CisAjsCktAltDC[i][1], _CisAjsCktAltDC[i][2], _CisAjsCktAltDC[i][3],
              _CisAjsCktAltDC[i][4], _CisAjsCktAltDC[i][5], _CisAjsCktAltDC[i][6], _CisAjsCktAltDC[i][7],
              LS_CisAjsCktAltDC_real[i]);
    }
    fprintf(_FileLSCisAjsCktAltDC, "\n");

    CalculatePhysVal_real(_QQQQ_real[2], _QQQQ_real[3],
              alpha, _QCisAjsCktAltQ_real, _nCisAjsCktAlt,
              _nLSHam, LS_CisAjsCktAlt_real);
    /* zvo_ls_cisajscktaltex.dat */
    for (i = 0; i < _nCisAjsCktAlt; i++) {
      fprintf(_FileLSCisAjsCktAlt, "% .18e 0.0 ", LS_CisAjsCktAlt_real[i]); // LS_CisAjsCktAlt_real
    }
    fprintf(_FileLSCisAjsCktAlt, "\n");
  }
  free(LS_CisAjs_real);
  free(LS_CisAjsCktAlt_real);
  free(LS_CisAjsCktAltDC_real);
  return 0;
}

int PhysCalLanczos_fcmp(
 double complex* _QQQQ,
 double complex* _QCisAjsQ,
 double complex* _QCisAjsCktAltQ,
 double complex* _QCisAjsCktAltQDC,
 const int _nLSHam,
 const int _Ns,
 const int _nCisAjs,
 const int _nCisAjsLz,
 int **_iOneBodyGIdx,
 int **_CisAjsLzIdx,
 const int _nCisAjsCktAlt,
 const int _nCisAjsCktAltDC,
 int **_CisAjsCktAltDC,
 const int _NLanczosmode,
 size_t sampleCount,
 FILE *_FileLS,
 FILE *_FileLSQQQQ,
 FILE *_FileLSQCisAjsQ,
 FILE *_FileLSQCisAjsCktAltQ,
 FILE *_FileLSCisAjs,
 FILE *_FileLSCisAjsCktAlt,
 FILE *_FileLSCisAjsCktAltDC
)
{
  int i=0;
  int idx=0;
  double alpha, ene, ene_v;
  double alpha_p, ene_p, ene_vp;
  double alpha_m, ene_m, ene_vm;
  double complex*LS_CisAjs;
  double complex*LS_CisAjsCktAlt;
  double complex*LS_CisAjsCktAltDC;
  LS_CisAjs = (double complex*)malloc(sizeof(double complex)*_nCisAjs);
  LS_CisAjsCktAlt = (double complex*)malloc(sizeof(double complex)*_nCisAjsCktAlt);
  LS_CisAjsCktAltDC = (double complex*)malloc(sizeof(double complex)*_nCisAjsCktAltDC);

  /* zvo_ls.dat */
  if (CalculateEne(creal(_QQQQ[2]), creal(_QQQQ[3]),
                   creal(_QQQQ[10]), creal(_QQQQ[11]),
                   creal(_QQQQ[15]), sampleCount, &alpha_p, &ene_p, &ene_vp, &alpha_m,
                   &ene_m, &ene_vm) != 0) {
    fprintf(stderr,"Error: singular or non-finite first-step Lanczos solution.\n");
    free(LS_CisAjs);
    free(LS_CisAjsCktAlt);
    free(LS_CisAjsCktAltDC);
    return -1;
  }
  //determine alpha
  if (ene_p > ene_m) {
    alpha = alpha_m;
    ene = ene_m;
    ene_v = ene_vm;
  }else {
    alpha = alpha_p;
    ene = ene_p;
    ene_v = ene_vp;
  }

  fprintf(_FileLS, "% .18e  ", ene);
  fprintf(_FileLS, "% .18e  ", ene_v);
  fprintf(_FileLS, "% .18e  ", alpha);

  /* zvo_ls_qqqq.dat */
  for (i = 0; i < _nLSHam * _nLSHam * _nLSHam * _nLSHam; i++) {
    fprintf(_FileLSQQQQ, "% .18e  ", creal(_QQQQ[i]));
  }
  fprintf(_FileLSQQQQ, "\n");

  if (_NLanczosmode > 1) {

#ifdef _DEBUG
    /* zvo_ls_qcisajsq.dat */
    for (i = 0; i < _nLSHam * _nLSHam * _nCisAjs; i++) {
      fprintf(_FileLSQCisAjsQ, "% .18e  ", creal(_QCisAjsQ[i]));
    }
    fprintf(_FileLSQCisAjsQ, "\n");
#endif

    CalculatePhysVal_fcmp(_QQQQ[2], _QQQQ[3],
              alpha, _QCisAjsQ, _nCisAjs,
              _nLSHam, LS_CisAjs);
    /* zvo_ls_cisajs.dat */

    for (i = 0; i < _nCisAjsLz; i++) {
      idx=_iOneBodyGIdx[_CisAjsLzIdx[i][0]+_CisAjsLzIdx[i][1]*_Ns][_CisAjsLzIdx[i][2]+_CisAjsLzIdx[i][3]*_Ns];
      fprintf(_FileLSCisAjs, "%d %d %d %d % .18e % .18e \n", _CisAjsLzIdx[idx][0], _CisAjsLzIdx[idx][1], _CisAjsLzIdx[idx][2], _CisAjsLzIdx[idx][3], creal(LS_CisAjs[idx]), cimag(LS_CisAjs[idx]));
    }
    fprintf(_FileLSCisAjs, "\n");

    /* zvo_ls_qcisajscktaltq.dat */
#ifdef _DEBUG
    for (i = 0; i < _nLSHam * _nLSHam * _nCisAjsCktAlt; i++) {
      fprintf(_FileLSQCisAjsCktAltQ, "% .18e ",  creal(_QCisAjsCktAltQ[i]));
    }
    fprintf(_FileLSQCisAjsCktAltQ, "\n");
#endif
    CalculatePhysVal_fcmp(_QQQQ[2], _QQQQ[3],
              alpha, _QCisAjsCktAltQDC, _nCisAjsCktAltDC,
              _nLSHam, LS_CisAjsCktAltDC);
    /* zvo_ls_cisajscktalt.dat */
    for (i = 0; i < _nCisAjsCktAltDC; i++) {
      fprintf(_FileLSCisAjsCktAltDC, "%d %d %d %d %d %d %d %d % .18e % .18e\n",
              _CisAjsCktAltDC[i][0], _CisAjsCktAltDC[i][1], _CisAjsCktAltDC[i][2], _CisAjsCktAltDC[i][3],
              _CisAjsCktAltDC[i][4], _CisAjsCktAltDC[i][5], _CisAjsCktAltDC[i][6], _CisAjsCktAltDC[i][7],
              creal(LS_CisAjsCktAltDC[i]), cimag(LS_CisAjsCktAltDC[i]));
    }
    fprintf(_FileLSCisAjsCktAltDC, "\n");

    CalculatePhysVal_fcmp(_QQQQ[2], _QQQQ[3],
              alpha, _QCisAjsCktAltQ, _nCisAjsCktAlt,
              _nLSHam, LS_CisAjsCktAlt);
    /* zvo_ls_cisajscktaltex.dat */
    for (i = 0; i < _nCisAjsCktAlt; i++) {
      fprintf(_FileLSCisAjsCktAlt, "% .18e % .18e ", creal(LS_CisAjsCktAlt[i]), cimag(LS_CisAjsCktAlt[i])); // LS_CisAjsCktAlt_real
    }
    fprintf(_FileLSCisAjsCktAlt, "\n");
  }
  free(LS_CisAjs);
  free(LS_CisAjsCktAlt);
  free(LS_CisAjsCktAltDC);
  return 0;
}

int CalculateEne(double H1, double H2_1, double H2_2, double H3, double H4,
                 size_t sampleCount, double *alpha_p, double *ene_p,
                 double *ene_vp, double *alpha_m, double *ene_m, double *ene_vm) {
  double tmp_AA, tmp_BB, tmp_CC, tmp_xp, tmp_xm;
  /* det S for (psi, H psi) is the sampled local-energy variance.
     Summation error in each moment grows with the number of samples.
     Use the valid rows from all chains in this bin, including split chains. */
  const double relativeTolerance = fmax(64.0, 4.0 * (double)sampleCount) * DBL_EPSILON;
  const double h1Squared = H1 * H1;
  if (!isfinite(H1) || !isfinite(H2_1) || !isfinite(H2_2) ||
      !isfinite(H3) || !isfinite(H4) ||
      H2_1 - h1Squared <= relativeTolerance * fmax(H2_1, h1Squared)) {
    return -1;
  }

  tmp_AA  = H2_1*(H2_1+H2_2)-2*H1*H3;
  tmp_BB  = -H1*H2_1+H3;
  tmp_CC  = H2_1*pow(H2_1+H2_2,2)-pow(H1,2)*H2_1*(H2_1+2.0*H2_2)+4*pow(H1, 3)*H3-2.0*H1*(2*H2_1+H2_2)*H3+H3*H3;
  if (!isfinite(tmp_AA) || tmp_AA == 0.0 || !isfinite(tmp_CC) || tmp_CC < 0.0) {
    return -1;
  }
  tmp_xp  = (tmp_BB+sqrt(tmp_CC))/tmp_AA;
  tmp_xm  = (tmp_BB-sqrt(tmp_CC))/tmp_AA;
  *alpha_p = tmp_xp ;
  *alpha_m = tmp_xm ;

  //calculate energy
  if (CalculateEneByAlpha(H1, H2_1, H2_2, H3, H4, *alpha_p, ene_p, ene_vp) != 0) {
    return -1;
  }
  if (CalculateEneByAlpha(H1, H2_1, H2_2, H3, H4, *alpha_m, ene_m, ene_vm) != 0) {
    return -1;
  }

  /* Only the lower-energy candidate is written. Its relative variance
     must also be defined when the other candidate has zero energy. */
  return isfinite(*ene_p > *ene_m ? *ene_vm : *ene_vp) ? 0 : -1;
}

int CalculateEneByAlpha(
    double H1, double H2_1, double H2_2,
    double H3, double H4, double alpha,
    double *ene, double*ene_V
){
  double tmp_ene, dnorm, tmp_ene_V;
  tmp_ene        = H1+alpha*(H2_1+H2_2)+alpha*alpha*H3;
  dnorm        = 1.0+2*alpha*H1+alpha*alpha*H2_1;
  tmp_ene_V        = H2_1+2*alpha*H3+alpha*alpha*H4;
  if (!isfinite(alpha) || !isfinite(dnorm) || dnorm <= 0.0 ||
      dnorm < 1e-12 * fabs(H1)) return -1;
  *ene = tmp_ene / dnorm;
  if (!isfinite(*ene)) return -1;
  const double e2 = (*ene) * (*ene);
  *ene_V = e2 > 0.0 && isfinite(e2) ? (tmp_ene_V / dnorm - e2) / e2 : NAN;
  return 0;
}

int CalculatePhysVal_real(
    double H1, double H2_1, double alpha,
    double *_QPhysQ_real, int _NPhys,  int _NLSHam,
    double *_Phys_LS_real
)
{
  int i;
  double A0;
  double A1_01;
  double A1_10;
  double A2_11;
  double dnorm = 1.0+2*alpha*H1+alpha*alpha*H2_1;
  for (i = 0; i < _NPhys; i++) {
    A0                = _QPhysQ_real[i];//00
    A1_01             = _QPhysQ_real[_NPhys+i];//01
    A1_10             = _QPhysQ_real[_NLSHam*_NPhys+i];//10
    A2_11             = _QPhysQ_real[_NLSHam*_NPhys+_NPhys+i];//11
    _Phys_LS_real[i]  = A0+alpha*(A1_01+A1_10)+alpha*alpha*A2_11;
    _Phys_LS_real[i] /= dnorm;
  }
  return 0;
}

int CalculatePhysVal_fcmp(
    double complex H1, double complex H2_1, double complex alpha,
    double complex *_QPhysQ, int _NPhys,  int _NLSHam,
    double complex *_Phys_LS
)
{
  int i;
  double complex A0;
  double complex A1_01;
  double complex A1_10;
  double complex A2_11;
  double dnorm = creal(1.0+2*alpha*H1+alpha*alpha*H2_1);
  for (i = 0; i < _NPhys; i++) {
    A0           = _QPhysQ[i];
    A1_01        = _QPhysQ[_NPhys+i];//01
    A1_10        = _QPhysQ[_NLSHam*_NPhys+i];//10
    A2_11        = _QPhysQ[_NLSHam*_NPhys+_NPhys+i];//11
    _Phys_LS[i]  = A0+alpha*(A1_01+A1_10)+alpha*alpha*A2_11;
    _Phys_LS[i] /= dnorm;
    //fprintf(stdout, "Debug: i=%d, A0=%lf, A1_01=%lf, A1_10=%lf, A2_11=%lf, _Phys_LS=%lf.\n", i, creal(A0), creal(A1_01), creal(A1_10), creal(A2_11), creal(_Phys_LS[i]));
  }
  return 0;
}
