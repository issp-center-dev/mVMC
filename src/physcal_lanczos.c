#include "physcal_lanczos.h"
#include "math.h"
#include "stdlib.h"

int CalculateEne(
		double H1, double H2_1, double H2_2, double H3, double H4,
		double *alpha_p, double *ene_p, double *ene_vp,
		double *alpha_m, double *ene_m, double *ene_vm
);
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
 const int _nLSHam,
 const int _Ns,
 const int _nCisAjs,
 const int _nCisAjsLz,
 int **_iOneBodyGIdx,
 int **_CisAjsLzIdx,
 const int _nCisAjsCktAlt,
 int **_CisAjsCktAlt,
 const int _NLanczosmode,
 FILE *_FileLS,
 FILE *_FileLSQQQQ,
 FILE *_FileLSQCisAjsQ,
 FILE *_FileLSQCisAjsCktAltQ,
 FILE *_FileLSCisAjs,
 FILE *_FileLSCisAjsCktAlt
 )
{
  int i=0;
	int idx=0;
	double alpha, ene, ene_v;
	double alpha_p, ene_p, ene_vp;
	double alpha_m, ene_m, ene_vm;
	double *LS_CisAjs_real;
	double *LS_CisAjsCktAlt_real;
	LS_CisAjs_real = (double*)malloc(sizeof(double)*_nCisAjs);
	LS_CisAjsCktAlt_real = (double*)malloc(sizeof(double)*_nCisAjsCktAlt);

	/* zvo_ls.dat */
	/*
  fprintf(_FileLS, "% .18e  ", _QQQQ_real[2]);  // H * I = QQQQ[1],[2],[4],[8]       //TBC
  fprintf(_FileLS, "% .18e  ", _QQQQ_real[3]);  // H * H = QQQQ[3],[6],[9],[12]      //TBC
  fprintf(_FileLS, "% .18e  ", _QQQQ_real[10]); // H^2 * I = QQQQ[5],[10]            //TBC
  fprintf(_FileLS, "% .18e  ", _QQQQ_real[11]); // H^2 * H = QQQQ[7],[11],[13],[14]  //TBC
  fprintf(_FileLS, "% .18e\n", _QQQQ_real[15]); // H^2 * H^2 = QQQQ[15]              //TBC
	*/

	CalculateEne(_QQQQ_real[2], _QQQQ_real[3], _QQQQ_real[10], _QQQQ_real[11], _QQQQ_real[15],
	 &alpha_p,  &ene_p,  &ene_vp, &alpha_m,  &ene_m,  &ene_vm);

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

/*
	fprintf(_FileLS, "% .18e  ", alpha_p);
	fprintf(_FileLS, "% .18e  ", ene_p);
	fprintf(_FileLS, "% .18e  ", ene_vp);
	fprintf(_FileLS, "% .18e  ", alpha_m);
	fprintf(_FileLS, "% .18e  ", ene_m);
	fprintf(_FileLS, "% .18e \n", ene_vm);
*/
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
      //fprintf(_FileLSCisAjs, "% .18e 0.0 ", LS_CisAjs_real[idx]);
		  //fprintf(_FileLSCisAjs, "%d % .18e  \n", i, LS_CisAjs_real[i]);
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
							alpha, _QCisAjsCktAltQ_real, _nCisAjsCktAlt,
							_nLSHam, LS_CisAjsCktAlt_real);
	  /* zvo_ls_cisajscktalt.dat */
	  for (i = 0; i < _nCisAjsCktAlt; i++) {
		  //fprintf(_FileLSCisAjsCktAlt, "% .18e 0.0 ", LS_CisAjsCktAlt_real[i]);
      fprintf(_FileLSCisAjsCktAlt, "%d %d %d %d %d %d %d %d % .18e 0.0\n",
              _CisAjsCktAlt[i][0], _CisAjsCktAlt[i][1], _CisAjsCktAlt[i][2], _CisAjsCktAlt[i][3],
              _CisAjsCktAlt[i][4], _CisAjsCktAlt[i][5], _CisAjsCktAlt[i][6], _CisAjsCktAlt[i][7],
              LS_CisAjsCktAlt_real[i]);

    }
	  fprintf(_FileLSCisAjsCktAlt, "\n");
  }
	free(LS_CisAjs_real);
	free(LS_CisAjsCktAlt_real);
	return 0;
}

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
 int **_CisAjsCktAlt,
 const int _NLanczosmode,
 FILE *_FileLS,
 FILE *_FileLSQQQQ,
 FILE *_FileLSQCisAjsQ,
 FILE *_FileLSQCisAjsCktAltQ,
 FILE *_FileLSCisAjs,
 FILE *_FileLSCisAjsCktAlt
)
{
  int i=0;
	int idx=0;
	double alpha, ene, ene_v;
	double alpha_p, ene_p, ene_vp;
	double alpha_m, ene_m, ene_vm;
	double complex*LS_CisAjs;
	double complex*LS_CisAjsCktAlt;
	LS_CisAjs = (double complex*)malloc(sizeof(double complex)*_nCisAjs);
	LS_CisAjsCktAlt = (double complex*)malloc(sizeof(double complex)*_nCisAjsCktAlt);

	/* zvo_ls.dat */
	/*
  fprintf(_FileLS, "% .18e  ", creal(_QQQQ[2]));  // H * I = _QQQQ[1],[2],[4],[8]       //TBC
  fprintf(_FileLS, "% .18e  ", creal(_QQQQ[3]));  // H * H = _QQQQ[3],[6],[9],[12]      //TBC
  fprintf(_FileLS, "% .18e  ", creal(_QQQQ[10])); // H^2 * I = _QQQQ[5],[10]            //TBC
  fprintf(_FileLS, "% .18e  ", creal(_QQQQ[11])); // H^2 * H = _QQQQ[7],[11],[13],[14]  //TBC
  fprintf(_FileLS, "% .18e\n", creal(_QQQQ[15])); // H^2 * H^2 = _QQQQ[15]              //TBC
*/

	if(!CalculateEne(creal(_QQQQ[2]),creal(_QQQQ[3]),
					 creal(_QQQQ[10]), creal(_QQQQ[11]), creal(_QQQQ[15]),
				 &alpha_p,  &ene_p,  &ene_vp, &alpha_m,  &ene_m,  &ene_vm)==0){
    fprintf(stderr,"Error: Lanczos method is failed due to illegal value of alpha.\n");
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

/*
	fprintf(_FileLS, "% .18e  ", alpha_p);
	fprintf(_FileLS, "% .18e  ", ene_p);
	fprintf(_FileLS, "% .18e  ", ene_vp);
	fprintf(_FileLS, "% .18e  ", alpha_m);
	fprintf(_FileLS, "% .18e  ", ene_m);
	fprintf(_FileLS, "% .18e \n", ene_vm);
*/
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
							alpha, _QCisAjsCktAltQ, _nCisAjsCktAlt,
							_nLSHam, LS_CisAjsCktAlt);
	  /* zvo_ls_cisajs.dat */
	  for (i = 0; i < _nCisAjsCktAlt; i++) {
      fprintf(_FileLSCisAjsCktAlt, "%d %d %d %d %d %d %d %d % .18e % .18e\n",
              _CisAjsCktAlt[i][0], _CisAjsCktAlt[i][1], _CisAjsCktAlt[i][2], _CisAjsCktAlt[i][3],
              _CisAjsCktAlt[i][4], _CisAjsCktAlt[i][5], _CisAjsCktAlt[i][6], _CisAjsCktAlt[i][7],
              creal(LS_CisAjsCktAlt[i]), cimag(LS_CisAjsCktAlt[i]));
		  //fprintf(_FileLSCisAjsCktAlt, "% .18e % .18e ", creal(LS_CisAjsCktAlt[i]), cimag(LS_CisAjsCktAlt[i]));
	  }
	  fprintf(_FileLSCisAjsCktAlt, "\n");
  }
	free(LS_CisAjs);
	free(LS_CisAjsCktAlt);
	return 0;
}

int CalculateEne(double H1, double H2_1, double H2_2, double H3, double H4,
				 double *alpha_p, double *ene_p, double *ene_vp,
				 double *alpha_m, double *ene_m, double *ene_vm
){
	double tmp_AA, tmp_BB, tmp_CC, tmp_xp, tmp_xm;
	//determine alpha
	tmp_AA  = H2_1*(H2_1+H2_2)-2*H1*H3;
	tmp_BB  = -H1*H2_1+H3;
	tmp_CC  = H2_1*pow(H2_1+H2_2,2)-pow(H1,2)*H2_1*(H2_1+2.0*H2_2)+4*pow(H1, 3)*H3-2.0*H1*(2*H2_1+H2_2)*H3+H3*H3;
	if(tmp_CC < 0){
		return -1;
	}
	tmp_xp  = (tmp_BB+sqrt(tmp_CC))/tmp_AA;
	tmp_xm  = (tmp_BB-sqrt(tmp_CC))/tmp_AA;
	*alpha_p = tmp_xp ;
	*alpha_m = tmp_xm ;

	//calculate energy
	if(!CalculateEneByAlpha(H1, H2_1, H2_2, H3, H4, *alpha_p, ene_p, ene_vp)==0){
		return -1;
	}
	if(!CalculateEneByAlpha(H1, H2_1, H2_2, H3, H4, *alpha_m, ene_m, ene_vm)==0){
		return -1;
	}

	return 0;
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
	if(fabs(dnorm/H1) < pow(10.0, -12)) return -1;
	*ene_V        =  ((tmp_ene_V/dnorm)-pow((tmp_ene/dnorm), 2))/pow((tmp_ene/dnorm),2);
	*ene		= tmp_ene/dnorm;
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
      A0=_QPhysQ_real[i];//00
		A1_01=_QPhysQ_real[_NPhys+i];//01
		A1_10=_QPhysQ_real[_NLSHam*_NPhys+i];//10
		A2_11=_QPhysQ_real[_NLSHam*_NPhys+_NPhys+i];//11
		_Phys_LS_real[i]=A0+alpha*(A1_01+A1_10)+alpha*alpha*A2_11;
		_Phys_LS_real[i]/=dnorm;
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
		A0=_QPhysQ[i];
		A1_01=_QPhysQ[_NPhys+i];//01
		A1_10=_QPhysQ[_NLSHam*_NPhys+i];//10
		A2_11=_QPhysQ[_NLSHam*_NPhys+_NPhys+i];//11
		_Phys_LS[i]=A0+alpha*(A1_01+A1_10)+alpha*alpha*A2_11;
		_Phys_LS[i]/=dnorm;
    //fprintf(stdout, "Debug: i=%d, A0=%lf, A1_01=%lf, A1_10=%lf, A2_11=%lf, _Phys_LS=%lf.\n", i, creal(A0), creal(A1_01), creal(A1_10), creal(A2_11), creal(_Phys_LS[i]));
	}
	return 0;
}
