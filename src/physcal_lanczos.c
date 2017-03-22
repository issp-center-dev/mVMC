#include "physcal_lanczos.h"
#include "math.h"
int CalculateEne(
		double H1, double H2_1, double H2_2, double H3, double H4,
		double *alpha_p, double *ene_p, double *ene_vp,
		double *alpha_m, double *ene_m, double *ene_vm
);
int CalculateEneByAlpha(double H1, double H2_1, double H2_2,
						double H3, double H4, double alpha,
						double *ene, double*ene_V);

int CalculateEne_fcmp(
		double complex H1, double complex H2_1, double complex H2_2,
		double complex H3, double complex H4,
		double *alpha_p, double *ene_p, double *ene_vp,
		double *alpha_m, double *ene_m, double *ene_vm
);

int CalculateEneByAlpha_fcmp(
		double complex H1, double complex H2_1, double complex H2_2,
		double complex H3, double complex H4, double alpha,
		double *ene, double*ene_V
);

int PhysCalLanczos_real
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
	double alpha_p, ene_p, ene_vp;
	double alpha_m, ene_m, ene_vm;

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
	fprintf(_FileLS, "% .18e  ", alpha_p);
	fprintf(_FileLS, "% .18e  ", ene_p);
	fprintf(_FileLS, "% .18e  ", ene_vp);
	fprintf(_FileLS, "% .18e  ", alpha_m);
	fprintf(_FileLS, "% .18e  ", ene_m);
	fprintf(_FileLS, "% .18e \n", ene_vm);

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
	return 0;
}

int PhysCalLanczos_fcmp(
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
	double alpha_p, ene_p, ene_vp;
	double alpha_m, ene_m, ene_vm;

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
		return -1;
	}

	/*
	if(!CalculateEne_fcmp(_QQQQ[2],_QQQQ[3],
					 _QQQQ[10], _QQQQ[11], _QQQQ[15],
					 &alpha_p,  &ene_p,  &ene_vp, &alpha_m,  &ene_m,  &ene_vm)==0){
		return -1;
	}
	 */

	fprintf(_FileLS, "% .18e  ", alpha_p);
	fprintf(_FileLS, "% .18e  ", ene_p);
	fprintf(_FileLS, "% .18e  ", ene_vp);
	fprintf(_FileLS, "% .18e  ", alpha_m);
	fprintf(_FileLS, "% .18e  ", ene_m);
	fprintf(_FileLS, "% .18e \n", ene_vm);

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

	return 0;
}

int CalculateEne(double H1, double H2_1, double H2_2, double H3, double H4,
				 double *alpha_p, double *ene_p, double *ene_vp,
				 double *alpha_m, double *ene_m, double *ene_vm
){
	double tmp_AA, tmp_BB, tmp_CC, tmp_xp, tmp_xm;
	double alpha;
	double tmp_1, tmp_2, tmp_3, tmp_V, tmp_ene;
	//determine alpha
	tmp_AA  = H2_1*(H2_1+H2_2)-2*H1*H3;
	tmp_BB  = -H1*H2_1+H3;
	tmp_CC  = H2_1*pow(H2_1+H2_2,2)-pow(H1,2)*H2_1*(H2_1+2.0*H2_2)+4*pow(H1, 3)*H3-2.0*H1*(2*H2_1+H2_2)*H3+H3*H3;
	if(tmp_CC < 0){
		tmp_xp=0;
		tmp_xm=0;
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
	double tmp_1, tmp_2, tmp_3;
	tmp_1        = H1+alpha*(H2_1+H2_2)+alpha*alpha*H3;
	tmp_2        = 1.0+2*alpha*H1+alpha*alpha*H2_1;
	tmp_3        = H2_1+2*alpha*H3+alpha*alpha*H4;
	if(fabs(tmp_2/H1) < pow(10.0, -12)) return -1;
	*ene_V        =  ((tmp_3/tmp_2)-pow((tmp_1/tmp_2), 2))/pow((tmp_1/tmp_2),2);
	*ene		= tmp_1/tmp_2;
	return 0;
}

int CalculateEne_fcmp(
		double complex H1, double complex H2_1, double complex H2_2,
		double complex H3, double complex H4,
		double *alpha_p, double *ene_p, double *ene_vp,
		double *alpha_m, double *ene_m, double *ene_vm
){
	double tmp_AA, tmp_BB, tmp_CC, tmp_xp, tmp_xm;
	//determine alpha
	tmp_AA  = creal(H2_1*(H2_1+H2_2)-2*H1*H3);
	tmp_BB  = creal(H1*H2_1+H3);
	tmp_CC  = creal(H2_1*pow(H2_1+H2_2,2)-pow(H1,2)*H2_1*(H2_1+2.0*H2_2)+4*pow(H1, 3)*H3-2.0*H1*(2*H2_1+H2_2)*H3+H3*H3);
	if(tmp_CC < 0){
		return -1;
	}
	tmp_xp  = (tmp_BB+sqrt(tmp_CC))/tmp_AA;
	tmp_xm  = (tmp_BB-sqrt(tmp_CC))/tmp_AA;
	*alpha_p = creal(tmp_xp);
	*alpha_m = creal(tmp_xm);

	//calculate energy
	if(!CalculateEneByAlpha_fcmp(H1, H2_1, H2_2, H3, H4, *alpha_p, ene_p, ene_vp)==0){
		return -1;
	}
	if(!CalculateEneByAlpha_fcmp(H1, H2_1, H2_2, H3, H4, *alpha_m, ene_m, ene_vm)==0){
		return -1;
	}

	return 0;
}

int CalculateEneByAlpha_fcmp(
		double complex H1, double complex H2_1, double complex H2_2,
		double complex H3, double complex H4, double alpha,
		double *ene, double*ene_V
){
	double complex tmp_1, tmp_2, tmp_3;
	tmp_1        = H1+alpha*(H2_1+H2_2)+alpha*alpha*H3;
	tmp_2        = 1.0+2*alpha*H1+alpha*alpha*H2_1;
	tmp_3        = H2_1+2*alpha*H3+alpha*alpha*H4;
	if(fabs(creal(tmp_2/H1)) < pow(10.0, -12)) return -1;
	*ene_V        =  creal(((tmp_3/tmp_2)-pow((tmp_1/tmp_2), 2))/pow((tmp_1/tmp_2),2));
	*ene		= creal(tmp_1/tmp_2);
	return 0;
}
