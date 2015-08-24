/*-------------------------------------------------------------
 *[ver.2009.05.25]
 *-------------------------------------------------------------
 * Copyright (C) 2007-2009 Daisuke Tahara. All rights reserved.
 * Copyright (C) 2009-     Takahiro Misawa. All rights reserved.
 * Some functions are added by TM.
 *-------------------------------------------------------------*/


/*=================================================================================================*/
typedef struct{
  double real;
  double img;
}doublecomplex;

int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
int dgemm_(char *jobz, char *uplo, int *m,int *n,int *k,double *alpha,  double *a, int *lda, double *b, int *ldb, double *beta,double *c,int *ldc);
int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
int zheev_(char *jobz, char *uplo, int *n, doublecomplex *a, int *lda, double *w, doublecomplex *work, int *lwork, double *rwork, int *info);

double dlamch_(char *cmach);
int dsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, 
	int *il, int *iu, double *abstol, int *m, double *w, double *z__, int *ldz, 
	double *work, int *lwork, int *iwork, int *ifail, int *info);

int MMProd(int Ns, int Ne, double **Mat_1, double **Mat_2,double **Mat_3){

	int i,j,k;
	int m,n,lda,ldb,ldc,info,lwork;
	char jobz, uplo;
	double *a,*b,*c;
  double alpha,beta;

  alpha = 1.0;
  beta  = 0.0;

	m=n=lda=lwork=Ns;

	a = (double*)malloc(Ns*Ne*sizeof(double));
	b = (double*)malloc(Ns*Ne*sizeof(double));
	c = (double*)malloc(Ns*Ns*sizeof(double));

	k=0;
	for(j=0;j<Ne;j++){
		for(i=0;i<Ns;i++){
			a[k] = Mat_1[i][j];
			k++;
		}
	}

	k=0;
	for(j=0;j<Ns;j++){
		for(i=0;i<Ne;i++){
			b[k] = Mat_2[i][j];
			k++;
		}
	}

	k=0;
	for(j=0;j<Ns;j++){
		for(i=0;i<Ns;i++){
			c[k] = Mat_3[i][j];
			k++;
		}
	}

	jobz = 'N';
	uplo = 'N';
	dgemm_(&jobz,&uplo,&Ns,&Ns,&Ne,&alpha,a,&Ns,b,&Ne,&beta,c,&Ns);

	k=0;
  for(j=0;j<Ns;j++){
  	for(i=0;i<Ns;i++){
			Mat_3[i][j] = c[k];
			k++;
		}
	}

	free(a);
	free(b);
	free(c);

	return 1;
}




int DInv(int xNsize, double **xM, double **xIM){

	int i,j,k;
	int m,n,lda,info,*piv,lwork;
	double *work;
	double *a;

	m=n=lda=lwork=xNsize;

	a = (double*)malloc(xNsize*xNsize*sizeof(double));
	work = (double*)malloc(xNsize*sizeof(double));
	piv = (int*)malloc(xNsize*sizeof(int));

	k=0;
	for(j=0;j<xNsize;j++){
		for(i=0;i<xNsize;i++){
			a[k] = xM[i][j];
			k++;
		}
	}

	dgetrf_(&m, &n, a, &lda, piv, &info);
	dgetri_(&n, a, &lda, piv, work, &lwork, &info);

	if(info != 0){
		free(a);
		free(work);
		free(piv);
		return 0;
	}

	for(k=0;k<xNsize*xNsize;k++){
		xIM[k%xNsize][k/xNsize] = a[k];
	}
	free(a);
	free(work);
	free(piv);

	return 1;
}



//added by Misawa 090623
//For complex Hermite matrix
int ZHEEVvalue(int xNsize, doublecomplex **A, double *r){
	int i,j,k;
	char jobz, uplo;
	int n, lda, lwork, info;
  double *rwork;
	double *w;
	doublecomplex *a, *work;

	n = lda = xNsize;
	lwork = 4*xNsize; /* 3*xNsize OK?*/

	a = (doublecomplex*)malloc(xNsize*xNsize*sizeof(doublecomplex));
	w = (double*)malloc(xNsize*sizeof(double));
	work = (doublecomplex*)malloc(lwork*sizeof(doublecomplex));
	rwork = (double*)malloc(lwork*sizeof(double));

	k=0;
	for(j=0;j<xNsize;j++){
		for(i=0;i<xNsize;i++){
			a[k].real = A[i][j].real;
			a[k].img = A[i][j].img;
			k++;
		}
	}

	jobz = 'N';
	uplo = 'U';

	zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);

	if(info != 0){
		free(a);
		free(w);
		free(work);
		return 0;
	}

  
	for(k=0;k<xNsize;k++){
		r[k] = w[k];
	}

	free(a);
	free(w);
	free(work);

	return 1;
}

int DSEVvalue(int xNsize, double **A, double *r){
	int i,j,k;
	char jobz, uplo;
	int n, lda, lwork, info;
	double *a, *w, *work;

	n = lda = xNsize;
	lwork = 4*xNsize; /* 3*xNsize OK?*/

	a = (double*)malloc(xNsize*xNsize*sizeof(double));
	w = (double*)malloc(xNsize*sizeof(double));
	work = (double*)malloc(lwork*sizeof(double));

	k=0;
	for(j=0;j<xNsize;j++){
		for(i=0;i<xNsize;i++){
			a[k] = A[i][j];
			k++;
		}
	}

	jobz = 'N';
	uplo = 'U';

	dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

	if(info != 0){
		free(a);
		free(w);
		free(work);
		return 0;
	}

	for(k=0;k<xNsize;k++){
		r[k] = w[k];
	}

	free(a);
	free(w);
	free(work);

	return 1;
}
// added by Misawa 20090309 
// obtain eigen vectors 
int DSEVvector(int xNsize, double **A, double *r, double **vec ){

	int i,j,k;
	char jobz, uplo;
	int n, lda, lwork, info;
	double *a, *w, *work;

	n = lda = xNsize;
	lwork = 4*xNsize; /* 3*xNsize OK?*/

	a = (double*)malloc(xNsize*xNsize*sizeof(double));
	w = (double*)malloc(xNsize*sizeof(double));
	work = (double*)malloc(lwork*sizeof(double));

	k=0;
	for(j=0;j<xNsize;j++){
		for(i=0;i<xNsize;i++){
			a[k] = A[i][j];
			k++;
		}
	}

	jobz = 'V';
	uplo = 'U';

	dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);


	if(info != 0){
		free(a);
		free(w);
		free(work);
		return 0;
	}

	k=0;
	for(i=0;i<xNsize;i++){
		for(j=0;j<xNsize;j++){
			vec[i][j]=a[k];
			k++;
		}
	}

	for(k=0;k<xNsize;k++){
		r[k] = w[k];
	}

	free(a);
	free(w);
	free(work);

	return 1;
}

int DSEVXU(int xNsize, double **A, double *r, double **X, int xNev){

	int i,j,k;
	char jobz, range, uplo;
	int n, lda, il, iu, m, ldz, lwork, *iwork, *ifail, info;
	double *a, *w, *work, vl, vu, abstol, *z;

	n = lda = ldz = xNsize;
	lwork = 8*xNsize;

	a = (double*)malloc(xNsize*xNsize*sizeof(double));
	w = (double*)malloc(xNsize*sizeof(double));
	z = (double*)malloc(xNsize*xNsize*sizeof(double));
	work  = (double*)malloc(lwork*sizeof(double));
	iwork = (int*)malloc(5*xNsize*sizeof(int));
	ifail = (int*)malloc(xNsize*sizeof(int));

	abstol = 2.0*dlamch_("S");
	vl = vu = 0.0;

	k=0;
	for(j=0;j<xNsize;j++){
		for(i=0;i<xNsize;i++){
			a[k] = A[i][j];
			k++;
		}
	}

	jobz  = 'V';
	range = 'I';
	uplo  = 'U';

	il = xNsize-xNev+1;
	iu = xNsize;
	m = iu-il+1;

	dsyevx_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, 
			&m, w, z, &ldz, work, &lwork, iwork, ifail, &info);

	if(info != 0){
		free(a);
		free(w);
		free(z);
		free(work);
		free(iwork);
		free(ifail);
		return 0;
	}

	for(k=0;k<xNev;k++){
		r[k+xNsize-xNev] = w[k];
	}

	for(k=0;k<xNsize*xNev;k++){
		X[k%xNsize][k/xNsize+xNsize-xNev] = z[k];
	}
	free(a);
	free(w);
	free(z);
	free(work);
	free(iwork);
	free(ifail);

	return 1;
}
