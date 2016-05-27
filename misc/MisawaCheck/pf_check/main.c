/*-------------------------------------------------------------
 * Variational Monte Carlo
 * matrix Package (LAPACK and Pfapack)
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>


  #define M_DGETRF dgetrf_
  #define M_DGETRI dgetri_
  #define M_DSKPFA dskpfa_
  #define M_ZGETRF zgetrf_
  #define M_ZGETRI zgetri_
  #define M_ZSKPFA zskpfa_


int M_DGETRF(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
int M_DGETRI(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
int M_DSKPFA(const char *uplo, const char *mthd, const int *n,
             double *a, const int *lda, double *pfaff, int *iwork,
             double *work, const int *lwork, int *info);
int M_ZGETRF(int *m, int *n, double complex *a, int *lda, int *ipiv, int *info);
int M_ZGETRI(int *n, double complex *a, int *lda, int *ipiv, double complex *work, int *lwork, int *info);
int M_ZSKPFA(const char *uplo, const char *mthd, const int *n,
             double complex *a, const int *lda, double complex *pfaff, int *iwork,
             double complex *work, const int *lwork, double *rwork, int *info);
int M_ZSKPFA(const char *uplo, const char *mthd, const int *n,
             double complex *a, const int *lda, double complex *pfaff, int *iwork,
             double complex *work, const int *lwork, double *rwork, int *info);


int main(void){
  char uplo='U', mthd='H';
  int Nsize,n,lda;
  double complex *A;
  double complex pfaff;
  int *iwork;
  double complex  *work;
  int lwork;
  double *rwork;
  int info;
  int i,j;
  
  Nsize = 4;
  n = lda = Nsize;
  lwork = n;

  A     = (double complex*)malloc(sizeof(double complex )*Nsize*Nsize);
  iwork = (int *)malloc(sizeof(int)*lwork);
  work  = (double complex*)malloc(sizeof(double complex)*lwork);
  rwork = (double *)malloc(sizeof(double )*(Nsize-1));

  for(i=0;i<Nsize;i++){
    for(j=i;j<Nsize;j++){
      if(i==j){
         A[j+Nsize*i] = 0.0;
      }else{
         A[j+Nsize*i] = i+2*j+0.2*i*I;
         A[i+Nsize*j] = -(i+2*j)-0.2*i*I;
      }
    }
  }
  for(i=0;i<Nsize;i++){
    for(j=0;j<Nsize;j++){
      printf(" %lf + %lf j",creal(A[j+i*Nsize]),cimag(A[j+i*Nsize]));
    }  
    printf(" \n");
  }

  //M_DSKPFA(&uplo, &mthd, &n, &a, &lda, &pfaff, &iwork, &optSize2, &lwork, &info);
  M_ZSKPFA(&uplo, &mthd, &n, A, &lda, &pfaff, iwork, work, &lwork, rwork, &info);

  printf(" %lf %lf: infp = %d \n",creal(pfaff),cimag(pfaff),info);

  return 0;

}

