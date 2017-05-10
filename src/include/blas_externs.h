/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/

#ifndef _BLAS_EXTERNS_H
#define _BLAS_EXTERNS_H

// BLAS
#define M_DAXPY daxpy_
#define M_DGEMM  dgemm_
#define M_DGEMV  dgemv_
#define M_DGER  dger_
#define M_ZAXPY zaxpy_
#define M_ZGEMM zgemm_
#define M_ZGERC zgerc_

// LAPACK
#define M_DGETRF dgetrf_
#define M_DGETRI dgetri_
#define M_DPOSV  dposv_
#define M_ZGETRF zgetrf_
#define M_ZGETRI zgetri_
#define M_ZPOSV  zposv_

// pfaPACK
#define M_DSKPFA dskpfa_
#define M_ZSKPFA zskpfa_

// pBLAS
#define M_PDGEMV  pdgemv_

// ScaLAPACK
#define M_PDPOSV  pdposv_
#define M_PDSYEVD  pdsyevd_

// blacs
#define M_NUMROC   numroc_
#define M_DESCINIT descinit_

// BLAS

extern void M_DGEMV(const char *trans, const int *m, const int *n, const double *alpha,
                    const double *a, const int *lda, const double *x, const int *incx,
                    const double *beta, double *y, const int *incy);
extern void M_DGEMM(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                    const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
                    const double *beta, double *c, const int *ldc);
extern void M_ZGEMM(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                    const double complex *alpha, const double complex *a, const int *lda,
                    const double complex *b, const int *ldb, const double complex *beta,
                    double complex *c, const int *ldc);
extern void M_DGER(const int *m, const int *n, const double *alpha, const double *x, const int *incx,
                   const double *y, const int *incy, double *a, const int *lda);
extern void M_ZGERC(const int *m, const int *n, const double complex *alpha, const double complex *x, const int *incx,
                    const double complex *y, const int *incy, double complex *a, const int *lda);
extern void M_DAXPY(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern void M_ZAXPY(const int *n, const double complex *alpha, const double complex *x, const int *incx, double complex *y, const int *incy);

// LAPACK
extern void M_DPOSV(const char* uplo, const int* n, const int* nrhs, double* a,
                    const int* lda, double* b, const int* ldb, int* info );
extern void M_DGETRF(const int* m, const int* n, double* a, const int* lda,
                     int* ipiv, int* info );
extern void M_DGETRI(const int* n, double* a, const int* lda,
                     const int* ipiv, double* work, const int* lwork,
                     int* info );
extern void M_ZGETRF(const int* m, const int* n, double complex* a,
                     const int* lda, int* ipiv, int* info );
extern void M_ZGETRI(const int* n, double complex* a, const int* lda,
                     const int* ipiv, double complex* work, const int* lwork,
                     int* info );

// pfapack
extern int M_DSKPFA(const char *uplo, const char *mthd, const int *n,
                    double *a, const int *lda, double *pfaff, int *iwork,
                    double *work, const int *lwork, int *info);
extern int M_ZSKPFA(const char *uplo, const char *mthd, const int *n,
                    double complex *a, const int *lda, double complex *pfaff, int *iwork,
                    double complex *work, const int *lwork, double *rwork, int *info);

// pBLAS
extern void M_PDGEMV(const char *trans, const int *m, const int *n,
                     const double *alpha, const double *a, const int *ia, const int *ja, const int *desca,
                     const double *x, const int *ix, const int *jx, const int *descx, const int *incx,
                     const double *beta, double *y, const int *iy, const int *jy, const int *descy, const int *incy );

// ScaLAPACK
extern void	M_PDSYEVD(const char* jobz, const char* uplo,
                      const int* n, const double* a, const int* ia, const int* ja, const int* desca,
                      double* w, double* z, const int* iz, const int* jz, const int* descz,
                      double* work, const int* lwork, int* iwork, const int* liwork, int* info);
extern void	M_PDPOSV(const char* uplo, const int* n, const int* nrhs,
                     double* a, const int* ia, const int* ja, const int* desca,
                     double* b, const int* ib, const int* jb, const int* descb, int* info);


// blacs
extern int Csys2blacs_handle(MPI_Comm comm);
extern int Cblacs_gridinit(int *ictxt, char *order, int nprow, int npcol);
extern int Cblacs_gridexit(int ictxt);
extern int Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *myprow, int *mypcol);
extern int M_NUMROC(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int M_DESCINIT(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc,
                      int *icsrc, int *ictxt, int *lld, int *info);

#endif // _BLAS_EXTERNS_H
