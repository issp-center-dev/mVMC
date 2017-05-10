/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

his program is developed based on the mVMC-mini program
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

// #define _DEBUG_STCOPT_CG
#define _DEBUG_STCOPT_CG_GEMV
// #define _DEBUG_STCOPT_CG_LAPACK


#ifdef _SYSTEM_A
 #define M_ZPOSV  ZPOSV
 #define M_ZGETRS  ZGETRS
 #define M_ZGETRF  ZGETRF
 #define M_DGEMV  DGEMV
#ifdef _DEBUG_STCOPT_CG_GEMV
  #define M_DGEMM  DGEMM
  #define M_DGER  DGER
#endif
#elif _lapack_small_nounderscore
 #define M_ZPOSV  zposv
 #define M_ZGETRS  zgetrs
 #define M_ZGETRF  zgetrf
 #define M_DGEMV  dgemv
#ifdef _DEBUG_STCOPT_CG_GEMV
 #define M_DGEMM  dgemm
 #define M_DGER  dger
#endif
#else
 #define M_ZPOSV  zposv_
 #define M_ZGETRS  zgetrs_
 #define M_ZGETRF  zgetrf_
 #define M_DGEMV  dgemv_
 #define M_ZGEMV  zgemv_
#ifdef _DEBUG_STCOPT_CG_GEMV
// #define M_DGEMM  dgemm_
 #define M_ZGEMM  dzemm_
 #define M_DGER   dger_
// #define M_ZGERC  zgerc_
#endif
#endif

extern void M_DGEMV(const char *trans, const int *m, const int *n, const double *alpha,
                    const double *a, const int *lda, const double *x, const int *incx,
                    const double *beta, double *y, const int *incy);


extern void M_ZGEMV(const char *trans, const int *m, const int *n, const double complex *alpha,
                    const double complex *a, const int *lda, const double complex *x, const int *incx,
                    const double complex *beta, double complex *y, const int *incy);

#ifdef _DEBUG_STCOPT_CG_GEMV
/*
extern void M_DGEMM(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                    const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
                    const double *beta, double *c, const int *ldc);
                    */
extern void M_ZGEMM(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                    const double complex *alpha, const double complex *a, const int *lda,
                    const double complex *b, const int *ldb, const double complex *beta,
                    double complex *c, const int *ldc);
extern void M_DGER(const int *m, const int *n, const double *alpha, const double *x, const int *incx,
                   const double *y, const int *incy, double *a, const int *lda);
/*
extern void M_ZGERC(const int *m, const int *n, const double complex *alpha, const double complex *x, const int *incx,
                    const double complex *y, const int *incy, double complex *a, const int *lda);
                    */
#endif

#define MVMC_SRCG_REAL
#include "stcopt_cg_impl.c"
#undef MVMC_SRCG_REAL
#include "stcopt_cg_impl.c"

