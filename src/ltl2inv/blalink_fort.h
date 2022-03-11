/**
 * \file blalink_fort.hh
 * Wrapper for external BLAS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <stdlib.h>
#include "blis.h"


// converter from BLIS enum to BLAS char.
BLIS_INLINE char trans2char(trans_t t)
{
    switch (t) {
    case BLIS_NO_TRANSPOSE:
        return 'N';
    case BLIS_TRANSPOSE:
        return 'T';
    case BLIS_CONJ_NO_TRANSPOSE:
        return 'C';
    default:
        abort();
    }
}
BLIS_INLINE char uplo2char(uplo_t t)
{
    switch (t) {
    case BLIS_UPPER:
        return 'U';
    case BLIS_LOWER:
        return 'L';
    default:
        abort();
    }
}
BLIS_INLINE char side2char(side_t t)
{
    switch (t) {
    case BLIS_LEFT:
        return 'L';
    case BLIS_RIGHT:
        return 'R';
    default:
        abort();
    }
}
BLIS_INLINE char diag2char(diag_t t)
{
    switch (t) {
    case BLIS_UNIT_DIAG:
        return 'U';
    case BLIS_NONUNIT_DIAG:
        return 'N';
    }
}


#ifdef __cplusplus
extern "C"
{
#endif

// BLAS {
#ifdef BLAS_EXTERNAL

// gemm
void sgemm_(char *transa, char *transb, dim_t *m, dim_t *n, dim_t *k, float *alpha,
            float *a, dim_t *lda, float *b, dim_t *ldb, float *beta, float *c, dim_t *ldc);
void dgemm_(char *transa, char *transb, dim_t *m, dim_t *n, dim_t *k, double *alpha,
            double *a, dim_t *lda, double *b, dim_t *ldb, double *beta, double *c, dim_t *ldc);
void cgemm_(char *transa, char *transb, dim_t *m, dim_t *n, dim_t *k, void *alpha,
            void *a, dim_t *lda, void *b, dim_t *ldb, void *beta, void *c, dim_t *ldc);
void zgemm_(char *transa, char *transb, dim_t *m, dim_t *n, dim_t *k, void *alpha,
            void *a, dim_t *lda, void *b, dim_t *ldb, void *beta, void *c, dim_t *ldc);


// ger
void sger_(dim_t *m, dim_t *n, float *alpha,
           float *x, dim_t *incx, float *y, dim_t *incy, float *a, dim_t *lda);
void dger_(dim_t *m, dim_t *n, double *alpha,
           double *x, dim_t *incx, double *y, dim_t *incy, double *a, dim_t *lda);
void cgeru_(dim_t *m, dim_t *n, void *alpha,
            void *x, dim_t *incx, void *y, dim_t *incy, void *a, dim_t *lda);
void zgeru_(dim_t *m, dim_t *n, void *alpha,
            void *x, dim_t *incx, void *y, dim_t *incy, void *a, dim_t *lda);


// gemv
void sgemv_(char *trans, dim_t *m, dim_t *n, float *alpha, float *a, dim_t *lda,
            float *x, dim_t *incx, float *beta, float *y, dim_t *incy);
void dgemv_(char *trans, dim_t *m, dim_t *n, double *alpha, double *a, dim_t *lda,
            double *x, dim_t *incx, double *beta, double *y, dim_t *incy);
void cgemv_(char *trans, dim_t *m, dim_t *n, void *alpha, void *a, dim_t *lda,
            void *x, dim_t *incx, void *beta, void *y, dim_t *incy);
void zgemv_(char *trans, dim_t *m, dim_t *n, void *alpha, void *a, dim_t *lda,
            void *x, dim_t *incx, void *beta, void *y, dim_t *incy);


// trmm
void strmm_(char *side, char *uplo, char *transa, char *diag, dim_t *m, dim_t *n,
            float *alpha, float *a, inc_t *lda, float *b, inc_t *ldb);	
void dtrmm_(char *side, char *uplo, char *transa, char *diag, dim_t *m, dim_t *n,
            double *alpha, double *a, inc_t *lda, double *b, inc_t *ldb);	
void ctrmm_(char *side, char *uplo, char *transa, char *diag, dim_t *m, dim_t *n,
            void *alpha, void *a, inc_t *lda, void *b, inc_t *ldb);	
void ztrmm_(char *side, char *uplo, char *transa, char *diag, dim_t *m, dim_t *n,
            void *alpha, void *a, inc_t *lda, void *b, inc_t *ldb);	


// trmv
void strmv_(char *uplo, char *transa, char *diag, dim_t *n,
            float *a, inc_t *lda, float *b, inc_t *ldb);	
void dtrmv_(char *uplo, char *transa, char *diag, dim_t *n,
            double *a, inc_t *lda, double *b, inc_t *ldb);	
void ctrmv_(char *uplo, char *transa, char *diag, dim_t *n,
            void *a, inc_t *lda, void *b, inc_t *ldb);	
void ztrmv_(char *uplo, char *transa, char *diag, dim_t *n,
            void *a, inc_t *lda, void *b, inc_t *ldb);	

// swap
void sswap_(dim_t *n, float *sx, dim_t *incx, float *sy, dim_t *incy);
void dswap_(dim_t *n, double *sx, dim_t *incx, double *sy, dim_t *incy);
void cswap_(dim_t *n, void *sx, dim_t *incx, void *sy, dim_t *incy);
void zswap_(dim_t *n, void *sx, dim_t *incx, void *sy, dim_t *incy);

// axpy
void saxpy_(dim_t *n, float *alpha, float *sx, dim_t *incx, float *sy, dim_t *incy);
void daxpy_(dim_t *n, double *alpha, double *sx, dim_t *incx, double *sy, dim_t *incy);
void caxpy_(dim_t *n, void *alpha, void *sx, dim_t *incx, void *sy, dim_t *incy);
void zaxpy_(dim_t *n, void *alpha, void *sx, dim_t *incx, void *sy, dim_t *incy);

// dot
float sdot_(dim_t *n, float *sx, dim_t *incx, float *sy, dim_t *incy);
double ddot_(dim_t *n, double *sx, dim_t *incx, double *sy, dim_t *incy);
scomplex cdotc_(dim_t *n, void *sx, dim_t *incx, void *sy, dim_t *incy);
dcomplex zdotc_(dim_t *n, void *sx, dim_t *incx, void *sy, dim_t *incy);

// }
#endif

// LAPACK {

// lacpy
void slacpy_(char *uplo, dim_t *m, dim_t *n, float *a, inc_t *lda, float *b, inc_t *ldb);
void dlacpy_(char *uplo, dim_t *m, dim_t *n, double *a, inc_t *lda, double *b, inc_t *ldb);
void clacpy_(char *uplo, dim_t *m, dim_t *n, void *a, inc_t *lda, void *b, inc_t *ldb);
void zlacpy_(char *uplo, dim_t *m, dim_t *n, void *a, inc_t *lda, void *b, inc_t *ldb);

// trtri
void strtri_(char *uplo, char *diag, dim_t *n, float *a, inc_t *lda, int *info);
void dtrtri_(char *uplo, char *diag, dim_t *n, double *a, inc_t *lda, int *info);
void ctrtri_(char *uplo, char *diag, dim_t *n, void *a, inc_t *lda, int *info);
void ztrtri_(char *uplo, char *diag, dim_t *n, void *a, inc_t *lda, int *info);

// }

// PFAPACK77 {

void ssktrf_(char *uplo, char *mode, dim_t *n, float *a, dim_t *lda, int *ipiv, float *work, dim_t *lwork, int *info);
void dsktrf_(char *uplo, char *mode, dim_t *n, double *a, dim_t *lda, int *ipiv, double *work, dim_t *lwork, int *info);
void csktrf_(char *uplo, char *mode, dim_t *n, void *a, dim_t *lda, int *ipiv, void *work, dim_t *lwork, int *info);
void zsktrf_(char *uplo, char *mode, dim_t *n, void *a, dim_t *lda, int *ipiv, void *work, dim_t *lwork, int *info);

// }

#ifdef __cplusplus
}
#endif

