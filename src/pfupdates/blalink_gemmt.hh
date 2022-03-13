/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "blalink.hh"


// gemmt is not part of BLAS standard.
// It's currently know as exposed by BLIS and MKL.
template <typename T>
inline void gemmt(uplo_t uploc,
                  trans_t transa, trans_t transb,
                  dim_t m, dim_t k,
                  T alpha,
                  T *a, inc_t lda,
                  T *b, inc_t ldb,
                  T beta,
                  T *c, inc_t ldc);
#if !( defined(BLAS_EXTERNAL) && defined(MKL) )
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemmt<cctype>(uplo_t uploc, \
                                          trans_t transa, trans_t transb, \
                                          dim_t m, dim_t k, \
                                          cctype alpha, \
                                          cctype *a, inc_t lda, \
                                          cctype *b, inc_t ldb, \
                                          cctype beta, \
                                          cctype *c, inc_t ldc) \
    { \
        bli_##cchar##gemmt(uploc, \
                           transa, transb, \
                           m, k, \
                           (ctype *)&alpha, \
                           (ctype *)a, 1, lda, \
                           (ctype *)b, 1, ldb, \
                           (ctype *)&beta, \
                           (ctype *)c, 1, ldc); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemmt<cctype>(uplo_t uploc, \
                                          trans_t transa, trans_t transb, \
                                          dim_t m, dim_t k, \
                                          cctype alpha, \
                                          cctype *a, inc_t lda, \
                                          cctype *b, inc_t ldb, \
                                          cctype beta, \
                                          cctype *c, inc_t ldc) \
    { \
        char ul = uplo2char(uploc); \
        char ta = trans2char(transa); \
        char tb = trans2char(transb); \
        cchar##gemmt_(&ul, \
                      &ta, &tb, \
                      &m, &k, \
                      (ctype *)&alpha, \
                      (ctype *)a, &lda, \
                      (ctype *)b, &ldb, \
                      (ctype *)&beta, \
                      (ctype *)c, &ldc); \
    }
void sgemmt_(char *uplo, char *transa, char *transb, dim_t *m, dim_t *k, float *alpha,
             float *a, dim_t *lda, float *b, dim_t *ldb, float *beta, float *c, dim_t *ldc);
void dgemmt_(char *uplo, char *transa, char *transb, dim_t *m, dim_t *k, double *alpha,
             double *a, dim_t *lda, double *b, dim_t *ldb, double *beta, double *c, dim_t *ldc);
void cgemmt_(char *uplo, char *transa, char *transb, dim_t *m, dim_t *k, void *alpha,
             void *a, dim_t *lda, void *b, dim_t *ldb, void *beta, void *c, dim_t *ldc);
void zgemmt_(char *uplo, char *transa, char *transb, dim_t *m, dim_t *k, void *alpha,
             void *a, dim_t *lda, void *b, dim_t *ldb, void *beta, void *c, dim_t *ldc);
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC

