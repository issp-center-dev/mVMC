/**
 * \file blalink.hh
 * C++ to C type wrapper for BLIS type-apis.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <complex>
#include <cassert>
typedef std::complex<float>  ccscmplx;
typedef std::complex<double> ccdcmplx;
// BLIS definitions.
#include "blis.h"
// BLAS definitions
#include "blalink_fort.h"

// error info
typedef enum {
    Pfaffine_SIGN_ERR = 1,
    Pfaffine_OUT_OF_BOUND,
    Pfaffine_BAD_SCRATCHPAD,
    Pfaffine_BAD_REPRESENTATION,
    Pfaffine_INTEGER_OVERFLOW,
    Pfaffine_DOUBLE_NAN_DETECTED,
    Pfaffine_NOT_IMPLEMNTED,
    Pfaffine_NUM_ERROR_TYPE
} skpfa_error_t;
inline signed err_info(signed type, signed pos)
{ return type * Pfaffine_NUM_ERROR_TYPE + pos; }

// gemm
template <typename T>
inline void gemm(trans_t transa, trans_t transb,
                 dim_t m, dim_t n, dim_t k,
                 T alpha,
                 T *a, inc_t lda,
                 T *b, inc_t ldb,
                 T beta,
                 T *c, inc_t ldc);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemm<cctype>(trans_t transa, trans_t transb, \
                                         dim_t m, dim_t n, dim_t k, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb, \
                                         cctype beta, \
                                         cctype *c, inc_t ldc) \
    { \
        bli_##cchar##gemm(transa, transb, \
                          m, n, k, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)b, 1, ldb, \
                          (ctype *)&beta, \
                          (ctype *)c, 1, ldc); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemm<cctype>(trans_t transa, trans_t transb, \
                                         dim_t m, dim_t n, dim_t k, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb, \
                                         cctype beta, \
                                         cctype *c, inc_t ldc) \
    { \
        char ta = trans2char(transa), \
             tb = trans2char(transb); \
        cchar##gemm_(&ta, &tb, &m, &n, &k, \
                     (ctype *)&alpha, \
                     (ctype *)a, &lda, \
                     (ctype *)b, &ldb, \
                     (ctype *)&beta, \
                     (ctype *)c, &ldc); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// ger
template <typename T>
inline void ger(dim_t m, dim_t n,
                T alpha,
                T *x, inc_t incx,
                T *y, inc_t incy,
                T *a, inc_t lda);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline void ger<cctype>(dim_t m, dim_t n, \
                                        cctype alpha, \
                                        cctype *x, inc_t incx, \
                                        cctype *y, inc_t incy, \
                                        cctype *a, inc_t lda) \
    { \
        bli_##cchar##ger(BLIS_NO_CONJUGATE, \
                         BLIS_NO_CONJUGATE, \
                         m, n, \
                         (ctype *)&alpha, \
                         (ctype *)x, incx, \
                         (ctype *)y, incy, \
                         (ctype *)a, 1, lda); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline void ger<cctype>(dim_t m, dim_t n, \
                                        cctype alpha, \
                                        cctype *x, inc_t incx, \
                                        cctype *y, inc_t incy, \
                                        cctype *a, inc_t lda) \
    { \
        cchar##cfunc##_(&m, &n, \
                        (ctype *)&alpha, \
                        (ctype *)x, &incx, \
                        (ctype *)y, &incy, \
                        (ctype *)a, &lda); \
    }
#endif
BLALINK_MAC( float,    float,    s, ger  )
BLALINK_MAC( double,   double,   d, ger  )
BLALINK_MAC( ccscmplx, scomplex, c, geru )
BLALINK_MAC( ccdcmplx, dcomplex, z, geru )
#undef BLALINK_MAC


// gemv
template <typename T>
inline void gemv(trans_t trans,
                 dim_t m, dim_t n,
                 T alpha,
                 T *a, inc_t lda,
                 T *x, inc_t incx,
                 T beta,
                 T *y, inc_t incy);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemv<cctype>(trans_t trans, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *x, inc_t incx, \
                                         cctype beta, \
                                         cctype *y, inc_t incy) \
    { \
        bli_##cchar##gemv(trans, \
                          BLIS_NO_CONJUGATE, \
                          m, n, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)x, incx, \
                          (ctype *)&beta, \
                          (ctype *)y, incy); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemv<cctype>(trans_t trans, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *x, inc_t incx, \
                                         cctype beta, \
                                         cctype *y, inc_t incy) \
    { \
        char t = trans2char(trans); \
        cchar##gemv_(&t, &m, &n, \
                     (ctype *)&alpha, \
                     (ctype *)a, &lda, \
                     (ctype *)x, &incx, \
                     (ctype *)&beta, \
                     (ctype *)y, &incy); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// trmm
template <typename T>
inline void trmm(side_t sidea, 
                 uplo_t uploa, 
                 trans_t transa, 
                 diag_t diaga,
                 dim_t m, dim_t n,
                 T alpha, 
                 T *a, inc_t lda, 
                 T *b, inc_t ldb);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void trmm<cctype>(side_t sidea, uplo_t uploa, trans_t transa, diag_t diaga, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb) \
    { \
        bli_##cchar##trmm(sidea, uploa, transa, diaga, \
                          m, n, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)b, 1, ldb); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void trmm<cctype>(side_t sidea, uplo_t uploa, trans_t transa, diag_t diaga, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb) \
    { \
        char ul = uplo2char(uploa), \
             si = side2char(sidea), \
             tr = trans2char(transa), \
             dg = diag2char(diaga); \
        cchar##trmm_(&si, &ul, &tr, &dg, \
                     &m, &n, \
                     (ctype *)&alpha, \
                     (ctype *)a, &lda, \
                     (ctype *)b, &ldb); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// trmv
template <typename T>
inline void trmv(uplo_t uploa,
                 trans_t transa,
                 dim_t m,
                 T alpha,
                 T *a, inc_t lda,
                 T *x, inc_t incx);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void trmv<cctype>(uplo_t uploa, trans_t transa, dim_t m, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *x, inc_t incx) \
    { \
        bli_##cchar##trmv(uploa, transa, \
                          BLIS_NONUNIT_DIAG, m, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)x, incx); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void trmv<cctype>(uplo_t uploa, trans_t transa, dim_t m, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *x, inc_t incx) \
    { \
        char ul = uplo2char(uploa), \
             tr = trans2char(transa), \
             dg = 'N'; \
        /*
         * TRMV has alpha in its templated interface, which is absent in BLAS.
         */ \
        assert(alpha == cctype(1.0)); \
        cchar##trmv_(&ul, &tr, &dg, &m, \
                     (ctype *)a, &lda, \
                     (ctype *)x, &incx); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// swap
template <typename T>
inline void swap(dim_t n,
                 T *x, inc_t incx,
                 T *y, inc_t incy);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void swap<cctype>(dim_t n, \
                                         cctype *x, inc_t incx, \
                                         cctype *y, inc_t incy) \
    { \
        bli_##cchar##swapv(n, \
                           (ctype *)x, incx, \
                           (ctype *)y, incy); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void swap<cctype>(dim_t n, \
                                         cctype *x, inc_t incx, \
                                         cctype *y, inc_t incy) \
    { \
        cchar##swap_(&n, \
                     (ctype *)x, &incx, \
                     (ctype *)y, &incy); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// axpy
template <typename T>
inline void axpy(dim_t n,
                 T alpha,
                 T *x, inc_t incx,
                 T *y, inc_t incy);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void axpy<cctype>(dim_t n, \
                                         cctype alpha, \
                                         cctype *x, inc_t incx, \
                                         cctype *y, inc_t incy) \
    { \
        bli_##cchar##axpyv(BLIS_NO_CONJUGATE, n, \
                           (ctype *)&alpha, \
                           (ctype *)x, incx, \
                           (ctype *)y, incy); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void axpy<cctype>(dim_t n, \
                                         cctype alpha, \
                                         cctype *x, inc_t incx, \
                                         cctype *y, inc_t incy) \
    { \
        cchar##axpy_(&n, \
                     (ctype *)&alpha, \
                     (ctype *)x, &incx, \
                     (ctype *)y, &incy); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC

// dot
template <typename T>
inline T dot(dim_t n,
             T *sx, inc_t incx,
             T *sy, inc_t incy);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline cctype dot<cctype>(dim_t n, \
                                          cctype *sx, inc_t incx, \
                                          cctype *sy, inc_t incy) \
    { \
        cctype rho; \
        bli_##cchar##dotv(BLIS_NO_CONJUGATE, \
                          BLIS_NO_CONJUGATE, \
                          n, \
                          (ctype *)sx, incx, \
                          (ctype *)sy, incy, \
                          (ctype *)&rho); \
        return rho; \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline cctype dot<cctype>(dim_t n, \
                                          cctype *sx, inc_t incx, \
                                          cctype *sy, inc_t incy) \
    { \
        ctype rho = cchar##cfunc##_(&n, \
                                    (ctype *)sx, &incx, \
                                    (ctype *)sy, &incy); \
        return *((cctype *)&rho); \
    }
#endif
BLALINK_MAC( float,    float,    s, dot  )
BLALINK_MAC( double,   double,   d, dot  )
#if defined(BLAS_EXTERNAL) && defined(F77_COMPLEX_RET_INTEL)
// Intel style complex return.
#undef BLALINK_MAC
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline cctype dot<cctype>(dim_t n, \
                                          cctype *sx, inc_t incx, \
                                          cctype *sy, inc_t incy) \
    { \
        cctype rho; \
        cchar##cfunc##_((ctype *)&rho, &n, \
                        (ctype *)sx, &incx, \
                        (ctype *)sy, &incy); \
        return rho; \
    }
#endif
BLALINK_MAC( ccscmplx, scomplex, c, dotu )
BLALINK_MAC( ccdcmplx, dcomplex, z, dotu )
#undef BLALINK_MAC


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
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// [LAPACK] lacpy
template <typename T>
inline void lacpy(uplo_t uploa,
                  dim_t m, dim_t n,
                  T *a, inc_t lda,
                  T *b, inc_t ldb);
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void lacpy<cctype>(uplo_t uploa, \
                                          dim_t m, dim_t n, \
                                          cctype *a, inc_t lda, \
                                          cctype *b, inc_t ldb) \
    { \
        char ul = uplo2char(uploa); \
        cchar##lacpy_(&ul, &m, &n, a, &lda, b, &ldb); \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC

// [LAPACK] trtri
template <typename T>
inline int trtri(uplo_t uploa,
                 diag_t diaga,
                 dim_t n,
                 T *a, inc_t lda);
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline int trtri(uplo_t uploa, \
                                 diag_t diaga, \
                                 dim_t n, \
                                 cctype *a, inc_t lda) \
    { \
        char ul = uplo2char(uploa); \
        char dg = diag2char(diaga); \
        int info; \
        cchar##trtri_(&ul, &dg, &n, a, &lda, &info); \
        return info; \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC

// [PFAPACK] sktrf
template <typename T>
inline int sktrf(uplo_t uploa,
                 dim_t n,
                 T *a, inc_t lda,
                 int *ipiv,
                 T *work, dim_t lwork);
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline int sktrf(uplo_t uploa, \
                                 dim_t n, \
                                 cctype *a, inc_t lda, \
                                 int *ipiv, \
                                 cctype *work, dim_t lwork) \
    { \
        char ul = uplo2char(uploa); \
        char mode = 'N'; \
        int info; \
        cchar##sktrf_(&ul, &mode, &n, a, &lda, ipiv, work, &lwork, &info); \
        return info; \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC

// [PFAPACK] skpfa
template <typename T>
inline int skpfa(uplo_t uploa,
                 dim_t n,
                 T *a, inc_t lda,
                 T *pfa, int *ipiv,
                 T *work, dim_t lwork);
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline int skpfa(uplo_t uploa, \
                                 dim_t n, \
                                 cctype *a, inc_t lda, \
                                 cctype *pfa, int *ipiv, \
                                 cctype *work, dim_t lwork) \
    { \
        char ul = uplo2char(uploa); \
        char mthd = 'P'; \
        int info; \
        cchar##skpfa_(&ul, &mthd, &n, a, &lda, pfa, ipiv, work, &lwork, &info); \
        return info; \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
#undef BLALINK_MAC
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline int skpfa(uplo_t uploa, \
                                 dim_t n, \
                                 cctype *a, inc_t lda, \
                                 cctype *pfa, int *ipiv, \
                                 cctype *work, dim_t lwork) \
    { \
        char ul = uplo2char(uploa); \
        char mthd = 'P'; \
        int info; \
        cchar##skpfa_(&ul, &mthd, &n, a, &lda, pfa, ipiv, work, &lwork, 0, &info); \
        return info; \
    }
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC

