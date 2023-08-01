/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <Eigen/Dense>
#include <cstdint>
#ifndef EIGEN_USE_LAPACKE // This lets Eigen/Dense include LAPACK definitions.
#include <Eigen/src/misc/lapacke.h>
#endif
#ifdef _BLIS
#define BLIS_DISABLE_BLAS // Screen BLAS definition even if BLIS' built with them.
#include "blis/blis.h"
#else
using Eigen::scomplex;
using Eigen::dcomplex;
#endif
#include "error.hh"
#include "fortran_pfapack.h"
using f77_int = lapack_int;
using ccscmplx = std::complex<float>;
using ccdcmplx = std::complex<double>;

namespace l2e
{

template <typename T>
struct le_mat_t
{
    const f77_int rows, cols; ///< Matrix size.
    const f77_int rs, cs; ///< Strides.

    T *const data;

protected:
    template <typename Mat>
    static f77_int _rs(const Mat &mat)
    {
        f77_int rs = 0;
        if (mat.rows() > 1)
            rs = &(mat(1, 0)) - &(mat(0, 0));
        return rs;
    }
    template <typename Mat>
    static f77_int _cs(const Mat &mat)
    {
        f77_int cs = 0;
        if (mat.cols() > 1)
            cs = &(mat(0, 1)) - &(mat(0, 0));
        return cs;
    }

public:
    le_mat_t(const f77_int &rows_, const f77_int &cols_,
             const f77_int &rs_, const f77_int &cs_, T *const data_)
    : rows(rows_), cols(cols_), rs(rs_), cs(cs_), data(data_) { }

    template <int RowsAtCompilation, int ColsAtCompilation>
    le_mat_t(Eigen::Matrix<T, RowsAtCompilation, ColsAtCompilation> &mat)
    : rows(mat.rows()), cols(mat.cols()), rs(_rs(mat)), cs(_cs(mat)), data(&(mat(0, 0))) { }

    template <typename Mat, int BlockRowsAtCompil, int BlockColsAtCompil, bool BlockProps>
    le_mat_t(Eigen::Block<Mat, BlockRowsAtCompil, BlockColsAtCompil, BlockProps> mat) // TODO: Check for data copying.
    : rows(mat.rows()), cols(mat.cols()), rs(_rs(mat)), cs(_cs(mat)), data(&(mat(0, 0))) { }

    template <int RowsAtCompilation, int ColsAtCompilation, int Alignment, typename Stride>
    le_mat_t(Eigen::Map<Eigen::Matrix<T, RowsAtCompilation, ColsAtCompilation>, Alignment, Stride> &mat)
    : rows(mat.rows()), cols(mat.cols()), rs(_rs(mat)), cs(_cs(mat)), data(&(mat(0, 0))) { }

    template <int Alignment = 0>
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Alignment,
               Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > to_eigen() const
    {
        using matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
        using stride_t = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
        return Eigen::Map<matrix_t, Alignment, stride_t>(data, rows, cols, stride_t(cs, rs));
    }

    le_mat_t<T> transpose() const
    { return le_mat_t<T>(cols, rows, cs, rs, data); }
};

}


// ---- BLAS/LAPACK routines used in our settings. -----


// [LAPACK] lacpy
#define BLALINK_PARAMS(T) const char &uploA, \
                          const l2e::le_mat_t<T> &A, \
                                l2e::le_mat_t<T> &B
template <typename T>
inline void lacpy( BLALINK_PARAMS(T) );
#define BLALINK_MAC(T, ctype, cchar) \
    template <> inline void lacpy<T>( BLALINK_PARAMS(T) ) \
    { \
        assert_(A.rs == 1 && B.rs == 1, "BLAS requires column-major."); \
        /* Eigen-vendored lapack.h misses const qualifier. Do force conversion here. */ \
        cchar##lacpy_((char *)&uploA, (f77_int *)&A.rows, (f77_int *)&A.cols, \
                      (ctype *)A.data, (f77_int *)&A.cs, \
                      (ctype *)B.data, (f77_int *)&B.cs); \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, _Complex float,  c )
BLALINK_MAC( ccdcmplx, _Complex double, z )
#undef BLALINK_MAC
#undef BLALINK_PARAMS


// [LAPACK] trtri
#define BLALINK_PARAMS(T) const char &uploA, \
                          const char &diagA, \
                          l2e::le_mat_t<T> &A
template <typename T>
inline f77_int trtri( BLALINK_PARAMS(T) );
#define BLALINK_MAC(T, ctype, cchar) \
    template <> inline f77_int trtri( BLALINK_PARAMS(T) ) \
    { \
        assert_(A.rs == 1, "BLAS requires column-major."); \
        f77_int info; \
        /* Eigen-vendored lapack.h misses const qualifier. Do force conversion here. */ \
        cchar##trtri_((char *)&uploA, (char *)&diagA, (f77_int *)&A.rows, \
                      (ctype *)A.data, (f77_int *)&A.cs, &info); \
        return info; \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, _Complex float,  c )
BLALINK_MAC( ccdcmplx, _Complex double, z )
#undef BLALINK_MAC
#undef BLALINK_PARAMS


// [PFAPACK] sktrf
#define BLALINK_PARAMS(T) const char &uplo, \
                          const char &mode, \
                          l2e::le_mat_t<T> &A, \
                          f77_int *ipiv, \
                          T *work, f77_int lwork
template <typename T>
inline f77_int sktrf( BLALINK_PARAMS(T) );
#define BLALINK_MAC(T, ctype, cchar) \
    template <> inline f77_int sktrf( BLALINK_PARAMS(T) ) \
    { \
        assert_(A.rs == 1, "BLAS requires column-major."); \
        f77_int info; \
        cchar##sktrf_(&uplo, &mode, &A.rows, (ctype *)A.data, &A.cs, ipiv, work, &lwork, &info); \
        return info; \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, floatcmplx,  c )
BLALINK_MAC( ccdcmplx, doublecmplx, z )
#undef BLALINK_MAC
#undef BLALINK_PARAMS


// [PFAPACK] skpfa
#define BLALINK_PARAMS(T) const char &uplo, \
                          const char &mthd, \
                          l2e::le_mat_t<T> &A, \
                          T *pfa, f77_int *ipiv, \
                          T *work, f77_int lwork
template <typename T>
inline f77_int skpfa( BLALINK_PARAMS(T) );
#define BLALINK_MAC(T, ctype, cchar) \
    template <> inline f77_int skpfa( BLALINK_PARAMS(T) ) \
    { \
        assert_(A.rs == 1, "BLAS requires column-major."); \
        int info; \
        cchar##skpfa_(&uplo, &mthd, &A.rows, (ctype *)A.data, &A.cs, pfa, ipiv, work, &lwork, &info); \
        return info; \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
#undef BLALINK_MAC
#define BLALINK_MAC(T, ctype, rtype, cchar) \
    template <> inline f77_int skpfa( BLALINK_PARAMS(T) ) \
    { \
        assert_(A.rs == 1, "BLAS requires column-major."); \
        rtype *rwork = 0; \
        if (mthd != 'p' && mthd != 'P') \
            rwork = new rtype[lwork]; \
        int info; \
        cchar##skpfa_(&uplo, &mthd, &A.rows, A.data, &A.cs, pfa, ipiv, work, &lwork, rwork, &info); \
        delete[] rwork; \
        return info; \
    }
BLALINK_MAC( ccscmplx, floatcmplx,  float,  c )
BLALINK_MAC( ccdcmplx, doublecmplx, double, z )
#undef BLALINK_MAC
#undef BLALINK_PARAMS


// [Non-standard] gemmt
#define BLALINK_PARAMS(T) const char &uploC, \
                          const char &transA, \
                          const char &transB, \
                          const T &alpha, \
                          const l2e::le_mat_t<T> &A, \
                          const l2e::le_mat_t<T> &B, \
                          const T &beta, \
                                l2e::le_mat_t<T> &C
template <typename T>
inline void gemmt( BLALINK_PARAMS(T) );
#if !defined(_BLIS) && !defined(_BLAS_HAS_GEMMT)
  #include "eigen_gemmt.tcc"
#else
#ifdef _BLIS
#define BLALINK_MAC(T, ctype, cchar) \
    template <> inline void gemmt<T>( BLALINK_PARAMS(T) ) \
    { \
        uplo_t bli_uploC; \
        trans_t bli_transA, bli_transB; \
        bli_param_map_netlib_to_blis_uplo( uploC, &bli_uploC ); \
        bli_param_map_netlib_to_blis_trans( transA, &bli_transA ); \
        bli_param_map_netlib_to_blis_trans( transB, &bli_transB ); \
        bli_##cchar##gemmt(bli_uploC, \
                           bli_transA, \
                           bli_transB, \
                           (transA == 'T' || transA == 't') ? A.cols : A.rows, \
                           (transA == 'T' || transA == 't') ? A.rows : A.cols, \
                           (const ctype *)&alpha, \
                           (const ctype *)A.data, A.rs, A.cs, \
                           (const ctype *)B.data, B.rs, B.cs, \
                           (const ctype *)&beta, \
                           (ctype       *)C.data, C.rs, C.cs); \
    }
#else
  // ifdef _BLAS_HAS_GEMMT
#define BLALINK_MAC(T, ctype, cchar) \
    template <> inline void gemmt<T>( BLALINK_PARAMS(T) ) \
    { \
        assert_(A.rs == 1 && B.rs == 1 && C.rs == 1, "BLAS requires column-major."); \
        cchar##gemmt_(&uploC, \
                      &transA, &transB, \
                      (transA == 'T' || transA == 't') ? &A.cols : &A.rows, \
                      (transA == 'T' || transA == 't') ? &A.rows : &A.cols, \
                      (const ctype *)&alpha, \
                      (const ctype *)A.data, &A.cs, \
                      (const ctype *)B.data, &B.cs, \
                      (const ctype *)&beta, \
                      (ctype       *)C.data, &C.cs); \
    }
extern "C" {
#define GEMMT_F77_PARAMS(T) const char *uploC, \
                            const char *transA, \
                            const char *transB, \
                            const f77_int *m, \
                            const f77_int *k, \
                            const T *alpha, \
                            const T *A, const f77_int *ldA, \
                            const T *B, const f77_int *ldB, \
                            const T *beta, \
                                  T *C, const f77_int *ldC
void sgemmt_( GEMMT_F77_PARAMS(float)    );
void dgemmt_( GEMMT_F77_PARAMS(double)   );
void cgemmt_( GEMMT_F77_PARAMS(scomplex) );
void zgemmt_( GEMMT_F77_PARAMS(dcomplex) );
#undef GEMMT_F77_PARAMS
}
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC
#endif
#undef BLALINK_PARAMS
