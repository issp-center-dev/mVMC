/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "lapack2eigen.hh"


#if !defined(_BLAS_HAS_GEMMT)
// Import the customized GEMMT impl. to the BLAS API.
// Required by Pfapack77.
#define BLAGEN_MAC(cctype, ctype, cchar) \
    extern "C" \
    void cchar##gemmt_(const char *uplo_, \
                       const char *transa_, \
                       const char *transb_, \
                       const int *m, const int *k, \
                       const ctype *alpha, \
                       ctype *a, const int *lda, \
                       ctype *b, const int *ldb, \
                       const ctype *beta, ctype *c, const int *ldc) \
    { \
        bool tA = (*transa_ == 'T' || *transa_ == 't'); \
        bool tB = (*transb_ == 'T' || *transb_ == 't'); \
        l2e::le_mat_t<cctype> A(tA ? *k : *m, tA ? *m : *k, 1, *lda, (cctype *)a); \
        l2e::le_mat_t<cctype> B(tB ? *k : *m, tB ? *m : *k, 1, *ldb, (cctype *)b); \
        l2e::le_mat_t<cctype> C(*m, *m, 1, *ldc, (cctype *)c); \
        gemmt<cctype>(*uplo_, *transa_, *transb_, *(cctype *)alpha, A, B, *(cctype *)beta, C); \
    }
BLAGEN_MAC( float,    float,    s )
BLAGEN_MAC( double,   double,   d )
BLAGEN_MAC( ccscmplx, scomplex, c )
BLAGEN_MAC( ccdcmplx, dcomplex, z )
#undef BLAGEN_MAC
#endif
