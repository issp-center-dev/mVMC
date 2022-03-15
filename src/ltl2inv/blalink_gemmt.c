/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blis.h"
#include "blalink_fort.h"

#if !( defined(BLAS_EXTERNAL) && defined(MKL) )
// Redefine GEMMT Fortran from BLIS (which is built without BLAS API).
// Required by Pfapack77.
#define BLAGEN_MAC(cctype, ctype, cchar) \
    void cchar##gemmt_(char *uplo_, char *transa_, char *transb_, dim_t *m, dim_t *k, ctype *alpha, \
                       ctype *a, dim_t *lda, ctype *b, dim_t *ldb, ctype *beta, ctype *c, dim_t *ldc) \
    { \
        uplo_t uploc; \
        trans_t transa; \
        trans_t transb; \
        bli_param_map_netlib_to_blis_uplo(*uplo_, &uploc) ; \
        bli_param_map_netlib_to_blis_trans(*transa_, &transa) ; \
        bli_param_map_netlib_to_blis_trans(*transb_, &transb) ; \
        bli_##cchar##gemmt(uploc, \
                           transa, transb, \
                           *m, *k, \
                           alpha, \
                           a, 1, &lda, \
                           b, 1, &ldb, \
                           beta, \
                           c, 1, &ldc); \
    }
BLAGEN_MAC( float,    float,    s )
BLAGEN_MAC( double,   double,   d )
BLAGEN_MAC( ccscmplx, scomplex, c )
BLAGEN_MAC( ccdcmplx, dcomplex, z )
#undef BLAGEN_MAC
#endif
