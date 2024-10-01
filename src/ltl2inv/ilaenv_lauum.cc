/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blalink.hh"
#include "ilaenv.h"
#include "ilaenv_lauum.hh"

/*
 * sloppy version exploits that 
 * ILAENV(ispec=1, xLAUUM, "", dummy, dummy, dummy, dummy)
 * returns 64 for x=S,D,C,Z .
 */

#undef SLOPPY_ILAENV
// #define SLOPPY_ILAENV

/*
#define EXPANDMAC(cctype, name) \
template <> int ilaenv_lauum<cctype>(uplo_t uplo, int n) \
{ \
    char uplo_ = uplo2char(uplo); \
    int ispec = 1; \
    int dummy = 0; \
    return ilaenv_(&ispec, #name, &uplo_, &n, &dummy, &dummy, &dummy); \
}
*/

#ifndef SLOPPY_ILAENV
#define EXPANDMAC(cctype, name) \
template <> int ilaenv_lauum<cctype>(uplo_t uplo, int n) \
{ \
    char uplo_ = uplo2char(uplo); \
    int ispec = 1; \
    int dummy = -1; \
    return ilaenv_wrap(ispec, #name, &uplo_, n, dummy, dummy, dummy); \
}
#else
#define EXPANDMAC(cctype, name) \
template <> int ilaenv_lauum<cctype>(uplo_t uplo, int n) \
{ \
    return 64; \
}
#endif

EXPANDMAC( float,    SLAUUM )
EXPANDMAC( double,   DLAUUM )
EXPANDMAC( ccscmplx, CLAUUM )
EXPANDMAC( ccdcmplx, ZLAUUM )
#undef EXPANDMAC

