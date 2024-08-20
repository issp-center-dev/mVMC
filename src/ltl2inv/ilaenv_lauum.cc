/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blalink.hh"
#include "ilaenv.h"
#include "ilaenv_lauum.hh"

#define EXPANDMAC(cctype, name) \
template <> int ilaenv_lauum<cctype>(uplo_t uplo, int n) \
{ \
    char uplo_ = uplo2char(uplo); \
    int ispec = 1; \
    int dummy = 0; \
    int n_ = n; \
    return ilaenv_(&ispec, #name, &uplo_, &n_, &dummy, &dummy, &dummy); \
}
EXPANDMAC( float,    SLAUUM )
EXPANDMAC( double,   DLAUUM )
EXPANDMAC( ccscmplx, CLAUUM )
EXPANDMAC( ccdcmplx, ZLAUUM )
#undef EXPANDMAC

