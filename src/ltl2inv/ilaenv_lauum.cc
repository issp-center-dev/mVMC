/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "ilaenv.h"
#include "ilaenv_lauum.hh"
#include <complex>
using ccscmplx = std::complex<float>;
using ccdcmplx = std::complex<double>;


#define EXPANDMAC(cctype, name) \
template <> int ilaenv_lauum<cctype>(const char uplo, const int n) \
{ \
    int ispec = 1; \
    int dummy = 0; \
    return ilaenv_(&ispec, #name, &uplo, &n, &dummy, &dummy, &dummy); \
}
EXPANDMAC( float,    SLAUUM )
EXPANDMAC( double,   DLAUUM )
EXPANDMAC( ccscmplx, CLAUUM )
EXPANDMAC( ccdcmplx, ZLAUUM )
#undef EXPANDMAC

