/**
 * \file skslc.tcc
 * Get a slice from antisymmetric matrix.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "colmaj.hh"

template<typename T>
inline signed skslc(uplo_t uploA,
                    unsigned n,
                    unsigned i,
                    T *x,
                    T *_A, unsigned ldA)
{
    colmaj<T> A(_A, ldA);
    x[i] = 0.0;
    switch (uploA) {
    case BLIS_UPPER:
        for (unsigned j = 0; j < i; ++j)
            x[j] =  A(j, i);
        for (unsigned j = i+1; j < n; ++j)
            x[j] = -A(i, j);
        return 0;

    default:
        std::cerr << "SKSLC: Lower triangular storage not implemented. Sorry." << std::endl;
        return err_info(Pfaffine_NOT_IMPLEMNTED, 0);
    }
}

