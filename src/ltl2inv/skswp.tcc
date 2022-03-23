/**
 * \file skswp.tcc
 * Swap 2 row-columns of an antisymmetric matrix.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "colmaj.hh"
#include "blalink.hh"
#include <iostream>

template <typename T> 
inline signed skswp(uplo_t uploa,
                    dim_t m,
                    colmaj<T> &A,
                    dim_t s, dim_t t) 
{
    // Swap s and t if t<s.
    if (s > t) {
        dim_t r = s;
        s = t;
        t = r;
    }

    switch (uploa) {
    case BLIS_UPPER:
        if (s)
            swap<T>(s, &A(0, s), 1, &A(0, t), 1);
        A(s, t) *= -1.0;
        if (t > s+1) {
            // swap<T>(t-s-1, &A(s+1, t), 1, &A(s, s+1), A.ld);
            for (dim_t j = s+1; j < t; ++j) {
                // A(j, t) *= -1;
                // A(s, j) *= -1;
                T r = -A(j, t);
                A(j, t) = -A(s, j);
                A(s, j) = r;
            }
        }
        if (t+1 < m)
            swap<T>(m-1-t, &A(s, t+1), A.ld, &A(t, t+1), A.ld);
        return 0;

    case BLIS_LOWER:
        if (s)
            swap<T>(s, &A(s, 0), A.ld, &A(t, 0), A.ld);
        A(t, s) *= -1.0;
        if (t > s+1) {
            for (dim_t j = s+1; j < t; ++j) {
                T r = -A(t, j);
                A(t, j) = -A(j, s);
                A(j, s) = r;
            }
        }
        if (t+1 < m)
            swap<T>(m-1-t, &A(t+1, s), 1, &A(t+1, t), 1);
        return 0;

    default:
        std::cerr << "SKSWP: Unsupported UpLo input." << std::endl;
        return err_info(Pfaffine_NOT_IMPLEMNTED, 0);
    }
}

