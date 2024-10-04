/**
 * \copyright Copyright (c) Dept. Phys., Univ. Tokyo
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "updated_tdi.tcc"
#include <cstdio>
#include <chrono>
#include <vector>

struct info_t {
    double err;
    long long tim;
};

template <typename T>
info_t test_tdi_update(dim_t nsite, std::vector<dim_t> cfg, std::vector<dim_t> from, std::vector<dim_t> to)
{
    dim_t nelec = cfg.size();
    dim_t mmax = from.size();
    T err = 0;
    long long tim;

    T *orb__ = new T[nsite * nsite];
    colmaj<T> orb_(orb__, nsite);
    orbital_mat<T> orb(BLIS_UPPER, nsite, orb_);
    orb.randomize(0.6, 1234);

    T *mat_ = new T[nelec * nelec];
    updated_tdi<T> tdi(orb, cfg, mat_, nelec, mmax);

    for (dim_t i = 0; i < mmax; ++i) {
        if (from.at(i) < 0)
            tdi.pop_update(true);
        else
            tdi.push_update_safe(to.at(i), from[i], true);

        std::vector<dim_t> cfg2(tdi.elem_cfg);
        // Apply hopping.
        for (dim_t j = 0; j < tdi.from_idx.size(); ++j)
          cfg2.at(tdi.from_idx.at(j)) = tdi.to_site.at(j);

        T *mat2_ = new T[nelec * nelec];
        T pfa2 = updated_tdi<T>(orb, cfg2, mat2_, nelec, 1).Pfa;

        double err_ = std::abs(pfa2 - tdi.Pfa * tdi.PfaRatio);
        printf("Pfa = %e; Pfa diff = %e\n", tdi.Pfa, err_);
        err += err_;

        delete[] mat2_;
    }

    delete[] mat_;
    delete[] orb__;

    return info_t{ err, tim };
}

int main(void)
{
    test_tdi_update<double>(200,
        std::vector<dim_t>{ 8, 63, 20, 53, 57, 14, 7, 22, 32, 67, 56, 66, 46, 16, 18, 26, 70, 28, 51, 64, 17, 68, 49, 4, 48, 41, 43, 21, 19, 45, 58, 11, 1, 24, 2, 69, 12, 23, 30, 10, 65, 3, 31, 55, 52, 37, 42, 34, 25, 60, 59, 50, 40, 9, 44, 54, 6, 36, 13, 29, 27, 38, 61, 33, 47, 5, 15, 39, 62, 35 },
        std::vector<dim_t>{  0,  39,  49,  10,  60,   3 },
        std::vector<dim_t>{ 99, 101, 120, 199, 133, 177 });

    return 0;
}

