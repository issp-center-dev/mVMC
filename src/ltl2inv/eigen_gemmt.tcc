/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "lapack2eigen.hh"
#include "ilaenv_lauum.hh"


template <typename T>
inline void gemmt(const char &uploC,
                  const char &transA,
                  const char &transB,
                  const T &alpha,
                  const l2e::le_mat_t<T> &A_,
                  const l2e::le_mat_t<T> &B_,
                  const T &beta,
                        l2e::le_mat_t<T> &C_)
{
    using LProxy = l2e::le_mat_t<T>;
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    const int npanel = ilaenv_lauum<T>(uploC, C_.rows);

    bool conjA = transA == 'C' || transA == 'c';
    bool conjB = transB == 'C' || transB == 'c';

    const auto &A = (transA == 'T' || transA == 't') ? A_.transpose().to_eigen() : A_.to_eigen();
    const auto &B = (transB == 'T' || transB == 't') ? B_.transpose().to_eigen() : B_.to_eigen();
    auto        C =                                                                C_.to_eigen();

    if (uploC == 'L' || uploC == 'l') {

        for (int ipanel = 0; ipanel < A.rows(); ipanel += npanel) {
            int sz_panel = std::min<int>(npanel, A.rows() - ipanel);

            // Diagonal.
            const auto &Apanel = A(Eigen::seq(ipanel, ipanel+sz_panel-1), Eigen::all);
            const auto &Bpanel = B(Eigen::all, Eigen::seq(ipanel, ipanel+sz_panel-1));
            const auto &Cpanel = C(Eigen::seq(ipanel, ipanel+sz_panel-1),
                                   Eigen::seq(ipanel, ipanel+sz_panel-1));
            Matrix Cupdated = Cpanel * beta + Apanel * Bpanel * alpha;
            LProxy Cpanel_(Cpanel);
            lacpy(uploC, LProxy(Cupdated), Cpanel_);

            // Off-diagonal.
            const auto &Aremain = A(Eigen::seq(ipanel+sz_panel, A.rows()-1), Eigen::all);
            const auto &Bremain = B(Eigen::all, Eigen::seq(ipanel+sz_panel, B.cols()-1));
            auto        Cremain = C(Eigen::seq(ipanel+sz_panel, C.rows()-1),
                                    Eigen::seq(ipanel, ipanel+sz_panel-1));
            Cremain = Cremain * beta + Aremain * Bpanel * alpha;
        }

    } else { // (uploC == 'U' || uploC == 'u')

        for (int ipanel = (A.rows()-1) / npanel * npanel; ipanel >= 0; ipanel -= npanel) {
            int sz_panel = std::min<int>(npanel, A.rows() - ipanel);

            // Diagonal.
            const auto &Apanel = A(Eigen::seq(ipanel, ipanel+sz_panel-1), Eigen::all);
            const auto &Bpanel = B(Eigen::all, Eigen::seq(ipanel, ipanel+sz_panel-1));
            const auto &Cpanel = C(Eigen::seq(ipanel, ipanel+sz_panel-1),
                                   Eigen::seq(ipanel, ipanel+sz_panel-1));
            Matrix Cupdated = Cpanel * beta + Apanel * Bpanel * alpha;
            LProxy Cpanel_(Cpanel);
            lacpy(uploC, LProxy(Cupdated), Cpanel_);

            // Off-diagonal.
            const auto &Aremain = A(Eigen::seq(0, ipanel-1), Eigen::all);
            const auto &Bremain = B(Eigen::all, Eigen::seq(0, ipanel-1));
            auto        Cremain = C(Eigen::seq(0, ipanel-1),
                                    Eigen::seq(ipanel, ipanel+sz_panel-1));
            Cremain = Cremain * beta + Aremain * Bpanel * alpha;
        }

    }
}
