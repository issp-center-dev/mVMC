/**
 * \copyright Copyright (c) Dept. Phys., Univ. Tokyo
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "config_manager.tcc"
#include "pfaffian.tcc"
#include "invert.tcc"
#include "skr2k.tcc"
#include "skslc.tcc"
#include "skmv.tcc"
#include <iostream>
#include <vector>

namespace vmc
{
namespace orbital
{
template <typename orbital_t_,
          typename amp_t_ = typename orbital_t_::T>
struct updated_tdi {
  using orbital_t = orbital_t_;
  using index_t = vmc::config_manager::index_t;
  using amp_t = amp_t_;
  using T = typename orbital_t::T;
  using invmat_t = typename orbital_t::matrix_t; ///< Just for mVMC's C wrapper. Mathematically not so good.
  using matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using submat_t = Eigen::Block<matrix_t, Eigen::Dynamic, Eigen::Dynamic, true>;
  using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  using idxvec_t = Eigen::VectorXi;

  orbital_t &Xij;
  const index_t nelec; ///< Total Fermionic number / Matrix size.
  const index_t mmax;  ///< Max number of updates contained.
  const bool single_hop_alpha; ///< Whether to simplify k=1 situations.

private:
  char uplo; ///< Can be switched only with update_uplo().
  index_t nq_updated; // Update of Q(:, i) can be delayed against U.

  // Matrix M.
  invmat_t &M; ///< (nelec, nelec)

  // Updator blocks U, inv(M)U and inv(M)V.
  // Note that V is not stored.
  matrix_t U; ///< Serves also as B*inv
  matrix_t Q; ///< Contains the whole B buffer. Assert![ Q = M U ]
  // This member causes pointer problem on mixed precision.
  // (Which was fixed by constructing Block objects on-the-fly.)
  // TODO: Check its reason.
  // submat_t P; ///< &Q(0, mmax).

  // Central update buffer and its blocks.
  matrix_t W;
  matrix_t UMU, UMV, VMV;
  matrix_t Cp;       ///< C+BMB = [ W -I; I 0 ] + [ UMU UMV; -UMV VMV ]
  matrix_t Gc;       ///< Used to be Gaussian vectors. Now just scratchpads.
  idxvec_t cPiv;     ///< Pivot when tri-diagonalizing C.
  amp_t Pfa;
  amp_t PfaRatio; //< Updated Pfaffian / Base Pfaffian.
  // TODO: Maybe one should log partially the computed Pfaffian history.

  // Configuration.
  vmc::config_manager elem_cfg;

public:

  void initialize() {
    using namespace std;
    elem_cfg.merge_config();
    elem_cfg.reserve_space(mmax);
    const auto &cfg = elem_cfg.config_base();
    matrix_t G(nelec, nelec);
    nq_updated = 0;

    switch (uplo) {
    case 'U':
    case 'u':
      for (index_t j = 0; j < nelec; ++j) {
        for (index_t i = 0; i < j; ++i) {
          M(i, j) = Xij(cfg.at(i), cfg.at(j));
        }
        M(j, j) = T(0.0);
      }
      break;

    default:
      cerr << "updated_tdi: Lower triangular is not implemented yet."
           << endl;
    }

    // Allocate scratchpad.
    vector_t vT(nelec - 1);
    idxvec_t iPiv(nelec);
    l2e::le_mat_t<T> M_(M);

    // Perform LTL factorization regardless of `uplo`.
    signed info = sktrf(uplo, 'N', M_, &iPiv(0), &G(0, 0), G.size());
    Pfa = ltl2pfa<matrix_t, idxvec_t, amp_t>(uplo, M, iPiv);
    PfaRatio = 1.0; // Reset Pfaffian ratio.
    ltl2inv(uplo, M, iPiv, vT, G);
  }

  updated_tdi(orbital_t &Xij_, invmat_t &M_, index_t nelec_, index_t mmax_, std::vector<index_t> cfg)
      : Xij(Xij_), nelec(nelec_), mmax(mmax_), nq_updated(0), single_hop_alpha(true),
        M(M_), U(nelec, mmax * 2), Q(nelec, mmax * 2),
        W(mmax, mmax), UMU(mmax, mmax), UMV(mmax, mmax), VMV(mmax, mmax),
        // Cp, Gc and cPiv are dynamically allocated.
        Pfa(0.0), PfaRatio(1.0), elem_cfg(Xij_.norb()/2, cfg), uplo('U')
  { initialize(); }

  // Unsafe construction without initializaion.
  updated_tdi(orbital_t &Xij_, invmat_t &M_, index_t nelec_, index_t mmax_)
      : Xij(Xij_), nelec(nelec_), mmax(mmax_), nq_updated(0), single_hop_alpha(true),
        M(M_), U(nelec, mmax * 2), Q(nelec, mmax * 2),
        W(mmax, mmax), UMU(mmax, mmax), UMV(mmax, mmax), VMV(mmax, mmax),
        Pfa(0.0), PfaRatio(1.0), elem_cfg(Xij_.norb()/2), uplo('U') { }

  // Copy construction could be unsafe since it does not always contain initializaion.
  // M is automatically allocated in this case.
  updated_tdi(const updated_tdi<orbital_t> &tdi_) = delete;
    /*
      : Xij(tdi_.Xij), nelec(tdi_.nelec), mmax(tdi_.mmax), nq_updated(0), single_hop_alpha(true),
        M(new invmat_t(...)), U(nelec, mmax * 2), Q(nelec, mmax * 2),
        W(mmax, mmax), UMU(mmax, mmax), UMV(mmax, mmax), VMV(mmax, mmax),
        Pfa(0.0), PfaRatio(1.0), elem_cfg(tdi_.Xij.norb()/2, tdi_.get_config()), uplo('U') {
    if (elem_cfg.config().size() == nelec &&
        elem_cfg.config(0) != elem_cfg.config(1) && // Short-circuit complicated check.
        is_perm(elem_cfg.config()))
      initialize();
  } */

  bool is_perm(const vmc::config_manager::base_t &cfg) const
  {
    std::vector<bool> occupied(Xij.norb(), false);
    for (const auto &xi : cfg) {
      if (xi > Xij.norb())
        return false;
      if (occupied.at(xi))
        return false;
      occupied.at(xi) = true;
    }
    return true;
  }

  amp_t get_amplitude() { return Pfa * get_amplitude_ratio(); }
  amp_t get_amplitude_ratio() {
    if (PfaRatio == 0.0)
      update_pfaffian_ratio();
    return PfaRatio;
  }
  int max_updates() const { return mmax; }
  int num_updates() const { return elem_cfg.from_idx().size(); }

  const vmc::config_manager::base_t get_config() const
  { return elem_cfg.config(); }

  const vmc::config_manager::base_t &get_config_base() const
  { return elem_cfg.config_base(); }

  const vmc::config_manager &get_config_manager() const
  { return elem_cfg; }

  void attach_config(const vmc::config_manager::base_t &cfg)
  {
    assert_(cfg.size() == nelec, typeid(*this).name(), "Invalid config feeded.");
    elem_cfg.attach_config(cfg);
  }

  void assemble_C_BMB() {
    using namespace std;
    using namespace Eigen;
    index_t k = elem_cfg.from_idx().size();
    if (!k)
      return;

    // Assemble (half of) C+BMB buffer.
    Cp = matrix_t(2 * k, 2 * k);
    switch (uplo) {
    case 'U':
    case 'u':
      // Left-lower block skipped due to uplo='U'.
      Cp << W(seq(0, k-1), seq(0, k-1)) + UMU(seq(0, k-1), seq(0, k-1)), UMV(seq(0, k-1), seq(0, k-1)) - matrix_t::Identity(k, k),
            matrix_t::Zero(k, k),                                        VMV(seq(0, k-1), seq(0, k-1));
      break;

    default:
      cerr << "updated_tdi<T>::assemble_C_BMB:"
           << " Only upper-triangular storage is supported." << endl;
    }
  }

  // Update Q rows buffer.
  bool require_Q(bool require_all) {
    using namespace Eigen;
    using namespace l2e;

    const index_t k = elem_cfg.from_idx().size();
    const index_t n = nelec;
    index_t k_cal;
    // Require all K columns instead of first K-1.
    if (require_all)
      k_cal = k;
    else
      k_cal = k - 1;

    if (nq_updated >= k_cal)
      // Do nothing if nothing is to be updated.
      return false;

    for (; nq_updated < k_cal; ++nq_updated) {
      // Update single column. Use SKMV.
      Q(all, nq_updated) = skmv(uplo, T(1.0), M, U(all, nq_updated));
    } /* else {
      // Update multiple columns. Use SKMM.
      skmm('L', uplo, 'N', 'N',
           T(1.0), le_mat_t<T>(M), le_mat_t<T>U(all, seq(nq_updated, k_cal-1)),
           T(0.0), le_mat_t<T>(Q(all, seq(nq_updated, k_cal-1))));
    }
    nq_updated = k_cal;
    */
    return true;
  }

  // UMU needs to be recalculated for hopping change.
  // TODO: Find some way to avoid recalculating UQ?
  bool require_UMU() {
    using namespace Eigen;

    const index_t k = elem_cfg.from_idx().size();
    const index_t n = nelec;

    for (index_t o = 0; o < k; ++o)
      for (index_t l = 0; l < o; ++l)
        switch (uplo) {
        case 'U':
        case 'u':
          UMU(l, o) = -(U(all, o).transpose() * Q(all, l))(0, 0);
          break;
        case 'L':
        case 'l':
          // Always assume l < o to calculate one less column of Q.
          UMU(o, l) = (U(all, o).transpose() * Q(all, l))(0, 0);
          break;
        default:
          break;
        }
    return true;
  }

  // Update osi <- os[msj].
  // i.e. c+_i c_{x_j}.
  void push_update(index_t osi, index_t msj, bool compute_pfa) {
    using namespace std;
    using namespace Eigen;

    // This is the k-th hopping.
    index_t n = nelec;
    index_t k = elem_cfg.from_idx().size();
    index_t osj = elem_cfg.config_base(msj);

    elem_cfg.push_update(osi, msj);

    // TODO: Check for hopping-backs.
    // This can only be handled by cancellation.
    // Singularity will emerge otherwise.

    // TODO: Speed-up lookup of Xij.
    for (index_t i = 0; i < n; ++i)
      U(i, k) = Xij(elem_cfg.config_base(i), osi) - Xij(elem_cfg.config_base(i), osj);
    // for (index_t i = 0; i < Xij.norb(); ++i)
    //   if (elem_cfg.inv(i) >= 0)
    //     U(elem_cfg.inv(i), k) = Xij(i, osi) - Xij(i, osj);

    // P matrix from latter half of Q (or rigorously speaking, B matrix).
    submat_t P(Q(Eigen::all, Eigen::seq(mmax, mmax * 2 - 1)));

    P(all, k) = skslc(uplo, msj, M);

    // Updated the already logged U.
    for (index_t l = 0; l < k; ++l) {
      index_t msl = elem_cfg.from_idx().at(l);
      T U_jl_ = U(msj, l); ///< Backup this value before change.

      // Write updates.
      U(msl, k) = Xij(elem_cfg.to_orb().at(l), osi) - Xij(elem_cfg.config_base(msl), osj);
      U(msj, l) = Xij(osi, elem_cfg.to_orb().at(l)) - Xij(osj, elem_cfg.config_base(msl));

      // Change of Q[:, l] induced by change of U[msj, l];
      // Utilize that P[:, k] = M[:, msj] ///< M is skew-symmetric.
      Q(all, l) += P(all, k) * (U(msj, l) - U_jl_);

      // Change of UMV[l, k] induced by change of U[msj, l].
      for (index_t o = 0; o < k; ++o)
        UMV(l, o) += (U(msj, l) - U_jl_) * P(msj, o);
    }

    // Update UMV and VMV. UMU is updated elsewhere.
    switch (uplo) {
    case 'U':
    case 'u':
      for (index_t l = 0; l < k; ++l)
        W(l, k) = -Xij(elem_cfg.to_orb(l), osi) + Xij(elem_cfg.config_base(elem_cfg.from_idx(l)), osj);

      for (index_t l = 0; l < k; ++l) {
        UMV(l, k) = (U(all, l).transpose() * P(all, k))(0, 0);
        UMV(k, l) = (U(all, k).transpose() * P(all, l))(0, 0);
      }
      UMV(k, k) = (U(all, k).transpose() * P(all, k))(0, 0);
      for (index_t l = 0; l < k; ++l)
        VMV(l, k) /* -VMV(k, l) */ = -P(msj, l); ///< P(from_idx.at(l), k);
      break;

    default:
      cerr << "updated_tdi::push_update:"
           << " Only upper-triangular storage is supported." << endl;
    }

    if (compute_pfa) {
      // NOTE: Update k to be new size?
      // k += 1;

      update_pfaffian_ratio();
    } else
      // Set to 0.0 to denote dirty.
      PfaRatio = 0.0;
  }

  void update_pfaffian_ratio(void)
  {
    index_t k = elem_cfg.from_idx().size();
    submat_t P(Q(Eigen::all, Eigen::seq(mmax, mmax * 2 - 1)));
    {
      // All compute_pfa requires first K-1 of Q.
      require_Q(false);

      // Recompute UMU.
      require_UMU();
      // Reassemble C and scratchpads.
      assemble_C_BMB();

      // If it's the first update Pfafian can be directly read out.
      if (k == 1) {
        PfaRatio = -UMV(0, 0) + T(1.0);
        return;
      }

      // Compute pfaffian.
      l2e::le_mat_t<T> Cp_(Cp);
      Gc = matrix_t::Zero(2*k, 2*k);
      cPiv = idxvec_t::Zero(2*k);
      // Use not skpfa<T> for mixed precision.
      signed info = sktrf(uplo, 'P', Cp_, &cPiv(0), &Gc(0, 0), Gc.size());
      PfaRatio = ltl2pfa<matrix_t, idxvec_t, amp_t>(uplo, Cp, cPiv);

      // Pfaffian of C = [ W -I; I 0 ].
      PfaRatio *= pow(-1.0, k * (k + 1) / 2);
    }
  }

  bool check_update_safety(const std::vector<index_t> &ms, bool merge_error = false)
  {
    if (elem_cfg.from_idx().size() + ms.size() > mmax)
      return false;

    // Check if alrady hopped out.
    for (index_t l = 0; l < elem_cfg.from_idx().size(); ++l)
      for (index_t j = 0; j < ms.size(); ++j)
        if (ms.at(j) == elem_cfg.from_idx().at(l)) {
          assert_(!merge_error, typeid(*this).name(), "Conflict on non-mergable push.");
          return false;
        }
    return true;
  }

  void push_update(index_t osi, index_t msj) { push_update(osi, msj, true); }

  // Check duplicate and push.
  void push_update_safe(index_t osi, index_t msj, bool compute_pfa, bool merge_error = false) {
    if (!check_update_safety({ msj }))
      merge_updates();

    push_update(osi, msj, compute_pfa);
  }

  void push_update_safe(index_t osi, index_t msj) {
    push_update_safe(osi, msj, true);
  }

  void pop_update(bool compute_pfa) {
    using namespace std;
    using namespace Eigen;

    // Pop out.
    index_t msj_ = elem_cfg.from_idx().back();
    index_t osi_ = elem_cfg.to_orb().back();
    index_t osj_ = elem_cfg.config_base(msj_);
    elem_cfg.pop_update();

    // New update size.
    index_t k = elem_cfg.from_idx().size();
    submat_t P(Q(Eigen::all, Eigen::seq(mmax, mmax * 2 - 1)));

    // Revert U[:, 1:k-1] update.
    for (index_t l = 0; l < k; ++l) {
      index_t msl = elem_cfg.from_idx(l);
      T U_jl_ = U(msj_, l);

      // Revert U rows. (final column is popped out).
      U(msj_, l) = Xij(osj_, elem_cfg.to_orb(l)) - Xij(osj_, elem_cfg.config_base(msl));

      // Change of Q[:, l] induced by change of U[msj, l];
      T U_diff_msj = U(msj_, l) - U_jl_;
      switch (uplo) {
      case 'U':
      case 'u':
        Q(seq(0, msj_-1), l) += M(seq(0, msj_-1), msj_) * U_diff_msj;
        Q(seq(msj_+1, nelec-1), l) -= M(msj_, seq(msj_+1, nelec-1)) * U_diff_msj;
        break;

      default:
        cerr << "updated_tdi<T>::pop_update:"
             << " Only upper-triangular storage is supported." << endl;
      }

      // Change of UMV[l, k] induced by change of U[msj, l].
      for (index_t o = 0; o < k; ++o)
        UMV(l, o) += (U(msj_, l) - U_jl_) * P(msj_, o);
    }
    // Pop out outdated Q.
    if (nq_updated > k)
      nq_updated = k;

    // Compute new (previous, in fact) Pfaffian.
    if (k == 0)
      PfaRatio = 1.0;
    else if (compute_pfa) {
      update_pfaffian_ratio();

    } else
      // Set to 0.0 to denote dirty.
      PfaRatio = 0.0;
  }

  void merge_updates() {
    using namespace std;

    index_t n = nelec;
    index_t k = elem_cfg.from_idx().size();
    if (k == 0)
      return;

    // Allocate scratchpad.
    vector_t vT(2 * k - 1);

    // Update whole Q.
    require_Q(true);

    if (k == 1) {
      // M is dirty.
      if (PfaRatio == 0.0) {
        require_UMU();
        assemble_C_BMB();
        PfaRatio = -Cp(0, 1);
      }
      // Trivial inverse.
      Cp(0, 1) = T(-1.0) / Cp(0, 1);
      Cp(1, 0) = T(-1.0) / Cp(1, 0);
    } else {
      // Redo the tridiagonal factorization.
      require_UMU();
      assemble_C_BMB();

      l2e::le_mat_t<T> Cp_(Cp);
      Gc = matrix_t::Zero(2*k, 2*k);
      cPiv = idxvec_t::Zero(2*k);
      signed info = sktrf(uplo, 'N', Cp_, &cPiv(0), &Gc(0, 0), 2*k * 2*k);
      // PfaRatio is dirty.
      if (PfaRatio == 0.0)
        PfaRatio = ltl2pfa<matrix_t, idxvec_t, amp_t>(uplo, Cp, cPiv) *
          T(pow(-1.0, k * (k + 1) / 2));
      ltl2inv(uplo, Cp, cPiv, vT, Gc);
    }
    assert_(Cp.rows() == 2 * k, typeid(*this).name(), "Cp is bad.");
    inv_update(Cp);

    // Apply hopping.
    elem_cfg.merge_config();
    nq_updated = 0;
    Pfa *= PfaRatio;
    PfaRatio = 1.0;
  }

  void inv_update(matrix_t &C) {
    using namespace Eigen;
    const index_t k = C.rows() / 2;
    submat_t P(Q(Eigen::all, Eigen::seq(mmax, mmax * 2 - 1)));

    auto ABC = U(all, seq(0, 2 * k - 1)); ///< Use U as inv(A)*B*upper(C) buffer.
    auto AB  = Q(all, seq(0, 2 * k - 1));

    if (k == 1 && single_hop_alpha) {
      // k == 1 requires no copying
      skr2k(uplo, 'N', C(0, 1), Q(all, 0), P(all, 0), T(1.0), M);
      // See below for reason of calling this procedule.
      // update_uplo(uplo);
      return;
    }

    // Close empty space between Q and P.
    if (k != mmax)
      for (index_t j = 0; j < k; ++j)
        AB(all, k + j) = P(all, j);

    l2e::le_mat_t<T> M_(M);
    ABC = AB * C;
    gemmt<T>(uplo, 'N', 'T', T(1.0),
             l2e::le_mat_t<T>(ABC),
             l2e::le_mat_t<T>(AB),
             T(1.0), M_);

    // Identity update to complete antisymmetric matrix.
    // This called due to lack of skmm support at the moment.
    // update_uplo(uplo);
  }

  /**
   * Complete antisymmetric matrix.
   */
  static void skcomplete(const char &uplo_, matrix_t &A) {
    const index_t n = A.rows();

    for (index_t j = 0; j < n; ++j) {
      for (index_t i = 0; i < j; ++i) {
        switch (uplo_) {
        case 'U':
        case 'u':
          A(j, i) = -A(i, j);
          break;

        case 'L':
        case 'l':
          A(i, j) = -A(j, i);
          break;

        default:
          break;
        }
      }
      A(j, j) = T(0.0);
    }
  }

  void update_uplo(const char &uplo_new) {
    using namespace Eigen;

    skcomplete(uplo /* NOTE: Old uplo. */, M);
    const int k = elem_cfg.from_idx().size();
    if (k) {
      skcomplete(uplo, UMU(seq(0, k-1), seq(0, k-1)));
      skcomplete(uplo, VMV(seq(0, k-1), seq(0, k-1)));
      skcomplete(uplo,   W(seq(0, k-1), seq(0, k-1)));
    }

    uplo = uplo_new;
  }

  Eigen::Vector<amp_t, Eigen::Dynamic> batch_query_amplitudes(int N, idxvec_t to_orbs, idxvec_t from_ids) {
    using namespace Eigen;
    using ampvec_t = Eigen::Vector<amp_t, Eigen::Dynamic>;

    bool upper = (uplo == 'U' || uplo == 'u');
    int k0 = from_ids.size();
    int k1 = k0 / N;
    assert_(k0 % N == 0, typeid(*this).name(), "Misaligned batch query.");
    ampvec_t result(k1);

    // Works w/ merged M only.
    merge_updates();
    // Operates on complete M.
    skcomplete(uplo, M);

    matrix_t UbB;
    if (N > 1)
      UbB = matrix_t::Zero(nelec, k1 * (N-1));
    matrix_t UbF= matrix_t::Zero(nelec, k1);
    matrix_t Pb = matrix_t::Zero(nelec, k0);

    // Grouped comput. of Ub.
    for (index_t l1 = 0; l1 < k1; ++l1) {
      for (index_t l = l1 * N; l < (l1 + 1) * N; ++l)
        elem_cfg.push_update(to_orbs(l), from_ids(l));

      int osi, osj;

      for (index_t l2 = 0; l2 < N-1; ++l2) {
        index_t l = l1 * N + l2;
        index_t lB = l1 * (N-1) + l2;

        osi = to_orbs(l);
        osj = elem_cfg.config_base(from_ids(l));

        for (index_t i = 0; i < nelec; ++i) {
          index_t osx = i == from_ids(l) ? elem_cfg.config_base(i) : elem_cfg.config(i);
          UbB(i, lB) = Xij(osx, osi) - Xij(elem_cfg.config_base(i), osj);
        }
      }

      // Final column in the grouped U.
      {
        index_t l = l1 * N + N - 1;
        osi = to_orbs(l);
        osj = elem_cfg.config_base(from_ids(l));

        for (index_t i = 0; i < nelec; ++i) {
          index_t osx = i == from_ids(l) ? elem_cfg.config_base(i) : elem_cfg.config(i);
          UbF(i, l1) = Xij(osx, osi) - Xij(elem_cfg.config_base(i), osj);
        }
      }

      for (index_t ll = 0; ll < N; ++ll)
        elem_cfg.pop_update();
    }

    // Pb needs no grouping.
    for (index_t l = 0; l < k0; ++l)
      Pb(all, l) = M(all, from_ids(l));

    matrix_t QbB;
    if (N > 1)
      // Require Q.
      QbB = M * UbB;

    // Batch assemble C.
    matrix_t Cb = matrix_t::Zero(2*N, 2*N * k1);
    for (index_t l1 = 0; l1 < k1; ++l1) {
      auto CCur = Cb(all, seq(2*N * l1, 2*N * (l1+1) - 1));
      auto UCur = UbB(all, seq((N-1) * l1, (N-1) * (l1+1) - 1));
      auto QCur = QbB(all, seq((N-1) * l1, (N-1) * (l1+1) - 1));
      auto PCur = Pb(all, seq(N * l1, N * (l1+1) - 1));

      // UMU + W.
      for (index_t o = 0; o < N-1; ++o)
        for (index_t l = 0; l < o; ++l)
          if (upper)
            CCur(l, o) = -UCur(all, o).transpose() * QCur(all, l) -
              Xij(to_orbs(l1 * N + l), to_orbs(l1 * N + o)) +
              Xij(elem_cfg.config_base(from_ids(l1 * N + l)),
                  elem_cfg.config_base(from_ids(l1 * N + o)));
          else
            CCur(o, l) =  UCur(all, o).transpose() * QCur(all, l) -
              Xij(to_orbs(l1 * N + l), to_orbs(l1 * N + o)) +
              Xij(elem_cfg.config_base(from_ids(l1 * N + l)),
                  elem_cfg.config_base(from_ids(l1 * N + o)));
      {
        for (index_t l = 0; l < N-1; ++l)
          if (upper)
            CCur(l, N-1) = -UbF(all, l1).transpose() * QCur(all, l) -
              Xij(to_orbs(l1 * N + l), to_orbs(l1 * N + N-1)) +
              Xij(elem_cfg.config_base(from_ids(l1 * N + l)),
                  elem_cfg.config_base(from_ids(l1 * N + N-1)));
          else
            CCur(N-1, l) =  UbF(all, l1).transpose() * QCur(all, l) -
              Xij(to_orbs(l1 * N + l), to_orbs(l1 * N + N-1)) +
              Xij(elem_cfg.config_base(from_ids(l1 * N + l)),
                  elem_cfg.config_base(from_ids(l1 * N + N-1)));
      }

      // UMV
      if (upper) {
        if (N > 1)
          CCur(seq(0, N - 2), seq(N, 2*N - 1)) = UCur.transpose() * PCur;
        CCur(N - 1, seq(N, 2*N - 1)) = UbF(all, l1).transpose() * PCur;
      } else
        assert_(false, typeid(*this).name(), "Lower diag not implemented.");
      // -I
      CCur(seq(0, N - 1), seq(N, 2*N - 1)) -= matrix_t::Identity(N, N);

      // VMV
      for (index_t o = 0; o < N; ++o)
        for (index_t l = 0; l < o; ++l)
          if (upper)
            CCur(N + l, N + o) = -PCur(from_ids(l1 * N + o), l);
          else
            CCur(N + o, N + l) =  PCur(from_ids(l1 * N + o), l);
    }

    // Batch Pfaffian call.
    matrix_t Gc = matrix_t::Zero(2*N, 2*N);
    idxvec_t cPiv = idxvec_t::Zero(2*N);
    for (index_t l1 = 0; l1 < k1; ++l1) {
      auto CCur = Cb(all, seq(2*N * l1, 2*N * (l1+1) - 1));
      if (N == 1) {
        result(l1) = -CCur(0, 1) + T(1.0);
        continue;
      }
      l2e::le_mat_t<T> C_(CCur);
      signed info = sktrf(uplo, 'P', C_, &cPiv(0), &Gc(0, 0), Gc.size());
      result(l1) = ltl2pfa<matrix_t, idxvec_t, amp_t>(uplo, CCur, cPiv);
    }
    if (N > 1)
      result = result * pow(-1.0, N * (N + 1) / 2);

    return result;
  }

};
}
}
