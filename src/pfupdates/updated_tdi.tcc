/**
 * \copyright Copyright (c) Dept. Phys., Univ. Tokyo
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "orbital_mat.tcc"
#include "blalink.hh"
#include "skpfa.hh"
#include "sktdi.hh"
#include "skr2k.tcc"
#include "skslc.tcc"
#include "skmv.tcc"
#include "optpanel.hh"
#include <iostream>
#include <vector>

template <typename T> struct updated_tdi {
  orbital_mat<T> &Xij;
  uplo_t uplo; ///< Can be switched only with update_uplo().
  const dim_t nelec;
  const dim_t mmax;
  dim_t nq_updated; // Update of Q(:, i) can be delayed against U.

  // Performance tuning constants.
  const bool single_hop_alpha; ///< Whether to simplify k=1 situations.
  const dim_t npanel_big;
  const dim_t npanel_sub;

  // Matrix M.
  matrix_t<T> M;

  // Updator blocks U, inv(M)U and inv(M)V.
  // Note that V is not stored.
  matrix_t<T> U; ///< Serves also as B*inv
  matrix_t<T> Q; ///< Contains the whole B buffer. Assert![ Q = M U ]
  matrix_t<T> P; ///< Q(0, mmax)

  // Central update buffer and its blocks.
  matrix_t<T> W;
  matrix_t<T> UMU, UMV, VMV;
  matrix_t<T> Cp;       ///< C+BMB = [ W -I; I 0 ] + [ UMU UMV; -UMV VMV ]
  matrix_t<T> Gc;       ///< Gaussian vectors when tri-diagonalizing C.
  signed *const cPov; ///< Pivot when tri-diagonalizing C.
  T Pfa;
  T PfaRatio; //< Updated Pfaffian / Base Pfaffian.
  // TODO: Maybe one should log all Pfaffian histories.

  std::vector<dim_t> elem_cfg;
  std::vector<dim_t> from_idx;
  std::vector<dim_t> to_site;

  void initialize() {
    using namespace std;
    auto &cfg = elem_cfg;
    #ifdef UseBoost
    matrix_t<T> G(nelec, nelec);
    #else
    matrix_t<T> G(new T[nelec * nelec], nelec);
    #endif
    from_idx.clear();
    to_site.clear();
    from_idx.reserve(mmax * sizeof(dim_t));
    to_site.reserve(mmax * sizeof(dim_t));

    switch (uplo) {
    case BLIS_UPPER:
      for (dim_t j = 0; j < nelec; ++j) {
        for (dim_t i = 0; i < j; ++i) {
          M(i, j) = Xij(cfg.at(i), cfg.at(j));
        }
        M(j, j) = T(0.0);
      }
      break;

    default:
      cerr << "updated_tdi<T>: BLIS_LOWER is not supported. The error is fetal."
           << endl;
    }

    // Allocate scratchpad.
    signed *iPivFull = new signed[nelec + 1];
    dim_t lwork = nelec * npanel_big;
    T *pfwork = new T[lwork];

    #ifdef UseBoost
    signed info = skpfa(uplo, nelec, &M(0, 0), M.size1(), &G(0, 0), G.size1(), iPivFull,
                        true, &Pfa, pfwork, lwork);
    #else
    signed info = skpfa(uplo, nelec, &M(0, 0), M.ld, &G(0, 0), G.ld, iPivFull,
                        true, &Pfa, pfwork, lwork);
    #endif
#ifdef _DEBUG
    cerr << "SKPFA+INV: n=" << nelec << " info=" << info << endl;
#endif

    delete[] iPivFull;
    delete[] pfwork;
    #ifndef UseBoost
    delete[](&G(0, 0));
    #endif
  }

  ~updated_tdi() {
    #ifndef UseBoost
    delete[](&U(0, 0));
    delete[](&Q(0, 0));

    delete[](&Cp(0, 0));
    delete[](&Gc(0, 0));

    delete[](&W(0, 0));
    delete[](&UMU(0, 0));
    delete[](&UMV(0, 0));
    delete[](&VMV(0, 0));
    #endif

    delete[] cPov;
  }

  updated_tdi(orbital_mat<T> &Xij_, std::vector<dim_t> &cfg, T *M_, inc_t ldM,
              dim_t mmax_)
      : Xij(Xij_), nelec(cfg.size()), mmax(mmax_), nq_updated(0),
        single_hop_alpha(true), npanel_big(optpanel(nelec, 4)), npanel_sub(4),
    #ifdef UseBoost
        M(nelec, nelec),
        U(nelec, mmax * 2), Q(nelec, mmax * 2),
        P(nelec, mmax), W(mmax, mmax),
        UMU(mmax, mmax), UMV(mmax, mmax),
        VMV(mmax, mmax), Cp(2 * mmax, 2 * mmax),
        Gc(2 * mmax, 2 * mmax),
    #else
        M(M_, ldM),
        U(new T[nelec * mmax * 2], nelec), Q(new T[nelec * mmax * 2], nelec),
        P(&Q(0, mmax), Q.ld), W(new T[mmax * mmax], mmax),
        UMU(new T[mmax * mmax], mmax), UMV(new T[mmax * mmax], mmax),
        VMV(new T[mmax * mmax], mmax), Cp(new T[2 * mmax * 2 * mmax], 2 * mmax),
        Gc(new T[2 * mmax * 2 * mmax], 2 * mmax),
    #endif
        cPov(new signed[2 * mmax + 1]), Pfa(0.0), PfaRatio(1.0), elem_cfg(cfg),
        from_idx(0), to_site(0), uplo(BLIS_UPPER) {
    #ifdef UseBoost
    colmaj<T> M_tmp(M_, ldM);
    for (dim_t j = 0; j < nelec; ++j)
      for (dim_t i = 0; i < nelec; ++i)
        M(i, j) = M_tmp(i, j);
    #endif
    initialize();
  }

  /*
  updated_tdi(orbital_mat<T> &Xij_, std::vector<dim_t> &cfg, matrix_t<T> &M_,
              dim_t mmax_)
      : Xij(Xij_), nelec(cfg.size()), mmax(mmax_), nq_updated(0),
        single_hop_alpha(true), npanel_big(optpanel(nelec, 4)), npanel_sub(4),
        M(M_),
        U(new T[nelec * mmax * 2], nelec), Q(new T[nelec * mmax * 2], nelec),
        P(&Q(0, mmax), Q.ld), W(new T[mmax * mmax], mmax),
        UMU(new T[mmax * mmax], mmax), UMV(new T[mmax * mmax], mmax),
        VMV(new T[mmax * mmax], mmax), Cp(new T[2 * mmax * 2 * mmax], 2 * mmax),
        Gc(new T[2 * mmax * 2 * mmax], 2 * mmax),
        cPov(new signed[2 * mmax + 1]), Pfa(0.0), PfaRatio(1.0), elem_cfg(cfg),
        from_idx(0), to_site(0), uplo(BLIS_UPPER) {
    initialize();
  }*/

  // Unsafe construction without initializaion.
  updated_tdi(orbital_mat<T> &Xij_, dim_t nelec_, T *M_, inc_t ldM, dim_t mmax_) 
      : Xij(Xij_), nelec(nelec_), mmax(mmax_), nq_updated(0),
        single_hop_alpha(true), npanel_big(optpanel(nelec, 4)), npanel_sub(4),
    #ifdef UseBoost
        M(nelec, nelec),
        U(nelec, mmax * 2), Q(nelec, mmax * 2),
        P(nelec, mmax), W(mmax, mmax),
        UMU(mmax, mmax), UMV(mmax, mmax),
        VMV(mmax, mmax), Cp(2 * mmax, 2 * mmax),
        Gc(2 * mmax, 2 * mmax),
    #else
        M(M_, ldM),
        U(new T[nelec * mmax * 2], nelec), Q(new T[nelec * mmax * 2], nelec),
        P(&Q(0, mmax), Q.ld), W(new T[mmax * mmax], mmax),
        UMU(new T[mmax * mmax], mmax), UMV(new T[mmax * mmax], mmax),
        VMV(new T[mmax * mmax], mmax), Cp(new T[2 * mmax * 2 * mmax], 2 * mmax),
        Gc(new T[2 * mmax * 2 * mmax], 2 * mmax),
    #endif
        cPov(new signed[2 * mmax + 1]), Pfa(0.0), PfaRatio(1.0), elem_cfg(nelec, 0),
        from_idx(0), to_site(0), uplo(BLIS_UPPER) { 
    #ifdef UseBoost
    colmaj<T> M_tmp(M_, ldM);
    for (dim_t j = 0; j < nelec; ++j)
      for (dim_t i = 0; i < nelec; ++i)
        M(i, j) = M_tmp(i, j);
    #endif
  }

  T get_Pfa() { return Pfa * PfaRatio; }

  void assemble_C_BMB() {
    using namespace std;
    dim_t k = from_idx.size();

    // Assemble (half of) C+BMB buffer.
    switch (uplo) {
    case BLIS_UPPER:
      for (dim_t j = 0; j < k; ++j)
        for (dim_t i = 0; i < j; ++i)
          Cp(i, j) = W(i, j) + UMU(i, j);
      for (dim_t j = 0; j < k; ++j) {
        for (dim_t i = 0; i < k; ++i)
          Cp(i, j + k) = +UMV(i, j);

        Cp(j, j + k) -= T(1.0);
      }
      for (dim_t j = 0; j < k; ++j)
        for (dim_t i = 0; i < j; ++i)
          Cp(i + k, j + k) = VMV(i, j);
      break;

    default:
      cerr << "updated_tdi<T>::assemble_C_BMB:"
           << " Only upper-triangular storage is supported." << endl;
    }
  }

  // Update Q rows buffer.
  bool require_Q(bool all) {
    const dim_t k = from_idx.size();
    const dim_t n = nelec;
    dim_t k_cal;
    // Require all K columns instead of first K-1.
    if (all)
      k_cal = k;
    else
      k_cal = k - 1;

    if (nq_updated >= k_cal)
      // Do nothing if nothing is to be updated.
      return false;

    for (; nq_updated < k_cal; ++nq_updated) {
      // Update single column. Use SKMV.
      #ifdef UseBoost
      skmv(uplo, n, T(1.0), &M(0, 0), M.size1(), &U(0, nq_updated), &Q(0, nq_updated));
      #else
      skmv(uplo, n, T(1.0), &M(0, 0), M.ld, &U(0, nq_updated), &Q(0, nq_updated));
      #endif
    } /* else {
      // Update multiple columns. Use SKMM.
      #ifdef UseBoost
      skmm(BLIS_LEFT, uplo, BLIS_NO_CONJUGATE, BLIS_NO_TRANSPOSE, n, k_cal - nq_updated,
           T(1.0), &M(0, 0), M.size1(), &U(0, nq_updated), U.size1(),
           T(0.0), &Q(0, nq_updated), Q.size1());
      #else
      skmm(BLIS_LEFT, uplo, BLIS_NO_CONJUGATE, BLIS_NO_TRANSPOSE, n, k_cal - nq_updated,
           T(1.0), &M(0, 0), M.ld, &U(0, nq_updated), U.ld,
           T(0.0), &Q(0, nq_updated), Q.ld);
      #endif
    }
    nq_updated = k_cal;
    */
    return true;
  }

  // UMU needs to be recalculated for hopping change.
  // TODO: Find some way to avoid recalculating UQ?
  bool require_UMU() {
    const dim_t k = from_idx.size();
    const dim_t n = nelec;

    for (dim_t o = 0; o < k; ++o)
      for (dim_t l = 0; l < o; ++l)
        switch (uplo) {
        case BLIS_UPPER:
          UMU(l, o) = -dot(n, &U(0, o), 1, &Q(0, l), 1);
          break;
        case BLIS_LOWER:
          // Always assume l < o to calculate one less column of Q.
          UMU(o, l) = dot(n, &U(0, o), 1, &Q(0, l), 1);
          break;
        }
    return true;
  }

  // Inverts a configuration from fermion index to site index.
  void config_to_site(std::vector<dim_t> &cfg, std::vector<signed> &site, bool init) {
    using namespace std;

    if (site.size() != Xij.nsite || cfg.size() != nelec) {
      cerr << "updated_tdi<T>::inv_config:"
           << " Invalid config / invert config buffer." << endl;
      return;
    }

    if (init)
      for (dim_t i = 0; i < site.size(); ++i)
        site.at(i) = -1;
    for (dim_t i = 0; i < cfg.size(); ++i)
      site.at(cfg.at(i)) = i;

    return;
  }

  // Update osi <- os[msj].
  // i.e. c+_i c_{x_j}.
  void push_update(dim_t osi, dim_t msj, bool compute_pfa) {
    using namespace std;

    // This is the k-th hopping.
    dim_t n = nelec;
    dim_t k = from_idx.size();
    dim_t osj = elem_cfg.at(msj);

    // TODO: Check bounds, duplicates, etc.
    from_idx.push_back(msj);
    to_site.push_back(osi);

    // TODO: Check for hopping-backs. 
    // This can only be handled by cancellation.
    // Singularity will emerge otherwise.

    // Speed-up lookup of Xij.
    std::vector<signed> sitecfg(Xij.nsite, -1);
    config_to_site(elem_cfg, sitecfg, false);
    for (dim_t i = 0; i < Xij.nsite; ++i)
      // U(i, k) = Xij(elem_cfg.at(i), osi) - Xij(elem_cfg.at(i), osj);
      if (sitecfg.at(i) >= 0)
        // Direct access: special case when installed into VMC.
        U(sitecfg.at(i), k) = Xij.X(i, osi) - Xij.X(i, osj);

    skslc(uplo, n, msj, &P(0, k), &M(0, 0), M.ld);

    // Updated the already logged U.
    for (dim_t l = 0; l < k; ++l) {
      dim_t msl = from_idx.at(l);
      T U_jl_ = U(msj, l); ///< Backup this value before change.

      // Write updates.
      U(msl, k) = Xij(to_site.at(l), osi) - Xij(elem_cfg.at(msl), osj);
      U(msj, l) = Xij(osi, to_site.at(l)) - Xij(osj, elem_cfg.at(msl));

      // Change of Q[:, l] induced by change of U[msj, l];
      // Utilize that P[:, k] = M[:, msj] ///< M is skew-symmetric.
      axpy(n, U(msj, l) - U_jl_, &P(0, k), 1, &Q(0, l), 1);

      // Change of UMV[l, k] induced by change of U[msj, l].
      for (dim_t o = 0; o < k; ++o)
        UMV(l, o) += (U(msj, l) - U_jl_) * P(msj, o);
    }

    // Update UMV and VMV. UMU is updated elsewhere.
    switch (uplo) {
    case BLIS_UPPER:
      for (dim_t l = 0; l < k; ++l)
        W(l, k) = -Xij(to_site.at(l), osi) + Xij(elem_cfg.at(from_idx.at(l)), osj);

      for (dim_t l = 0; l < k; ++l) {
        UMV(l, k) = dot(n, &U(0, l), 1, &P(0, k), 1);
        UMV(k, l) = dot(n, &U(0, k), 1, &P(0, l), 1);
      }
      UMV(k, k) = dot(n, &U(0, k), 1, &P(0, k), 1);
      for (dim_t l = 0; l < k; ++l)
        VMV(l, k) /* -VMV(k, l) */ = -P(msj, l); ///< P(from_idx.at(l), k);
      break;

    default:
      cerr << "updated_tdi<T>::push_update:"
           << " Only upper-triangular storage is supported." << endl;
    }

    if (compute_pfa) {
      // NOTE: Update k to be new size.
      k += 1;
      // Allocate scratchpad.
      dim_t lwork = 2 * k * npanel_sub;
      T *pfwork = new T[lwork];

      // Calculate unupdated columns of U.
      require_Q(false);

      require_UMU();
      // Assemble (half of) C+BMB buffer.
      assemble_C_BMB();

      // If it's the first update Pfafian can be directly read out.
      if (k == 1) {
        PfaRatio = -UMV(0, 0) + T(1.0);
        delete[] pfwork;
        return;
      }

      // Compute pfaffian.
      #ifdef UseBoost
      signed info = skpfa(uplo, 2 * k, &Cp(0, 0), Cp.size1(), &Gc(0, 0), Gc.size1(), cPov,
                          false, &PfaRatio, pfwork, lwork);
      #else
      signed info = skpfa(uplo, 2 * k, &Cp(0, 0), Cp.ld, &Gc(0, 0), Gc.ld, cPov,
                          false, &PfaRatio, pfwork, lwork);
      #endif
#ifdef _DEBUG
      cerr << "SKPFA: info=" << info << endl;
#endif
      // Pfaffian of C = [ W -I; I 0 ].
      PfaRatio *= pow(-1.0, k * (k + 1) / 2);

      delete[] pfwork;
    } else
      // Set to 0.0 to denote dirty.
      PfaRatio = 0.0;
  }

  void push_update(dim_t osi, dim_t msj) { push_update(osi, msj, true); }

  // Check duplicate and push.
  void push_update_safe(dim_t osi, dim_t msj, bool compute_pfa) {
    if (from_idx.size() >= mmax)
      merge_updates();

    // Check if already hopped out.
    for (dim_t l = 0; l < from_idx.size(); ++l)
      if (msj == from_idx.at(l)) {
        merge_updates(); 
        break;
      }
    push_update(osi, msj, compute_pfa);
  }

  void push_update_safe(dim_t osi, dim_t msj) {
    push_update_safe(osi, msj, true);
  }

  void pop_update(bool compute_pfa) {
    using namespace std;

    // Pop out.
    dim_t msj_ = from_idx.at(from_idx.size() - 1);
    dim_t osi_ = to_site.at(to_site.size() - 1);
    dim_t osj_ = elem_cfg.at(msj_);
    from_idx.pop_back();
    to_site.pop_back();

    // New update size.
    dim_t k = from_idx.size();

    // Revert U[:, 1:k-1] update.
    for (dim_t l = 0; l < k; ++l) {
      dim_t msl = from_idx.at(l);
      T U_jl_ = U(msj_, l);

      // Revert U rows. (final column is popped out).
      U(msj_, l) = Xij(osj_, to_site.at(l)) - Xij(osj_, elem_cfg.at(msl));

      // Change of Q[:, l] induced by change of U[msj, l];
      T U_diff_msj = U(msj_, l) - U_jl_;
      switch (uplo) {
      case BLIS_UPPER:
        for (dim_t i = 0; i < msj_; ++i)
          Q(i, l) += M(i, msj_) * U_diff_msj;
        for (dim_t i = msj_ + 1; i < nelec; ++i)
          Q(i, l) -= M(msj_, i) * U_diff_msj;
        break;

      default:
        cerr << "updated_tdi<T>::pop_update:"
             << " Only upper-triangular storage is supported." << endl;
      }

      // Change of UMV[l, k] induced by change of U[msj, l].
      for (dim_t o = 0; o < k; ++o)
        UMV(l, o) += (U(msj_, l) - U_jl_) * P(msj_, o);
    }
    // Pop out outdated Q.
    if (nq_updated > k)
      nq_updated = k;

    // Compute new (previous, in fact) Pfaffian.
    if (compute_pfa) {
      // All compute_pfa requires first K-1 of Q.
      require_Q(false);

      // Recompute UMU.
      require_UMU();
      // Reassemble C and scratchpads.
      assemble_C_BMB();
      dim_t lwork = 2 * k * npanel_sub;
      T *pfwork = new T[lwork];

      #ifdef UseBoost
      signed info = skpfa(uplo, 2 * k, &Cp(0, 0), Cp.size1(), &Gc(0, 0), Gc.size1(), cPov,
                          false, &PfaRatio, pfwork, lwork);
      #else
      signed info = skpfa(uplo, 2 * k, &Cp(0, 0), Cp.ld, &Gc(0, 0), Gc.ld, cPov,
                          false, &PfaRatio, pfwork, lwork);
      #endif
      PfaRatio *= pow(-1.0, k * (k + 1) / 2);

      delete[] pfwork;
    } else
      // Set to 0.0 to denote dirty.
      PfaRatio = 0.0;
  }

  void merge_updates() {
    using namespace std;

    dim_t n = nelec;
    dim_t k = from_idx.size();
    if (k == 0)
      return;

    // Allocate scratchpad.
    dim_t lwork = 2 * k * npanel_sub;
    T *pfwork = new T[lwork];
    
    // Update whole Q.
    require_Q(true);

    // If M is dirty, redo the tridiagonal factorization.
    if (PfaRatio == T(0.0)) {
      require_UMU();
      assemble_C_BMB();
      #ifdef UseBoost
      signed info = skpfa(uplo, 2 * k, &Cp(0, 0), Cp.size1(), &Gc(0, 0), Gc.size1(), cPov,
                          false, &PfaRatio, pfwork, lwork);
      #else
      signed info = skpfa(uplo, 2 * k, &Cp(0, 0), Cp.ld, &Gc(0, 0), Gc.ld, cPov,
                          false, &PfaRatio, pfwork, lwork);
      #endif
      PfaRatio *= pow(-1.0, k * (k + 1) / 2);
    }

    if (k == 1) {
      // Trivial inverse.
      Cp(0, 1) = T(-1.0) / Cp(0, 1);
      Cp(1, 0) = T(-1.0) / Cp(1, 0);
    } else
      #ifdef UseBoost
      signed info = sktdi(uplo, 2 * k, &Cp(0, 0), Cp.size1(), &Gc(0, 0), Gc.size1(), cPov,
                          pfwork, lwork);
      #else
      signed info = sktdi(uplo, 2 * k, &Cp(0, 0), Cp.ld, &Gc(0, 0), Gc.ld, cPov,
                          pfwork, lwork);
      #endif
    inv_update(k, Cp);

    // Apply hopping.
    for (int j = 0; j < k; ++j)
      elem_cfg.at(from_idx.at(j)) = to_site.at(j);
    from_idx.clear();
    to_site.clear();
    nq_updated = 0;
    Pfa *= PfaRatio;
    PfaRatio = 1.0;

    delete[] pfwork;
  }

  void inv_update(dim_t k, matrix_t<T> &C) {
    #ifdef UseBoost
    matrix_t<T> &ABC = U; ///< Use U as inv(A)*B*upper(C) buffer.
    matrix_t<T> &AB = Q;
    #else
    matrix_t<T> ABC(&U(0, 0), U.ld); ///< Use U as inv(A)*B*upper(C) buffer.
    matrix_t<T> AB(&Q(0, 0), Q.ld);
    #endif

    if (k == 1 && single_hop_alpha) {
      // k == 1 requires no copying
      #ifdef UseBoost
      skr2k(uplo, BLIS_NO_TRANSPOSE, nelec, 1, C(0, 1), &Q(0, 0), Q.size1(),
            &P(0, 0), P.size1(), T(1.0), &M(0, 0), M.size1());
      #else
      skr2k(uplo, BLIS_NO_TRANSPOSE, nelec, 1, C(0, 1), &Q(0, 0), Q.ld,
            &P(0, 0), P.ld, T(1.0), &M(0, 0), M.ld);
      #endif
      // See below for reason of calling this procedule.
      // update_uplo(uplo);
      return;
    }

    // Close empty space between Q and P.
    #ifndef UseBoost
    if (k != mmax)
    #endif
      for (dim_t j = 0; j < k; ++j)
        memcpy(&AB(0, k + j), &P(0, j), nelec * sizeof(T));

    # if 0
    // Copy AB to ABC for TRMM interface.
    for (dim_t j = 0; j < 2 * k - 1; ++j)
      // inv(A)*U  [ 0 + + +
      //             0 0 + +
      //             0 0 0 +
      //             0 0 0 0 ] => AB[:, 0:2] -> ABC[:, 1:3]
      memcpy(&ABC(0, j + 1), &AB(0, j), nelec * sizeof(T));
    #ifdef UseBoost
    trmm(BLIS_RIGHT, BLIS_UPPER, BLIS_NO_TRANSPOSE, nelec, 2 * k - 1, T(1.0),
         &C(0, 1), C.size1(), &ABC(0, 1), ABC.size1());
    #else
    trmm(BLIS_RIGHT, BLIS_UPPER, BLIS_NO_TRANSPOSE, nelec, 2 * k - 1, T(1.0),
         &C(0, 1), C.ld, &ABC(0, 1), ABC.ld);
    #endif

    // Update: write to M.
    #ifdef UseBoost
    skr2k(uplo, BLIS_NO_TRANSPOSE, nelec, 2 * k - 1, T(1.0), &ABC(0, 1), ABC.size1(),
          &AB(0, 1), AB.size1(), T(1.0), &M(0, 0), M.size1());
    #else
    skr2k(uplo, BLIS_NO_TRANSPOSE, nelec, 2 * k - 1, T(1.0), &ABC(0, 1), ABC.ld,
          &AB(0, 1), AB.ld, T(1.0), &M(0, 0), M.ld);
    #endif
    #else
    gemm(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
         nelec, 2 * k, 2 * k,
         T(1.0),
         &AB(), AB.ld,
         &C(), C.ld,
         T(0.0),
         &ABC(), ABC.ld);
    gemmt(uplo, BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE,
          nelec, 2 * k,
          T(1.0),
          &ABC(), ABC.ld,
          &AB(), AB.ld,
          T(1.0),
          &M(), M.ld);
    #endif

    // Identity update to complete antisymmetric matrix.
    // This called due to lack of skmm support at the moment.
    // update_uplo(uplo);
  }

  /**
   * Complete antisymmetric matrix.
   */
  void skcomplete(uplo_t uplo_, dim_t n, matrix_t<T> &A) {
    for (dim_t j = 0; j < n; ++j) {
      for (dim_t i = 0; i < j; ++i) {
        switch (uplo_) {
        case BLIS_UPPER:
          A(j, i) = -A(i, j);
          break;

        case BLIS_LOWER:
          A(i, j) = -A(j, i);
          break;

        default:
          break;
        }
      }
      A(j, j) = T(0.0);
    }
  }

  void update_uplo(uplo_t uplo_new) {
    skcomplete(uplo /* NOTE: old uplo */, nelec, M);
    if (from_idx.size()) {
      skcomplete(uplo, from_idx.size(), UMU);
      skcomplete(uplo, from_idx.size(), VMV);
      skcomplete(uplo, from_idx.size(), W);
    }

    uplo = uplo_new;
  }

}; 
