/**
 * \copyright Copyright (c) Dept. Phys., Univ. Tokyo
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "error.hh"
#include <random>
#include <vector>
#include <cstdint>

namespace vmc
{
struct config_manager
{
  using index_t = int32_t;
  using base_t = std::vector<index_t>;

  const index_t nsites;
  index_t norbs()
  { return nsites * 2; }

  const  base_t &config_base()          const { return _cfg_base; }
  const index_t  config_base(index_t i) const { return _cfg_base.at(i); }
  const  base_t &   from_idx()          const { return _from_idx; }
  const index_t     from_idx(index_t i) const { return _from_idx.at(i); }
  const  base_t &     to_orb()          const { return _to_orb; }
  const index_t       to_orb(index_t i) const { return _to_orb.at(i); }
  const  base_t &     config()          const { return _cfg; }
  const index_t       config(index_t i) const { return _cfg.at(i); }
  const  base_t &        num()          const { return _num; }
  const index_t          num(index_t i) const { return _num.at(i); }
  const  base_t &        inv()          const { return _inv; }
  const index_t          inv(index_t i) const { return _inv.at(i); }
  const index_t         spin(index_t i) const { return _num.at(i) == 1 ? 1 : _num.at(i) == 2 ? -1 : 0; }

  void reserve_space(index_t mmax)
  {
    _from_idx.reserve(mmax * sizeof(index_t));
    _to_orb.reserve(mmax * sizeof(index_t));
  }

  void attach_config(const base_t &cfg)
  {
    _cfg_base = cfg;
    _cfg = cfg;
    _from_idx.clear();
    _to_orb.clear();

    // Reset _num.
    for (index_t i = 0; i < nsites; ++i)
      _num.at(i) = 0;

    // Generate num.
    for (index_t i = 0; i < cfg.size(); ++i) {
      // 0b01 for up, 0b10 for down, 0b11 for double occupation.
      index_t xi = cfg.at(i) % nsites;
      index_t si = cfg.at(i) / nsites;
      _num.at(xi) += (si + 1);
    }

    // Generate inv.
    for (index_t i = 0; i < _inv.size(); ++i)
      _inv.at(i) = -1;
    for (index_t i = 0; i < cfg.size(); ++i)
      _inv.at(cfg.at(i)) = i;
  }

  void merge_config()
  {
    // OR: cfg_base = cfg;
    for (index_t j = 0; j < _from_idx.size(); ++j)
      _cfg_base.at(_from_idx.at(j)) = _to_orb.at(j);
    _from_idx.clear();
    _to_orb.clear();
  }

  void push_update(index_t osi, index_t msj)
  {
    for (auto i = _from_idx.begin(); i != _from_idx.end(); ++i)
      assert_(*i != msj, typeid(*this).name(), "Consecutive hopping unsupported. Please merge manually before proceeding.");
    index_t osj = _cfg.at(msj);
    index_t xj = osj % nsites;
    index_t sj = osj / nsites;
    index_t xi = osi % nsites;
    index_t si = osi / nsites;
    assert_(_inv.at(osj) == msj && _inv.at(osi) == -1, typeid(*this).name(), "Inv error @ push.");

    _from_idx.push_back(msj);
    _to_orb.push_back(osi);
    _cfg.at(msj) = osi;
    // Update num.
    _num.at(xj) -= sj + 1;
    _num.at(xi) += si + 1;
    // Update inv.
    _inv.at(osj) = -1;
    _inv.at(osi) = msj;
  }

  void pop_update()
  {
    assert_(_from_idx.size(), typeid(*this).name(), "Trying to pop empty updates.");
    index_t msj = _from_idx.back();
    index_t osj = _cfg_base.at(msj);
    index_t xj  = osj % nsites;
    index_t sj  = osj / nsites;
    index_t osi = _cfg.at(msj);
    index_t xi  = osi % nsites;
    index_t si  = osi / nsites;
    assert_(_inv.at(osi) == msj && _inv.at(osj) == -1, typeid(*this).name(), "Inv error @ pop.");

    _cfg.at(msj) = osj;
    // Update num.
    _num.at(xi) -= si + 1;
    _num.at(xj) += sj + 1;
    // Update inv.
    _inv.at(osi) = -1;
    _inv.at(osj) = msj;
    _from_idx.pop_back();
    _to_orb.pop_back();
  }

  void check()
  {
    assert_(_cfg.size() <= norbs(), typeid(*this).name(), "More Fermions than orbitals.");

    base_t num(nsites, 0);
    for (index_t i = 0; i < num.size(); ++i) {
      index_t xi = _cfg.at(i) % nsites;
      index_t si = _cfg.at(i) / nsites;

      assert_(!(num.at(xi) & (si + 1)), typeid(*this).name(), "More than two Fermions occupying the same orbital.");
      num.at(xi) += (si + 1);
    }
    for (index_t i = 0; i < nsites; ++i)
      assert_(num.at(i) == _num.at(i), typeid(*this).name(), "Found broken member: _num.");

    base_t cfg_curr = _cfg_base;
    for (index_t j = 0; j < _from_idx.size(); ++j)
      cfg_curr.at(_from_idx.at(j)) = _to_orb.at(j);
    for (index_t i = 0; i < cfg_curr.size(); ++i)
      assert_(cfg_curr.at(i) == _cfg.at(i), typeid(*this).name(), "Found broken member: _cfg.");
    // Check inv.
  }

  config_manager(index_t nsites_)
  : nsites(nsites_), _num(nsites_), _inv(norbs()) {  }

  config_manager(index_t nsites_, base_t cfg)
  : nsites(nsites_), _num(nsites_), _inv(norbs())
  { attach_config(cfg); }

  template <typename rng_t>
  static base_t singlet_config(unsigned nsite, unsigned nsinglet, rng_t &rng)
  {
    using dist_t = std::uniform_int_distribution<unsigned>;

    bool allow_doublon = nsinglet * 2 > nsite;
    base_t marks(nsite, 0);
    base_t cfg(0);

    // Spin up.
    for (int i = 0; i < nsinglet; ++i) {
      int iseq = dist_t(0, nsite - i - 1)(rng);

      int isite;
      int isite_c = 0;
      do {
        // Find the next available site.
        while (marks.at(isite_c)) { isite_c++; }
        isite = isite_c;
        isite_c++;
      } while (iseq-- > 0);

      marks.at(isite) = 1;
      cfg.push_back(isite);
    }

    // Spin down.
    for (int i = 0; i < nsinglet; ++i) {
      int iseq;
      if (allow_doublon)
        iseq = dist_t(0, nsite - i - 1)(rng);
      else
        iseq = dist_t(0, nsite - i - 1 - nsinglet)(rng);

      int isite;
      int isite_c = 0;
      do {
        // Find the next unoccupied site / orb.
        if (allow_doublon)
          while (marks.at(isite_c) > 1) { isite_c++; }
        else
          while (marks.at(isite_c) > 0) { isite_c++; }
        isite = isite_c;
        isite_c++;
      } while (iseq-- > 0);

      marks.at(isite) = 2;
      cfg.push_back(isite + nsite);
    }
    return cfg;
  }

  template <typename rng_t>
  std::pair<int, int> propose_hop(rng_t &rng) const
  {
    int loop_breaker = 0;
    int xi;
    int msj = std::uniform_int_distribution<unsigned>(0, _cfg.size() - 1)(rng);
    std::uniform_int_distribution<unsigned> select_site(0, nsites - 1);
    int sj = msj / nsites;

    do {
      assert_(loop_breaker++ < 10000000, typeid(*this).name(), "Too many loops.");
      xi = select_site(rng);
    } while (_num.at(xi) & (sj + 1));

    return std::make_pair(xi + sj * nsites, msj);
  }

  template <typename rng_t>
  std::tuple<int, int, int, int> propose_exc(rng_t &rng) const
  {
    using dist_t = std::uniform_int_distribution<unsigned>;

    int nup = 0;
    int ndn = 0;
    int loop_breaker = 0;
    int msj, xj = 0;
    int msl, xl = 0;

    // Count # of spin-up/down sites.
    for (int i = 0; i < _num.size(); ++i)
      if (num(i) == 1)
        nup++;
      else if (num(i) == 2)
        ndn++;

    // Pick a singly occupied spin-up site.
    int iseq_up = dist_t(0, nup - 1)(rng);
    while (true) {
      while (num(xj) != 1) ++xj;
      if (--iseq_up < 0)
        break;
      else
        ++xj;
    }

    // Pick a singly occupied spin-down site.
    int iseq_dn = dist_t(0, ndn - 1)(rng);
    while (true) {
      while (num(xl) != 2) ++xl;
      if (--iseq_dn < 0)
        break;
      else
        ++xl;
    }

    // Reverse lookup config indices.
    msj = inv(xj);
    msl = inv(xl + nsites);
    assert_(msj >= 0 && msl >= 0, typeid(*this).name(), "ExchangeConfig error.");

    return std::make_tuple(xl, msj, xj + nsites, msl);
  }

private:
  base_t _cfg_base;
  base_t _from_idx;
  base_t _to_orb;
  base_t _cfg;
  base_t _num;
  base_t _inv;
};
}

