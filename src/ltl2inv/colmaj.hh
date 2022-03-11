/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once

template <typename T>
struct colmaj
{
  T *dat;
  int ld;

  colmaj() = delete;
  colmaj(T *dat_, int ld_)
  : dat(dat_), ld(ld_) { }

  T &operator()(int i, int j)
  {
    return dat[ i + j * ld ];
  }
  T &operator()()
  {
    return dat[ 0 ];
  }
};
