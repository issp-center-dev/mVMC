/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <cstdlib>
#include <iostream>


inline void abort(const char *msg)
{
  std::cout << msg << std::endl;
  abort();
}

inline void abort(const char *msg0, const char *msg1)
{
  std::cout << msg0 << ": " << msg1 << std::endl;
  abort();
}

#if (defined(_DEBUG) || !defined(SKIP_CHECKS))

inline void assert_(const bool cond, const char *msg)
{
  if (!cond) abort(msg);
}

inline void assert_(const bool cond, const char *msg0, const char *msg1)
{
  if (!cond) abort(msg0, msg1);
}

#else

inline void assert_(const bool cond, const char *msg) { }

inline void assert_(const bool cond, const char *msg0, const char *msg1) { }

#endif

