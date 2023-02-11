/**
 * To link to this module, the following C/C++-compiler feature is required
 *  BEYOND C99/C++14 standard:
 *  - all pointers has uint64 underlying datatype.
 *  This is usually satisfied if one uses 64-bit compilers with compatible ABI.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "stdint.h"

// Redefine double complex for (future) non-C99 compatibility.
#ifndef __cplusplus
typedef _Complex double ccdcmplx;
typedef _Complex float  ccscmplx;
#else
extern "C" {
#endif

// Function naming scheme:
// updated_tdi_[v means ``operate on vector'']_[operation]_[datatype].

#define EXPANDNAME( funcname, cblachar ) funcname##_##cblachar

#define GENDEF( ctype, cblachar ) \
  void EXPANDNAME( updated_tdi_v_seq_init_precomp, cblachar ) \
    ( uint64_t  num_qp, \
      uint64_t  nsite, \
      uint64_t  norbs, \
      uint64_t  nelec, \
      ctype    *orbmat_base, \
      int64_t   orbmat_stride, \
      ctype    *invmat_base, \
      int64_t   invmat_stride, \
      int32_t  *eleidx, \
      int32_t  *elespn, \
      uint64_t  mmax, \
      ctype     pfav[], \
      void     *objv[], \
      void     *orbv[], \
      void     *matv[], \
      void     *mapv[] );
// GENDEF( float,    s )
GENDEF( double,   d )
// GENDEF( ccscmplx, c )
GENDEF( ccdcmplx, z )
#undef GENDEF


#define GENDEF( ctype, cblachar ) \
   void EXPANDNAME( updated_tdi_v_init, cblachar ) \
    ( uint64_t  num_qp, \
      uint64_t  nsite, \
      uint64_t  norbs, \
      uint64_t  nelec, \
      ctype    *orbmat_base, \
      int64_t   orbmat_stride, \
      ctype    *invmat_base, \
      int64_t   invmat_stride, \
      int32_t  *eleidx, \
      int32_t  *elespn, \
      uint64_t  mmax, \
      void     *objv[], \
      void     *orbv[], \
      void     *matv[], \
      void     *mapv[] );

// GENDEF( float,    s )
GENDEF( double,   d )
// GENDEF( ccscmplx, c )
GENDEF( ccdcmplx, z )
#undef GENDEF

#define GENDEF( ctype, cblachar ) \
   void EXPANDNAME( updated_tdi_v_free, cblachar ) \
    ( uint64_t  num_qp, \
      void     *objv[], \
      void     *orbv[], \
      void     *matv[], \
      void     *mapv[] );

// GENDEF( float,    s )
GENDEF( double,   d )
// GENDEF( ccscmplx, c )
GENDEF( ccdcmplx, z )
#undef GENDEF

#define GENDEF( ctype, cblachar ) \
   void EXPANDNAME( updated_tdi_v_get_pfa, cblachar ) \
    ( uint64_t  num_qp, \
      ctype     pfav[], \
      void     *objv[] );

// GENDEF( float,    s )
GENDEF( double,   d )
// GENDEF( ccscmplx, c )
GENDEF( ccdcmplx, z )
#undef GENDEF

#define GENDEF( ctype, cblachar ) \
   void EXPANDNAME( updated_tdi_v_push, cblachar ) \
    ( uint64_t  num_qp, \
      int64_t   osi, \
      int64_t   msj, \
      int64_t   cal_pfa, \
      void     *objv[] );

// GENDEF( float,    s )
GENDEF( double,   d )
// GENDEF( ccscmplx, c )
GENDEF( ccdcmplx, z )
#undef GENDEF

#define GENDEF( ctype, cblachar ) \
   void EXPANDNAME( updated_tdi_v_push_pair, cblachar ) \
    ( uint64_t  num_qp, \
      int64_t   osi, \
      int64_t   msj, \
      int64_t   osk, \
      int64_t   msl, \
      int64_t   cal_pfa, \
      void     *objv[] );

// GENDEF( float,    s )
GENDEF( double,   d )
// GENDEF( ccscmplx, c )
GENDEF( ccdcmplx, z )
#undef GENDEF

#define GENDEF( ctype, cblachar ) \
   void EXPANDNAME( updated_tdi_v_pop, cblachar ) \
    ( uint64_t  num_qp, \
      int64_t   cal_pfa, \
      void     *objv[] );

// GENDEF( float,    s )
GENDEF( double,   d )
// GENDEF( ccscmplx, c )
GENDEF( ccdcmplx, z )
#undef GENDEF

#define GENDEF( ctype, cblachar ) \
  void EXPANDNAME( updated_tdi_v_omp_proc_batch_greentwo, cblachar ) \
    ( uint64_t  num_qp, \
      int       num_gf, \
      int       needs_comput[], \
      int       unpack_idx[], \
      int       to_orbs[], \
      int       from_ids[], \
      ctype     pfav[], \
      void     *objv[], \
      void     *orbv[], \
      void     *matv[], \
      void     *mapv[] );
// GENDEF( float,    s )
GENDEF( double,   d )
// GENDEF( ccscmplx, c )
GENDEF( ccdcmplx, z )
#undef GENDEF

#undef EXPANDNAME

#ifdef __cplusplus
}
#endif

