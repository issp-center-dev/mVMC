/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "updated_tdi.tcc"
// In implementation, this should appear AFTER blis.h.
#include "pf_interface.h"

// Non-# Pragma support differs from compiler to compiler
#if defined(__INTEL_COMPILER)
#define OMP_PARALLEL_FOR_SHARED __pragma(omp parallel for default(shared))
#elif defined(__GNUC__)
#define OMP_PARALLEL_FOR_SHARED _Pragma("omp parallel for default(shared)")
#else
#error "Valid non-preprocessor _pragma() not found."
#endif

#define orbv( i, ctype ) ( (orbital_mat<ctype> *)orbv[i] )
#define objv( i, ctype ) ( (updated_tdi<ctype> *)objv[i] )

#define EXPANDNAME( funcname, cblachar ) funcname##_##cblachar

#define GENIMPL( ctype, cblachar ) \
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
      void     *orbv[] ) \
{ \
  OMP_PARALLEL_FOR_SHARED \
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    orbv[iqp] = new orbital_mat<ctype>( \
        BLIS_UPPER, norbs, orbmat_base + iqp * orbmat_stride, norbs); \
    objv[iqp] = new updated_tdi<ctype>( \
        *orbv(iqp, ctype), nelec, invmat_base + iqp * invmat_stride, nelec, mmax); \
\
    /* Compute initial. */ \
    auto &cfg_i = objv(iqp, ctype)->elem_cfg; \
    for (int msi = 0; msi < nelec; ++msi) \
      cfg_i.at(msi) = eleidx[msi] + elespn[msi]*nsite; \
    objv(iqp, ctype)->initialize(); \
  } \
}
GENIMPL( float,    s )
GENIMPL( double,   d )
GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#define GENIMPL( ctype, cblachar ) \
  void EXPANDNAME( updated_tdi_v_free, cblachar ) \
    ( uint64_t  num_qp, \
      void     *objv[], \
      void     *orbv[] ) \
{ \
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    delete orbv(iqp, ctype); \
    delete objv(iqp, ctype); \
  } \
}
GENIMPL( float,    s )
GENIMPL( double,   d )
GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#define GENIMPL( ctype, cblachar ) \
  void EXPANDNAME( updated_tdi_v_get_pfa, cblachar ) \
    ( uint64_t  num_qp, \
      ctype     pfav[], \
      void     *objv[] ) \
{ \
  /* Minus sign appears from that SlaterElm is in fact interpreted as row-major. */ \
  for (int iqp = 0; iqp < num_qp; ++iqp) \
    pfav[iqp] = -objv(iqp, ctype)->get_Pfa(); \
}
GENIMPL( float,    s )
GENIMPL( double,   d )
GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#define GENIMPL( ctype, cblachar ) \
  void EXPANDNAME( updated_tdi_v_push, cblachar ) \
    ( uint64_t  num_qp, \
      int64_t   osi, \
      int64_t   msj, \
      int64_t   cal_pfa, \
      void     *objv[] ) \
{ \
  OMP_PARALLEL_FOR_SHARED \
  for (int iqp = 0; iqp < num_qp; ++iqp) \
    objv(iqp, ctype)->push_update_safe(osi, msj, cal_pfa!=0); \
}
GENIMPL( float,    s )
GENIMPL( double,   d )
GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#define GENIMPL( ctype, cblachar ) \
  void EXPANDNAME( updated_tdi_v_push_pair, cblachar ) \
    ( uint64_t  num_qp, \
      int64_t   osi, \
      int64_t   msj, \
      int64_t   osk, \
      int64_t   msl, \
      int64_t   cal_pfa, \
      void     *objv[] ) \
{ \
  OMP_PARALLEL_FOR_SHARED \
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    auto &from_i = objv(iqp, ctype)->from_idx; \
    if (from_i.size() > objv(iqp, ctype)->mmax - 2) \
      objv(iqp, ctype)->merge_updates(); \
    else \
      for (int ik = 0; ik < from_i.size(); ++ik) \
        if (from_i.at(ik) == msj || \
            from_i.at(ik) == msl) { \
          objv(iqp, ctype)->merge_updates(); \
          break; \
        } \
    objv(iqp, ctype)->push_update(osi, msj, false); \
    objv(iqp, ctype)->push_update(osk, msl, cal_pfa!=0); \
  } \
}
GENIMPL( float,    s )
GENIMPL( double,   d )
GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#define GENIMPL( ctype, cblachar ) \
  void EXPANDNAME( updated_tdi_v_pop, cblachar ) \
    ( uint64_t  num_qp, \
      int64_t   cal_pfa, \
      void     *objv[] ) \
{ \
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    objv(iqp, ctype)->pop_update(cal_pfa!=0); \
  } \
}
GENIMPL( float,    s )
GENIMPL( double,   d )
GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#undef EXPANDNAME

