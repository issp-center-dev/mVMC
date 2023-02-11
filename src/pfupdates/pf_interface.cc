/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "updated_tdi.tcc"
#include "orbital_mat.tcc"
// In implementation, this should appear AFTER blis.h.
#include "pf_interface.h"
#include <omp.h>

// Non-# Pragma support differs from compiler to compiler
#if defined(__INTEL_COMPILER)
#define _ccPragma(_1) __pragma(_1)
#elif defined(__GNUC__)
#define _ccPragma(_1) _Pragma(#_1)
#else
#error "Valid non-preprocessor _Pragma() not found."
#endif

using namespace Eigen;
using namespace vmc::orbital;
template <typename T>
using matrix_t = Map<Matrix<T, Dynamic, Dynamic>, 0, OuterStride<> >;

#define matv( i, ctype ) ( (                        matrix_t<ctype>     *)matv[i] )
#define mapv( i, ctype ) ( (                        matrix_t<ctype>     *)mapv[i] )
#define orbv( i, ctype ) ( (            orbital_mat<matrix_t<ctype> >   *)orbv[i] )
#define objv( i, ctype ) ( (updated_tdi<orbital_mat<matrix_t<ctype> > > *)objv[i] )

#define EXPANDNAME( funcname, cblachar ) funcname##_##cblachar

#define GENIMPL( ctype, cblachar ) \
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
      void     *mapv[] ) \
{ \
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    mapv[iqp] = new matrix_t<ctype>(orbmat_base + iqp * orbmat_stride, norbs, norbs, OuterStride<>(norbs)); \
    matv[iqp] = new matrix_t<ctype>(invmat_base + iqp * invmat_stride, nelec, nelec, OuterStride<>(nelec)); \
    orbv[iqp] = new orbital_mat<matrix_t<ctype>>('U', norbs, *mapv(iqp, ctype) ); \
    objv[iqp] = new updated_tdi<orbital_mat<matrix_t<ctype>>> \
        ( *orbv(iqp, ctype), *matv(iqp, ctype), nelec, mmax); \
\
    /* Compute initial. */ \
    vmc::config_manager::base_t cfg_i(nelec); \
    for (int msi = 0; msi < nelec; ++msi) \
      cfg_i.at(msi) = eleidx[msi] + elespn[msi]*nsite; \
    objv(iqp, ctype)->attach_config(cfg_i); \
    objv(iqp, ctype)->initialize_precomputed(pfav[iqp]); \
  } \
}
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

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
      void     *orbv[], \
      void     *matv[], \
      void     *mapv[] ) \
{ \
  _ccPragma(omp parallel for default(shared)) \
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    mapv[iqp] = new matrix_t<ctype>(orbmat_base + iqp * orbmat_stride, norbs, norbs, OuterStride<>(norbs)); \
    matv[iqp] = new matrix_t<ctype>(invmat_base + iqp * invmat_stride, nelec, nelec, OuterStride<>(nelec)); \
    orbv[iqp] = new orbital_mat<matrix_t<ctype>>('U', norbs, *mapv(iqp, ctype) ); \
    objv[iqp] = new updated_tdi<orbital_mat<matrix_t<ctype>>> \
        ( *orbv(iqp, ctype), *matv(iqp, ctype), nelec, mmax); \
\
    /* Compute initial. */ \
    vmc::config_manager::base_t cfg_i(nelec); \
    for (int msi = 0; msi < nelec; ++msi) \
      cfg_i.at(msi) = eleidx[msi] + elespn[msi]*nsite; \
    objv(iqp, ctype)->attach_config(cfg_i); \
    objv(iqp, ctype)->initialize(); \
  } \
}
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#define GENIMPL( ctype, cblachar ) \
  void EXPANDNAME( updated_tdi_v_free, cblachar ) \
    ( uint64_t  num_qp, \
      void     *objv[], \
      void     *orbv[], \
      void     *matv[], \
      void     *mapv[] ) \
{ \
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    delete orbv(iqp, ctype); \
    delete objv(iqp, ctype); \
    delete matv(iqp, ctype); \
    delete mapv(iqp, ctype); \
  } \
}
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
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
    pfav[iqp] = -objv(iqp, ctype)->get_amplitude(); \
}
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
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
  _ccPragma(omp parallel for default(shared)) \
  for (int iqp = 0; iqp < num_qp; ++iqp) \
    objv(iqp, ctype)->push_update_safe(osi, msj, cal_pfa!=0); \
}
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
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
  _ccPragma(omp parallel for default(shared)) \
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    const auto &from_i = objv(iqp, ctype)->get_config_manager().from_idx(); \
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
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
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
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#define GENIMPL( ctype, cblachar ) \
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
      void     *mapv[] ) \
{ \
  int ith = omp_get_thread_num(); \
  int nth = omp_get_num_threads(); \
\
  int gf_count = 0; \
  for (int i = 0; i < num_gf; ++i) \
    if ( needs_comput[i] ) { \
      if ( ith == 0 ) { \
        /* Pack info space. */ \
        to_orbs[gf_count * 2    ] = to_orbs[i * 2]; \
        to_orbs[gf_count * 2 + 1] = to_orbs[i * 2 + 1]; \
        from_ids[gf_count * 2    ] = from_ids[i * 2]; \
        from_ids[gf_count * 2 + 1] = from_ids[i * 2 + 1]; \
        unpack_idx[gf_count] = i; \
      } \
      /* Count Green's funcs. */ \
      gf_count++; \
    } \
\
  int gf_cnt_thread = (gf_count + nth - 1) / nth; \
  int gfStart = gf_cnt_thread * ith; \
  int gfEnd = std::min(gf_count, gfStart + gf_cnt_thread); \
\
  Eigen::Map<VectorXi>  to_orbs_thread( to_orbs + gfStart * 2, (gfEnd - gfStart) * 2); \
  Eigen::Map<VectorXi> from_ids_thread(from_ids + gfStart * 2, (gfEnd - gfStart) * 2); \
\
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    Eigen::Vector<ctype, Eigen::Dynamic> pfaBatch = \
      objv(iqp, ctype)->batch_query_amplitudes(2, to_orbs_thread, from_ids_thread); \
\
    /* Unpack computed Pfaffians. */ \
    for (int igf = gfStart; igf < gfEnd; ++igf) \
      pfav[unpack_idx[igf] * num_qp + iqp] = pfaBatch[igf - gfStart]; \
  } \
}
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#undef EXPANDNAME

