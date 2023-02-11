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
  const int ith = omp_get_thread_num(); \
  const int nth = omp_get_num_threads(); \
\
  /* Count Green's funcs. */ \
  int gf_count = 0; \
  for (int i = 0; i < num_gf; ++i) \
    if ( needs_comput[i] ) \
      gf_count++; \
\
  /* Partitioning. */ \
  int gfCountThread = (gf_count + nth - 1) / nth; \
  const int gfStart = gfCountThread * ith; \
  gfCountThread = std::min(gf_count - gfStart, gfCountThread); \
\
  /* Pack info space. */ \
  gf_count = -gfStart; \
  int *to_orbs_thread_, *from_ids_thread_, *unpack_idx_thread; \
  for (int i = 0; i < num_gf && gf_count < gfCountThread; ++i) \
    if ( needs_comput[i] ) { \
      if ( gf_count >= 0 ) { \
        if ( gf_count == 0 ) { \
          /* Pin down scratchpad spots. */ \
          to_orbs_thread_   = to_orbs  + i * 2; \
          from_ids_thread_  = from_ids + i * 2; \
          unpack_idx_thread = unpack_idx + i; \
          unpack_idx_thread[0] = i; \
        } else { \
          /* 0th idx already in the place. Pack from 1st. */ \
          to_orbs_thread_[gf_count * 2    ] = to_orbs[i * 2]; \
          to_orbs_thread_[gf_count * 2 + 1] = to_orbs[i * 2 + 1]; \
          from_ids_thread_[gf_count * 2    ] = from_ids[i * 2]; \
          from_ids_thread_[gf_count * 2 + 1] = from_ids[i * 2 + 1]; \
          unpack_idx_thread[gf_count] = i; \
        } \
      } \
      gf_count++; \
    } \
\
  Eigen::Map<VectorXi>  to_orbs_thread( to_orbs_thread_, gfCountThread * 2); \
  Eigen::Map<VectorXi> from_ids_thread(from_ids_thread_, gfCountThread * 2); \
\
  for (int iqp = 0; iqp < num_qp; ++iqp) { \
    Eigen::Vector<ctype, Eigen::Dynamic> pfaBatch = \
      objv(iqp, ctype)->batch_query_amplitudes(2, to_orbs_thread, from_ids_thread); \
    pfaBatch *= objv(iqp, ctype)->get_amplitude(); \
\
    /* Unpack computed Pfaffians. */ \
    for (int igf = 0; igf < gfCountThread; ++igf) \
      pfav[unpack_idx_thread[igf] * num_qp + iqp] = pfaBatch[igf]; \
  } \
}
// GENIMPL( float,    s )
GENIMPL( double,   d )
// GENIMPL( ccscmplx, c )
GENIMPL( ccdcmplx, z )
#undef GENIMPL

#undef EXPANDNAME

