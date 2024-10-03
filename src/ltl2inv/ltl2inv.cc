/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "pfaffian.tcc"
#include "invert.tcc"



extern "C" {

#define GENDEF(schema, cchar, ctype, cctype) \
void schema ##2inv_## cchar(int n, cctype *A_, int ldA, int *iPiv, cctype *vT, cctype *M_, int ldM) \
{ \
  using namespace Eigen; \
  Map<Matrix<cctype, Dynamic, Dynamic>, 0, OuterStride<> > A(A_, n, n, OuterStride<>(ldA)); \
  Map<Matrix<cctype, Dynamic, Dynamic>, 0, OuterStride<> > M(M_, n, n, OuterStride<>(ldM)); \
  Map<VectorXi> Piv(iPiv, n); \
  Map<Vector<cctype, Dynamic>> T(vT, n); \
  schema ##2inv(A, Piv, vT, M); \
}
GENDEF(ltl, s, float,    float)
GENDEF(ltl, d, double,   double)
GENDEF(ltl, c, ccscmplx, ccscmplx)
GENDEF(ltl, z, ccdcmplx, ccdcmplx)
GENDEF(utu, s, float,    float)
GENDEF(utu, d, double,   double)
GENDEF(utu, c, ccscmplx, ccscmplx)
GENDEF(utu, z, ccdcmplx, ccdcmplx)

#undef GENDEF
#define GENDEF(schema, cchar, ctype, cctype) \
void schema ##2pfa_## cchar(int n, cctype *A_, int ldA, int *iPiv, cctype *Pfa) \
{ \
  using namespace Eigen; \
  Map<Matrix<cctype, Dynamic, Dynamic>, 0, OuterStride<> > A(A_, n, n, OuterStride<>(ldA)); \
  Map<VectorXi> Piv(iPiv, n); \
  *Pfa = schema ##2pfa(A, Piv); \
}
GENDEF(ltl, s, float,    float)
GENDEF(ltl, d, double,   double)
GENDEF(ltl, c, ccscmplx, ccscmplx)
GENDEF(ltl, z, ccdcmplx, ccdcmplx)
GENDEF(utu, s, float,    float)
GENDEF(utu, d, double,   double)
GENDEF(utu, c, ccscmplx, ccscmplx)
GENDEF(utu, z, ccdcmplx, ccdcmplx)

}

