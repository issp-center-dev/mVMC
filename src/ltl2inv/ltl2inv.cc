/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "pfaffian.tcc"
#include "invert.tcc"



extern "C" {

void ltl2inv_s(int n, float    *A_, int ldA, int *iPiv, float    *vT, float    *M_, int ldM) { ltl2inv<float   >(n, A_, ldA, iPiv, vT, M_, ldM); }
void ltl2inv_d(int n, double   *A_, int ldA, int *iPiv, double   *vT, double   *M_, int ldM) { ltl2inv<double  >(n, A_, ldA, iPiv, vT, M_, ldM); }
void ltl2inv_c(int n, ccscmplx *A_, int ldA, int *iPiv, ccscmplx *vT, ccscmplx *M_, int ldM) { ltl2inv<ccscmplx>(n, A_, ldA, iPiv, vT, M_, ldM); }
void ltl2inv_z(int n, ccdcmplx *A_, int ldA, int *iPiv, ccdcmplx *vT, ccdcmplx *M_, int ldM) { ltl2inv<ccdcmplx>(n, A_, ldA, iPiv, vT, M_, ldM); }
void utu2inv_s(int n, float    *A_, int ldA, int *iPiv, float    *vT, float    *M_, int ldM) { utu2inv<float   >(n, A_, ldA, iPiv, vT, M_, ldM); }
void utu2inv_d(int n, double   *A_, int ldA, int *iPiv, double   *vT, double   *M_, int ldM) { utu2inv<double  >(n, A_, ldA, iPiv, vT, M_, ldM); }
void utu2inv_c(int n, ccscmplx *A_, int ldA, int *iPiv, ccscmplx *vT, ccscmplx *M_, int ldM) { utu2inv<ccscmplx>(n, A_, ldA, iPiv, vT, M_, ldM); }
void utu2inv_z(int n, ccdcmplx *A_, int ldA, int *iPiv, ccdcmplx *vT, ccdcmplx *M_, int ldM) { utu2inv<ccdcmplx>(n, A_, ldA, iPiv, vT, M_, ldM); }


void ltl2pfa_s(int n, float    *A_, int ldA, int *iPiv, float    *Pfa) { *Pfa = ltl2pfa<float   >(n, A_, ldA, iPiv); }
void ltl2pfa_d(int n, double   *A_, int ldA, int *iPiv, double   *Pfa) { *Pfa = ltl2pfa<double  >(n, A_, ldA, iPiv); }
void ltl2pfa_c(int n, ccscmplx *A_, int ldA, int *iPiv, ccscmplx *Pfa) { *Pfa = ltl2pfa<ccscmplx>(n, A_, ldA, iPiv); }
void ltl2pfa_z(int n, ccdcmplx *A_, int ldA, int *iPiv, ccdcmplx *Pfa) { *Pfa = ltl2pfa<ccdcmplx>(n, A_, ldA, iPiv); }
void utu2pfa_s(int n, float    *A_, int ldA, int *iPiv, float    *Pfa) { *Pfa = utu2pfa<float   >(n, A_, ldA, iPiv); }
void utu2pfa_d(int n, double   *A_, int ldA, int *iPiv, double   *Pfa) { *Pfa = utu2pfa<double  >(n, A_, ldA, iPiv); }
void utu2pfa_c(int n, ccscmplx *A_, int ldA, int *iPiv, ccscmplx *Pfa) { *Pfa = utu2pfa<ccscmplx>(n, A_, ldA, iPiv); }
void utu2pfa_z(int n, ccdcmplx *A_, int ldA, int *iPiv, ccdcmplx *Pfa) { *Pfa = utu2pfa<ccdcmplx>(n, A_, ldA, iPiv); }

}

