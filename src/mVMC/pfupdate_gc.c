#include <complex.h>
#include <math.h>
#include <stddef.h>

#include "include/blas_externs.h"
#include "include/gc_size.h"
#include "include/global.h"
#include "include/pfupdate_gc.h"
#include "include/workspace.h"

#ifndef _SRC_PFUPDATE_GC
#define _SRC_PFUPDATE_GC

static size_t GCInvStride(void) {
  return GCCheckedMulSize((size_t)NsizeMax, (size_t)NsizeMax);
}

static size_t GCSlaterStride(void) {
  return GCCheckedMulSize((size_t)Nsite2, (size_t)Nsite2);
}

/* The proposal has already replaced eleIdx[ma]. InvM/PfM still describe the
 * old configuration. The active inverse uses leading dimension NsizeMax. */
void CalculateNewPfMHopGC(const int ma, double complex *pfMNew,
                          const int *eleIdx, const int ncur,
                          const int qpStart, const int qpEnd) {
  const size_t invStride = GCInvStride();
  const size_t slaterStride = GCSlaterStride();
  int qpidx;
  for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
    const int rsa = eleIdx[ma];
    const double complex *slater =
        SlaterElm + (size_t)(qpStart + qpidx) * slaterStride;
    const double complex *slaterRow =
        slater + (size_t)rsa * (size_t)Nsite2;
    const double complex *invRow =
        InvM + (size_t)qpidx * invStride +
        (size_t)ma * (size_t)NsizeMax;
    double complex ratio = 0.0;
    int msj;
    for (msj = 0; msj < ncur; msj++) {
      ratio += invRow[msj] * slaterRow[eleIdx[msj]];
    }
    pfMNew[qpidx] = ratio * PfM[qpidx];
  }
}

static void UpdateMAllHopGCChild(
    const int ma, const int *eleIdx, const int ncur, const int qpStart,
    const int qpidx, double complex *vec1, double complex *vec2) {
  const size_t invStride = GCInvStride();
  const size_t slaterStride = GCSlaterStride();
  const int rsa = eleIdx[ma];
  const double complex *slater =
      SlaterElm + (size_t)(qpStart + qpidx) * slaterStride;
  const double complex *slaterRow =
      slater + (size_t)rsa * (size_t)Nsite2;
  double complex *invM = InvM + (size_t)qpidx * invStride;
  double complex *invRow = invM + (size_t)ma * (size_t)NsizeMax;
  double complex pivot;
  int msi;
  int msj;

  for (msi = 0; msi < ncur; msi++) vec1[msi] = 0.0;
  for (msj = 0; msj < ncur; msj++) {
    const double complex value = slaterRow[eleIdx[msj]];
    const double complex *invColumnAsRow =
        invM + (size_t)msj * (size_t)NsizeMax;
    for (msi = 0; msi < ncur; msi++) {
      vec1[msi] += -invColumnAsRow[msi] * value;
    }
  }

  pivot = vec1[ma];
  PfM[qpidx] *= pivot;
  for (msi = 0; msi < ncur; msi++) vec2[msi] = invRow[msi];
  for (msi = 0; msi < ncur; msi++) {
    double complex *invMRow;
    const double complex oldColumn = -vec2[msi];
    if (msi == ma) continue;
    invMRow = invM + (size_t)msi * (size_t)NsizeMax;
    for (msj = 0; msj < ncur; msj++) {
      if (msj == ma) continue;
      invMRow[msj] -=
          (vec1[msj] * oldColumn + vec1[msi] * vec2[msj]) / pivot;
    }
  }
  invRow[ma] = 0.0;
  for (msi = 0; msi < ncur; msi++) {
    double complex *invMRow =
        invM + (size_t)msi * (size_t)NsizeMax;
    if (msi == ma) continue;
    invRow[msi] = vec2[msi] / pivot;
    invMRow[ma] = -invRow[msi];
  }
}

void UpdateMAllHopGC(const int ma, const int *eleIdx, const int ncur,
                     const int qpStart, const int qpEnd) {
  const int workspaceCount = GCCheckedSizeToInt(
      GCCheckedMulSize(2, (size_t)ncur), "GC hop update workspace");
  int qpidx;
  RequestWorkSpaceThreadComplex(workspaceCount);
#pragma omp parallel default(shared)
  {
    double complex *vec1 = GetWorkSpaceThreadComplex(ncur);
    double complex *vec2 = GetWorkSpaceThreadComplex(ncur);
#pragma omp for
    for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
      UpdateMAllHopGCChild(ma, eleIdx, ncur, qpStart, qpidx, vec1, vec2);
    }
  }
  ReleaseWorkSpaceThreadComplex();
}

static double complex CalculateNewPfMTwoHopGCChild(
    const int ma, const int mb, const int *eleIdx, const int ncur,
    const int qpStart, const int qpidx, double complex *vecA,
    double complex *vecB) {
  const size_t invStride = GCInvStride();
  const size_t slaterStride = GCSlaterStride();
  const double complex *slater =
      SlaterElm + (size_t)(qpStart + qpidx) * slaterStride;
  const double complex *slaterA =
      slater + (size_t)eleIdx[ma] * (size_t)Nsite2;
  const double complex *slaterB =
      slater + (size_t)eleIdx[mb] * (size_t)Nsite2;
  const double complex *invM = InvM + (size_t)qpidx * invStride;
  const double complex *invA =
      invM + (size_t)ma * (size_t)NsizeMax;
  const double complex *invB =
      invM + (size_t)mb * (size_t)NsizeMax;
  double complex pA = 0.0;
  double complex pB = 0.0;
  double complex qA = 0.0;
  double complex qB = 0.0;
  double complex bMa = 0.0;
  double complex ratio;
  int i;

  for (i = 0; i < ncur; i++) {
    vecA[i] = slaterA[eleIdx[i]];
    vecB[i] = slaterB[eleIdx[i]];
    pA += -invA[i] * vecA[i];
    pB += -invB[i] * vecA[i];
    qA += -invA[i] * vecB[i];
    qB += -invB[i] * vecB[i];
  }
  for (i = 0; i < ncur; i++) {
    const double complex *invRow =
        invM + (size_t)i * (size_t)NsizeMax;
    double complex inner = 0.0;
    int j;
    for (j = 0; j < ncur; j++) inner += -invRow[j] * vecA[j];
    bMa += vecB[i] * inner;
  }
  ratio = -invA[mb] * vecB[ma] - invA[mb] * bMa +
          pA * qB - pB * qA;
  return ratio * PfM[qpidx];
}

void CalculateNewPfMTwoHopGC(const int ma, const int mb,
                             double complex *pfMNew, const int *eleIdx,
                             const int ncur, const int qpStart,
                             const int qpEnd) {
  const int workspaceCount = GCCheckedSizeToInt(
      GCCheckedMulSize(2, (size_t)ncur), "GC two-hop workspace");
  int qpidx;
  if (ma == mb) {
    CalculateNewPfMHopGC(ma, pfMNew, eleIdx, ncur, qpStart, qpEnd);
    return;
  }
  RequestWorkSpaceComplex(workspaceCount);
  {
    double complex *vecA = GetWorkSpaceComplex(ncur);
    double complex *vecB = GetWorkSpaceComplex(ncur);
    for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
      pfMNew[qpidx] = CalculateNewPfMTwoHopGCChild(
          ma, mb, eleIdx, ncur, qpStart, qpidx, vecA, vecB);
    }
  }
  ReleaseWorkSpaceComplex();
}

double complex CalculateNewPfMNGC(
    const int qpidx, const int n, const int *msa, const int *rsa,
    const int *eleIdx, const int ncur, double complex *vec,
    double complex *w, double complex *smallMat, int *pfIWork,
    double complex *pfWork, const int pfLWork, double *pfRWork) {
  size_t smallCount;
  const size_t invStride = GCInvStride();
  const size_t slaterStride = GCSlaterStride();
  const double complex *slater;
  const double complex *invM;
  double complex pfaffian;
  double sign;
  char uplo = 'U';
  char method = 'P';
  int n2;
  int lda;
  int lwork = pfLWork;
  int info = 0;
  int k;

  if (n <= 0 || ncur <= 0 || ncur > NsizeMax || qpidx < 0 ||
      qpidx >= NQPFull || msa == NULL || rsa == NULL || eleIdx == NULL ||
      vec == NULL || w == NULL || smallMat == NULL || pfIWork == NULL ||
      pfWork == NULL || pfRWork == NULL) {
    return NAN + NAN * I;
  }
  n2 = GCCheckedSizeToInt(GCCheckedMulSize(2, (size_t)n),
                          "GC n-row dimension");
  smallCount = GCCheckedMulSize((size_t)n2, (size_t)n2);
  if ((size_t)pfLWork < smallCount) return NAN + NAN * I;
  for (k = 0; k < n; k++) {
    if (msa[k] < 0 || msa[k] >= ncur || rsa[k] < 0 || rsa[k] >= Nsite2) {
      return NAN + NAN * I;
    }
  }

  slater = SlaterElm + (size_t)qpidx * slaterStride;
  invM = InvM + (size_t)qpidx * invStride;
  for (k = 0; k < n; k++) {
    double complex *vecK = vec + (size_t)k * (size_t)NsizeMax;
    const double complex *slaterK =
        slater + (size_t)rsa[k] * (size_t)Nsite2;
    int i;
    for (i = 0; i < ncur; i++) vecK[i] = slaterK[eleIdx[i]];
  }

  for (k = 0; k < n; k++) {
    double complex *columnK = smallMat + (size_t)n2 * (size_t)k;
    double complex *vecK = vec + (size_t)k * (size_t)NsizeMax;
    int l;
    for (l = k + 1; l < n; l++) {
      const double complex *vecL =
          vec + (size_t)l * (size_t)NsizeMax;
      double complex value = vecK[msa[l]];
      int i;
      for (i = 0; i < ncur; i++) {
        const double complex *invRow =
            invM + (size_t)i * (size_t)NsizeMax;
        double complex inner = 0.0;
        int j;
        for (j = 0; j < ncur; j++) inner += -invRow[j] * vecL[j];
        w[i] = inner;
      }
      for (i = 0; i < ncur; i++) value += w[i] * vecK[i];
      columnK[l] = value;
    }
  }
  for (k = 0; k < n; k++) {
    double complex *columnK =
        smallMat + (size_t)n2 * (size_t)k + (size_t)n;
    const double complex *vecK =
        vec + (size_t)k * (size_t)NsizeMax;
    int l;
    for (l = 0; l < n; l++) {
      const double complex *invL =
          invM + (size_t)msa[l] * (size_t)NsizeMax;
      double complex value = 0.0;
      int i;
      for (i = 0; i < ncur; i++) value += vecK[i] * (-invL[i]);
      columnK[l] = value;
    }
  }
  for (k = 0; k < n; k++) {
    double complex *columnK =
        smallMat + (size_t)n2 * (size_t)(k + n) + (size_t)n;
    const double complex *invK =
        invM + (size_t)msa[k] * (size_t)NsizeMax;
    int l;
    for (l = k + 1; l < n; l++) columnK[l] = -invK[msa[l]];
  }
  for (k = 0; k < n2; k++) {
    int l;
    for (l = 0; l < k; l++) {
      smallMat[(size_t)n2 * (size_t)k + (size_t)l] =
          -smallMat[(size_t)n2 * (size_t)l + (size_t)k];
    }
    smallMat[(size_t)n2 * (size_t)k + (size_t)k] = 0.0;
  }

  lda = n2;
  M_ZSKPFA(&uplo, &method, &n2, smallMat, &lda, &pfaffian, pfIWork,
           pfWork, &lwork, pfRWork, &info);
  if (info != 0 || !isfinite(creal(pfaffian)) ||
      !isfinite(cimag(pfaffian))) {
    return NAN + NAN * I;
  }
  sign = (((n * (n - 1) / 2) % 2) == 0) ? 1.0 : -1.0;
  return sign * pfaffian * PfM[qpidx];
}

#endif
