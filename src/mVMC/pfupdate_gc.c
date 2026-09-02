#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "include/blas_externs.h"
#include "include/gc_config.h"
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

void CalculateNewPfMTwoHopGCWorkspace(
    const int ma, const int mb, double complex *pfMNew,
    const int *eleIdx, const int ncur, const int qpStart, const int qpEnd,
    double complex *vecA, double complex *vecB) {
  int qpidx;
  if (ma == mb) {
    CalculateNewPfMHopGC(ma, pfMNew, eleIdx, ncur, qpStart, qpEnd);
    return;
  }
  for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
    pfMNew[qpidx] = CalculateNewPfMTwoHopGCChild(
        ma, mb, eleIdx, ncur, qpStart, qpidx, vecA, vecB);
  }
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

static double complex GCAddSchur(
    const int rsa, const int rsb, const int *eleIdx, const int ncurOld,
    const int qpStart, const int qpidx, double complex *y0,
    double complex *y1) {
  const size_t invStride = GCInvStride();
  const size_t slaterStride = GCSlaterStride();
  const double complex *slater =
      SlaterElm + (size_t)(qpStart + qpidx) * slaterStride;
  const double complex *slaterA =
      slater + (size_t)rsa * (size_t)Nsite2;
  const double complex *slaterB =
      slater + (size_t)rsb * (size_t)Nsite2;
  const double complex *invM = InvM + (size_t)qpidx * invStride;
  double complex d01 = -slaterA[rsb];
  int i;
  for (i = 0; i < ncurOld; i++) {
    const double complex *invRow =
        invM + (size_t)i * (size_t)NsizeMax;
    double complex value0 = 0.0;
    double complex value1 = 0.0;
    int j;
    for (j = 0; j < ncurOld; j++) {
      const int rsj = eleIdx[j];
      value0 += invRow[j] * slaterA[rsj];
      value1 += invRow[j] * slaterB[rsj];
    }
    y0[i] = value0;
    y1[i] = value1;
  }
  for (i = 0; i < ncurOld; i++) d01 += slaterA[eleIdx[i]] * y1[i];
  return d01;
}

void CalculateNewPfMAddGC(const int rsa, const int rsb,
                          double complex *pfMNew, const int *eleIdx,
                          const int ncurOld, const int qpStart,
                          const int qpEnd) {
  const int workspaceCount = GCCheckedSizeToInt(
      GCCheckedMulSize(2, (size_t)ncurOld), "GC add candidate workspace");
  RequestWorkSpaceComplex(workspaceCount);
  {
    double complex *y0 = GetWorkSpaceComplex(ncurOld);
    double complex *y1 = GetWorkSpaceComplex(ncurOld);
    CalculateNewPfMAddGCWorkspace(rsa, rsb, pfMNew, eleIdx, ncurOld,
                                  qpStart, qpEnd, y0, y1);
  }
  ReleaseWorkSpaceComplex();
}

void CalculateNewPfMAddGCWorkspace(
    const int rsa, const int rsb, double complex *pfMNew,
    const int *eleIdx, const int ncurOld, const int qpStart,
    const int qpEnd, double complex *y0, double complex *y1) {
  int qpidx;
  for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
    const double complex d01 =
        GCAddSchur(rsa, rsb, eleIdx, ncurOld, qpStart, qpidx, y0, y1);
    pfMNew[qpidx] = d01 * PfM[qpidx];
  }
}

static void UpdateMAllAddGCChild(
    const int rsa, const int rsb, const int *eleIdx, const int ncurOld,
    const int qpStart, const int qpidx, double complex *y0,
    double complex *y1) {
  const size_t invStride = GCInvStride();
  double complex *invM = InvM + (size_t)qpidx * invStride;
  const double complex d01 =
      GCAddSchur(rsa, rsb, eleIdx, ncurOld, qpStart, qpidx, y0, y1);
  int i;
  PfM[qpidx] *= d01;
  for (i = 0; i < ncurOld; i++) {
    double complex *invRow =
        invM + (size_t)i * (size_t)NsizeMax;
    int j;
    for (j = 0; j < ncurOld; j++) {
      invRow[j] += (-y0[i] * y1[j] + y1[i] * y0[j]) / d01;
    }
    invRow[ncurOld] = -y1[i] / d01;
    invRow[ncurOld + 1] = y0[i] / d01;
  }
  {
    double complex *row0 =
        invM + (size_t)ncurOld * (size_t)NsizeMax;
    double complex *row1 =
        invM + (size_t)(ncurOld + 1) * (size_t)NsizeMax;
    for (i = 0; i < ncurOld; i++) {
      row0[i] = y1[i] / d01;
      row1[i] = -y0[i] / d01;
    }
    row0[ncurOld] = 0.0;
    row0[ncurOld + 1] = -1.0 / d01;
    row1[ncurOld] = 1.0 / d01;
    row1[ncurOld + 1] = 0.0;
  }
}

void UpdateMAllAddGC(const int rsa, const int rsb, const int *eleIdx,
                     const int ncurOld, const int qpStart,
                     const int qpEnd) {
  const int workspaceCount = GCCheckedSizeToInt(
      GCCheckedMulSize(2, (size_t)ncurOld), "GC add update workspace");
  int qpidx;
  RequestWorkSpaceThreadComplex(workspaceCount);
#pragma omp parallel default(shared)
  {
    double complex *y0 = GetWorkSpaceThreadComplex(ncurOld);
    double complex *y1 = GetWorkSpaceThreadComplex(ncurOld);
#pragma omp for
    for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
      UpdateMAllAddGCChild(rsa, rsb, eleIdx, ncurOld, qpStart, qpidx, y0,
                           y1);
    }
  }
  ReleaseWorkSpaceThreadComplex();
}

static int GCOldPositionAfterRemoveSwap(const int newPosition,
                                        const int pos0, const int pos1,
                                        const int ncurOld) {
  int oldPosition = newPosition;
  if (oldPosition == pos0) {
    oldPosition = ncurOld - 2;
  } else if (oldPosition == ncurOld - 2) {
    oldPosition = pos0;
  }
  if (oldPosition == pos1) {
    oldPosition = ncurOld - 1;
  } else if (oldPosition == ncurOld - 1) {
    oldPosition = pos1;
  }
  return oldPosition;
}

static double complex GCRemoveRatio(const int pos0, const int pos1,
                                    const int ncurOld, const int qpidx) {
  const size_t invStride = GCInvStride();
  const double complex *invM = InvM + (size_t)qpidx * invStride;
  const int tail0 = GCOldPositionAfterRemoveSwap(
      ncurOld - 2, pos0, pos1, ncurOld);
  const int tail1 = GCOldPositionAfterRemoveSwap(
      ncurOld - 1, pos0, pos1, ncurOld);
  const double complex r =
      invM[(size_t)tail0 * (size_t)NsizeMax + (size_t)tail1];
  return -(double)GCRemoveParitySign(pos0, pos1, ncurOld) * r;
}

void CalculateNewPfMRemoveGC(const int pos0, const int pos1,
                             double complex *pfMNew, const int *eleIdx,
                             const int ncurOld, const int qpStart,
                             const int qpEnd) {
  int qpidx;
  (void)eleIdx;
  (void)qpStart;
  for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
    pfMNew[qpidx] =
        GCRemoveRatio(pos0, pos1, ncurOld, qpidx) * PfM[qpidx];
  }
}

static void UpdateMAllRemoveGCChild(
    const int pos0, const int pos1, const int ncurOld, const int qpidx,
    double complex *oldInv) {
  const size_t invStride = GCInvStride();
  const int survivor = ncurOld - 2;
  double complex *invM = InvM + (size_t)qpidx * invStride;
  const int tail0 = GCOldPositionAfterRemoveSwap(
      survivor, pos0, pos1, ncurOld);
  const int tail1 = GCOldPositionAfterRemoveSwap(
      survivor + 1, pos0, pos1, ncurOld);
  double complex r;
  int i;
  for (i = 0; i < ncurOld; i++) {
    memcpy(oldInv + (size_t)i * (size_t)ncurOld,
           invM + (size_t)i * (size_t)NsizeMax,
           (size_t)ncurOld * sizeof(*oldInv));
  }
  r = oldInv[(size_t)tail0 * (size_t)ncurOld + (size_t)tail1];
  PfM[qpidx] *= -(double)GCRemoveParitySign(pos0, pos1, ncurOld) * r;
  for (i = 0; i < survivor; i++) {
    const int oldI =
        GCOldPositionAfterRemoveSwap(i, pos0, pos1, ncurOld);
    double complex *newRow =
        invM + (size_t)i * (size_t)NsizeMax;
    int j;
    for (j = 0; j < survivor; j++) {
      const int oldJ =
          GCOldPositionAfterRemoveSwap(j, pos0, pos1, ncurOld);
      const double complex p =
          oldInv[(size_t)oldI * (size_t)ncurOld + (size_t)oldJ];
      const double complex qi0 =
          oldInv[(size_t)oldI * (size_t)ncurOld + (size_t)tail0];
      const double complex qi1 =
          oldInv[(size_t)oldI * (size_t)ncurOld + (size_t)tail1];
      const double complex u0j =
          oldInv[(size_t)tail0 * (size_t)ncurOld + (size_t)oldJ];
      const double complex u1j =
          oldInv[(size_t)tail1 * (size_t)ncurOld + (size_t)oldJ];
      newRow[j] = p + (qi0 * u1j - qi1 * u0j) / r;
    }
  }
}

void UpdateMAllRemoveGC(const int pos0, const int pos1,
                        const int *eleIdx, const int ncurOld,
                        const int qpStart, const int qpEnd) {
  const size_t matrixCount =
      GCCheckedMulSize((size_t)ncurOld, (size_t)ncurOld);
  const int workspaceCount =
      GCCheckedSizeToInt(matrixCount, "GC remove update workspace");
  int qpidx;
  (void)eleIdx;
  (void)qpStart;
  RequestWorkSpaceThreadComplex(workspaceCount);
#pragma omp parallel default(shared)
  {
    double complex *oldInv = GetWorkSpaceThreadComplex(workspaceCount);
#pragma omp for
    for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
      UpdateMAllRemoveGCChild(pos0, pos1, ncurOld, qpidx, oldInv);
    }
  }
  ReleaseWorkSpaceThreadComplex();
}

#endif
