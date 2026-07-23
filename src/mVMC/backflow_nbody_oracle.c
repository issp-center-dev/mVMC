#include "./include/backflow_nbody.h"

#include "./include/global.h"
#include "./include/matrix.h"
#include "./include/projection.h"
#include "./include/qp.h"
#include "./include/slater.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int enabled;
  int useFsz;
  FILE *stream;
  BFNBodyScratchSizes sizes;
  BFNBodyScratch scratch;
  int *intBase;
  double complex *complexBase;
  double *doubleBase;
  double complex *slaterBefore;
  double *slaterRealBefore;
  double complex *invBefore;
  double complex *pfBefore;
  double complex *etaBefore;
  int *etaFlagBefore;
  size_t slaterCount;
  size_t invCount;
  size_t etaCount;
} BFNBodyOracle;

static int BFNBodyOracleSizeMul(size_t left, size_t right, size_t *value) {
  if(value == NULL || (right != 0 && left > SIZE_MAX/right)) return -1;
  *value = left*right;
  return 0;
}

static void BFNBodyOracleClose(BFNBodyOracle *oracle) {
  if(oracle == NULL) return;
  if(oracle->stream != NULL) fclose(oracle->stream);
  free(oracle->intBase);
  free(oracle->complexBase);
  free(oracle->doubleBase);
  free(oracle->slaterBefore);
  free(oracle->slaterRealBefore);
  free(oracle->invBefore);
  free(oracle->pfBefore);
  free(oracle->etaBefore);
  free(oracle->etaFlagBefore);
  memset(oracle, 0, sizeof(*oracle));
}

static int BFNBodyOracleOpen(BFNBodyOracle *oracle, int parentRank,
                             int parentSize, int useFsz) {
  const char *dumpValue;
  char dumpPath[1024];
  int maxOrder;
  int pathLength;
  size_t nsizeSquare;

  if(oracle == NULL || (useFsz != 0 && useFsz != 1)) return -1;
  memset(oracle, 0, sizeof(*oracle));
  oracle->useFsz = useFsz;
  dumpValue = getenv("MVMC_BF_NBODY_ORACLE_DUMP");
  if(dumpValue == NULL || dumpValue[0] == '\0'
     || strcmp(dumpValue, "0") == 0) {
    return 0;
  }
  if(NVMCCalMode != 1 || AllComplexFlag == 0 || NBackFlowIdx <= 0
     || parentRank < 0 || parentSize < 1 || parentRank >= parentSize) {
    return -1;
  }

  maxOrder = NBodyGMaxN > NBodyInterAllMaxN
           ? NBodyGMaxN : NBodyInterAllMaxN;
  if(maxOrder < 1
     || GetBFNBodyScratchSizes(maxOrder, useFsz, &oracle->sizes)
          != BF_NBODY_OK
     || BFNBodyOracleSizeMul((size_t)Nsize, (size_t)Nsize,
                             &nsizeSquare) != 0
     || BFNBodyOracleSizeMul((size_t)NQPFull, nsizeSquare,
                             &oracle->invCount) != 0
     || BFNBodyOracleSizeMul((size_t)Nsite, (size_t)Nsite,
                             &oracle->etaCount) != 0) {
    return -1;
  }
  oracle->slaterCount = oracle->sizes.slaterCount;

  oracle->intBase =
      (int *)malloc(oracle->sizes.intCount*sizeof(int));
  oracle->complexBase = (double complex *)malloc(
      oracle->sizes.complexCount*sizeof(double complex));
  oracle->doubleBase =
      (double *)malloc(oracle->sizes.doubleCount*sizeof(double));
  oracle->slaterBefore = (double complex *)malloc(
      oracle->slaterCount*sizeof(double complex));
  oracle->slaterRealBefore =
      (double *)malloc(oracle->slaterCount*sizeof(double));
  oracle->invBefore =
      (double complex *)malloc(oracle->invCount*sizeof(double complex));
  oracle->pfBefore =
      (double complex *)malloc((size_t)NQPFull*sizeof(double complex));
  oracle->etaBefore =
      (double complex *)malloc(oracle->etaCount*sizeof(double complex));
  oracle->etaFlagBefore =
      (int *)malloc(oracle->etaCount*sizeof(int));
  if(oracle->intBase == NULL || oracle->complexBase == NULL
     || oracle->doubleBase == NULL || oracle->slaterBefore == NULL
     || oracle->slaterRealBefore == NULL || oracle->invBefore == NULL
     || oracle->pfBefore == NULL || oracle->etaBefore == NULL
     || oracle->etaFlagBefore == NULL
     || BindBFNBodyScratch(
            &oracle->sizes, oracle->intBase, oracle->sizes.intCount,
            oracle->complexBase, oracle->sizes.complexCount,
            oracle->doubleBase, oracle->sizes.doubleCount,
            &oracle->scratch) != BF_NBODY_OK) {
    BFNBodyOracleClose(oracle);
    return -1;
  }

  if(parentSize > 1) {
    pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s.rank%04d",
                          dumpValue, parentRank);
  } else {
    pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s", dumpValue);
  }
  if(pathLength < 0 || (size_t)pathLength >= sizeof(dumpPath)) {
    BFNBodyOracleClose(oracle);
    return -1;
  }
  oracle->stream = fopen(dumpPath, "w");
  if(oracle->stream == NULL) {
    BFNBodyOracleClose(oracle);
    return -1;
  }
  oracle->enabled = 1;
  return 0;
}

static void BFNBodyOracleSaveGlobals(BFNBodyOracle *oracle) {
  int site;

  memcpy(oracle->slaterBefore, SlaterElmBF,
         oracle->slaterCount*sizeof(double complex));
  memcpy(oracle->slaterRealBefore, SlaterElmBF_real,
         oracle->slaterCount*sizeof(double));
  memcpy(oracle->invBefore, InvM,
         oracle->invCount*sizeof(double complex));
  memcpy(oracle->pfBefore, PfM,
         (size_t)NQPFull*sizeof(double complex));
  for(site=0;site<Nsite;site++) {
    memcpy(oracle->etaBefore+(size_t)site*(size_t)Nsite, eta[site],
           (size_t)Nsite*sizeof(double complex));
    memcpy(oracle->etaFlagBefore+(size_t)site*(size_t)Nsite, etaFlag[site],
           (size_t)Nsite*sizeof(int));
  }
}

static void BFNBodyOracleRestoreGlobals(const BFNBodyOracle *oracle) {
  int site;

  memcpy(SlaterElmBF, oracle->slaterBefore,
         oracle->slaterCount*sizeof(double complex));
  memcpy(SlaterElmBF_real, oracle->slaterRealBefore,
         oracle->slaterCount*sizeof(double));
  memcpy(InvM, oracle->invBefore,
         oracle->invCount*sizeof(double complex));
  memcpy(PfM, oracle->pfBefore,
         (size_t)NQPFull*sizeof(double complex));
  for(site=0;site<Nsite;site++) {
    memcpy(eta[site], oracle->etaBefore+(size_t)site*(size_t)Nsite,
           (size_t)Nsite*sizeof(double complex));
    memcpy(etaFlag[site],
           oracle->etaFlagBefore+(size_t)site*(size_t)Nsite,
           (size_t)Nsite*sizeof(int));
  }
}

static int BFNBodyOracleGlobalsDiffer(const BFNBodyOracle *oracle) {
  int site;

  if(memcmp(SlaterElmBF, oracle->slaterBefore,
            oracle->slaterCount*sizeof(double complex)) != 0
     || memcmp(SlaterElmBF_real, oracle->slaterRealBefore,
               oracle->slaterCount*sizeof(double)) != 0
     || memcmp(InvM, oracle->invBefore,
               oracle->invCount*sizeof(double complex)) != 0
     || memcmp(PfM, oracle->pfBefore,
               (size_t)NQPFull*sizeof(double complex)) != 0) {
    return 1;
  }
  for(site=0;site<Nsite;site++) {
    if(memcmp(eta[site],
              oracle->etaBefore+(size_t)site*(size_t)Nsite,
              (size_t)Nsite*sizeof(double complex)) != 0
       || memcmp(etaFlag[site],
                 oracle->etaFlagBefore+(size_t)site*(size_t)Nsite,
                 (size_t)Nsite*sizeof(int)) != 0) {
      return 1;
    }
  }
  return 0;
}

static int BFNBodyOracleCanonicalEleIdx(int *canonical,
                                       const int *eleNum) {
  int spin;

  if(canonical == NULL || eleNum == NULL) return -1;
  for(spin=0;spin<2;spin++) {
    int label = 0;
    int site;
    for(site=0;site<Nsite;site++) {
      if(eleNum[site+spin*Nsite] == 1) {
        if(label >= Ne) return -1;
        canonical[label+spin*Ne] = site;
        label++;
      }
    }
    if(label != Ne) return -1;
  }
  return 0;
}

static int BFNBodyOracleFSZCanonicalParity(
    const int *eleIdx, const int *eleSpn, const int *eleNum,
    int *parity) {
  int value = 1;
  int particle;

  if(eleIdx == NULL || eleSpn == NULL || eleNum == NULL
     || parity == NULL) {
    return -1;
  }
  for(particle=0;particle<Nsize;particle++) {
    const int orbital =
        eleIdx[particle]+eleSpn[particle]*Nsite;
    int earlier;
    if(eleIdx[particle] < 0 || eleIdx[particle] >= Nsite
       || eleSpn[particle] < 0 || eleSpn[particle] > 1
       || orbital < 0 || orbital >= Nsite2
       || eleNum[orbital] != 1) {
      return -1;
    }
    for(earlier=0;earlier<particle;earlier++) {
      const int earlierOrbital =
          eleIdx[earlier]+eleSpn[earlier]*Nsite;
      if(earlierOrbital == orbital) return -1;
      if(earlierOrbital > orbital) value = -value;
    }
  }
  *parity = value;
  return 0;
}

static int BFNBodyOracleProjectedAmplitude(
    BFNBodyOracle *oracle, const int *eleIdx, const int *eleSpn,
    const int *eleNum, double complex *psi) {
  double complex ipFull;
  double logProj;
  int canonicalParity = 1;

  if(oracle == NULL || eleNum == NULL || psi == NULL
     || (oracle->useFsz
         && (eleIdx == NULL || eleSpn == NULL))) {
    return -1;
  }
  MakeProjCnt(oracle->scratch.projCnt, eleNum);
  MakeProjBFCnt(oracle->scratch.projBFCnt, eleNum);
  if(oracle->useFsz) {
    MakeSlaterElmBF_fsz(eleNum, oracle->scratch.projBFCnt);
    if(CalculateMAll_BF_fsz(eleIdx, eleSpn, 0, NQPFull) != 0
       || BFNBodyOracleFSZCanonicalParity(
              eleIdx, eleSpn, eleNum, &canonicalParity) != 0) {
      return -1;
    }
  } else {
    if(BFNBodyOracleCanonicalEleIdx(
           oracle->scratch.eleSpn, eleNum) != 0) {
      return -1;
    }
    MakeSlaterElmBF_fcmp(eleNum, oracle->scratch.projBFCnt);
    if(CalculateMAll_BF_fcmp(
           oracle->scratch.eleSpn, 0, NQPFull) != 0) {
      return -1;
    }
  }
  ipFull = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
  logProj = LogProjVal(oracle->scratch.projCnt);
  *psi = (double)canonicalParity*exp(logProj)*ipFull;
  return isfinite(creal(*psi)) && isfinite(cimag(*psi)) ? 0 : -1;
}

/*
 * Apply the original factor pairs, not their reduced representation.
 * The Python oracle independently supplies the fermionic sign.
 */
static int BFNBodyOracleApplyOriginal(
    BFNBodyOracle *oracle, int n, const int *rsi, const int *rsj,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleSpn) {
  int k;

  memcpy(oracle->scratch.eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(oracle->scratch.eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(oracle->scratch.eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  if(oracle->useFsz) {
    if(eleSpn == NULL) return -1;
    memcpy(oracle->scratch.eleSpn, eleSpn, (size_t)Nsize*sizeof(int));
  }

  for(k=n-1;k>=0;k--) {
    const int target = rsi[k];
    const int source = rsj[k];
    const int sourceSpin = source/Nsite;
    const int targetSpin = target/Nsite;
    int localLabel;
    int particle;

    if(source < 0 || source >= Nsite2 || target < 0 || target >= Nsite2
       || (!oracle->useFsz && sourceSpin != targetSpin)) {
      return -1;
    }
    if(oracle->scratch.eleNum[source] == 0) return 0;
    localLabel = oracle->scratch.eleCfg[source];
    particle = oracle->useFsz
             ? localLabel : localLabel+sourceSpin*Ne;
    if((oracle->useFsz
        ? (localLabel < 0 || localLabel >= Nsize)
        : (localLabel < 0 || localLabel >= Ne))
       || particle < 0 || particle >= Nsize
       || oracle->scratch.eleIdx[particle] != source%Nsite) {
      return -1;
    }
    if(oracle->useFsz
       && oracle->scratch.eleSpn[particle] != sourceSpin) {
      return -1;
    }

    oracle->scratch.eleCfg[source] = -1;
    oracle->scratch.eleNum[source] = 0;
    if(oracle->scratch.eleNum[target] != 0) return 0;
    oracle->scratch.eleCfg[target] =
        oracle->useFsz ? particle : localLabel;
    oracle->scratch.eleNum[target] = 1;
    oracle->scratch.eleIdx[particle] = target%Nsite;
    if(oracle->useFsz) {
      oracle->scratch.eleSpn[particle] = targetSpin;
    }
  }
  return 1;
}

static int BFNBodyOracleWriteRow(
    BFNBodyOracle *oracle, int sample, char sourceKind, int term, int n,
    const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt,
    double complex basePsi, double complex coefficient) {
  BFNBodyResult result;
  double complex targetPsi = 0.0+0.0*I;
  int targetValid;
  int orbital;

  if(oracle->useFsz) {
    result = GreenFuncNBF_fsz(
        n, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
        eleProjCnt, eleSpn, eleProjBFCnt, &oracle->scratch);
  } else {
    result = GreenFuncNBF(
        n, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
        eleProjCnt, eleProjBFCnt, &oracle->scratch);
  }
  if(result.status != BF_NBODY_OK
     && result.status != BF_NBODY_PHYSICAL_ZERO) {
    return -1;
  }

  targetValid = BFNBodyOracleApplyOriginal(
      oracle, n, rsi, rsj, eleIdx, eleCfg, eleNum, eleSpn);
  if(targetValid < 0) return -1;
  if(targetValid == 0) {
    if(result.status != BF_NBODY_PHYSICAL_ZERO
       || cabs(result.value) != 0.0) {
      return -1;
    }
  } else {
    if(BFNBodyOracleProjectedAmplitude(
           oracle, oracle->scratch.eleIdx, oracle->scratch.eleSpn,
           oracle->scratch.eleNum, &targetPsi) != 0) {
      BFNBodyOracleRestoreGlobals(oracle);
      return -1;
    }
    BFNBodyOracleRestoreGlobals(oracle);
  }

  fprintf(oracle->stream,
          "sample %d source %c term %d n %d status %d base_occ",
          sample, sourceKind, term, n, (int)result.status);
  for(orbital=0;orbital<Nsite2;orbital++) {
    fprintf(oracle->stream, " %d", eleNum[orbital]);
  }
  fprintf(oracle->stream, " target_valid %d target_occ", targetValid);
  for(orbital=0;orbital<Nsite2;orbital++) {
    fprintf(oracle->stream, " %d", oracle->scratch.eleNum[orbital]);
  }
  fprintf(oracle->stream,
          " base_psi %.17e %.17e target_psi %.17e %.17e"
          " value %.17e %.17e coeff %.17e %.17e\n",
          creal(basePsi), cimag(basePsi),
          creal(targetPsi), cimag(targetPsi),
          creal(result.value), cimag(result.value),
          creal(coefficient), cimag(coefficient));
  return ferror(oracle->stream) ? -1 : 0;
}

static int BFNBodyOracleDumpSample(
    BFNBodyOracle *oracle, int sample, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleSpn,
    const int *eleProjBFCnt) {
  double complex basePsi;
  int term;
  int k;

  if(oracle == NULL || !oracle->enabled || oracle->stream == NULL
     || eleIdx == NULL || eleCfg == NULL || eleNum == NULL
     || (NProj > 0 && eleProjCnt == NULL) || eleProjBFCnt == NULL
     || (oracle != NULL && oracle->useFsz && eleSpn == NULL)
     || !isfinite(creal(ip)) || !isfinite(cimag(ip)) || cabs(ip) == 0.0) {
    return -1;
  }

  BFNBodyOracleSaveGlobals(oracle);
  if(BFNBodyOracleProjectedAmplitude(
         oracle, eleIdx, eleSpn, eleNum, &basePsi) != 0
     || cabs(basePsi) == 0.0) {
    BFNBodyOracleRestoreGlobals(oracle);
    return -1;
  }
  BFNBodyOracleRestoreGlobals(oracle);

  for(term=0;term<NNBodyG;term++) {
    const int n = NBodyGN[term];
    const int offset = NBodyGOffset[term];
    for(k=0;k<n;k++) {
      oracle->scratch.inputRsi[k] =
          NBodyGIdx[offset+k][0]+NBodyGIdx[offset+k][1]*Nsite;
      oracle->scratch.inputRsj[k] =
          NBodyGIdx[offset+k][2]+NBodyGIdx[offset+k][3]*Nsite;
    }
    if(BFNBodyOracleWriteRow(
           oracle, sample, 'g', term, n,
           oracle->scratch.inputRsi, oracle->scratch.inputRsj,
           ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
           eleProjBFCnt,
           basePsi, 1.0+0.0*I) != 0) {
      BFNBodyOracleRestoreGlobals(oracle);
      return -1;
    }
  }
  for(term=0;term<NNBodyInterAll;term++) {
    const int n = NBodyInterAllN[term];
    const int offset = NBodyInterAllOffset[term];
    for(k=0;k<n;k++) {
      oracle->scratch.inputRsi[k] =
          NBodyInterAllIdx[offset+k][0]
          +NBodyInterAllIdx[offset+k][1]*Nsite;
      oracle->scratch.inputRsj[k] =
          NBodyInterAllIdx[offset+k][2]
          +NBodyInterAllIdx[offset+k][3]*Nsite;
    }
    if(BFNBodyOracleWriteRow(
           oracle, sample, 'h', term, n,
           oracle->scratch.inputRsi, oracle->scratch.inputRsj,
           ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
           eleProjBFCnt,
           basePsi, ParaNBodyInterAll[term]) != 0) {
      BFNBodyOracleRestoreGlobals(oracle);
      return -1;
    }
  }

  BFNBodyOracleRestoreGlobals(oracle);
  if(BFNBodyOracleGlobalsDiffer(oracle) || fflush(oracle->stream) != 0) {
    return -1;
  }
  return 0;
}
