/*-------------------------------------------------------------
 * Variational Monte Carlo
 * fast Pfaffian update
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

void CalculateNewPfM(const int mi, const int s, double *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd);
void CalculateNewPfM2(const int mi, const int s, double *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd);
void UpdateMAll(const int mi, const int s, const int *eleIdx,
                const int qpStart, const int qpEnd);
void updateMAll_child(const int ma, const int s, const int *eleIdx,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double *vec1, double *vec2);

/* Calculate new pfaffian. The ma-th electron with spin s hops. */
void CalculateNewPfM(const int ma, const int s, double *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd) {
  #pragma procedure serial
  const int qpNum = qpEnd-qpStart;
  const int msa = ma+s*Ne;
  const int rsa = eleIdx[msa] + s*Nsite;

  int qpidx;
  int msj,rsj;
  const double *sltE_a; /* update elements of msa-th row */
  const double *invM_a;
  double ratio;

  /* optimization for Kei */
  const int nsize = Nsize;
  const int ne = Ne;

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    sltE_a = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;
    invM_a = InvM + qpidx*Nsize*Nsize + msa*Nsize;

    ratio = 0.0;
    for(msj=0;msj<ne;msj++) {
      rsj = eleIdx[msj];
      ratio += invM_a[msj] * sltE_a[rsj];
    }
    for(msj=ne;msj<nsize;msj++) {
      rsj = eleIdx[msj] + Nsite;
      ratio += invM_a[msj] * sltE_a[rsj];
    }

    pfMNew[qpidx] = -ratio*PfM[qpidx];
  }

  return;
}

/* thread parallel version of CalculateNewPfM */
void CalculateNewPfM2(const int ma, const int s, double *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  const int msa = ma+s*Ne;
  const int rsa = eleIdx[msa] + s*Nsite;

  int qpidx;
  int msj,rsj;
  const double *sltE_a; /* update elements of msa-th row */
  const double *invM_a;
  double ratio;

  /* optimization for Kei */
  const int nsize = Nsize;
  const int ne = Ne;

  #pragma omp parallel for default(shared)        \
    private(qpidx,msj,sltE_a,invM_a,ratio,rsj)
  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    sltE_a = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;
    invM_a = InvM + qpidx*Nsize*Nsize + msa*Nsize;

    ratio = 0.0;
    for(msj=0;msj<ne;msj++) {
      rsj = eleIdx[msj];
      ratio += invM_a[msj] * sltE_a[rsj];
    }
    for(msj=ne;msj<nsize;msj++) {
      rsj = eleIdx[msj] + Nsite;
      ratio += invM_a[msj] * sltE_a[rsj];
    }

    pfMNew[qpidx] = -ratio*PfM[qpidx];
  }

  return;
}

/* Update PfM and InvM. The ma-th electron with spin s hops to site ra=eleIdx[msi] */
void UpdateMAll(const int ma, const int s, const int *eleIdx,
                const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  double *vec1,*vec2;

  RequestWorkSpaceThreadDouble(2*Nsize);

  #pragma omp parallel default(shared) private(vec1,vec2)
  {
    vec1 = GetWorkSpaceThreadDouble(Nsize);
    vec2 = GetWorkSpaceThreadDouble(Nsize);
   
    #pragma omp for private(qpidx)
    #pragma loop nounroll
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      updateMAll_child(ma, s, eleIdx, qpStart, qpEnd, qpidx, vec1, vec2);
    }
  }

  ReleaseWorkSpaceThreadDouble();
  return;
}

void updateMAll_child(const int ma, const int s, const int *eleIdx,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double *vec1, double *vec2) {
  #pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  const int msa = ma+s*Ne;
  const int rsa = eleIdx[msa] + s*Nsite;
  const int nsize = Nsize; /* optimization for Kei */

  int msi,msj,rsj;

  const double *sltE_a; /* update elements of msa-th row */
  double sltE_aj;
  double *invM;
  double *invM_i,*invM_j,*invM_a;

  double vec1_i,vec2_i;
  double invVec1_a;
  double tmp;

  sltE_a = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;

  invM = InvM + qpidx*Nsize*Nsize;
  invM_a = invM + msa*Nsize;

  for(msi=0;msi<nsize;msi++) vec1[msi] = 0.0;

  /* Calculate vec1[i] = sum_j invM[i][j] sltE[a][j] */
  /* Note tah invM[i][j] = -invM[j][i] */
  #pragma loop noalias
  for(msj=0;msj<nsize;msj++) {
    rsj = eleIdx[msj] + (msj/Ne)*Nsite;
    sltE_aj = sltE_a[rsj];
    invM_j = invM + msj*Nsize;

    for(msi=0;msi<nsize;msi++) {
      vec1[msi] += -invM_j[msi] * sltE_aj;
    }
  }

  /* Update Pfaffian */
  /* Calculate -1.0/bufV_a to reduce devision */
  tmp = vec1[msa];
  PfM[qpidx] *= -tmp;
  invVec1_a = -1.0/tmp;

  /* Calculate vec2[i] = -InvM[a][i]/vec1[a] */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    vec2[msi] = invM_a[msi] * invVec1_a;
  }

  /* Update InvM */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    vec1_i = vec1[msi];
    vec2_i = vec2[msi];

    for(msj=0;msj<nsize;msj++) {
      invM_i[msj] += vec1_i * vec2[msj] - vec1[msj] * vec2_i;
    }

    invM_i[msa] -= vec2_i;
  }

  #pragma loop noalias
  for(msj=0;msj<nsize;msj++) {
    invM_a[msj] += vec2[msj];
  }
  /* end of update invM */

  return;
}
