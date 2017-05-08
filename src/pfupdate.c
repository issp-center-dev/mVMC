/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * fast Pfaffian update
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "pfupdate.h"
#ifndef _PFUDATE_SRC
#define _PFUDATE_SRC

double complex updateMAll_BF_fcmp_child(
        const int qpidx,
        const int n, const int *msa,
        const int *eleIdx);


/* Calculate new pfaffian. The ma-th electron with spin s hops. */
void CalculateNewPfM(const int ma, const int s, double complex *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd) {
  #pragma procedure serial
  const int qpNum = qpEnd-qpStart;
  const int msa = ma+s*Ne;
  const int rsa = eleIdx[msa] + s*Nsite;

  int qpidx;
  int msj,rsj;
  const double complex *sltE_a; /* update elements of msa-th row */
  const double complex *invM_a;
  double complex ratio;

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
void CalculateNewPfM2(const int ma, const int s, double complex *pfMNew, const int *eleIdx,
                     const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  const int msa = ma+s*Ne;
  const int rsa = eleIdx[msa] + s*Nsite;

  int qpidx;
  int msj,rsj;
  const double complex *sltE_a; /* update elements of msa-th row */
  const double complex *invM_a;
  double complex ratio;

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
      //printf("DEBUG:msj=%d rsj=%d: invM=%lf %lf : slt=%lf %lf \n",msj,rsj,creal(invM_a[msj]),cimag(invM_a[msj]),creal(sltE_a[rsj]),cimag(sltE_a[rsj]));
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
  double complex *vec1,*vec2;

  RequestWorkSpaceThreadComplex(2*Nsize);

  #pragma omp parallel default(shared) private(vec1,vec2)
  {
    vec1 = GetWorkSpaceThreadComplex(Nsize);
    vec2 = GetWorkSpaceThreadComplex(Nsize);
   
    #pragma omp for private(qpidx)
    #pragma loop nounroll
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      updateMAll_child(ma, s, eleIdx, qpStart, qpEnd, qpidx, vec1, vec2);
    }
  }

  ReleaseWorkSpaceThreadComplex();
  return;
}

void updateMAll_child(const int ma, const int s, const int *eleIdx,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double complex *vec1, double complex *vec2) {
  #pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  const int msa = ma+s*Ne;
  const int rsa = eleIdx[msa] + s*Nsite;
  const int nsize = Nsize; /* optimization for Kei */

  int msi,msj,rsj;

  const double complex *sltE_a; /* update elements of msa-th row */
  double complex sltE_aj;
  double complex *invM;
  double complex *invM_i,*invM_j,*invM_a;

  double complex vec1_i,vec2_i;
  double complex invVec1_a;
  double complex tmp;

  sltE_a = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;

  invM = InvM + qpidx*Nsize*Nsize;
  invM_a = invM + msa*Nsize;

  for(msi=0;msi<nsize;msi++) vec1[msi] = 0.0+0.0*I; //TBC

  /* Calculate vec1[i] = sum_j invM[i][j] sltE[a][j] */
  /* Note that invM[i][j] = -invM[j][i] */
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


/* Calculate new pfaffian with Backflow effects.
   The ma-th electron with spin s hops from ra to rb */
void CalculateNewPfMBF(const int *icount, const int *msaTmp,
                       double complex* pfMNew, const int *eleIdx,
                       const int qpStart, const int qpEnd, const double complex* bufM) {
  //#pragma procedure serial
  int i;
  const int nsize = Nsize;
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  int *msa;

  for(qpidx=0;qpidx<qpNum;qpidx++) {
    //Store msa//
    msa=(int *)malloc(sizeof(int)*icount[qpidx]);
    //printf("Total=%d\n",icount[qpidx]);
    for(i=0;i<icount[qpidx];i++){
      msa[i] = msaTmp[i+qpidx*Nsite];
      //printf("hop[%d]=%d\n",i,msa[i]);
    }

    /* calculateNewPfM */
    pfMNew[qpidx] = calculateNewPfMBFN4_child(qpidx,icount[qpidx],msa,eleIdx,bufM);

    free(msa);
  }

  return;
}

double complex calculateNewPfMBFN4_child(const int qpidx, const int n, const int *msa,
                                         const int *eleIdx, const double complex* bufM)
{
  const int nsize = Nsize;
  const int n2 = 2*n;
  const double complex *sltE;
  const double complex *sltE_k;
  double complex *invM;
  double complex *invM_i, *invM_k, *invM_l;

  double complex *vec; /* vec[n][nsize] */
  //double complex vec[n*nsize]; /* vec[n][nsize] */
  double complex *vec_k, *vec_l;
  //double complex mat[n2*n2]; /* mat[n2][n2] */
  double complex *mat; /* mat[n2][n2] */
  double complex *mat_k;
  double complex *invMat; /* mat[n2][n2] */
  //double complex matUV[n2*nsize]; /* mat[n2][nsize] */
  double complex *matUV; /* mat[n2][nsize] */
  double complex *invMat_k, *matUV_i, *matUV_j;
  double sgn;

  int rsi,rsk,msi,msj,k,l;
  double complex val,tmp;

  /* for ZSKPFA */
  char uplo='U', mthd='P';
  int m,nn,lda,info=0;
  double complex pfaff;
  int iwork[n2];
  int iwork2[n2];
  //double complex work[n2*n2]; /* [n2][n2] */
  double complex *work; /* [n2][n2] */
  double rwork[n2*n2]; /* [n2][n2] */
  int lwork = n2*n2;
  nn=lda=n2;

  sltE = bufM + qpidx*Nsite2*Nsite2;
  invM = InvM + qpidx*Nsize*Nsize;

  vec = (double complex*)malloc(sizeof(double complex)*n*nsize);
  //invMat = (double *)malloc(sizeof(double)*n2*n2);
  mat = (double complex*)malloc(sizeof(double complex)*n2*n2);
  //matUV = (double *)malloc(sizeof(double)*n2*nsize);
  work = (double complex*)malloc(sizeof(double complex)*n2*n2);

  //#pragma loop noalias
  for(k=0;k<n;k++) {
    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
    sltE_k = sltE + rsk*Nsite2;
    vec_k = vec + k*nsize;
    //#pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) {
      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
      vec_k[msi] = sltE_k[rsi];
    }
  }

  /* X_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k;
    vec_k = vec + k*nsize;
    for(l=k+1;l<n;l++) {
      vec_l = vec + l*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        invM_i = invM + msi*nsize;
        tmp = 0.0;
        for(msj=0;msj<nsize;msj++) {
          tmp += invM_i[msj] * vec_l[msj];
        }
        val += tmp * vec_k[msi];
      }
      mat_k[l] = val + vec_k[msa[l]];
    }
  }

  /* Y_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k + n;
    vec_k = vec + k*nsize;
    for(l=0;l<n;l++) {
      invM_l = invM + msa[l]*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        val += vec_k[msi] * invM_l[msi];
      }
      mat_k[l] = val;
    }
  }

  /* Z_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*(k+n) + n;
    invM_k = invM + msa[k]*nsize;
    for(l=k+1;l<n;l++) {
      mat_k[l] = invM_k[msa[l]];
    }
  }

  //#pragma loop noalias
  for(k=0;k<n2;k++) {
    //#pragma loop norecurrence
    for(l=0;l<k;l++) {
      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
    }
    mat[n2*k + k] = 0.0; /* diagonal elements */
  }

  /* calculate Pf M */
  //TODO: CHECK rwork is real
  M_ZSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
  //M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);

  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;

  //free(matUV);
  free(mat);
  free(vec);
  //free(invMat);
  free(work);
  return sgn * pfaff * PfM[qpidx];
}

void UpdateMAll_BF_fcmp(const int *icount, const int *msaTmp,
                        double complex* pfMNew, const int *eleIdx,
                        const int qpStart, const int qpEnd)
{
#pragma procedure serial
  const int nsize = Nsize;
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  //double complex *sltE;
  //double complex *sltE_i;
  int *msa;
  int msi,msj,rsi,rsj,i;
  //int *hop;
  //double complex diff;

  for(qpidx=0;qpidx<qpNum;qpidx++) {
    //Store msa//
    msa=(int *)malloc(sizeof(int)*icount[qpidx]);
    for(i=0;i<icount[qpidx];i++){
      msa[i] = msaTmp[i+qpidx*Nsite];
    }

    /* calculateNewPfM */
    pfMNew[qpidx] = updateMAll_BF_fcmp_child(qpidx,icount[qpidx],msa,eleIdx);

    free(msa);
  }

  return;
}

/* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
/* buffer size = n*Nsize */
//void updateMAllBF3_child(const int qpidx, const int n, const int *msa,
double complex updateMAll_BF_fcmp_child(
        const int qpidx,
        const int n, const int *msa,
        const int *eleIdx)
{
  const int nsize = Nsize;
  const int n2 = 2*n;
  const double complex*sltE;
  const double complex*sltE_k;
  double complex*invM;
  double complex *invM_i, *invM_k, *invM_l;

  double complex*vec; /* vec[n][nsize] */
  double complex*vec_k, *vec_l;
  //double complex mat[n2*n2]; /* mat[n2][n2] */
  double complex*mat; /* mat[n2][n2] */
  double complex*mat_k;
  //double complex invMat[n2*n2]; /* mat[n2][n2] */
  double complex*invMat; /* mat[n2][n2] */
  //double complex matUV[n2*nsize]; /* mat[n2][nsize] */
  double complex*matUV; /* mat[n2][nsize] */
  double complex*invMat_k, *matUV_i, *matUV_j;
  double complex sgn;

  int rsi,rsk,msi,msj,k,l;
  double complex val,tmp;

  /* for DSKPFA */
  char uplo='U', mthd='P';
  int m,nn,lda,info=0;
  double complex pfaff;
  int iwork[n2];
  int iwork2[n2];
  //double complex work[n2*n2]; /* [n2][n2] */
  double rwork[n2*n2]; /* [n2][n2] */
  double complex*work; /* [n2][n2] */
  double complex*work2; /* [n2][n2] */
  int lwork = n2*n2;
  int lwork2 = n2;
  m=nn=lda=n2;

  sltE = SlaterElmBF + qpidx*Nsite2*Nsite2;
  invM = InvM + qpidx*Nsize*Nsize;

  //vec = bufferc; /* n*nsize */
  vec = (double complex*)malloc(sizeof(double complex)*n*nsize);
  invMat = (double complex*)malloc(sizeof(double complex)*n2*n2);
  mat = (double complex*)malloc(sizeof(double complex)*n2*n2);
  matUV = (double complex*)malloc(sizeof(double complex)*n2*nsize);
  work = (double complex*)malloc(sizeof(double complex)*n2*n2);
  work2 = (double complex*)malloc(sizeof(double complex)*n2);

#pragma loop noalias
  for(k=0;k<n;k++) {
    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
    sltE_k = sltE + rsk*Nsite2;
    vec_k = vec + k*nsize;
#pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) {
      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
      vec_k[msi] = sltE_k[rsi];
    }
  }

  /* X_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k;
    vec_k = vec + k*nsize;
    for(l=k+1;l<n;l++) {
      vec_l = vec + l*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        invM_i = invM + msi*nsize;
        tmp = 0.0;
        for(msj=0;msj<nsize;msj++) {
          tmp += invM_i[msj] * vec_l[msj];
        }
        val += tmp * vec_k[msi];
      }
      mat_k[l] = val + vec_k[msa[l]];
    }
  }

  /* Y_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k + n;
    vec_k = vec + k*nsize;
    for(l=0;l<n;l++) {
      invM_l = invM + msa[l]*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        val += vec_k[msi] * invM_l[msi];
      }
      mat_k[l] = val;
    }
  }

  /* Z_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*(k+n) + n;
    invM_k = invM + msa[k]*nsize;
    for(l=k+1;l<n;l++) {
      mat_k[l] = invM_k[msa[l]];
    }
  }

#pragma loop noalias
  for(k=0;k<n2;k++) {
#pragma loop norecurrence
    for(l=0;l<k;l++) {
      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
    }
    mat[n2*k + k] = 0.0; /* diagonal elements */
  }

  //#pragma loop noalias
  for(k=0;k<n2;k++) {
    //#pragma loop norecurrence
    for(l=0;l<n2;l++) {
      invMat[n2*k + l] = mat[n2*l + k];
    }
  }

  /* matUV_ik */
  for(msi=0;msi<nsize;msi++) {
    matUV_i = matUV + n2*msi;
    invM_i = invM + msi*nsize;
    for(k=0;k<n;k++) {
      vec_k = vec + k*nsize;
      val = 0.0;
      for(msj=0;msj<nsize;msj++) {
        val -= invM_i[msj] * vec_k[msj];
      }
      if(msa[k]==msi){val -= 1.0;}
      matUV_i[k] = val;
      matUV_i[k+n] = invM_i[msa[k]];
    }
  }
  /* DInv */
  m=nn=lda=n2;
  M_ZGETRF(&m, &nn, invMat, &lda, iwork2, &info); /* ipiv = iwork */
  //M_DGETRF(&m, &nn, invMat, &lda, iwork2, &info); /* ipiv = iwork */
  if(info!=0){
    printf("Error M_ZGETRF %d\n",info);
    //return;
  }

  M_ZGETRI(&nn, invMat, &lda, iwork2, work2, &lwork2, &info);
  //M_DGETRI(&nn, invMat, &lda, iwork2, work2, &lwork2, &info);
  if(info!=0){
    printf("Error M_ZGETRI %d\n",info);
    //return;
  }

  for(msi=0;msi<nsize;msi++) {
    matUV_i = matUV + n2*msi;
    invM_i = invM + msi*nsize;
    for(msj=0;msj<nsize;msj++) {
      matUV_j = matUV + n2*msj;
      val = 0.0;
      for(k=0;k<n2;k++) {
        invMat_k = invMat + k*n2;
        tmp=0.0;
        for(l=0;l<n2;l++) {
          tmp += invMat_k[l]*matUV_j[l];
        }
        val += matUV_i[k]*tmp;
      }
      invM_i[msj] -= val;
    }
  }

  for(msi=0;msi<Ne;msi++) {
    //invM_i = invM + msi*nsize;
    for(msj=Ne;msj<Nsize;msj++) {
      InvM[msj*nsize+msi] = -InvM[msi*nsize+msj];
    }
  }
  for(msi=Ne;msi<Nsize;msi++) {
    for(msj=Ne;msj<Nsize;msj++) {
      InvM[msi*nsize+msj] = 0.0;
    }
  }
  for(msi=0;msi<Ne;msi++) {
    for(msj=0;msj<Ne;msj++) {
      InvM[msi*nsize+msj] = 0.0;
    }
  }

  /* calculate Pf M */
  M_ZSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
  //M_DSKPFA(&uplo, &mthd, &nn, invM, &lda, &pfaff, iwork, work, &lwork, &info);
  //M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);

  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;

  free(matUV);
  free(mat);
  free(vec);
  free(invMat);
  free(work);
  free(work2);

  return sgn * pfaff * PfM[qpidx];
}


#endif
