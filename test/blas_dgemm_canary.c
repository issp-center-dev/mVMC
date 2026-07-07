#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern void dgemm_(char *, char *, int *, int *, int *, const double *,
                   const double *, int *, const double *, int *,
                   const double *, double *, int *);

static double a_value(int row, int col) {
  return 0.125 + 0.03125 * (double)(row + 1) - 0.0175 * (double)(col + 1);
}

static double b_value(int row, int col) {
  return -0.0625 + 0.02125 * (double)(row + 2) + 0.01375 * (double)(col + 1);
}

static int check_tn_shape(int m, int n, int k) {
  int lda = k;
  int ldb = k;
  int ldc = m;
  const double alpha = 1.0;
  const double beta = 0.0;
  char transa = 'T';
  char transb = 'N';
  double *a = (double *)malloc(sizeof(double) * (size_t)lda * (size_t)m);
  double *b = (double *)malloc(sizeof(double) * (size_t)ldb * (size_t)n);
  double *c = (double *)malloc(sizeof(double) * (size_t)ldc * (size_t)n);
  double maxdiff = 0.0;
  int i, j, l;

  if(a == NULL || b == NULL || c == NULL) {
    fprintf(stderr, "BLAS_DGEMM_Canary allocation failed for M=%d N=%d K=%d\n", m, n, k);
    free(a);
    free(b);
    free(c);
    return 1;
  }

  for(j=0;j<m;j++) {
    for(i=0;i<k;i++) {
      a[i + j*lda] = a_value(i, j);
    }
  }
  for(j=0;j<n;j++) {
    for(i=0;i<k;i++) {
      b[i + j*ldb] = b_value(i, j);
    }
  }
  for(j=0;j<n;j++) {
    for(i=0;i<m;i++) {
      c[i + j*ldc] = -12345.0;
    }
  }

  dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);

  for(j=0;j<n;j++) {
    for(i=0;i<m;i++) {
      double ref = 0.0;
      double diff;
      for(l=0;l<k;l++) {
        ref += a[l + i*lda] * b[l + j*ldb];
      }
      diff = fabs(c[i + j*ldc] - ref);
      if(diff > maxdiff) maxdiff = diff;
    }
  }

  free(a);
  free(b);
  free(c);

  if(!isfinite(maxdiff) || maxdiff > 1.0e-11) {
    fprintf(stderr, "BLAS_DGEMM_Canary FAIL M=%d N=%d K=%d maxdiff=%.17e\n", m, n, k, maxdiff);
    return 1;
  }

  return 0;
}

int main(void) {
  const int ms[] = {8, 16, 32};
  const int ns[] = {4, 16};
  const int ks[] = {2, 8, 32};
  int im, in, ik;

  for(im=0;im<(int)(sizeof(ms)/sizeof(ms[0]));im++) {
    for(in=0;in<(int)(sizeof(ns)/sizeof(ns[0]));in++) {
      for(ik=0;ik<(int)(sizeof(ks)/sizeof(ks[0]));ik++) {
        if(check_tn_shape(ms[im], ns[in], ks[ik]) != 0) return 1;
      }
    }
  }

  printf("BLAS_DGEMM_Canary PASS\n");
  return 0;
}
