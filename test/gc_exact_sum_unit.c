#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gc_exact_sum.h"

static int failures;
static double complex make_complex(double re, double im) {
  double parts[2] = {re, im};
  double complex value;
  memcpy(&value, parts, sizeof value);
  return value;
}

static void check(const char *name, const double complex *v, int n,
                  int expected) {
  int zero = -1;
  int status = GCExactComplexSumIsZero(v, n, &zero);
  if (status != 0 || zero != expected) {
    fprintf(stderr, "%s: status=%d zero=%d expected=%d\n", name, status,
            zero, expected);
    ++failures;
  }
}

static void permutations(double complex *v, int start, int n, int expected) {
  int i;
  if (start == n) {
    check("permutation", v, n, expected);
    return;
  }
  for (i = start; i < n; ++i) {
    double complex t = v[start];
    v[start] = v[i]; v[i] = t;
    permutations(v, start + 1, n, expected);
    t = v[start]; v[start] = v[i]; v[i] = t;
  }
}

int main(void) {
  double complex phases[] = {0.25, -0.25 * I, -0.25, 0.25 * I};
  double complex cancellation[] = {1.0, 0x1p-53, -1.0};
  double complex extreme[] = {DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX,
                              0x1p-1074};
  double complex signedZero[] = {0.0, make_complex(-0.0, -0.0)};
  double complex invalid[1];
  int zero, i;
  if (GCExactComplexSumSelfTest() != 0) ++failures;
  permutations(phases, 0, 4, 1);
  permutations(cancellation, 0, 3, 0);
  permutations(extreme, 0, 5, 0);
  extreme[4] = 0.0;
  permutations(extreme, 0, 5, 1);
  check("signed zeros", signedZero, 2, 1);
  for (i = -1074; i <= 1023; i += 17) {
    double x = scalbn(1.0, i);
    double complex v[] = {make_complex(x, -x), make_complex(-x, x)};
    check("exact cancellation across exponent range", v, 2, 1);
    v[1] = make_complex(-x, nextafter(x, INFINITY));
    check("one ulp imaginary perturbation", v, 2, 0);
  }
  for (i = 0; i < 6; ++i) {
    const double nonfinite[] = {NAN, INFINITY, -INFINITY};
    invalid[0] = i < 3 ? make_complex(nonfinite[i], 0.0)
                      : make_complex(0.0, nonfinite[i - 3]);
    zero = 7;
    if (GCExactComplexSumIsZero(invalid, 1, &zero) == 0 || zero != 0)
      ++failures;
  }
  zero = 7;
  if (GCExactComplexSumIsZero(NULL, 1, &zero) == 0 || zero != 0)
    ++failures;
  if (GCExactComplexSumIsZero(phases, 0, &zero) == 0 || zero != 0)
    ++failures;
  if (GCExactComplexSumIsZero(phases, -1, &zero) == 0 || zero != 0)
    ++failures;
  if (GCExactComplexSumIsZero(phases, 4, NULL) == 0) ++failures;
  if (failures) {
    fprintf(stderr, "GC exact-sum: %d failures\n", failures);
    return 1;
  }
  puts("GC exact-sum: PASS (permutations, extremes, subnormal, nonfinite)");
  return 0;
}
