#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "update_weight.h"

static int failures = 0;

static void Expect(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "FAIL: %s\n", message);
    failures++;
  }
}

static FILE *OpenText(const char *text) {
  FILE *fp = tmpfile();
  if (fp == NULL) return NULL;
  fputs(text, fp);
  rewind(fp);
  return fp;
}

static void TestParserAndSelector(void) {
  const char *text =
      "================================\n"
      "NUpdateWeight 3\n"
      "================================\n"
      "========UpdateType Weight=======\n"
      "================================\n"
      "Exchange 2.0\n"
      "LocalSpinFlip 1.0\n"
      "PairSpinFlip 1.0\n";
  char error[256] = "";
  double weights[UpdateWeightCount];
  FILE *fp = OpenText(text);
  Expect(fp != NULL, "tmpfile opens");
  if (fp == NULL) return;
  Expect(ReadUpdateWeight(fp, "unit", weights, error, sizeof(error)) == 0,
         error);
  fclose(fp);
  Expect(fabs(weights[UpdateWeightExchange] - 0.5) < 1.0e-15,
         "exchange weight normalizes");
  Expect(fabs(weights[UpdateWeightLocalSpinFlip] - 0.25) < 1.0e-15,
         "single-flip weight normalizes");
  Expect(fabs(weights[UpdateWeightPairSpinFlip] - 0.25) < 1.0e-15,
         "pair-flip weight normalizes");
  Expect(SelectUpdateWeight(weights, 0.0) == UpdateWeightExchange,
         "lower interval selects exchange");
  Expect(SelectUpdateWeight(weights, 0.5) == UpdateWeightLocalSpinFlip,
         "middle interval selects local spin flip");
  Expect(SelectUpdateWeight(weights, 0.75) == UpdateWeightPairSpinFlip,
         "upper interval selects pair spin flip");
}

static void TestMissingKernelDefaultsToZero(void) {
  const char *text =
      "================================\n"
      "NUpdateWeight 1\n"
      "================================\n"
      "========UpdateType Weight=======\n"
      "================================\n"
      "Exchange 4.0\n";
  char error[256] = "";
  double weights[UpdateWeightCount];
  FILE *fp = OpenText(text);
  Expect(ReadUpdateWeight(fp, "unit", weights, error, sizeof(error)) == 0,
         error);
  fclose(fp);
  Expect(weights[UpdateWeightExchange] == 1.0,
         "single named kernel normalizes to one");
  Expect(weights[UpdateWeightLocalSpinFlip] == 0.0,
         "missing local spin flip defaults to zero");
  Expect(weights[UpdateWeightPairSpinFlip] == 0.0,
         "missing pair spin flip defaults to zero");
}

static void TestParserFailures(void) {
  const char *cases[] = {
      "=\nNUpdateWeight 2\n=\n=\n=\nExchange 1\nExchange 2\n",
      "=\nNUpdateWeight 1\n=\n=\n=\nUnknown 1\n",
      "=\nNUpdateWeight 1\n=\n=\n=\nExchange -1\n",
      "=\nNUpdateWeight 2\n=\n=\n=\nExchange 1\n"};
  int index;
  for (index = 0; index < 4; index++) {
    char error[256] = "";
    double weights[UpdateWeightCount];
    FILE *fp = OpenText(cases[index]);
    Expect(ReadUpdateWeight(fp, "invalid", weights, error, sizeof(error)) != 0,
           "invalid parser case fails");
    Expect(error[0] != '\0', "invalid parser case reports diagnostic");
    fclose(fp);
  }
}

static void TestContract(void) {
  char error[256] = "";
  const double pairWeights[UpdateWeightCount] = {1.0 / 3.0, 1.0 / 3.0,
                                                  1.0 / 3.0};
  const double legacyWeights[UpdateWeightCount] = {0.5, 0.5, 0.0};
  const int allLocal[4] = {1, 1, 1, 1};
  const int mixed[4] = {1, 1, 0, 1};

  Expect(ValidateUpdateWeightContract(1, 2, -1, 4, 4, 2, 1,
                                      pairWeights, error, sizeof(error)) == 0,
         "pure-spin pair-flip contract passes");
  Expect(ValidateUpdateWeightLocSpin(1, 4, 4, allLocal, pairWeights,
                                     error, sizeof(error)) == 0,
         "all-local-spin detail contract passes");
  Expect(ValidateUpdateWeightContract(1, 1, -1, 4, 4, 2, 1,
                                      legacyWeights, error, sizeof(error)) != 0,
         "weight file rejects non-path-2 mode");
  Expect(ValidateUpdateWeightContract(1, 2, 0, 4, 4, 2, 1,
                                      pairWeights, error, sizeof(error)) != 0,
         "pair flip rejects fixed Sz");
  Expect(ValidateUpdateWeightLocSpin(1, 4, 4, mixed, pairWeights,
                                     error, sizeof(error)) != 0,
         "pair flip rejects non-local site");
  Expect(ValidateUpdateWeightContract(0, 0, 0, 0, 0, 0, 0, NULL,
                                      error, sizeof(error)) == 0,
         "disabled contract is a no-op");
}

int main(void) {
  TestParserAndSelector();
  TestMissingKernelDefaultsToZero();
  TestParserFailures();
  TestContract();
  if (failures != 0) {
    fprintf(stderr, "%d UpdateWeight unit checks failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("UpdateWeight unit checks: PASS\n");
  return EXIT_SUCCESS;
}
