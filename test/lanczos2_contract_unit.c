#include <stdio.h>
#include <string.h>

#include "lanczos2_contract.h"

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "Lanczos2Contract_Unit FAIL: ");                         \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                    \
      failures++;                                                               \
    }                                                                           \
  } while (0)

static Lanczos2Contract SupportedContract(void) {
  Lanczos2Contract contract = {
      2, /* step */
      1, /* lanczosMode */
      1, /* vmcCalMode */
      0, /* orbitalGeneral */
      0, /* nProjBF */
      0, /* flagRBM */
      0, /* exUpdatePath */
      0, /* nPairHopping */
      0, /* nExchangeCoupling */
      0, /* nInterAll */
      0, /* nNBodyInterAll */
      0, /* nNBodyG */
      0  /* nSpinFlipTransfer */
  };
  return contract;
}

static void CheckStatus(const Lanczos2Contract *contract,
                        Lanczos2ContractStatus expected,
                        const char *label) {
  const Lanczos2ContractStatus actual = ValidateLanczos2Contract(contract);
  CHECK(actual == expected, "%s: status=%d expected=%d", label, (int)actual,
        (int)expected);
  if (actual != LANCZOS2_CONTRACT_OK) {
    CHECK(Lanczos2ContractError(actual)[0] != '\0',
          "%s: empty error message", label);
  }
}

static void TestAllowedInputs(void) {
  Lanczos2Contract contract = SupportedContract();
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK, "supported step 2");

  /*
   * Step 1 is the compatibility path: the new step-2 capability gate must not
   * reject combinations that remain owned by the existing Lanczos code.
   */
  memset(&contract, 0x5a, sizeof(contract));
  contract.step = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK, "step 1 compatibility");
}

static void TestStepDomain(void) {
  const int badSteps[] = {-2, -1, 0, 3, 17};
  size_t i;
  for (i = 0; i < sizeof(badSteps) / sizeof(badSteps[0]); i++) {
    Lanczos2Contract contract = SupportedContract();
    contract.step = badSteps[i];
    CheckStatus(&contract, LANCZOS2_CONTRACT_INVALID_STEP,
                "invalid step domain");
  }
  CheckStatus(NULL, LANCZOS2_CONTRACT_INVALID_STEP, "NULL contract");
}

static void TestEachUnsupportedCapability(void) {
  Lanczos2Contract contract;

  contract = SupportedContract();
  contract.lanczosMode = 0;
  CheckStatus(&contract, LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1,
              "Lanczos mode");

  contract = SupportedContract();
  contract.vmcCalMode = 0;
  CheckStatus(&contract, LANCZOS2_CONTRACT_REQUIRES_VMC_CAL_MODE_1,
              "VMC calculation mode");

  contract = SupportedContract();
  contract.orbitalGeneral = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL,
              "orbital-general");

  contract = SupportedContract();
  contract.nProjBF = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_BACKFLOW, "BackFlow");

  contract = SupportedContract();
  contract.flagRBM = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_RBM, "RBM");

  contract = SupportedContract();
  contract.exUpdatePath = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH,
              "update path");

  contract = SupportedContract();
  contract.nPairHopping = 2;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_PAIR_HOPPING,
              "PairHopping");

  contract = SupportedContract();
  contract.nExchangeCoupling = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_EXCHANGE, "Exchange");

  contract = SupportedContract();
  contract.nInterAll = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_INTER_ALL, "InterAll");

  contract = SupportedContract();
  contract.nNBodyInterAll = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL,
              "NBodyInterAll");

  contract = SupportedContract();
  contract.nNBodyG = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_NBODY_G, "NBodyG");

  contract = SupportedContract();
  contract.nSpinFlipTransfer = 1;
  CheckStatus(
      &contract, LANCZOS2_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER,
      "spin-flip Transfer");
}

static void TestFirstFailureIsStable(void) {
  Lanczos2Contract contract = SupportedContract();
  contract.lanczosMode = 0;
  contract.nNBodyG = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1,
              "first-failure order");
}

int main(void) {
  TestAllowedInputs();
  TestStepDomain();
  TestEachUnsupportedCapability();
  TestFirstFailureIsStable();

  if (failures != 0) {
    fprintf(stderr, "Lanczos2Contract_Unit: %d failure(s)\n", failures);
    return 1;
  }
  printf("Lanczos2Contract_Unit: PASS\n");
  return 0;
}
