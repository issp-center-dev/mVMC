#include <limits.h>
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
  Lanczos2Contract contract;
  memset(&contract, 0, sizeof(contract));
  contract.step = 2;
  contract.lanczosMode = 1;
  contract.vmcCalMode = 1;
  contract.nsite = 4;
  contract.ne = 2;
  contract.nTransfer = 8;
  contract.nQPFull = 1;
  return contract;
}

static Lanczos2Contract SupportedHeisenbergContract(void) {
  Lanczos2Contract contract = SupportedContract();
  contract.exUpdatePath = 2;
  contract.nExchangeCoupling = 4;
  contract.nLocSpn = 4;
  contract.nTransfer = 0;
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

static void CheckLegacyExecutionStatus(const Lanczos2Contract *contract,
                                       Lanczos2ContractStatus expected,
                                       const char *label) {
  const Lanczos2ContractStatus actual =
      ValidateLegacyLanczos2ExecutionContract(contract);
  CHECK(actual == expected, "%s: legacy status=%d expected=%d", label,
        (int)actual, (int)expected);
  if (actual != LANCZOS2_CONTRACT_OK) {
    CHECK(Lanczos2ContractError(actual)[0] != '\0',
          "%s: empty legacy error message", label);
  }
}

static void TestAllowedInputs(void) {
  Lanczos2Contract contract = SupportedContract();
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK, "supported step 2");

  contract = SupportedHeisenbergContract();
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK,
              "supported pure-spin Heisenberg step 2");
  CHECK(ClassifyLanczos2Model(&contract) ==
            LANCZOS2_MODEL_LOCAL_SPIN_EXCHANGE,
        "pure-spin model classification");

  contract.nExchangeCoupling = 0;
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK,
              "supported pure-spin Ising step 2");

  contract = SupportedContract();
  CHECK(ClassifyLanczos2Model(&contract) == LANCZOS2_MODEL_ELECTRONIC_VK,
        "electronic model classification");

  /* Step 1 remains the legacy compatibility path, but not an unconditional OK. */
  contract = SupportedContract();
  contract.step = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK, "step 1 compatibility");

  contract = SupportedContract();
  contract.step = 1;
  contract.lanczosMode = 2;
  contract.nProjBF = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK,
              "step 1 legacy capability compatibility");
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
  contract.twoSz = 2;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_TWO_SZ, "2Sz");

  contract = SupportedContract();
  contract.nProjBF = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_BACKFLOW, "BackFlow");

  contract = SupportedContract();
  contract.flagRBM = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_RBM, "RBM");

  contract = SupportedContract();
  contract.reweight = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_REWEIGHT,
              "reweight");

  contract = SupportedContract();
  contract.exUpdatePath = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH,
              "update path");

  contract = SupportedContract();
  contract.nTransfer = 0;
  CheckStatus(&contract,
              LANCZOS2_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER,
              "electronic Transfer");

  contract = SupportedHeisenbergContract();
  contract.nLocSpn = 3;
  CheckStatus(&contract, LANCZOS2_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN,
              "mixed localized and itinerant model");

  contract = SupportedHeisenbergContract();
  contract.nsite = 5;
  CheckStatus(&contract, LANCZOS2_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN,
              "localized spin does not cover every site");

  contract = SupportedHeisenbergContract();
  contract.nTransfer = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_TRANSFER,
              "pure-spin Transfer");

  contract = SupportedHeisenbergContract();
  contract.nQPFull = 2;
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK,
              "pure-spin multi-QP production contract");

  contract = SupportedHeisenbergContract();
  contract.ne = INT_MAX;
  contract.nsite = INT_MAX;
  contract.nLocSpn = INT_MAX;
  CheckStatus(&contract, LANCZOS2_CONTRACT_SIZE_OVERFLOW,
              "pure-spin electron-count overflow");

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

  {
    const int unsupportedPaths[] = {1, 3, 4, 5, 6};
    size_t i;
    for (i = 0;
         i < sizeof(unsupportedPaths) / sizeof(unsupportedPaths[0]); i++) {
      contract = SupportedHeisenbergContract();
      contract.exUpdatePath = unsupportedPaths[i];
      CheckStatus(&contract, LANCZOS2_CONTRACT_UNSUPPORTED_UPDATE_PATH,
                  "pure-spin unsupported update path");
    }
  }
}

static void TestCommonDomainAndOverflow(void) {
  Lanczos2Contract contract = SupportedContract();
  contract.step = 1;
  contract.ne = 0;
  CheckStatus(&contract, LANCZOS2_CONTRACT_INVALID_COUNT,
              "step 1 zero Ne rejected");

  contract = SupportedContract();
  contract.step = 1;
  contract.nPairHopping = -1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_INVALID_COUNT,
              "step 1 negative count rejected");

  contract = SupportedContract();
  contract.step = 1;
  contract.nsite = INT_MAX;
  CheckStatus(&contract, LANCZOS2_CONTRACT_SIZE_OVERFLOW,
              "step 1 Nsite2 overflow rejected");

  contract = SupportedContract();
  contract.nQPFull = 0;
  CheckStatus(&contract, LANCZOS2_CONTRACT_INVALID_COUNT,
              "step 2 zero NQPFull rejected");

  contract = SupportedContract();
  contract.ne = INT_MAX;
  CheckStatus(&contract, LANCZOS2_CONTRACT_SIZE_OVERFLOW,
              "step 2 Nsize overflow rejected");
}

static void TestFirstFailureIsStable(void) {
  Lanczos2Contract contract = SupportedContract();
  contract.lanczosMode = 0;
  contract.nNBodyG = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1,
              "first-failure order");

  contract = SupportedHeisenbergContract();
  contract.lanczosMode = 2;
  contract.nTransfer = 1;
  CheckStatus(&contract, LANCZOS2_CONTRACT_REQUIRES_LANCZOS_MODE_1,
              "pure-spin first-failure order");
}

static void TestLegacyPureSpinExecutionLimit(void) {
  Lanczos2Contract contract = SupportedHeisenbergContract();

  contract.nQPFull = 2;
  CheckStatus(&contract, LANCZOS2_CONTRACT_OK,
              "pure-spin multi-QP remains a production capability");
  CheckLegacyExecutionStatus(
      &contract, LANCZOS2_CONTRACT_UNSUPPORTED_QUANTUM_PROJECTION,
      "legacy pure-spin multi-QP execution");

  contract.nQPFull = 1;
  CheckLegacyExecutionStatus(&contract, LANCZOS2_CONTRACT_OK,
                             "legacy pure-spin single-QP execution");

  contract = SupportedContract();
  contract.nQPFull = 4;
  CheckLegacyExecutionStatus(&contract, LANCZOS2_CONTRACT_OK,
                             "legacy electronic multi-QP execution");

  contract = SupportedHeisenbergContract();
  contract.step = 1;
  contract.nQPFull = 4;
  CheckLegacyExecutionStatus(&contract, LANCZOS2_CONTRACT_OK,
                             "legacy pure-spin step-one multi-QP execution");

  contract = SupportedHeisenbergContract();
  contract.vmcCalMode = 0;
  contract.nQPFull = 4;
  CheckLegacyExecutionStatus(
      &contract, LANCZOS2_CONTRACT_REQUIRES_VMC_CAL_MODE_1,
      "legacy validation preserves general failure precedence");
}

int main(void) {
  TestAllowedInputs();
  TestStepDomain();
  TestEachUnsupportedCapability();
  TestCommonDomainAndOverflow();
  TestFirstFailureIsStable();
  TestLegacyPureSpinExecutionLimit();

  if (failures != 0) {
    fprintf(stderr, "Lanczos2Contract_Unit: %d failure(s)\n", failures);
    return 1;
  }
  printf("Lanczos2Contract_Unit: PASS\n");
  return 0;
}
