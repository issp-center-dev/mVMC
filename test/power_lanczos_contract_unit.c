#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "power_lanczos_contract.h"

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "PowerLanczosContract_Unit FAIL: ");                     \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                    \
      failures++;                                                               \
    }                                                                           \
  } while (0)

static MVMCPowerLanczosContract ElectronicContract(int order) {
  MVMCPowerLanczosContract contract;
  memset(&contract, 0, sizeof(contract));
  contract.order = order;
  contract.lanczosMode = 1;
  contract.vmcCalMode = 1;
  contract.nsite = 6;
  contract.ne = 3;
  contract.nTransfer = 24;
  contract.nQPFull = 4;
  return contract;
}

static MVMCPowerLanczosContract PureSpinContract(int order) {
  MVMCPowerLanczosContract contract = ElectronicContract(order);
  contract.exUpdatePath = 2;
  contract.nsite = 4;
  contract.ne = 2;
  contract.nLocSpn = 4;
  contract.nTransfer = 0;
  contract.nExchangeCoupling = 4;
  contract.nQPFull = 8;
  return contract;
}

static void CheckStatus(const MVMCPowerLanczosContract *contract,
                        MVMCPowerLanczosContractStatus expected,
                        const char *label) {
  const MVMCPowerLanczosContractStatus actual =
      mvmc_power_lanczos_validate_production_contract(contract);
  CHECK(actual == expected, "%s: status=%d expected=%d", label, (int)actual,
        (int)expected);
  if (actual != MVMC_POWER_LANCZOS_CONTRACT_OK) {
    CHECK(mvmc_power_lanczos_contract_error(actual)[0] != '\0',
          "%s: empty error message", label);
  }
}

static void TestAllowedScope(void) {
  MVMCPowerLanczosDerivedSizes sizes;
  MVMCPowerLanczosContract contract = ElectronicContract(1);

  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_OK,
              "electronic order 1");
  CHECK(mvmc_power_lanczos_classify_model(&contract) ==
            MVMC_POWER_LANCZOS_MODEL_ELECTRONIC_VK,
        "electronic classification");
  CHECK(mvmc_power_lanczos_compute_derived_sizes(&contract, &sizes) ==
            MVMC_POWER_LANCZOS_CONTRACT_OK,
        "derived-size computation");
  CHECK(sizes.nsite2 == 12u && sizes.nsize == 6u &&
            sizes.coefficient_basis_dimension == 2u &&
            sizes.anchor_count == 3u,
        "order 1 derived sizes");

  contract = ElectronicContract(2);
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_OK,
              "electronic order 2");

  contract = PureSpinContract(1);
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_OK,
              "pure spin order 1 multi-QP");
  CHECK(mvmc_power_lanczos_classify_model(&contract) ==
            MVMC_POWER_LANCZOS_MODEL_LOCAL_SPIN_EXCHANGE,
        "pure-spin classification");

  contract = PureSpinContract(2);
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_OK,
              "pure spin order 2 multi-QP");
}

static void TestCommonDomainAndOverflow(void) {
  MVMCPowerLanczosContract contract = ElectronicContract(1);
  CheckStatus(NULL, MVMC_POWER_LANCZOS_CONTRACT_INVALID_ARGUMENT,
              "NULL contract");

  contract.order = 0;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_INVALID_ORDER,
              "invalid low order");
  contract = ElectronicContract(3);
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_INVALID_ORDER,
              "invalid high order");

  contract = ElectronicContract(1);
  contract.nsite = 0;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_INVALID_COUNT,
              "zero site count");

  contract = ElectronicContract(1);
  contract.ne = -1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_INVALID_COUNT,
              "negative electron count");

  contract = ElectronicContract(1);
  contract.nQPFull = 0;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_INVALID_COUNT,
              "zero NQPFull");

  contract = ElectronicContract(1);
  contract.nPairHopping = -1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_INVALID_COUNT,
              "negative capability count");

  contract = ElectronicContract(1);
  contract.ne = INT_MAX / 2 + 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_SIZE_OVERFLOW,
              "Nsize int overflow");

  contract = ElectronicContract(1);
  contract.nsite = INT_MAX / 2 + 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_SIZE_OVERFLOW,
              "Nsite2 int overflow");
}

static void TestRejectedScope(void) {
  MVMCPowerLanczosContract contract = ElectronicContract(1);

  contract.lanczosMode = 2;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_LANCZOS_MODE_1,
              "legacy observable mode");

  contract = ElectronicContract(1);
  contract.vmcCalMode = 0;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_VMC_CAL_MODE_1,
              "VMC calculation mode");

  contract = ElectronicContract(1);
  contract.orbitalGeneral = 1;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_ORBITAL_GENERAL,
              "FSZ/orbital-general");

  contract = ElectronicContract(1);
  contract.twoSz = 2;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TWO_SZ,
              "nonzero 2Sz");

  contract = ElectronicContract(1);
  contract.nProjBF = 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_BACKFLOW,
              "BackFlow");

  contract = ElectronicContract(1);
  contract.flagRBM = 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_RBM,
              "RBM");

  contract = ElectronicContract(1);
  contract.reweight = 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_REWEIGHT,
              "reweight");

  contract = ElectronicContract(1);
  contract.exUpdatePath = 4;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_UPDATE_PATH,
              "t-J update");

  contract = ElectronicContract(1);
  contract.nTransfer = 0;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_SPIN_CONSERVING_TRANSFER,
              "electronic transfer");

  contract = ElectronicContract(1);
  contract.nPairHopping = 1;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_PAIR_HOPPING,
              "PairHopping");

  contract = ElectronicContract(1);
  contract.nExchangeCoupling = 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_EXCHANGE,
              "electronic Exchange");

  contract = ElectronicContract(1);
  contract.nInterAll = 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_INTER_ALL,
              "InterAll");

  contract = ElectronicContract(1);
  contract.nNBodyInterAll = 1;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_INTER_ALL,
              "NBodyInterAll");

  contract = ElectronicContract(1);
  contract.nNBodyG = 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_NBODY_G,
              "NBodyG");

  contract = ElectronicContract(1);
  contract.nSpinFlipTransfer = 1;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_SPIN_FLIP_TRANSFER,
              "spin-flip Transfer");

  contract = PureSpinContract(1);
  contract.nTransfer = 1;
  CheckStatus(&contract, MVMC_POWER_LANCZOS_CONTRACT_UNSUPPORTED_TRANSFER,
              "pure-spin Transfer");

  contract = PureSpinContract(1);
  contract.nLocSpn = 3;
  CheckStatus(&contract,
              MVMC_POWER_LANCZOS_CONTRACT_REQUIRES_PURE_LOCALIZED_SPIN,
              "mixed localized spin");
}

int main(void) {
  CHECK(strcmp(mvmc_power_lanczos_contract_id(),
               MVMC_POWER_LANCZOS_PRODUCTION_CONTRACT_ID) == 0,
        "contract id mismatch");
  TestAllowedScope();
  TestCommonDomainAndOverflow();
  TestRejectedScope();

  if (failures != 0) {
    fprintf(stderr, "PowerLanczosContract_Unit: %d failure(s)\n", failures);
    return 1;
  }
  printf("PowerLanczosContract_Unit: PASS\n");
  return 0;
}
