/*
 * Deterministic oracle tests for the RBM primitives used by the real VMC
 * path.  The oracle below evaluates the RBM amplitude directly from the
 * occupation features and does not use rbmCnt.
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double complex LogWeightRBM(const double complex *rbmCnt);
double complex RBMRatio(const double complex *rbmCntNew,
                        const double complex *rbmCntOld);
double complex LogRBMRatio(const double complex *rbmCntNew,
                           const double complex *rbmCntOld);
void MakeRBMCnt(double complex *rbmCnt, const int *eleNum);
void UpdateRBMCnt(int ri, int rj, int spin, double complex *rbmCntNew,
                  const double complex *rbmCntOld, const int *eleNum);
void RBMDiff(double complex *srOptO, const double complex *rbmCnt,
             const int *eleNum);

/* rbm.c owns these globals through global.h. */
extern int Nsite;
extern int Nsite2;
extern int Nneuron;
extern int NneuronGeneral;
extern int NneuronCharge;
extern int NneuronSpin;
extern int NRBM;
extern int NRBM_HiddenLayerIdx;
extern int NRBM_PhysLayerIdx;
extern int NRBM_PhysHiddenIdx;
extern int NGeneralRBM_HiddenLayerIdx;
extern int *GeneralRBM_HiddenLayerIdx;
extern int NGeneralRBM_PhysLayerIdx;
extern int *GeneralRBM_PhysLayerIdx;
extern int NGeneralRBM_PhysHiddenIdx;
extern int **GeneralRBM_PhysHiddenIdx;
extern int NChargeRBM_HiddenLayerIdx;
extern int *ChargeRBM_HiddenLayerIdx;
extern int NChargeRBM_PhysLayerIdx;
extern int *ChargeRBM_PhysLayerIdx;
extern int NChargeRBM_PhysHiddenIdx;
extern int **ChargeRBM_PhysHiddenIdx;
extern int NSpinRBM_HiddenLayerIdx;
extern int *SpinRBM_HiddenLayerIdx;
extern int NSpinRBM_PhysLayerIdx;
extern int *SpinRBM_PhysLayerIdx;
extern int NSpinRBM_PhysHiddenIdx;
extern int **SpinRBM_PhysHiddenIdx;
extern int NBlockSize_RBMRatio;
extern double complex *RBM;

enum {
  SITE_COUNT = 2,
  SPIN_ORBITAL_COUNT = 4,
  PHYS_COUNT = 8,
  HIDDEN_COUNT = 3,
  COUPLING_COUNT = 8,
  PARAM_COUNT = PHYS_COUNT + HIDDEN_COUNT + COUPLING_COUNT,
  COUNT_SIZE = PHYS_COUNT + HIDDEN_COUNT
};

static int charge_phys[] = {0, 1};
static int spin_phys[] = {0, 1};
static int general_phys[] = {0, 1, 2, 3};
static int charge_hidden[] = {0};
static int spin_hidden[] = {0};
static int general_hidden[] = {0};
static int charge_coupling_storage[][1] = {{0}, {1}};
static int spin_coupling_storage[][1] = {{0}, {1}};
static int general_coupling_storage[][1] = {{0}, {1}, {2}, {3}};
static int *charge_coupling[] = {
    charge_coupling_storage[0], charge_coupling_storage[1]};
static int *spin_coupling[] = {
    spin_coupling_storage[0], spin_coupling_storage[1]};
static int *general_coupling[] = {
    general_coupling_storage[0], general_coupling_storage[1],
    general_coupling_storage[2], general_coupling_storage[3]};
static double complex parameters[PARAM_COUNT];

static int failures = 0;

static void check_close(const char *label, double complex actual,
                        double complex expected, double tolerance) {
  const double error = cabs(actual - expected);
  if (!isfinite(creal(actual)) || !isfinite(cimag(actual)) ||
      error > tolerance) {
    fprintf(stderr,
            "FAIL %s: actual=(%.17g, %.17g) expected=(%.17g, %.17g) "
            "error=%.3e tolerance=%.3e\n",
            label, creal(actual), cimag(actual), creal(expected),
            cimag(expected), error, tolerance);
    failures++;
  }
}

static void configure_fixture(void) {
  int i;

  Nsite = SITE_COUNT;
  Nsite2 = SPIN_ORBITAL_COUNT;
  NneuronCharge = 1;
  NneuronSpin = 1;
  NneuronGeneral = 1;
  Nneuron = HIDDEN_COUNT;

  NChargeRBM_PhysLayerIdx = 2;
  NSpinRBM_PhysLayerIdx = 2;
  NGeneralRBM_PhysLayerIdx = 4;
  NRBM_PhysLayerIdx = PHYS_COUNT;
  ChargeRBM_PhysLayerIdx = charge_phys;
  SpinRBM_PhysLayerIdx = spin_phys;
  GeneralRBM_PhysLayerIdx = general_phys;

  NChargeRBM_HiddenLayerIdx = 1;
  NSpinRBM_HiddenLayerIdx = 1;
  NGeneralRBM_HiddenLayerIdx = 1;
  NRBM_HiddenLayerIdx = HIDDEN_COUNT;
  ChargeRBM_HiddenLayerIdx = charge_hidden;
  SpinRBM_HiddenLayerIdx = spin_hidden;
  GeneralRBM_HiddenLayerIdx = general_hidden;

  NChargeRBM_PhysHiddenIdx = 2;
  NSpinRBM_PhysHiddenIdx = 2;
  NGeneralRBM_PhysHiddenIdx = 4;
  NRBM_PhysHiddenIdx = COUPLING_COUNT;
  ChargeRBM_PhysHiddenIdx = charge_coupling;
  SpinRBM_PhysHiddenIdx = spin_coupling;
  GeneralRBM_PhysHiddenIdx = general_coupling;

  NRBM = PARAM_COUNT;
  NBlockSize_RBMRatio = 8;
  RBM = parameters;
  for (i = 0; i < PARAM_COUNT; ++i) {
    parameters[i] = -0.19 + 0.03 * i;
  }
}

static void make_features(const int *occupation, int *charge, int *spin,
                          int *general) {
  int site;
  for (site = 0; site < SITE_COUNT; ++site) {
    charge[site] = occupation[site] + occupation[site + SITE_COUNT] - 1;
    spin[site] = occupation[site] - occupation[site + SITE_COUNT];
  }
  for (site = 0; site < SPIN_ORBITAL_COUNT; ++site) {
    general[site] = 2 * occupation[site] - 1;
  }
}

static double complex direct_log_amplitude(const int *occupation) {
  int charge[SITE_COUNT];
  int spin[SITE_COUNT];
  int general[SPIN_ORBITAL_COUNT];
  double complex visible = 0.0;
  double complex theta_charge = parameters[PHYS_COUNT];
  double complex theta_spin = parameters[PHYS_COUNT + 1];
  double complex theta_general = parameters[PHYS_COUNT + 2];
  int i;

  make_features(occupation, charge, spin, general);
  for (i = 0; i < SITE_COUNT; ++i) {
    visible += parameters[i] * charge[i];
    visible += parameters[2 + i] * spin[i];
    theta_charge += parameters[PHYS_COUNT + HIDDEN_COUNT + i] * charge[i];
    theta_spin +=
        parameters[PHYS_COUNT + HIDDEN_COUNT + 2 + i] * spin[i];
  }
  for (i = 0; i < SPIN_ORBITAL_COUNT; ++i) {
    visible += parameters[4 + i] * general[i];
    theta_general +=
        parameters[PHYS_COUNT + HIDDEN_COUNT + 4 + i] * general[i];
  }
  return visible + clog(ccosh(theta_charge)) + clog(ccosh(theta_spin)) +
         clog(ccosh(theta_general));
}

static void check_counts_and_amplitude(const int *occupation,
                                       const char *case_name) {
  double complex count[COUNT_SIZE];
  double complex expected_count[COUNT_SIZE];
  int charge[SITE_COUNT];
  int spin[SITE_COUNT];
  int general[SPIN_ORBITAL_COUNT];
  char label[128];
  int i;

  make_features(occupation, charge, spin, general);
  for (i = 0; i < COUNT_SIZE; ++i) expected_count[i] = 0.0;
  for (i = 0; i < SITE_COUNT; ++i) {
    expected_count[i] = charge[i];
    expected_count[2 + i] = spin[i];
  }
  for (i = 0; i < SPIN_ORBITAL_COUNT; ++i) {
    expected_count[4 + i] = general[i];
  }
  expected_count[PHYS_COUNT] = parameters[PHYS_COUNT];
  expected_count[PHYS_COUNT + 1] = parameters[PHYS_COUNT + 1];
  expected_count[PHYS_COUNT + 2] = parameters[PHYS_COUNT + 2];
  for (i = 0; i < SITE_COUNT; ++i) {
    expected_count[PHYS_COUNT] +=
        parameters[PHYS_COUNT + HIDDEN_COUNT + i] * charge[i];
    expected_count[PHYS_COUNT + 1] +=
        parameters[PHYS_COUNT + HIDDEN_COUNT + 2 + i] * spin[i];
  }
  for (i = 0; i < SPIN_ORBITAL_COUNT; ++i) {
    expected_count[PHYS_COUNT + 2] +=
        parameters[PHYS_COUNT + HIDDEN_COUNT + 4 + i] * general[i];
  }

  MakeRBMCnt(count, occupation);
  for (i = 0; i < COUNT_SIZE; ++i) {
    snprintf(label, sizeof(label), "%s count[%d]", case_name, i);
    check_close(label, count[i], expected_count[i], 2e-14);
  }
  snprintf(label, sizeof(label), "%s log amplitude", case_name);
  check_close(label, LogWeightRBM(count), direct_log_amplitude(occupation),
              3e-14);
}

static void check_all_hop_updates(void) {
  int old_occupation[SPIN_ORBITAL_COUNT];
  int new_occupation[SPIN_ORBITAL_COUNT];
  double complex old_count[COUNT_SIZE];
  double complex incremental_count[COUNT_SIZE];
  double complex rebuilt_count[COUNT_SIZE];
  char label[160];
  int mask, spin_value, from, to, i;

  for (mask = 0; mask < (1 << SPIN_ORBITAL_COUNT); ++mask) {
    for (i = 0; i < SPIN_ORBITAL_COUNT; ++i) {
      old_occupation[i] = (mask >> i) & 1;
    }
    MakeRBMCnt(old_count, old_occupation);
    for (spin_value = 0; spin_value < 2; ++spin_value) {
      for (from = 0; from < SITE_COUNT; ++from) {
        for (to = 0; to < SITE_COUNT; ++to) {
          const int source = from + spin_value * SITE_COUNT;
          const int destination = to + spin_value * SITE_COUNT;
          double complex log_ratio;
          double complex ratio;
          double complex oracle_log_ratio;

          if (from == to || !old_occupation[source] ||
              old_occupation[destination]) {
            continue;
          }
          for (i = 0; i < SPIN_ORBITAL_COUNT; ++i) {
            new_occupation[i] = old_occupation[i];
          }
          new_occupation[source] = 0;
          new_occupation[destination] = 1;
          UpdateRBMCnt(from, to, spin_value, incremental_count, old_count,
                       new_occupation);
          MakeRBMCnt(rebuilt_count, new_occupation);
          for (i = 0; i < COUNT_SIZE; ++i) {
            snprintf(label, sizeof(label),
                     "mask=%d hop=%d:%d->%d count[%d]", mask, spin_value,
                     from, to, i);
            check_close(label, incremental_count[i], rebuilt_count[i], 3e-14);
          }

          oracle_log_ratio = direct_log_amplitude(new_occupation) -
                             direct_log_amplitude(old_occupation);
          log_ratio = LogRBMRatio(incremental_count, old_count);
          ratio = RBMRatio(incremental_count, old_count);
          snprintf(label, sizeof(label), "mask=%d hop=%d:%d->%d log ratio",
                   mask, spin_value, from, to);
          check_close(label, log_ratio, oracle_log_ratio, 8e-14);
          snprintf(label, sizeof(label), "mask=%d hop=%d:%d->%d ratio", mask,
                   spin_value, from, to);
          check_close(label, ratio, cexp(oracle_log_ratio), 8e-14);
          snprintf(label, sizeof(label), "mask=%d hop=%d:%d->%d log imag",
                   mask, spin_value, from, to);
          check_close(label, I * cimag(log_ratio), 0.0, 8e-14);
          snprintf(label, sizeof(label), "mask=%d hop=%d:%d->%d ratio imag",
                   mask, spin_value, from, to);
          check_close(label, I * cimag(ratio), 0.0, 8e-14);
        }
      }
    }
  }
}

static void check_exchange_update(void) {
  const int old_occupation[SPIN_ORBITAL_COUNT] = {1, 0, 0, 1};
  int intermediate[SPIN_ORBITAL_COUNT] = {0, 1, 0, 1};
  const int final_occupation[SPIN_ORBITAL_COUNT] = {0, 1, 1, 0};
  double complex old_count[COUNT_SIZE];
  double complex incremental_count[COUNT_SIZE];
  double complex rebuilt_count[COUNT_SIZE];
  int i;

  MakeRBMCnt(old_count, old_occupation);
  UpdateRBMCnt(0, 1, 0, incremental_count, old_count, intermediate);
  intermediate[3] = 0;
  intermediate[2] = 1;
  UpdateRBMCnt(1, 0, 1, incremental_count, incremental_count, intermediate);
  MakeRBMCnt(rebuilt_count, final_occupation);
  for (i = 0; i < COUNT_SIZE; ++i) {
    char label[80];
    snprintf(label, sizeof(label), "exchange count[%d]", i);
    check_close(label, incremental_count[i], rebuilt_count[i], 3e-14);
  }
  check_close("exchange log ratio",
              LogRBMRatio(incremental_count, old_count),
              direct_log_amplitude(final_occupation) -
                  direct_log_amplitude(old_occupation),
              8e-14);
}

static void check_derivative(void) {
  const int occupation[SPIN_ORBITAL_COUNT] = {1, 0, 0, 1};
  double complex count[COUNT_SIZE];
  double complex derivative[2 * PARAM_COUNT];
  const double epsilon = 1e-6;
  int i;

  MakeRBMCnt(count, occupation);
  RBMDiff(derivative, count, occupation);
  for (i = 0; i < PARAM_COUNT; ++i) {
    const double complex saved = parameters[i];
    double complex plus;
    double complex minus;
    double complex finite_difference;
    char label[96];

    parameters[i] = saved + epsilon;
    plus = direct_log_amplitude(occupation);
    parameters[i] = saved - epsilon;
    minus = direct_log_amplitude(occupation);
    parameters[i] = saved;
    finite_difference = (plus - minus) / (2.0 * epsilon);

    snprintf(label, sizeof(label), "RBMDiff real parameter[%d]", i);
    check_close(label, derivative[2 * i], finite_difference, 2e-9);
    snprintf(label, sizeof(label), "RBMDiff imaginary slot[%d]", i);
    check_close(label, derivative[2 * i + 1], I * derivative[2 * i],
                3e-14);
  }
}

int main(void) {
  int occupation[SPIN_ORBITAL_COUNT];
  int mask, i;

  configure_fixture();
  for (mask = 0; mask < (1 << SPIN_ORBITAL_COUNT); ++mask) {
    char case_name[32];
    for (i = 0; i < SPIN_ORBITAL_COUNT; ++i) {
      occupation[i] = (mask >> i) & 1;
    }
    snprintf(case_name, sizeof(case_name), "occupation-%02d", mask);
    check_counts_and_amplitude(occupation, case_name);
  }
  check_all_hop_updates();
  check_exchange_update();
  check_derivative();

  if (failures != 0) {
    fprintf(stderr, "real RBM oracle: %d failure(s)\n", failures);
    return EXIT_FAILURE;
  }
  printf("real RBM primitive oracle passed\n");
  return EXIT_SUCCESS;
}
