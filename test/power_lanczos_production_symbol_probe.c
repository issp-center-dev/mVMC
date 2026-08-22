int oracle_power_lanczos_forbidden_probe(void) { return 17; }

int mvmc_power_lanczos_symbol_policy_probe(void) {
  return oracle_power_lanczos_forbidden_probe();
}

int main(void) {
  return mvmc_power_lanczos_symbol_policy_probe() == 17 ? 0 : 1;
}
