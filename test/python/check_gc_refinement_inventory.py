#!/usr/bin/env python3
"""Fail CI if the exact-input backend silently omits its regression tests."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inventory", type=Path)
parser.add_argument("--mpi", choices=("on", "off"), required=True)
args = parser.parse_args()
expected = {"GCLanczos_Refinement_Unit", "GCLanczos_Refinement_FinalSteps"}
expected |= {"GCLanczos_FiniteSpin_" + name for name in
             ("Oracle", "Convergence", "Replay", "Initial", "Fault", "Refined", "Contract", "Sample")}
expected |= {"GCLanczos_FiniteSpin_Statistics_" + name for name in
             ("general", "triplet", "maximal", "multi", "pbc_k1", "apbc_k1")}
if args.mpi == "on":
    expected |= {"GCLanczos_FiniteSpin_" + name + "_mpi"
                 for name in ("Layouts", "Fault", "Refined")}
    expected |= {f"GCLanczos_Refinement_Unit_np{n}_mpi" for n in (2, 3)}
    expected |= {f"GCLanczos_Refinement_FinalSteps_np{n}_mpi"
                 for n in (2, 3, 4)}
tests = json.loads(args.inventory.read_text())["tests"]
actual = {t["name"]: t for t in tests
          if t["name"].startswith(("GCLanczos_Refinement_", "GCLanczos_FiniteSpin_"))}
if set(actual) != expected:
    raise SystemExit(f"Refinement inventory mismatch: {sorted(actual)}; "
                     f"expected {sorted(expected)}")
for name, test in actual.items():
    props = {p["name"]: p["value"] for p in test.get("properties", [])}
    if props.get("DISABLED", False) or not test.get("command"):
        raise SystemExit(f"Refinement test is disabled or not executable: {name}")
    if name.startswith("GCLanczos_FiniteSpin_"):
        labels = props.get("LABELS", [])
        command = test["command"]
        if name.endswith("_mpi"):
            if "requires_mpi" not in labels or props.get("PROCESSORS") != 4:
                raise SystemExit(f"Missing MPI scheduling properties: {name}")
        elif name.endswith(("_Fault", "_Refined")) and "--no-mpi" not in command:
            raise SystemExit(f"Serial test may launch multiple MPI ranks: {name}")
        if "_Statistics_" in name and not (
                "long" in labels and "gc_finite_spin_statistics" in labels and
                props.get("RUN_SERIAL") and props.get("TIMEOUT", 0) >= 1800):
            raise SystemExit(f"Missing statistics scheduling properties: {name}")
print(f"Refinement inventory: {len(actual)} enabled tests")
