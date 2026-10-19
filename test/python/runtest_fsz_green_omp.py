"""Mixed spin indices must be private in the non-BF FSZ Green loops."""
import math
import os
from pathlib import Path

from runtest_bf_fsz import (prepare_case, update_modpara,
    write_spin_changing_c2_defs, write_two_body_rows, run_vmc, read_float_rows)


def main():
    root = os.getcwd()
    snapshots = {}
    for threads in (1, 4):
        work = prepare_case(root, "FSZ_Green_SpinIndices_omp{}".format(threads), False)
        update_modpara(work, {"2Sz": "-1", "NVMCSample": "16"})
        terms = write_spin_changing_c2_defs(work)
        # Interleave spin classes rather than grouping equal spins in chunks.
        # Repeated observables also let the test check scheduling independence
        # within one run, without statistical tolerances.
        write_two_body_rows(work, terms * 32)
        one = Path(work) / "greenone.def"
        lines = one.read_text().splitlines()
        terms1 = ["{} {} {} {}".format(i % 4, i // 4, j % 4, j // 4)
                  for i in range(8) for j in range(8)]
        lines[1] = "NCisAjs {}".format(len(terms1))
        one.write_text("\n".join(lines[:5] + terms1) + "\n")
        proc = run_vmc(root, work, extra_env={"OMP_NUM_THREADS": str(threads),
            "OMP_NESTED": "true", "OMP_MAX_ACTIVE_LEVELS": "2"})
        if proc.returncode:
            raise RuntimeError(proc.stdout)
        values = {}
        for name, period in (("zvo_cisajs_001.dat", len(terms1)),
                             ("zvo_cisajscktalt_001.dat", len(terms))):
            rows = read_float_rows(str(Path(work) / "output" / name))
            assert len(rows) == period * (32 if name == "zvo_cisajscktalt_001.dat" else 1), (name, len(rows))
            assert all(math.isfinite(v) for row in rows for v in row), name
            for i, row in enumerate(rows):
                reference = rows[i % period]
                assert row[:-2] == reference[:-2], (name, i)
                assert max(abs(a-b) for a, b in zip(row[-2:], reference[-2:])) <= 1e-12, (threads, name, i, "duplicate operator changed")
                if name == "zvo_cisajscktalt_001.dat" and i % period in (2, 5, 6, 7, 12):
                    assert max(abs(v) for v in row[-2:]) <= 1e-12, (threads, i, "structural zero")
            values[name] = rows
        snapshots[threads] = values
    for name in snapshots[1]:
        for a, b in zip(snapshots[1][name], snapshots[4][name]):
            assert a[:-2] == b[:-2]
            assert max(abs(x-y) for x, y in zip(a[-2:], b[-2:])) <= 1e-11, (name, "OMP mismatch")
    print("FSZ mixed-spin Green OMP=1/4: duplicate operators, structural zeros, and outputs agree")


if __name__ == "__main__":
    main()
