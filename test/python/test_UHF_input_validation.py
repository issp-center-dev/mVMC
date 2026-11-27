"""Input-validation tests for the ComplexUHF tool (UHF).

Each case copies the UHF_HubbardSquare input, changes one definition file,
runs UHF and checks that the run is either accepted, or rejected with a
clean nonzero exit status (not a signal) and an error message.
"""
from __future__ import print_function

import os
import re
import shutil
import subprocess
import unittest

import numpy as np

ROOT = os.getcwd()
DATA = os.path.join(ROOT, "data", "UHF_HubbardSquare")
WORK = os.path.join(ROOT, "work", "UHF_InputValidation")
UHF = os.path.join(ROOT, "..", "..", "src", "ComplexUHF", "UHF")
TIMEOUT = 300
HEADER_LINES = 5


def prepare(name):
    d = os.path.join(WORK, name)
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    for f in os.listdir(DATA):
        if f.endswith(".def"):
            shutil.copyfile(os.path.join(DATA, f), os.path.join(d, f))
    return d


def edit(d, fn, func):
    p = os.path.join(d, fn)
    with open(p) as f:
        text = f.read()
    with open(p, "w") as f:
        f.write(func(text))


def set_record(d, fn, row, new_line):
    """Replace data record `row` (0-based, after the header) of file `fn`."""
    def f(text):
        lines = text.splitlines(True)
        lines[HEADER_LINES + row] = new_line + "\n"
        return "".join(lines)
    edit(d, fn, f)


def append_record(d, fn, new_line):
    edit(d, fn, lambda t: t.rstrip("\n") + "\n" + new_line + "\n")


def run(d):
    return subprocess.run(
        [UHF, "namelist.def"],
        cwd=d,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=TIMEOUT,
    )


class TestUHFInputValidation(unittest.TestCase):
    def assertAccepted(self, res):
        self.assertEqual(res.returncode, 0, msg=res.stdout + res.stderr)

    def assertRejected(self, res):
        # A positive status is a clean error exit; a negative one is a signal.
        self.assertGreater(res.returncode, 0, msg=res.stdout + res.stderr)
        self.assertIn("Error", res.stdout + res.stderr)

    def test_reference_input_is_accepted(self):
        self.assertAccepted(run(prepare("ok")))

    def test_ne_gives_the_same_result_as_ncond(self):
        ref = prepare("ncond")
        self.assertAccepted(run(ref))
        d = prepare("ne")
        edit(d, "modpara.def",
             lambda t: re.sub(r"^Ncond\s+8.*$", "Ne             4", t, flags=re.M))
        self.assertAccepted(run(d))
        for fn in ("zvo_eigen.dat", "zqp_APOrbital_opt.dat", "zvo_UHF_cisajs.dat"):
            a = np.loadtxt(os.path.join(ref, fn), skiprows=0 if fn.startswith("zvo") else 5)
            b = np.loadtxt(os.path.join(d, fn), skiprows=0 if fn.startswith("zvo") else 5)
            np.testing.assert_allclose(a, b, rtol=0.0, atol=1e-12, err_msg=fn)

    def test_missing_electron_number_is_rejected(self):
        d = prepare("no_ncond")
        edit(d, "modpara.def", lambda t: re.sub(r"^Ncond.*\n", "", t, flags=re.M))
        self.assertRejected(run(d))

    def test_electron_number_outside_range_is_rejected(self):
        for name, value in (("ncond_zero", 0), ("ncond_full", 16), ("ncond_negative", -2)):
            d = prepare(name)
            edit(d, "modpara.def",
                 lambda t, v=value: re.sub(r"^Ncond\s+8.*$", "Ncond          %d" % v, t, flags=re.M))
            self.assertRejected(run(d))

    def test_missing_definition_file_is_rejected(self):
        d = prepare("missing_trans")
        os.remove(os.path.join(d, "trans.def"))
        self.assertRejected(run(d))

    def test_unknown_namelist_keyword_is_skipped(self):
        d = prepare("unknown_keyword")
        append_record(d, "namelist.def", "UnknownKeyword  unknown.def")
        res = run(d)
        self.assertAccepted(res)
        self.assertIn("Warning", res.stdout + res.stderr)

    def test_extra_records_are_rejected(self):
        d = prepare("extra_trans")
        append_record(d, "trans.def", "    0     0     1     0         1.000000000000000         0.000000000000000")
        self.assertRejected(run(d))

    def test_site_index_out_of_range_is_rejected(self):
        for fn, line in (
            ("trans.def", "    8     0     0     0         1.000000000000000         0.000000000000000"),
            ("coulombintra.def", "    8   4.000000000000000"),
            ("locspn.def", "    8      0"),
            ("orbitalidx.def", "    8      0      0"),
            ("initial.def", "    8     0     0     0   0.500000000000000   0.000000000000000"),
        ):
            d = prepare("site_" + fn.replace(".", "_"))
            set_record(d, fn, 0, line)
            self.assertRejected(run(d))

    def test_spin_index_out_of_range_is_rejected(self):
        d = prepare("spin_trans")
        set_record(d, "trans.def", 0, "    1     2     0     0         1.000000000000000         0.000000000000000")
        self.assertRejected(run(d))

    def test_orbital_index_out_of_range_is_rejected(self):
        d = prepare("orbital_index")
        set_record(d, "orbitalidx.def", 0, "    0      0     32")
        self.assertRejected(run(d))

    def test_short_record_is_rejected(self):
        for fn, line in (
            ("orbitalidx.def", "    0      1"),
            ("trans.def", "    1     0     0"),
            ("coulombintra.def", "    0"),
        ):
            d = prepare("short_" + fn.replace(".", "_"))
            set_record(d, fn, 0, line)
            self.assertRejected(run(d))

    def test_non_numeric_record_is_rejected(self):
        d = prepare("non_numeric")
        set_record(d, "coulombintra.def", 0, "    a   4.000000000000000")
        self.assertRejected(run(d))

    def test_odd_local_spin_count_with_ne_is_accepted(self):
        # Site 0 is a localized spin (NlocalSpin = 1) and modpara gives Ne, not
        # Ncond. 2*Ne = NlocalSpin + Ncond is even, so the derived odd Ncond is
        # accepted even though an orbital file is present (review R-1).
        d = prepare("odd_locspn_ne")
        set_record(d, "locspn.def", 0, "    0      1")
        edit(d, "locspn.def",
             lambda t: re.sub(r"^NlocalSpin\s+0", "NlocalSpin     1", t, flags=re.M))
        edit(d, "modpara.def",
             lambda t: re.sub(r"^Ncond\s+8.*$", "Ne             4", t, flags=re.M))
        self.assertAccepted(run(d))

    def test_local_spin_header_mismatch_is_rejected(self):
        # Only the header NlocalSpin enters the electron count, so it must
        # agree with the flags in the body (review R-2).
        d = prepare("locspn_header")
        edit(d, "locspn.def",
             lambda t: re.sub(r"^NlocalSpin\s+0", "NlocalSpin     2", t, flags=re.M))
        self.assertRejected(run(d))

    def test_duplicate_local_spin_site_is_rejected(self):
        d = prepare("locspn_duplicate")
        set_record(d, "locspn.def", 1, "    0      0")
        self.assertRejected(run(d))

    def test_negative_record_count_is_rejected(self):
        d = prepare("negative_count")
        edit(d, "trans.def", lambda t: re.sub(r"^NTransfer\s+48", "NTransfer      -1", t, flags=re.M))
        self.assertRejected(run(d))


if __name__ == "__main__":
    unittest.main()
