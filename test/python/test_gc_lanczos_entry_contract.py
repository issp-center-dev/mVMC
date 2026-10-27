"""Reject crashes, unrelated failures and stale output in finite-S entry tests."""
import pathlib
import subprocess
import tempfile
import unittest

from runtest_gc_lanczos import CheckFailure, require_finite_spin_entry_outcome


class EntryContract(unittest.TestCase):
    def setUp(self):
        # CTest's working directory is outside the source tree.
        self.temp = tempfile.TemporaryDirectory(dir=".")
        self.addCleanup(self.temp.cleanup)
        self.work = pathlib.Path(self.temp.name)
        self.error = (
            "Error: stabilized power-Lanczos failed: amplitude failure.\n"
            "GC Lanczos failure: stage=final-sample chain=0 N=4 "
            "status=amplitude failure amplitude_state=numeric_zero "
            "log_abs=-30 log_error_bound=-inf retry_reason=numeric_zero\n")

    def check(self, code=1, error=None, backend="off"):
        done = subprocess.CompletedProcess([], code, "", self.error if error is None else error)
        require_finite_spin_entry_outcome("contract", done, self.work, 1, backend)

    def test_allowed_off_stops(self):
        self.check()
        self.check(error=self.error.replace("amplitude_state=numeric_zero", "amplitude_state=finite")
                   .replace("retry_reason=numeric_zero", "retry_reason=local_energy"))

    def test_exit_codes_and_on_failure(self):
        for code in (0, 2, 99, -11, 139):
            with self.subTest(code=code), self.assertRaises(CheckFailure):
                self.check(code=code)
        with self.assertRaises(CheckFailure):
            self.check(backend="on")

    def test_failure_classification(self):
        mutations = (
            ("amplitude failure", "nonfinite value"),
            ("amplitude failure", "resource limit"),
            ("final-sample", "final-bootstrap"),
            ("final-sample", "pilot-sample"),
            ("retry_reason=numeric_zero", "retry_reason=no_retry"),
            ("retry_reason=numeric_zero", "retry_reason=local_energy"),
            ("amplitude_state=numeric_zero", "amplitude_state=nonfinite"),
            ("chain=0", "chain=1"), ("log_abs=-30", "log_abs=nan"))
        for old, new in mutations:
            with self.subTest(change=new), self.assertRaises(CheckFailure):
                self.check(error=self.error.replace(old, new))
        with self.assertRaises(CheckFailure):
            self.check(error=self.error + self.error)

    def test_failed_run_cannot_publish_results(self):
        output = self.work / "output"
        output.mkdir()
        (output / "zvo_gc_pl_number_001.dat").write_text("stale result\n")
        with self.assertRaises(CheckFailure):
            self.check()


if __name__ == "__main__":
    unittest.main()
