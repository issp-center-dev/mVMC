"""Regression checks for the OptTrans optimization test's output validation."""

import unittest

from runtest_opttrans import check_weights


def valid_rows():
    return [["0"] * 6 + ["1", "0", "0", str(step / 100.0), "0", "0"]
            for step in range(100)]


class OptTransOutputTest(unittest.TestCase):
    def test_valid_optimization(self):
        self.assertEqual(check_weights(valid_rows()), (1.0, 0.99))

    def test_nonfinite_weights(self):
        for value in ("nan", "inf", "-inf"):
            for step in (0, 49, 99):
                for column in (-6, -5, -3, -2):
                    with self.subTest(value=value, step=step, column=column):
                        rows = valid_rows()
                        rows[step][column] = value
                        with self.assertRaisesRegex(ValueError, "non-finite"):
                            check_weights(rows)

    def test_all_nan_weights(self):
        rows = valid_rows()
        for row in rows:
            row[-6] = row[-3] = "nan"
        with self.assertRaisesRegex(ValueError, "non-finite"):
            check_weights(rows)

    def test_incomplete_output(self):
        for count in (0, 1, 99, 101):
            with self.subTest(count=count):
                rows = (valid_rows() * 2)[:count]
                with self.assertRaisesRegex(ValueError, "expected 100 SR steps"):
                    check_weights(rows)

    def test_malformed_row(self):
        rows = valid_rows()
        rows[50] = rows[50][:-1]
        with self.assertRaisesRegex(ValueError, "malformed"):
            check_weights(rows)

    def test_fixed_weight(self):
        rows = valid_rows()
        for row in rows:
            row[-3] = "0"
        with self.assertRaisesRegex(ValueError, "did not grow"):
            check_weights(rows)

    def test_invalid_normalization(self):
        rows = valid_rows()
        rows[50][-6] = "2"
        with self.assertRaisesRegex(ValueError, "not normalized"):
            check_weights(rows)

    def test_imaginary_weight(self):
        rows = valid_rows()
        rows[50][-2] = "0.1"
        with self.assertRaisesRegex(ValueError, "not real"):
            check_weights(rows)


if __name__ == "__main__":
    unittest.main()
