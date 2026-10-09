"""Focused regression for the retained CleanBench 7 sweep plot."""
import os
os.environ.setdefault("MPLBACKEND", "Agg")

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "work"))
import plot_fields


class CleanBench7SweepTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(plt.close, "all")
        self.table = Path(self.temp.name) / "sweep.csv"
        self.h = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05])
        self.approximate = np.full_like(self.h, 10.0)
        self.exact = np.array([10.0, 10.1, 10.4, 10.9, 11.6, 12.5])
        self.write_table(np.column_stack((self.h, self.approximate, self.exact)))

    def write_table(self, values):
        np.savetxt(self.table, values, delimiter=";")

    def test_rendered_curves_normalization_quadratic_guide_and_zero_mask(self):
        args = SimpleNamespace(table=str(self.table), output=None)
        with patch.object(plot_fields, "finish"):
            plot_fields.sweep_plot(args)

        left, right = plt.gcf().axes
        np.testing.assert_array_equal(left.lines[0].get_xdata(), self.h)
        np.testing.assert_array_equal(left.lines[0].get_ydata(), self.exact)
        np.testing.assert_array_equal(left.lines[1].get_xdata(), self.h)
        np.testing.assert_array_equal(left.lines[1].get_ydata(), self.approximate)
        self.assertEqual(left.get_ylabel(), "Centre potential (m² s⁻²)")

        expected_error = np.abs(self.exact - self.approximate) / np.max(np.abs(self.approximate))
        expected_guide = expected_error[3] * (self.h / self.h[3]) ** 2
        valid = (self.h > 0) & (expected_error > 0)
        np.testing.assert_array_equal(right.lines[0].get_xdata(), self.h[valid])
        np.testing.assert_allclose(right.lines[0].get_ydata(), expected_error[valid] * 100.0)
        np.testing.assert_array_equal(right.lines[1].get_xdata(), self.h[self.h > 0])
        np.testing.assert_allclose(right.lines[1].get_ydata(), expected_guide[self.h > 0] * 100.0)
        self.assertFalse(np.any(right.lines[0].get_xdata() == 0.0))
        self.assertEqual(left.lines[0].get_xdata()[0], 0.0)
        self.assertEqual(right.get_yscale(), "log")
        self.assertEqual(right.get_xscale(), "log")

    def test_rejects_wrong_shape_nonfinite_zero_scale_and_invalid_anchor(self):
        invalid_tables = [
            (np.ones((4, 2)), "exactly three columns"),
            (np.array([[0, 1, 1], [1, 1, np.nan], [2, 1, 2], [3, 1, 4]]), "non-finite"),
            (np.column_stack((self.h, np.zeros_like(self.h), self.exact)), "nonzero scale"),
            (np.column_stack((self.h, self.approximate, self.approximate)), "fourth sample"),
        ]
        for values, message in invalid_tables:
            with self.subTest(message=message):
                self.write_table(values)
                with self.assertRaisesRegex(ValueError, message):
                    plot_fields.load_sweep_table(self.table)


if __name__ == "__main__":
    unittest.main()
