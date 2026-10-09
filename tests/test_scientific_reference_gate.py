#!/usr/bin/env python3
"""Negative checks for scientific-reference coverage and physical limits."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "scientific_reference_gate", Path(__file__).with_name("scientific_reference_gate.py"))
gate_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate_module)


class ScientificReferenceGateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.reference = self.root / "reference.csv"
        self.actual = self.root / "actual.csv"
        self.rows = [[f"homogeneous/solution/0/0/{i}", "0", "0"]
                     for i in range(9)]
        self.rows += [["homogeneous/source/0", "1", "0"]]
        self.write(self.reference, self.rows)
        self.write(self.actual, self.rows)

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def write(path, rows):
        with path.open("w", newline="") as stream:
            csv.writer(stream).writerows(rows)

    def test_matching_retained_records_pass(self):
        gate_module.gate(self.reference, self.actual, "stage00")

    def assert_gate_fails(self, rows, fragment):
        self.write(self.actual, rows)
        with self.assertRaisesRegex(ValueError, fragment):
            gate_module.gate(self.reference, self.actual, "stage00")

    def test_missing_record_fails(self):
        self.assert_gate_fails(self.rows[:-1], "key coverage failure")

    def test_extra_record_fails(self):
        self.assert_gate_fails(self.rows + [["homogeneous/source/extra", "0"]],
                               "key coverage failure")

    def test_nonfinite_record_fails(self):
        rows = list(self.rows)
        rows[-1] = ["homogeneous/source/0", "nan", "0"]
        self.assert_gate_fails(rows, "nonfinite")

    def test_excessive_action_error_fails(self):
        rows = list(self.rows)
        rows[-1] = ["homogeneous/source/0", "2", "0"]
        self.assert_gate_fails(rows, "algebra/action norm failed")

    def test_excessive_physical_potential_error_fails(self):
        rows = list(self.rows)
        rows[0] = ["homogeneous/solution/0/0/0", "1e-8", "0"]
        self.assert_gate_fails(rows, "physical potential bound failed")


if __name__ == "__main__":
    unittest.main()
