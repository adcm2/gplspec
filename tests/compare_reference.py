#!/usr/bin/env python3
"""Compare indexed stage-00 numerical records with frozen baseline data."""
import csv
import math
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit("usage: compare_reference.py BASELINE.csv ACTUAL.csv")


def read(path):
    records = {}
    with Path(path).open(newline="") as stream:
        for row in csv.reader(stream):
            if len(row) not in (2, 3):
                raise SystemExit(f"{path}: malformed record: {row[:4]}")
            key = row[0]
            if key in records:
                raise SystemExit(f"{path}: duplicate record {key}")
            try:
                values = tuple(float(item) for item in row[1:])
            except ValueError as exc:
                raise SystemExit(f"{path}: nonnumeric record {key}: {exc}")
            if not all(math.isfinite(item) for item in values):
                raise SystemExit(f"{path}: nonfinite record {key}")
            records[key] = values
    return records


reference, actual = read(sys.argv[1]), read(sys.argv[2])
if reference.keys() != actual.keys():
    missing, extra = reference.keys() - actual.keys(), actual.keys() - reference.keys()
    raise SystemExit(f"record mismatch: missing={len(missing)} extra={len(extra)}")
rtol, atol = 1e-13, 1e-15
max_abs = 0.0
max_scaled = 0.0
changed = 0
for key, expected in reference.items():
    observed = actual[key]
    if len(expected) != len(observed):
        raise SystemExit(f"component mismatch: {key}")
    for lhs, rhs in zip(expected, observed):
        delta = abs(lhs - rhs)
        max_abs = max(max_abs, delta)
        max_scaled = max(max_scaled, delta / max(atol, abs(lhs)))
        if delta > atol + rtol * abs(lhs):
            raise SystemExit(f"tolerance failed: {key}: expected={lhs:.17g} actual={rhs:.17g}")
        changed += delta != 0.0
print(f"PASS: {len(reference)} records; changed components={changed}; max_abs={max_abs:.3g}; max_scaled={max_scaled:.3g}; rtol={rtol:g}; atol={atol:g}")
