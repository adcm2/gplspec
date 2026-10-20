#!/usr/bin/env python3
"""Report frozen Stage 05 parity after retiring temporary mapping buffers."""
import csv
import math
import sys
from pathlib import Path


def read(path):
    records = {}
    preamble = None
    with Path(path).open(newline="") as stream:
        for line, row in enumerate(csv.reader(stream), 1):
            if line == 1 and row == ["Finding a"]:
                preamble = row[0]
                continue
            if len(row) != 3:
                raise SystemExit(f"{path}:{line}: malformed numeric record {row[:4]}")
            key = row[0]
            if key in records:
                raise SystemExit(f"{path}:{line}: duplicate key {key}")
            try:
                values = tuple(float(value) for value in row[1:])
            except ValueError as exc:
                raise SystemExit(f"{path}:{line}: invalid value for {key}: {exc}")
            if not all(math.isfinite(value) for value in values):
                raise SystemExit(f"{path}:{line}: nonfinite value for {key}")
            records[key] = values
    return preamble, records


if len(sys.argv) != 3:
    raise SystemExit("usage: compare_stage05_file_perturbation.py REFERENCE.csv ACTUAL.csv")

reference_preamble, reference = read(sys.argv[1])
actual_preamble, actual = read(sys.argv[2])
retired = {"dxi", "dxilm", "df"}


def category(key):
    fields = key.split("/")
    if len(fields) < 2 or fields[0] != "file_perturbation":
        raise SystemExit(f"unexpected Stage 05 record key: {key}")
    return fields[1]


expected = {key for key in reference if category(key) not in retired}
if reference_preamble != actual_preamble:
    raise SystemExit("frozen and candidate preambles differ")
if expected != actual.keys():
    missing, extra = expected - actual.keys(), actual.keys() - expected
    raise SystemExit(f"retained key mismatch: missing={len(missing)} extra={len(extra)}")

changed_components = 0
for key in sorted(expected):
    lhs, rhs = reference[key], actual[key]
    if len(lhs) != len(rhs):
        raise SystemExit(f"component-count mismatch: {key}")
    for expected_value, actual_value in zip(lhs, rhs):
        difference = abs(expected_value - actual_value)
        changed_components += difference != 0.0
        if expected_value != actual_value:
            raise SystemExit(
                f"strict frozen numeric comparison differs at {key}: "
                f"expected={expected_value:.17g} actual={actual_value:.17g}"
            )

print(
    f"PASS: retained {len(expected)} frozen records; retired "
    f"{sum(category(key) in retired for key in reference)} intermediate records; "
    f"changed components={changed_components}; retained components are numerically exact."
)
