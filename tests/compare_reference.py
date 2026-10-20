#!/usr/bin/env python3
"""Compare retained stage-00 records with frozen baseline data.

MappingPerturbation no longer exposes its temporary displacement, harmonic,
and gradient buffers. Their old rows remain in the frozen CSV as evidence;
the focused 3E.2 finite-difference tests replace that coverage.
"""
import csv
import math
import sys
from pathlib import Path

retire_mapping_intermediates = (
    len(sys.argv) == 4 and sys.argv[1] == "--retire-mapping-intermediates"
)
if len(sys.argv) not in (3, 4) or (
    len(sys.argv) == 4 and not retire_mapping_intermediates
):
    raise SystemExit(
        "usage: compare_reference.py [--retire-mapping-intermediates] BASELINE.csv ACTUAL.csv"
    )
paths = sys.argv[2:] if retire_mapping_intermediates else sys.argv[1:]


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


reference, actual = read(paths[0]), read(paths[1])
retired_prefixes = (
    "mapping_perturbation/dxi/",
    "mapping_perturbation/dxilm/",
    "mapping_perturbation/df/",
)
retired = (
    {key for key in reference if key.startswith(retired_prefixes)}
    if retire_mapping_intermediates
    else set()
)
expected_keys = reference.keys() - retired
if expected_keys != actual.keys():
    missing, extra = expected_keys - actual.keys(), actual.keys() - expected_keys
    raise SystemExit(f"record mismatch: missing={len(missing)} extra={len(extra)}")
rtol, atol = 1e-13, 1e-15
max_abs = 0.0
max_scaled = 0.0
changed = 0
for key, expected in reference.items():
    if key in retired:
        continue
    observed = actual[key]
    if len(expected) != len(observed):
        raise SystemExit(f"component mismatch: {key}")
    for lhs, rhs in zip(expected, observed):
        delta = abs(lhs - rhs)
        max_abs = max(max_abs, delta)
        max_scaled = max(max_scaled, delta / max(atol, abs(lhs)))
        if delta > atol + rtol * abs(lhs):
            raise SystemExit(
                f"strict tolerance failed: {key}: expected={lhs:.17g} actual={rhs:.17g}"
            )
        changed += delta != 0.0
print(f"PASS: retained {len(expected_keys)} of {len(reference)} frozen records; "
      f"retired intermediate records={len(retired)}; changed components={changed}; "
      f"max_abs={max_abs:.3g}; max_scaled={max_scaled:.3g}; "
      f"rtol={rtol:g}; atol={atol:g}")
