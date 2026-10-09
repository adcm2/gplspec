#!/usr/bin/env python3
"""Summarize frozen/candidate indexed CSVs without editing either input."""
import csv
import math
import sys
from collections import defaultdict

PREAMBLE = ("Finding a",)
INDEXED_LABELS = {"Ax", "boundary_delta", "force", "interface", "origin"}

if len(sys.argv) != 3:
    raise SystemExit("usage: summarize.py FROZEN.csv ACTUAL.csv")

def read(path):
    records = {}
    preambles = 0
    with open(path, newline="") as stream:
        for line_number, row in enumerate(csv.reader(stream), 1):
            if tuple(row) == PREAMBLE:
                preambles += 1
                continue
            if not row:
                raise SystemExit(f"{path}:{line_number}: empty record")
            if row[0] in INDEXED_LABELS and len(row) >= 3:
                key, fields = f"{row[0]}/{row[1]}", row[2:]
            else:
                key, fields = row[0], row[1:]
            if key in records:
                raise SystemExit(f"{path}:{line_number}: duplicate key {key}")
            try:
                values = tuple(float(value) for value in fields)
            except ValueError as exc:
                raise SystemExit(f"{path}:{line_number}: invalid numeric record {key}: {exc}")
            if not values or not all(math.isfinite(value) for value in values):
                raise SystemExit(f"{path}:{line_number}: empty or non-finite numeric record {key}")
            records[key] = values
    return records, preambles

reference, reference_preambles = read(sys.argv[1])
actual, actual_preambles = read(sys.argv[2])
if reference.keys() != actual.keys():
    missing = sorted(reference.keys() - actual.keys())
    extra = sorted(actual.keys() - reference.keys())
    raise SystemExit(f"record-key mismatch: missing={missing[:3]} extra={extra[:3]}")
if reference_preambles != actual_preambles:
    raise SystemExit(f"preamble-count mismatch: {reference_preambles} != {actual_preambles}")

changed = defaultdict(int)
max_abs = defaultdict(float)
category_changed = defaultdict(int)
changed_records = defaultdict(int)
source_scaled = defaultdict(float)
radial_squared = defaultdict(float)
for key, old in reference.items():
    new = actual[key]
    if len(old) != len(new):
        raise SystemExit(f"tuple-count mismatch for {key}: {len(old)} != {len(new)}")
    group = key.split("/", 1)[0]
    category = ("source" if group == "force" else
                key.split("/", 2)[1] if "/" in key else group)
    record_changed = False
    for x, y in zip(old, new):
        delta = abs(x - y)
        record_changed |= x != y
        changed[group] += x != y
        category_changed[(group, category)] += x != y
        max_abs[group] = max(max_abs[group], delta)
        if category == "source":
            source_scaled[group] = max(source_scaled[group],
                                       delta / max(1e-13, abs(x)))
        if category == "solution":
            radial_squared[(group, key.rsplit("/", 1)[0])] += delta * delta
    if record_changed:
        changed_records[(group, category)] += 1

print(f"numeric_records={len(reference)} preamble_records={reference_preambles} "
      f"changed_numeric_components={sum(changed.values())}")
for group in sorted(changed):
    if changed[group]:
        print(f"{group}: changed={changed[group]} max_abs={max_abs[group]:.17g}")
for (group, category), count in sorted(category_changed.items()):
    if count:
        print(f"{group}/{category}: changed_components={count} "
              f"changed_records={changed_records[(group, category)]}")
for group, scaled in sorted(source_scaled.items()):
    if scaled:
        print(f"{group}/source_scaled_error_at_1e-13_floor={scaled:.17g}")
radial_norms = defaultdict(list)
for (group, _radial_row), squared_norm in radial_squared.items():
    radial_norms[group].append(math.sqrt(squared_norm))
potential_norm = (6371000.0 / 3600.0) ** 2
for group, norms in sorted(radial_norms.items()):
    coefficient_norm = max(norms)
    bound = 3.0 / math.sqrt(4.0 * math.pi) * coefficient_norm * potential_norm
    print(f"{group}/max_radial_solution_coefficient_l2={coefficient_norm:.17g}")
    print(f"{group}/lmax2_pointwise_potential_bound_m2s2={bound:.17g}")
