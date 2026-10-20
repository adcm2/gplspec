#!/usr/bin/env python3
"""Gate retained Stage 00/05 records by quantity, preserving frozen inputs."""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

# Only intermediate arrays removed by the tensor-only perturbation API are
# retired. The frozen files remain unchanged and all other records are exact
# coverage requirements.
RETIRED = {
    "stage00": ("mapping_perturbation/dxi/", "mapping_perturbation/dxilm/",
                "mapping_perturbation/df/"),
    "stage05_file": ("file_perturbation/dxi/", "file_perturbation/dxilm/",
                     "file_perturbation/df/"),
    "stage05_1d": (),
}
# These unit-radius fixtures have density at least 0.8*rho_norm over the body:
# homogeneous/layered density is >=1, the lateral fractional variation is
# within +/-0.2, and the smooth map preserves total mass through its Jacobian.
# Thus M >= 0.8*(4*pi/3)*Mnorm > 1.2*Mnorm while Rout <= 1.2*Rnorm, so
# GMnorm/Rnorm is a conservative lower bound on GM/R for the field scale.
FIELD_FLOOR = 1e-12 * 6.67430e-11 * 5.972e24 / 6_371_000.0
POTENTIAL_NORM = (6_371_000.0 / 3600.0) ** 2
L_MAX = 2


def read_stage00(path):
    rows = {}
    with Path(path).open(newline="") as stream:
        for line, row in enumerate(csv.reader(stream), 1):
            if len(row) not in (2, 3):
                raise ValueError(f"{path}:{line}: malformed Stage 00 row")
            key = row[0]
            if key in rows:
                raise ValueError(f"{path}:{line}: duplicate key {key}")
            values = tuple(float(value) for value in row[1:])
            if not values or not all(math.isfinite(value) for value in values):
                raise ValueError(f"{path}:{line}: nonfinite/empty row {key}")
            rows[key] = values
    return rows, 0


def read_stage05(path):
    rows = {}
    preamble = 0
    with Path(path).open(newline="") as stream:
        for line, row in enumerate(csv.reader(stream), 1):
            if line == 1 and row == ["Finding a"]:
                preamble = 1
                continue
            if len(row) not in (2, 3, 4, 5):
                raise ValueError(f"{path}:{line}: malformed Stage 05 row")
            if len(row) == 2:
                key, fields = row[0], row[1:]
            elif row[0].startswith("file_perturbation/"):
                key, fields = row[0], row[1:]
            elif row[0] in {"rows", "xAy", "Ax_y"}:
                key, fields = row[0], row[1:]
            else:
                key, fields = "/".join(row[:2]), row[2:]
            if key in rows:
                raise ValueError(f"{path}:{line}: duplicate key {key}")
            values = tuple(float(value) for value in fields)
            if not values or not all(math.isfinite(value) for value in values):
                raise ValueError(f"{path}:{line}: nonfinite/empty row {key}")
            rows[key] = values
    return rows, preamble


def category(key, kind):
    parts = key.split("/")
    if kind == "stage00":
        if len(parts) < 2 or parts[0] not in {
            "homogeneous", "layered", "lateral_density", "smooth_mapping",
            "mapping_perturbation"
        }:
            raise ValueError(f"unknown Stage 00 key: {key}")
        return parts[1]
    if kind == "stage05_file":
        if len(parts) < 2 or parts[0] != "file_perturbation":
            raise ValueError(f"unknown Stage 05 file key: {key}")
        return parts[1]
    return parts[0]


def complex_value(values):
    if len(values) == 1:
        return complex(values[0], 0.0)
    if len(values) == 2:
        return complex(*values)
    raise ValueError(f"unsupported record tuple length {len(values)}")


def relative_norm_bound(reference, actual, relative=1e-11, absolute_scale=1e-13):
    delta = math.sqrt(sum(abs(a - b) ** 2 for a, b in zip(reference, actual)))
    scale = math.sqrt(sum(abs(a) ** 2 for a in reference))
    limit = relative * scale + absolute_scale * max(1.0, scale)
    return delta, limit


def gate(reference_path, actual_path, kind):
    reader = read_stage00 if kind == "stage00" else read_stage05
    reference, reference_preamble = reader(reference_path)
    actual, actual_preamble = reader(actual_path)
    retired = {key for key in reference if key.startswith(RETIRED[kind])}
    expected = reference.keys() - retired
    if expected != actual.keys():
        missing = sorted(expected - actual.keys())
        extra = sorted(actual.keys() - expected)
        raise ValueError(f"key coverage failure: missing={missing[:3]} extra={extra[:3]}")
    if reference_preamble != actual_preamble:
        raise ValueError("preamble count differs from frozen reference")
    for key, old in reference.items():
        if key in retired:
            continue
        if len(old) != len(actual[key]):
            raise ValueError(f"tuple-count mismatch for {key}")

    # Compare complex solution harmonics together at each radius. The addition
    # theorem bounds pointwise physical error, and therefore any positive-
    # weight physical-volume RMS. Requiring the conservative bound below the
    # smallest predeclared 1e-12 GM/R floor is sufficient to pass the full gate.
    solutions = defaultdict(dict)
    vectors = defaultdict(list)
    scalar_bounds = []
    for key in expected:
        cat = category(key, kind)
        old = reference[key]
        new = actual[key]
        if kind == "stage00":
            parts = key.split("/")
            if cat == "solution":
                radial = "/".join(parts[:-1])
                harmonic = parts[-1]
                solutions[radial][int(harmonic)] = (
                    complex_value(old), complex_value(new))
            elif cat in {"source", "operator_y", "da"}:
                vectors[parts[0] + "/" + cat].append(
                    (complex_value(old), complex_value(new)))
            elif cat == "operator_x":
                if old != new:
                    raise ValueError(f"deterministic action input changed at {key}")
            elif cat in {"density", "mapping", "jacobian", "inverse_f",
                         "laplace_tensor"}:
                scale = max(1.0, *(abs(value) for value in old))
                for lhs, rhs in zip(old, new):
                    if abs(lhs - rhs) > 1e-12 * abs(lhs) + 1e-14 * scale:
                        raise ValueError(f"geometry/material criterion failed at {key}")
            else:
                raise ValueError(f"unclassified Stage 00 quantity {key}")
        elif kind == "stage05_file":
            parts = key.split("/")
            if cat == "solution":
                radial = "/".join(parts[:-1])
                harmonic = parts[-1]
                solutions[radial][int(harmonic)] = (
                    complex_value(old), complex_value(new))
            elif cat in {"source", "operator_y", "operator_ay", "da",
                         "Ax_y", "xAy"}:
                vectors["file/" + cat].append((complex_value(old), complex_value(new)))
            elif cat == "operator_x":
                if old != new:
                    raise ValueError(f"deterministic action input changed at {key}")
            else:
                raise ValueError(f"unclassified Stage 05 quantity {key}")
        else:
            parts = key.split("/")
            if cat in {"Ax", "boundary_delta", "xAy", "Ax_y"}:
                vectors["1d/action/" + cat].append(
                    (complex_value(old[-2:]), complex_value(new[-2:])))
            elif cat == "force":
                vectors["1d/source"].append(
                    (complex_value(old[-2:]), complex_value(new[-2:])))
            elif cat in {"rows", "interface", "origin"}:
                if cat == "rows":
                    if old != new:
                        raise ValueError("Stage 05 1-D record count changed")
                else:
                    scale = max(1.0, *(abs(value) for value in old))
                    if any(abs(lhs - rhs) > 1e-12 * abs(lhs) + 1e-14 * scale
                           for lhs, rhs in zip(old, new)):
                        raise ValueError(f"interface/origin criterion failed at {key}")
            else:
                raise ValueError(f"unclassified Stage 05 1-D quantity {key}")

    for label, values in vectors.items():
        old = [pair[0] for pair in values]
        new = [pair[1] for pair in values]
        delta, limit = relative_norm_bound(old, new)
        if delta > limit:
            raise ValueError(f"algebra/action norm failed for {label}: {delta:g} > {limit:g}")
    if kind != "stage05_1d" and not solutions:
        raise ValueError("no potential solution records were compared")
    for label, harmonics in sorted(solutions.items()):
        if set(harmonics) != set(range((L_MAX + 1) ** 2)):
            raise ValueError(f"solution harmonic coverage is incomplete for {label}")
        coeff_delta = math.sqrt(sum(abs(new - old) ** 2
                                    for old, new in harmonics.values()))
        # All Stage 00 and file-perturbation solution fixtures use lmax=2.
        bound = (L_MAX + 1) / math.sqrt(4.0 * math.pi) * coeff_delta * POTENTIAL_NORM
        scalar_bounds.append((label, bound))
        if bound > FIELD_FLOOR:
            raise ValueError(f"physical potential bound failed for {label}: "
                             f"{bound:g} > {FIELD_FLOOR:g} m^2/s^2")
    print(f"PASS {kind}: retained {len(expected)} records; retired {len(retired)}; "
          f"physical potential floor={FIELD_FLOOR:.9g} m^2/s^2")
    if scalar_bounds:
        label, bound = max(scalar_bounds, key=lambda item: item[1])
        print(f"worst pointwise/RMS potential bound ({label})="
              f"{bound:.9g} m^2/s^2")
    for label, values in sorted(vectors.items()):
        delta, limit = relative_norm_bound([p[0] for p in values],
                                           [p[1] for p in values])
        print(f"{label}: norm difference={delta:.9g}, limit={limit:.9g}")


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in {
        "stage00", "stage05_1d", "stage05_file"
    }:
        raise SystemExit("usage: scientific_reference_gate.py stage00|stage05_1d|stage05_file REFERENCE.csv ACTUAL.csv")
    try:
        gate(sys.argv[2], sys.argv[3], sys.argv[1])
    except (ValueError, OSError) as exc:
        raise SystemExit(f"SCIENTIFIC GATE FAIL: {exc}")
