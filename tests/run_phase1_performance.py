#!/usr/bin/env python3
"""Run five controlled Stage 1D processes and aggregate their CSV measurements."""

import csv
import json
import os
import pathlib
import statistics
import subprocess
import sys


def summarize(values):
    median = statistics.median(values)
    mad = statistics.median(abs(value - median) for value in values)
    mean = statistics.mean(values)
    result = {
        "samples": len(values),
        "median": median,
        "min": min(values),
        "max": max(values),
        "coefficient_of_variation": statistics.pstdev(values) / mean if mean else 0.0,
        "median_absolute_deviation": mad,
        "median_plus_3_mad": median + 3.0 * mad,
    }
    if values and max(values) > 10000:
        result["max_plus_10_percent"] = max(values) * 1.10
    return result


if len(sys.argv) != 3:
    raise SystemExit("usage: run_phase1_performance.py EXECUTABLE OUTPUT_DIRECTORY")

executable = pathlib.Path(sys.argv[1]).resolve()
output = pathlib.Path(sys.argv[2]).resolve()
repo = pathlib.Path(__file__).resolve().parents[1]
if not executable.is_file():
    raise SystemExit(f"benchmark executable not found: {executable}")
output.mkdir(parents=True, exist_ok=True)

environment = os.environ.copy()
environment.update({
    "OMP_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
})

rows = []
for process in range(1, 6):
    name = f"process-{process:02d}.csv"
    result = subprocess.run(
        [str(executable), str(process)],
        cwd=repo,
        env=environment,
        check=True,
        text=True,
        capture_output=True,
    )
    (output / name).write_text(result.stdout)
    if result.stderr:
        (output / f"process-{process:02d}.stderr.log").write_text(result.stderr)
    rows.extend(csv.DictReader(result.stdout.splitlines()))

if not rows:
    raise RuntimeError("benchmark runs emitted no CSV records")

for row in rows:
    row["value"] = float(row["value"])
iterations = [row["value"] for row in rows
              if row["metric"] == "public_solver_iterations"]
diagnostic_iterations = [row["value"] for row in rows
                         if row["metric"] == "diagnostic_iterations"]
relative_residuals = [row["value"] for row in rows
                      if row["metric"] == "independent_relative_residual"]
if len(iterations) != 30 or any(value != 3 for value in iterations):
    raise RuntimeError("public solve iterations changed from the Stage 1C three-iteration fixture")
if len(diagnostic_iterations) != 5 or any(value != 3 for value in diagnostic_iterations):
    raise RuntimeError("diagnostic solver iteration count changed")
if len(relative_residuals) != 5 or any(value > 1e-6 for value in relative_residuals):
    raise RuntimeError("independent residual exceeded the existing 1e-6 tolerance")

groups = {}
for row in rows:
    if row["unit"] in ("s", "KiB"):
        key = f"{row['scenario']}/{row['metric']}"
        groups.setdefault(key, []).append(row["value"])

summary = {
    "processes": 5,
    "warm_repeats_per_process": 5,
    "thread_environment": {
        key: environment[key]
        for key in ("OMP_NUM_THREADS", "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    },
    "metrics": {key: summarize(values) for key, values in sorted(groups.items())},
    "correctness": {
        "public_solve_iteration_count": 3,
        "public_solve_samples": len(iterations),
        "independent_relative_residual_max": max(relative_residuals),
    },
}
(output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(f"wrote {len(rows)} measurements from five fresh processes to {output}")
