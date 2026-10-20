#!/usr/bin/env python3
"""Run controlled Phase 3A baseline performance samples in isolated processes."""
import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

BASELINE = "15582edab8112dbb63e0101be5b26f5135298644"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(source: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--label", default="baseline")
    parser.add_argument("--expected-commit", default=BASELINE)
    parser.add_argument("--samples", type=int, default=5)
    parser.add_argument("--batch-seconds", type=float, default=0.25)
    parser.add_argument("--scenario", choices=("mapped", "file_backed", "both"), default="both")
    args = parser.parse_args()

    executable = args.executable.resolve()
    source = args.source_root.resolve()
    output = args.output_dir.resolve()
    if not executable.is_file() or not source.is_dir():
        parser.error("executable and source root must exist")
    if args.samples < 3 or args.batch_seconds <= 0:
        parser.error("use at least three samples and a positive batch duration")
    commit = git(source, "rev-parse", "HEAD")
    if commit != args.expected_commit:
        parser.error(f"source HEAD {commit} does not match expected {args.expected_commit}")
    harness_root = Path(__file__).resolve().parent.parent
    harness_commit = git(harness_root, "rev-parse", "HEAD")
    production_changes = git(
        harness_root, "diff", "--name-only", args.expected_commit, "--",
        "gplspec", "examples/clean_bench_3.cpp", "examples/clean_bench_6.cpp",
    ).splitlines()
    if production_changes:
        parser.error("production sources or benchmark fixture settings differ from the baseline: "
                     + ", ".join(production_changes))
    time_bin = shutil.which("/usr/bin/time")
    if not time_bin:
        parser.error("/usr/bin/time is required to record per-process peak RSS")

    output.mkdir(parents=True, exist_ok=False)
    scenarios = ["mapped", "file_backed"] if args.scenario == "both" else [args.scenario]
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    measurement_rows = []
    raw_files = []
    sample_exit_codes = []
    residual_gate_failures = 0
    for scenario in scenarios:
        scenario_dir = output / scenario
        scenario_dir.mkdir()
        for sample in range(1, args.samples + 1):
            stem = f"sample-{sample:02d}"
            stdout_path = scenario_dir / f"{stem}.csv"
            stderr_path = scenario_dir / f"{stem}.stderr.log"
            rss_path = scenario_dir / f"{stem}.max-rss-kb.txt"
            command = [time_bin, "-f", "%M", "-o", str(rss_path),
                       str(executable), scenario, str(args.batch_seconds)]
            with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
                result = subprocess.run(command, cwd=source, env=env,
                                        stdout=stdout, stderr=stderr, check=False)
            sample_exit_codes.append({"scenario": scenario, "sample": sample,
                                      "exit_code": result.returncode})
            if result.returncode not in (0, 3):
                print(f"failed sample: {scenario}/{stem}; see {stderr_path}",
                      file=sys.stderr)
                return result.returncode
            if result.returncode == 3:
                residual_gate_failures += 1
            rss_kb = int(rss_path.read_text().strip().splitlines()[-1])
            with stdout_path.open(newline="") as stream:
                for row in csv.DictReader(stream):
                    row["sample"] = str(sample)
                    row["process_exit_code"] = str(result.returncode)
                    row["peak_rss_kb"] = str(rss_kb)
                    measurement_rows.append(row)
            for path in (stdout_path, stderr_path, rss_path):
                raw_files.append({"path": str(path.relative_to(output)),
                                  "sha256": sha256(path),
                                  "bytes": path.stat().st_size})

    columns = ["sample", "process_exit_code", "peak_rss_kb", "scenario", "metric", "index", "repeats",
               "total_seconds", "mean_seconds", "min_seconds", "median_seconds",
               "max_seconds", "solver_info", "iterations", "requested_tolerance",
               "solver_error", "absolute_residual", "relative_residual",
               "independent_gate_pass", "checksum",
               "target_met", "solver_setup_seconds", "initial_guess_seconds",
               "iteration_seconds"]
    results_path = output / "measurements.csv"
    with results_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(measurement_rows)

    inputs = ["modeldata/prem.200", "modeldata/S40RTS_dvs.nc"]
    manifest = {
        "label": args.label,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": commit,
        "expected_commit": args.expected_commit,
        "source_root": str(source),
        "harness_checkout": str(harness_root),
        "harness_commit": harness_commit,
        "production_source_changes_vs_baseline": production_changes,
        "executable": str(executable),
        "executable_sha256": sha256(executable),
        "runner_sha256": sha256(Path(__file__).resolve()),
        "harness_source_sha256": sha256(Path(__file__).with_name("phase3a_performance.cpp")),
        "host": platform.uname()._asdict(),
        "environment": {key: env[key] for key in
                        ("OMP_NUM_THREADS", "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "scenarios": scenarios,
        "samples_per_scenario": args.samples,
        "target_batch_seconds_per_timed_component": args.batch_seconds,
        "command_template": [time_bin, "-f", "%M", "-o", "<sample>.max-rss-kb.txt",
                             str(executable), "<scenario>", str(args.batch_seconds)],
        "fixture_settings": {
            "mapped": {"source_case": "clean_bench_3", "npoly": 5, "lmax": 2,
                       "max_step": 0.1, "ball_radius": 1.4, "map": "h*r, h=0.2",
                       "tolerance": 1e-12, "input_files": []},
            "file_backed": {"source_case": "clean_bench_6", "npoly": 12, "lmax": 20,
                            "max_step": 0.05, "ball_radius": 1.2, "tolerance": 1e-14,
                            "input_files": inputs},
        },
        "input_hashes": {item: sha256(source / item) for item in inputs},
        "sample_exit_codes": sample_exit_codes,
        "residual_gate_failures": residual_gate_failures,
        "scientific_gate_pass": residual_gate_failures == 0,
        "raw_samples": raw_files,
        "measurements_csv_sha256": sha256(results_path),
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"wrote {len(measurement_rows)} measurement rows to {results_path}")
    return 3 if residual_gate_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
