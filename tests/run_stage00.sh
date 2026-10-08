#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "usage: run_stage00.sh PATH_TO_STAGE00_REFERENCE_EXECUTABLE" >&2
  exit 2
fi
runner=$(realpath "$1")
root=$(cd "$(dirname "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
for run in 1 2; do
  "$runner" "$work/output-$run" "$work/actual-$run.csv" >"$work/run-$run.log"
  python3 "$root/tests/compare_reference.py" \
    "$root/tests/reference/stage00-baseline.csv" "$work/actual-$run.csv"
  cmp "$root/tests/reference/homogeneous-MatrixSolution.out" \
      "$work/output-$run/homogeneous/MatrixSolution.out"
done
python3 "$root/tests/compare_reference.py" \
  "$work/actual-1.csv" "$work/actual-2.csv"
echo "PASS: representative MatrixSolution.out is byte-identical to the frozen baseline"
