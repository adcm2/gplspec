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
  if python3 "$root/tests/compare_reference.py" --retire-mapping-intermediates \
      "$root/tests/reference/stage00-baseline.csv" "$work/actual-$run.csv" \
      >"$work/strict-$run.log" 2>&1; then
    echo "DIAGNOSTIC strict Stage 00 comparison run $run: PASS"
  else
    echo "DIAGNOSTIC strict Stage 00 comparison run $run: DIFFERENCES (see details below)"
    sed -n '1,8p' "$work/strict-$run.log"
  fi
  python3 "$root/tests/scientific_reference_gate.py" stage00 \
    "$root/tests/reference/stage00-baseline.csv" "$work/actual-$run.csv"
  if cmp -s "$root/tests/reference/homogeneous-MatrixSolution.out" \
      "$work/output-$run/homogeneous/MatrixSolution.out"; then
    echo "DIAGNOSTIC frozen MatrixSolution.out comparison run $run: byte-identical"
  else
    echo "DIAGNOSTIC frozen MatrixSolution.out comparison run $run: differs"
  fi
done
if python3 "$root/tests/compare_reference.py" \
    "$work/actual-1.csv" "$work/actual-2.csv" >"$work/repeat-strict.log" 2>&1; then
  echo "DIAGNOSTIC repeated Stage 00 strict comparison: PASS"
else
  echo "DIAGNOSTIC repeated Stage 00 strict comparison: DIFFERENCES"
  sed -n '1,8p' "$work/repeat-strict.log"
fi
python3 "$root/tests/scientific_reference_gate.py" stage00 \
  "$work/actual-1.csv" "$work/actual-2.csv"
echo "PASS: Stage 00 scientific quantity gate passed; strict and bytewise comparisons above are diagnostics."
