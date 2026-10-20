#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 stage05_1d|stage05_file EXECUTABLE ORIGINAL_REFERENCE_CSV" >&2
  exit 2
fi
kind=$1
exe=$2
reference=$3
root=$(cd "$(dirname "$0")/.." && pwd)
actual=$(mktemp)
strict=$(mktemp)
trap 'rm -f "$actual" "$strict"' EXIT
"$exe" > "$actual"

case "$kind" in
  stage05_1d)
    if diff -u "$reference" "$actual" >"$strict" 2>&1; then
      echo "DIAGNOSTIC strict Stage 05 1-D comparison: PASS"
    else
      echo "DIAGNOSTIC strict Stage 05 1-D comparison: DIFFERENCES"
      sed -n '1,8p' "$strict"
    fi
    ;;
  stage05_file)
    if python3 "$root/tests/compare_stage05_file_perturbation.py" \
        "$reference" "$actual" >"$strict" 2>&1; then
      echo "DIAGNOSTIC strict Stage 05 file comparison: PASS"
    else
      echo "DIAGNOSTIC strict Stage 05 file comparison: DIFFERENCES"
      sed -n '1,8p' "$strict"
    fi
    ;;
  *)
    echo "unknown scientific reference kind: $kind" >&2
    exit 2
    ;;
esac

python3 "$root/tests/scientific_reference_gate.py" \
  "$kind" "$reference" "$actual"
