#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 STAGE05_OPERATOR_EXECUTABLE ORIGINAL_REFERENCE_CSV" >&2
  exit 2
fi
exe=$1
reference=$2
actual=$(mktemp)
trap 'rm -f "$actual"' EXIT
"$exe" > "$actual"
diff -u "$reference" "$actual"
