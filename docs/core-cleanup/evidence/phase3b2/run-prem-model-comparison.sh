#!/usr/bin/env bash
set -euo pipefail

EVIDENCE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$EVIDENCE/../../../.." && pwd)"
BASE=/tmp/gplspec-phase3-baseline
BASE_BUILD=/tmp/gplspec-phase3-baseline-release
CANDIDATE_BUILD=/tmp/gplspec-phase3-candidate-release
OUT=/tmp/gplspec-phase3b2-prem-repro
mkdir -p "$OUT"

DEPS=(
  -I_DEP_/gaussquad-src
  -I_DEP_/numericconcepts-src/include
  -I_DEP_/gshtrans-src
  -I_DEP_/fftwpp-src
  -I_DEP_/planetarymodel-src
  -I_DEP_/interpolation-src
  -I_DEP_/tomographymodels-src
  -isystem_DEP_/eigen3-src
)
compile_probe() {
  local source="$1" root="$2" build="$3" output="$4"
  local dep_args=()
  for arg in "${DEPS[@]}"; do
    dep_args+=("${arg//_DEP_/$build/_deps}")
  done
  /usr/bin/c++ -O0 -DNDEBUG -std=gnu++23 -fopenmp -I"$root" \
    "${dep_args[@]}" "$source" -o "$output" -lfftw3 -lfftw3f \
    -lfftw3l -lnetcdf
}

compile_probe "$ROOT/docs/core-cleanup/evidence/phase3b1/prem-model-probe.cpp" \
  "$BASE" "$BASE_BUILD" "$OUT/prem-old"
compile_probe "$EVIDENCE/prem-model-probe.cpp" "$ROOT" \
  "$CANDIDATE_BUILD" "$OUT/prem-new"
(
  cd "$BASE"
  OMP_NUM_THREADS=1 "$OUT/prem-old" modeldata/PREM500 > "$OUT/PREM500-old.csv"
  OMP_NUM_THREADS=1 "$OUT/prem-old" modeldata/prem.200 > "$OUT/prem.200-old.csv"
)
(
  cd "$ROOT"
  OMP_NUM_THREADS=1 "$OUT/prem-new" modeldata/PREM500 > "$OUT/PREM500-new.csv"
  OMP_NUM_THREADS=1 "$OUT/prem-new" modeldata/prem.200 > "$OUT/prem.200-new.csv"
)
cmp "$OUT/PREM500-old.csv" "$OUT/PREM500-new.csv"
cmp "$OUT/prem.200-old.csv" "$OUT/prem.200-new.csv"
sha256sum "$OUT/PREM500-old.csv" "$OUT/PREM500-new.csv" \
  "$OUT/prem.200-old.csv" "$OUT/prem.200-new.csv"
