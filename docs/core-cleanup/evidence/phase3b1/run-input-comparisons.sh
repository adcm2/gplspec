#!/usr/bin/env bash
set -euo pipefail

BASE=/tmp/gplspec-phase3-baseline
BASE_BUILD=/tmp/gplspec-phase3-baseline-release
CANDIDATE=/home/adcm2/Documents/Research/Projects/gplspec
CANDIDATE_BUILD=/tmp/gplspec-phase3-candidate-release
EVIDENCE="$(cd "$(dirname "$0")" && pwd)"
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
  local source="$1" root="$2" build="$3" output="$4" typed="$5"
  local dep_args=()
  for arg in "${DEPS[@]}"; do
    dep_args+=("${arg//_DEP_/$build/_deps}")
  done
  local typed_arg=()
  if [[ "$typed" == yes ]]; then
    typed_arg=(-DGPLSPEC_TYPED_INPUT)
  fi
  /usr/bin/c++ -O0 -DNDEBUG -std=gnu++23 -fopenmp "${typed_arg[@]}" \
    -I"$root" "${dep_args[@]}" "$EVIDENCE/$source" -o "$output" \
    -lfftw3 -lfftw3f -lfftw3l -lnetcdf
}

compile_probe prem-model-probe.cpp "$BASE" "$BASE_BUILD" \
  /tmp/gplspec-phase3b-model-probe-old no
compile_probe prem-model-probe.cpp "$CANDIDATE" "$CANDIDATE_BUILD" \
  /tmp/gplspec-phase3b-model-probe-new no
compile_probe phobos-reader-probe.cpp "$BASE" "$BASE_BUILD" \
  /tmp/gplspec-phase3b-phobos-probe-old no
compile_probe phobos-reader-probe.cpp "$CANDIDATE" "$CANDIDATE_BUILD" \
  /tmp/gplspec-phase3b-phobos-probe-new yes

(
  cd "$BASE"
  OMP_NUM_THREADS=1 /tmp/gplspec-phase3b-model-probe-old modeldata/PREM500 \
    > /tmp/gplspec-phase3b-PREM500-old.csv
  OMP_NUM_THREADS=1 /tmp/gplspec-phase3b-model-probe-old modeldata/prem.200 \
    > /tmp/gplspec-phase3b-prem200-old.csv
  OMP_NUM_THREADS=1 /tmp/gplspec-phase3b-phobos-probe-old \
    > /tmp/gplspec-phase3b-phobos-old.csv
)
(
  cd "$CANDIDATE"
  OMP_NUM_THREADS=1 /tmp/gplspec-phase3b-model-probe-new modeldata/PREM500 \
    > /tmp/gplspec-phase3b-PREM500-new.csv
  OMP_NUM_THREADS=1 /tmp/gplspec-phase3b-model-probe-new modeldata/prem.200 \
    > /tmp/gplspec-phase3b-prem200-new.csv
  OMP_NUM_THREADS=1 /tmp/gplspec-phase3b-phobos-probe-new \
    > /tmp/gplspec-phase3b-phobos-new.csv
)

cmp /tmp/gplspec-phase3b-PREM500-old.csv /tmp/gplspec-phase3b-PREM500-new.csv
cmp /tmp/gplspec-phase3b-prem200-old.csv /tmp/gplspec-phase3b-prem200-new.csv
cmp /tmp/gplspec-phase3b-phobos-old.csv /tmp/gplspec-phase3b-phobos-new.csv
sha256sum /tmp/gplspec-phase3b-PREM500-old.csv \
  /tmp/gplspec-phase3b-prem200-old.csv /tmp/gplspec-phase3b-phobos-old.csv
