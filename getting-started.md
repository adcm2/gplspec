---
title: "Getting Started"
permalink: /getting-started/
---

# Getting started with gplspec

gplspec is a header-only C++ library. CMake supplies the pinned header dependencies and propagates the include and link requirements through the `gplspec` interface target. A supported CMake consumer must use C++23.

## Requirements

- A C++23 compiler. The cleanup validation used GNU C++ 13.3.0.
- CMake 3.28.3 for the recorded build and test workflow (the project declares a lower minimum, but older CMake versions are not covered by this validation).
- Git for CMake FetchContent source dependencies.
- OpenMP development support.
- FFTW3 development libraries, including `fftw3` and `fftw3f`, plus headers (`fftw3l` is detected and linked when available).
- NetCDF development headers and library, required by the pinned model dependencies.
- Python 3 and Bash for the registered regression tests.

The stage-00 baseline's recorded system packages were Eigen 3.4.0, FFTW 3.3.10, and NetCDF 4.9.2. Stage 07 used Eigen 3.4.0 from the SHA-256-pinned source archive; a parent project may instead supply `Eigen3::Eigen`, which GPLSpec reuses. The recorded FFTW and NetCDF versions are test provenance, not CMake-enforced requirements. The seven Git dependencies are pinned to commits in [the baseline manifest](docs/cleanup/baseline-manifest.md).

## Build the examples

```sh
git clone https://github.com/adcm2/gplspec.git
cd gplspec
git checkout b813285fc34b13e9d4868860ae0651eb405457f4
cmake -S . -B build -DMY_PROJECT_BUILD_EXAMPLES=ON
cmake --build build -j2
```

The executable is written to `build/bin/clean_bench_1`. Its working directory matters: this example writes under `./work/Bench1`. Run it from a disposable directory so its scientific output does not replace repository data:

```sh
mkdir -p /tmp/gplspec-bench1-run/work/Bench1
cd /tmp/gplspec-bench1-run
/absolute/path/to/gplspec/build/bin/clean_bench_1
```

## Run the cleanup regression suite

The opt-in baseline harness enables the tests and their registered CTest checks:

```sh
cd /absolute/path/to/gplspec
cmake -S . -B /tmp/gplspec-build \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON
cmake --build /tmp/gplspec-build -j2
ctest --test-dir /tmp/gplspec-build --output-on-failure
tests/run_stage00.sh /tmp/gplspec-build/bin/stage00_reference
/tmp/gplspec-build/bin/stage02_header_link
```

The frozen stage-00 runner executes the candidate twice. It compares each 20,700-record output with the checked-in original-source reference using `rtol=1e-13` and `atol=1e-15`, compares the two runs to each other, and checks each representative output byte-for-byte. Never regenerate the frozen references as part of ordinary validation.

## Use gplspec from another CMake project

Link the `gplspec` interface target so its C++ dependencies and include requirements propagate. Set C++23 on the consumer target. The current project calls `find_package(FFTW)` from its top-level CMake file, so a FetchContent parent must also expose gplspec's CMake module directory before making the dependency available. A parent-provided `Eigen3::Eigen` target is reused; otherwise gplspec fetches Eigen 3.4.0.

```cmake
cmake_minimum_required(VERSION 3.28)
project(my_consumer LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(gplspec
  GIT_REPOSITORY https://github.com/adcm2/gplspec.git
  GIT_TAG b813285fc34b13e9d4868860ae0651eb405457f4)
FetchContent_GetProperties(gplspec)
if(NOT gplspec_POPULATED)
  FetchContent_Populate(gplspec)
endif()
list(APPEND CMAKE_MODULE_PATH "${gplspec_SOURCE_DIR}/cmake")
if(NOT TARGET gplspec)
  add_subdirectory("${gplspec_SOURCE_DIR}" "${gplspec_BINARY_DIR}")
endif()

add_executable(my_consumer main.cpp)
target_link_libraries(my_consumer PRIVATE gplspec)
target_compile_features(my_consumer PRIVATE cxx_std_23)
```

The example pins the accepted stage-06 base used for stage 07, since the normal `main` branch does not yet contain this cleanup baseline. The parent remains responsible for providing the system FFTW and NetCDF development packages. The compiler include paths for any additional public headers should come from linking the `gplspec` target rather than manually listing dependency directories.

## Known scope

The cleanup preserves existing angle and density conventions, including the slice-rotation sign and one-sided epsilon checks. It also preserves the small Hermitian-form discrepancies recorded in the deferred-issues file. Validation covers finite deterministic fixtures and the recorded dependency/compiler stack; it does not establish behavior for every model, degree, or platform. Experimental sources under `experimental/` are preserved, excluded from production targets, and numerically unvalidated.
