---
title: "Getting Started"
permalink: /getting-started/
---

# Getting started

`gplspec` is a header-only C++ library. Link the CMake interface target named `gplspec`; it supplies the public include paths, pinned header dependencies, and C++23 compile feature. The project currently declares CMake 3.20 as its minimum. There is no install/export package yet, so add the checkout as a subdirectory.

## Build examples and regression checks

Requirements are CMake 3.20 or newer, a C++23 compiler, Git, OpenMP development support, FFTW3 development files (`fftw3` and `fftw3f`; `fftw3l` is used when available), NetCDF development files, Python 3, and Bash. CMake fetches the Git-pinned header dependencies. Eigen 3.4.0 is fetched unless the parent project already defines `Eigen3::Eigen`.

The recorded Phase 2 validation used CMake 3.28.3, GNU C/C++ 13.3, Eigen 3.4.0, FFTW 3.3.10, and NetCDF 4.9.2. Those are evidence provenance, not enforced exact versions beyond the pinned source revisions.

These commands describe the current working tree. Until the Phase 2 changes are published, a fresh clone of the public default branch may not contain this target-level build contract.

```sh
cd /path/to/current/gplspec-checkout
cmake -S . -B build -DMY_PROJECT_BUILD_EXAMPLES=ON
cmake --build build -j2
```

To configure the registered test harness as well:

```sh
cmake -S . -B build-tests \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON
cmake --build build-tests -j2
ctest --test-dir build-tests --output-on-failure
```

The fixed Stage 00 runner and header-link smoke check are separate commands:

```sh
tests/run_stage00.sh /absolute/path/to/gplspec/build-tests/bin/stage00_reference
/absolute/path/to/gplspec/build-tests/bin/stage02_header_link
```

The Stage 00 runner compares three 20,700-record outputs with frozen references, using `rtol=1e-13` and `atol=1e-15`, and checks representative output bytes. Do not regenerate those references during ordinary validation. The regression suite is bounded evidence for its listed fixtures; it does not validate every example, model, or platform. See [workflow and output contracts](benchmark.md).

## Consume the current checkout from another CMake project

This recipe uses the current checkout as a subdirectory. The consumer links only `gplspec`; the target propagates C++23 and GPLSpec's dependency include/link requirements. GPLSpec's top-level CMake handles its own FFTW module lookup, so the parent does not edit `CMAKE_MODULE_PATH`.

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_consumer LANGUAGES CXX)

add_subdirectory("/absolute/path/to/gplspec" gplspec-build)
add_executable(my_consumer main.cpp)
target_link_libraries(my_consumer PRIVATE gplspec)
```

CMake fetches GPLSpec's pinned sources unless compatible dependency targets/sources are already provided through the supported parent configuration. The tracked [external consumer check](tests/external_consumer/CMakeLists.txt) exercises subdirectory integration with parent-provided Eigen and checks that the target supplies language and module requirements.

The old website recipe pinned immutable PR #7 commit `a77742f60fee305a1e058663e55e5ca62c2f3739`; that revision predates the current target propagation. It is historical and is not a recipe for this working tree.

## Run an example safely

Executables use paths relative to their process working directory. For example, `clean_bench_1` writes to `work/Bench1`; create that directory and run from a disposable working directory:

```sh
mkdir -p /tmp/gplspec-run/work/Bench1
cd /tmp/gplspec-run
/absolute/path/to/gplspec/build/bin/clean_bench_1
```

The examples do not create their output directories. File-backed examples also expect the tracked `modeldata/` paths relative to the working directory; arrange those inputs in an isolated run directory before launching. Some examples read values from stdin. The [workflow map](benchmark.md) lists each program's prompts, inputs, output paths, and plot consumers. Building an executable does not run or validate its scientific workflow.
