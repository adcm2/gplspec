---
title: "Getting Started"
permalink: /getting-started/
---

# Getting started

`gplspec` is a header-only C++ library. Link the CMake interface target named `gplspec`; it supplies the public include paths, pinned header dependencies, and C++23 compile feature. The project currently declares CMake 3.20 as its minimum. There is no install/export package yet, so add the checkout as a subdirectory.

## Build examples and regression checks

Requirements are CMake 3.20 or newer, a C++23 compiler, Git, OpenMP development support, FFTW3 development files (`fftw3` and `fftw3f`; `fftw3l` is used when available), NetCDF development files, Python 3, and Bash. CMake fetches the Git-pinned header dependencies. Eigen 3.4.0 is fetched unless the parent project already defines `Eigen3::Eigen`. NumPy and Matplotlib are needed for the registered field-bundle and plotting tests and for `work/plot_fields.py`; they are Python workflow dependencies and are not linked by the header-only C++ library.

The recorded Phase 2 validation used CMake 3.28.3, GNU C/C++ 13.3, Eigen 3.4.0, FFTW 3.3.10, and NetCDF 4.9.2. Those are evidence provenance, not enforced exact versions beyond the pinned source revisions.

Use these commands to configure and build the checkout at its current revision. For workflow settings and expected outputs, see [scientific workflows](benchmark.md); for isolated producer and plotting commands, see the [field-bundle guide](work/plots.md).

```sh
cd /path/to/current/gplspec-checkout
cmake -S . -B build -DMY_PROJECT_BUILD_EXAMPLES=ON
cmake --build build -j2
```

To configure the registered test harness, install its plotting/test dependencies into a virtual environment and point CMake at that same interpreter. This also supplies the interpreter CMake uses for its Python-based harness scripts:

```sh
python3.12 -m venv build/python-test-env
build/python-test-env/bin/python -m pip install -r tests/requirements.txt
cmake -S . -B build-tests \
  -DMY_PROJECT_BUILD_EXAMPLES=ON \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON \
  -DPYTHON_EXECUTABLE="$PWD/build/python-test-env/bin/python"
cmake --build build-tests -j2
ctest --test-dir build-tests --output-on-failure
```

The same Python 3.12 virtual environment can run the documented plotting commands. C++ consumers that only add GPLSpec as a CMake subdirectory and link `gplspec` do not need NumPy or Matplotlib.

The Stage 00 scientific runner and standalone header-link smoke check can also be run directly:

```sh
tests/run_stage00.sh /absolute/path/to/gplspec/build-tests/bin/stage00_reference
/absolute/path/to/gplspec/build-tests/bin/stage02_header_link
```

The Stage 00 runner takes two fresh captures and gates all 18,648 retained scientific records against the unchanged 20,700-record frozen reference. It explicitly retires only the 2,052 intermediate `dxi`, `dxilm`, and `df` records. The runner also prints the original strict exact comparison and representative byte comparison as diagnostics; those diagnostics can report floating-point or serialized-output differences while the separately defined scientific criteria pass. Stage 05 source/action and file-backed field gates run under CTest. Do not regenerate frozen references during ordinary validation. These checks are bounded evidence for their listed fixtures; they do not validate every example, model, or platform. See [workflow and output contracts](benchmark.md).

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

## Run an example safely

Executables use paths relative to their process working directory. For example, `clean_bench_1` writes to `work/Bench1`; create that directory and run from a disposable working directory:

```sh
mkdir -p /tmp/gplspec-run/work/Bench1
cd /tmp/gplspec-run
/absolute/path/to/gplspec/build/bin/clean_bench_1
```

Field-bundle writers create their output directories. File-backed examples expect the tracked `modeldata/` paths relative to the working directory; link or copy those inputs into an isolated run directory before launching. Some examples read values from stdin. The [workflow map](benchmark.md) lists each program's prompts, inputs, and effective settings. The [bundle plotting guide](work/plots.md) documents current bundle paths and plot commands. Building an executable does not run or validate its scientific workflow.
