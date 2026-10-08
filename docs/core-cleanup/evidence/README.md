# Stage 1C evidence

Candidate evidence from the fresh Release build `/tmp/gplspec-core-phase1-build` on Linux x86_64. Toolchain: CMake 3.28.3, GNU C/C++ 13.3.0, Eigen 3.4.0 (fetched pinned archive), FFTW 3.3.10, NetCDF 4.9.2, OpenMP 4.5, Python 3.12.3. Configuration used Release, Ninja, examples enabled, baseline harness enabled, and `-j2` builds. The configure log records the existing CMP0148 developer warning from `FindPythonInterp`.

Commands and saved outputs:

- Configure: `cmake -S . -B /tmp/gplspec-core-phase1-build -G Ninja -DCMAKE_BUILD_TYPE=Release -DMY_PROJECT_BUILD_EXAMPLES=ON -DGPLSPEC_BUILD_BASELINE_HARNESS=ON`; `stage1c-configure.log`.
- Initial full build: `cmake --build /tmp/gplspec-core-phase1-build -j2`; `stage1c-initial-full-build.log` (64 build steps, all examples and then-existing tests). The two new test targets were subsequently built; the final full build reported no remaining work.
- Final tests: `ctest --test-dir /tmp/gplspec-core-phase1-build --output-on-failure`; `stage1c-ctest.log` (14/14 passed, including the standalone external consumer).
- Standalone public umbrella compile checks: the existing `Core`, `GeneralModels`, `SimpleModels`, and `All` targets plus new `Test` and `Timer` targets. The initial Luna review found the latter two missing. After adding them, both were force-recompiled with `cmake --build /tmp/gplspec-core-phase1-build --target stage02_test_smoke stage02_timer_smoke -j2`; `stage1c-umbrella-followup.log` records both successful compilations.
- Frozen numerical references: `tests/run_stage00.sh /tmp/gplspec-core-phase1-build/bin/stage00_reference`; `stage1c-stage00.log` (three 20,700-record comparisons, zero changes, fixed `rtol=1e-13`, `atol=1e-15`, output byte-identical).
- Solver diagnostic: `/tmp/gplspec-core-phase1-build/bin/stage1c_file_tomography_solver`; `stage1c-solver-diagnostic.log` (3 iterations at `1e-6`; solver error `1.62767e-7`; absolute residual `3.04753e-5`; relative residual `1.62767e-7`; public formatted result difference `0`).
- External consumer direct run: `python3 tests/run_external_consumer.py /home/adcm2/Documents/Research/Projects/gplspec /tmp/gplspec-core-phase1-build/tests/external-consumer-build /tmp/gplspec-core-phase1-build/_deps /tmp/gplspec-core-phase1-build/_deps/eigen3-src`; `stage1c-external-consumer.log` (configure, build, run passed).

The direct inputs and unchanged frozen reference hashes are:

| File | SHA-256 |
| --- | --- |
| `modeldata/PREM500` | `ec9869508847eb6753cf0165f6176dd3eea1f03ac5717533ffb6dcfe09489685` |
| `modeldata/S40RTS_dvs.nc` | `4ebd0ca09731a593c7960ced6a8ddde0abbe03fa91e49e70cf897dd41019c73e` |
| `tests/reference/stage00-baseline.csv` | `ef20bd7e3016661c60903c290c74c599151fdfc0417443bcf81c4215401f105e` |
| `tests/reference/homogeneous-MatrixSolution.out` | `a902254d3d416c55f103058c045a30036e7b5edd8b56166e086d43de4de13c3f` |
| `tests/reference/stage05-1d-original.csv` | `a9e89b9ac1250f3f01b695ec5be0fcb79990a85da42e6eb80ae001f5f88a43d6` |
| `tests/reference/stage05-file-perturbation-original.csv` | `08cf04f1243fa128645490c4eba13eed5f7b5fceff76b032f2f5d8141efe3263` |
| `tests/reference/stage06-output-baseline.sha256` | `f58a90f7f05c1f6e88712a8d6b20f4202284cb2bb266782961e1a018c8b22236` |

The root CMake dependency pins were not modified: NumericConcepts `888126b44a979fda5dcc9d0b11d80cb8b3f61b52`; FFTWpp `06cc1fb04c4398407839e637cb5414b179b0e9c3`; GaussQuad `fbe37c7eef93695317dbf9c48715110bf64c800f`; Interpolation `1557aad571ceac1899c2d3795696aa4047f41d8a`; PlanetaryModel `a801b9e9c2205c3b32daa745fd523d596a9d2ede`; GSHTrans `f0a0e24a3579ceb7c99f8782328cae14bef195de`; TomographyModels `97869aeea3f3d87901705db59b5de975960bd9e1`. Eigen is the existing 3.4.0 archive with the root-pinned SHA-256 recorded in `CMakeLists.txt`.

`stage1c-content-manifest.md` lists 96 candidate source/build/CI content hashes against the exact Stage 1C parent and separately hashes tracked test references plus the two file-backed tomography inputs. It excludes status/review records and itself to avoid recursive hashes.

The saved configure transcripts `stage1c-configure.log` and `stage1c-external-consumer.log` have trailing spaces removed from tool-generated report lines to satisfy the staged whitespace check. This formatting-only cleanup changes no command, version, exit status, or numerical result.

The CI workflow has been added but has not run on GitHub. No performance measurements were made; those belong to Stage 1D. No solver algorithm, requested tolerance, reference, input or dependency was changed.
