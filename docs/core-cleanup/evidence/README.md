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

## Stage 1D performance baseline

The opt-in `phase1_performance` executable uses `std::chrono::steady_clock` and is not registered as a timing CTest. It measures a mapped 3-iteration solve fixture and the file-backed `PREM500`/`S40RTS_dvs.nc` tomography path. The file case separately times Earth-model input read, NetCDF read, model/discretization construction and one complex matrix-vector action. The mapped case separately times base-model setup, model/discretization construction, first and repeated matrix-vector actions, then first and repeated public `FindGravitationalPotential` calls. Each full public-solver timing includes operator, force, preconditioner, solve, result-format setup, harness `CoutCapture` construction, diagnostic-text capture and string copy; each call creates a fresh solver. The first public solve occurs after model/operator setup and six matrix-vector actions, so it is not a cold end-to-end or cold-OS solve. “Repeat” means same-process repeated calls, not reuse of solver state. The harness independently recomputes residual and captures reported diagnostics once per process after the timed public calls; residual is not recomputed for every timed call. The harness checks all public solver calls retain the three-iteration result.

The five-process sample set is generated by `python3 tests/run_phase1_performance.py /tmp/gplspec-core-phase1-build/bin/phase1_performance docs/core-cleanup/evidence/stage1d-raw`. The driver serially sets `OMP_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, `OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`, saves each process's CSV, checks the 30 public solver calls and five independent residuals, and writes `summary.json` with median, min/max, CV, MAD and `median + 3 MAD` values. Each process performed five repeated matvecs and five repeated complete public solves after the first call. Inputs, resolutions, and paths are fixed to the Stage 1C case; their hashes are in `stage1c-content-manifest.md`. Machine: Linux x86_64, Intel Core Ultra 5 235U (14 logical CPUs reported), kernel `7.0.0-38-generic`. Toolchain: CMake 3.28.3, GNU C++ 13.3.0, Eigen 3.4.0, FFTW 3.3.10, NetCDF 4.9.2, OpenMP 4.5. CMake Release was used (`-O3 -DNDEBUG`), with the pinned dependency sources from Stage 1C.

| Metric | Samples | Median | Observed range | CV | Provisional local threshold |
| --- | ---: | ---: | ---: | ---: | --- |
| EarthModel input read | 5 processes | 0.001496 s | 0.001112–0.004241 s | 59.6% | unavailable: high variability |
| NetCDF input read | 5 processes | 0.13575 s | 0.12113–0.15152 s | 7.7% | 0.16787 s (`median + 3 MAD`); local reference only |
| File-model discretization/setup | 5 processes | 0.000705 s | 0.000668–0.000802 s | 6.5% | unavailable: max exceeds `median + 3 MAD` |
| File-tomography first matvec | 5 processes | 0.000594 s | 0.000561–0.000728 s | 9.8% | unavailable: max exceeds `median + 3 MAD` |
| Mapped model discretization/setup | 5 processes | 0.001396 s | 0.001088–0.002269 s | 27.4% | unavailable: high variability |
| Mapped repeated matvec | 25 calls | 0.000204 s | 0.000114–0.000348 s | 33.2% | unavailable: observed timing noise is too large |
| First mapped full public solve | 5 processes | 0.001043 s | 0.000657–0.001578 s | 30.5% | unavailable: observed timing noise is too large |
| Repeated mapped full public solve | 25 calls | 0.000795 s | 0.000558–0.001269 s | 23.0% | unavailable: observed timing noise is too large |
| Process peak RSS | 5 processes | 250,420 KiB | 250,140–250,532 KiB | 0.054% | 275,585 KiB (10% above largest sample); same host only |

Every measured mapped public solve used 3 iterations. Independent residual stayed at relative `1.62767e-7` and absolute `3.04753e-5`; Eigen-reported error was `1.62767e-7`. These are correctness diagnostics, not timing outcomes. Fine-grained mapped timings are too noisy for thresholds on this host; do not gate CI on them. NetCDF read and same-host RSS thresholds are provisional references only, not cross-machine guarantees.

One clean-first Release build measured all configured targets (70 build steps) with `-j2` in the already-configured build tree and populated dependency directories. Command: `OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/time -v -o docs/core-cleanup/evidence/stage1d-clean-full-time.txt cmake --build /tmp/gplspec-core-phase1-build --clean-first -j2 > docs/core-cleanup/evidence/stage1d-clean-full-build.log 2>&1`. `/usr/bin/time -v` recorded 3:41.47 wall time and 980,376 KiB peak RSS; all targets built successfully. This clean includes examples and test targets, and excludes configure and dependency fetch. It is one observation with no provisional regression threshold. Raw logs: `stage1d-clean-full-build.log` and `stage1d-clean-full-time.txt`.

Harness SHA-256: `eab5227292ec89f3edcc642c9698daab0c62ccb442fb8daa06b85c6c0a214334`; five-process driver SHA-256: `177bb3e947a9e9a4a3e0431a509e9737a25db1f8cb38fe615bbc9d87362f4b0c`. Benchmark target, raw samples, generated aggregate, build-timing logs, and the measured methodology are candidate evidence; no algorithm, input, tolerance, frozen reference, or dependency was modified.
