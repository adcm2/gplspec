# Stage-00 baseline manifest

## Source and toolchain

- Repository: `https://github.com/adcm2/gplspec.git`; origin `main` and `develop` were unchanged.
- Source baseline: `4ef3a66c62d64408d99989dd51c3ccbdc46d0b0c`; annotated tag `gplspec-cleanup-start` already pointed to this commit.
- Compiler: `/usr/bin/c++`, GNU 13.3.0 (Ubuntu `13.3.0-6ubuntu2~24.04.1`). CMake 3.28.3. OpenMP 4.5.
- Exact system packages: Eigen 3.4.0 (`libeigen3-dev 3.4.0-4build0.1`), FFTW 3.3.10 (`libfftw3-dev 3.3.10-1ubuntu3`), NetCDF 4.9.2 (`libnetcdf-dev 1:4.9.2-5ubuntu4`).

## Locked source dependencies

All local source checkouts below were clean at the recorded revision. Root CMake declares the complete closure before population so transitive moving-branch declarations cannot win.

| Dependency | Repository | Commit |
|---|---|---|
| GSHTrans fork | `https://github.com/adcm2/GSHTrans.git` | `f0a0e24a3579ceb7c99f8782328cae14bef195de` |
| NumericConcepts | `https://github.com/da380/NumericConcepts.git` | `888126b44a979fda5dcc9d0b11d80cb8b3f61b52` |
| GaussQuad | `https://github.com/da380/GaussQuad.git` | `fbe37c7eef93695317dbf9c48715110bf64c800f` |
| FFTWpp | `https://github.com/da380/FFTWpp.git` | `06cc1fb04c4398407839e637cb5414b179b0e9c3` |
| Interpolation | `https://github.com/da380/Interpolation.git` | `1557aad571ceac1899c2d3795696aa4047f41d8a` |
| PlanetaryModel | `https://github.com/da380/PlanetaryModel.git` | `a801b9e9c2205c3b32daa745fd523d596a9d2ede` |
| TomographyModels | `https://github.com/adcm2/TomographyModels.git` | `97869aeea3f3d87901705db59b5de975960bd9e1` |

The Eigen source clone was not used: a fresh 3.4.1 GitLab clone did not complete in a reasonable time. The validated stack uses the exact installed Eigen 3.4.0 package; `find_package(Eigen3 3.4.0 EXACT)` and a small FetchContent shim ensure nested dependencies use the same package.

## Build and regression commands

```sh
cmake -S . -B /tmp/gplspec-cleanup-baseline \
  -DMY_PROJECT_BUILD_EXAMPLES=OFF \
  -DGPLSPEC_BUILD_BASELINE_HARNESS=ON
cmake --build /tmp/gplspec-cleanup-baseline --target stage00_reference -j2
tests/run_stage00.sh /tmp/gplspec-cleanup-baseline/bin/stage00_reference
```

These exact commands were run with CMake FetchContent downloading the locked SHA revisions (no external source overrides). Configure, build, and both regression runs passed. The checker requires matching keys and finite values and uses `rtol=1e-13`, `atol=1e-15`. The output artifact is checked byte-for-byte.

## Immutable reference provenance

- Reference source worktree: `/tmp/gplspec-stage00-base-ref`, detached at the original source commit above.
- Reference executable: `/tmp/gplspec-stage00-base-build/bin/stage00_reference`; SHA-256 `54bb65502536691aa2ceb7414ba8a1acbc3add8bf1122229f38878b8b58fd3ba`.
- The reference worktree received only temporary CMake glue to link GSHTrans and add the baseline harness target. No numerical source files were changed. It used the exact clean dependency revisions and toolchain above.
- Frozen CSV: `tests/reference/stage00-baseline.csv`, 20,700 records, SHA-256 `ef20bd7e3016661c60903c290c74c599151fdfc0417443bcf81c4215401f105e`.
- Frozen existing-format output: `tests/reference/homogeneous-MatrixSolution.out`, 329 bytes, SHA-256 `a902254d3d416c55f103058c045a30036e7b5edd8b56166e086d43de4de13c3f`.
- Two original-source runs and one candidate-source run were byte-identical: max absolute difference 0 across all recorded values.
