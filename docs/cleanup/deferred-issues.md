# Deferred numerical and build observations

These are recorded observations or validation boundaries. Stage 07 documents them; it does not investigate or repair them.

## Preserved numerical behavior

- The original-source stage-00 baseline reports zero iterations for the homogeneous sphere, layered sphere, and lateral-density solve at relative tolerance `1e-6`; the smooth mapping solve reports one iteration, the mapping-perturbation base solve one, and the perturbation solve zero. Fixture force/operator vectors and returned values are finite. These iteration counts are preserved observations, not claims about arbitrary inputs.
- Original-source stage-05 probes retain small nonzero differences between complex sesquilinear forms. In 1D, `<x,Ay>` is `302694.90083236271 - 302420.33613110462i` and `<Ax,y>` is `302694.86333236267 - 302420.37363110413i`. For the 3D file perturbation, `<x,Ay>` is `1967.8406259289259 - 1923.06562591585i` and `<Ax,y>` is `1967.8406259037456 - 1923.0656259410546i`. The cleanup records these values and imposes no symmetry assumption or correction.
- Output-angle conventions remain caller-specific. In particular, `RotateSliceToEquator` keeps its opposite `tmp4` sign and its one-sided epsilon checks. Existing density-normalization behavior, including caller-specific physical versus referential density selection, is also preserved. These quirks are not reinterpreted or repaired by the cleanup.

## Coverage limits

The baseline and stage-05 operator evidence use finite deterministic fixtures, including a degree-2 coefficient perturbation. They do not establish exhaustive behavior for arbitrary high degrees, every model, or every platform. Stage-06 output comparisons cover the named writers and fixtures on the recorded toolchain. They are behavior-preservation evidence within that coverage, not a general scientific validation of GPLSpec.

All files under `experimental/`, including the historical ellipticity calculations, remain outside production targets and have not been automatically built or numerically validated.

## Dependency and consumer provenance

The baseline records GNU C++ 13.3.0, CMake 3.28.3, Eigen 3.4.0, FFTW 3.3.10, and NetCDF 4.9.2. Eigen 3.4.0 is fetched from a SHA-256-pinned archive unless the parent project supplies `Eigen3::Eigen`. The seven Git source dependencies are pinned in [baseline-manifest.md](baseline-manifest.md). CMake does not enforce the recorded system FFTW and NetCDF versions; other versions may configure but are outside the recorded validation stack.

An embedded CMake consumer must set C++23 on its target, provide the FFTW and NetCDF system development dependencies, and add the fetched GPLSpec `cmake/` directory to `CMAKE_MODULE_PATH` before configuring GPLSpec so `find_package(FFTW)` resolves the repository's `FindFFTW.cmake`. Linking the `gplspec` interface target supplies its transitive library/header dependencies. A parent-provided Eigen target is reused; otherwise GPLSpec fetches its pinned Eigen archive. These consumer constraints are documented rather than changed in build code.
